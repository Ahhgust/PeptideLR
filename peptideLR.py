#!/usr/bin/python3
# Written by August Woerner

# MIT License

# Copyright (c) [2019] [August E. Woerner]

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


from math import log, exp, factorial
from collections import namedtuple
from collections import defaultdict
import sys
import csv
import gzip
from itertools import combinations, combinations_with_replacement, accumulate, product, chain
import random
import argparse
import os.path
import suffer
import numpy
import glob
# 0.0011 -> added cache for peptide hits (keep last query)
# 0.0012 -> added cache for peptide hits (keep all queries) AND added cache for keeping likelihoods (-N argument)
# 0.0013 -> added weights (-W) and population filtering
# 0.0014 -> added RMP
# 0.0015 -> modified to accommodate new suffix array (suffer.py) interface.
# 0.0016 -> decoupled likelihoods in the LR presentation 
VERSION_NUM = "0.0016"

# A glorified named tuple
# records the program state 
class State(namedtuple("State", ["experimentType", "Population", "randomSeed"])):
  __slots__ = ()

  def __str__(self):

    global LOGLIKES

    likeStr="\tNone"
    if LOGLIKES:
      likeStr="\tLoglikelihood"
    pop = self.Population
    if pop == "":
      pop= "Pooled"
      
    global IS_RANDOMIZED
    if not IS_RANDOMIZED:
      return self.experimentType + "\t" + pop + "\t" + VERSION_NUM + likeStr
    return self.experimentType + "\t" + pop + "\t" + VERSION_NUM + "\t" + str(self.randomSeed) + likeStr
  
  def keyHeader(self):

    likeStr="\tLikelihoodTransform"

    
    global IS_RANDOMIZED
    if not IS_RANDOMIZED:
      return "experimentType\tPopulation\tversionNumber"+likeStr
    return "experimentType\tPopulation\tversionNumber\trandomSeed"+likeStr



  
# suffixes added to diploid IDs (e.g., NA12887)
# to make them haploid (into NA12887_1)
MATERNAL_SUFFIX = "_1"
PATERNAL_SUFFIX = "_2"

def unknownHaplotypeGenerator(nHaps, nPeeps, nMc=0, cumWeights=None):
  """
  We need to walk through every possible diploid combination of haplotypes
  from some fixed number of individuals
  nHaps is the number of haplotypes possible
  nPeeps is the number of (unknown) individuals to consider
  *
  alternatively, generate Monte Carlo samples
  (weighted by the probability of the haplotypes assuming HWE)

  returns a flattened list of haplotype identifiers; ints in range(nHaps)
  consecutive pairs represent the haplotype-pair for an unknown individual
  """

  if nMc < 1:
    
    # possible pairs of haplotypes (1 person)
    diplotypes = combinations_with_replacement(range(nHaps), 2)
    # if we were more careful with our counting we could make the below work...
    #allPeeps = combinations_with_replacement(diplotypes, nPeeps)

    # but the cartesian product is what we actually want.
    allPeeps = product(diplotypes, repeat=nPeeps)

    for mix in allPeeps:
      # and flatten the list...
      # (and make it into a list as opposed to an iterable)
      yield list(chain.from_iterable(mix))
      
  else:
    if cumWeights is None:
      print("I need cumulative weights to work!", file=sys.stderr)
      return None
    
    intList = list(range(nHaps))
    n = nPeeps + nPeeps # this is how many haplotypes we want to sample
    for i in range(nMc):
      yield random.choices(intList, cum_weights=cumWeights, k=n)
      
def diploidIDsToHaploid(ids):
  """
  Convenience function
  converts diploid ids (of some person)
  and gives us the IDs associated with the maternal and paternal chromosomes
  """

  # scalar implementation
  if type(ids) == str:
    return [ids + MATERNAL_SUFFIX, ids + PATERNAL_SUFFIX ]

  # and vector (or any other iterable other than strings...)
  ret = []
  for who in ids:
    ret.append(who + MATERNAL_SUFFIX)
    ret.append(who + PATERNAL_SUFFIX)
    
  return ret


def parseSamps2Pops(filename, sampleCol, popCol):
  """
  This parses a TSV file of samples associated with their respective population and 
  returns 2 dictionaries; one associating a sample -> their pop (samps2pops)
  and another associating the pop -> the list of samples (pops2samps)
  e.g.,
  a sample ID (NA00015) is they to their corresponding population affiliation (dictionary)
  and it returns a population (TSI) to a list of IDs associated with that population
  
  The arguments needed is the filename to parse
  The file needs to be a TSV file needs with a header

  and columns are the *names* of the columns (args2 and 3)
  associated with the sample ID and the population IDs, respectively
  """

  samps2pops = {}
  pops2samps = defaultdict(list)
  with open(filename) as tsv:
    reader = csv.DictReader(tsv, delimiter="\t")
    for row in reader:
      if sampleCol not in row or popCol not in row:
        print(filename , " does not have the appropriate column names in it. expecting " , sampleCol , "\t" , popCol, file=sys.stderr)
        sys.exit()
      # associate a sample ID to their population
      # and a population to a list of sample IDs
      samps2pops[ row[sampleCol] ] = row[ popCol ]
      pops2samps[ row[ popCol ] ].append( row[ sampleCol])


  return (samps2pops, pops2samps)
      
def isoleucine2leucine(s):
  """
  converts all isoleucines to leucines in a string...
  this is done to remove equivalent weight amino acids
  it serves to reduce the alphabet from 20 aa to 19.
  The I->L reduction is pre-performed in the database (suffix array)
  """

  ss = list(s.upper())

  for i in range(len(ss)):
    if ss[i] == 'I':
      ss[i] = 'L'

  return "".join(ss)


def homozygousProbWithTheta(theta, pu):
  """
  Notation from Baldings and Nichols (1994)
  taken from Buckleton FSI Gen (23) 2016; 91-100
  https://doi.org/10.1016/j.fsigen.2016.03.004
  """
  # 2theta + (1-theta)*pu
  numerator1 = (theta+theta) + (1-theta)*pu
  # 3theta + (1-theta)*pu == theta + numerator1
  numerator2 = theta + numerator1
  
  # 1 + theta TIMES 1 + 2theta
  denom = (1+theta) * (1 + theta + theta)

  return (numerator1*numerator2)/denom

def heterozygousProbWithTheta(theta, pu, pv):
  """
  Notation from Baldings and Nichols (1994)
  taken from Buckleton FSI Gen (23) 2016; 91-100
  https://doi.org/10.1016/j.fsigen.2016.03.004
  """

  # theta + (1-theta)pu
  numerator1 = (theta + (1-theta)*pu)
  numerator2 = (theta + (1-theta)*pv)

  denom = (1+theta)*(1+theta+theta)

  return (2*numerator1*numerator2)/denom


def proteoGenomicLookup(peptides, suffixarrayObject, samps2pops):

  peptides = list(sorted(peptides))
  noI= [ isoleucine2leucine(s) for s in peptides ]
  whos = suffixarrayObject.listAllPeople()
  
  search = suffer.ProteosPeptideSearch(suffixarrayObject, noI, suffixarrayObject.tmpdir)
  curMat = search.getNextMatch()

  print("Peptide", "Peptide20AA", "EnsemblID", "Chromosome", "Start", "Stop", "Strand", sep="\t")
  # walk through the peptide detections (once)
  # report the pooled allele frequency
  i=0
  while curMat is not None:
    curPep = curMat.sequence
    pepSought = peptides[i] # may be in 20AA alphabet
    i += 1

    prots = curMat.protData
    d = set()
    for prot in prots:

      if len(prot) < 7:
        continue
      
      protID = prot[1]
      chrom = prot[3]
      strand = prot[6]

      # lower and upper genomic coordinates
      lb = int(prot[4])
      ub = int(prot[5])
      if lb > ub:
        (lb, ub) = (ub, lb)
      
      if protID not in d:
        j = 7
        while j + 4 < len(prot):
          a = int(prot[j+1])
          b = int(prot[j+2])
          if a < lb:
            lb = a
          elif a > ub:
            ub = a
            
          if b < lb:
            lb = b
          elif b > ub:
            ub = b
          j += 4

        d.add(protID)
        print(curPep, pepSought, protID, chrom, lb, ub, strand, sep="\t")

    
    
    curMat = search.getNextMatch()
    
def summarizePanel(peptides, suffixarrayObject, samps2pops, verbose=False):
  """
  Computes some quick allele summary information on a panel...
  Note: Peptides have ISOLEUCINES
  """

  peptides = list(sorted(peptides))
  
  noI= [ isoleucine2leucine(s) for s in peptides ]

  search = suffer.ProteosPeptideSearch(suffixarrayObject, noI, suffixarrayObject.tmpdir)
  whos = suffixarrayObject.listAllPeople()

  popTots = {} # pop -> total sample size of population
  pepCounts = {} # peptide -> popuation -> count 
  if verbose and samps2pops is not None:
    # who is a haploid ID (ID + _MATERNAL_SUFFIX)
    for who in whos:
      dipID = who[:-2]
      if dipID in samps2pops:
        pop = samps2pops[dipID]
        if pop not in popTots:
          popTots[pop] = 1
        else:
          popTots[pop] += 1
  
  curMat = search.getNextMatch()
  print("Population", "Peptide", "Peptide20AA", "AlleleCount", "AlleleFrequency", "Total", sep="\t")

  popsOrdered = sorted(popTots.keys()) # fix the population ordering...
  
  # walk through the peptide detections (once)
  # report the pooled allele frequency
  i=0
  while curMat is not None:
    curPep = curMat.sequence
    pepSought = peptides[i] # may be in 20AA alphabet
    i += 1
    curOwn = curMat.getFullPeptideOwnership()
    print("Total", curPep, pepSought, len(curOwn), len(whos), (len(curOwn)+0.5)/(len(whos)+1), sep="\t")

    if verbose and samps2pops is not None:
      # initialize all population counts to 0
      pepCounts[pepSought]={}
      for pop in popsOrdered:
        pepCounts[pop]=0
      
      for own in curOwn:
        dipID = own[:-2]
        if dipID in samps2pops:
          pop = samps2pops[dipID]
          pepCounts[pop] += 1
          
      for pop in popsOrdered:
        print(pop, curPep, pepSought, pepCounts[pop], popTots[pop], (pepCounts[pop]+0.5)/(popTots[pop]+1) ,sep="\t")
        
    curMat = search.getNextMatch()
    
  
def computeTheta(peptides, suffixarrayObject, pops2samps, ALT_MIN_THETA=1e-9, infiniteDB=True):
  """
  This is a custom theta computation
  performed at the level of the locus.
  The locus is defined at the level of the chromosome;
  ie, considering all SAPs (translated into SNPs; or vice versa)
  but using the EXOME data.
  haplotypes that span the chromosome are at these sites are constructed
  and "theta" (Bi from Buckleton 2016; https://doi.org/10.1016/j.fsigen.2016.03.004)
  is computed and returned.
  This function is computed anew for each set of peptides for every chromosome detected
  """

  allInds = set([s[:-2] for s in suffixarrayObject.listAllPeople() ])
  nPopGot = 0
  for pop, samples in pops2samps.items():
    for s in samples:
      if s in allInds:
        nPopGot += 1
        break

  if nPopGot < 2:
    print("To estimate theta I need individuals in the suffix array from at least two populations...", \
          "for now 1e-6 is being used...", file=sys.stderr)
    return 1e-6
  

  popCounts = {} # population -> alleles -> counts (dictionary of dictionaries)
  popTots = {} # population -> diploidSamplesize (dictioary of string/int pairs)
  popHomozygosities = {} # this is M_i in Buckelton 2016
  # from this paper l is omitted as this computation is
  # done per locus on the fly.
  for pop, samples in pops2samps.items():

    # records, for a POPULATION, what the allele counts are
    hapCounts = defaultdict(int)
    (hapCounts, nSamp) = peptides2hot(peptides, suffixarrayObject, \
                                        len(samples)*2,\
                                        include=diploidIDsToHaploid(samples))
    
    popCounts[pop] = hapCounts

    homozygosities = 0.0
    if infiniteDB:
      divisor = (nSamp)*( nSamp) # approximation in paper
    else:
      divisor = (nSamp)*( nSamp-1) # exact computation (also in paper)
    #however,
    # this can create problems as the number of comparisons can be very
    # different between populations as opposed to within.
    # they also give the approximation (equivalent to this)
    # which also makes negative Fst a lot more difficult!
    
    popTots[pop] = nSamp
    # count is the number of times the haplotype was seen
    # the homozygosities is, given all pairs of individuals, the number of times
    # you see the same genotype. this is: count(count-1)/n(n-1)
    # changed to count^2/n^2
    # (equivalent to using allele frequencies as opposed to allele counts...)
    for count in hapCounts.values():
      if infiniteDB:
        homozygosities += (count*(count)) /divisor
      else:
        homozygosities += (count*(count-1)) /divisor
        
    popHomozygosities[pop] = homozygosities # the haplotype homozygosity for population pop

  pops = list( pops2samps.keys() )
  betweenPopHomo=0. # This is Mb from Buckleton 2016
  # the between-population homozygosities
  
  for (pop1, pop2) in combinations(pops, 2):
    pop1Alleles = popCounts[pop1]
    pop2Alleles = popCounts[pop2]
    denom = float(popTots[pop1]*popTots[pop2])
    # this computes Mij (Buckleton)
    for a1 in pop1Alleles:
      # if not in pop2Alleles then c2 is 0, and numerator is 0 => no addition needed!
      if a1 in pop2Alleles:
        c1 = pop1Alleles[a1]
        c2 = pop2Alleles[a1]
        numerator = c1*c2
        
        betweenPopHomo += numerator / denom

  # in the paper they divide by |r|*(|r|-1)
  # which is the number of pairs of populations evaluated
  # BUT
  # they evaluate all pairs redundantly (as ordered pairs, CEU and CHB, then CHB and CEU)
  # our computation uses unordered pairs (CEU,CHB but not CHB,CEU)
  # we need to adjust by /2
  mbDivisor = (len(pops)*(len(pops)-1))/2
  Mb = betweenPopHomo / mbDivisor

  Bw = 0 # Bw from Buckleton
  for pop in pops:
    Mi = popHomozygosities[pop]
    Bw += (Mi-Mb) /(1-Mb)
    
  # Taken from Appendix A of Buckleton  
  Bw /= len(pops)
    
  return max(Bw, ALT_MIN_THETA)

def haplotypeEventCalculator(observed, comparisonTo, comparisonToAlso, D, C, weights=None, logtrans=False ):
  """
  count up the number of dropin, contamination,
  and NOT dropin and NOT contamination "events" to explain the observed
  alleles assuming that the set of haplotypes in comparisonTo is the origin
  of said contamination.
  The number of each type of events is computed
  and the probability of the sequence of events (assuming independence)
  is returned
  """
  cont=drop=notdrop=notcont=0

  nAlleles = len(observed)

  contaminationWeight=1

  for i in range(nAlleles):
    total=0
    # in the pseudocode the two sets of haplotypes (known and unknown) are joined
    # it's more convenient to keep them separate (computationally)
    # one or the other vectors can be empty.
    for hap in comparisonTo:
      if hap[i] == "1":
        total += 1

    for hap in comparisonToAlso:
      if hap[i] == "1":
        total += 1
    

    if observed[i] == "0": # allele not observed in evidence
      if total > 0:
        drop += total # thus everyone that has it must've dropped it.
        notcont += 1  # and no, this isn't contamination
    else:
      if total==0: # observed, but nobody has it!
        cont += 1  # must be contamination
        
        if weights is not None: # weights[i] == the probabiltity of allele i dropping in. Typically this is set to the allele frequency.
          contaminationWeight *= weights[i] # when weights is none then the allele frequency is modeled as 1. This is the approach of Mitchell in their
          # admittedly flawed FST program, though I would argue the most major flaw therein was how the system was implemented.
        
      else:
        notdrop += total
        notcont += 1

  # equivalent to log of numpy.prod( ... )
  if logtrans:
    return( sum( [numpy.log(pow(C, cont)), numpy.log(pow(1-C, notcont)), numpy.log(pow(D, drop)), numpy.log(pow(1-D, notdrop)), numpy.log(contaminationWeight) ] ) )
    
  return numpy.prod( [pow(C, cont), pow(1-C, notcont), pow(D, drop), pow(1-D, notdrop), contaminationWeight])

def makeHotVector(peptidesInPanel, peptidesDetected):
  """
  Converts a set (datatype) of peptides detected
  and the panel of peptides (list)
  and returns a 1-hot encoded vector (as a string)
  The vector is an indicator variable; parallel to the order of the peptides in the panel
  1 means that the peptide is in the sample
  """
  det = ["0"] * len(peptidesInPanel)
  i=0
  for pep in peptidesInPanel:
    if pep in peptidesDetected:
      det[i] = "1"
    i += 1

  return "".join(det)
  

def findOwners(suffixArrayObject, peptides):
  """
  Newest iteration of suffer retooled how peptide search works
  and it changed the interface
  This is a retrofit (old interface output from new interface)
  returns a list of lists, parallel to peptides
  each inner list is a list of ints; ints are indexes in whos (list of haploid identifiers)
  s.t. that individual has that peptide
  """

  if not type(peptides) is list:
    print("Peptides type error!", type(peptides), file=sys.stderr)
    exit(1)

  if suffixArrayObject.detectsList is not None:
    search = suffer.ProteosPeptideSearch(suffixArrayObject, suffixArrayObject.detectsList, suffixArrayObject.tmpdir)
  else:
    search = suffer.ProteosPeptideSearch(suffixArrayObject, peptides, suffixArrayObject.tmpdir)
    
  whos = suffixArrayObject.listAllPeople()

  i=0
  w2i = {}
  for who in whos:
    w2i[who]=i
    i += 1

  out = []
  for p in peptides:
    out.append([])

  curMat = search.getNextMatch()
  i=0
  while curMat is not None:
    
    curPep = curMat.sequence
  
    # optimization: query is run once 
    if curPep not in peptides:
      curMat = search.getNextMatch()
      continue
    
    curOwn = curMat.getFullPeptideOwnership()
    inner = out[i]
    i += 1
    for who in curOwn:
      if who in w2i:
        inner.append( w2i[who] )

    curMat = search.getNextMatch()
    
  return (out, whos)

def peptides2hot(peptides, suffixarrayObject, nSamps, include=(), exclude=()):
  """
  This takes in a set of peptides (vector of strings)
  and a protein suffix array object (tells me which individuals have which peptides)
  and returns a dictionary associating a unique haplotype (binary, 1-hot encoded) with its count
  
  include and exclude lists will (only) include people in the inclusion list
  or exclude anybody in the exclusion list...
  inclusion/exclusion lists must refer to HAPLOID ids

  nSamps is the haploid sample size (which I need to get which people have the 0-vector haplotype)
  """

  # first, figure out who has which peptide
  #note: if an individual has none of the peptides, they are NOT in ww
  ww = dict() # who has what. associates an individual with the peptides they have

  global SUFFIX_QUERY_RESULTS
  global SUFFIX_QUERY
  global SUFFIX_QUERY_ARRAY
  global SUFFIX_CACHE
  # optimization: I cache queries. (all of them)
  # this lets me keep the previous query, and recycle the search results.
  #  if SUFFIX_QUERY is None or peptides != SUFFIX_QUERY or SUFFIX_QUERY_ARRAY != suffixarrayObject:
  # at most, the suffix cache be of size 2 (only 2 suffix arrays are supported, and all queries should involve the same peptides)
  i=0
  for tup in SUFFIX_CACHE:
    if tup[0] == suffixarrayObject and tup[1] == peptides:
      SUFFIX_QUERY_RESULTS=tup[2]
      SUFFIX_QUERY=peptides
      SUFFIX_QUERY_ARRAY=suffixarrayObject
      break
    i+=1

  # cache miss. 
  if i == len(SUFFIX_CACHE):
    SUFFIX_QUERY_RESULTS = findOwners(suffixarrayObject,peptides) # TODO: revamp. redirect to give an equivalent encoding.
    SUFFIX_QUERY = peptides
    SUFFIX_QUERY_ARRAY = suffixarrayObject
    SUFFIX_CACHE.append( (suffixarrayObject, peptides, SUFFIX_QUERY_RESULTS) )

  
  # idx is a list of lists (integers).
  # whos is a list of strings (sample IDs)
  # the inner list's integers are indexes into whos
#  res = suffixarrayObject.findOwners(peptides)
  res = SUFFIX_QUERY_RESULTS
  idx = res[0]
  whos = res[1]
  
  nPeps = len(peptides)

  global USE_WHITELIST
  global WHITELIST
  global SAPOINTER
  useWhitelist = USE_WHITELIST and SAPOINTER==suffixarrayObject
  
  if len(include):
    nSamps = len(include)

  if useWhitelist:
    nSamps = len(WHITELIST)

  for i in range(nPeps):
    for index in idx[i]:
      who = whos[index]

      if useWhitelist:
        if who not in WHITELIST:
          continue 
      
      if len(include) == 0 or who in include:
        if who not in ww:
          ww[who] = set()
          
        ww[who].add(peptides[i]) # haploid identifier (who)-> set of peptides associated with that (haploid) person
    
  
  if len(exclude):
    for e in exclude:
      if e in ww:
        del ww[e]
      # regardless if they have any of the (variant) peptides, the sample size is 1-less  
      nSamps -= 1
  
  if nSamps < 2:
    print("Not good. After: ", exclude , " your sample size is not possible (<=1) ", nSamps, file=sys.stderr)
    sys.exit(1)

  hotDict = {}
  for (who, alleles) in ww.items():
    listy = ["0"] * len(peptides)
    i = 0
    for pep in peptides:
      if pep in alleles:
        listy[i] = "1"
      i += 1

    # associate the haplotype (as a string) with its count
    s = "".join(listy)
    if s not in hotDict:
      hotDict[s] = 1
    else:
      hotDict[s] += 1

  # how many individuals have a vector of all 0s?
  # ie, none of the requested peptides?
  leftovers = nSamps - len(ww)
  if leftovers < 0:
    print("Cannot be. A negative number of haplotypes inferred!!", nSamps, len(ww), file=sys.stderr)
    sys.exit(1)

  if leftovers > 0:
    s = "0" * len(peptides)
    hotDict[s] = leftovers

  return(hotDict, nSamps)

def readSingleColumnFile(f, column2grab=1, sep="\t", header=False, callback=None):
  """
  Takes in a filename (string)
  and a column to grab (1-based index)
  and a separator
  the first line can be skipped (header=TRUE)
  and a callback function can be applied to each element
  and returns the column vector referred to (str type)
  None is the type contained in the vector element if the row index not referred to is not present.
  """

  if f is None or not os.path.exists(f):
    print("Failed to find file: ", f , "", sep="\n", file=sys.stderr)
    return [None]
  
  if f.endswith(".gz"):
    h = gzip.open(f)
  else:
    h = open(f)

  column2grab -= 1 # convert to 0-based indexing
  
  col = []
  for line in h:
    if header: # optionally skip the 1st line...
      header=False
      continue
    
    row = line.rstrip().split(sep)
    if column2grab < len(row):
      if callback is None:
        col.append(row[column2grab])
      else:
        col.append( callback(row[column2grab]) )
    else:
      col.append(None)

  h.close()
  return(col)

def readTwoColumnFile(f, column2grab=1, nextColumn2grab=2, sep="\t", header=False, callback=None, otherCallback=None):
  """
  Takes in a filename (string)
  and a column to grab (1-based index)
  and ANOTHER column to grab
  and a separator
  the first line can be skipped (header=TRUE)
  and a callback function can be applied to each element
  and returns the column vector referred to (str type)
  None is the type contained in the vector element if the row index not referred to is not present.

  Note that calling readSingleColumnFile twice would accomplish the same thing
  EXCEPT
  such an approach does not work with a stream.... this does.
  """

  if f is None or not os.path.exists(f):
    print("Failed to find file: ", f , "", sep="\n", file=sys.stderr)
    return [None]
  
  if f.endswith(".gz"):
    h = gzip.open(f)
  else:
    h = open(f)

  column2grab -= 1 # convert to 0-based indexing
  nextColumn2grab -= 1
  
  col = []
  otherCol = []
  
  for line in h:
    if header: # optionally skip the 1st line...
      header=False
      continue
    
    row = line.rstrip().split(sep)
    
    if column2grab < len(row):
      if callback is None:
        col.append(row[column2grab])
      else:
        col.append( callback(row[column2grab]) )
    else:
      col.append(None)

    if nextColumn2grab < len(row):
      if otherCallback is None:
        otherCol.append(row[nextColumn2grab])
      else:
        otherCol.append( otherCallback(row[nextColumn2grab]) )
    else:
      otherCol.append(None)

  h.close()
  return((col, otherCol))


def initSuffixArrayOld(binary, array, tmpdir, cutfile):
  """
  Initializes the suffix array object...
  (this is for the V2 implementation of profinman. now retired.
  """
  names = glob.glob(array + "/*names")
  if len(names)==0:
    print("No names file found for: ", array , "Not good!", file=sys.stderr, sep="\n")
    return None
  
  return suffer.SuffixArraySearch(binary, array, names[0], tmpdir, cutfile)

def initSuffixArray(binary, array, tmpdir):
  """
  Initializes the suffix array object...
  """
  
  sa = suffer.ProteosSuffixWrapper(binary, array)
  sa.tmpdir=tmpdir # kludge. Keep the directory as a cache
  sa.detectsList=None # When we cache hits (using suffer.py) the search needs to remain the same
  # the initial search is often genomic; this lets me dynamically subset
  return sa

def computeLRWithKnowns(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, state, saKnown, knownIDs, nUnknown):
  """
  This computes a nested likelihood ratio
  some number of known individuals (haploid ids: knownIDS) + 
  some number of unknown individuals (nUnknown; diploid counter) is assumed
  and the likelihood of the alleles | this hypothesis is assessed. (this is H1; aka Hp)

  then the known individuals are set to the empty set, and the number of
  unknowns is increased (s.t. the number of individuals, known+unknown remains the same)
  and the likelihood of this hypothsis (H2; aka Hd) is assessed.

  The arguments to this function are the same as crossValidate
  SAVE
  the list of known individuals (knownIDs, haploid IDs; consecutive pairs -> 1 diploid individual
  (and saKnown, which is the corresponding suffix array)
  and the number of unknown individuals (0+), counting diploid individuals
  """


  knownHaps = list()
  
  for i in range(0, len(knownIDs), 2):
    (haplotypes1, nSampK) = peptides2hot(peptides, saKnown, 2, include= knownIDs[i:(i+2)])
    
    theseHaps = list( haplotypes1.keys() )
    if len(theseHaps) == 1: # homozygote.
      knownHaps.append(theseHaps[0])
      knownHaps.append(theseHaps[0])
    else:
      knownHaps.append(theseHaps[0])
      knownHaps.append(theseHaps[1])

  dipIDs = sorted(list(set([who[:-2] for who in knownIDs ]))) 
  numer=computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=nUnknown, knownHaps = knownHaps, excludeList= [])
  denom=computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=nUnknown + len(dipIDs), knownHaps = [], excludeList= [])

  print("ID", "Dropin_Rate", "Dropout_Rate", "Denom", "Num", "Theta", state.keyHeader(), sep="\t")
    
  for tup in denom.keys():
    print(",".join(dipIDs), tup[0], tup[1], denom[tup], numer[tup], theta, str(state), sep="\t")
  

def allNestedLRs(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, state, saKnown, numKnown):
  """
  Computes all nested likelihood ratios where the total sample size is numKnown (numerator and denominator)
  It does this for all individuals in some suffix array (saKnown)

  For a description of the arguments see allSingleSourceLR

  """
  
  allHaps = getAllDiplotypes(saKnown, peptides)
  diploidIds = sorted(list(allHaps.keys()))

  print("ID", "NKnown", "NTotal", "Dropin_Rate", "Dropout_Rate", "Likelihood", "Theta", state.keyHeader(), sep="\t")

  # save on redundant computations; if we specify the same set of knowns and same number of unknowns, keep the likelihood
  denomCache= {}
  
  for i in range(1, (numKnown+1)):

    totalSS = i
    # consider every person (i==1),
    # every pair of people (i==2),
    # ...
    for peeps in combinations(diploidIds, i):
      kHaps = []
      peeps = list(peeps)
      numID = ",".join(peeps)

      for person in peeps:
        kHaps.extend(allHaps[person])

      # and compute the likelihood of the combination
      num = computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=0, knownHaps = kHaps, excludeList=[])
      for tup in num.keys():
        print(numID , i, i, tup[0], tup[1], num[tup], theta, str(state), sep="\t")
          
      # and consider all proper subsets (ie, the powerset save the identitical set)
      for subset in chain.from_iterable(combinations(peeps,n) for n in range(totalSS)):
        nKnown=len(subset)
        subsetID = ",".join(subset)
        if nKnown==0:
          subsetID = "Random"
        kHaps = []
        for person in subset:
          kHaps.extend(allHaps[person])
          
        # key is: the known individuals , the # of unknowns
        khapStr = ";".join(kHaps)
        cacheKey = (khapStr, totalSS-nKnown)
        if cacheKey not in denomCache:
          denom = computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=totalSS-nKnown, knownHaps = kHaps, excludeList=[])
          denomCache[ cacheKey ] = denom
        else: # redundant computation. skip it.
          continue

          #print("DenominatorID", "NKnownDenom", "Dropin_Rate", "Dropout_Rate", "Denom", "Num", "Theta", state.keyHeader(), sep="\t")
        for tup in denom.keys():
          print(subsetID, len(subset), i, tup[0], tup[1], denom[tup], theta, str(state), sep="\t")
        
      
def allSingleSourceLR(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, state, saKnown):
  """
  This computes a single-source LR
  specifically, it takes every person in a suffix array (saKnown)
  computes the p(alleles from some panel (peptides) |person)
  and then reports the p(alleles | random person ). The latter is wrt to a population (saObject)
  """
  
  denom=computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=1, knownHaps = [], excludeList= [])
  for tup in denom.keys():
    print("RandomPerson", tup[0], tup[1], denom[tup], theta, str(state), sep="\t")
    
  peeps = saKnown.listAllPeople() # these are haploid, but I can't guarantee the order
  diploidIds = sorted(list(set( [f[:-2] for f in peeps] ))) # now diploid 
  hapIds = diploidIDsToHaploid(diploidIds) # and ordered, haploid

  print("ID", "Dropin_Rate", "Dropout_Rate", "Likelihood", "Theta", state.keyHeader(), sep="\t")
  
  for i in range(0, len(hapIds), 2):
    (haplotypes1, nSampK) = peptides2hot(peptides, saKnown, 2, include= hapIds[i:(i+2)])

    knownHaps=[]
    theseHaps = list( haplotypes1.keys() )
    if len(theseHaps) == 1: # homozygote.
      knownHaps.append(theseHaps[0])
      knownHaps.append(theseHaps[0])
    else:
      knownHaps.append(theseHaps[0])
      knownHaps.append(theseHaps[1])

    numer=computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=0, knownHaps = knownHaps, excludeList= [])
    for tup in denom.keys():
      print(hapIds[i][:-2], tup[0], tup[1], numer[tup], theta, str(state), sep="\t")

def getAllDiplotypes(saObject, peptides):
  """
  Convenience function
  Takes in a suffix array, and a peptide panel and creates a dictionary associating an individual (diploid)
  with a pair of haplotypes 
  This calls a suffix array lookup N times, for N individuals
  (that part could be smarter)
  But it is nice when we want to look at pairs, triples, ... of known individuals...

  returns a dictionary associating the ID of the person to the pair of haplotypes..
  """
  d = {}
  
  peeps = saObject.listAllPeople()
  # remove the _1/_2 identifiers; get me the diploid sample IDs
  diploidIds = sorted(list(set( [f[:-2] for f in peeps] )))
  for person in diploidIds:
    ashap = diploidIDsToHaploid(person) # and remake the haploid ids so we know to include them.
    (haplotypes1, nSampK) = peptides2hot(peptides, saObject, len(peeps), include= ashap)

    kHaps = list( haplotypes1.keys() )
    if len(kHaps) == 1: # homozygote.
      kHaps = [kHaps[0], kHaps[0]]
    d[person] = kHaps

  return d

def crossValidate(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, state, nCross=1, computeDenominator=True):
  """
  Runs leave-one-out cross-validation on a suffix array (internal)
  Take a list of *alleles* detected from some set of peptides (a panel)
  assume some dropoutRate and some dropinRate

  then, take each diploid person from the saObject
  put them in the numerator of the LR
  *and*
  remove them from the suffix array

  then compute the denominator on the modified suffix array
  (note removal is implicit)
  
  other variables (theta, nMC) are used in the LR computation...
  and state is a program state variable (reproducibility)

  """
  
  allHaps = getAllDiplotypes(saObject, peptides)

  # remove the _1/_2 identifiers; get me the diploid sample IDs
  diploidIds = allHaps.keys()

  print("ID", "Dropin_Rate", "Dropout_Rate", "Denom", "Num", "Theta", state.keyHeader(), sep="\t")
  for peeps in combinations(diploidIds, nCross):
    kHaps = []
    ashap = []

    for person in peeps:
      ashap.extend( diploidIDsToHaploid(person) ) # and remake the haploid ids so we know to exclude them.
      kHaps.extend(allHaps[person])



    if computeDenominator:
      e = computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=nCross, knownHaps = [], excludeList= ashap)
      s = computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=0, knownHaps = kHaps, excludeList=[])

      for tup in e.keys():
        print(",".join(peeps), tup[0], tup[1], e[tup], s[tup], theta, str(state), sep="\t")

    else:
      s = computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC, nUnknown=0, knownHaps = kHaps, excludeList=[])
      for tup in s.keys():
        print(",".join(peeps), tup[0], tup[1], -1, s[tup], theta, str(state), sep="\t")



def computeRMP(saObject, peptides, theta, numIts, dropin, state, altMinCount=5.0):
  """
  This is an adaptation of the RMP of Woerner et al. 2020 https://doi.org/10.1016/j.fsigen.2020.102295
  It takes in a suffix array object, and list of peptides detected
  A theta correction value 
  And it further needs a vector/list of dropin rates and the program state,
  as well as an alternative minimum allele count (wherein all counts < this are set to this)
  RMPs are used to express the rarity of some set of peptides detected under the assumption
  that the data are single source, with the variety provided permitting some global dropin (false +)
  rate.
  """
  
  peeps = saObject.listAllPeople()
  (haplotypes, nSamp) = peptides2hot(peptides, saObject, len(peeps) )
  
  nPeps = len(peptides)

  hItems = [ (h[0], h[1]) for h in haplotypes.items() ]

  altMin = altMinCount/nSamp
  altRmp = heterozygousProbWithTheta(theta, altMin, 1-altMin)


  # for a given mask (1 means peptide dropped out, 0 means not)
  # the RMP is deterministic
  # when drop-in is small, the same mask will be seen often. This recycles...
  cache = {}
  # stores marginal theta-corrected RMPS and Naive rmps
  posteriorT = [-1] * numIts
  posteriorN = [-1] * numIts


  print("Dropin", "Npeps", "RMP_theta", "RMP_theta_0.025", "RMP_theta_0.975", "RMP_naive", "RMP_naive_0.025", "RMP_naive_0.975", state.keyHeader(), sep="\t")
  
  for d in dropin:

    # Monte carlo iterations
    for pI in range(numIts):
      
      mask = ['0'] * nPeps
      for i in range(nPeps):
        if random.random() < d:
          mask[i] = '1'
        
      mask = "".join(mask)

      # cache miss
      if mask not in cache:
        # all of the heterozygous cases...
        sumRMP = 0.
        naiveRMP = 0.
        for ((hap1, count1), (hap2, count2)) in combinations(hItems, 2):
      
          broken=False
          for i in range(nPeps):
            if mask[i] == '1': # masked base; simulates that this base is dropin
              continue
            if hap1[i] == '0' and hap2[i]=='0': # niether haplotype has the peptide. Cannot have left evidence.
              broken=True
              break
        
          if not broken:
            sumRMP += heterozygousProbWithTheta(theta, max(altMinCount, count1)/nSamp, max(altMinCount, count2)/nSamp)
            # 2pq
            naiveRMP += (count1/nSamp) * (count2/nSamp) * 2
            
        # homozygous case
        for (hap, count) in hItems:
          broken=False
          for i in range(nPeps):
            if mask[i] == '1':
              continue
            if hap[i] == '0':
              broken=True
              break

          if not broken:
            sumRMP += homozygousProbWithTheta(theta, max(altMinCount, count)/nSamp)
            # p^2
            naiveRMP += (count/nSamp) * (count/nSamp)
            
      
        if sumRMP < altRmp:
          sumRMP = altRmp
        elif sumRMP > 1: # many 5/2n can add up to more than one... as can the homo/hetero prob computation as it approaches a frequency of 1
          sumRMP = 1.
        
        cache[mask] = (sumRMP, naiveRMP)

      (rmpT, rmpN) = cache[mask]
      posteriorT[pI] = rmpT
      posteriorN[pI] = rmpN
      
    eRmpT = sum(posteriorT)/numIts
    eRmpN = sum(posteriorN)/numIts
    ordsT = numpy.quantile(posteriorT, [0.025, 0.975])
    ordsN = numpy.quantile(posteriorN, [0.025, 0.975])

    print(d, nPeps, eRmpT, ordsT[0], ordsT[1], eRmpN, ordsN[0], ordsN[1], state, sep="\t")
    
      
def computeEverything(saObject, alleles, peptides, theta, dropoutRate, dropinRate, nMC=-1, nUnknown=1, knownHaps = [], excludeList=[]):
  """
  Needs: suffix array object
  alleles detected (set)
  peptides in panel (vector)
  estimate of theta (from computeTheta, or from the command line)
  a dropin (error rate, false positive rate; as a vector)
  a dropout (error rate, false negative rate; as a vector)
  the number of Monte Carlo samples to take (-1 -> exhaustive search)
  the number of DIPLOID individuals to sample from the database
  a list of known haplotypes 
  """
  if type(peptides) != list:
    print("Your arguments are messed. Check your peptides data-type", file=sys.stderr)
    return None
  if len(alleles) > len(peptides):
    print("Your arguments are messed. Your panel is smaller than the number of alleles you've detected?!", file=sys.stderr)
    return None
  if nUnknown == 0 and len(knownHaps)==0:
    print("Cannot happen. No unknown or knowns?!", file=sys.stderr)
    return None
  
  observed = makeHotVector(peptides, alleles)

  global LOGLIKES
  logTransform=LOGLIKES

  alleleWeights=None
  global WEIGHT_ALLELES
  global WEIGHTS
  if WEIGHT_ALLELES:
    alleleWeights=WEIGHTS

  finalOutput = {}
  # handle the case of NO unknown haplotypes separately
  if nUnknown == 0:
    for c,d in product(dropinRate, dropoutRate):
      # can decouple events from probability (more efficient)
      eventProb = haplotypeEventCalculator(observed,\
                                           [], knownHaps,\
                                           d, c, alleleWeights, logTransform)

      finalOutput[(c,d)] = eventProb

      
    return finalOutput

  # from this point there is some number of unknown diploid individuals
  # first: get the haplotypes
  peeps = saObject.listAllPeople()
  (haplotypes, nSamp) = peptides2hot(peptides, saObject, len(peeps), exclude=excludeList )

  # haplotypes are now ordered
  hapVectors = list(haplotypes.keys())
  hapCounts = [ haplotypes[i] for i in hapVectors ]
  # make the ecd
  hapCumulative = list(accumulate(hapCounts))

  nHaps = len(hapVectors)
    
  intList = list(range(nHaps))


  for unkInts in unknownHaplotypeGenerator(nHaps, nUnknown, nMC, hapCumulative):
          
    knownProb = 1.0

    haps = [ hapVectors[j] for j in unkInts ]

    
    if nMC > 0:
      knownProb = 1.0/nMC
    else:
      # the probability of the unknown haplotypes
      for j in range(0, len(haps), 2):
        ## TODO: handle case for MC sample.## (prob is wrong)
        if unkInts[j] == unkInts[j+1]: # min b/c the adjusted probability computation is statistically biased ( may exceed 1)
          knownProb *= min(1.0, homozygousProbWithTheta(theta, hapCounts[ unkInts[j] ] / hapCumulative[-1]))
        else:
          knownProb *= min(1.0, heterozygousProbWithTheta(theta, hapCounts[ unkInts[j] ] / hapCumulative[-1],hapCounts[ unkInts[j+1] ] / hapCumulative[-1]) )

    # all combinations (cartesian product) of drop-in and dropout rates
    for c,d in product(dropinRate, dropoutRate):

      # can decouple events from probability (more efficient)
      eventProb = haplotypeEventCalculator(observed,\
                                           haps, knownHaps,\
                                           d, c, alleleWeights, logTransform)

      # use the log-sum-exp trick to prevent underflow...
      if logTransform:
        if (c,d) not in finalOutput:
          finalOutput[(c,d)] = []
        
        finalOutput[(c,d)].append( numpy.log( knownProb ) + eventProb )

        
      else:
        if (c,d) not in finalOutput:
          finalOutput[(c,d)] = knownProb * eventProb
        else:
          finalOutput[(c,d)] += knownProb * eventProb

  # we don't keep the sum in this case, but an array of log likelihoods
  # we want the sum of the likelihoods
  # uses the log-sum-exp trick to get the log of this sum
  # convert to a likelihood (optionally) with numpy.exp( finalOutput[k] )
  if logTransform:
    for k in finalOutput:
      finalOutput[k] = numpy.logaddexp.reduce( finalOutput[k ] )

  return finalOutput


def preprocessArgv(argv):
  """
  This lets command lines be specified with -@ / --aptetail FILE
  The file can have arbitrary numbers of command line arguments specified in it
  from this a new argv is constructed (missing arg[0] and missing the -@ FILE)
  """
  newargv=[]

  i=0
  stop = len(argv)-1
  while i < stop:
    i += 1
    
    if argv[i] == '-@' or argv[i] == '--apetail':
      i += 1
      foo = readSingleColumnFile(argv[i], 1, sep="\n", header=False)
      for line in foo:
        for word in line.split():
          newargv.append(word)
    else:
      newargv.append(argv[i])

  return newargv



def peptide_main(argv):
  parser = argparse.ArgumentParser(description="Let's compute some likelihood ratios!\n")


  parser.add_argument('-p', '--peptide_hits',           dest='P', help="single-column file; no headers; peptides detected", type=str, default='')
  parser.add_argument('-a', '--peptide_panel',          dest='A', help="single-column file; no headers; peptides in panel", type=str, default='')
  parser.add_argument('-c', '--contamination_rate',     dest='C', help="Contamination rate(s) (i.e., drop-in rates) used in the LR calculation.",type=float, default=[1e-6, 0.03, 0.05], nargs="+")
  parser.add_argument('-d', '--dropout_rate',           dest='D', help="Dropout rate(s) (i.e., false negative rates) used in the LR calculation.",type=float, default=[1e-6, 0.05, 0.10], nargs="+")
  parser.add_argument('-S', '--suffix_array',           dest='Array', help="the suffix array (directory) itself; used to generate *unknown* individuals; for single source, trypically used in denominator", type=str, default='/home/becrloon/ProtengineR2/Prot_1000_AC')
  parser.add_argument('-L', '--known_suffix_array',     dest='KnownArray', help="the suffix array (directory) itself; used to generate *known* individuals; for single-source, typically used in the numerator (and the denominator in more complex cases)", type=str, default='/home/becrloon/ProtengineR2/Prot_25_AC')

  
  parser.add_argument('-V', '--cross_validate',         dest='V', help="Runs cross-validation; -V 1 runs leave-one out cross-validation, -V 2 is leave two-out, ...; combine with -S", type=int, default=0)
  parser.add_argument('-X', '--cross_validate_nodenom', dest='X', help="Runs cross-validation without the denominator; -X 1 runs leave-one out cross-validation, -X 2 is leave two-out, ...; combine with -S", type=int, default=0)
  parser.add_argument('-1', '--all_single_source',      dest='SingleSource', help='compute all single source LRs; computes the likelihood for each person in one suffix array (-L, numerator) relative to the haplotypes generated from another (-S, denominator)', action="store_true")
  parser.add_argument('-W', '--weights',                dest='W', help='Weights contamination events by the allele frequency', action="store_true")
  parser.add_argument('-l', '--loglikes',               dest='LL', help='Returns log likelihoods instead of raw likelihoods. Use when/if underflow occurs.', action="store_true")
  parser.add_argument('-N', '--all_nested_lrs',         dest='Nested', help='compute all nested LRs; -N 1 equivalent to -1 ; -N 2 computes all nested LRs that involve a 2-person mixture, -N 3 for 3-person...', default=0, type=int)
  parser.add_argument('-r', '--rmp',                    dest='R', help="Computes the RMP using -R Monte Carlo simulations as per Woerner et al. 2019. Respects -P", default=-1, type=int)
  
  
  
  # manually compute the likelihood function; specifying particular people... (fast, but you need to manually specify your hypotheses)
  parser.add_argument('-K', '--known_individuals',      dest='KNOWN', help="Named known individuals under H1 (aka Hp).",type=str, default=[], nargs="+")
  parser.add_argument('-U', '--n_unknown',              dest='U', help="the number of (additional) unknown individuals", type=int, default=0)

  parser.add_argument('-P', '--population',              dest='Pops', help="Restrict to the population (refers to -q)", type=str, default="")
  
  
  ###
  ### Only for use when run with Monte Carlo sampling... (not recommended; no theta correction.)
  ###
  parser.add_argument('-s', '--seed',                   dest='S', help="This sets the seed on the random number generator; helps with reproducibility", type=int, default=1)
  parser.add_argument('-i', '--iterations',             dest='I', help="The number of Monte Carlo *I*terations", type=int, default=-1)

  ###
  ### Used for the theta;
  ###
  # manually set theta:
  parser.add_argument('-n', '--no_theta',               dest='N', help="Turns off the theta computation (theta=0)", action="store_true")
  parser.add_argument('-t', '--theta',                  dest='T', help="Sets theta to a fixed constant", type=float, default=-1)
  # or use a point-estimate of theta:
  parser.add_argument('-q', '--population_to_samples',  dest='Q', help="This is a file gives the population IDs of (some of the) samples specified in the suffix array (-S)", type=str, default="")


  ## uncommon parameters
  ## for the suffix array
  parser.add_argument('-T', '--tempDir',                dest='Tmp', help="temporary directory for suffix array lookup (pooled)", type=str, default='TMP')
  parser.add_argument('-Y', '--tempNum',                dest='Tmp_N', help="temporary directory for suffix array lookup (numerator, results cached)", type=str, default='TMP_N')
  parser.add_argument('-Z', '--tempDenom',                dest='Tmp_D', help="temporary directory for suffix array lookup (denominator, results cached)", type=str, default='TMP_D')
  parser.add_argument('-B', '--suffix_array_binary',    dest='Binary', help="the binary (executable) suffix array (include the path)", type=str, default="/home/becrloon/ProtengineR3/zproj_profinman/builds/bin_x64_linux/profinman") #default='/home/becrloon/ProtengineR2/Proteos/proteos/workflows/suffix_array/builds/bin_x64_linux')
#  parser.add_argument('-C', '--cut_file',               dest='Cut', help="suffix array: how cuts (e.g., trypsin digest) are defined", type=str, default="/home/becrloon/ProtengineR2/Digests/trypsin_only.dig")

  parser.add_argument('-@', '--apetail',                dest='APE', help="Additional arguments can be specified in the file", type=str, default='')
  parser.add_argument('-Q', '--quick_summaries',        dest='Summaries', help='Computes summary statistics on a peptide panel', action="store_true", default=False)
  parser.add_argument('-F', '--full_summaries',         dest='Full', help="Same as -Q, but it also prints who has which alleles (may be big!!)", action='store_true', default=False)
  parser.add_argument('-G', '--proteogenomic_summaries', dest='ProtGenom', help="Provides genomic information of peptides", action='store_true', default=False)


  argv = preprocessArgv(argv)

  results = parser.parse_known_args(argv)[0]
  args = parser.parse_known_args(argv)[1]

  global IS_RANDOMIZED
  IS_RANDOMIZED=True
  random.seed(results.S)
  if results.I < 1:
    IS_RANDOMIZED = False

  global SUFFIX_QUERY_RESULTS
  global SUFFIX_QUERY
  global SUFFIX_QUERY_ARRAY
  global SUFFIX_CACHE
  SUFFIX_QUERY_ARRAY = SUFFIX_QUERY_RESULTS = SUFFIX_QUERY = None
  SUFFIX_CACHE = []

  global WEIGHT_ALLELES
  WEIGHT_ALLELES=False
  global WEIGHTS
  WEIGHTS=None
  
  global LOGLIKES
  LOGLIKES=False
  if results.LL:
    LOGLIKES=True

  global USE_WHITELIST
  global WHITELIST
  global SAPOINTER
  USE_WHITELIST=False
  WHITELIST=SAPOINTER=None
  
    
  # programatically, treat cross validation (no denom) as cross-validation with a negative number (-1 == loocv with w. no denom)
  if results.X:
    results.V = -results.X
  

  tmpdir = results.Tmp
  if results.Tmp_D != "TMP_D": # default
    tmpdir = results.Tmp_D
    
      # init the suffix array
  saObject = initSuffixArray(results.Binary, results.Array, tmpdir)
  if saObject is None:
    parser.print_help()
    sys.exit(1)

  errors=0
  if results.Summaries or results.Full or results.ProtGenom:
    peptidePanel = readSingleColumnFile(results.A, callback=None)
  elif results.W:
    (peptidePanel, weights) = readTwoColumnFile(results.A, callback=isoleucine2leucine)
    if None in weights:
      print("Problem parsing the weights from the peptide panel file. I need a two-column file of peptides; the first is the peptide, the second is the weight.", file=sys.stderr)
      errors+=1
    elif len(weights) == 0:
      print("Your peptide panel is ... blank?", file=sys.stderr)
      sys.exit(1)
    try:
      weights = [ float(w) for w in weights ]
    except ValueError:
      print("At least one of the float values is not a number", file=sys.stderr)
      errors += 1

    # set up a global variable to house the weights
    # everything is set up as a parallel array
    # so just referring to the index is sufficient.
    WEIGHT_ALLELES=True
    WEIGHTS=weights

    if len(weights) and (min(weights) < 0 or max(weights)>1):
      print("Proper weights must be frequencies between 0 and 1...", file=sys.stderr)
      errors += 1
      
  elif results.R < 1: # no panel for RMPs...
    peptidePanel = readSingleColumnFile(results.A, callback=isoleucine2leucine)
  else:
    peptidePanel = []
    
  if None in peptidePanel and not( results.Summaries or results.Full or results.ProtGenom ):
    print("Problems parsing the peptide panel file: ", results.A, " at least one of the rows is... blank? irregular?", file=sys.stderr)
    parser.print_help()
    sys.exit(1)
  elif errors:
    parser.print_help()
    sys.exit(1)
    
  samps2pops=None
  pops2samps=None
    
  if results.Summaries or results.Full or results.ProtGenom:
    if samps2pops is None and not results.Summaries:
      (samps2pops, pops2samps) = parseSamps2Pops(results.Q, "Individual ID", "Population")
    if results.ProtGenom:
      proteoGenomicLookup(peptidePanel, saObject, samps2pops)
    else:
      summarizePanel(peptidePanel, saObject, samps2pops, results.Full)
    exit(0)

  if results.Q != "" and not os.path.isfile(results.Q):
    print("There is no file: " , results.Q , file=sys.stderr, sep="\n")
    parser.print_help()
    sys.exit(1)  
    

  peptidesDetected = set( readSingleColumnFile(results.P, callback=isoleucine2leucine) )
  saObject.detectsList= list( sorted( peptidesDetected))
  
  if None in peptidesDetected:
    print("Problems parsing the peptide-hits file: ", results.P, " at least one of the rows is... blank? irregular?", file=sys.stderr)
    parser.print_help()
    sys.exit(1)
    
  # forcibly take an intersection of the panel and the detects (let's the panel vary and the detects stay samsies)
  if results.R < 1: # but only if it's an LR. If it's an RMP there's no panel 
    peptidesDetected = set([ p for p in peptidesDetected if p in peptidePanel] )
  else: #case of RMP. we need the detects to be ordered and indexable

    if len(peptidePanel):
      peptidesDetected = list([ p for p in peptidesDetected if p in peptidePanel])
    else:
      peptidesDetected = list( sorted(peptidesDetected ))

    

  theta = results.T
  if results.N or results.Q == "":
    theta = 0
  elif theta < 0:
    (samps2pops, pops2samps) = parseSamps2Pops(results.Q, "Individual ID", "Population")
    theta = computeTheta(peptidesDetected, saObject, pops2samps, 0)    
    
  if len(peptidePanel) == 0 and results.R<1:
    print("Your panel is blank. Nothing to compute", file=sys.stderr)
    sys.exit(0)

  
  if results.U < 0:
    print("The number of unknowns must be a non-negative integer! Not , ", results.U, file=sys.stderr)
    parser.print_help()
    sys.exit(1)

  if results.Pops != "":
    if pops2samps is None:
      (samps2pops, pops2samps) = parseSamps2Pops(results.Q, "Individual ID", "Population")
      if pops2samps is None:
        print("Cannot filter to population: " , results.Pops, " no populations have been defined", file=sys.stderr)
        exit(1)
      elif results.Pops not in pops2samps:
        print("Cannot filter population: ", results.Pops, "It's not in the file: ", results.Q, pops2samps.keys(), sep="\n", file=sys.stderr)
        exit(1)

    # alll people in some population
    wl= diploidIDsToHaploid( pops2samps[ results.Pops ] ) 
    WHITELIST = set()
    USE_WHITELIST=True
    SAPOINTER=saObject
    allpeeps = set(saObject.listAllPeople()) #all people in all populations
    for w in wl:# intersect (don't ask me why you have to do this, but you do. All files are coming from 1KG project...)
      if w in allpeeps:
        WHITELIST.add(w)
        
    if len(WHITELIST)==0:
      print("Should never happen. The intersection of the database and the population: ", results.Pops, " is empty...", file=sys.stderr)
      exit(1)

  if results.R > 0:

    IS_RANDOMIZED = True
    
    state =State("RMP" + str(results.R), results.Pops, results.S)
    
    e = computeRMP(saObject, peptidesDetected, theta, results.R, results.D, state)
    
    
  elif results.V != 0:

    state =State("CV" + str(results.V), results.Pops, results.S)
    computeDenom=True
    if results.V < 0:
      computeDenom=False
      results.V = -results.V

    e = crossValidate(saObject, peptidesDetected, peptidePanel, theta, results.D, results.C, results.I, state, results.V, computeDenom)
  elif results.KNOWN:
    #results.KNOWN  are diploid ids (eg, SA001)
    # they need to be mapped into haploid ids (eg, SA001_1, SA001_2) (below)
    known = []
    
    saObjectKnown = initSuffixArray(results.Binary, results.KnownArray, results.Tmp_N)
    saObjectKnown.detectsList = saObject.detectsList
    
    allKnowns = saObjectKnown.listAllPeople()
    nFound = 0
    for k in results.KNOWN:
      if k + MATERNAL_SUFFIX in allKnowns:
        nFound += 1
        known.append( k + MATERNAL_SUFFIX )
        known.append( k + PATERNAL_SUFFIX )
      else:
        print("Failed to find known individual: ", k, file=sys.stderr)
    
    if nFound == 0 or nFound < len(results.KNOWN):
      t = sorted(list(set([who[:-2] for who in allKnowns ]))) # strip out haploid identifiers..., uniquify, and sort...
      print("Problems finding the known individuals. Their IDs must match (exactly).",
            " ".join(t),
            "is a list of them all...", "", sep="\n", file=sys.stderr)
      sys.exit(1)

    state =State("LR_Known_" + ";".join(results.KNOWN) + ":" + str(results.U) + "_unknown", results.Pops, results.S)
    e = computeLRWithKnowns(saObject, peptidesDetected, peptidePanel, theta, results.D, results.C, results.I, state, saObjectKnown, known, results.U)
  elif results.SingleSource:
    saObjectKnown = initSuffixArray(results.Binary, results.KnownArray, results.Tmp_N)
    saObjectKnown.detectsList = saObject.detectsList # needed for v2 of suffer.py. (caching)
    
    state =State("LR_Known_AllSingleSource", results.Pops, results.S)
    e = allSingleSourceLR(saObject, peptidesDetected, peptidePanel, theta, results.D, results.C, results.I, state, saObjectKnown)  
  elif results.Nested:
    saObjectKnown = initSuffixArray(results.Binary, results.KnownArray, results.Tmp_N)
    saObjectKnown.detectsList = saObject.detectsList
    
    state =State("LR_Nested_" + str(results.Nested), results.Pops, results.S)
    e = allNestedLRs(saObject, peptidesDetected, peptidePanel, theta, results.D, results.C, results.I, state, saObjectKnown, results.Nested)
  else:
    print("\n\nI don't know what to do!\n\n", file=sys.stderr)
    parser.print_help()
    sys.exit(1)

    
if __name__ == "__main__":
  peptide_main(sys.argv)

    
