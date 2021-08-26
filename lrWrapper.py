#!/usr/bin/python3

import csv
import sys
import glob
import os
import argparse
import gzip
import pparall

MIN_THETA=1e-6
DEFAULT_PROTENGINE_DIR = 'ProtengineR3' # must be the RELATIVE path. by default, this needs to be in the src/ directory for the main script
DEFAULT_MAX_PROCS=10 # regardless of the number of CPUs on the system, this is the max asked for by default

MIN_CONTAMINATION=0.05

def runOrDie(command):
  """
  Homage to perl's die command
  this exec's the command in the shell
  AND it captures the return/exit code
  and tests for 0 (clean exit)
  """  
  if os.system(command) != 0:
    print("Some problem with command:" , command , file=sys.stderr, sep="\n")
    exit(1)
    
  return True

def getArrays(callingFrom):
  """
  This takes in sys.argv[0]
  and it finds the directory in which PeptideLR was originally installed
  and it grabs all directories in DEFAULT_PROTENGINE_DIR (as a relative path in that dir)
  this corresponds to all arrays ALSO
  it corresponds to the binary (profinman/protengine)
  """
  path = os.path.realpath(callingFrom)
  dirname = os.path.join( os.path.dirname(path) , DEFAULT_PROTENGINE_DIR )
  out = {}
  for dirp in os.listdir(dirname):
    fullpath = os.path.join(dirname, dirp)
    if os.path.isdir( fullpath ):
      out[dirp] = fullpath
      
  return out
  

def getPepLR(callingFrom):
  path = os.path.dirname( os.path.realpath(callingFrom) )
  lr = os.path.join(path,  "peptideLR.py")
  combiner = os.path.join(path,  "combineLikelihoods.R -c " + str(MIN_CONTAMINATION))
  return (lr, combiner)

def getBinary(zprojDir):
  """
  Takes in a path (zproj directory)
  and returns a full path to the binary for the suffix array
  """
  bins = glob.glob(zprojDir + "/**/profinman", recursive=True)
  if len(bins):
    return bins[0]
  return None


def makeCommands(profinman, pepLR, combiner, args, detects, outdir, nullArray, altArray=None, pepFreqs=None):


  catCommand = 'cat' # cat (or catenate) ; unix utility
  if os.name != 'posix':
    catCommand = 'type' # cat equivalent on Windows...

  pepLRCommand = pepLR + " -B " + profinman + " -T " + os.path.join(outdir, args.N) + " -S " + nullArray 

  if os.path.isfile( os.path.join(nullArray, "samples2populations.tsv")  ):
    pepLRCommand += " -q " + os.path.join(nullArray, "samples2populations.tsv")

    
  # peptide query
  outfile = os.path.join(outdir, "peptideQuery." + args.N + ".tsv")
  detectsFile = os.path.join(outdir, "rawDetects.tsv")
  queryResultsFile =  os.path.join(outdir, "queryResults." + args.N + ".tsv")

  if not os.path.isfile(detectsFile):
    with open(detectsFile, "w") as fh:
      for peps in detects.values():
        for pep in peps:
          print(pep, file=fh)
          
  if not os.path.isfile(outfile):
    with open(outfile, "w") as fh:
      command = pepLRCommand
      if False and os.path.isfile( os.path.join(nullArray, "samples2populations.tsv")  ):
        command += " -F"
      else:
        command += " -Q"
      
      command += " -a " + detectsFile + " > " + queryResultsFile
      print(command, file=fh)

    # this performs the suffix array search (once)
    # wrt to the null suffix array (1000 Genomes and the like)
    if not os.path.isfile(queryResultsFile):
      runOrDie(command)

    if not os.path.isfile(queryResultsFile):
      print("Cannot find file: " , queryResultsFile, "Must be some file permission issues...", file=sys.stderr)
      exit(1)

  genomicFile = os.path.join(outdir, "genomicInformation.tsv")
  if args.G and not os.path.isfile(genomicFile):
    command = pepLRCommand + " -G -a " + detectsFile + " > " + genomicFile
    with open(outfile, "a") as fh:
      print(command, fh)
      
    runOrDie(command)


  # we query the suffix array once (regardless of the type)
  # this serves to do the peptide lookup once
  # and the downstream rmp/lr calculations can recycle the original query
  if altArray is not None:
    outfile = os.path.join(outdir, "peptideQuery." + args.A + ".tsv")
    queryResultsFileAlt =  os.path.join(outdir, "queryResults." + args.A + ".tsv")
    if not os.path.isfile(queryResultsFileAlt):
      # Note that the TMP dir changes 
      pepLRCommand = pepLR + " -F -B " + profinman + " -T " + os.path.join(outdir, args.A) + " -S " + altArray 

      with open(outfile, "w") as fh:
        command = pepLRCommand
        command += " -a " + detectsFile + " > " + queryResultsFileAlt
        print(command, file=fh)
      
      runOrDie(command)
      
  # invert the detects (was chrom -> peptides. make peptides -> chrom)
  pep2chrom = dict()
  # panel file handles
  handles = {}

  makePanels=False
  panelFiles = {}

  if not os.path.isdir( os.path.join(outdir, "Panels")):
    os.mkdir(os.path.join(outdir, "Panels"))
  
  for chrom in detects:
    for pep in detects[chrom]:
      pep2chrom[pep] = chrom

    for population in [args.P]:
      filename = os.path.join(outdir, "Panels", "panel." + population + "." + chrom)
      if population == 'Total':
        panelFiles[(population, chrom)] = "-W -a " + filename
      else:
        panelFiles[(population, chrom)] = "-W -a " + filename + " -P " + population + " -q " + os.path.join(nullArray, "samples2populations.tsv")
      
      if not population in handles:
        handles[population]={}
      if not os.path.isfile(filename):
        handles[population][chrom] = open(filename, "w")
        makePanels=True

        
  # create seperate PANEL files
  # one for each chromosome / population pair.
  if makePanels:
    if pepFreqs is None:

      with open(queryResultsFile) as fh:
        first=True
        for line in fh:
          if first:
            first=False
            continue
          sp = line.rstrip().split("\t")
          pop = sp[0]
          pep = sp[2] # 20AA peptide sequence
          freq = sp[-1]

          if pep not in pep2chrom:
            print("Should never happen! No partition for pep:" , pep, file=sys.stderr)
            exit(1)
        
          chrom = pep2chrom[pep]
          # population we care about
          if pop in handles and chrom in handles[pop]:
            fh = handles[pop][chrom]
            print(pep, freq, sep="\t", file=fh)
          
    else: # user-supplied allele weights (frequencies)
      for (pep, freq) in pepFreqs.items():
        chrom = pep2chrom[pep]

        for (pop, inner) in handles.items():
          print(pep, freq, sep="\t", file=inner[chrom])
        
  # close the barn door!
  for pop in handles:
    for chrom in handles[pop]:
      handles[pop][chrom].close()

  # -Y is for the numerator/population suffix array
  pepLRCommand = pepLR + " -B " + profinman + " -Z " + os.path.join(outdir, args.N) + " -S " + nullArray + " -t " + str(args.T)

  commands = []

  logFH = open(os.path.join(outdir,"logFile.txt"), "a+")

  rmpFiles = []
  # make per-chromosome annotations
  if args.R:
    
    if not os.path.isdir( os.path.join(outdir, "RMPs")):
      os.mkdir(os.path.join(outdir, "RMPs"))
      
    command = pepLRCommand + " -p " + detectsFile
    
    outfileBase = os.path.join(outdir, "RMPs", "RMP")
    command += " -r 1000"
    
    for block in detects:

      for pop in [args.P]:
        outfile = outfileBase + "." + pop + "." + block
        append = ""
        rmpFiles.append(outfile)
        
        if not os.path.isfile(outfile):
          #print(command , append , panelFiles[(pop, block)] ,  ">", outfile, sep=' ')
          finalCommand = " ".join( (command,  append , panelFiles[(pop, block)] ,  ">", outfile))
          print(finalCommand, file=logFH)
          if args.C < 2: # just one core. avoid fancy parallelism and just run it synchronously
            os.system(finalCommand)
          else:
            commands.append(finalCommand)

  # block-level (ie, chromosome level) RMPs/LRs get written to directories of that name

  # LR computation.
  # add in MC sampling (if that's requested)
  if args.M>0:
    pepLRCommand += " -i " + str(args.M)


  likeFiles = []
  if args.L:

    if not os.path.isdir( os.path.join(outdir, "LRs")):
      os.mkdir(os.path.join(outdir, "LRs"))
    
    command = pepLRCommand + " -p " + detectsFile + " -Y " + os.path.join(outdir, args.A) + " -L " + altArray 
    outfileBase = os.path.join(outdir, "LRs", "LR."+ str(args.L))

    if args.L==1: # -1 implementation is faster. Let's use that if we can.
      command += " -1"
    else:
      command += " -N " + str(args.L)
    
    for block in detects:
      for pop in [args.P]:
        outfile = outfileBase + "." + pop + "." + block
        append = ""
        likeFiles.append(outfile)
        
        if not os.path.isfile(outfile):
          #print(command , append , panelFiles[(pop, block)] ,  ">", outfile, sep=' ')
          finalCommand = " ".join( (command , append , panelFiles[(pop, block)] ,  ">", outfile))
          print(finalCommand, file=logFH)
          if args.C < 2:
            os.system(finalCommand)
          else:
            commands.append(finalCommand)

  if args.W:
    if args.W > 1:
      print("Only single-contributor estimates of the LR are possible for the W-estimator: ", args.W, file=sys.stderr)
      
    if not os.path.isdir( os.path.join(outdir, "LRs")):
      os.mkdir(os.path.join(outdir, "LRs"))

    command = pepLRCommand + " -p " + detectsFile + " -Y " + os.path.join(outdir, args.A) + " -W -a - -L " + altArray + " -g " + str(args.W)
    outfileBase = os.path.join(outdir, "LRs", "LR."+ str(args.W))
    for pop in [args.P]:
      outfile = outfileBase + "." + pop + ".genomic"
      append = ""

      catty = catCommand + " " + os.path.join(outdir, "Panels", "panel." + pop + ".* | ")
      
      if not os.path.isfile(outfile):

        finalCommand = " ".join( (catty, command ,  ">", outfile))
        print(finalCommand, file=logFH)
        if args.C < 2:
          os.system(finalCommand)
        else:
          commands.append(finalCommand)

    
  if len(commands)>0:
    if pparall.inParallel(commands, args.C, check=True):
      print("At least one command failed!", file=sys.stderr)
      exit(1)

  if len(rmpFiles):
    outfile = os.path.join(outdir, "RMPs", "CombinedRMP.tsv")
    if not os.path.isfile(outfile):
      rmpCommand = "Rscript " + combiner + " " + " ".join(rmpFiles) + " > " + outfile
      print(rmpCommand, file=logFH)
      #print(rmpCommand)
      os.system(rmpCommand)

  if len(likeFiles):
    outfile = os.path.join(outdir, "LRs", "CombinedLR." + str(args.L) + ".tsv")
    if not os.path.isfile(outfile):
      lrCommand = "Rscript " + combiner + " " + " ".join(likeFiles) + " > " + outfile
      print(lrCommand, file=logFH)

      os.system(lrCommand)
    
  logFH.close()  
      
def buildArgvParser(parser):
  """
  Name says it all. Processes command line arguments"
  """
  parser.add_argument('-r', '--rmp', dest='R', help="Computes the RMP", action='store_true')
  parser.add_argument('-l', '--likelihoods', dest='L', help="Computes likelihoods for 1..L contributors", type=int, default=0)
  parser.add_argument('-w', '--w_likelihoods', dest='W', help="Computes the likelihood of Woerner et al.", type=int, default=0)
  parser.add_argument('-p', '--population', dest='P', help="Sets the reference population (P) in the likelihood estimation. See -W. Defaults to pooled frequencies (Total)", type=str, default="Total")
  parser.add_argument('-t', '--theta', dest='T', help="Turns on the theta-correction", default=MIN_THETA)
  parser.add_argument('-d', '--detects', dest='D', help="A file with the peptide detections...", default="-")
  parser.add_argument('-P', '--detects_peptide_colname', dest="pepcol", help="In -D, the column name for the peptide detections", default='Peptide')
  parser.add_argument('-C', '--detects_chromosome_colname', dest="chromcol", help="In -C, the column name for the chromosome (or any categorical variable used to partition the detections)", default='Chromosome')
  parser.add_argument('-F', '--detects_frequency_colname', dest="freqcol", help="Rather than using population-specific frequencies use those defined in the column specified in the detects file ", default='')
  parser.add_argument('-M', '--monte_carlo_sims', dest="M", help="Number of Monte Carlo simulations", default=0, type=int)
  

  parser.add_argument('-q', '--query_allele_frequencies', dest='Q', help="Computes allele frequencies on --detects", action='store_true')
  parser.add_argument('-g', '--genomic', dest='G', help="Generates genomic information on peptides", action='store_true')
  parser.add_argument('-n', '--null_array', dest='N', help="The null/reference array. Default: HG38_Clean", default="HG38_Clean")
  parser.add_argument('-a', '--alt_array', dest='A', help="The comparison array", default="")
  parser.add_argument('-o', '--output_directory', dest='O', help="The directory where the analysis is conducted", default='')
  parser.add_argument('-c', '--n_cpus', dest='C', help='The number of CPUs (degree of multi-processing) used. Only applies to RMP/LR calculation', default=min(DEFAULT_MAX_PROCS, os.cpu_count()))
  parser.add_argument('-L', '--ls', dest='ls', help="Lists the suffix arrays in the default directory", action='store_true')
  parser.add_argument('-W', '--which_pops', dest='WhichPops', help="Lists the populations available in the suffix array", action='store_true')
  
def parser_main(argv):

  parser = argparse.ArgumentParser(description="Let's generate some command files!")
  buildArgvParser(parser)
  results = parser.parse_known_args(argv[1:])[0]
  args = parser.parse_known_args(argv[1:])[1]
  if len(args):
    print("Should never be! Unused arguments: ", args,file=sys.stderr)
    return 1


  suffixArraysAndMore = getArrays(argv[0])
  binary = ""
  for k in suffixArraysAndMore.keys():
    if "profinman" in k:
      if binary =="":
        binary = k
      else:
        print("Installation problem. Multiple binaries detected...", binary, k, file=sys.stderr)
        return 1
      
  if binary == "":
    print("Installation problem! No profinman detected!", suffixArraysAndMore.keys(), file=sys.stderr)
    return 1

  
  if results.ls:
    print("The current arrays are: ")
    print("\n".join( [ k for k in suffixArraysAndMore.keys() if k != binary ] ) )
    return 1

  profinman = getBinary( suffixArraysAndMore[binary])
  (lr, combiner) = getPepLR(argv[0])

  if profinman is None:
    print("Problem finding profinman...", file=sys.stderr)
    return 1

  if lr is None:
    print("Problem finding the LR calculator...", file=sys.stderr)
    return 1
  
  if results.D != "-" and not os.path.isfile(results.D):
    print("File: " , results.D , " cannot be found. Maybe there's a typo in the name?", file=sys.stderr)
    exit(1)
    
      
  if results.D == '-':
    fh = sys.stdin
  elif results.D.endswith(".gz"):
    fh= gzip.open(results.D, "rt")
    outdir = os.path.splitext(os.path.basename(results.D))[0]
    if outdir.find(".") != -1:
      outdir = os.path.splitext(outdir)[0]
    else:
      outdir += "_out"
      
  else:
    fh = open(results.D, "r")
    outdir = os.path.splitext(os.path.basename(results.D))[0]
    if outdir.find(".") == -1: # hoo-ray for the case of no file extension..
      outdir += "_out"
    

  if results.O != '':
    outdir = results.O

    
  if results.N not in suffixArraysAndMore:
    print("Cannot find the null/reference suffix array: ", results.N, " Run\n", argv[0] , " --ls \nto list all suffix arrays", file=sys.stderr)
    exit(1)

  altArray = None
  if results.A != "":
    if results.A not in suffixArraysAndMore:
      print("Cannot find the alt suffix array: ", results.A , " Run\n", argv[0] , " --ls \nto list all suffix arrays", file=sys.stderr)
      exit(1)
    altArray = suffixArraysAndMore[ results.A ]


  if results.L > 0 and altArray is None:
    print("Cannot be: You asked to compute a LR but you did not specify an alternative (alt) suffix array!", file=sys.stderr)
    exit(1)

  if results.W > 0 and altArray is None:
    print("Cannot be: You asked to compute a LR (W estimate) but you did not specify an alternative (alt) suffix array!", file=sys.stderr)
    exit(1)
    
    
  refArray = suffixArraysAndMore[ results.N ]
  otherPops = False
  if results.P != 'Total':
    otherPops=True
    
  if results.WhichPops or otherPops:
    popsFile = os.path.join(refArray, "samples2populations.tsv")
    # print the unique populations...
    d = set()
    with open(popsFile) as popsfh:

      first=True
      if not otherPops:
        print("Available populations:")
        print("Total")
        
      for line in popsfh:
        if first:
          first=False
          continue
        s = line.rstrip().split("\t")
        if not otherPops and s[-1] not in d:
          print(s[-1])
        d.add(s[-1])
    if not otherPops:
      exit(0)
    else:

      if results.P != "Total" and results.P not in d:
        print("Cannot find population: ", pop, file=sys.stderr)
        print("One or more errors detected... the available populations are" , "\n".join(d), file=sys.stderr)
        exit(1)
          
  outdir = os.path.abspath(outdir)
  
  if not os.path.exists(outdir):
    os.mkdir(outdir)

  detects = {}
  reader=csv.DictReader(fh, delimiter="\t")
  if results.chromcol not in reader.fieldnames and results.chromcol != "":
    print("Your detections file has no column: ", results.chromcol, file=sys.stderr)
    return 1
    
  if results.pepcol not in reader.fieldnames:
    print("Your detections file has no column: ", results.pepcol, file=sys.stderr)
    exit(1)

  pepFreqs = None
  if results.freqcol != "":
    if results.freqcol not in reader.fieldnames:
      print("Your detections file has no column: ", results.freqcol, file=sys.stderr)
      exit(1)
    pepFreqs = {} # peptide -> peptide frequency
    
  for row in reader:
    if results.chromcol != "":
      chrom = row[ results.chromcol ]
    else:
      chrom='pooled'
      
    pep   = row[ results.pepcol ]

    if chrom not in detects:
      detects[chrom] = set()
    detects[chrom].add(pep)

    if results.freqcol != "":
      freq = float( row[ results.freqcol ] )
      pepFreqs[pep]=freq
    
  if fh != sys.stdin:
    fh.close()

    
  if len(detects)>0:
    makeCommands(profinman, lr, combiner, results, detects, outdir, refArray, altArray, pepFreqs)

if __name__ == "__main__":
  parser_main(sys.argv)






