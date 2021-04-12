#!/usr/bin/python3

import csv
import sys
import glob
import os
import argparse
import gzip

MIN_THETA=1e-6
DEFAULT_PROTENGINE_DIR = 'ProtengineR3' # must be the RELATIVE path. by default, this needs to be in the src/ directory for the main script

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
  return lr

def getBinary(zprojDir):
  """
  Takes in a path (zproj directory)
  and returns a full path to the binary for the suffix array
  """
  bins = glob.glob(zprojDir + "/**/profinman", recursive=True)
  if len(bins):
    return bins[0]
  return None


def makeCommands(profinman, pepLR, args, detects, outdir, nullArray, altArray=None):


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
      if os.path.isfile( os.path.join(nullArray, "samples2populations.tsv")  ):
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
      pepLRCommand = pepLR + " -Q -B " + profinman + " -T " + os.path.join(outdir, args.A) + " -S " + altArray 

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

    for population in args.P:
      filename = os.path.join(outdir, "Panels", "panel." + population + "." + chrom)
      if population == 'Total':
        panelFiles[(population, chrom)] = "-W -a " + filename
      else:
        panelFiles[(population, chrom)] = "-W -a " + filename + " -P " + population
      
      if not population in handles:
        handles[population]={}
      if not os.path.isfile(filename):
        handles[population][chrom] = open(filename, "w")
        makePanels=True

        
  # create seperate PANEL files
  # one for each chromosome / population pair.
  if makePanels:
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

  # close the barn door!
  for pop in handles:
    for chrom in handles[pop]:
      handles[pop][chrom].close()

  # -Y is for the numberator/population suffix array
  pepLRCommand = pepLR + " -B " + profinman + " -Z " + os.path.join(outdir, args.N) + " -S " + nullArray 
  # make per-chromosome annotations
  if args.R:
    
    if not os.path.isdir( os.path.join(outdir, "RMPs")):
      os.mkdir(os.path.join(outdir, "RMPs"))
      
    command = pepLRCommand + " -p " + detectsFile
    
    outfileBase = os.path.join(outdir, "RMPs", "RMP")
    command += " -r 1000"
       
    for block in detects:

      for pop in args.P:
        outfile = outfileBase + "." + pop + "." + block
        append = ""

        if not os.path.isfile(outfile):
          print(command , append , panelFiles[(pop, block)] ,  ">", outfile, sep=' ')

  # block-level (ie, chromosome level) RMPs/LRs get written to directories of that name  
  if args.L:

    if not os.path.isdir( os.path.join(outdir, "LRs")):
      os.mkdir(os.path.join(outdir, "LRs"))
    
    command = pepLRCommand + " -p " + detectsFile + " -Y " + os.path.join(outdir, args.A) + " -L " + altArray
    outfileBase = os.path.join(outdir, "LRs", "LR")

    if args.L==1: # -1 implementation is faster. Let's use that if we can.
      command += " -1"
    else:
      command += " -N " + str(args.L)
      
    for block in detects:
      for pop in args.P:
        outfile = outfileBase + "." + pop + "." + block
        append = ""

        if not os.path.isfile(outfile):
          print(command , append , panelFiles[(pop, block)] ,  ">", outfile, sep=' ')

      
def buildArgvParser(parser):
  """
  Name says it all. Processes command line arguments"
  """
  parser.add_argument('-r', '--rmp', dest='R', help="Computes the RMP", action='store_true')
  parser.add_argument('-l', '--likelihoods', dest='L', help="Computes likelihoods for L contributors", type=int, default=0)
  parser.add_argument('-p', '--population', dest='P', help="Sets the reference population (P) in the likelihood estimation. Defaults to pooled frequencies (Total)", type=str, nargs="+", default=["Total"])
  parser.add_argument('-t', '--theta', dest='T', help="Turns on the theta-correction", default=MIN_THETA)
  parser.add_argument('-d', '--detects', dest='D', help="A file with the peptide detections...", default="-")
  parser.add_argument('-P', '--detects_peptide_colname', dest="pepcol", help="In -D, the column name for the peptide detections", default='Peptide')
  parser.add_argument('-C', '--detects_chromosome_colname', dest="chromcol", help="In -C, the column name for the chromosome (or any categorical variable used to partition the detections)", default='Chromosome')

  parser.add_argument('-q', '--query_allele_frequencies', dest='Q', help="Computes allele frequencies on --detects", action='store_true')
  parser.add_argument('-g', '--genomic', dest='G', help="Generates genomic information on peptides", action='store_true')
  parser.add_argument('-n', '--null_array', dest='N', help="The null/reference array. Default: HG38_Clean", default="HG38_Clean")
  parser.add_argument('-a', '--alt_array', dest='A', help="The comparison array", default="")
  parser.add_argument('-o', '--output_directory', dest='O', help="The directory where the analysis is conducted", default='')
  parser.add_argument('-c', '--n_cpus', dest='C', help='The number of CPUs (degree of multi-processing) used. Only applies to RMP/LR calculation', default=1)
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
  lr = getPepLR(argv[0])

  
  if profinman is None:
    print("Problem finding profinman...", file=sys.stderr)
    return 1

  if lr is None:
    print("Problem finding the LR calculator...", file=sys.stderr)
    return 1
  
  outdir = results.O
  if outdir == '':
    if results.D=='-':
      outdir='out'

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
    fh = open(results.D, "r")
    outdir = os.path.splitext(os.path.basename(results.D))[0]

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
    
    
  refArray = suffixArraysAndMore[ results.N ]
  if results.WhichPops:
    popsFile = os.path.join(refArray, "samples2populations.tsv")
    # print the unique populations...
    with open(popsFile) as fh:
      d = set()
      first=True
      print("Available populations:")
      print("Total")
      for line in fh:
        if first:
          first=False
          continue
        s = line.rstrip().split("\t")
        if s[-1] not in d:
          print(s[-1])
        d.add(s[-1])
    exit(0)
    
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
    
  for row in reader:
    if results.chromcol != "":
      chrom = row[ results.chromcol ]
    else:
      chrom='pooled'
      
    pep   = row[ results.pepcol ]

    if chrom not in detects:
      detects[chrom] = set()
    detects[chrom].add(pep)

  if fh != sys.stdin:
    fh.close()

  makeCommands(profinman, lr, results, detects, outdir, refArray, altArray)

if __name__ == "__main__":
  parser_main(sys.argv)






