>#

The package herein can be used to compute a semi-continuous likelihood ratio [(to know more about what that is and how we do it, click here!)](about.md) using peptides detected by [liquid chromatography tandem mass spectrometry](https://en.wikipedia.org/wiki/Liquid_chromatography%E2%80%93mass_spectrometry) (LC-MS/MS). \
It can also be used to compute some genomic properties of peptides. Namely, where in the genome they might be found. Of importance, this query is done on a population database; so long as the allele is present (in at least one individual in the database) it is query-able. 

# Dependencies (pretty standard)
* For Windows users
  * [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) (WSL)
  * Note: WSL is not a strict requirement for the code to run, but the installation instructions need to be made Windows-compatible 
* g++
  * pthreads 
  * zlib1g-dev (zlib, development version)
* Python 3.*
  * Numpy
  * joblib
* R
  * tidyverse

* Compute recommendations
  * ~65 Gb of storage (to store the ~4000 whole proteomes)
  * 8 cores (soft requirement; the more the better)
  * 16 Gb of memory (less important, but the more the better)

<br>

# Quick start for *nix including [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) (WSL)
User-level installation (installs to $HOME/bin/, which is assumed to be in your $PATH): \
Note this includes downloading ~4000 whole proteomes (twice), which takes a while!\
Maybe it's time to grab a coffee?
<br>

## Installation
```
git clone --recursive https://github.com/Ahhgust/PeptideLR.git # download this project and Dr. Crysup's suffix array project (protengine)
cd PeptideLR/protengine # and build
chmod +x Buildit_x64_linux.sh && ./Buildit_x64_linux.sh && ln -s $PWD/pywrap/suffer.py ..
# and grab the public datasets (forgive the dropbox link but these files are BIG!)
cd ..
mkdir ProtengineR3 && cd ProtengineR3 && wget -O ProtengineR3.zip 'https://www.dropbox.com/sh/xp3wzs5fy9taqvl/AADjmPYTT201_MPtNXjVEZuaa?dl=1'
unzip ProtengineR3.zip && rm ProtengineR3.zip
for file in *.tbz; do tar -xf $file && rm $file & done
wait
cd ..
ln -s $PWD/protengine ProtengineR3/zproj_profinman
if [ -d $HOME/bin ]; then ln -s $PWD/lrWrapper.py $HOME/bin; fi
```
<br>

Note, this program can be installed at the system level, but the appropriate (user/group/other-readable) location must be used to install the code, and the lrWrapper.py in the bin (e.g., in /usr/bin or /usr/local/bin) **needs to be a symlink**!

## Quick start
Using the Example/ directory as, well, an example: <br>
* First, test that the program was successfully installed. Type: `lrWrapper.py --ls`. This *should* list the currently available suffix arrays. You should see
  * GRCh38Ref
  * HG38_Clean
  * HG38
  * (any additional suffix arrays that you've built; if not, check your warning/error messages from the installation)
* Second, test the calculator itself. From the `PeptideLR/` directory, type
  * `cd Example`
  * `lrWrapper.py -n HG38_Clean -a GRCh38Ref -l 1 -r -g -q -p NFE -d peptides.tsv `
    * (the above should take a minute or so to run; please forgive the stray "warnings" from the tidyverse)

The above lrWrapper command will compute single-source likelihoods (-l 1; -l 2 would consider all single source and two-person mixtures), the RMP (-r), give genomic summary stats (-g) and query the allele frequencies (-q) for the peptides given by -d (peptides.tsv). The results will be written to: peptides/. The "null" hypothesis (-n) in this case is used to generate the haplotypes for the random man (HG38_Clean, which gives the clean/passed SNP calls from the [1000 Genomes + HGDP](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg) dataset), and the "alternative" hypothesis (-a) is that the contributor is the reference sequence (GRCh38Ref). This is obviously not an appropriate question to consider (for one, the reference is haploid and a composite of several individuals), but it should serve as a jumping off point.
<br>

## Walk-through
The code-base is designed to work on peptides, specifically on peptides that are polymorphic in populations. To make the code easier to use a wrapper script (`lrWrapper.py`) was made. The wrapper can be used to compute a [semi-continuous likelihood](about.md), [random match probability](https://doi.org/10.1016/j.fsigen.2020.102295), as well as basic assessments of peptide allele frequency and it's proteo-genomic location(s).

The code expects tabular data (tab-delimited). These data need to be written as a regular file (no streaming). The file itself needs to have two named columns (see below; additional columns are okay and the order of the columns is arbitrary). The first column (Peptide) needs to give the peptide, and the second column gives the Chromosome associated with a given peptide. Technically the "chromosome" is used to group peptides; that is, the likelihood is assessed on each set of peptides that share the same group. Chromosomes are (biologically) independent (in which case the label "chromosome" needs to stay), but in principle protein identifiers (or chromosome + megabase identifiers) can be used if there is a reasonable belief that such groupings make sense (knowing that the groups need to represent statistically independent pieces of biological information). Note that a peptide may not be associated with multiple chromosomes. Let's say we want to look up the information on four peptides:
<br>

| Peptide | Chromosome |
| :---------: | :--------: |
| AASSQTPTMCTTTVTIK  | 18 |
| AASSQTPTMCTTTVTVK  | 18 |
| ADFSGMSAEK         | 18 |
| ADFSGMSTEK         | 18 |

<br>
Let's see what sort of genomic annotations are associated with these peptides. To do so, let's make a directory to work in (let's call it foo/. Pick some other name if that name has already been chosen.)
<br>

`mkdir foo && cd foo`

<br>

Now, select the above table (including the headers!) and copy it (ctrl+c). Let's make a file with this contents: Type:
<br>

`tr -d ' ' | cat > peptides.tsv`

<br>

and paste the content of the copied table (e.g., the top-left icon in putty) and type ctrl+d (which marks the end of file).

Let's begin with pulling out some genomic annotations associated with these peptides. To do so type:
<br>

`lrWrapper.py -d peptides.tsv -n HG38_Clean -g`

<br>
This command takes a little while to complete. What it is doing is searching every protein sequence in every individual from the 1000 Genomes Project (HG38_Clean), and checking to see which individuals have a tryptic peptide that is equivalent (remembering that Is and Ls are mass-equivalent) to the ones queried. Of those individuals, it then pulls out the information on which proteins correspond to these peptides.

<br>

When the above completes it will create a directory ( `peptides/` ). In it will be a file ( `genomicInformation.tsv` ) that looks like:

| Peptide            | Peptide20AA          |   EnsemblID       | Chromosome  |   Start   |  Stop      |  Strand |
| :----------------: | :------------------: |   :-------------: | :---------: | :-------: | :--------: | :-----: |
| AASSQTPTMCTTTVTLK  |AASSQTPTMCTTTVTIK     |   ENSP00000257197 | chr18       | 31130802  |  31162594  |      -  |
| AASSQTPTMCTTTVTLK  |AASSQTPTMCTTTVTIK     |   ENSP00000257198 | chr18       | 31130513  |  31162594  |      -  |
| AASSQTPTMCTTTVTVK  |AASSQTPTMCTTTVTVK     |   ENSP00000257197 | chr18       | 31130802  |  31162594  |      -  |
| AASSQTPTMCTTTVTVK  |AASSQTPTMCTTTVTVK     |   ENSP00000257198 | chr18       | 31130513  |  31162594  |      -  |
| ADFSGMSAEK         |ADFSGMSAEK            |   ENSP00000331368 | chr18       | 63978308  |  63987278  |      +  |
| ADFSGMSAEK         |ADFSGMSAEK            |   ENSP00000438328 | chr18       | 63983700  |  63987278  |      +  |
| ADFSGMSAEK         |ADFSGMSAEK            |   ENSP00000381072 | chr18       | 63978308  |  63987278  |      +  |
| ADFSGMSTEK         |ADFSGMSTEK            |   ENSP00000331368 | chr18       | 63978308  |  63987278  |      +  |
| ADFSGMSTEK         |ADFSGMSTEK            |   ENSP00000438328 | chr18       | 63983700  |  63987278  |      +  |
| ADFSGMSTEK         |ADFSGMSTEK            |   ENSP00000381072 | chr18       | 63978308  |  63987278  |      +  |

<br>

This file gives the peptide sequence (in a 19 amino acid alphabet, with **I**sos -> **L**eus), the 20 amino-acid version (what was queried), and the protein ID (EnsemblID), as well as the genomic coordinates of the protein (GRCh38/hg38, chromosome, start, stop, strand) of the protein (not the peptide). What these annotations mean is that there is at least one individual in the [1000 Genomes + HGDP](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg) dataset that has these tryptic peptides, and the protein records associated with these peptides are as described above.

<br>

The data in this table are presented in a "[tidy](https://r4ds.had.co.nz/tidy-data.html)" format. If a peptide is associated with multiple proteins (as is expected), the same **Peptide** be repeated across multiple rows.
It is a reasonable bet that many of these proteins correspond to alternative transcripts of the same gene. They are different proteins, but proteins that likely correspond to the same gene, a fact that is alluded to by the shared start/stop coordinates.

<br>

This utility can also be used to estimate the peptide allele frequency. Even though the GVPs that we are searching for have a basis in the genome (they overlap some number of SNPs), estimating the GVP allele frequency can be tricky. For instance, if you have a biallelic non-synonymous SNP, any other protein-altering SNP that lands in the same region (regardless if you're looking for that SNP) will impact the frequency (generally lower) of one or both of the GVPs sought. Hard examples are when the GVP is multi-copy (within the same chromosome), in which case there is no way to estimate the allele frequency unless you know which SNPs cause the formation of the GVP *and* you know the co-variance (that is, which haploid individuals have both alleles). 
Instead, we estimate the allele frequency of peptides (GVPs) directly in the protein sequences to circumvent these issues. Namely, we compute the allele frequency as:

( | ![interval-union](https://latex.codecogs.com/svg.latex?\cup) (all haploid individuals that have some GVP)  | + 0.5) / ( 2*number of diploid individuals + 1)

where *have* means that the GVP is a substring of at least one of their proteins and it is cleavable by trypsin.

<br>

The population genetic data we provide is from gnomAD, and we use the population labels available [here](https://gnomad.broadinstitute.org/news/2018-10-gnomad-v2-1/). Note that we are restricted to those populations that are publically available (so not all populations are available)


<br>

# RMP

  `lrWrapper.py` uses caching to speed up the computation. In short, the suffix array is large, and while the search times are fast, it is faster still to keep a local copy of the query (and result). e.g., the directory `pep/HG38_Clean` has the results of search the peptides in `pep/rawDetects.tsv` in the HG38_Clean suffix array. To see the benefits of using a local cache, try running <br>

`lrWrapper.py -d peptides.tsv -n HG38_Clean -r`

Specifying -r computes the [random match probability](https://doi.org/10.1016/j.fsigen.2020.102295) <br>
Wherein it computes the RMP at the level of the chromosome, and combines these probabilities by taking products. RMPs (and likelihoods) <br>
are computed in parallel (control with -c) and then combined (RMPs/CombinedRMP.tsv).
It is important to remember what the RMP means in this context-- it is useful for saying how rare a collection of GVPs are in a single-source sample. It cannot be used, however, to associate a person to some piece of evidence. For that, a likelihood ratio is more appropriate.

# Likelihood

`lrWrapper.py` can be used to estimate a semi-continuous likeihood ratio [(what is this? click here!)](about.md). While the RMP requires a pair of information (the GVPs and some population databaes), the likelihood ratio requires three pieces of information: the GVPs, a population database, and a database of proposed contributors. Additional parameters include the hypothesis, as well as the proposed rates of drop-in and drop-out. The code-base provides an (over)simplified way of considering a hypothesis. Namely, two suffix arrays are given, and the number of individuals in the mixture is proposed (as an upper bound). Then all ways of constructing hypotheses involving that many individuals (or fewer) in constructed. In the simplest case (single source sample, specified with `-l 1`), each person in the contributor suffix array database (specified with `-a`) has the likelihood assessed, as well as the likelihood from the population database (roughly speaking, the mean is taken over all individuals in the database). If a two-person mixture is proposed, all single source likelihoods are assessed, as well as all ways of constructing a 2-person mixture (involving all pairs of known contributors (specified in the array `-a`)), as well as for each contributor plus a random individual, and as well for two random individuals.
<br>

The likelihoods get written to files in the `LRs/` directory (e.g., in the Quick Start, peptides/LRs/ is the directory referred to). These files include the likelihoods estimated within a chromosome (using haplotypes), as well as the combined genomic/proteomic LR (CombinedLR.?.tsv, where ? is the number of individuals proposed), wherein the final likelihood is combined by products over chromosomes (ie, by the product rule). 



## All options

  -h, --help            show this help message and exit <br>
  -r, --rmp             Computes the [RMP](https://doi.org/10.1016/j.fsigen.2020.102295) <br>
  -l L, --likelihoods L
                        Computes likelihoods for 1..L contributors   <br>
  -p P, --population P  Sets the reference population (P) in the likelihood
                        estimation. See -W. Defaults to pooled frequencies (Total) <br>
  -t T, --theta T       Turns on the theta-correction <br>
  -d D, --detects D     A file with the peptide detections... <br>
  -P PEPCOL, --detects_peptide_colname PEPCOL
                        In -D, the column name for the peptide detections <br>
  -C CHROMCOL, --detects_chromosome_colname CHROMCOL
                        In -C, the column name for the chromosome (or any
                        categorical variable used to partition the detections) <br>
  -F FREQCOL, --detects_frequency_colname FREQCOL
                        Rather than using population-specific frequencies use
                        those defined in the column specified in the detects
                        file <br>
  -M M, --monte_carlo_sims M
                        Number of Monte Carlo simulations <br>
  -q, --query_allele_frequencies
                        Computes allele frequencies on --detects <br>
  -g, --genomic         Generates genomic information on peptides <br>
  -n N, --null_array N  The null/reference array. Default: HG38_Clean <br>
  -a A, --alt_array A   The comparison array <br>
  -o O, --output_directory O
                        The directory where the analysis is conducted <br>
  -c C, --n_cpus C      The number of CPUs (degree of multi-processing) used.
                        Only applies to RMP/LR calculation <br>
  -L, --ls              Lists the suffix arrays in the default directory <br>
  -W, --which_pops      Lists the populations available in the suffix array <br>






## Installation notes and other gotchas
The code base is designed to work on a single file system on a single operating system (ie, at most, one binary, Windows or *nix is supported). This is typically the case, but some more exotic situations (two operating systems mounting the same file system) are not handled<br>
The code base is composed of three components:
* The data (suffix arrays). Data must be stored in a directory call ProtengineR3/ in the same directory as the installation (PeptideLR/). Symbolic links are okay.
* The suffix array source code. The LR calculator uses underlying routines (written in C/C++) from Protengine (the suffix array and related code). Protengine must be installed/compiled exactly as it is described in the *Setup* section.

## Funding

This research is based upon work supported in part by the Office of the Director of National Intelligence (ODNI), Intelligence Advanced Research Projects Activity (IARPA), via contract number 2018-18041000003. The views and conclusions contained herein are those of the authors and should not be interpreted as necessarily representing the official policies, either expressed or implied, of ODNI, IARPA, or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprints for governmental purposes notwithstanding any copyright annotation therein.

