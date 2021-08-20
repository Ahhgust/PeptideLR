#

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
unzip ProtengineR3.zip
for file in *.tbz; do tar -xf $file && rm $file & done
wait
cd ..
if [ -d $HOME/bin ]; then ln -s $PWD/lrWrapper.py $HOME/bin; fi
```
<br>

Note, this program can be installed at the system level, but the appropriate (user/group/other-readable) location must be used to install the code, and the lrWrapper.py in the bin (e.g., in /usr/bin or /usr/local/bin) **needs to be a symlink**!

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

`lrWrapper.py -d peptides.tsv -l HG38_Clean -g`

<br>
This command takes a little while to complete. What it is doing is searching every protein sequence in every individual from the 1000 Genomes Project (HG38_Clean), and checking to see which individuals have a tryptic peptide that is equivalent (remembering that Is and Ls are mass-equivalent) to the ones queried. Of those individuals, it then pulls out the information on which proteins correspond to these peptides.

<br>
When the above completes it will create a directory (`peptides/`). In it will be a file (`genomicInformation.tsv`) that looks like:

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






## Installation notes and other gotchas
The code base is designed to work on a single file system on a single operating system (ie, at most, one binary, Windows or *nix is supported). This is typically the case, but some more exotic situations (two operating systems mounting the same file system) are not handled<br>
The code base is composed of three components:
* The data (suffix arrays). Data must be stored in a directory call ProtengineR3/ in the same directory as the installation (PeptideLR/). Symbolic links are okay.
* The suffix array source code. The LR calculator uses underlying routines (written in C/C++) from Protengine (the suffix array and related code). Protengine must be installed/compiled exactly as it is described in the *Setup* section.
* The LR calculator itself. The code for this is written in Python with minimal dependencies (v3.*, though note it does need numpy)
## Funding

This research is based upon work supported in part by the Office of the Director of National Intelligence (ODNI), Intelligence Advanced Research Projects Activity (IARPA), via contract number 2018-18041000003. The views and conclusions contained herein are those of the authors and should not be interpreted as necessarily representing the official policies, either expressed or implied, of ODNI, IARPA, or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprints for governmental purposes notwithstanding any copyright annotation therein.

