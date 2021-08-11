#

peptideLR.py can be used to estimate a semicontinuous likelihood ratio from peptide alleles. 
It is designed based on [LoComatioN](https://doi.org/10.1016/j.forsciint.2006.04.016) of Gill et al., 
which in turn is based on the set theoretic (semicontinuous) likelihood ratio as per [Curran et al](https://doi.org/10.1016/j.forsciint.2004.04.077) <br>

The basic premise is that there is evidence (*E*), and the probability of this evidence can be considered given some hypothesis (*Pr(E|H)*).
For example, one hypothesis (*H<sub>1</sub>*) might be that the evidence is from some individual, while an alternative hypothesis (*H<sub>2</sub>*) may be that instead some random person deposited the evidence. The relative strength of the evidence given these two competing hypotheses can be evaluated with a likelihood ratio (*LR*) taken as:


![\Large LR=\frac{Pr(E|H1}{2a}](https://latex.codecogs.com/svg.latex?LR%3D%5Cfrac%7BPr(E|H_1)%7D%7BPr(E|H_2)%7D)

The evidence (*E*) is, in this case, some set of peptide markers. <br><br>
In this approach all peptide variations that are considered are those that derive from changes that are genetic. Namely, that there is some set of alleles (in DNA) that give rise to some specific peptide form. Even with this constraint, applying a LR to peptides requires some modifications and extensions from that of Curran et al. <br> Namely
* Peptides need not be independent
  * The same peptide may overlap multiple SNP sites
  * Two peptides may be adjacent, or more generally may be in genetic linkage, thus they may convey redundant information
  

The peptides in question are assumed to have been derived from some set of proteins, proteins that have been in turn digested with some enzyme like [trypsin](https://en.wikipedia.org/wiki/Trypsin). <br>
The above formulation is assessed by considering alleles

# Dependencies (pretty standard)
* For Windows users
  * [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) (WSL)
  * Note: WSL is not a strict requirement for the code to run, but the installation instructions need to be made Windows-compatible 
* g++
  * pthreads (*nix OR WSL)
  * zlib1g-dev (zlib, development version)
  * mthreads (Windows only: ie, win32 threads. This should be standard)
* Python 3.*
  * Numpy
  * joblib

<br>

# Quick start for *nix including [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about) (WSL)
User-level installation:
<br>
```
git clone --recursive https://github.com/Ahhgust/PeptideLR.git # download this project and Dr. Crysup's suffix array project
cd PeptideLR/protengine # and build
chmod +x Buildit_x64_linux.sh && ./Buildit_x64_linux.sh && ln -s $PWD/pywrap/suffer.py ..
# and grab the public datasets (forgive the dropbox link but these files are BIG!)
cd ..
mkdir ProtengineR3 && cd ProtengineR3 && wget -O ProtengineR3.zip 'https://www.dropbox.com/sh/xp3wzs5fy9taqvl/AADjmPYTT201_MPtNXjVEZuaa?dl=1'
unzip ProtengineR3.zip
for file in *.tbz; do tar -xf $file & done
wait
cd ..
if [ -d $HOME/bin ]; then ln -s $PWD/lrWrapper.py $HOME/bin; fi
```
<br>
Note, this program can be installed at the system level, but the appropriate (user/group/other-readable) location must be used to install the code, and the lrWrapper.py in the bin (e.g., /usr/local/bin) *needs to be a symlink* to that same file in the same directory!


## Installation notes
The code base is designed to work on a single file system on a single operating system (ie, at most, one binary, Windows or *nix is supported). This is typically the case, but some more exotic situations (two operating systems mounting the same file system) are not handled<br>
The code base is composed of three components:
* The data (suffix arrays). Data must be stored in a directory call ProtengineR3/ in the same directory as the installation (PeptideLR/). Symbolic links are okay.
* The suffix array source code. The LR calculator uses underlying routines (written in C/C++) from Protengine (the suffix array and related code). Protengine must be installed/compiled exactly as it is described in the *Setup* section.
* The LR calculator itself. The code for this is written in Python with minimal dependencies (v3.*, though note it does need numpy)
## Funding

This research is based upon work supported in part by the Office of the Director of National Intelligence (ODNI), Intelligence Advanced Research Projects Activity (IARPA), via contract number 2018-18041000003. The views and conclusions contained herein are those of the authors and should not be interpreted as necessarily representing the official policies, either expressed or implied, of ODNI, IARPA, or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprints for governmental purposes notwithstanding any copyright annotation therein.

