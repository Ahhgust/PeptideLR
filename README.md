#

The package herein can be used to compute a [semi-continuous likelihood ratio](about.md) \
It can also be used to compute some genomic properties of peptides. Namely, where in the genome they might be found. Of importance, this query is done on a population database; so long as the allele is present (in at least one individual in the database) it is query-able.\

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
User-level installation (installs to $HOME/bin/, which is assumed to be in your $PATH): \
Note this includes downloading ~4000 whole proteomes (twice), which takes a while!\
Maybe it's time to grab a coffee?
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

Note, this program can be installed at the system level, but the appropriate (user/group/other-readable) location must be used to install the code, and the lrWrapper.py in the bin (e.g., /usr/local/bin) **needs to be a symlink** to that same file in the same directory!


## Installation notes
The code base is designed to work on a single file system on a single operating system (ie, at most, one binary, Windows or *nix is supported). This is typically the case, but some more exotic situations (two operating systems mounting the same file system) are not handled<br>
The code base is composed of three components:
* The data (suffix arrays). Data must be stored in a directory call ProtengineR3/ in the same directory as the installation (PeptideLR/). Symbolic links are okay.
* The suffix array source code. The LR calculator uses underlying routines (written in C/C++) from Protengine (the suffix array and related code). Protengine must be installed/compiled exactly as it is described in the *Setup* section.
* The LR calculator itself. The code for this is written in Python with minimal dependencies (v3.*, though note it does need numpy)
## Funding

This research is based upon work supported in part by the Office of the Director of National Intelligence (ODNI), Intelligence Advanced Research Projects Activity (IARPA), via contract number 2018-18041000003. The views and conclusions contained herein are those of the authors and should not be interpreted as necessarily representing the official policies, either expressed or implied, of ODNI, IARPA, or the U.S. Government. The U.S. Government is authorized to reproduce and distribute reprints for governmental purposes notwithstanding any copyright annotation therein.

