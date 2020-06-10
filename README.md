#

peptideLR.py can be used to estimate a semicontinuous likelihood ratio from peptide alleles. 
It is designed based on [LoComatioN](https://doi.org/10.1016/j.forsciint.2006.04.016) of Gill et al., 
which in turn is based on the set theoretic (semicontinuous) likelihood ratio as per [Curran et al](https://doi.org/10.1016/j.forsciint.2004.04.077) <br>

The basic premise is that there is evidence (*E*), and the probability of this evidence can be considered given some hypothesis (*Pr(E|H)*).
For example, one hypothesis (*H<sub>1</sub>*) might be that the evidence is from some individual, while an alternative hypothesis (*H<sub>2</sub>*) may be that instead some random person deposited the evidence. The relative strength of the evidence given these two competing hypotheses can be evaluated with a likelihood ratio (*LR*) taken as:


![\Large LR=\frac{Pr(E|H1}{2a}](https://latex.codecogs.com/svg.latex?LR%3D%5Cfrac%7BPr(E|H_1)%7D%7BPr(E|H_2)%7D)

The evidence (*E*) is, in this case, some set of peptide markers. <br>
The peptides in question are assumed to have been derived from some set of proteins, proteins that have been in turn digested with some enzyme like [trypsin](https://en.wikipedia.org/wiki/Trypsin)



# Quick start
<br>
git clone --recursive https://github.com/Ahhgust/PeptideLR.git
<br>

