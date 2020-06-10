#

peptideLR.py can be used to estimate a semicontinuous likelihood ratio from peptide alleles. 
It is designed based on [LoComatioN](https://doi.org/10.1016/j.forsciint.2006.04.016) of Gill et al., 
which in turn is based on the set theoretic (semicontinuous) likelihood ratio as per [Curran et al](https://doi.org/10.1016/j.forsciint.2004.04.077) <br>

The basic premise is that there is evidence (*E*), and the probability of this evidence can be considered given some hypothesis (*Pr(E|H)*).<br>
For example, one hypothesis might be that the evidence is from some individual, while an alternative hypothesis is that instead some random person deposited the evidence


![\Large LR=\frac{Pr(E|H1}{2a}](https://latex.codecogs.com/svg.latex?LR%3D%5Cfrac%7BPr(E|H_1)%7D%7BPr(E|H_2)%7D)





# Quick start
<br>
git clone --recursive https://github.com/Ahhgust/PeptideLR.git
<br>

