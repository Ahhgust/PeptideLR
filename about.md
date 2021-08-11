#

The code herein can be used to estimate a semicontinuous likelihood ratio from peptide alleles. 
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