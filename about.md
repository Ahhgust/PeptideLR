#

The code herein can be used to estimate a semicontinuous likelihood ratio from peptide alleles. 
It is directly derived from the forensic literature. In particular, it is based on [LoComatioN](https://doi.org/10.1016/j.forsciint.2006.04.016) of Gill et al., 
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

We need three pieces of information to compute a likelihood: \
* A set of peptides
* A reference profile (from some hypothesized contributor; to ask if the contributor "matches" the evidence")
* A population database (to estimate the rarity of such a "match" in populations)

Starting with the first bullet, the peptides should be variable in the population (genetically variable peptides or GVPs). The code itself is designed to work solely on autosomal GVPs. GVPs may be linked (provide redundant information). This code permits linkage (up to the level of the chromosome). However, you must tell it which chromosome each GVP comes from. In practice you can choose a finer level of linkage (units of Morgans would also make sense), however for this presentation we control for linkage at the level of the chromosome. \
<br>

The reference profile is also needed. This is the full proteome of some individual as estimated from whole genome/exome sequencing. Click [here](creation.md) to see how to make such a file. \
<br>

The last piece of information is the population database. This database was created in the same way as the reference profile, however it was done so on population database. In this case it is the publically available genotypes from the [gnomAD project](https://gnomad.broadinstitute.org/about) [(raw data)](https://gnomad.broadinstitute.org/downloads#v3-hgdp-1kg) which includes the [1000 Genomes Project](https://en.wikipedia.org/wiki/1000_Genomes_Project) and samples from the [Human Genome Diversity Panel](https://en.wikipedia.org/wiki/Human_Genome_Diversity_Project)
