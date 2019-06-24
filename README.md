# ness
R package to compute hypergeometric standardization and probabilistic measures of similarity for community ecology

The ness package provides three functions to compute NNESS and CNESS, two similarity measures based on the expected species shared between random draws of individuals without replacement.
    The package is a tentative to translate into R the functions provided by Eugene D. Gallagher in his fortran program COMPAH96. The theoretical background and analytical solutions can be found in the must-read COMPAH documentation (Gallagher, 1996).
    The NNESS function computes the New Normalized Expected Specis Shared (NNESS) similarity index from a random draw of NESSm individuals and provides the similarity matrix between samples.
    The CNESS function converts raw community data to hypergeometric probabilities from a random draw of NESSm individuals followed by station normalization. The Chord-distance Normalized Expected Specis Shared (CNESS) can then be obtained by computing the Euclidean distance on the transformed data (Trueblood et al., 1994, Gallagher, 1996).
    The NESSm function computes the sample size that provides the best trade-off between NESSm = 1 (giving high weight to abundant species) and NESSm= "minimum sample total" (giving high weight to rare species)
    
# References
Gallagher, E.D., 1996. COMPAH documentation. p. 65. http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.9.1334&rep=rep1&type=pdf

Trueblood, D.D., Gallagher, E.D., Gould, D.M., 1994. Three stages of seasonal succession on the Savin Hill Cove mudflat, Boston Harbor. Limnology and Oceanography 39, 1440-1454.

# Installation
You can install the package in R from github using devtools.

devtools::install_github('Lenaick-Menot/ness')

library(ness)
