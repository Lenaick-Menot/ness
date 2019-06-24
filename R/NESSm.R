#'
#' NESSm sample size
#'
#' Computes the sample size that provides the best trade-off between NESSm = 1 (giving high weight to abundant species) and NESSm= "minimum sample total" (giving high weight to rare species), after Trueblood et al. (1994) and Gallagher (1996).
#'
#'
#'
#'
#' @param X Community data matrix, samples as rows, species as column
#'
#' @return sample size to be used in CNESS or NNESS to produce an index sensitive to both rare and abundant species
#'
#' @examples
#'require(vegan)
#'data(varespec)
#'m <- NESSm(varespec) # Be patient, computation is slow
#'
#'### Note that the function doesn't work with large numbers nor data matrix with no variance
#' largenumber <- matrix(rnorm(100, 1000, 100), 10 , 10)
#' m <- NESSm(largenumber)
#'
#' @references
#' Gallagher, E.D., 1996. COMPAH documentation. 65 p. http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.9.1334&rep=rep1&type=pdf.
#' Trueblood, D.D., Gallagher, E.D., Gould, D.M., 1994. Three stages of seasonal succession on the Savin Hill Cove mudflat, Boston Harbor. Limnology and Oceanography 39, 1440-1454.
#'
#' @importFrom vegan decostand vegdist
#'
#' @importFrom pcaPP cor.fk
#'
#'@export

NESSm <- function(X){

# The function works only with integer
X <- round(X)

##Define the range of m values (i.e. from 1 to mmax = the minimu sample total)
mmax <- Inf
for (a in 1:nrow(X)) {
  if (sum(X[a, ]) < mmax) mmax <- sum(X[a, ])
  }
m <- c(1:mmax)

# compute CNESS for all values of m
H <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
nbligne <- (nrow(X)^2 - nrow(X))/2
test.tau <- matrix(data=NA, nr=length(m), nc=nbligne)
for (E in 1:length(m)) {
    for (i in 1:nrow(X)) {
        for (k in 1:ncol(X)) {
      if ((sum(X[i, ]) - X[i, k] - m[E]) < 0) H[i, k] <- 1
      else if (is.finite(choose((sum(X[i,])-X[i,k]), m[E])/choose(sum(X[i,]), m[E])) == FALSE) H[i, k] <- 1
      else H[i, k] <- 1 -(choose((sum(X[i,])-X[i,k]), m[E])/choose(sum(X[i,]), m[E]))
          }
      }
  normH <- decostand(H, "normalize")
  distH <- vegdist(normH, method="euclidean", diag=TRUE, upper=FALSE)
  test.tau[E, ] <- as.vector(distH)
  }

# Kendall tau correlations between all pairs of distance matrix for m=1 to m=mmax
kendall <- cor.fk(t(test.tau))

# Which value of m has a correlation that is roughly the same between CNESS(NESSm=1) and CNESS(NESSm=mmax)?
mintau <- numeric(length(m))
for (j in 1:length(m)) {
  mintau[j] <- abs(kendall[j, 1] - kendall[length(m), j])
}
N <- m[which.min(mintau)]
return(N)
}

