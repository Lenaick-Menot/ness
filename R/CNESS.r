#'
#'The Chord-distance Normalized Expected Specis Shared (CNESS)
#'
#'
#'The function converts raw community data to hypergeometric probabilities from a random draw of NESSm individuals followed by station normalization. The CNESS distance can then be obtained by computing the Euclidean distance on the transformed data (see Trueblood et al., 1994; Gallagher 1996).
#'
#'
#' @param X Community data matrix, samples as rows, species as column
#'
#' @param m NESSm value, can range from 1 (high weight to adundant species) to the minimum sample total (high weight to rare species). The NESSm function compute the best trade-off.
#'
#' @return Standardized data frame
#'
#' @examples
#'library(vegan)
#'data(varespec)
#'m <- NESSm(varespec)
#'cness.mat <- CNESS(varespec, m)
#'cness.dist <- dist(cness.mat, "euclidean")
#'
#' @references
#' Gallagher, E.D., 1996. COMPAH documentation. 65 p. http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.9.1334&rep=rep1&type=pdf.
#' Trueblood, D.D., Gallagher, E.D., Gould, D.M., 1994. Three stages of seasonal succession on the Savin Hill Cove mudflat, Boston Harbor. Limnology and Oceanography 39, 1440-1454.
#'
#' @importFrom vegan decostand vegdist
#'
#' @export

CNESS <- function(X, m){
X <- round(X)
H <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
for (i in 1:nrow(X)) {
for (k in 1:ncol(X)) {
if ((sum(X[i, ]) - X[i, k] - m) < 0) H[i, k] <- 1
else if (is.finite(choose((sum(X[i,])-X[i,k]), m)/choose(sum(X[i,]), m)) == FALSE) H[i, k] <- 1
else H[i, k] <- 1 -(choose((sum(X[i,])-X[i,k]), m)/choose(sum(X[i,]), m))
k <- k+1} # fin de boucle pour esp?ce
i <- i+1} # fin de boucle pour sample
normH <- decostand(H, "normalize")
row.names(normH) <- row.names(X)
colnames(normH) <- colnames(X)
return(normH)
}




