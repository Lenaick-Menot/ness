#'
#' New Normalized Expected Specis Shared (NNESS)
#'
#'
#' This functions computes the NNESS similarity index from a random draw of NESSm individuals (see Gallagher, 1996).
#'
#'
#'
#' @param X Community data matrix, samples as rows, species as column
#'
#' @param m NESSm value, can range from 1 (high weight to adundant species) to the minimum sample total (high weight to rare species). The NESSm function compute the best trade-off.
#'
#' @return distance object, the similarity matrix between samples
#'
#'  @examples
#'library(vegan)
#'data(varespec)
#'m <- NESSm(varespec)
#'nness.dist <- NNESS(varespec, m)
#'
#' @references
#' Gallagher, E.D., 1996. COMPAH documentation. 65 p. http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.9.1334&rep=rep1&type=pdf.
#'
#' @importFrom stats as.dist
#'
#' @export

NNESS <- function(X, m){
X = round(X)
H <- matrix(NA, nrow = nrow(X), ncol = ncol(X))
row.names(H) = row.names(X)
i <- 1
for (i in 1:nrow(X)) {
for (k in 1:ncol(X)) {
if ((sum(X[i, ]) - X[i, k] - m) < 0) H[i, k] <- 1
else H[i, k] <- 1 -(choose((sum(X[i,])-X[i,k]), m)/choose(sum(X[i,]), m))
k <- k+1}
i <- i+1}
ESS = H %*% aperm(H)
NewNESS = ESS
for (i in 1:nrow(ESS)){
  for (j in 1:ncol(ESS)){
NewNESS[i,j]= ESS[i,j]/(0.5*(ESS[i,i] + ESS[j,j]))
}}
res <- as.dist(NewNESS)
return(res)
}




