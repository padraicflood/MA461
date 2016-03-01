require(expm)

#input sequences and split them for looping.
x <- unlist(strsplit('AGTCCATGAAT',split=''))
y <- unlist(strsplit('ACGTCGTGCCA',split=''))

genetic_distance_mle <- function(s1, s2, upper = 10){
  #input: two aligned sequences. 
  #       (optional) upper bound for search
  #output: Maximum likelihood estimate of genetic distance
  G = diag(x = -1 -1/3, 4,4) + matrix(1/3,4,4) #generator matrix for JC model.
  P <- function(t){
    Pt <- expm(G*t)
    rownames(Pt) <- c('A', 'G', 'C', 'T')
    colnames(Pt) <- c('A', 'G', 'C', 'T')
    product <- 1
    for (i in 1:length(s1)){
      product <- product * Pt[s1[i], s2[i]]
    }
    product
  }
  optimise(P, c(0, upper), maximum = T)$maximum
}

#Or using formula:
p <- sum(x != y)/length(x) #p-distance
d <- -(3/4)*log(1-(4/3)*p) #genetic distance
