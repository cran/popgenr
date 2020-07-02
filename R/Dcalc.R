#' @title LD calculation
#'
#' @description Function for calculating the degree of linkage disequilibirum between two biallelic diploid loci
#'
#' @param AABB
#'
#' @param AaBB
#'
#' @param aaBB
#'
#' @param AABb
#'
#' @param AaBb
#'
#' @param aaBb
#'
#' @param AAbb
#'
#' @param Aabb
#'
#' @param aabb
#'
#' @return NULL
#'
#' @examples Dcalc(AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb)
#'
#' @export Dcalc

Dcalc <- function(AABB, AaBB, aaBB, AABb, AaBb, aaBb, AAbb, Aabb, aabb){
  
  # total observations
  sum <- AABB+AaBB+aaBB+AABb+AaBb+aaBb+AAbb+Aabb+aabb
  
  # guess initial haplotype frequencies
  AB <- 1/4
  Ab <- 1/4
  aB <- 1/4
  ab <- 1/4
  
  # initial probabilities of the genotypes
  pAABB <- AB*AB
  pAaBB <- 2*AB*aB
  paaBB <- aB*aB
  pAABb <- 2*AB*Ab
  pAaBb <- 2*AB*ab+2*Ab*aB
  paaBb <- 2*aB*ab
  pAAbb <- Ab*Ab
  pAabb <- 2*Ab*ab
  paabb <- ab*ab
  
  # vector of the probabilities
  x <- c(pAABB,pAaBB,paaBB,pAABb,pAaBb,paaBb,pAAbb,pAabb,paabb)
  
  # object to hold iteration information
  Diter <- matrix(nrow=10,ncol=8)
  colnames(Diter) <- c('Iter', 'AB','Ab','aB','ab','AB_phase','D','like')
  # iterations for plotting
  g <- c(1:10)
  plotAB <- list()
  plotAb <- list()
  plotaB <- list()
  plotab <- list()
  plotD <- list()
  plotlike <- list()
  
  for(i in 1:10){
    # estimate the fraction of double heterozygotes in the 'AB/ab' phase
    ABphase <- AB*ab/(AB*ab+aB*Ab)
    # calculate D
    D <- abs(AB*ab-Ab*aB) 
    # calculate the likelihood of the data
    like <- dmultinom(c(AABB,AaBB,aaBB,AABb,AaBb,aaBb,AAbb,Aabb,aabb),prob=x)
    
    # store current iteration results
    plotAB[i] <- AB
    plotAb[i] <- Ab
    plotaB[i] <- aB
    plotab[i] <- ab
    plotD[i] <- D
    plotlike[i] <- like
    
    # save current iteration results
    Diter[i,] <- c(i, AB, Ab, aB, ab, ABphase, D, like)
    
    # update haplotype frequency estimates
    AB <- (2*AABB+AABb+AaBB+AaBb*ABphase)/(2*sum)
    Ab <- (2*AAbb+AABb++Aabb+AaBb*(1-ABphase))/(2*sum)
    aB <- (2*aaBB+aaBb+AaBB+AaBb*(1-ABphase))/(2*sum)
    ab <- (2*aabb+aaBb+Aabb+AaBb*ABphase)/(2*sum)
    
    # update genotype frequency estimates
    pAABB <- AB*AB
    pAaBB <- 2*AB*aB
    paaBB <- aB*aB
    pAABb <- 2*AB*Ab
    pAaBb <- 2*AB*ab+2*Ab*aB
    paaBb <- 2*aB*ab
    pAAbb <- Ab*Ab
    pAabb <- 2*Ab*ab
    paabb <- ab*ab
    # save genotype frequencies in a vector for the multinomial
    x <- c(pAABB,pAaBB,paaBB,pAABb,pAaBb,paaBb,pAAbb,pAabb,paabb)
  }
  
  # plot AB frequency
  plot(g, plotAB, ylim=c(0,1), type='b', col='blue', lwd=3, xlab="iterations", 
       ylab="Haplotype frequency or D")
  # plot Ab
  par(new=TRUE)# save the last plot to overlay
  plot(g, plotAb, type="b", lwd=3, ylim=c(0,1), col="red", xlab=" ", ylab=" ")
  # plot aB
  par(new=TRUE)# save the last plot to overlay
  plot(g, plotaB, type="b", lwd=3, ylim=c(0,1), col="green", xlab=" ", ylab=" ")
  # plot ab
  par(new=TRUE)# save the last plot to overlay
  plot(g, plotab, type="b", lwd=3, ylim=c(0,1), col="orange", xlab=" ", ylab=" ")
  # plot D
  par(new=TRUE)# save the last plot to overlay
  plot(g, plotD, ylim=c(0,1), type='b', col='gray', lwd=3, xlab=" ", ylab=" ")
  ## Legend
  legend('topleft', bg='white', legend=c('AB','Ab','aB','ab','D'),
         col=c('blue','red','green','orange','gray'), lty=1, lwd=3)
  
  #print(Diter)
}