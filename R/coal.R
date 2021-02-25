#' @title Simple coalescent simulator
#'
#' @description Function that provides a simple starting off point to simulate a coalescent process
#'
#' @param length
#'
#' @param number
#'
#' @param muscale
#'
#' @param reps
#'
#' @param prnt
#'
#' @return NULL
#'
#' @examples coal(length, number, muscale, reps, prnt)
#'
#' @export coal

coal <- function(length, number, muscale, reps, prnt=0){

  # create lists to keep track of variables
  thetawlist<-numeric()
  pilist<-numeric()


  for(r in 1:reps){
    seq=matrix(0,number,length) # need to reset each rep
    # pick first coalescent time of two lineages
    n=2
    t1=rexp(1,1/2)

    # pick number of mutations on each lineage
    mut1=rpois(2,t1*muscale)

    # add mutations
    for(i in 1:mut1[1]){
      if(mut1[1]==0){break}
      pos=sample(1:length, 1)
      seq[1,pos]=1
    }
    for(i in 1:mut1[2]){
      if(mut1[2]==0){break}
      pos=sample(1:length, 1)
      seq[2,pos]=1
    }

    # make a loop for three and greater lineages
    for(n in 3:number){
      pick=sample(n-1,1)

      for(i in 1:length){
        seq[n,i]=seq[pick,i]
      }

      # times during segment
      ways=n*(n-1)/2
      t2=rexp(1,ways/2)

      # mutation number
      mut2=rpois(n,t2*muscale)

      # add mutations
      for(j in 1:n){
        for(i in 1:mut2[j]){
          if(mut2[j]==0){break}
          pos=sample(1:length, 1)
          seq[j,pos]=1
        }
      }

    }

    #if(prnt){print(seq)}

    # calculate S
    S=0
    for(i in 1:length){
      dif=0
      ref=seq[1,i]
      for(j in 1:number){
        if(seq[j,i]!=ref){
          dif=1
        }
      }
      if(dif==1){S=S+1}
    }
    if(prnt){print(paste("S: ",S))}
    # convert to Watterson's estimator of theta
    # calculate harmonic number
    harm=0
    for(i in 1:(number-1)){
      harm=harm+1/i
    }
    thetaw=S/harm
    if(prnt){print(paste("theta_w: ",thetaw))}
    # record the value of this repitition
    thetawlist[r]=thetaw

    # calculate pi (average nucleotide heterozygosity, H)
    # a.k.a. Tajima's estimator of 4Nemu, sometimes represented by pi
    # need frequency of 1 at each site
    # 1 - p^2 - (1-p)^2
    # average over all sites
    hetsum=0
    for(i in 1:length){
      count=0
      for(j in 1:number){
        if(seq[j,i]==1){
          count=count+1
        }
      }
      p=count/number
      #print(c("freq",p))
      het=1-(p*p+(1-p)*(1-p))
      #apply small sample-size correction (Nei and Roychoudhury 1974)
      het=het*number/(number-1)
      #print(c("het",het))
      hetsum=hetsum+het
    }
    pi=hetsum
    if(prnt){print(paste("pi: ",pi))}
    # record current value
    pilist[r]=pi

  }

  # plot the two theta estimates for each tree
  maximum <- max(c(thetawlist,pilist))
  plot(thetawlist,pilist,xlab="4Neu, Watterson's",ylab="4Neu, Tajima's",xlim=c(0,maximum),ylim=c(0,maximum))
  abline(0,1)

  if(prnt==1){
    # print variances and covariance
    print(paste("Variance, Taj.: ", var(pilist)))
    print(paste("Variance, Wat.:  ", var(thetawlist)))
    print(paste("Covariance: ", cov(thetawlist, pilist)))
  }
}
