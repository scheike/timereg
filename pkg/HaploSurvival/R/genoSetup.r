geno.setup <- function(G,returnNames=TRUE,haplo.baseline=NULL)
{
  orderedAlleles <- apply(G,2,unique)
  unorderedAlleles <- list()
  nLoci <- ncol(G)/2
  for(i in 1:nLoci){
    if(is.list(orderedAlleles)){
      unorderedAlleles[[i]] <- sort(unique(c(orderedAlleles[[2*i - 1]],orderedAlleles[[2*i]])))
    } else if(is.matrix(orderedAlleles)){
      unorderedAlleles[[i]] <- sort(unique(c(orderedAlleles[,2*i - 1],orderedAlleles[,2*i])))
    } else if(is.vector(orderedAlleles)){
      unorderedAlleles[[i]] <- sort(unique(c(orderedAlleles[2*i - 1],orderedAlleles[2*i])))
    } else {
      browser()
      stop('Do not know what orderedAlleles is...')
    }
  }
  nAllelesPerLocus <- sapply(unorderedAlleles,length)
  indexOfFirstMultiAllelicLocus <- min(which(nAllelesPerLocus>1))
  numberOfOrderedGenotypePairs <- c(nAllelesPerLocus[1:indexOfFirstMultiAllelicLocus]*(nAllelesPerLocus[1:indexOfFirstMultiAllelicLocus]+1)/2,
                                    nAllelesPerLocus[(indexOfFirstMultiAllelicLocus+1):nLoci]^2)
  if(prod(nAllelesPerLocus) > .Machine$integer.max){
    stop('Too many possible haplotypes')
  }

  Ggeneric <- matrix(0,nrow(G),ncol(G))
  for(i in 1:nLoci){
    Ggeneric[,2*i - 1] <- match(as.vector(G[,2*i - 1]),unorderedAlleles[[i]]) 
    Ggeneric[,2*i]     <- match(G[,2*i],    unorderedAlleles[[i]])
  }

  nPeople <- nrow(G)

  cumProdForOGPI <- cumprod(c(1,nAllelesPerLocus[(nLoci):2]))[nLoci:1]
  HPI <- list()
  HPN <- list()
  # Encode the possible genotype pairs for each person at each locus
  for(i in 1:nPeople){
    OGPI <- list()
    distinct <- FALSE
    for(j in 1:nLoci){
      Row <- min(Ggeneric[i,2*j+(-1:0)])
      Col <- max(Ggeneric[i,2*j+(-1:0)])
      if(Row==Col || distinct == FALSE){
##      if(Row==Col){
        OGPI[[j]] <- nAllelesPerLocus[j]*(Row-1)+Col
      } else {
        OGPI[[j]] <- c(nAllelesPerLocus[j]*(Row-1)+Col,
                       nAllelesPerLocus[j]*(Col-1)+Row)
      }
      if(distinct == FALSE){
        distinct <- Row!=Col
      }        
    }
    # Make a matrix of possible ordered genotype pairs
    M <- as.matrix(expand.grid(OGPI))
    OGPI <- list()
    OGPN <- list()
    # Encode this as haplotypes - indexed as integers
    for(j in 1:nrow(M)){
      M2 <- matrix(0,2,nLoci)
      for(k in 1:nLoci){
        M2[1,k] <- ((M[j,k] - 1) %% nAllelesPerLocus[k]) + 1
        M2[2,k] <- ((M[j,k] - 1) %/% nAllelesPerLocus[k]) + 1
      }
      OGPI[[j]] <- as.vector(M2%*%cumProdForOGPI)
      if(returnNames==TRUE){
        OGPN[[j]] <- apply(sapply(1:nLoci,
                     function(i,newNames,oldNames){unorderedAlleles[[i]][newNames[,i]]},
                     newNames=M2,oldNames=unorderedAlleles)
                     ,1, function(x){paste(x,collapse = ',')})
      }
    }
    HPI[[i]] <- OGPI
    HPN[[i]] <- OGPN    
  }

  # Record all unique possible haplotypes  
  uniqueHaplos <- sort(unique(unlist(HPI)))
  if(returnNames == TRUE){
###    uniqueHaploNames <- unique(unlist(HPN))
    uniqueHaploNames <- sort(unique(unlist(HPN)))
  } else {
    uniqueHaploNames <- NULL
  }
  
  # check for specification of baseline haplotype 
  if (is.null(haplo.baseline)==FALSE) {
     index<-(haplo.baseline==uniqueHaploNames)
     if (sum(index)==1) {
        base.index<-(1:length(uniqueHaplos))[index]
        uniqueHaplos<-c(uniqueHaplos[-base.index],uniqueHaplos[base.index])
        uniqueHaploNames<-c(uniqueHaploNames[-base.index],haplo.baseline)
     }
     else 
     print("Haplo-baseline not consistent with uniqueHaploNames")
  }

  # Relabel possible haplotypes according to this ordering
###print(HPI)
###print(HPN)
###print(uniqueHaploNames)
###print(uniqueHaplos)

  HPIordered <- HPN
  for(i in 1:length(HPIordered)){
    for(j in 1:length(HPIordered[[i]])){
      HPIordered[[i]][[j]] <- match(HPIordered[[i]][[j]],uniqueHaploNames)
    }
  }

  # simple guess on frequencies, counting the occurrences
  ini.freqs<-as.vector(table(unlist(HPIordered))/sum(table(unlist(HPIordered))))

  # return all the parts that are required to calculate the likelihood and derivatives
  out <- list(HPIordered = HPIordered,
              uniqueHaploNames = uniqueHaploNames,
              nAllelesPerLocus = nAllelesPerLocus,
              unorderedAlleles = unorderedAlleles,
              nPeople = nPeople,
              nPossHaps=sapply(HPIordered,length),
              nLoci = nLoci,
              initial.freqs=ini.freqs)

  class(out) <- "geno.setup"
  
  return(out)
}


print.geno.setup <- function(x,...){

  cat("A geno.setup object:",fill = TRUE)
  cat("  Data for ",x$nPeople," people",sep = "",fill = TRUE)
  cat("  Data for ",x$nLoci," loci",sep = "",fill = TRUE)
  cat("  Number of alleles per loci",x$nAllel,sep=" ",fill = TRUE)
  cat("  ",length(x$uniqueHaploNames)," possible haplotypes",sep = "",fill = TRUE)
  cat("    for example: ",paste(head(x$uniqueHaploNames),collapse=' ')," ...",sep = "",fill = TRUE)
  cat("  Initial estimate of frequencies: ",paste(round(head(x$initial.freqs),3),collapse=' ')," ...",sep = "",fill = TRUE)
  cat("  Number of possible haplotype pairs per person: ",sep = "",fill = TRUE)
  cat("    ",paste(head(x$nPossHaps),collapse=' ')," ...",sep = "",fill = TRUE)

}



summary.geno.setup <- function(object,...){
  print.geno.setup(object,...)
}

