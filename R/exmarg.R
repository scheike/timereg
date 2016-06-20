
ExMarg <- function(Y0,XX0,W0,dS0,midx1=seq(ncol(XX0)/2),midx2=seq(ncol(XX0)/2)+ncol(XX0)/2,eqmarg=TRUE,allmarg=FALSE,Z0=NULL,id=NULL,...) {  
  ii1 <- which(is.na(Y0[,2]) & !is.na(Y0[,1]))
  ii2 <- which(is.na(Y0[,1]) & !is.na(Y0[,2]))
  ii0 <- which(is.na(Y0[,1]) & is.na(Y0[,2]))
  margidx <- c(ii1,ii2)
  id1 <- id2 <-  NULL
  both <- setdiff(seq(nrow(Y0)),c(ii1,ii2,ii0))
  idB <- seq_len(length(both))
  if (allmarg) {
    id1 <- c(seq_len(length(ii1))+length(idB), idB)
    id2 <- c(seq_len(length(ii2))+length(idB)+length(id1), idB)
    ii1 <- c(ii1,both)
    ii2 <- c(ii2,both)
  }
  Y0_marg <- XX0_marg <- X0_marg1 <- X0_marg2 <- dS0_marg <- W0_marg <- NULL
  id0 <- id
  idmarg0 <- NULL
  if (length(margidx)>0) {
      Y0_marg <- cbind(c(Y0[ii1,1],Y0[ii2,2]))
      idmarg0 <- c(id[ii1],id[ii2])     
      X0_marg1 <- XX0[ii1,midx1,drop=FALSE]
      X0_marg2 <- XX0[ii2,midx2,drop=FALSE]
      dS0_marg <- dS0[,1,drop=FALSE]
      if (eqmarg) {
          XX0_marg <- rbind(X0_marg1,X0_marg2)
      } else {
          XX0_marg <- XX0[c(ii1,ii2),,drop=FALSE]
      }    
      if (!is.null(W0)) {
      W0_marg <- cbind(c(W0[ii1,1],W0[ii2,2]))
      W0 <- W0[-c(margidx,ii0),,drop=FALSE]
      }
      id0 <- id[-c(margidx,ii0)]
      Y0 <- Y0[-c(margidx,ii0),,drop=FALSE]
      if (!is.null(Z0))  Z0 <- Z0[-c(margidx,ii0),,drop=FALSE]
      XX0 <- XX0[-c(margidx,ii0),,drop=FALSE]    
  }
  
  res <- list(Y0=Y0,XX0=XX0,W0=W0,
              Y0_marg=Y0_marg, XX0_marg=XX0_marg,
              X0_marg1=X0_marg1, X0_marg2=X0_marg2,
              dS0_marg=dS0_marg, W0_marg=W0_marg,              
              id=idB, idmarg=c(id1,id2),
              ii1=ii1,
              id0=id0,idmarg0=idmarg0, ## Original id's
              margidx=margidx,
              Z0=Z0)
}
