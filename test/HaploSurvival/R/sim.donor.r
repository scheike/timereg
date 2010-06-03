sim.donor<-function(gp,gs,haps)
{
  ##print(setup)

  n<-nrow(gp); nPossHaps <- sapply(gs$HPIordered,length);
  inames<-sapply(gs$uniqueHaploNames,function(x) as.numeric(gsub(",",".",x)))

  Zh<-c(); hapD<-rep(0,4); 
  for (i in 1:n) 
    {
      if (nPossHaps[i]==1)  Zh<-c(Zh,1); 
      if (nPossHaps[i]>1)  {
        ppp<-rep(0.25*0.25,nPossHaps[i])
        pgeno<-sum(ppp) 
        ri<-sample(1:nPossHaps[i],1,prob= ppp/pgeno)

        hapo<- gs$HPIordered[[i]][[ri]]
        hapo1<-gs$uniqueHaploNames[hapo[1]]
        hapo2<-gs$uniqueHaploNames[hapo[2]]
                                        #print(setup$HPIordered[[i]][[ri]])
        hapD[c(1,3)]<-haps[hapo[1],]
        hapD[c(2,4)]<-haps[hapo[2],]
        match11<-sum(gp[i,c(1,3)]==hapD[c(1,3)])==2 
        match12<-sum(gp[i,c(1,3)]==hapD[c(2,4)])==2 
        match21<-sum(gp[i,c(2,4)]==hapD[c(1,3)])==2 
        match22<-sum(gp[i,c(2,4)]==hapD[c(2,4)])==2 
        match<-match11*match22+match12*match21; 
        if (match) Zh<-c(Zh,1) else Zh<-c(Zh,0)
      }
      ##print("==============================")
    } ## }}}
  return(Zh)
}
