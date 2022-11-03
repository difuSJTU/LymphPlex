
#Matrix with predefined clusters, for example Lymphogen as in example data
#Each column represents a gene but for the last column as cluster
#Each row represents a sample

GetCenter=function(predefined){
  prdf=predefined[,-1]
  rownames(prdf)=predefined[,1]
  colnames(prdf)[ncol(prdf)]="Cluster"
  
  prdf=predefined[,-1]
  rownames(prdf)=predefined[,1]
  colnames(prdf)[ncol(prdf)]="Cluster"
  
  Type=unique(prdf$Cluster)
  c=ncol(prdf)-1
  r=length(Type)
  center=matrix(ncol=c,nrow=r)
  
  colnames(center)=colnames(prdf)[1:c]
  rownames(center)=Type
  
  for (i in 1:r) {
    use.i=subset(prdf,Cluster==Type[i])[,1:c]
    center[i,]=colMeans(use.i)
  }
  center=cbind(rownames(center),center)
  return(center)
}

#PAM to identify specific cluster for each sample

GetCluster=function(test,center,seeds){
  tst=test[,-1]
  rownames(tst)=test[,1]
  
  cent=center[,-1]
  rownames(cent)=center[,1]
  c=nrow(cent)
   
  colnames(test)[1]="Sample"
  colnames(seeds)=c("Cluster","Seed")
 
  r=nrow(test)
  
  sds=unique(seeds[,1])
  s=length(sds)
  
  for (i in 1:r) {
    use.i=tst[i,]
    me=melt(use.i)
    colnames(me)=c("Seed","IF")
    s.i=left_join(seeds,me,by="Seed")
    judge=c()
    for (j in 1:s) {
      s.j=subset(s.i,Cluster==paste0(sds[j]))$IF
      judge[j]=sum(s.j)-length(s.j)
    }
    zero=which(judge==0)
    
    if(length(zero)==1){
      test$res[i]=sds[zero]
    } else {
      use=rbind(cent,use.i)
      pF=pam(use,c,medoids=c(1:c),stand = FALSE, cluster.only = FALSE,do.swap = F)
      clustF=data.frame(pF$clustering)
      test$res[i]=rownames(clustF)[which(clustF[-(c+1),1]==clustF[(c+1),1])]
    }
  }
  clust=cbind(test$Sample,test$res)
  return(clust)
}
