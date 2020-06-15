RC.p<-function(comm,method="bray",rand=1000,nworker=8)
{
  # by Daliang Ning (ningdaliang@gmail.com) on 2015.2.12 #
  # revised by Daliang according to James Stegen's code on 2015.8.5
  ## cite the original reference: Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ et al. (2013). Quantifying community assembly processes and identifying features that impose them. Isme Journal 7: 2069-2079.##
  # comm: otu table, rownames are sample names, colnames are OTU names
  # nworker is the thread number for parallel computation.
  
  library(vegan)
  library(parallel)
  com<-comm[,colSums(comm)>0]
  BC.obs<-as.matrix(vegdist(com,method=method))
  com.rd0=com
  com.rd0[]=0
  id<-(1:ncol(com))
  prob.sp<-colSums(com>0)
  prob.ab<-colSums(com)
  Si<-rowSums(com>0)
  Ni<-rowSums(com)
  samp.num=nrow(com)
  
  BC.rand<-function(j,com.rd0,samp.num,id,prob.sp,prob.ab,Si,Ni,method)
  {
    library(vegan)
    com.rd=com.rd0
    for(i in 1:samp.num)
    {
      id.sp<-sample(id,Si[i],replace=FALSE,prob=prob.sp)
      if(length(id.sp)==1){count=rep(id.sp,Ni[i])}else{
        count<-sample(id.sp,(Ni[i]-Si[i]),replace=TRUE,prob=prob.ab[id.sp])
      }
      table<-table(count)
      com.rd[i,as.numeric(names(table))]=as.vector(table)
      com.rd[i,id.sp]=com.rd[i,id.sp]+1
    }
    BCrand=as.matrix(vegdist(com.rd,method=method))
    BCrand
  }
  
  c1<-makeCluster(nworker,type="PSOCK")
  message("Now parallel computing. begin at ", date(),". Please wait...")
  BC.rd<-parLapply(c1,1:rand,BC.rand,com.rd0=com.rd0,samp.num=samp.num,id=id,prob.sp=prob.sp,prob.ab=prob.ab,Si=Si,Ni=Ni,method=method)
  stopCluster(c1)
  
  BC.rd=array(unlist(BC.rd),dim=c(nrow(BC.rd[[1]]),ncol(BC.rd[[1]]),length(BC.rd)))
  gc()
  
  comp<-function(x,c){(x<c)+0.5*(x==c)}
  message("----now calculating rc at ",date(),"----")
  alpha=matrix(rowSums(apply(BC.rd,3,comp,c=BC.obs)),nrow=nrow(BC.obs))/rand
  rc=(alpha-0.5)*2
  rownames(rc)=rownames(BC.obs)
  colnames(rc)=colnames(BC.obs)
  rc
}

setwd("/Volumes/Jianshu/Downloads/Submit_fungal_ecology/Fungi_paddy/Paddy_fungi_macqiime_new")
fungi_all <- read.table("otu_table_mc2_w_tax_no_singletons_nolowcoverage_2500_rarefied_2500_181_1.csv", header = T, row.names=1, sep=",")
dim(fungi_all)
otu_fungi=fungi_all[,7:ncol(fungi_all)]
otu_fungi=fungi_all
comm=t(otu_fungi)


RC.p(comm = comm,method="bray",rand=1000)
