#!/usr/bin/env Rscript

### for RC-bray calculation using mutiple processors, developed by Jianshu Zhao (jianshu.zhao@gatech.edu) and Daliang Ning (ningdaliang@gmail.com)
## Collect arguments
args <- commandArgs(TRUE)

## Print help message when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      RC_bray_multip.r
      
      Arguments:
      --input_file=path/to/file    - (relative) path to file containing OTU table, must be in TAB seperated format, txt or tsv, columns are samples, rows are OTUs/ASVs, no taxonomy
      --output_file=path/to/file    - (relative) path to output file [default='output.txt']
      --processors                - numerical, number of processors to use for parallel computing
      --help                      - print this text
      
      Example:
      ./RC_bray_multip.r --input_file=OTU_table.txt --output_file=output.txt --processors=8 \n\n")
  
  q(save="no")
}
### check and install pacakges needed
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')
library(BiocManager,quietly=TRUE)
requiredPackages = c('vegan','parallel')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE, quietly=TRUE)) BiocManager::install(p,update = FALSE,quietly=TRUE)
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
args.df <- as.data.frame(do.call("rbind", parseArgs(args)))

args.list <- as.list(as.character(args.df$V2))
names(args.list) <- args.df$V1

num_processor <- as.numeric(args.list$processors)

## Arg default
if(is.null(args.list$input_file)) {
  stop("input file must be supplied (OTU table).\n", call.=FALSE)
}
if(is.null(args.list$output_file)) {
  stop("Output file must be supplied (output file).\n", call.=FALSE)
}
# Print args list to STDOUT
if(length(args) > 1) {
  for( i in names(args.list) ) {
    cat( i, "\t", args.list[[i]], "\n")
  }
}

# these variables are passed to the workflow
input.path <- normalizePath( args.list$input_file )
output.file <- ifelse( is.null(args.list$output_file), "output", args.list$output_file )

if (!file.exists(output.file)){
cat( output.file , "will be created", "\n")
} else {
    stop("file already exists! If you do not want to overwrite it, please use a new file path.\n", call.=FALSE)
}


RC.p<-function(comm,method="bray",rand=1000,nworker=num_processor)
{
  # by Daliang Ning (ningdaliang@gmail.com) on 2015.2.12 #
  # revised by Daliang according to James Stegen's code on 2015.8.5
  ## cite the original reference: Stegen JC, Lin X, Fredrickson JK, Chen X, Kennedy DW, Murray CJ et al. (2013). Quantifying community assembly processes and identifying features that impose them. Isme Journal 7: 2069-2079.##
  # comm: otu table, rownames are sample names, colnames are OTU names
  # nworker is the number of processors to use for parallel computation.
  # Note that the memory requirement will increase exponetially as you have more samples more species in you input
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
### load data and run the function
otu_all <- read.table(input.path, header = T, row.names=1, sep="\t")
dim(otu_all)
comm=t(otu_all)
RC_bray <- RC.p(comm = comm,method = "bray",rand = 1000, nworker = num_processor)
### write output files
write.table(file=output.file, RC_bray,sep = "\t")
cat("Raup-Crick metric calculation based on Bray-Curtis distance done")