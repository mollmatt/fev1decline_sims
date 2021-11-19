library(tidyverse)
library(data.table)

iter=1:1000
snpnumber=1:2

simPath<-"/proj/regeps/regep00/studies/COPDGene/analyses/remol/Fev1decline/sims/"

res<-fread(paste0(simPath,"sim","1","_snp","1","_lmm.csv"),data.table=F)

res<-res[0,]

for(i in 1:length(iter)){
  
  for(j in 1:length(snpnumber)){

    res2<-fread(paste0(simPath,"sim",i,"_snp",j,"_lmm.csv"),data.table=F)
    
    res <- rbind(res,res2)
    
    
  }
}

res

write.csv(res,file=paste0(simPath,"sims_lmm_compiled.csv"),row.names = F)

sessionInfo()

