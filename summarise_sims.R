library(tidyverse)
library(data.table)
library(optparse)

# option_list = list(
#   make_option(c("-i", "--infile"), type = "character", default = NULL,
#               help = "Path to file containing input data",
#               metavar = "character"),
#   make_option(c("-o", "--outfile"), type = "character", default = NULL,
#               help = "Path to output file ", metavar = "character"),
#   make_option(c("-r", "--chr_range"),type = "integer", default = NULL,
#               help = " chromosomes", metavar = "number"),
#   make_option(c("-p", "--phenovar"), type = "character", default = NULL,
#               help = "Phenovar", metavar = "character"),
#   make_option(c("-v", "--covariates"), type = "character", default = NULL,
#               help = "Covariates", metavar = "character")
# )
# 
# # parse arguments
# opt_parser = OptionParser(option_list = option_list)
# opt = parse_args(opt_parser)


simPath<-"/proj/regeps/regep00/studies/COPDGene/analyses/remol/Fev1decline/sims/"

sims <- fread(paste0(simPath,"sims_lmm_compiled.csv"),data.table=F)

# head(sims %>% filter(modelnum==1),20)

snpnumber<-snakemake@wildcards[["snpid"]]
print(snpnumber)
# snpnumber<-opt$snpid


if(snpnumber==1){
  snpname="snp_rd"
} else {
  snpname="snp_lb"
}



# for(snpname in c("snp_rd","snp_lb")){
  
sims <- sims %>% 
  filter(snp==snpname) 

  simSum<-sims %>% 
    dplyr::group_by(modelnum) %>% 
    summarise(avgAIC=mean(aic_qic,na.rm=T),
              minAIC=min(aic_qic,na.rm=T),
              pctConverged=sum(converged)/length(converged)*100)
  
  simSum<-as.data.frame(simSum)
  # kable(simSum)
  
  
  # get best AIC models and attach
  sim.best <- sims[0,]
  
  for(i in 1:length(simSum[,1])){
    
    bestAIC<-as.numeric(simSum[i,"minAIC"])
    simres<-sims %>% 
      filter(modelnum==i) %>% 
      filter(aic_qic==bestAIC)
    
    sim.best<-rbind(sim.best,simres)
  }
  
  sim.best.merge<-merge(sim.best,
                        simSum,
                        by="modelnum",all.x=TRUE)
  
  sim.best.merge <- sim.best.merge %>% 
    dplyr::select(modelnum,iteration,formula,random_slope,interaction,snp,converged,pctConverged,beta,p,avgAIC,minAIC)
  
  # Now arrange by model but report the proportion of times a model had the best AIC and the percentage of times it converged
  bestmodlist<-list()
  
  for(j in 1:length(unique(sims$iteration))){
    
    iter=unique(sims$iteration)[j]
    
    aicTab <- sims %>% 
      filter(iteration==iter) 
    
    bestmodlist[j]<-as.numeric(unique(aicTab[which(aicTab$aic_qic==min(aicTab$aic_qic)),"modelnum"]))
  }
  
  bestmodlist<-as.numeric(as.character(bestmodlist))
  
  ## add to table
  pctBestAIC<-list()
  
  for(k in 1:length(unique(sim.best.merge$modelnum))){
    
    mnum <- unique(sim.best.merge$modelnum)[k]
    pctBestAIC[k]<-sum(mnum==bestmodlist)/length(bestmodlist)*100
  }
  
  pctBestAIC<-as.numeric(as.character(pctBestAIC))
  
  # check it adds to 100%
  sum(pctBestAIC)
  
  bestAicTab=data.frame(
    modelnum=1:length(unique(sim.best.merge$modelnum)),
    pctBestAIC=pctBestAIC
  )
  
  sim.best.merge <- merge(sim.best.merge,
                          bestAicTab,
                          by="modelnum",all.x = TRUE)
  
  # kable(sim.best.merge)
  
  write.csv(sim.best.merge,file=paste0(simPath,"sim_summary_snp",snpnumber,".csv"),row.names = F)
  
  
# }

sessionInfo()
