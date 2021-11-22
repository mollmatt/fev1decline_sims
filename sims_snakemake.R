library(tidyverse)
library(data.table)
library(lmerTest)
# library(geepack)
# library(knitr)

# iter=1
iter=snakemake@wildcards[["iter"]]
set.seed(iter)

simPath<-"/proj/regeps/regep00/studies/COPDGene/analyses/remol/Fev1decline/sims/"

# load data
covariates<-fread(paste0(simPath,"covariates3.csv"),data.table=F)
dfsim<-fread(cmd=paste0("zcat ",simPath,"lungFx10kLong.csv.gz"),data.table=F)

# load functions
source(paste0(simPath,"simfunctions.R"))

snpnumber<-snakemake@wildcards[["snpid"]]
print(snpnumber)


if(snpnumber==1){
  snpname="snp_rd"
} else {
  snpname="snp_lb"
}


dfsim <- dfsim %>%
  mutate(snp_rd=ifelse(grp %in% c("lb","normal"),0,1),
         snp_lb=ifelse(grp %in% c("normal","rd"),0,1)) %>%
  mutate(snp=get(snpname))





# make model formulas
for(i in 1:length(covariates[,1])){

  covars<-names(covariates)[!is.na(covariates[i,])]
  covars<-covars[2:length(covars)]
  covars<-covars[!covars %in% c("randomint")]
  assign(paste0("modform",i),
         paste0("fev1~snp+",paste0(covars,collapse="+")),
         envir = .GlobalEnv)

}




# sample data
df<-sampledata(df=dfsim,n=3000,s=8,total_time = 15,agemin=20,agemax=90,proportions = c(0.35/3, 0.35/3, 0.65, 0.35/3))

df <- df %>%
  mutate(age_c2=scale(age)^2,
         time_c2=scale(time)^2,
         age_c2_baseline=scale(age_baseline)^2)


## Run models
sumTab<-data.frame(
  iteration=character(),
  modelnum=numeric(),
  model=character(),
  random_slope=character(),
  formula=character(),
  snp=character(),
  aic_qic=numeric(),
  converged=character(),
  interaction=character(),
  beta=character(),
  p=numeric()
)


for(i in 1:length(covariates[,1])){

  modform<-get(paste0("modform",i))
  slopevar<-covariates[i,"randomint"]

  lmmSlope<-lmer(as.formula(paste0(modform,"+(",slopevar,"|id)")),data = df)
  # print(summary(lmmSlope))

  convcode.lmmSlope<-ifelse(is.null(lmmSlope@optinfo$conv$lme4$code),0,1)

  if(convcode.lmmSlope==0){

    lmmSlope2<-update(lmmSlope,  start=getME(lmmSlope, c("theta","fixef")),  control=lmerControl(optCtrl=list(maxfun=2e4)))
    convcode.lmmSlope<-ifelse(is.null(lmmSlope2@optinfo$conv$lme4$code),0,1)

  }

  aic.lmmSlope<-AIC(lmmSlope)

  res.lmmSlope<-as.data.frame(summary(lmmSlope)$coefficients)
  varstokeep<-rownames(res.lmmSlope)[grep(":",rownames(res.lmmSlope))]

  res.lmmSlope<-res.lmmSlope[which(rownames(res.lmmSlope) %in% varstokeep),] %>%
    mutate(beta=paste0(
      signif(Estimate,2)," (",
      signif(Estimate-1.96*`Std. Error`,2)," - ",
      signif(Estimate+1.96*`Std. Error`,2),")"
    ),
    p=signif(`Pr(>|t|)`,2))



  res<-data.frame(
    iteration=c(rep(iter,length(varstokeep))),
    modelnum=c(rep(i,length(varstokeep))),
    model=c(rep("LMM random intercept + slope",length(varstokeep))),
    random_slope=c(rep(slopevar,length(varstokeep))),
    formula=c(rep(modform,length(varstokeep))),
    snp=c(rep(snpname,length(varstokeep))),
    aic_qic=c(rep(aic.lmmSlope,length(varstokeep))),
    converged=c(rep(convcode.lmmSlope,length(varstokeep))),
    interaction=c(varstokeep),
    beta=res.lmmSlope$beta,
    p=res.lmmSlope$p
  )


  sumTab=rbind(sumTab,res)


}

sumTab <- sumTab %>% arrange(modelnum)

write.csv(sumTab,file=paste0(simPath,"sim",iter,"_snp",snpnumber,"_lmm.csv"),row.names = F)

sessionInfo()


