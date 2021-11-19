
sampledata<-function(df,n,s,total_time,agemin,agemax,proportions){
  
  
  # df = data frame
  # n = number of unique individuals
  # s = number of time points to sample
  # total_time = total time of the study
  # agemin = minimum age for study entry
  # agemax = max age for study entry
  # proportions = proportions of individuals from each group; order is lb,lbrd,normal,rd
  
  # Parameters to vary: proportions of groups, baseline fev1 
  # set.seed(1010)
  
  # config
  # n=1800
  # s=6 # number time points to sample
  # total_time=10 # number of years in study - must be at least as many as the number of time points
  # agemin=40
  # agemax=80
  # 
  # proportions<-c(0.5,0.25,0.125,0.125) # order is lb,lbrd,normal,rd
  
  df.select<-df %>% 
    filter(age >= agemin) %>% 
    filter(age <= agemax)
  
  
  grpnums<-n*proportions
  grps=table(df.select$grp)
  names(grps)
  
  if(any(length(unique(df.select[which(df.select$grp==names(grps)[1]),"id"]))<grpnums)){
    
    maxsample<-length(unique(df.select[which(df.select$grp==names(grps)[1]),"id"]))
    grpnums<-floor(grpnums/max(grpnums)*maxsample)
    
  } 
  
  idx1<-sample(unique(df.select[which(df.select$grp==names(grps)[1]),"id"]),grpnums[1])
  idx2<-sample(unique(df.select[which(df.select$grp==names(grps)[2]),"id"]),grpnums[2])
  idx3<-sample(unique(df.select[which(df.select$grp==names(grps)[3]),"id"]),grpnums[3])
  idx4<-sample(unique(df.select[which(df.select$grp==names(grps)[4]),"id"]),grpnums[4])
  
  
  
  
  ec<-df.select[c(idx1,idx2,idx3,idx4),]
  
  
  idlist<-ec$id
  
  ec.long<-ec[0,]
  
  # generate age distribution based on weibull sampling
  ages <- unique(df.select[which(df.select$age <= (agemax-(s+1))),"age"])
  wbprobs <- rweibull(ages,shape=2,scale=1)
  agedist=sample(ages,length(ages),replace=TRUE,prob=wbprobs) 
  
  for(i in 1:length(idlist)){
    
    idnum=idlist[i]
    

    age_baseline=sample(agedist,1,replace = TRUE)
      
    df.filt<-df.select %>% 
      filter(id==idnum) %>% 
      filter(age >= age_baseline) %>% 
      filter(age <= (age_baseline + total_time))
    
    # ec.bl<-ec[which(ec$id==idnum),]
    
    idx.ind<-sample(1:nrow(df.filt),s,replace = FALSE)
    
    ec.new<-df.filt[idx.ind,] %>% arrange(age)
    
    # ec.new<-rbind(ec.bl,
    #               df.filt)
    
    ec.new$time=seq(1,s,by=1)
    ec.new$age_baseline=age_baseline
    
    ec.long<-rbind(ec.long,ec.new)
    
  }
  
  ec.long
}



rangestd <- function(x){(x-min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm = T))} # range standardization
