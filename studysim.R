# start of code

##### R libraries required
library(truncnorm)
library(doParallel) ##if you can run in parallel (i.e. your machine has enough cores, 
                    ## it is set up to use windows - sorry if you use something else!)
library(data.table) ## I use data.table when I can with some tidyverse where it makes sense.
                    ## feel free to adapt to what you code with
library(lme4)
library(tidyverse)
library(metafor)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

################################################################################
############################## Data Preparation ###############################
################################################################################

#tru = true effect as a log response ratio
#noise = natural interannual variation that study designs must cope with to find signal amongst noise
#cht = proportional change in treatment sites
#chc = proportional change in control sites

#set true effect sizes ranging from 4 times decrease to 4 times increase
#i.e., a log response ratio of -4 to 4.
#use a uniform distribution for true effect to see how results vary

#this requires generating a range of cht and chc values to get a uniform distribution for
#tru. I have done this below, but there may be an easier or better way of doing this.


#need to think about whether we would reasonably expect control sites to change as much?

#here we set these as response ratios (not log response ratios)
#we use these to multiply by the population abundance in the simulation later on
valcht <- seq(0.25,4,((4-0.25)/999)) #treatment sites can change by 0.25 to 4 times

valchc <- seq(0.25,4,((4-0.25)/999)) #control sites can change by 0.25 to 4 times

tru_eff <- round(seq(log(1/4),log(4/1),(log(4/1)-log(1/4))/999),digits=3)
summary(tru_eff)


#find all possible combinations and then just randomly select only those that will produce
#a given true effect size, such that the final distribution of the true effect 
#log(cht/chc) (equivalent to log(t-c after/ t-c before))
#is a roughly uniform distribution of n=100 true effect sizes

allcombs <- round(outer(X=valcht, Y=valchc, FUN=function(X,Y){log(X/Y)}),digits=3)
rightcombs <- lapply(tru_eff,function(x){which(allcombs==x, arr.ind = T)})

rowstopick <- unlist(lapply(1:1000,function(x){sample(1:nrow(rightcombs[[x]]),1)}))
rowspicked <- lapply(1:1000,function(x){
  rowsandcols <- as.numeric(rightcombs[[x]][rowstopick[x],])
  return(data.table(tru=allcombs[rowsandcols[1],rowsandcols[2]],
                    cht=valcht[rowsandcols[1]],
                    chc=valchc[rowsandcols[2]]))
})

tru_eff_runs <- do.call(rbind,rowspicked[seq(1,length(rowspicked),length(rowspicked)/100)])

#check distributions of chc, cht, and tru. tru should have roughly uniform distribution
hist(tru_eff_runs$chc)
hist(tru_eff_runs$cht)
hist(tru_eff_runs$tru,n=20)
summary(tru_eff_runs$tru)
hist(tru_eff,n=20)

#set given amount of noise or interannual variation
#there may be a better way of doing this
#this ranges from 10% of pop.abundance to 50% of pop.abundance
#set just 3 different values of noise to keep simulation 'small'
noise <- c(0.1,0.25,0.5)
summary(noise)


### put parameters together in data.table ready for simulation

params = do.call(rbind,lapply(1:length(noise),function(x){
  return(data.table(tru_eff_runs,noise=noise[x]))}))
rep.sim <- 1000 #repeat simulation 1000 times for each row of params 
                #total runs= 100 tru effect values * 3 noise values * 1000 repititions)
params.rep <- params[rep(params[, .I], rep.sim)]


################################################################################
########################## Simulation function #################################
################################################################################

#studysim is the simulation function
#zeta=num1 #this is to debug simulation, each value of zeta is a unique run of the simulation

studysim <- function(zeta){
  
  ### set lambda to 50 for true population abundance before the intervention/treatment
  ### we assume this is the same for all sites beforehand
  ### this is completely arbitrary
  startvalt = 50
  
  ### set values for change in treatment, control, and noise
  ### we get these from the dataset params.rep that we made earlier in data preparation
  ### it's indexed by the value zeta for each unique run
  cht1 = params.rep$cht[zeta]
  chc1 = params.rep$chc[zeta]
  noise1 = params.rep$noise[zeta]
  
  ### we set some possible numbers of treatment AND control sites sampled (i.e., 5 and 5 sites)
  ### we could add more values but keep to three to keep simulation 'small'
  sitesno = c(5,20,40) 
  
  ### create list to capture estimates for each design later on
  results = list()
  
  #set counter to index results (we will cycle through different numbers of sites)
  z=1
  
  ###here we assume the total number of time steps is 2 years before and after
  ###doing this allows for the time series design called the 'after' design 
  ###(simplest one taking difference between first and last time step after treatment/intervention)
  ntime=2
  
  ###here we generate some true mean values for the population abundance before the
  ###intervention/treatment and introduce some noise
  ###these true values for each time step are drawn from a truncated normal distribution
  ###abundance must be greater than 0. We could use poisson but then would not be able
  ###to control noise? But maybe poisson would be more appropriate here.
  
  before.means <- rtruncnorm(2, 0,Inf,startvalt, sd = startvalt*noise1)
  
  ###we then generate some true mean values for the population abundance after the
  ###intervention/treatment and introduce noise here too.
  ###this is different for control versus treatment sites.
  ###treatment sites change by cht1 (a proportional change or response ratio)
  ###control sites change by chc1 (a separate proportional change or response ratio)
  after.treatment.means <- rtruncnorm(2, 0,Inf,startvalt*cht1, sd = startvalt*cht1*noise1)
  after.control.means <- rtruncnorm(2, 0,Inf,startvalt*chc1, sd = startvalt*chc1*noise1)
  
  ###for every parameter value of sites (treatment and control), we loop through and generate
  ###some site data to sample from
  
  #debugging variables - st sets the index for the number of sites (either 1, 2 or 3 corresponding to 5, 20 or 40 sites)
  # - tt sets the number of time steps (either 1 or 2)
  #st=5
  #tt=1
  
  
  for(st in 1:length(sitesno)){ #loop through each number of sites
    
    s1=sitesno[st] #set s1 as the no. treatment and control sites to generate
    sitedata <- list() #create list to put sitedata in for analysis later
    
    for(tt in 1:ntime){ #loop through each timestep i.e., year 1 and year 2
      t1=tt #set t1 as time step to sample
      
      #### generate 1000 sites to sample from based on the tru mean abundance for a 
      #### given year
      #### we use poisson here given that is the most appropriate for abundance
      after.treatment<- rpois(1000, after.treatment.means[t1]) #sites after treatment in year t1 for treatment sites
      after.control<- rpois(1000, after.control.means[t1]) #sites after treatment in year t1 for control sites
      before.treatment<- rpois(1000, before.means[t1]) #treatment sites before treatment in year t1 
      before.control<-  before.treatment #control sites before treatment in year t1 (same pool as treatment sites)
      
      ### for randomised designs, sample the 1000 sites above using random draws
      ### number of random draws = number of sites
      ### do this for each type of site (treatment and control) 
      ### and for each time period (before and after)
      ### create a data.table to hold this so we can run glms later on easily
      rand.after.treatment = data.table(randomised="randomised",before.after="after",
                                     year=t1+ntime,control.treatment="treatment",
                                     data=sample(after.treatment,s1)) #After treatment
      
      rand.after.control = data.table(randomised="randomised",before.after="after",
                                      year=t1+ntime,control.treatment="control",
                                      data=sample(after.control,s1)) #After control
      
      rand.before.treatment = data.table(randomised="randomised",before.after="before",
                                      year=t1,control.treatment="treatment",
                                      data=sample(before.treatment,s1)) #before treatment
      
      rand.before.control = data.table(randomised="randomised",before.after="before",
                                       year=t1,control.treatment="control",
                                       data=sample(before.control,s1)) #before control
      
      ###now for non-randomised designs we do things differently
      ###first I alternate the direction of bias in which 
      ###we allocate sites to treatment versus control
      ###this is alternated each run so as not to introduce systematic bias into simulation
      ###do this based on whether zeta is odd or even 
      ### if zeta is even, allocate sites above average to treatment and below average to control
      ### if zeta is odd, do the opposite
      ### we could do this differently. In my previous paper, I generated a parameter
      ### that determined the difference between treatment and control sites
      ### before the intervention/treatment, so we could revert to that.
      ### this would have a similar effect of introducing site selection bias
      
      if((zeta %% 2) ==0){
        
        nonrand.after.treatment = data.table(randomised="non-randomised",before.after="after",
                                          year=t1+ntime,control.treatment="treatment",
                                          data=sample(after.treatment[after.treatment>after.treatment.means[t1]],s1)) #After treatment
        
        nonrand.after.control = data.table(randomised="non-randomised",before.after="after",
                                           year=t1+ntime,control.treatment="control",
                                           data=sample(after.control[after.control<after.control.means[t1]],s1)) #After control
        
        
        nonrand.before.treatment = data.table(randomised="non-randomised",before.after="before",
                                           year=t1,control.treatment="treatment",
                                           data=sample(before.treatment[before.treatment>before.means[t1]],s1))#Before treatment
        
        nonrand.before.control = data.table(randomised="non-randomised",before.after="before",
                                            year=t1,control.treatment="control",
                                            data=sample(before.control[before.control<before.means[t1]],s1)) #Before control
      } else{
        
        nonrand.after.treatment = data.table(randomised="non-randomised",before.after="after",
                                          year=t1+ntime,control.treatment="treatment",
                                          data=sample(after.treatment[after.treatment<after.treatment.means[t1]],s1)) #After treatment
        
        nonrand.after.control = data.table(randomised="non-randomised",before.after="after",
                                           year=t1+ntime,control.treatment="control",
                                           data=sample(after.control[after.control>after.control.means[t1]],s1)) #After control
        
        
        nonrand.before.treatment = data.table(randomised="non-randomised",before.after="before",
                                           year=t1,control.treatment="treatment",
                                           data=sample(before.treatment[before.treatment<before.means[t1]],s1))#Before treatment
        
        nonrand.before.control = data.table(randomised="non-randomised",before.after="before",
                                            year=t1,control.treatment="control",
                                            data=sample(before.control[before.control>before.means[t1]],s1)) #Before control
      }
      
      ####once we have generated all the sites we need, we can put them in a list ready
      ####for analysis using glms in the next stage
      sitedata[[t1]]<- rbind(rand.after.control,rand.before.control,rand.after.treatment,rand.before.treatment,
                             nonrand.after.control,nonrand.before.control,nonrand.after.treatment,nonrand.before.treatment)
    }
    
    ####bring all the site data together ready for analysis
    allsitedata <- do.call(rbind,sitedata)
    
    ####prepare allsitedata variables for glms 
    allsitedata$before.after <- factor(allsitedata$before.after,levels=c("before","after"))
    allsitedata$control.treatment <- factor(allsitedata$control.treatment,levels=c("control","treatment"))
    allsitedata$yearfactor<- factor(allsitedata$year,levels=c((2*ntime):1))
    allsitedata$year <- as.numeric(allsitedata$year)
    
    ####run poisson glms on the site data generated earlier
    #### each glm is set up to analyse the data with a different study design
    #### rct (randomsed controlled trial), rbaci (randomised DiD), baci (DiD),
    #### ci (control-impact), ba (before-after), or after (time series)
    rct.glm <- glm(data~control.treatment, data=allsitedata[randomised=="randomised"&before.after=="after"], family="poisson") ###RCT
    rbaci.glm <- glm(data~before.after*control.treatment, data=allsitedata[randomised=="randomised"], family="poisson") ###RBACI
    baci.glm <- glm(data~before.after*control.treatment, data=allsitedata[randomised=="non-randomised"], family="poisson") ###BACI
    ci.glm <- glm(data~control.treatment, data=allsitedata[randomised=="non-randomised"&before.after=="after"], family="poisson") ###CI
    ba.glm <- glm(data~before.after, data=allsitedata[randomised=="non-randomised"&control.treatment=="treatment"], family="poisson") ###BA
    after.glm <- glm(data~yearfactor, data=allsitedata[randomised=="non-randomised"&before.after=="after"&control.treatment=="treatment"], family="poisson") ###AFTER
    
    #### extract coefficients we need from the glms
    rct.glm.coef <- cbind(est=coef(summary(rct.glm))[2,1],err=coef(summary(rct.glm))[2,2],design="rct")
    rbaci.glm.coef  <- cbind(est=coef(summary(rbaci.glm))[4,1],err=coef(summary(rbaci.glm))[4,2],design="rbaci")
    baci.glm.coef  <- cbind(est=coef(summary(baci.glm))[4,1],err=coef(summary(baci.glm))[4,2],design="baci")
    ci.glm.coef  <- cbind(est=coef(summary(ci.glm))[2,1],err=coef(summary(ci.glm))[2,2],design="ci")
    ba.glm.coef  <- cbind(est=coef(summary(ba.glm))[2,1],err=coef(summary(ba.glm))[2,2],design="ba")
    after.glm.coef <- cbind(est=coef(summary(after.glm))[2,1],err=coef(summary(after.glm))[2,2],design="after")
    
    #### put these estimates for each design into a master list that is outputted at the 
    #### end of the simulation.  
    results[[z]] <- cbind(rbind(rct.glm.coef,rbaci.glm.coef,baci.glm.coef,ci.glm.coef,ba.glm.coef,after.glm.coef),sites=s1,tru=params.rep$tru[zeta],noise=params.rep$noise[zeta]) ###number of treatment and control sites sampled and true effect size
    
    ### Update the index for the next loop - remember we need to repeat this for a 
    ### different number of sites
    z=z+1
  }
  
  ###package up results ready for output
  return(do.call(rbind,results))
}

####debugchecks if required
#num1<-sample(1:nrow(params.rep),1)
#studysim(num1)

##########################################################################
############################## Windows parallelisation ###################
##########################################################################

####set number of cores for your machine
####mine has 8 but I find using only 7 is necessary for performance
cl <- makeCluster(7)

####notify cluster of which libraries it needs to load for simulation
clusterEvalQ(cl,c(library(data.table),library(lme4),library(truncnorm)))

####export the data.tables that the simulation needs to run to the cluster
clusterExport(cl,c('params.rep'))

#### run the simulation for every row of params.rep (each unique run)
system.time({output1 <- parLapplyLB(cl,1:nrow(params.rep),studysim)})

#### make sure we end cluster
stopCluster(cl)

#### pull data together after outputting from simulation
output <- data.table(do.call(rbind,output1))
str(output)

#### save the output so we can reload again later
save(output, file = "studysim_output.rda")

#### load it back in if needed
#studysimfilepath <- "YOUR-FILEPATH-TO-studysim_output.rda"
#load(studysimfilepath)

##########################################################################
############################## Check simulation results ##################
##########################################################################

#### set variables 
output$est<-as.numeric(output$est)
output$err<-as.numeric(output$err)
output$sites<-as.numeric(output$sites)
output$tru<-as.numeric(output$tru)
output$deltatru <- abs(output$tru - output$est) #this is the absolute difference between
                                                #true effect and each estimate
                                                #i.e., a measure of bias
output$ID <- as.factor(1:nrow(output)) #put in an index

#### below we can plot the data just to check we are getting the results we'd expect
#### from the simulation

#### set parameter values for numbers of sites and types of design
site.options <- unique(output$sites)
design.options <- unique(output$design)

### get data in a form to plot it via ggplot2
output.results.sum <- data.table(output %>% 
                                   group_by(design,sites)%>%
                                   summarise(deltatru.abs.mean=
                                               list(mean_cl_normal(deltatru)%>%
                                                      rename(deltatru.abs.mean.y=y,
                                                             deltatru.abs.mean.ymin=ymin,
                                                             deltatru.abs.mean.ymax=ymax)),
                                             err=list(mean_cl_normal(err)%>%
                                                        rename(err.y=y,
                                                               err.ymin=ymin,
                                                               err.ymax=ymax)))%>% 
                                   unnest(cols=c(deltatru.abs.mean,err)))

#check it's ready for plotting
str(output.results.sum)

#set number of sites as a factor
output.results.sum$sites <- as.factor(output.results.sum$sites)

##plot in ggplot2
ggplot(aes(x=sites,y=deltatru.abs.mean.y),data=output.results.sum)+
  geom_line(aes(group=sites))+
  geom_pointrange(aes(ymin=deltatru.abs.mean.ymin,ymax=deltatru.abs.mean.ymax,fill=sites),shape=21)+
  facet_wrap(~design,nrow=1)+
  xlab("Number of sites")+
  ylab("True - estimated effect size")+
  scale_fill_manual(values=brewer.pal(5,"Oranges"))


##########################################################################
######## Data preparation for Meta-analysis simulation ###################
##########################################################################

#### set a parameter dataset with all combinations of numbers of sites,
#### study designs, and numbers of studies

#how many studies we add to the meta-analysis
num.marg.studies.options <- c(1,5,10)
#the design of studies already in the meta-analysis
design.options <- unique(output$design)
#the number of studies already in the meta-analysis
num.studies.options <- c(2,5,10) 
#get the data from the simulation on the different runs
#(each one has an associated pool of studies remember (1000 reps = studies))
trueff.noise.des.sites.runs <- unique(output[,list(tru,noise,design,sites)])

#find all combinations of these parameters and put into a dataset ready to load into metasim
params.meta.marg.dat <- data.table(crossing(trueff.noise.des.sites.runs,
                                            num.marg.studies=num.marg.studies.options,
                                            marg.design = design.options,
                                            num.studies=num.studies.options))


##########################################################################
############################# Meta-analysis simulation ###################
##########################################################################

#### function to run meta-analyses for different combinations of studies and designs

metasim.marg <- function(subset.dat, marg.study.design, num.studies, num.marg.studies){

  #function is given the estimates for all studies for that run of simulation
  #also information on the marginal study design to add
  #number of studies to add
  #and number of studies that should already be in the meta-analysis
  
  #first check if the design to add to meta-analysis is the same as one in original meta-analysis
  if(length(unique(subset.dat$design))==1){ #make sure to select right studies in initial MA
    subset.design <- subset.dat[design==marg.study.design,]
  } else{
    subset.design <- subset.dat[design!=marg.study.design,]
  }
  
  #now sample some studies of given design to put in original meta-analysis
  studies.sample <- sample(1:nrow(subset.design),num.studies)
  subset.design.studies <- subset.design[studies.sample,] 
  
  #run a simple random effects meta-analysis using inverse variance weighting
  meta.mod <- rma(subset.design.studies$est,subset.design.studies$err^2, data=subset.design.studies) #simple fixed effects currently, should it be random?
  
  
  #now let's rerun meta-analysis but this time adding in extra studies of a given design
  subset.design.marg <- subset.dat[design==marg.study.design,]#now select studies to add of particular design and rerun MA 
  studies.sample.marg <- sample(1:nrow(subset.design.marg),num.marg.studies)
  subset.design.studies.marg <- rbind(subset.design.studies,subset.design.marg[studies.sample.marg,])
  
  #now rerun meta-analysis with extra studies in of given design (meta-analysis unchanged otherwise)
  meta.mod.marg <- rma(subset.design.studies.marg$est,subset.design.studies.marg$err^2, data=subset.design.studies.marg) #simple fixed effects currently, should it be random?
  
  #package up original and new meta-analysis results
  return(data.table(rbind(coef(summary(meta.mod))[,c(1,2,4)],coef(summary(meta.mod.marg))[,c(1,2,4)]),marginal=c(0,1)))
}

#for debugging
#dat <- params.meta.marg.dat[7,]

#use another function to deliver the data that metasim.marg needs to run meta-analyses
#it takes a row of a dataset containing the information on which unique run of the simulation
#to sample study design estimates from
runmetasim.marg <- function(dat){
  #here we take the data we need for the run of the simulation we're interested in
  output.subset <- output[tru==dat$tru & noise==dat$noise &
                            (design==dat$design | design==dat$marg.design)&
                            sites==dat$sites,]
  
  #make a list to take the results we'll generate
  metasim.output.marg <- list()
  #set an index
  j=1
  #loop through this 100 times (could go up to 1000 but trying to keep simulation 'small')
  for(j in 1:100){
    #put in a tryCatch in case of errors.
    #run metasim.marg function with data we give it
    #marginal refers to whether it is the original meta-analysis or the one with added studies
    metasim.output.marg[[j]] <- tryCatch(cbind(metasim.marg(output.subset,
                                                   marg.study.design=dat$marg.design,
                                                   num.studies = dat$num.studies,
                                                   num.marg.studies = dat$num.marg.studies)),
                                         error = function(e){data.table(estimate=c(NA,NA),
                                                                        se=c(NA,NA),
                                                                        zval=c(NA,NA),
                                                                        pval=c(NA,NA),
                                                                        ci.lb=c(NA,NA),
                                                                        ci.ub=c(NA,NA),
                                                                        marginal=c(0,1))})
  }
  #package up data that we get back
  return(data.table(do.call(rbind,metasim.output.marg),dat))
}

#debugging and testing time it takes to run
#system.time({runmetasim.marg(params.meta.marg.dat[77777,])})


##########################################################################
############################## Windows parallelisation ###################
##########################################################################
#same as before - make a cluster, give it the libraries and data it needs
cl <- makeCluster(7)
clusterEvalQ(cl,c(library(data.table),library(lme4),library(truncnorm),library(metafor)))
clusterExport(cl,c('params.meta.marg.dat','metasim.marg','runmetasim.marg','output'))

#run metasim.results.marg for all parameter combinations
#add tryCatch in case of errors
metasim.results.marg <- parLapplyLB(cl,1:nrow(params.meta.marg.dat),function(x){
  tryCatch(runmetasim.marg(params.meta.marg.dat[x,]),
            error=function(e){data.table(estimate=c(NA,NA),
                                         se=c(NA,NA),
                                         zval=c(NA,NA),
                                         pval=c(NA,NA),
                                         ci.lb=c(NA,NA),
                                         ci.ub=c(NA,NA),
                                         marginal=c(0,1),params.meta.marg.dat)})
           }
)

#remember to stop cluster
stopCluster(cl)

#package up the meta-analysis results
metasim.results.all.marg <- do.call(function(x){rbind(x,fill=TRUE)},metasim.results.marg) 

#save these data
save(metasim.results.all.marg,metasim.results.marg, file = "metasim.rda")

#### load it back in if needed
#metasimsimfilepath <- "YOUR-FILEPATH-TO-metasimsim.rda"
#load(metasimsimfilepath)


##########################################################################
############################## Check meta-simulation results #############
##########################################################################

#set variables required for plotting
metasim.results.all.marg$deltatru.abs <- abs(metasim.results.all.marg$tru - metasim.results.all.marg$estimate)
metasim.results.all.marg$study.design <- factor(toupper(metasim.results.all.marg$design),levels=c("RCT","RBACI","BACI","CI","BA","AFTER"))
metasim.results.all.marg$marg.design <- factor(paste0("Add ",toupper(metasim.results.all.marg$marg.design)),levels=c("Add RCT","Add RBACI","Add BACI","Add CI","Add BA","Add AFTER"))

#calculate marginal changes in difference between true effect and meta-analytic estimate from adding studies of different designs
marg.diff <- data.table(metasim.results.all.marg %>%
                           group_by(study.design,sites, num.studies,marg.design,num.marg.studies)%>%
                           mutate(deltatru.abs.marg = c(0,diff(deltatru.abs)),
                                  se.marg = c(0,diff(se))))
marg.diff <- marg.diff[marginal==1] #only keep differences between marginal and original results, i.e., rows with marginal ==1

#save these data
save(marg.diff, file = "marg.diff.rda")
#### load it back in if needed
#metasimsimfilepath <- "YOUR-FILEPATH-TO-metasimsim.rda"
#load(metasimsimfilepath)

#summarise to find average change in absolute difference from adding a study across all runs and effect sizes
metasim.results.sum.marg <- data.table(marg.diff %>% 
                                    group_by(study.design,sites, num.studies,marg.design,num.marg.studies)%>%
                                    summarise(deltatru.abs.marg.mean=list(mean_cl_normal(deltatru.abs.marg)%>%
                                                                            rename(deltatru.abs.marg.mean.y=y,
                                                                                   deltatru.abs.marg.mean.ymin=ymin,
                                                                                   deltatru.abs.marg.mean.ymax=ymax)),
                                              se.marg=list(mean_cl_normal(se.marg)%>%
                                                          rename(se.marg.y=y,
                                                                 se.marg.ymin=ymin,
                                                                 se.marg.ymax=ymax)))%>%
                                      unnest(cols=c(deltatru.abs.marg.mean,se.marg)))

#set variables for plotting
metasim.results.sum.marg$sample.size <- as.factor(metasim.results.sum.marg$sites)
metasim.results.sum.marg$num.studies <- as.factor(metasim.results.sum.marg$num.studies)

#plot data on absolute difference 
ggplot(aes(x=num.marg.studies,y=deltatru.abs.marg.mean.y),data=metasim.results.sum.marg)+
  geom_line(aes(group=interaction(sample.size,num.studies),colour=num.studies))+
  geom_pointrange(aes(y=deltatru.abs.marg.mean.y, ymin = deltatru.abs.marg.mean.ymin,
                      ymax=deltatru.abs.marg.mean.ymax,fill=sample.size),shape=21)+
  geom_hline(aes(yintercept=0),linetype=2)+
  scale_x_continuous(breaks=c(1,5,10))+
  facet_grid(marg.design~study.design)+
  xlab("Number of studies added")+
  ylab("Change in absolute deviation from true effect")+
  scale_fill_manual(values=brewer.pal(5,"Oranges"))+
  scale_colour_manual(values=brewer.pal(5,"Oranges")[c(2,4,5)])+
  theme_cowplot()

### plot data on standard error
ggplot(aes(x=num.marg.studies,y=err.marg.y),data=metasim.results.sum.marg)+
  geom_line(aes(group=interaction(sample.size,num.studies),colour=num.studies))+
  geom_pointrange(aes(y=err.marg.y, ymin = err.marg.ymin,
                      ymax=err.marg.ymax,fill=sample.size),shape=21)+
  geom_hline(aes(yintercept=0),linetype=2)+
  scale_x_continuous(breaks=c(1,5,10))+
  facet_grid(marg.design~study.design)+
  xlab("Number of studies added")+
  ylab("Change in standard error")+
  scale_fill_manual(values=brewer.pal(5,"Oranges"))+
  scale_colour_manual(values=brewer.pal(5,"Oranges")[c(2,4,5)])+
  theme_cowplot()


### another possible measure in terms of statistical 
### significance of estimate - not sure this is right!
### sign and p-value errors

#find whether there was a p-value of <0.05 and the estimate's sign was 
#in the right direction/right sign
marg.diff.errrate <- data.table(metasim.results.all.marg %>%
                          group_by(study.design,sites, num.studies,marg.design,num.marg.studies,marginal)%>%
                          mutate(err.rate=1-mean(as.numeric(pval<0.05 & sign(estimate)==sign(trueffsize))))%>%

marg.diff.errrate <- marg.diff.errrate %>%
                        mutate(diff.errrate = c(0,diff(err.rate))))

marg.diff.errrate <- marg.diff.errrate[marginal==1] #only keep differences between marginal and original results, i.e., rows with marginal ==1

### summarise this data for plotting
metasim.results.sum.marg.errrate <- data.table(marg.diff.errrate %>% 
                                         group_by(study.design,sites, num.studies,marg.design,num.marg.studies)%>%
                                         summarise(err.rate=list(mean_cl_normal(diff.errrate)%>%
                                             rename(err.rate.y=y,
                                                    err.rate.ymin=ymin,
                                                    err.rate.ymax=ymax)
                                         ))%>%
                                         unnest(cols=err.rate))

#set variables for plotting
metasim.results.sum.marg.errrate$sample.size <- as.factor(metasim.results.sum.marg.errrate$sites)
metasim.results.sum.marg.errrate$num.studies <- as.factor(metasim.results.sum.marg.errrate$num.studies)

#plot data on statistical significance and sign
ggplot(aes(x=num.marg.studies,y=err.rate.y),data=metasim.results.sum.marg.errrate)+
  geom_line(aes(group=interaction(sample.size,num.studies),colour=num.studies))+
  geom_pointrange(aes(y=err.rate.y, ymin = err.rate.ymin,
                      ymax=err.rate.ymax,fill=sample.size),shape=21)+
  geom_hline(aes(yintercept=0),linetype=2)+
  scale_x_continuous(breaks=c(1,5,10))+
  facet_grid(marg.design~study.design)+
  xlab("Number of studies added")+
  ylab("Change in error rate")+
  scale_fill_manual(values=brewer.pal(5,"Oranges"))+
  scale_colour_manual(values=brewer.pal(5,"Oranges")[c(2,4,5)])+
  theme_cowplot()


# end of code