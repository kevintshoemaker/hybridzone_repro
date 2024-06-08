##################################################
# Wood Rat Reproduction                           #
# Reproductive success and reproductive 'choice' #
#K Shoemaker and E Hunter, 2021		               #
##################################################


#######################
# TODO

# KTS: remove strange male that skews the reproduction summaries: "S-659"
#   replace this individual in the parent dataset with  

# Incorporate environmental variables in rep prob and mate-pairing models? [NOPE]
# test to see if pure-pairings were more productive than out-cross pairings

#######################

rm(list=ls())

source("Repro_functions_v6.R")

##################################
##########  LOAD PACKAGES

library(lubridate)
# library(rgeos)
# library(rgdal)
library(raster)
library(sp)
library(edfun)
library(tidyverse)
library(EMT)
library(glmmTMB)

###################################
##########   GLOBAL VARS

Species <- c("Fuscipes","Fus_Back","Hybrid","Mac_Back","Macrotis")
Species2 <- c("N. fuscipes","BC-fuscipes","F1","BC-macrotis","N. macrotis")
Species3 <- c(expression(italic("N. fuscipes")),expression(italic("BC-fuscipes")),"F1",expression(italic("BC-macrotis")),expression(italic("N. macrotis")))
speccol <- c("dodgerblue4",rgb(70, 150, 170, maxColorValue=255),"palegreen3",rgb(210, 210, 60, maxColorValue=255),"gold3")
names(speccol) <- Species

pure <- c("Macrotis", "Fuscipes")

###################################
##########   READ IN DATA

######
# capture data

temp <- readcapdata()              # read in and process capture data
rat <- temp$rat
rat.clip <- temp$rat.clip
names(rat)


###########
# double check summary stats

length(sort(unique(rat$PJM.ID)))  # 2026 individuals captured between 2008 and 2013
sort(unique(rat$Year))

temp <- subset(rat,Year!=2013&!is.na(qFu))

length(sort(unique(temp$PJM.ID)))    # 1683 have associated genotypes

temp2 <- subset(temp, AgeFI=="A"&SexFI=="F")

length(sort(unique(temp2$PJM.ID)))    # 519 unique females

length(which(sapply(sort(unique(temp2$PJM.ID)),function(t) {t2 = subset(temp,PJM.ID==t); length(unique(t2$AgeFI)  )}  )>=2))

temp2 <- subset(temp, AgeFI=="A"&SexFI=="M")

length(sort(unique(temp2$PJM.ID)))    # 493 unique males

length(which(sapply(sort(unique(temp2$PJM.ID)),function(t) {t2 = subset(temp,PJM.ID==t); length(unique(t2$AgeFI)  )}  )>=2))

temp2 <- subset(temp, AgeFI=="J")

length(sort(unique(temp2$PJM.ID)))    # 719 unique juveniles

temp3 = sapply(sort(unique(rat$Year)), function(t){t2 = subset(temp,Year==t);length(sort(unique(t2$PJM.ID))) } )  # 191
names(temp3) = sort(unique(rat$Year))

temp3
######
# parentage data

temp <- readparentagedata()                # read in and process offspring parentage data
parent <- temp$parent                                 # all offspring parentage data
parent_knowndad <- temp$parent_knowndad               # all offspring parentage data with known dads or moms
parent_knownmom <- temp$parent_knownmom   
parent_pairings <- temp$parent_pairings               # one unique record for each pairing each year
parent_bothknown <- temp$parent_bothknown

allyears <- sort(unique(parent$Year))
rm(temp)


##### summary stats

nrow(parent_pairings)     # 579 unique pairings
nrow(parent)   # 823 offspring
nrow(parent_bothknown)  # 782 offspring with both parents known

nrow(parent_knownmom)

mean(parent_pairings$MateDist,na.rm=T)   # mean pairing dist of 62.56 m

temp2 <- subset(parent_pairings,Mom.spec=="Macrotis"&Dad.spec=="Macrotis")
nrow(temp2)  # 180 macrotis*macrotis

temp2 <- subset(parent_pairings,Mom.spec=="Fuscipes"&Dad.spec=="Fuscipes")
nrow(temp2)  # 296 macrotis*macrotis


###########################
#  EXPLORATORY VISUALS AND DESCRIPTIVE STATS
###########################

########################
### descriptive: size difference by species/genotype

# rat %>% group_by(SexFI,Gen2) %>% summarize(meanwt = mean(WtMax_g,na.rm=T))
plotsizedif()     # plot size difference by species for both sexes

make_summarytable1()

plotmatedist()     # plot mate distances (and write to svg)
 
reprostats <- doReproStats()   # summarize reproductive fraction, offspring production etc by year and genotype

reprostats$meanOPF

reprostats$meanOPF_gt

############
# GENERATE NULL-MODEL PAIRINGS
############

# runnullmods <- T
runnullmods <- F

if(runnullmods){       # NOTE: this can take a long time!
  exp.out.f <- generateNullPairings(sex="F",nsims=50)
  exp.out.m <- generateNullPairings(sex="M",nsims=50)
  
  write.csv(exp.out.f,"NullPairings_F.csv",row.names = F)
  write.csv(exp.out.m,"NullPairings_M.csv",row.names = F)
}else{
  exp.out.f <-read.csv("NullPairings_F.csv")
  exp.out.m <-read.csv("NullPairings_M.csv")
}


############
# STATISTICAL TESTS AND VISUALIZATIONS: 
############

# makeELIZboxplots("F","all",TRUE)
# makeELIZboxplots("F","all",FALSE)

                   # random pairings generated for each female
mergeddf_f <- makemergeddf(parent=parent_pairings,zone="all",sex="F")   # make merged dataset for running tests
nrow(mergeddf_f)

                   # random pairings generated for each male
mergeddf_m <- makemergeddf(parent=parent_pairings,zone="all",sex="M")   # make merged dataset for running tests
nrow(mergeddf_m)

                   # random pairings generated for both females and males together
mergeddf_all <- makemergeddf(parent=parent_pairings,zone="all",sex="all")   # make merged dataset for running tests
nrow(mergeddf_all)

nulldists_f <- makenulldists(mergeddf_f,sex='F')
nulldists_m <- makenulldists(mergeddf_m,sex='M')
nulldists_all <- makenulldists(mergeddf_all,sex='all')

#########
# apparent mate choice  [test WHO the different genotypes actually paired with and see if it's different 
#                 from who whey would have been expected to pair with under the null model]

testmateselection(sex="F",zone="all")    # for all known reproductive females, see if mate choice followed expectations by genotype

testmateselection(sex="M",zone="all")    # for all known reproductive males, see if mate choice followed expectations by genotype

##############
# test if some sex/genotype pairing types produce more or fewer offspring than expected...

pairings <- makepairingsdf()

######
# test: are some pairings more or less frequently observed than expected [Answer: YES]

chisq.test(pairings$offspring.obs,p=pairings$exp.freq,simulate.p.value = T)
test <- multinomial.test(pairings$offspring.obs,p=pairings$exp.freq,MonteCarlo = T, useChisq = T, ntrial=100000)
test$p.value


chisq.test(pairings$pairs.obs,p=pairings$exp.freq,simulate.p.value = T)
test <- multinomial.test(pairings$pairs.obs,p=pairings$exp.freq,MonteCarlo = T, useChisq = T, ntrial=100000)
test$p.value      ## Strong evidence that some pairings result in more offspring than other pairings

#####
# test: are successful pairings are more likely to have at least one pure parent 
#               than expected from the null model?

test_obs_vs_exp_purity(var="pureparent",mergeddf=mergeddf_all,nulldists=nulldists_all,hyp="greater",lab="frequency at least one pure parent")

#####
# test: are successful pairings more likely to have both pure parents?

test_obs_vs_exp_purity("bothpure",mergeddf_all,nulldists_all,"greater","frequency both pure parents")

#####
# test: are successful pairings less likely to have both pure parents of opposite types

test_obs_vs_exp_purity("bothpurex",mergeddf_all,nulldists_all,"less","frequency both pure parents opposite species")

#####
# test: are successful pairings more likely to have both pure parents of the same types?

test_obs_vs_exp_purity("bothpures",mergeddf_all,nulldists_all,"greater","frequency both pure parents same species")

#####
# test: are successful pairings less likely to have both hybrid parents (note-hybrid here refers to backcrosses as well)

test_obs_vs_exp_purity(var="bothhyb",mergeddf=mergeddf_all,nulldists=nulldists_all,hyp="less",lab="frequency both hybrid parents")


####
# test: are successful pairings with hybrid fathers less likely?

test_obs_vs_exp_purity("hybdad",mergeddf_f,nulldists_f,"less","frequency hybrid father")

####
# test: are successful pairings with hybrid mothers less likely?

test_obs_vs_exp_purity("hybmom",mergeddf_m,nulldists_m,"less","frequency hybrid mother")

mean(abs(parent_pairings$diff),na.rm=T)   # mean of 0.087 difference in q val

temp <- subset(mergeddf_f,truemat==0)
mean(abs(temp$diff),na.rm=T)   # mean of 0.087 difference in q val

######
# visualize some of the above tests 

make_purepar_all()     # save svg file

##############
## test if hybrid males and females are more likely to have pure partner in successful pairings

make_hybrid_puremat()   # save svg file

########
# test if cross-pairings in the direction of macrotis female with fuscipes male
#   more rare than cross-pairings of fuscipes female with macrotis male

testcrosspairings(sex="F")   # test: are female fuscipes more likely than female macrotis to cross-breed?
testcrosspairings(sex="M")   # test: are male fuscipes less likely than male macrotis to cross-breed?
testcrosspairings(sex="all")   

###############
#  Model prob of reproduction
###############

############
# Make master dataframe of reproductive success by individual/year
############

fuavl <- make_avlvar("Fuscipes")
maavl <- make_avlvar("Macrotis")

repdf <- make_repdf()
names(repdf)
repdf_F <- subset(repdf,sex=="F")
repdf_M <- subset(repdf,sex=="M")

plot_specbysize()    # make SVG file

plot_avlbytime()    # make SVG file

# Marjorie: is there a way to separate the species-sex combos to see who had the greatest loss over the course of the study?

rep_onlypure <- subset(repdf,species%in%pure)
table(rep_onlypure$didrep,rep_onlypure$year)


############################
# Model factors influencing probability of reproduction
############################

model_rep_prob()   # run candidate models and save to workspace

####### Visualize probability of reproduction

plot_probrepro()


###############
#  Model mate 'selection' (what determines 
#            who a female ends up successfully mating with given she is reproductive)
##############


#########
# pick a sex and species to model

sex <- "F"
thisspec <- Species[1]      #  "Fuscipes" "Fus_Back" "Hybrid"   "Mac_Back" "Macrotis"

#######################
# prepare a dataset of known pairings and associated random pairings 
#              for the conditional logistic regression models 

df_FUfem <- preparedata_forclr("F", "Fuscipes")
nrow(df_FUfem)

df_MAfem <- preparedata_forclr("F", "Macrotis")
nrow(df_MAfem)

df_FUmale <- preparedata_forclr("M", "Fuscipes")
nrow(df_FUmale)

df_MAmale <- preparedata_forclr("M", "Macrotis")
nrow(df_MAmale)

########################
# Fit models

bestmod_FUfem <- fit_mateselection(df_FUfem)
bestmod_MAfem <- fit_mateselection(df_MAfem)
bestmod_FUmale <- fit_mateselection(df_FUmale)
bestmod_MAmale <- fit_mateselection(df_MAmale)


#########
# Visualize relationships

visualize_mateselection(bestmod_FUfem,df_FUfem,"Fuscipes","F")

visualize_mateselection(bestmod_MAfem,df_MAfem,"Macrotis","F")

visualize_mateselection(bestmod_FUmale,df_FUmale,"Fuscipes","M")

visualize_mateselection(bestmod_MAmale,df_MAmale,"Macrotis","M")











#########
# START HERE
#########


# hist(sub0$dadsize,freq=F,breaks=8,xlim=c(100,600))
# hist(sub1$dadsize,add=T,col=rgb(0,0.2,0.9,0.5),freq=F,breaks=8)    # macrotis selects larger males to breed with, fuscipes not so much
# 
# mean(sub0$dadsize,na.rm=T)
# mean(sub1$dadsize,na.rm=T)
# 
# 
# allfudad <- unique(df2$FatherID[df2$Dad.spec=="Fuscipes"])
# allmacdad <- unique(df2$FatherID[df2$Dad.spec=="Macrotis"])
# 
# intersect(allfudad,allmacdad)
# setdiff(allfudad,allmacdad)   # debug- okay now
# 
# subset(df2,FatherID=="S-515")
# 
# names(df2)
# df3 <- subset(df2,FatherID%in%allfudad&truemat==0)    
# mean(df3$dadsize,na.rm=T)
# df3 <- subset(df2,FatherID%in%allmacdad&truemat==0)
# mean(df3$dadsize,na.rm=T)
# 
# 
# 
# 
# 
# 
# 
# ################
# # TRY GLOBAL MODELS
# ################
# 
# df3 <- df2
# 
# # Prepare covariates
# df3$qdiff_s <- scale(df3$diff)
# df3$sizediff_s <- scale(df3$sizediff)
# df3$stratum <- df3$mid_year
# df3$matepair <- as.factor(df3$truemat)
# df3$qFu.dad_s <- scale(df3$qFu.dad) 
# df3$dadsize_s <- scale(df3$dadsize)
# df3$hzs_s <- scale(df3$hzs)
# df3$hmf_s <- scale(df3$hzs)
# df3$fsurv_s <- scale(df3$fsurv)
# 
# Species2 <- c("Fuscipes","Hybrid","Macrotis")
# pure <- c("Fuscipes","Macrotis")
# 
# df3$Mom.spec2 <- "Hybrid"
# df3$Mom.spec2[df3$Mom.spec=="Macrotis"] <- "Macrotis"
# df3$Mom.spec2[df3$Mom.spec=="Fuscipes"] <- "Fuscipes"
# df3$puremom <- ifelse(df3$Mom.spec%in%pure,1,0)
# df3$same <- ifelse((df3$Mom.spec==df3$Dad.spec),1,0 )
# df3$ismac <- ifelse(df3$Mom.spec=="Macrotis",1,0)
# df3$ishyb <- ifelse(df3$Mom.spec2=="Hybrid",1,0)
# 
# df3$Avail_fu <- as.numeric(df3$Avail_fu)
# df3$Avail_hy <- as.numeric(df3$Avail_hy)
# df3$Avail_ma <- as.numeric(df3$Avail_ma)
# 
# ######
# # Remove any NA observations
# 
# ndx <- complete.cases(df3[,c("qdiff_s","sizediff_s","stratum","qFu.dad_s","dadsize_s")])
# 
# df3 <- df3[ndx,]
# nrow(df3)         #27154
# 
# df3$weights <- ifelse(df3$truemat==1,1,1000)
# 
# 
# 
# 
# mod1f = glmmTMB(matepair ~ qFu.dad_s + qFu.dad_s:ismac + qFu.dad_s:ishyb +
#                   dadsize_s + dadsize_s:ismac + dadsize_s:ishyb +        #  
#                     (1|stratum),
#                  family=binomial,
#                  weights=weights,
#                  data=df3
# )
# summary(mod1f)
# AIC(mod1f)
# 
# 
# names(df3)
# library(survival)
# df3$matepair2 <- as.numeric(df3$matepair)
# mod1f = survival::clogit(matepair2 ~ qFu.dad_s + qFu.dad_s:ismac + qFu.dad_s:ishyb +  # MateDist + 
#                   sizediff_s + sizediff_s:ismac + sizediff_s:ishyb + 
#                   strata(stratum),        #  
#                 data=df3
# )
# summary(mod1f)
# AIC(mod1f)
# 
# 
# mod1f = survival::clogit(matepair2 ~ qFu.dad_s + qFu.dad_s:ismac + qFu.dad_s:ishyb + MateDist + 
#                            sizediff_s + poly(sizediff_s,2) +
#                            sizediff_s:ismac + sizediff_s:ishyb + 
#                            poly(sizediff_s,2):ismac + poly(sizediff_s,2):ishyb +
#                            strata(stratum),        #  
#                          data=df3
# )
# summary(mod1f)
# AIC(mod1f)
# 
# 
# mod1f = survival::clogit(matepair2 ~ qFu.dad_s + qFu.dad_s:ismac + qFu.dad_s:ishyb + MateDist + 
#                            poly(sizediff_s,2)[,1] +poly(sizediff_s,2)[,2] +
#                            poly(sizediff_s,2)[,1]:ismac +  poly(sizediff_s,2)[,1]:ishyb +
#                            poly(sizediff_s,2)[,2]:ismac +  poly(sizediff_s,2)[,2]:ishyb +
#                            strata(stratum),        #  
#                          data=df3
# )
# summary(mod1f)
# AIC(mod1f)
# 
# mod1g = survival::clogit(matepair2 ~ qFu.dad_s + qFu.dad_s:ismac + qFu.dad_s:ishyb + MateDist + poly(sizediff_s,2) +
#                            poly(sizediff_s,2) : Mom.spec2  +
#                            strata(stratum),        #  
#                          data=df3
# )
# summary(mod1g)
# AIC(mod1g)
# 
# names(df3)
#   cor(df3$Avail_fu,df3$Avail_ma,use="complete.obs")
# #### TRY AVAILABILITY
# mod1x = glmmTMB(matepair ~ qFu.dad_s + qFu.dad_s:ismac + qFu.dad_s:ishyb +
#                   dadsize_s + dadsize_s:ismac + dadsize_s:ishyb +        #
#                   qFu.dad_s:Avail_fu + qFu.dad_s:Avail_fu:ismac + qFu.dad_s:Avail_fu:ishyb +
#                   (1|stratum),
#                 family=binomial,
#                 weights=weights,
#                 data=df3
# )
# summary(mod1x)
# AIC(mod1x)
# 
# library(survival)
# df3$matepair2 <- as.numeric(df3$matepair)
# mod1x = survival::clogit(matepair2 ~ qFu.dad_s + qFu.dad_s:ismac + qFu.dad_s:ishyb +  # MateDist + 
#                            sizediff_s + sizediff_s:ismac + sizediff_s:ishyb + 
#                            qFu.dad_s:Avail_fu + qFu.dad_s:Avail_fu:ismac + qFu.dad_s:Avail_fu:ishyb +
#                            strata(stratum),        #  
#                          data=df3
# )
# summary(mod1x)
# AIC(mod1x)
# 
# # try cubic splines?  splines::bs()
# # try gam?
# 
# 
# 
# 
# ################
# # Ancillary code (not totally relevant but possibly worth hanging on to)
# 
# ###########
# # QUESTION: does male fuscipes size change over time?  [answer: no]
# 
# sub1 <- subset(rat,SexFI=="M"&Gen2=="Fu"&AgeFI=="A")
# 
# plot(sub1$WtMax_g~sub1$Year)
# 
# sub2 <- sub1 %>% group_by(Year,PJM.ID) %>% summarize(meanwt = mean(WtMax_g,na.rm=T))
# 
# plot(sub2$meanwt~sub2$Year)
# 
# summary(lm(sub2$meanwt~sub2$Year))
# 
# 
# ## model linking mother and father q value to determine offspring q
# mod_offspringq <- lm(parent$qFu.offspring~parent$qFu.mom+parent$qFu.dad-1)
# summary(mod_offspringq)
# 












########
# NOTES

# model reproductive probability as a function of percent availability of pure parents?

# separate concepts: 
# (1) apparent mate selection (of those individuals we know reproduced, who did they mate with)
# (2a) probability of mating (of those captured individuals, what factors determine who successfully reproduced)
#    (2b) reproductive output (related to 2a: of those captured individuals, what factors determine how many offspring successfully produced?

# test if successful pairings are more likely to have at least one pure parent 
#   than would expected based on availability (global test)  [TRY IT] [DONE- PURE PARENT MORE LIKELY]
#   test if F1 and F2 are disproportionately mating with pure parent [TRY IT]  [DONE- YES THEY ARE ]

# are matings in the direction of macrotis female with fuscipes male
#   more rare than matings of fuscipes female with macrotis male [TRY IT] [DONE- NOT DETECTABLY DIFFERENT]





#########
# FUTURE TODO

# TODO: incorporate availability explicitly...  [NOT FOR FIRST PAPER] 
# Build conditional logistic regression 'step selection' model of mate pairings
# the linear predictor represents a weighting scheme.
# use the linear predictor to weight probabilities of mate selection after accounting for distance. 
# NOTE: forester 2009 and avgar 2015 suggests including inter-pair distance (step length) as a covariate to unbias the fact that empirical distances already reflect selection patterns. 

# try using amt package (animal movement toolbox) to run the conditional logistic regression models...? 

