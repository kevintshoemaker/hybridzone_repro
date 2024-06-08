#################################
# Functions to support woodrat reproduction paper
#################################


#################################
# K Shoemaker and E Hunter 2021
#################################



#################################
# OVERVIEW OF SIMULATION SCRIPT
#################################

# First, read in master data and offspring/parent data / compute summary stats
# For every observed mated female, determine potential mate pairings based only on availability
#   STATISTICS: are females 'selecting' mates based on availability? [Repeat for males]
# For every observed adult female in the master dataset, determine 100+ potential mate pairings based on availability


# Combine with 'true' (observed) pairings to create a 'presence/background' dataset
#   STATISTICS: are 'real' pairings (vs. available pairings) more likely to have at least one pure parent?
#   STATISTICS: are matings between fuscipes female and macrotis male more rare than matings between macrotis female and fuscipes male?
#   STATISTICS: logistic regression: what is the probability of a successful pairing for macrotis females? For fuscipes females? Hybrids?
#   STATISTICS: global model of mate-pairing success?



##################
# FUNCTION: add figure labels (e.g., A B C labels to figure panels)

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}



##################
# FUNCTION: plot size by species for each sex

plotsizedif <- function(){
  svg("SizebySpec_MandF.svg",4,6)
  par(mfrow=c(2,1))
  thissex="M"
  sub1 <- subset(rat,SexFI==thissex)
  sub2 <- subset(sub1,Gen2=="Fu"&AgeFI=="A")
  hist(sub2$WtMax_g,ylim=c(0,0.02),main=paste0("size dif by species, ",thissex),
       col=adjustcolor(speccol["Fuscipes"],alpha.f = 0.4),freq=F,xlab="Size, g")
  cat(sprintf("Fuscipes %s: %#.3f\n",thissex,mean(sub2$WtMax_g,na.rm=T)))
  
  sub2 <- subset(sub1,Gen2=="Ma"&AgeFI=="A")
  hist(sub2$WtMax_g,add=T,
       col=adjustcolor(speccol["Macrotis"],alpha.f = 0.4),freq=F)
  cat(sprintf("Macrotis %s: %#.3f\n",thissex,mean(sub2$WtMax_g,na.rm=T)))
  #mean(sub2$WtMax_g,na.rm=T)
  
  legend("topleft",bty="n",fill=c(adjustcolor(speccol["Fuscipes"],alpha.f = 0.4),adjustcolor(speccol["Macrotis"],alpha.f = 0.4)),
         legend=c("Fuscipes","Macrotis"))
  
  thissex="F"
  sub1 <- subset(rat,SexFI==thissex)
  sub2 <- subset(sub1,Gen2=="Fu"&AgeFI=="A")
  hist(sub2$WtMax_g,ylim=c(0,0.02),main=paste0("size dif by species, ",thissex),
       col=adjustcolor(speccol["Fuscipes"],alpha.f = 0.4),freq=F,xlab="Size, g")
  cat(sprintf("Fuscipes %s: %#.3f\n",thissex,mean(sub2$WtMax_g,na.rm=T)))
  #mean(sub2$WtMax_g,na.rm=T)
  
  sub2 <- subset(sub1,Gen2=="Ma"&AgeFI=="A")
  hist(sub2$WtMax_g,add=T,
       col=adjustcolor(speccol["Macrotis"],alpha.f = 0.4),freq=F)
  cat(sprintf("Macrotis %s: %#.3f\n",thissex,mean(sub2$WtMax_g,na.rm=T)))
  #mean(sub2$WtMax_g,na.rm=T)
  
  legend("topleft",bty="n",fill=c(adjustcolor(speccol["Fuscipes"],alpha.f = 0.4),adjustcolor(speccol["Macrotis"],alpha.f = 0.4)),
         legend=c("Fuscipes","Macrotis"))
  
  dev.off()
}




####################
# FUNCTION: read and process woodrat capture data

#### read in capture data
####      this dataset stores information on all known individuals captured at camp roberts
####          id, genotype, sex, date, body size, etc.
####          each row is capture event, each column is individual attribute
####          genotype is coded as 'qFu', or the fraction of the genome that is 'fuscipes-like'

readcapdata <- function(){
  
  rat <- read.csv("Data\\CR_all_12dec2016.csv", header=TRUE, na.strings = c("NA",".",""), stringsAsFactors = FALSE)   # all rat captures
  nrow(rat)     # 13825 captures
  
  rat$DATE <- dmy(rat$Date1900)
  rat <- rat[order(rat$DATE),]
  rat$qFu <- as.numeric(as.character(rat$qFu))
  rat$Gen2 <- ifelse(rat$qFu>0.1 & rat$qFu<0.4, "bcMa", ifelse(rat$qFu>0.6 & rat$qFu<0.89, "bcFu", rat$Gen))
  
  table(rat$Gen2)	 
  
  names(rat)
  
  rat$id_year <- paste0(rat$PJM.ID,rat$Year)
  rat$Species <- ifelse(rat$qFu<=0.1,"Macrotis",ifelse(rat$qFu>=0.9,"Fuscipes",ifelse(rat$qFu<0.4,"Mac_Back",ifelse(rat$qFu>0.6,"Fus_Back","Hybrid"))))
  
  ndx <- !grepl("unk",rat$PJM.ID)
  rat <- rat[ndx,]    # remove unknown individs
  
  #Within the heavily sampled central area, they sampled every year, only did the extremes in 2008 and 2010
  #Within that central area, they were not sampling every area every month
  #But they were sampling every area within the central area in both early/late seasons
  
  #Need to remove captures in those extremity areas 
  rat.coords <- coordinates(rat[,c(38,37)])
  rat.spdf <- SpatialPointsDataFrame(coords=rat.coords, data=rat)
  #e <- extent(698500, 700500, 3960000, 3963000)
  e <- extent(rat.spdf)  #LEAVING IN EXTRIMITY AREAS FOR NOW
  rat.clip <- crop(rat.spdf, e)
  
  ndx <- complete.cases(rat.clip@data[,c("Year", "SexFI", "AgeFI")])
  rat.clip <- rat.clip[ndx,]
  # plot(rat.clip)
  
  #Add in column for zone
  rat.clip@data$Zone <- ifelse(rat.clip@data$UTM_N<3960975, "South", ifelse(rat.clip@data$UTM_N>3961300, "North", "Hybrid"))
  
  rat <- rat.clip@data
  
  return(list(rat=rat,rat.clip=rat.clip))
}


################
# FUNCTION: read in and process offspring parentage data
################

###########
# Read in parentage data
#    This data set stores the putative parents assigned to offspring captured in the study
#    Assignment was done with genetic assignment tests (see text for detail)


readparentagedata <- function(){
  
  #setwd(DataDir)
  parent <- read.csv("Data\\parent5.csv", header=TRUE, stringsAsFactors = FALSE,na.strings = c("NA","","-9999"))
  
  #In 2010, IDs that start with S need a hyphen after the S
  parent$OffspringID <- ifelse(substr(parent$OffspringID, 1, 1) == "S" & !grepl("-", parent$OffspringID), paste("S-", substr(parent$OffspringID, 2, nchar(parent$OffspringID)), sep=""), parent$OffspringID)
  parent$MotherID <- ifelse(substr(parent$MotherID, 1, 1) == "S" & !grepl("-", parent$MotherID), paste("S-", substr(parent$MotherID, 2, nchar(parent$MotherID)), sep=""), parent$MotherID)
  parent$FatherID <- ifelse(substr(parent$FatherID, 1, 1) == "S" & !grepl("-", parent$FatherID), paste("S-", substr(parent$FatherID, 2, nchar(parent$FatherID)), sep=""), parent$FatherID)
  
  # table(parent$MotherID)
  # table(parent$FatherID)
  
  #Remove female 729...not in full dataset
  parent <- parent[parent$MotherID!="729",]
  
  strangemale <- "S-659"
  parent[parent$FatherID%in%strangemale,]$FatherID <- NA
  
  
  # #################
  # # correct data errors
  # 
  to.rm = c()
  to.rm = c(to.rm,rownames(subset(parent, OffspringID=="461"))[2])
  to.rm = c(to.rm,rownames(subset(parent, OffspringID=="S-434"))[2])
  to.rm = c(to.rm,rownames(subset(parent, OffspringID=="S-822"))[2])
  to.rm = c(to.rm,rownames(subset(parent, OffspringID=="S-936"))[2])
  
  parent <- parent[setdiff(rownames(parent),to.rm),]
  
  parent <- parent[,!grepl("X",names(parent))]
  
  #parent[parent$MotherID=="3552",]
  
  ###########
  # Further develop "parent" dataset for comparison with simulated data
  
  #Get q-values in the "parent" df
  lu <- rat[,c("PJM.ID", "qFu","SexFI")]      # make lookup table
  lu <- lu[!duplicated(lu$PJM.ID),]
  nrow(lu)   # 2032 unique individuals in lookup table
  
  # subset(lu,PJM.ID=="S-515")   #test
  
  # names(parent)
  
  parent$qFu.dad <- lu$qFu[match(parent$FatherID,lu$PJM.ID)]
  parent$qFu.mom <- lu$qFu[match(parent$MotherID,lu$PJM.ID)]
  parent$qFu.offspring <- lu$qFu[match(parent$OffspringID,lu$PJM.ID)]
  parent$sex.offspring <- lu$SexFI[match(parent$OffspringID,lu$PJM.ID)]
  
  tmp <- subset(parent,FatherID=="S-515")   #test
  tmp
  
  parent$diff <-  parent$qFu.mom - parent$qFu.dad
  parent$sqdiff <- parent$diff^2
  
  parent$mid_year <- paste0(parent$MotherID,parent$Year)    # unique id-year combos
  parent$pid_year <- paste0(parent$FatherID,parent$Year)
  
  parent$matpairid <- paste0(parent$MotherID,"_",parent$FatherID,"_",parent$Year)
  
  parent <- parent[!is.na(parent$OffspringID),]
  
  nrow(parent) # 823 offspring with known parents
  length(unique(parent$OffspringID)) 
  
  length(unique(parent$MotherID))   # 321 unique mother ids
  length(unique(paste0(parent$MotherID,parent$Year)))   # 423 unique mother/year combos
  length(unique(parent$mid_year))
  
  parent$Mom.spec <- ifelse(parent$qFu.mom<=0.1,"Macrotis",ifelse(parent$qFu.mom>=0.9,"Fuscipes",ifelse(parent$qFu.mom<0.4,"Mac_Back",ifelse(parent$qFu.mom>0.6,"Fus_Back","Hybrid"))))
  parent$Dad.spec <- ifelse(parent$qFu.dad<=0.1,"Macrotis",ifelse(parent$qFu.dad>=0.9,"Fuscipes",ifelse(parent$qFu.dad<0.4,"Mac_Back",ifelse(parent$qFu.dad>0.6,"Fus_Back","Hybrid"))))
  parent$Offspring.spec <- ifelse(parent$qFu.offspring<=0.1,"Macrotis",ifelse(parent$qFu.offspring>=0.9,"Fuscipes",ifelse(parent$qFu.offspring<0.4,"Mac_Back",ifelse(parent$qFu.offspring>0.6,"Fus_Back","Hybrid"))))
  
  ######### Add spatial info to parent dataset
  
  #Offspring and mother location are the same, so look at distance b/w offspring and father locations
  #Use first location for each within a year as "true" location?
  
  parent$MotherUTM_N <- NA
  parent$MotherUTM_E <- NA
  parent$FatherUTM_N <- NA
  parent$FatherUTM_E <- NA
  parent$MateDist <- NA
  parent$MotherZone <- NA
  parent$FatherZone <- NA
  
  for(i in 1:nrow(parent)){
    temp.off <- rat[rat$PJM.ID==parent$MotherID[i] & rat$Year==parent$Year[i],]
    parent$MotherUTM_N[i] <- temp.off$UTM_N[1]
    parent$MotherUTM_E[i] <- temp.off$UTM_E[1]
    temp.dad <- rat[rat$PJM.ID==parent$FatherID[i] & rat$Year==parent$Year[i],]
    parent$FatherUTM_N[i] <- temp.dad$UTM_N[1]
    parent$FatherUTM_E[i] <- temp.dad$UTM_E[1]
    parent$MateDist[i] <- sqrt((parent$MotherUTM_E[i] - parent$FatherUTM_E[i])^2 + (parent$MotherUTM_N[i] - parent$MotherUTM_N[i])^2)
    parent$MotherZone[i] <- temp.off$Zone[1]
    parent$FatherZone[i] <- temp.dad$Zone[1]
  }
  
  
  #############
  # Make separate datasets for offspring with known mom and known dads (and both known)
  
  parent_knownmom <- parent[!is.na(parent$MotherID),]
  parent_knowndad <- parent[!is.na(parent$FatherID),]
  
  ndx <- complete.cases(cbind(parent$MotherID,parent$FatherID))
  parent_bothknown <- parent[ndx,]
  
  # names(parent)
  # tmp <- subset(parent,FatherID=="S-515")   # debug check
  # tmp
  
  
  ###########################
  # Make dataset with only unique pairings (exclude pairings that resulted in multiple detected offspring)
  #    This data set stores all unique parent-pairs assigned to offspring captured in the study each year
  #    Assignment was done with genetic assigment tests (see text for detail)
  ###########################
  
  #names(parent)
  #parent$pairid <- paste0(parent$mid_year,parent$pid_year)
  # length(unique(parent$matpairid))   # 613 unique pairings
  # nrow(parent)     # 823 offspring
  
  parent_pairings <- parent_bothknown[!duplicated(parent_bothknown$matpairid),]
  nrow(parent_pairings)
  
  
  return(
    list(parent=parent,
         parent_knowndad=parent_knowndad,
         parent_knownmom=parent_knownmom,
         parent_pairings=parent_pairings,
         parent_bothknown=parent_bothknown)
    )
}


#################
# FUNCTION: plot mating distances


plotmatedist <- function(){
  #HISTOGRAM OF MATING DISTANCES
  # - to create expectation, use the 95%ile as the distance threshold
  
  # setwd(FiguresDir)
  
  svg("MatingDists2.svg", width=6.5, height=4)
  par(mfrow=c(1,2))
  hist(parent_pairings$MateDist, breaks=100, xlab="Mating Distance", main="All distances")
  hist(parent_pairings$MateDist, breaks=100, xlim=c(0,500), xlab="Mating Distance", main="Zoomed in to distances <500m")
  dev.off()
  
  #graphics.off()
  #plot mating distances by area
  # setwd(FiguresDir)
  svg("MatingDistsLocs2.svg", width=6.5, height=7)
  plot(parent_pairings$MotherUTM_N ~ parent_pairings$MotherUTM_E, col=heat.colors(parent_pairings$MateDist))
  legend("topleft", legend=c("close", "far"), fill=c("yellow", "red"))
  dev.off()
  
}


#######################
# FUNCTION: STANDARD ERROR OF WEIGHTED MEAN
#######################

weighted.var.se <- function(x, w, na.rm=FALSE)
  #  Computes the variance of a weighted mean following Cochran 1977 definition
{
  if (na.rm) { w <- w[i <- !is.na(x)]; x <- x[i] }
  n = length(w)
  xWbar = weighted.mean(x,w,na.rm=na.rm)
  wbar = mean(w)
  out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  return(out)
}

####################
# FUNCTION: REPRODUCTION STATS (by year, genotype etc)
####################

# compute offspring per year and fraction reproductive for all genotypes and both sexes

# OPF = offspring per female
# OPM = offspring per male
# FRF = fraction reproductive, female
# FRM = fraction reproductive, male
# PPF = unique pairings per female
# PPM = unique pairings per male
# PPF2 = unique pairings per successful female
# PPM2 = unique pairings per successful male

# PPF2 = unique pairings per successful female
# PPM2 = unique pairings per successful male

doReproStats <- function(){
  
  reprostats <- list()
  
  ## determine the average number of known offspring per female per year etc.
  
  temp <- subset(rat,(SexFI=="F")&(AgeFI=="A"))
  allmidyear <- sort(unique(temp$id_year))
  reprostats$meanOPF <-  mean(sapply(allmidyear,function(t) length(which(parent$mid_year==t))))   # 0.86
  reprostats$meanPPF <- mean(sapply(allmidyear,function(t) length(which(parent_pairings$mid_year==t))))
  reprostats$meanFRF <- mean(sapply(allmidyear,function(t) ifelse(t%in%parent$mid_year,1,0) ))
  
  temp <- subset(rat,(SexFI=="M")&(AgeFI=="A"))
  allpidyear <- sort(unique(temp$id_year))
  reprostats$meanOPM <-  mean(sapply(allpidyear,function(t) length(which(parent$pid_year==t))))    #0.97
  reprostats$meanPPM <- mean(sapply(allpidyear,function(t) length(which(parent_pairings$pid_year==t))))
  reprostats$meanFRM <- mean(sapply(allpidyear,function(t) ifelse(t%in%parent$pid_year,1,0) ))
  
  Species
  
  allmidyear_known <- sort(unique(parent$mid_year))
  reprostats$meanOPF2 <- mean(sapply(allmidyear_known,function(t) length(which(parent$mid_year==t))))  # 1.95
  
  allpidyear_known <- sort(unique(parent$pid_year))
  reprostats$meanOPM2 <- mean(sapply(allpidyear_known,function(t) length(which(parent$pid_year==t))))  # 2.61
  
  reprostats$meanPPF2 <- mean(sapply(allmidyear_known,function(t) length(which(parent_pairings$mid_year==t))))  # 1.36
  reprostats$meanPPM2 <- mean(sapply(allpidyear_known,function(t) length(which(parent_pairings$pid_year==t))))  # 1.83
  
  #########
  # First compute annual stats
  # then summarize overall by genotype
  
  
  ############
  # determine the fraction of repeat offspring that are from a different male
  
  
  head(parent)
  allreps <- c()
  
  
  i=1
  for(i in 1:length(allmidyear_known)){
    thisindyr <- allmidyear_known[i] 
    temp <- subset(parent,mid_year==thisindyr)
    if(nrow(temp)>1){
      j=2
      for(j in 2:nrow(temp)){
        allreps <- c(allreps,ifelse(temp$FatherID[j]!=temp$FatherID[j-1],1,0))
      }
    }
  }
  
  mean(allreps,na.rm=T)    # 0.5275591 is the fraction of repeat offspring that are from a different male
  
  allreps <- c()
  i=1
  for(i in 1:length(allpidyear_known)){
    thisindyr <- allpidyear_known[i] 
    temp <- subset(parent,pid_year==thisindyr)
    if(nrow(temp)>1){
      j=2
      for(j in 2:nrow(temp)){
        allreps <- c(allreps,ifelse(temp$MotherID[j]!=temp$MotherID[j-1],1,0))
      }
    }
  }
  
  mean(allreps,na.rm=T)    # 0.6003937 is the fraction of repeat offspring that are from a different male
  
   
  
  allyears = sort(unique(rat$Year))
  
  meanpairingsf_year <- numeric(length(allyears)-1)
  meanpairingsm_year <- numeric(length(allyears)-1)
  meanoffspringf_year <- numeric(length(allyears)-1)
  meanoffspringm_year <- numeric(length(allyears)-1)
  meanoffspringf2_year <- numeric(length(allyears)-1)    # offspring conditional on reproducing 
  meanoffspringm2_year <- numeric(length(allyears)-1)
  fracreprodf_year <- numeric(length(allyears)-1)
  fracreprodm_year <- numeric(length(allyears)-1)
    
  check <- logical(length(allyears))
  
  meanpairingsf_gtyear <- matrix(NA,length(Species),(length(allyears))-1)
  meanpairingsm_gtyear <- matrix(NA,length(Species),(length(allyears))-1)
  meanoffspringf_gtyear <- matrix(NA,length(Species),(length(allyears))-1)
  meanoffspringm_gtyear <- matrix(NA,length(Species),(length(allyears))-1)
  meanoffspringf2_gtyear <- matrix(NA,length(Species),(length(allyears))-1)
  meanoffspringm2_gtyear <- matrix(NA,length(Species),(length(allyears))-1)
  fracreprodf_gtyear <- matrix(NA,length(Species),(length(allyears))-1)
  fracreprodm_gtyear <- matrix(NA,length(Species),(length(allyears))-1)
  
  totobs_eachyear_f <- numeric(length(allyears)-1)
  totobs_eachyear_m <- numeric(length(allyears)-1)
  
  alloffspringdataf <- c()
  alloffspringdatam <- c()
  
  y=1
  for(y in 1:(length(allyears)-1)){
    if(y==1){
      ys <- allyears[1]
    }else{
      ys <- allyears[y]   # try limiting to just the focal year # c(y-1,y,y+1)
    }
    temp <- subset(rat,(Year%in%ys)&(SexFI=="F")&(AgeFI=="A"))
    allpos <- sort(unique(temp$PJM.ID))
    temp2 <- subset(parent_pairings,Year==allyears[y])
    temp3 <- subset(parent,Year==allyears[y])
    allpos2 <- sort(unique(temp2$MotherID))
    
    check[y] <- all(temp2$MotherID %in% allpos)  # check to make sure all the known-parent females that year are in the master dataset for that year
    
    meanpairingsf_year[y] <- mean(sapply(allpos,function(t) length(which(temp2$MotherID==t))))
    meanoffspringf_year[y] <- mean(sapply(allpos,function(t) length(which(temp3$MotherID==t))))
    meanoffspringf2_year[y] <- mean(sapply(allpos2,function(t) length(which(temp3$MotherID==t))))
    fracreprodf_year[y] <- mean(sapply(allpos,function(t) ifelse(t%in%temp3$MotherID,1,0)  ))
    
    alloffspringdataf <- c(alloffspringdataf,sapply(allpos,function(t) length(which(temp3$MotherID==t))))
    
    meanpairingsf_gtyear[,y] <- sapply(Species,function(z) {
      sub1 <- subset(temp,Species==z)
      allpos <- sort(unique(sub1$PJM.ID))
      mean(sapply(allpos,function(t) length(which(temp2$MotherID==t))))
    }   )
    
    meanoffspringf_gtyear[,y] <- sapply(Species,function(z) {
      sub1 <- subset(temp,Species==z)
      allpos <- sort(unique(sub1$PJM.ID))
      mean(sapply(allpos,function(t) length(which(temp3$MotherID==t))))
    }   )
    
    meanoffspringf2_gtyear[,y] <- unlist(sapply(Species,function(z) {
      sub1 <- subset(temp3,Mom.spec==z)
      allpos <- sort(unique(sub1$MotherID))
      if(length(allpos)>0){ mean(sapply(allpos,function(t) length(which(temp3$MotherID==t)))) }else{NA}
    }   ))
    
    fracreprodf_gtyear[,y] <- sapply(Species,function(z) {
      sub1 <- subset(temp,Species==z)
      allpos <- sort(unique(sub1$PJM.ID))
      mean(sapply(allpos,function(t) ifelse(t%in%temp3$MotherID,1,0)    ))
    }   )
    
    totobs_eachyear_f[y] <- length(allpos)
    
    ############################
    ##...and for males
    temp <- subset(rat,(Year%in%ys)&(SexFI=="M")&(AgeFI=="A"))
    allpos <- sort(unique(temp$PJM.ID))
    temp2 <- subset(parent_pairings,Year==allyears[y])
    temp2 <- temp2[!is.na(temp2$FatherID),]
    temp3 <- subset(parent,Year==allyears[y])
    allpos2 <- sort(unique(temp3$FatherID))
    
    check[y] <- all(temp2$FatherID %in% allpos)  # check to make sure all the known-parent females that year are in the master dataset for that year
    
    meanpairingsm_year[y] <- mean(sapply(allpos,function(t) length(which(temp2$FatherID==t))))
    meanoffspringm_year[y] <- mean(sapply(allpos,function(t) length(which(temp3$FatherID==t))))
    meanoffspringm2_year[y] <- mean(sapply(allpos2,function(t) length(which(temp3$FatherID==t))))
    fracreprodm_year[y] <- mean(sapply(allpos,function(t) ifelse(t%in%temp3$FatherID,1,0)  ))
    
    alloffspringdatam <- c(alloffspringdatam,sapply(allpos,function(t) length(which(temp3$FatherID==t))))
    
    meanpairingsm_gtyear[,y] <- sapply(Species,function(z) {
      sub1 <- subset(temp,Species==z)
      allpos <- sort(unique(sub1$PJM.ID))
      mean(sapply(allpos,function(t) length(which(temp2$FatherID==t))))
    }   )
    
    meanoffspringm_gtyear[,y] <- sapply(Species,function(z) {
      sub1 <- subset(temp,Species==z)
      allpos <- sort(unique(sub1$PJM.ID))
      mean(sapply(allpos,function(t) length(which(temp3$FatherID==t))))
    }   )
    
    meanoffspringm2_gtyear[,y] <- sapply(Species,function(z) {
      sub1 <- subset(temp3,Dad.spec==z)
      allpos <- sort(unique(sub1$FatherID))
      if(length(allpos)>0){ mean(sapply(allpos,function(t) length(which(temp3$FatherID==t)))) }else{NA}
    }   )
    
    fracreprodm_gtyear[,y] <- sapply(Species,function(z) {
      sub1 <- subset(temp,Species==z)
      allpos <- sort(unique(sub1$PJM.ID))
      mean(sapply(allpos,function(t) ifelse(t%in%temp3$FatherID,1,0)    ))
    }   )
    
    totobs_eachyear_m[y] <- length(allpos)   # use for weighted means
    
  }
  
  alloffspringdataf
  alloffspringdatam
  
  meanpairingsf_year   # by year
  meanoffspringf_year
  meanoffspringf2_year
  fracreprodf_year
  
  rownames(meanpairingsf_gtyear) <- Species
  colnames(meanpairingsf_gtyear) <- allyears[-length(allyears)]
  
  rownames(meanoffspringf_gtyear) <- Species
  colnames(meanoffspringf_gtyear) <- allyears[-length(allyears)]
  
  rownames(meanoffspringf2_gtyear) <- Species
  colnames(meanoffspringf2_gtyear) <- allyears[-length(allyears)]
  
  rownames(fracreprodf_gtyear) <- Species
  colnames(fracreprodf_gtyear) <- allyears[-length(allyears)]
  
  rownames(meanpairingsm_gtyear) <- Species
  colnames(meanpairingsm_gtyear) <- allyears[-length(allyears)]
  
  rownames(meanoffspringm_gtyear) <- Species
  colnames(meanoffspringm_gtyear) <- allyears[-length(allyears)]
  
  rownames(meanoffspringm2_gtyear) <- Species
  colnames(meanoffspringm2_gtyear) <- allyears[-length(allyears)]
  
  rownames(fracreprodm_gtyear) <- Species
  colnames(fracreprodm_gtyear) <- allyears[-length(allyears)]
  
  meanpairingsf_gtyear
  meanoffspringf_gtyear
  meanoffspringf2_gtyear
  fracreprodf_gtyear
  
  meanpairingsm_gtyear
  meanoffspringm_gtyear
  meanoffspringm2_gtyear
  fracreprodm_gtyear
  
  meanpairingsf_gt <- sapply(1:length(Species),function(t) sum(meanpairingsf_gtyear[Species[t],]*totobs_eachyear_f,na.rm = T)/sum(totobs_eachyear_f,na.rm = T)  )
  meanoffspringf_gt <- sapply(1:length(Species),function(t) sum(meanoffspringf_gtyear[Species[t],]*totobs_eachyear_f,na.rm = T)/sum(totobs_eachyear_f,na.rm = T)  )
  meanoffspringf2_gt <- sapply(1:length(Species),function(t) sum(meanoffspringf2_gtyear[Species[t],]*totobs_eachyear_f,na.rm = T)/sum(totobs_eachyear_f,na.rm = T)  )
  fracreprodf_gt <- sapply(1:length(Species),function(t) sum(fracreprodf_gtyear[Species[t],]*totobs_eachyear_f,na.rm = T)/sum(totobs_eachyear_f,na.rm = T)  )
  
  names(meanpairingsf_gt) <- Species
  names(meanoffspringf_gt) <- Species
  names(meanoffspringf2_gt) <- Species
  names(fracreprodf_gt) <- Species
  
  meanpairingsf_gt
  meanoffspringf_gt
  meanoffspringf2_gt
  fracreprodf_gt
  
  ####### COMPUTE STD ERROR (can do for the other statistics if needed- right now implemented only for offspring per parent)
  
  #t="549"
  thesespecf <- sapply(names(alloffspringdataf),function(t){temp1 <- subset(rat,PJM.ID==t);temp1$Species[1]  }  )
  thesespecm <- sapply(names(alloffspringdatam),function(t){temp1 <- subset(rat,PJM.ID==t);temp1$Species[1]  }  )
  
  t=1
  seoffspringf_gt <- sapply(1:length(Species),function(t) sd(alloffspringdataf[thesespecf==Species[t]],na.rm=T )/sqrt(length(which(!is.na(alloffspringdataf[thesespecf==Species[t]]))))  )
  seoffspringm_gt <- sapply(1:length(Species),function(t) sd(alloffspringdatam[thesespecm==Species[t]],na.rm=T )/sqrt(length(which(!is.na(alloffspringdatam[thesespecm==Species[t]]))))  )
  
  # note: could get se for conditional offspring easily by modifying the above lines
  
  names(seoffspringf_gt) <- Species
  names(seoffspringm_gt) <- Species
  cat(seoffspringf_gt)
  cat(seoffspringm_gt)
  
  # Species    # NOTE: can't compute SE this way- need to go back to the raw data
  # t=1    
  # sepairingsf_gt <- sapply(1:length(Species),function(t) sqrt(weighted.var.se(meanpairingsf_gtyear[Species[t],],totobs_eachyear_f,na.rm = T))  )
  # seoffspringf_gt <- sapply(1:length(Species),function(t) sqrt(weighted.var.se(meanoffspringf_gtyear[Species[t],],totobs_eachyear_f,na.rm = T))  )
  # sefracreprodf_gt <- sapply(1:length(Species),function(t) sum(fracreprodf_gtyear[Species[t],]*totobs_eachyear_f,na.rm = T)/sum(totobs_eachyear_f,na.rm = T)  )
  
  meanpairingsm_gt <- sapply(1:length(Species),function(t) sum(meanpairingsm_gtyear[Species[t],]*totobs_eachyear_m,na.rm = T)/sum(totobs_eachyear_m,na.rm = T)  )
  meanoffspringm_gt <- sapply(1:length(Species),function(t) sum(meanoffspringm_gtyear[Species[t],]*totobs_eachyear_m,na.rm = T)/sum(totobs_eachyear_m,na.rm = T)  )
  meanoffspringm2_gt <- sapply(1:length(Species),function(t) sum(meanoffspringm2_gtyear[Species[t],]*totobs_eachyear_m,na.rm = T)/sum(totobs_eachyear_m,na.rm = T)  )
  fracreprodm_gt <- sapply(1:length(Species),function(t) sum(fracreprodm_gtyear[Species[t],]*totobs_eachyear_m,na.rm = T)/sum(totobs_eachyear_m,na.rm = T)  )
  
  names(meanpairingsm_gt) <- Species
  names(meanoffspringm_gt) <- Species
  names(meanoffspringm2_gt) <- Species
  names(fracreprodm_gt) <- Species
  
  meanpairingsm_gt
  meanoffspringm_gt
  meanoffspringm2_gt
  fracreprodm_gt
  
  ### Save as CSV for table
  thistable <- rbind(fracreprodf_gt,fracreprodm_gt,meanpairingsf_gt,meanpairingsm_gt,meanoffspringf_gt,meanoffspringm_gt,
                     meanoffspringf2_gt,meanoffspringm2_gt)
  thistable <- round(thistable,2)
  rownames(thistable) <- c("Fraction successfully reproducing, F","Fraction successfully reproducing, M",
                           "Mean number of pairings, F","Mean number of pairings, M",
                           "Mean number of offspring, F","Mean number of offspring, M",
                           "Mean number of offspring2, F","Mean number of offspring2, M")
  write.csv(thistable,"summary_table1.csv")
  
  meanpairingsm_gtyear    # by year and genotype
  
  # OPF = offspring per female
  # OPM = offspring per male
  # FRF = fraction reproductive, female
  # FRM = fraction reproductive, male
  # PPF = unique pairings per female
  # PPM = unique pairings per male
  
  reprostats$meanOPF_year <- meanoffspringf_year
  reprostats$meanOPF_gt <- meanoffspringf_gt
  reprostats$meanOPF_gtyr <- meanoffspringf_gtyear
  reprostats$meanOPF2_gtyr <- meanoffspringf2_gtyear
  
  reprostats$meanPPF_year <- meanpairingsf_year
  reprostats$meanPPF_gt <- meanpairingsf_gt
  reprostats$meanPPF_gtyr <- meanpairingsf_gtyear
  
  reprostats$meanFRF_year <- fracreprodf_year
  reprostats$meanFRF_gt <- fracreprodf_gt
  reprostats$meanFRF_gtyr <- fracreprodf_gtyear
  
  reprostats$meanOPM_year <- meanoffspringm_year
  reprostats$meanOPM_gt <- meanoffspringm_gt
  reprostats$meanOPM_gtyr <- meanoffspringm_gtyear
  reprostats$meanOPM2_gtyr <- meanoffspringm2_gtyear
  
  reprostats$meanPPM_year <- meanpairingsm_year
  reprostats$meanPPM_gt <- meanpairingsm_gt
  reprostats$meanPPM_gtyr <- meanpairingsm_gtyear
  
  reprostats$meanFRM_year <- fracreprodm_year
  reprostats$meanFRM_gt <- fracreprodm_gt
  reprostats$meanFRM_gtyr <- fracreprodm_gtyear
  
  return(reprostats)
}


###########################
# FUNCTION: generate distribution of mating pairs under the null expectation


#For each year, and for each known mother/father, draw from potential mates based on empirical distance distribution of matings

generateNullPairings <- function(sex="M",nsims=2){
  
  if(sex=="F"){
    parent <- parent_knownmom
  }else{
    parent <- parent_knowndad
  }
  
  oppsex <- ifelse(sex=="F","M","F")
  
  #First, create empirical distribution of female-male pair distances
  temp <- edfun(na.omit(parent$MateDist))
  f1 <- temp$dfun
  
  svg(sprintf("EmpDist_matedist_%s.svg",sex),4,3)
  curve(f1, 0,200,ylab="Prob Density",xlab="Distance (m)")
  dev.off()
  
  #Select focal individuals (all adult individuals of the selected sex)
  
  temp <- subset(rat,AgeFI=="A"&SexFI==sex)
  focals <- sort(unique(temp$PJM.ID))
  
  exp.out <- NULL
  
  s=1
  for(s in 1:nsims){
    i=1
    for(i in 1:length(focals)){
      this <- focals[i]
      this.rat <- subset(rat,PJM.ID==this&AgeFI=="A") #rat[rat$PJM.ID==this,]
      t.yrs <- unique(this.rat$Year)
      
      #...loop through years
      j=1
      for(j in 1:length(t.yrs)){
        temp.yr <- t.yrs[j]
        
        its.loc <- this.rat[this.rat$Year==t.yrs[j],]
        its.zone <- its.loc$Zone[1] 
        its.loc <- SpatialPoints(its.loc[1,c("UTM_E", "UTM_N")])
        
        #Then, for all potential mates in that year, find their distances to that individual, and calculate probability of being selected
        #based on distance and empirical distribution
        t.mates <- rat.clip[rat.clip@data$Year==temp.yr & rat.clip@data$SexFI==oppsex & rat.clip@data$AgeFI=="A",]
        #Pare down to first location for males
        t.mates <- t.mates[!duplicated(t.mates@data$PJM.ID),]
        t.mates@data$dists <- gDistance(its.loc, t.mates, byid=T)
        t.mates@data$distprob <- f1(t.mates@data$dists)/f1(6.1)  #Max prob at dist 6.1    # KTS: changed 'f' to 'f1'
        #get some NAs here for very far away mates, but I think that's fine, we'll remove them
        t.mates <- t.mates[!is.na(t.mates@data$distprob),]
        #also remove mates with unknown q-values
        t.mates <- t.mates[!is.na(t.mates@data$qFu),]
        
        #Sample a mate with probability from above
        samp.mate <- sample(t.mates@data$PJM.ID, 1, prob=t.mates@data$distprob)
        
        #Find difference in q-value between female and male
        q.this <- rat$qFu[rat$PJM.ID==this][1]
        q.mate <- rat$qFu[rat$PJM.ID==samp.mate][1]
        
        temp.out <- c(s,this, samp.mate, temp.yr, q.this, q.mate, its.zone)
        
        exp.out <- rbind(exp.out, temp.out)
      }
    }
  }
  
  exp.out <- as.data.frame(exp.out)
  
  if(sex=="F"){
    names(exp.out) <- c("Simulation","MotherID", "FatherID", "Year", "qFu.mom", "qFu.dad", "MotherZone")
  }else{
    names(exp.out) <- c("Simulation","FatherID", "MotherID", "Year", "qFu.dad", "qFu.mom",  "FatherZone")
  }
  
  exp.out$qFu.mom <- as.numeric(exp.out$qFu.mom)
  exp.out$qFu.dad <- as.numeric(exp.out$qFu.dad)
  exp.out$diff <- exp.out$qFu.mom - exp.out$qFu.dad
  if(sex=="M")  exp.out$diff <- exp.out$qFu.dad - exp.out$qFu.mom
  exp.out$sqdiff <- exp.out$diff^2
  exp.out$mid_year <- paste0(exp.out$MotherID,exp.out$Year)
  exp.out$pid_year <- paste0(exp.out$FatherID,exp.out$Year)
  exp.out$matpairid <- paste0(exp.out$MotherID,"_",exp.out$FatherID,"_",exp.out$Year)
  
  # names(exp.out)
  # tmp <- subset(exp.out.f,FatherID=="S-515")   # debug
  # tmp
  # 
  exp.out$Mom.spec <- ifelse(exp.out$qFu.mom<=0.1,"Macrotis",ifelse(exp.out$qFu.mom>=0.9,"Fuscipes",ifelse(exp.out$qFu.mom<0.4,"Mac_Back",ifelse(exp.out$qFu.mom>0.6,"Fus_Back","Hybrid"))))
  exp.out$Dad.spec <- ifelse(exp.out$qFu.dad<=0.1,"Macrotis",ifelse(exp.out$qFu.dad>=0.9,"Fuscipes",ifelse(exp.out$qFu.dad<0.4,"Mac_Back",ifelse(exp.out$qFu.dad>0.6,"Fus_Back","Hybrid"))))
  
  # NotE: can integrate and remove this next time we rerun the above
  exp.out$matpairid <- paste0(exp.out$MotherID,"_",exp.out$FatherID,"_",exp.out$Year)
  #exp.out.m$matpairid <- paste0(exp.out.m$MotherID,"_",exp.out.m$FatherID,"_",exp.out.m$Year)
  # names(exp.out.f)[names(exp.out.f)=="MotherZone"] <- "Zone"
  # names(exp.out.m)[names(exp.out.m)=="FatherZone"] <- "Zone"
  
  if(sex=="F"){
    exp.out$Zone <- exp.out$MotherZone
  }else{
    exp.out$Zone <- exp.out$FatherZone
  }
  
  #exp.out.m$Zone <- exp.out.m$FatherZone
  exp.out <- exp.out[complete.cases(exp.out[,c("Mom.spec","Dad.spec")]),]
  #exp.out.m <- exp.out.m[complete.cases(exp.out.m[,c("Mom.spec","Dad.spec")]),]
  # sort(unique(exp.out$Year))
  # sort(unique(rat$Year))
  allyears <- intersect(sort(unique(exp.out$Year)),sort(unique(parent$Year)))
  exp.out <- subset(exp.out,Year%in%allyears) 
  #exp.out.m <- subset(exp.out.m,Year%in%allyears)
  
  return(exp.out)
  
}

####debug

# idyr <- 5582008 #"S-3722010"
# sex="M"
# if(sex=="F"){
#   exp.out.f[exp.out.f$mid_year==idyr,]
#   #exp.out.m[exp.out.m$mid_year==idyr,]
# }else{
#   #exp.out.f[exp.out.f$pid_year==idyr,]
#   exp.out.m[exp.out.m$pid_year==idyr,]
# }
# 
# rat[rat$id_year==idyr,]



#####################
# FUNCTION: apparent mate selection tests

testmateselection <- function(sex="F",zone="all"){
  
  exp.out <- thin_knownrep(sex)
  if(sex=="F"){
    idcol <- "mid_year"
  }else{
    idcol <- "pid_year"
  }
  
  temp <- subset(exp.out,Simulation==2)   
  
  nrow(temp)
  
  ndx <- exp.out[[idcol]]%in%parent[[idcol]]
  any(!ndx)   # all true [must be true]
  
  ndx <- parent[[idcol]]%in%exp.out[[idcol]]
  any(!ndx)   # not all true- these are the unknown sex individuals
  
  # parent[!ndx,]
  # parent[806,]
  # rat[rat$PJM.ID=="S-973",]
  
  parent2_1 <- parent[ndx,]
  
  length(unique(parent2_1[[idcol]]))     # 421 unique 
  length(unique(exp.out[[idcol]]))     # 421 unique -- it should match
  
  #######
  # subset by zone
  
  zonevar <- ifelse(sex=="F","MotherZone","FatherZone")
  allzones <- sort(unique(parent[[zonevar]]))
  thesezones <- zone
  if(zone=="all") thesezones <- allzones
  
  exp.out2 <- exp.out[exp.out[[zonevar]]%in%thesezones,]
  
  parent2 <- parent2_1[parent2_1[[zonevar]]%in%thesezones,]
  
  ######
  # null distribution: mean squared difference
  
  nulldist_sqdiff <- exp.out2 %>% group_by(Simulation) %>% 
    summarize(meansqdiff = mean(sqdiff))
  nsims <- nrow(nulldist_sqdiff)
  nrow(exp.out)
  
  ######
  # observed value: mean squared difference
  #idcolnum <- which(names(parent2)==idcol)
  obsdf <- parent2 %>% group_by_at(idcol) %>% 
    summarize(sqdiff = sqdiff[which(!is.na(sqdiff))[1]])       # KTS: keep only the first observation (this is a lazy shortcut but may be okay- we can talk about other ways to do this)
  nrow(obsdf)
  
  real_sqdiff <- mean(obsdf$sqdiff,na.rm=T)
  
  #######
  # statistical test: mean squared difference
  
  svg(sprintf("globaltest_nonrandom_mating_%s_%s.svg",sex,zone),4,4)
  hist(nulldist_sqdiff$meansqdiff,xlim=c(0,0.5),
       main = "",
       xlab="Mean squared dif in q values in matings")  # null distribution of squared differences in q val between mother and father
  abline(v=real_sqdiff,col="blue",lwd=2)
  pval <- length(which(nulldist_sqdiff$meansqdiff<real_sqdiff))/nsims
  legend("topright",bty="n",pch="",legend=sprintf("p = %s",round(pval,3)))
  dev.off()
  
  #######
  # Next: do the same thing but with the different 'species' (macrotis, fuscipes, hybrid, backcross)
  
  specvar <- ifelse(sex=="F","Mom.spec","Dad.spec")
  table(exp.out[[specvar]])
  
  matespecvar <- ifelse(sex=="F","Dad.spec","Mom.spec")
  
  ######
  # null distribution: mean squared difference and mean difference
  
  if(sex=="F"){
    nulldist_sqdiff <- exp.out2 %>% group_by(Simulation,Mom.spec) %>% 
      summarize(meansqdiff = mean(sqdiff),
                meandiff = mean(diff),
                meanmateq = mean(qFu.dad)
      )
  }else{
    nulldist_sqdiff <- exp.out2 %>% group_by(Simulation,Dad.spec) %>% 
      summarize(meansqdiff = mean(sqdiff),
                meandiff = mean(diff),
                meanmateq = mean(qFu.mom)
      )
  }
  
  nrow(nulldist_sqdiff)
  
  ######
  # observed value: mean squared difference and mean difference
  
  table(parent2$Mom.spec)
  
  if(sex=="F"){
    obsdf <- parent2 %>% group_by(mid_year,Mom.spec) %>% 
      summarize(sqdiff = sqdiff[which(!is.na(sqdiff))[1]],
                diff = diff[which(!is.na(diff))[1]],
                mateq = qFu.dad[which(!is.na(qFu.dad))[1]]
      )       # KTS: keep only the first observation (this is a lazy shortcut but may be okay- we can talk about other ways to do this)
  }else{
    obsdf <- parent2 %>% group_by(mid_year,Dad.spec) %>% 
      summarize(sqdiff = sqdiff[which(!is.na(sqdiff))[1]],
                diff = diff[which(!is.na(diff))[1]],
                mateq = qFu.dad[which(!is.na(qFu.mom))[1]]
      ) 
  }
 
  nrow(obsdf)
  
  real_sqdiff <- tapply(obsdf$sqdiff,obsdf[[specvar]],mean,na.rm=T)
  real_diff <- tapply(obsdf$diff,obsdf[[specvar]],mean,na.rm=T)
  real_mateq <- tapply(obsdf$mateq,obsdf[[specvar]],mean,na.rm=T)
  
  #######
  # statistical test: mean squared difference and mean difference
  
  #Species <- c("Fuscipes","Macrotis","Hybrid","Fus_Back","Mac_Back")
  
  names(speccol) <- Species
  
  # layout(matrix(1:6,ncol=2))
  
  # s=Species[3]
  # for(s in Species){
  #   nulldist <- subset(nulldist_sqdiff,Mom.spec==s)
  #   hist(nulldist$meansqdiff,xlim=c(0,0.5),
  #        main = s,
  #        xlab="Mean squared dif in q values in matings",freq=F)  # null distribution of squared differences in q val between mother and father
  #   abline(v=real_sqdiff[s],col="blue",lwd=2)
  #   pval <- length(which(nulldist$meansqdiff<real_sqdiff[s]))/nsims
  #   text(0.3,15,paste0("p = ",round(pval,3)))
  # }
  
  table(nulldist_sqdiff[[specvar]])
  
  svg(sprintf("ApparentMateSelTest_%s_%s_forinkscape.svg",sex,zone),6,3.5)
  
  layout(matrix(1:6,ncol=3,byrow=F))
  
  par(mai=c(0.7,0.7,0.05,0.1))
  s=Species[2]
  for(s in Species[1:5]){
    s3 <- Species3[Species==s]
    if(sex=="F"){
      nulldist <- subset(nulldist_sqdiff,Mom.spec==s)
    }else{
      nulldist <- subset(nulldist_sqdiff,Dad.spec==s)
    }
    
    hist(nulldist$meanmateq,breaks=8,freq=T,xlim=c(0,1),main = "",xaxt="n",yaxt="n",
         xlab="",ylab = "Frequency",col=gray(0.6),border=gray(0.6),ylim=c(0,15))
    axis(1,at=c(0,0.5,1),labels=c(0,0.5,1))
    axis(2,at=c(0,5,10,15),labels=c(0,5,10,15))
    mtext(expression('q'[italic("fuscipes")]), side=1, line=3.2, cex=1.2)
    # polygon(c(0,0.1,0.1,0),c(0,0,50,50),col=speccol[5],border=speccol[5])
    # polygon(c(0.1,0.4,0.4,0.1),c(0,0,50,50),col=speccol[4],border=speccol[4])
    # polygon(c(0.4,0.6,0.6,0.4),c(0,0,50,50),col=speccol[3],border=speccol[3])
    # polygon(c(0.6,0.9,0.9,0.6),c(0,0,50,50),col=speccol[2],border=speccol[2])
    # polygon(c(0.9,1,1,0.9),c(0,0,50,50),col=speccol[1],border=speccol[1])
    # hist(nulldist$meanmateq,
    #      add=T,
    #      breaks=8,
    #      freq=F,col=gray(0.6),border=gray(0.9))  # null distribution of differences in q val between mother and father
    abline(v=real_mateq[s],col="red3",lwd=3)
    legend("topleft",pch="",legend=s3)
    pval <- min(length(which(nulldist$meanmateq<real_mateq[s]))/nsims,length(which(nulldist$meanmateq<real_mateq[s]))/nsims)
    text(-0.3,3,paste0("p = ",round(pval,3)))
  }
  # par(mai=c(0,0,0,0))
  # plot(1,1,pch="",xaxt="n",xlab="",ylab="",yaxt="n",bty="n")
  # legend("top",fill=speccol,bty="n",legend=Species3,cex=1.5)
  
  dev.off()
  
  ########
  # look at fraction of available mate pairings by genotype
  
  ndx <- !is.na(exp.out2[[matespecvar]])
  table(ndx)
  exp.out2 <- exp.out2[ndx,]
  
  ndx <- !is.na(parent2[[matespecvar]])
  table(ndx)
  parent2 <- parent2[ndx,]
  
  if(sex=="F"){
    permatepair_sim <- sapply(Species,function(t){temp <- subset(exp.out2,Mom.spec==t);temp$Dad.spec <- factor(temp$Dad.spec,levels=Species) ;  as.numeric(table(temp$Dad.spec)/nrow(temp))  }  )
    rownames(permatepair_sim) = Species
    
    permatepair_obs <- sapply(Species,function(t){temp <- subset(parent2,Mom.spec==t);temp$Dad.spec <- factor(temp$Dad.spec,levels=Species) ;  as.numeric(table(temp$Dad.spec)/nrow(temp))  }  )
    rownames(permatepair_obs) = Species
  }else{
    permatepair_sim <- sapply(Species,function(t){temp <- subset(exp.out2,Dad.spec==t);temp$Mom.spec <- factor(temp$Mom.spec,levels=Species) ;  as.numeric(table(temp$Mom.spec)/nrow(temp))  }  )
    rownames(permatepair_sim) = Species
    
    permatepair_obs <- sapply(Species,function(t){temp <- subset(parent2,Dad.spec==t);temp$Mom.spec <- factor(temp$Mom.spec,levels=Species) ;  as.numeric(table(temp$Mom.spec)/nrow(temp))  }  )
    rownames(permatepair_obs) = Species
  }
  
  apply(permatepair_obs,2,sum)

  svg(sprintf("Avail_vs_obs_pairings_%s_%s.svg",sex,zone),6,7)
  layout(matrix(1:2,nrow=2))
  par(las=1)
  par(mai=c(1.5,0.8,0.4,0.8))
  x<-barplot(permatepair_sim*100, col=speccol, xlim=c(0,8),main="Available pairings",
             border="white", xlab="", #sprintf("%s type",ifelse(sex=="F","Maternal","Paternal")), 
             ylab="% available mates",las=3, names.arg=Species3)
  legend("topright",fill=speccol,legend = Species3)
  
  
  x<-barplot(permatepair_obs*100, col=speccol, xlim=c(0,8),main="Observed pairings",
             border="white", xlab="", #sprintf("%s type",ifelse(sex=="F","Maternal","Paternal")), 
             ylab="% observed mates",las=3,names.arg=Species3)
  legend("topright",fill=speccol,legend = Species3)
  dev.off()
  
}



####################
# FUNCTION: merge null hypothesis datasets

   # note that 'diff' (q difference) has different meaning for male and female simulated datasets- changed here

mergesimdata <- function(){
  temp1 <- exp.out.f         # first merge the two exp.out datasets (to do for all sexes)
  temp2 <- exp.out.m
  temp2$diff <- -1*temp2$diff   # put this variable on the same playing field (mom q minus dad q)
  cmn <- intersect(names(temp1),names(temp2))
  temp2 <- temp2[,cmn]
  exp.out <- rbind(temp1[,cmn],temp2[,cmn])
  return(exp.out)
}
  



######################
# FUNCTION: make 'pairings' dataset to look at all pairings types and their reproductive output/success

# NOTE: expected pairings will be filled in later after simulations are performed.
# expected offspring (offspring.exp) determined assuming number of offspring equal across all pairing types and proportional to availability of each pairing type

# NOTE: could also determine expected offspring in other ways. Another way is to determine the number of actual pairs in the parent dataset and multiply by the offspring per female
# RESULT: no pairings have more representation in the 

makepairingsdf <- function(){
  
  exp.out <- mergesimdata()
  nrow(exp.out)
  
  ## determine the frequency of different pairings
  pairings <- expand.grid(Species,Species)
  names(pairings) <- c("Mom.spec","Dad.spec")    # all possible pairings
  
  pairings$pairs.obs <- NA
  pairings$pairs.exp <- NA
  pairings$offspring.obs <- NA
  pairings$offspring.exp <- NA
  pairings$exp.freq <- NA
  pairings$offsp.fracfem <- NA
  # pairings$freq.obs <- NA
  # pairings$freq.exp <- NA
  # 
  
  exp.out.nona <- exp.out[complete.cases(exp.out[,c("Mom.spec","Dad.spec")]),]
  nrow(exp.out.nona)
  
  thispairingtype <- 1
  for(thispairingtype in 1:nrow(pairings)){
    thismom <- pairings$Mom.spec[thispairingtype]
    thisdad <- pairings$Dad.spec[thispairingtype]
    pairings$pairs.obs[thispairingtype] <- length(which((parent_pairings$Mom.spec==thismom) & (parent_pairings$Dad.spec==thisdad)))
    pairings$offspring.obs[thispairingtype] <- length(which((parent$Mom.spec==thismom) & (parent$Dad.spec==thisdad)))
    pairings$offsp.fracfem[thispairingtype] <- length(which((parent$Mom.spec==thismom) & (parent$Dad.spec==thisdad) & parent$sex.offspring=="F"))/pairings$offspring.obs[thispairingtype]
    #pairings$freq.obs[thispairingtype] <- 
    #pairings$obsfreq[thispairingtype] <- pairings$offspring.obs[thispairingtype]/nrow(parent)
  }
  
  
  thispairingtype <- 1
  for(thispairingtype in 1:nrow(pairings)){
    thisfrac <- length(which((exp.out.nona$Mom.spec==pairings$Mom.spec[thispairingtype])&(exp.out.nona$Dad.spec==pairings$Dad.spec[thispairingtype])))/nrow(exp.out.nona)
    pairings$pairs.exp[thispairingtype] <- round(thisfrac*sum(pairings$pairs.obs),1)
    pairings$offspring.exp[thispairingtype] <- round(thisfrac*sum(pairings$offspring.obs),1)
    pairings$exp.freq[thispairingtype] <- thisfrac
  }
  
  write.csv(pairings,"pairings_obs_exp3.csv",row.names = F)
  
  return(pairings)
}


###############
# Function: make merged df with both observed and expected pairings

makemergeddf <- function(parent=parent_pairings,zone="all",sex="all"){
  
  if(sex=="all"){
    exp.out <- mergesimdata()
  }else{
    if(sex=="F"){
      exp.out <- exp.out.f
    }else{
      exp.out <- exp.out.m
    }
  }
  
  # table(exp.out$Mom.spec)
  # table(exp.out$Dad.spec)
  
  # if(sex=="F"){
  #   allzones <- sort(unique(parent$MotherZone))
  # }else{
  #   allzones <- sort(unique(parent$FatherZone))
  # }
  
  allzones <- sort(unique(exp.out$Zone))
  
  if(zone=="all"){
    thesezones <- allzones #[1]  #[2]  #[2] #[1]   #[2]   # hybrid zone best for this test?
  }else{
    thesezones <- zone
  } 
  
  parent$truemat  <- 1
  exp.out$truemat <- 0
  
  names(exp.out)
  parent$Simulation <- 0
  
  parent$Zone <- parent$MotherZone
  
  commn <- intersect(names(parent),names(exp.out))
  
  mergeddf <- rbind(parent[,commn],exp.out[,commn])
  
  names(mergeddf)
  nrow(mergeddf)
  
  # if(sex=="F"){
  #   mergeddf2 <- subset(mergeddf,(MotherZone%in%thesezones))  #&((Mom.spec%in%thisspec)|(Dad.spec%in%thisspec))
  # }else{
  #   mergeddf2 <- subset(mergeddf,(FatherZone%in%thesezones))  #&((Mom.spec%in%thisspec)|(Dad.spec%in%thisspec))
  # }
  mergeddf2 <- subset(mergeddf,(Zone%in%thesezones))  #&((
  
  pure <- c("Macrotis","Fuscipes")
  
  mergeddf2$pureparent <- ifelse((mergeddf2$Mom.spec%in%pure)|(mergeddf2$Dad.spec%in%pure),1,0 )
  mergeddf2$bothpure <- ifelse((mergeddf2$Mom.spec%in%pure)&(mergeddf2$Dad.spec%in%pure),1,0 )
  mergeddf2$bothpurex <- ifelse((mergeddf2$Mom.spec%in%pure)&(mergeddf2$Dad.spec%in%pure)&(mergeddf2$Mom.spec!=mergeddf2$Dad.spec),1,0 )
  mergeddf2$bothpures <- ifelse((mergeddf2$Mom.spec%in%pure)&(mergeddf2$Dad.spec%in%pure)&(mergeddf2$Mom.spec==mergeddf2$Dad.spec),1,0 )
  mergeddf2$bothhyb <- ifelse((!mergeddf2$Mom.spec%in%pure)&(!mergeddf2$Dad.spec%in%pure),1,0 )
  
  if(sex=="all"){
    mergeddf2$hybdad <- ifelse((!mergeddf2$Dad.spec%in%pure),1,0 )
    mergeddf2$hybmom <- ifelse((!mergeddf2$Mom.spec%in%pure),1,0 )
  }else{
    if(sex=="F"){
      mergeddf2$hybdad <- ifelse((!mergeddf2$Dad.spec%in%pure),1,0 )
    }else{
      mergeddf2$hybmom <- ifelse((!mergeddf2$Mom.spec%in%pure),1,0 )    # can only compute this for the opposite sex
    }
  }
  
  return(mergeddf2)
}


#################
# Function: run tests to see if observed mate purity differs from expected 

test_obs_vs_exp_purity <- function(var="bothpure",mergeddf = mergeddf_f, nulldists=nulldists_f,
                                   hyp="less",lab="freq both pure of same species",main=""){
  
  tbl <- table(mergeddf[[var]],mergeddf$truemat)
  colnames(tbl) <- c("exp","obs")
  rownames(tbl) <- c("cond false","cond true")
  #print(tbl)
  print(chisq.test(tbl))
  cls <- apply(tbl,2,sum)
  obs.bothpure <- round(tbl[2,2]/cls[2],3)
  exp.bothpure <- round(tbl[2,1]/cls[1],3)
  cat(sprintf("observed: %s \nexpected: %s \n",obs.bothpure,exp.bothpure))
  if(hyp=="greater"){
    p1 = length(which(nulldists[[var]]>obs.bothpure))/nrow(nulldists) 
  }else{
    p1 = length(which(nulldists[[var]]<obs.bothpure))/nrow(nulldists)
  }
  alltemp <- c(nulldists[[var]],obs.bothpure)
  hist(nulldists[[var]],xlim=c(min(alltemp)*0.85,min(1,max(alltemp)*1.2)),breaks=8,main=main,xlab=lab)
  abline(v=obs.bothpure,col="blue",lwd=3)
  legend("topleft",pch="",legend=sprintf("p = %#.3f",p1),bty="n")
  
}


################
# Function: make null distributions

makenulldists <- function(mergeddf=mergeddf_all,sex="all"){
  null1 <- subset(mergeddf,truemat==0)
  
  if(sex=="all"){
    nulldist_pureparent <- null1 %>% group_by(Simulation) %>% 
      summarize(pureparent = mean(pureparent),
                bothpure = mean(bothpure),
                bothpurex = mean(bothpurex),
                bothpures = mean(bothpures),
                bothhyb = mean(bothhyb),
                hybmom = mean(hybmom),
                hybdad = mean(hybdad)
      )
  }else{
    if(sex=="F"){
      nulldist_pureparent <- null1 %>% group_by(Simulation) %>% 
        summarize(pureparent = mean(pureparent),
                  bothpure = mean(bothpure),
                  bothpurex = mean(bothpurex),
                  bothpures = mean(bothpures),
                  bothhyb = mean(bothhyb),
                  #hybmom = mean(hybmom),
                  hybdad = mean(hybdad)
        )
    }else{
      nulldist_pureparent <- null1 %>% group_by(Simulation) %>% 
        summarize(pureparent = mean(pureparent),
                  bothpure = mean(bothpure),
                  bothpurex = mean(bothpurex),
                  bothpures = mean(bothpures),
                  bothhyb = mean(bothhyb),
                  hybmom = mean(hybmom)
                  #hybdad = mean(hybdad)
        )
    }
  }
  nrow(nulldist_pureparent)
  return(nulldist_pureparent)
}


################
# Function: make figure "purepar_all.svg"

make_purepar_all <- function(){
  name <- "purepar_all.svg"
  svg(name,6,5)
  
  par(mai=c(.8,.8,.2,.1))
  par(mfcol=c(2,2))
  
  test_obs_vs_exp_purity("bothpures",mergeddf_all,nulldists_all,"greater","freq both pure parents same species")
  
  test_obs_vs_exp_purity("bothpurex",mergeddf_all,nulldists_all,"less","freq both pure parents dif species")
  
  test_obs_vs_exp_purity("bothpure",mergeddf_all,nulldists_all,"greater","freq at least one pure parent")
  
  test_obs_vs_exp_purity("bothhyb",mergeddf_all,nulldists_all,"less","freq both hybrid parents")
  
  dev.off()
}


###################
# Function: test if hybrid males and females tend to mate with pure individuals of opposite sex.

make_hybrid_puremat <- function(){
  
  # female focus: test if hybrid females are more likely to have pure male partner
  # male focus: test if hybrid males are more likely to have pure female partner
  
  mergeddf_f2 <- subset(mergeddf_f,!Mom.spec%in%pure)   # only non-pure moms
  nrow(mergeddf_f2)
  
  mergeddf_m2 <- subset(mergeddf_m,!Dad.spec%in%pure)   # only non-pure dads
  nrow(mergeddf_m2)
  
  nulldists_f2 <- makenulldists(mergeddf_f2,sex='F')
  nulldists_m2 <- makenulldists(mergeddf_m2,sex='M')
  
  name = "puremat_hybrid.svg"
  svg(name,5.5,4)
  
  par(mfcol = c(1,2))
  
  test_obs_vs_exp_purity("pureparent",mergeddf_f2,nulldists_f2,"greater","frequency pure mate","Non-pure females")
  test_obs_vs_exp_purity("pureparent",mergeddf_m2,nulldists_m2,"greater","frequency pure mate","Non-pure males")
  
  dev.off()
}



####################
# FUNCTION: test if cross-pairings with Fuscipes are more rare than cross pairings with macrotis

# strategy: combine the 'parent_pairings' and corresponding 'exp.out' datasets
#             then remove all hybrids, keeping only all-pure pairings
#             then run GLM to test strength of 'selection' one way or another 


testcrosspairings <- function(sex="F"){
  mergeddf <- makemergeddf(parent = parent_pairings, zone="all",sex = sex)   # use female only because we are just focusing on female "choice
  nrow(mergeddf)
  
  # names(mergeddf)
  # nrow(mergeddf)
  
  mergeddf2 <- subset(mergeddf,(Mom.spec%in%pure)&(Dad.spec%in%pure))
  nrow(mergeddf2)
  
  mergeddf2$cross <- as.factor(ifelse((mergeddf2$Mom.spec%in%pure)&(mergeddf2$Dad.spec%in%pure)&(mergeddf2$Mom.spec!=mergeddf2$Dad.spec),1,0 ))
  if(sex%in%c("F","all")){
    mergeddf2$fufoc <- as.factor(ifelse(mergeddf2$Mom.spec=="Fuscipes",1,0))
  }else{
    mergeddf2$fufoc <- as.factor(ifelse(mergeddf2$Dad.spec=="Fuscipes",1,0))
  }
  
  #mergeddf2$noncross <- as.factor(ifelse(mergeddf2$cross!=1,1,0))
  #mergeddf2$crossxspec <- as.factor(ifelse(mergeddf2$noncross==1,"same",ifelse(mergeddf2$)))
  
  mergeddf3 <- subset(mergeddf2,cross==1)
  
  thistab <- table(mergeddf3$fufoc,mergeddf3$truemat)
  colnames(thistab) <- c("Sim","Obs")
  rownames(thistab) <- c("MA","FU")
  print(thistab)
  
  colsums <- colSums(thistab)
  temp <- round(thistab[2,]/colsums,3)
  
  cat(sprintf("Observed: %s\nSimulated: %s\n",temp[2],temp[1]))
  
  fisher.test(thistab,alternative = "two.sided")
  
  
  # mod1 <- glm(truemat~cross,data=mergeddf2,family="binomial")
  # summary(mod1)
  # 
  # mod2 <- glm(truemat~cross*fufem,data=mergeddf2,family="binomial")   # weak evidence for stronger avoidance of fuscipes males by macrotis
  # summary(mod2)      # 
  # 
  # mod3 <- glm(truemat~fufem,data=mergeddf3,family="binomial")   # no evidence for stronger avoidance in one way than another
  # summary(mod3)
  # 
  # library(lme4)
  # mod3 <- glmer(truemat~cross*fufem+(1|strata),data=mergeddf2,family="binomial")  #no evidence for stronger avoidance in one way than another# KTS: does it matter that there are some strata with no true pairings?
  # summary(mod3)
  # 
  # nd <- data.frame(
  #   cross = as.factor(c(0,0,1,1)),
  #   fufem = as.factor(c(0,1,0,1))
  # )
  # 
  # allstrat <- levels(mergeddf2$strata)
  # nd$strata <- factor(allstrat[100],levels=allstrat)
  # 
  # predict(mod2,nd,type="response")
  
}



################
# FUNCTION: characterize availability

make_avlvar <- function(spec="Fuscipes"){
  #table(exp.out.f$mid_year)  # check
  allfems <- sort(unique(exp.out.f$mid_year))
  #table(exp.out.m$pid_year)
  allmales <- sort(unique(exp.out.m$pid_year))
  allinds <- c(allfems,allmales)
  allyears <- intersect(sort(unique(rat$Year)),sort(unique(parent$Year)))
  
  avl_f <- numeric(length(allfems))
  names(avl_f) <- allfems
  f=1
  for(f in 1:length(allfems)){
    temp <- subset(exp.out.f,mid_year==allfems[f])
    temp2 <- subset(temp,Dad.spec==spec)
    avl_f[f] <- nrow(temp2)/nrow(temp)
  }
  
  avl_m <- numeric(length(allmales))
  names(avl_m) <- allmales
  m=1
  for(m in 1:length(allmales)){
    temp <- subset(exp.out.m,pid_year==allmales[m])
    temp2 <- subset(temp,Mom.spec==spec)
    avl_m[m] <- nrow(temp2)/nrow(temp)
  }
  avl <- c(avl_f,avl_m)
  return(avl)
}



#########################
# FUNCTION: make data frame of reproductive success by individual

make_repdf <- function(){
  
  allfems <- sort(unique(exp.out.f$mid_year))
  #table(exp.out.m$pid_year)
  allmales <- sort(unique(exp.out.m$pid_year))
  allinds <- c(allfems,allmales)
  allyears <- intersect(sort(unique(rat$Year)),sort(unique(parent$Year)))
  
  repdf <- data.frame(id_year=allinds)
  rownames(repdf) <- allinds
  repdf$year <- NA
  repdf$sex <- NA
  repdf$fuavl <- NA
  repdf$maavl <- NA
  repdf$didrep <- NA
  repdf$species <- NA
  repdf$qFu <- NA
  repdf$size <- NA
  
  i=1
  for(i in 1:length(allinds)){
    thisind <- allinds[i]
    temp1 <- subset(rat,id_year==thisind)
    temp2 <- subset(parent,(mid_year==thisind)|(pid_year==thisind))
    repdf$year[i] <- temp1$Year[!is.na(temp1$Year)][1]
    repdf$sex[i] <- temp1$SexFI[!is.na(temp1$SexFI)][1]
    repdf$fuavl[i] <- fuavl[names(fuavl)==thisind][1]
    repdf$maavl[i] <- maavl[names(maavl)==thisind][1]
    repdf$didrep[i] <- ifelse(nrow(temp2)>0,1,0)
    repdf$species[i] <- temp1$Species[!is.na(temp1$Species)][1]
    repdf$qFu[i] <- temp1$qFu[!is.na(temp1$qFu)][1]
    repdf$size[i] <- ifelse(any(!is.na(temp1$WtMax_g)), max(temp1$WtMax_g[!is.na(temp1$WtMax_g)]),NA)#[1]
  }
  nrow(repdf)
  repdf$species2 <- repdf$species
  repdf$species2[!repdf$species%in%pure] <- 'Hybrid'
  
  repdf$species <- factor(repdf$species,levels = Species)
  repdf$species2 <- factor(repdf$species2,levels = c("Fuscipes","Hybrid","Macrotis"))
  
  repdf$nonpure <- 1-(repdf$fuavl+repdf$maavl)
  
  return(repdf)
}

#################
# FUNCTION: plot size by species

plot_specbysize <- function(){
  svg("specbysize_adults.svg",5.5,4)
  #repdf$didrep
  par(mfrow=c(1,2))     # look at sizes by species
  plot(repdf_M$size~repdf_M$species,ylab="size",xlab="",main="male",las=2,ylim=c(0,500))
  plot(repdf_F$size~repdf_F$species,ylab="size",xlab="",main="female",las=2,ylim=c(0,500))
  dev.off()
  
  Male <- round(tapply(repdf_M$size,repdf_M$species,median,na.rm=T),2)
  Female <- round(tapply(repdf_F$size,repdf_F$species,median,na.rm=T),2)
  write.csv(rbind(Male,Female),file="mediansizebygt.csv")
  # M2 <- round(tapply(repdf_M$size,repdf_M$species,mean,na.rm=T),2)
  # F2 <- round(tapply(repdf_F$size,repdf_F$species,mean,na.rm=T),2)
  # 
  # rbind(M2,F2)
  # rbind(Male,Female)
  
  Male2 <- tapply(repdf_M$size,repdf_M$species,function(x) round(c("min"=min(x,na.rm=T),quantile(x,c(0.25,0.5,0.75),na.rm=T),"max"=max(x,na.rm=T)),1)   )
  Male2 <- as.data.frame(do.call(rbind,Male2))
  Female2 <- tapply(repdf_F$size,repdf_F$species,function(x) round(c("min"=min(x,na.rm=T),quantile(x,c(0.25,0.5,0.75),na.rm=T),"max"=max(x,na.rm=T)),1)   )
  Female2 <- as.data.frame(do.call(rbind,Female2))
  Male2$Sex = "Male"
  Female2$Sex = "Female"
  nms <- names(Female2)[c(6,1:5)]
  Male2 <- Male2[,nms]
  Female2 <- Female2[,nms]
  towrite <- rbind(Male2,Female2)
  write.csv(towrite,file="size_data_table.csv")

}


#################
# FUNCTION: plot pure mate availability over time


plot_avlbytime <- function(){
  svg("mateavl_bytime.svg",6,2.5)
  par(mfrow=c(1,3))
  par(mai=c(0.8,0.6,0.5,0.1))
  sub1 <- subset(repdf,species2=="Fuscipes"&sex=="F")
  plot(jitter(sub1$fuavl)~jitter(sub1$year),main=expression(italic("N. fuscipes")),
       ylab="Fraction conspecific mates",xlab="Year")
  m1 <- lm(asin(sqrt(fuavl))~year,data=sub1)
  summary(m1)
  nd <- data.frame(year=seq(min(allyears)-1,max(allyears),length.out = 100))
  pred <- predict(m1,newdata=nd,interval="confidence")
  pred <- sin(pred)^2
  lines(nd[,1],pred[,1],lwd=2)
  lines(nd[,1],pred[,2],lwd=1,lty=2)
  lines(nd[,1],pred[,3],lwd=1,lty=2)
  fig_label("A", cex=1.5)
  
  sub1 <- subset(repdf,species2=="Macrotis"&sex=="F")
  plot(jitter(sub1$maavl)~jitter(sub1$year),main=expression(italic("N. macrotis")),
       ylab="Fraction conspecific mates",xlab="Year")
  m1 <- lm(asin(sqrt(maavl))~year,data=sub1)
  summary(m1)
  nd <- data.frame(year=seq(min(allyears)-1,max(allyears),length.out = 100))
  pred <- predict(m1,newdata=nd,interval="confidence")
  pred <- sin(pred)^2
  lines(nd[,1],pred[,1],lwd=2)
  lines(nd[,1],pred[,2],lwd=1,lty=2)
  lines(nd[,1],pred[,3],lwd=1,lty=2)
  fig_label("B", cex=1.5)
  
  sub1 <- subset(repdf,sex=="F")
  plot(jitter(sub1$nonpure)~jitter(sub1$year),main="All females",
       ylab="Fraction admixed mates available",xlab="Year")
  m1 <- lm(asin(sqrt(nonpure))~year,data=sub1)
  summary(m1)
  nd <- data.frame(year=seq(min(allyears)-1,max(allyears),length.out = 100))
  pred <- predict(m1,newdata=nd,interval="confidence")
  pred <- sin(pred)^2
  lines(nd[,1],pred[,1],lwd=2)
  lines(nd[,1],pred[,2],lwd=1,lty=2)
  lines(nd[,1],pred[,3],lwd=1,lty=2)
  fig_label("C", cex=1.5)
  
  dev.off()
}


#################
# FUNCTION: model prob of repro

model_rep_prob <- function(){
  ####FEMALES
  
  repdf_F <- repdf_F[complete.cases(repdf_F),]
  
  ## save this
  model1_F <<- glm(didrep~species2+size+nonpure, data=repdf_F,family="binomial")   # best model
  summary(model1_F)
  
  ## save this
  model2_F <<- glm(didrep~species2*size+nonpure, data=repdf_F,family="binomial")
  summary(model2_F)

  
  AIC(model1_F,model2_F)
  
  
  ####MALES
  
  repdf_M <- repdf_M[complete.cases(repdf_M),]
  
  # save this
  model1_M <<- glm(didrep~species2+size+nonpure, data=repdf_M,family="binomial")  
  summary(model1_M)
  
  # save this
  model2_M <<- glm(didrep~species2*size+nonpure, data=repdf_M,family="binomial")
  summary(model2_M)
  
  # save this
  model3_M <<- glm(didrep~species2*size, data=repdf_M,family="binomial")   # best model
  summary(model3_M)
  
  # save this
  model4_M <<- glm(didrep~species2, data=repdf_M,family="binomial")
  summary(model4_M)
  
  # save this
  model5_M <<- glm(didrep~species2+size, data=repdf_M,family="binomial")
  summary(model5_M)
  
  AIC(model1_M,model2_M,model3_M,model4_M,model5_M)  # model 5 is best
  
  
  ## supplemental analysis- which females are more affected by reproductive interference?
  
  rep_onlypuremacF <- subset(repdf_F,species=="Macrotis")
  rep_onlypurefusF <- subset(repdf_F,species=="Fuscipes")
  rep_onlypureF <- rbind(rep_onlypuremacF,rep_onlypurefusF)
  
  model1_Fp <<- glm(didrep~species2*nonpure+size, data=rep_onlypureF,family="binomial")
  summary(model1_Fp)
  
  # look at repro success for female pure macrotis
  
  model1_Fm <<- glm(didrep~size+nonpure, data=rep_onlypuremacF,family="binomial")   # best model
  summary(model1_Fm)
  
  # look at repro success for female pure fuscipes
  
  model1_Ff <<- glm(didrep~size+nonpure, data=rep_onlypurefusF,family="binomial")   # best model
  summary(model1_Ff)
  
  
}


#################
# TODO: visualize probability of reproduction as function of species, size, and nonpure

### TODO: rename species


###########
# function: pdp for rep prob

plotpdp_repprob <- function(mod=model1_F,sex="Fp",var="species2",xlab="Fraction non-pure"){
  
  if(sex=="F"){
    dat = repdf_F
  }else if(sex=="M"){
    dat = repdf_M
  }else{
    dat = rep_onlypureF
    dat$species <- factor(dat$species,levels=sort(unique(dat$species)))
    dat$species2 <- factor(dat$species2,levels=sort(unique(dat$species2)))
  }
  
  if(is.factor(dat[[var]])){
    thisseq = levels(dat[[var]])
  }else{
    thisseq = seq(min(dat[[var]],na.rm=T),max(dat[[var]],na.rm=T),length=100)
  }
  
  newdat <- as.data.frame(thisseq)
  names(newdat)[1] <- var
  
  if(sex=="Fp"&var=="nonpure"){
    var2 <- c("nonpure","species2")
    thisseq1 = seq(min(dat[[var]],na.rm=T),max(dat[[var]],na.rm=T),length=100)
    thisseq2 = levels(dat$species)
    newdat2 <- expand.grid(thisseq2,thisseq1)
    names(newdat2) <- c("species2","nonpure") 
    othervars <- setdiff(rownames(attr(mod$terms,"factors")),c("didrep",var2) )
    v="size"
    for(v in othervars){
      if(is.factor(dat[[v]])){
        newdat2[[v]] <- names(which.max(table(dat[[v]])))[1]
      }else{
        newdat2[[v]] <- mean(dat[[v]],na.rm=T) 
      }
    }
    
    pred2 <- predict(mod,newdata = newdat2,type="response",se.fit=T)
  }
  
  #newdat$species2 <- spec
  
  othervars <- setdiff(rownames(attr(mod$terms,"factors")),c("didrep",var) )
  
  v="species2"
  for(v in othervars){
    if(is.factor(dat[[v]])){
      newdat[[v]] <- names(which.max(table(dat[[v]])))[1]
    }else{
      newdat[[v]] <- mean(dat[[v]],na.rm=T) 
    }
  }
  
  pred <- predict(mod,newdata = newdat,type="response",se.fit=T)
  
  if(sex=="Fp"&var=="nonpure"){
    plot(newdat[[var]],pred$fit,type="l",ylim=c(max(0,min(pred$fit-3*pred$se.fit)),min(1,max(pred$fit+3*pred$se.fit))),lwd=0,
         xlab=xlab,ylab="Probability of reproduction")
    # lines(newdat[[var]],pred$fit+2*pred$se.fit,lty=2)
    # lines(newdat[[var]],pmax(0,pred$fit-2*pred$se.fit),lty=2)
    ndxf <- newdat2$species2=="Fuscipes"
    ndxm <- newdat2$species2=="Macrotis"
    # speccol
    lines(newdat[[var]],pred2$fit[ndxf]+0*pred2$se.fit[ndxf],lty=1,lwd=2,col=speccol["Fuscipes"])
    polygon(c(newdat[[var]],rev(newdat[[var]])),c(pred2$fit[ndxf]-2*pred2$se.fit[ndxf],rev(pred2$fit[ndxf]+2*pred2$se.fit[ndxf])),
            col=makeTransparent(speccol["Fuscipes"],100),border=NA )
    lines(newdat[[var]],pred2$fit[ndxm]+0*pred2$se.fit[ndxm],lty=1,lwd=2,col=speccol["Macrotis"])
    polygon(c(newdat[[var]],rev(newdat[[var]])),c(pred2$fit[ndxm]-2*pred2$se.fit[ndxm],rev(pred2$fit[ndxm]+2*pred2$se.fit[ndxm])),
            col=makeTransparent(speccol["Macrotis"],100),border=NA )
    legend("bottomleft", fill=makeTransparent(speccol[c("Fuscipes","Macrotis")],100),border=c(NA,NA),
           col=speccol[c("Fuscipes","Macrotis")], bty="n",legend= c(expression(italic("N. fuscipes")), expression(italic("N. macrotis")) )   ) # lwd=c(2,2), lty=c(1,1),
    
  }else{
    if(is.factor(dat[[var]])){
      # plot(as.factor(newdat[[var]]),pred$fit,type="l",ylim=c(max(0,min(pred$fit-3*pred$se.fit)),min(1,max(pred$fit+3*pred$se.fit))),lwd=3,
      #      xlab=xlab,ylab="Prob rep.")
      Hmisc::errbar(as.numeric(as.factor(newdat[[var]])),
                    pred$fit,pred$fit+2*pred$se.fit,pred$fit-2*pred$se.fit,
                    xlab="",xaxt="n",ylab="Probability of reproduction",xlim=c(0.5,3.5))
      axis(1,at=1:3,las=3,labels=c(expression(italic("N. fuscipes")), expression("BC and F1"),expression(italic("N. macrotis"))))      # levels(dat[[var]])
    }else{
      plot(newdat[[var]],pred$fit,type="l",ylim=c(max(0,min(pred$fit-3*pred$se.fit)),min(1,max(pred$fit+3*pred$se.fit))),lwd=3,
           xlab=xlab,ylab="Probability of reproduction")
      lines(newdat[[var]],pred$fit+2*pred$se.fit,lty=2)
      lines(newdat[[var]],pmax(0,pred$fit-2*pred$se.fit),lty=2)
    }
  }
  
  
}

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#################
# FUNCTION: visualize prob repro model

plot_probrepro <- function(){
  
  svg("probrep2.svg",6.7,5)
  par(mfcol=c(2,2))
  thissex <- "M"
  for(thissex in c("F","M")){
    pures <- list()
    sizes <- list()
    if(thissex=="M"){
      thismod <- model2_M
      for(sp in c("Fuscipes","Hybrid","Macrotis")){
        temp <- subset(repdf_M,species2==sp)
        sizes[[sp]] <- quantile(temp$size,probs=c(0.05,0.5,0.95),na.rm=T)
        pures[[sp]] <- quantile(temp$nonpure,probs=c(0.05,0.95),na.rm=T)
      }
      
    }else{
      thismod <- model1_F
      for(sp in c("Fuscipes","Hybrid","Macrotis")){
        temp <- subset(repdf_F,species2==sp)
        sizes[[sp]] <- quantile(temp$size,probs=c(0.05,0.5,0.95),na.rm=T)
        pures[[sp]] <- quantile(temp$nonpure,probs=c(0.05,0.95),na.rm=T)
      }
    }
    
    
    newdat <- data.frame(
      species2=rep(c("Fuscipes","Hybrid","Macrotis"),each=3),
      size=c(sizes$Fuscipes,sizes$Hybrid,sizes$Macrotis),
      nonpure=rep(c(pures$Fuscipes[1],pures$Hybrid[1],pures$Macrotis[1]),each=3)
    )
    
    pred <- predict(thismod,newdata = newdat,type="response",se.fit=T)
    
    pred2 <- matrix(pred$fit,nrow=3,ncol=3)
    rownames(pred2) <- c("sm","md","lg")
    colnames(pred2) <- c("Fuscipes","Hybrid","Macrotis")
    x<-barplot(pred2,beside = T,space=c(0,1),ylim=c(0,1),legend=F)
    arrows(x,pred$fit-pred$se.fit,x,pred$fit+pred$se.fit,angle=90,code=3,length=0.04,main=thissex)
    title(sprintf("%s, low hybrid avail",thissex))
    
    newdat <- data.frame(
      species2=rep(c("Fuscipes","Hybrid","Macrotis"),each=3),
      size=c(sizes$Fuscipes,sizes$Hybrid,sizes$Macrotis),
      nonpure=rep(c(pures$Fuscipes[2],pures$Hybrid[2],pures$Macrotis[2]),each=3)
    )
    
    pred <- predict(thismod,newdata = newdat,type="response",se.fit=T)
    
    pred2 <- matrix(pred$fit,nrow=3,ncol=3)
    rownames(pred2) <- c("sm","md","lg")
    colnames(pred2) <- c("Fuscipes","Hybrid","Macrotis")
    x<-barplot(pred2,beside = T,space=c(0,1),ylim=c(0,1),legend=F)
    arrows(x,pred$fit-pred$se.fit,x,pred$fit+pred$se.fit,angle=90,code=3,length=0.04,main=thissex)
    title(sprintf("%s, high hybrid avail",thissex))
  }
  dev.off()
  
  
  png("probrep_pdp1.png",5.5,3.1,units = "in",res=600)
    par(mfrow=c(1,3),mai=c(1.7,0.7,0.1,0.1))
    plotpdp_repprob(mod=model1_F,sex="F",var="size",xlab="Body size (g)")
    fig_label("A",cex=1.5)
    plotpdp_repprob(mod=model1_F,sex="F",var="species2",xlab="Species")
    fig_label("B",cex=1.5)
    plotpdp_repprob(mod=model1_F,sex="F",var="nonpure",xlab="Fraction admixed mates\navailable")
    fig_label("C",cex=1.5)
  dev.off()
  
  png("probrep_pdp2.png",5.5,3.1,units = "in",res=600)
    par(mfrow=c(1,3),mai=c(1.7,0.7,0.1,0.1))
    plotpdp_repprob(mod=model1_F,sex="F",var="size",xlab="Body size (g)")
    fig_label("A",cex=1.5)
    plotpdp_repprob(mod=model1_F,sex="F",var="species2",xlab="Species")
    fig_label("B",cex=1.5)
    plotpdp_repprob(mod=model1_Fp,sex="Fp",var="nonpure",xlab="Fraction admixed mates\navailable")
    fig_label("C",cex=1.5)
  dev.off()
  
}


#####################
# FUNCTION: grab only known-reproductive individuals from the null-simulated dataset


thin_knownrep <- function(sex="F"){
  if(sex=="F"){
    exp.out <- exp.out.f
    allinds <- sort(unique(parent$mid_year))
    exp.out <- exp.out[exp.out$mid_year%in%allinds,]
  }else{
    exp.out <- exp.out.m
    allinds <- sort(unique(parent$pid_year))
    exp.out <- exp.out[exp.out$pid_year%in%allinds,]
  }
  return(exp.out)
}

#exp.out <- exp.out.f 



###################
# FUNCTION: prepare data for conditional logistic regression models

preparedata_forclr <- function(sex="F",spec="Fuscipes"){
  exp.out <- thin_knownrep(sex)
  
  #######
  # load hybrid zone survival data
  #######
  
  hzs <- read.csv("Data\\survival_hybrid_zone.csv")
  # head(hzs)
  # summary(hzs)
  
  #########
  # merge parent and exp.out to generate a master dataset for logistic regression type analysis
  mergeddf <- makemergeddf(parent_pairings,"all",sex)
  nrow(mergeddf)
  names(mergeddf)
  
  if(sex=="F"){
    idvar <- "mid_year"
  }else{
    idvar <- "pid_year"
  }
  
  
  #########
  # run conditional logistic regression: covariates are size differential and q value differential. 
  # build separate model for each species/type?
  
  ### thin for testing
  
  nper <- 25    # random observations per stratum
  allstrata <- sort(unique(mergeddf[[idvar]]))
  
  df2 <- NULL
  
  i=50
  for(i in 1:length(allstrata)){
    temp <- mergeddf[mergeddf[[idvar]]==allstrata[i],]   #subset(mergeddf,mid_year==allstrata[i])
    temp_repdf <- subset(repdf,id_year==allstrata[i])
    if(nrow(temp_repdf)>0){
      temp$Avail_fu <- temp_repdf$fuavl[1]
      temp$Avail_hy <- temp_repdf$nonpure[1]
      temp$Avail_ma <- temp_repdf$maavl[1]
    } else{
      temp$Avail_fu <- NA
      temp$Avail_hy <- NA
      temp$Avail_ma <- NA
    }
    true <- subset(temp,truemat==1)
    rand <- subset(temp,truemat==0)
    rand <- rand[1:nper,]
    if(nrow(true)>0){
      new <- rbind(true,rand)
      df2 <- rbind(df2,new) 
    }
  }
  
  # nrow(df2)    # 10760 observations
  
  
  ###########
  # determine the body size differential 
  ###########
  
  ## make sure this is correct- talk with E and M!
  
  allmoms <- sort(unique(df2$MotherID))
  alldads <- sort(unique(df2$FatherID))
  
  df2$dadsize <- NA
  df2$momsize <- NA
  i=1
  for(i in 1:length(allmoms)){
    thismom <- allmoms[i]
    momdat <- subset(rat,PJM.ID==thismom)
    temp <- subset(df2,MotherID==thismom)
    allyears <- sort(unique(temp$Year))
    j=allyears[1]
    for(j in allyears){
      temp2 <- subset(momdat,Year==j)
      thiswt <- mean(as.numeric(temp2$WtMax_g),na.rm=T)
      ndx <- which((df2$MotherID==thismom)&(df2$Year==j))
      df2$momsize[ndx] <- thiswt
    }
  }
  
  i=1
  for(i in 1:length(alldads)){
    thisdad <- alldads[i]
    daddat <- subset(rat,PJM.ID==thisdad)
    temp <- subset(df2,FatherID==thisdad)
    allyears <- sort(unique(temp$Year))
    j=allyears[1]
    for(j in allyears){
      temp2 <- subset(daddat,Year==j)
      thiswt <- mean(as.numeric(temp2$WtMax_g),na.rm=T)
      ndx <- which((df2$FatherID==thisdad)&(df2$Year==j))
      df2$dadsize[ndx] <- thiswt
    }
  }
  
  # hist(df2$momsize)
  # hist(df2$dadsize)
  
  df2$sizediff <- df2$dadsize-df2$momsize 
  
  # hist(df2$sizediff)
  # hist(df2$sizediff[df2$truemat==1],freq=F,col="red",density=5)
  # hist(df2$sizediff[df2$truemat==0],freq=F,add=T)
  
  
  # summary(df2)
  
  #########
  # correlation between weight and q-value
  
  # names(df2)
  # plot(qFu.dad~dadsize,data=df2)
  # cor(df2[,c("qFu.dad","dadsize")],use="complete.obs")
  # 
  # plot(qFu.mom~momsize,data=df2)
  # cor(df2[,c("qFu.mom","momsize")],use="complete.obs")
  
  
  ########
  # Add hybrid survival
  
  names(df2)
  
  hzs2 <- tidyr::pivot_wider(hzs,names_from=genotype,values_from=median_survival)
  
  ndx <- match(df2$Year,hzs2$year)
  df2$hzs <- hzs2$hybrid[ndx]
  df2$hmf <- hzs2$hybrid[ndx]-hzs2$fuscipes[ndx]
  df2$fsurv <- hzs2$fuscipes[ndx]
  
  
  if(sex=="F"){
    df3 <-  subset(df2,Mom.spec==spec)
  }else{
    df3 <-  subset(df2,Dad.spec==spec)
  }
  
  nrow(df3)
  
  
  ## correlation between size difference and q difference
  
  # cor.test(df3$diff,df3$sizediff)   # negative correlation for fuscipes
  # 
  # plot(diff~sizediff,data=df3)
  # 
  # plot(df3$sizediff~df3$dadsize)
  # cor(df3$sizediff,df3$dadsize,use="complete.obs") # 0.87 correlation for fuscipes
  
  # Prepare covariates
  df3$qdiff_s <- scale(df3$diff)
  df3$sizediff_s <- scale(df3$sizediff)
  if(sex=="F"){
    df3$stratum <- df3$mid_year
    df3$qFu.mate_s <- scale(df3$qFu.dad) 
    df3$matesize <- df3$dadsize
    df3$matesize_s <- scale(df3$dadsize)
    df3$qFu.mate <- df3$qFu.dad
  }else{
    df3$stratum <- df3$pid_year
    df3$qFu.mate <- df3$qFu.mom
    df3$qFu.mate_s <- scale(df3$qFu.mom) 
    df3$matesize_s <- scale(df3$momsize)
    df3$matesize <- df3$momsize
  }
  
  df3$matepair <- as.factor(df3$truemat)   # response variable
  
  df3$hzs_s <- scale(df3$hzs)
  df3$hmf_s <- scale(df3$hzs)
  df3$fsurv_s <- scale(df3$fsurv)
  
  ######
  # Remove any NA observations
  
  
  ndx <- complete.cases(df3[,c("qdiff_s","sizediff_s","stratum","qFu.mate_s","matesize_s")])
  
  df3 <- df3[ndx,]
  nrow(df3)             # 11799
  
  df3$weights <- ifelse(df3$truemat==1,1,1000)
  df3$Year2 <- as.numeric(df3$Year)-min(as.numeric(df3$Year))+1 
  return(df3)
}






#####################
# FUNCTION: fit conditional logistic regression models
#####################


fit_mateselection <- function(data=df3){
  
  allmods <- list()
  
  # allmods[[1]] = glmmTMB(matepair ~ qdiff_s + sizediff_s +
  #                  (1|stratum),
  #                family=binomial,
  #                weights=weights,
  #                data=data
  # )
  
  # summary(allmods[[1]])
  # AIC(allmods[[1]])
  
 allmods[[1]] = glmmTMB(matepair ~ qFu.mate_s + matesize_s +     # BEST BASE MODEL FOR FUSCIPES  [NOW ALSO BEST FOR MACROTIS]
                     (1|stratum),
                   family=binomial,
                   weights=weights,
                   data=data
  )
  # summary (allmods[[1]])
  
  
  allmods[[2]] = glmmTMB(matepair ~ qFu.mate_s + poly(matesize_s,2) +      # BEST MODEL FOR MACROTIS  [NOT AFTER REMOVING DUPLICATED OBSERVATIONS IN PARENT DATASET]
                     (1|stratum),
                   family=binomial,
                   weights=weights,
                   data=data
  )
  # summary(allmods[[2]])
  
  allmods[[3]] = glmmTMB(matepair ~ qFu.mate_s + matesize_s + qFu.mate_s:matesize_s +    # BEST BASE MODEL FOR FUSCIPES  [NOW ALSO BEST FOR MACROTIS]
                           (1|stratum),
                         family=binomial,
                         weights=weights,
                         data=data
  )
  
  allmods[[4]] = glmmTMB(matepair ~ qFu.mate_s + poly(matesize_s,2) + qFu.mate_s:matesize_s +     # BEST MODEL FOR MACROTIS  [NOT AFTER REMOVING DUPLICATED OBSERVATIONS IN PARENT DATASET]
                           (1|stratum),
                         family=binomial,
                         weights=weights,
                         data=data
  )
  #summary(allmods[[4]])
  
  ## KTS: NOTE: what might explain the fact that the quadratic effect disappeared for macrotis after removing duplicates?? 
  # more offspring produced with smaller size disparity?
  
  
  # allmods[[3]] = glmmTMB(matepair ~ qFu.mate_s + matesize_s + qFu.mate_s:hzs_s +     # BEST MODEL FOR FUSCIPES? [NOT AFTER REMOVING DUPLICATED OBSERVATIONS IN PARENT DATASET]
  #                    (1|stratum),
  #                  family=binomial,
  #                  weights=weights,
  #                  data=data
  # )
  # summary(allmods[[3]])
  
  
  ## KTS: NOTE: what might explain the fact that the survival effect disappeared for fuscipes after removing duplicates?? 
  
  
                  
  # allmods[[4]] = glmmTMB(matepair ~ qFu.mate_s + poly(sizediff_s,2) + qFu.mate_s:Year2 +        # NOT USEFUL MODEL 
  #                    (1|stratum),
  #                  family=binomial,
  #                  weights=weights,
  #                  data=data
  # )
  # summary(allmods[[4]])
  
  
  allmods[[5]] = glmmTMB(matepair ~ qFu.mate_s + matesize_s + qFu.mate_s:Year2 +        # NOT USEFUL MODEL
                     (1|stratum),
                   family=binomial,
                   weights=weights,
                   data=data
  )
  
  aics <- sapply(allmods,AIC)
  which.min(aics)
  bestmod <- allmods[[which.min(aics)]]
  summary(bestmod)
  
  return(bestmod)
  

}


################
# FUNCTION: visualize mate pairing model results

visualize_mateselection <- function(model=bestmod_FUfem,data=df_FUfem,species="Fuscipes",sex="F"){
  # if(thisspec=="Macrotis"){   # model 1.1 best for fuscipes. mod1.1 best for macrotis
  #   model=mod1.1
  #   allvars = c("qFu.dad_s","dadsize_s") #"sizediff_s")
  # }else{
  #   model=mod1.1
  #   allvars = c("qFu.dad_s","dadsize_s") #,"hzs_s")
  # } 
  allvars_raw <- c("qFu.mate_s","matesize_s","hzs_s","Year2")
 
  response <- "matepair"
  #othervars <- c("stratum","(weights)")
  
  keep <- c()
  n=1
  for(n in 1:length(allvars_raw)){
    thisvar <- allvars_raw[n]
    if(any(grepl(allvars_raw[n], names(model$frame)))) keep = c(keep,n)
  }
  
  allvars <- allvars_raw[keep]
  
  svg(sprintf("MateSelectionFig_%s_%s.svg",species,sex),5,3)
  par(mfrow=c(1,2))
  
  predvars <- c("matesize_s","qFu.mate_s")
  varnames <- c("Mate Size","Mate genotype (frac FU)")
  for(p in 1:length(predvars)){
    predvar= predvars[p]  #"matesize_s" # "qFu.dad_s" # "qFu.dad_s" # "sizediff_s"
    varname= varnames[p]  #"Mate Size" # "Dad Genotype (qFu)" #  "Dad Genotype (qFu)" # "Dad Size" #"Size Diff" #    # NOTE: make size difference reversed (dad usually bigger)
    
    len <- 25
    
    dim <- data[,predvar]
    range <- seq(min(dim),max(dim),length=len)
    predvar2 <- gsub("_s","",predvar)
    realmean <- mean(data[[predvar2]])
    realsd <- sd(data[[predvar2]])
    newdata <- data.frame(temp=range)
    names(newdata) <- c(predvar)
    othervars <- allvars[!allvars%in%c(predvar,response)]
    
    for(var in othervars){
      newdata[,var] <- 0
    }
    
    # if(predvar=="sizediff_s") newdata$qFu.dad_s=min(df3$qFu.dad_s)
    
    newdata$stratum <- NA  #data$ID[1]
    newdata$weights = 1000
    
    pred <- predict(model,newdata,type="response",se.fit=T)   # ,type="conditional",se.fit=F
    
    col="black"
    # if(predvar=="sizediff_s") col = speccol[5]
    
    
    plot(range,pred$fit,xlab=varname,ylab="Relative Pairing Success",type="l",lwd=2,xaxt="n",col=col,
         ylim=c(0,max(pred$fit)*2),yaxt = "n")
    points(range,(pred$fit+pred$se.fit*2),type="l",lty=2,col=col)
    points(range,(pred$fit-pred$se.fit*2),type="l",lty=2,col=col)
    ats <- seq(min(range),max(range),length=6)
    axis(1,ats,labels = round(realmean+ats*realsd,2))
    rug(jitter(data[seq(1,nrow(data),3),][[predvar]]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
    ats <- seq(min(pred$fit),max(pred$fit)*1.5,length=4)
    axis(2,ats,labels=round(ats*1000,2))
  }
  dev.off()
    
  # 
  # if(predvar=="sizediff_s"){
  #   newdata$qFu.dad_s=max(df3$qFu.dad_s)
  #   pred <- predict(model,newdata,type="response",se.fit=T)   # ,type="conditional",se.fit=F
  #   lines(range,pred$fit,lwd=2,col=speccol[1])
  #   lines(range,(pred$fit+pred$se.fit*2),lty=2,col=speccol[1])
  #   lines(range,(pred$fit-pred$se.fit*2),lty=2,col=speccol[1])
  #   legend("topright",bty="n",col=c(speccol[5],speccol[1]),lwd=c(2,2),legend=c("Macrotis","Fuscipes"))
  # }
  
  
  ##################
  #### Plot interaction between hybrid zone survival and male q value
  
  
  # # c("qFu.dad_s","dadsize_s","hzs_s")
  # 
  # len=15
  # var1 = "qFu.dad_s"
  # var2 = "hzs_s"
  # dim1 <- df3[,var1]
  # dim2 <- df3[,var2]
  # range1 <- seq(min(dim1),max(dim1),length=len)
  # range2 <- seq(min(dim2),max(dim2),length=len)
  # var1_2 <- gsub("_s","",var1)
  # var2_2 <- gsub("_s","",var2)
  # realmean1 <- mean(df3[[var1_2]])
  # realmean2 <- mean(df3[[var2_2]])
  # realsd1 <- sd(df3[[var1_2]])
  # realsd2 <- sd(df3[[var2_2]])
  # 
  # nd = expand.grid(range1,range2)
  # names(nd) = c(var1,var2)
  # 
  # othervars <- allvars[!allvars%in%c(var1,var2,"matepair")]
  # 
  # for(var in othervars){
  #   nd[,var] <- 0
  # }
  # 
  # nd$stratum <- NA  #data$ID[1]
  # nd$weights = 1000
  # 
  # pred <- predict(mod2.1,nd,type="response",se.fit=F)   # ,type="conditional",se.fit=F
  # pred <- pred*1000
  # 
  # predmat = matrix(pred,nrow=length(range1),ncol=length(range2))
  # 
  # svg("3dplot_1.svg",5,5)
  # persp(range1,range2,predmat,theta=-25,phi=15,xlab=var1,ylab=var2,zlim=c(0,max(pred)))
  # dev.off()
  # 
  # 
  # 
  # ##################
  # ####
  # # Plot size differences (by species and by sex)
  # 
  # 
  # # boxplot figure for marjorie
  # par(mfrow=c(1,2))
  # for(i in 1:2){
  #   thisspec <- Species[i]      # "Fuscipes" "Macrotis" "Hybrid" "Fus_Back" "Mac_Back"
  #   df3 <-  subset(df2,Mom.spec==thisspec)
  #   
  #   thisspec <- Species[i]
  #   df4 <- subset(df3,Dad.spec==thisspec)
  #   
  #   sub1 <- subset(df4,truemat==1)
  #   sub0 <- subset(df4,truemat==0)
  #   
  #   sub2 <- rbind(sub1,sub0)
  #   
  #   plot(sub2$dadsize~factor(sub2$truemat,levels=c(0,1),labels=c("No","Yes")),main=thisspec,ylab="Dad weight", xlab="Successful pairing")
  #   
  # }
  
}



#####################
# FUNCTION: make summary table
#####################


make_summarytable1 <- function(){
  rat2 <- subset(rat,Year!=2013)
  rat2$Species <- factor(rat2$Species,levels=Species)
  temp1 <- subset(rat2,SexFI=="M"&AgeFI=="A")
  allmales <- sapply(Species,function(t){temp2 <- subset(temp1,Species==t);length(unique(temp2$PJM.ID))}   )
  temp1 <- subset(rat2,SexFI=="F"&AgeFI=="A")
  allfems <- sapply(Species,function(t){temp2 <- subset(temp1,Species==t);length(unique(temp2$PJM.ID))}   )
  parent2 <- parent
  parent2$Dad.spec <- factor(parent$Dad.spec,levels=Species)
  parent2$Mom.spec <- factor(parent$Mom.spec,levels=Species)
  parent_pairings2 <- parent_pairings
  parent_pairings2$Dad.spec <- factor(parent_pairings$Dad.spec,levels=Species)
  parent_pairings2$Mom.spec <- factor(parent_pairings$Mom.spec,levels=Species)
  tot_offsp_m <- table(parent2$Dad.spec)  # NOTE: there are 15 instances of a hybrid father siring an offspring and 51 instances of a hybrid mother
  tot_offsp_f <- table(parent2$Mom.spec)
  tot_pairings_m  <- table(parent_pairings2$Dad.spec)
  tot_pairings_f <- table(parent_pairings2$Mom.spec)
  tot_repind_m <- sapply(Species,function(t){temp <- subset(parent, Dad.spec==t & !is.na(FatherID)); length(unique(temp$FatherID))} )
  tot_repind_f <- sapply(Species,function(t){temp <- subset(parent, Mom.spec==t & !is.na(MotherID)); length(unique(temp$MotherID))} )
  towrite <- rbind(
    allmales,tot_offsp_m,tot_pairings_m,tot_repind_m,allfems,tot_offsp_f,tot_pairings_f,tot_repind_f
  )
  #rownames
  write.csv(towrite,"RawSummaryTable.csv",row.names = T)
}


####################
# FUNCTION: make elizabeth box plots



makeELIZboxplots <- function(sex="F", year='all',onlyexp = FALSE){
  
  ######
  #FIGURE FOR RAPID:
  
  lab3 <- ifelse(onlyexp,"background","withobs")
  filename=sprintf("ElizBoxPlot_%s_%s_%s.pdf",sex,year,lab3) 
  pdf(filename, width=6, height=4)
  par(mfrow=c(1,1), mai=c(1,3,0.1,0.1), las=1)
  
  if(year=="all"){
    years <- allyears
  }else{
    year <- year
  }
  
  if(sex=="F"){
    exp.out <- exp.out.f[exp.out.f$Year%in%years,]
    exp.out <- subset(exp.out,Simulation%in%(1:5))
    exp.out$Mom.num <- ifelse(exp.out$Mom.spec=="Macrotis", 0, ifelse(exp.out$Mom.spec=="Mac_Back", 0.25, ifelse(exp.out$Mom.spec=="Hybrid", 0.5, ifelse(exp.out$Mom.spec=="Fus_Back", 0.75, 1))))
    exp.out$Dad.num <- ifelse(exp.out$Dad.spec=="Macrotis", 0, ifelse(exp.out$Dad.spec=="Mac_Back", 0.25, ifelse(exp.out$Dad.spec=="Hybrid", 0.5, ifelse(exp.out$Dad.spec=="Fus_Back", 0.75, 1))))
    
    parent2 <- parent_pairings[parent_pairings$Year%in%years,]
    parent2$Mom.num <- ifelse(parent2$Mom.spec=="Macrotis", 0, ifelse(parent2$Mom.spec=="Mac_Back", 0.25, ifelse(parent2$Mom.spec=="Hybrid", 0.5, ifelse(parent2$Mom.spec=="Fus_Back", 0.75, 1))))
    parent2$Dad.num <- ifelse(parent2$Dad.spec=="Macrotis", 0, ifelse(parent2$Dad.spec=="Mac_Back", 0.25, ifelse(parent2$Dad.spec=="Hybrid", 0.5, ifelse(parent2$Dad.spec=="Fus_Back", 0.75, 1))))
    
    plot(jitter(exp.out$qFu.dad) ~ jitter(exp.out$qFu.mom), pch=20, cex=0.1, xlab="Mother genotype", 
         ylab="", xaxt="n", yaxt="n",col=rgb(0,0,0,0.5))
    axis(side=1, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("macrotis", "mac-back", "F1", "fus-back", "fuscipes"), cex=0.8)
    axis(side=2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c("macrotis", "mac-back", "F1", "fus-back", "fuscipes"), cex=0.8)
    title(ylab="Father genotype",line=6)
    
    if(!onlyexp) points(parent2$qFu.dad ~ jitter(parent2$qFu.mom), col=rgb(1,0,0,0.7), pch=16)
  }else{
    exp.out <- exp.out.m[exp.out.m$Year%in%years,]
    parent2 <- parent_pairings[parent_pairings$Year%in%years,]
    plot(exp.out$qFu.mom ~ jitter(exp.out$qFu.dad), pch=16, xlab="Male qFu", ylab="Female qFu")
    if(!onlyexp) points(parent2$qFu.mom ~ jitter(parent2$qFu.dad), col="red", pch=16)
  }
 
  dev.off()
}

