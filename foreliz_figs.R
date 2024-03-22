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




##############
# figure panels: change in admixed/pure mate availability over time

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



############
# body size figure (not sure exactly how marjorie wants this information presented yet...)

svg("specbysize_adults.svg",5.5,4)
#repdf$didrep
par(mfrow=c(1,2))     # look at sizes by species
plot(repdf_M$size~repdf_M$species,ylab="size",xlab="",main="male",las=2,ylim=c(0,500))
plot(repdf_F$size~repdf_F$species,ylab="size",xlab="",main="female",las=2,ylim=c(0,500))
dev.off()



############
# repro prob model

###########
# function: pdp for rep prob

plotpdp_repprob <- function(mod=model1_F,sex="F",var="species2",xlab="Fraction non-pure"){
  
  if(sex=="F"){
    dat = repdf_F
  }else{
    dat = repdf_M
  }
  
  if(is.factor(dat[[var]])){
    thisseq = levels(dat[[var]])
  }else{
    thisseq = seq(min(dat[[var]],na.rm=T),max(dat[[var]],na.rm=T),length=100)
  }
  
  
  newdat <- as.data.frame(thisseq)
  names(newdat)[1] <- var
  
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


svg("probrep_pdp1.svg",5.5,3.1)
par(mfrow=c(1,3),mai=c(1.7,0.7,0.1,0.1))
plotpdp_repprob(mod=model1_F,sex="F",var="size",xlab="Body size (g)")
fig_label("A",cex=1.5)
plotpdp_repprob(mod=model1_F,sex="F",var="species2",xlab="Species")
fig_label("B",cex=1.5)
plotpdp_repprob(mod=model1_F,sex="F",var="nonpure",xlab="Fraction admixed mates\navailable")
fig_label("C",cex=1.5)
dev.off()


#############
# package data for E

# save(repdf,allyears,repdf_M,repdf_F,file="dataforE.RData")

# save(model1_F,repdf_F,repdf_M, file="dataforE2.RData")
















