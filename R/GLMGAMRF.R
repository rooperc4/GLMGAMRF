#' A function to predict CPUE-index from a delta-lognormal GLM with delta method errors
#'
#' This function uses a binomial GLM and a log-normal GLM to compute a delta-lognormal
#' index of abundance with associated 95% confidence intervals estimated by the delta
#' method. A dummy variable for year must be included to get the annual abundance
#' index. The function back-transforms the index and the CI's.
#' @param pa.model binomial GLM predicting presence or absence
#' @param cpue.model log-normal GLM predicting abundance where catch is positive
#' @keywords delta-lognormal glm, survey abundance index, delta method
#' @export
#' @examples
#' glm.pa.xvars<-c("inverts","slope","btemp","bdepth")
#' glm.pa.yvar<-ifelse(PA.data[species.name]>0,1,0)
#' glm.pa.form <- as.formula(paste("glm.pa.yvar ~", paste(glm.pa.xvars,collapse="+"),"+as.factor(year)",sep=""))
#' pa.glm <- glm(glm.pa.form, family = binomial, data = PA.data)
#' 
#' CPUE.data<-subset(PA.data,PA.data[species.name]>0)
#' glm.cpue.yvar<-unlist(log(CPUE.data[species.name]))
#' glm.cpue.xvars<-c("inverts","slope","btemp","bdepth")
#' glm.cpue.form <- as.formula(paste("glm.cpue.yvar ~", paste(glm.cpue.xvars,collapse="+"),"+as.factor(year)",sep=""))
#' cpue.glm <- glm(glm.cpue.form, family = gaussian, data = CPUE.data)
#' 
#' predict.glm.index(pa.glm,cpue.glm)
#' 


###Delta-lognormal GLM with Delta method errors  
predict.glm.index<-function(pa.model,cpue.model){
  pa.data<-pa.model$data
  cpue.data<-cpue.model$data
  pa.meds<-t(data.frame(apply(pa.data,MARGIN=2,FUN=median)))
  pa.meds.names<-colnames(pa.meds)
  yrs<-unique(pa.data$year)
  
  pa.meds<-data.frame(rep.row(pa.meds,length(yrs)))
  colnames(pa.meds)<-pa.meds.names
  pa.meds$year<-yrs
  
  pa.predict<-predict(pa.model,newdata=pa.meds,type="response",se.fit=TRUE)
  cpue.predict<-predict(cpue.model,newdata=pa.meds,type="response",se.fit=TRUE)
  glm.raw.index<-pa.predict$fit*cpue.predict$fit

d<-aggregate(cpue.model$fitted.values,by=list(cpue.model$data$year),FUN=mean)
P<-aggregate(pa.model$fitted.values,by=list(pa.model$data$year),FUN=mean)
Vard<-cpue.predict$se.fit^2
VarP<-pa.predict$se.fit^2

pa.obsyr<-aggregate(pa.model$y,by=list(pa.data$year),FUN=mean)
cpue.obsyr<-aggregate(cpue.model$y,by=list(cpue.data$year),FUN=mean)

###COV term using observations
#rho<-cor(pa.obsyr$x,cpue.obsyr$x)
#SEP<-aggregate(pa.model$y,by=list(pa.data$year),FUN=sd)
#SEd<-aggregate(cpue.model$y,by=list(cpue.data$year),FUN=sd)
#SEPn<-aggregate(pa.model$y,by=list(pa.data$year),FUN=length)
#SEdn<-aggregate(cpue.model$y,by=list(cpue.data$year),FUN=length)
#SEP<-SEP$x/sqrt(SEPn$x)
#SEd<-SEd$x/sqrt(SEdn$x)
#Covdp<-cov(pa.obsyr$x,cpue.obsyr$x)

###COV term using predictions
#rho<-cor(pa.predict$fit,cpue.predict$fit)
#SEP<-aggregate(pa.model$fitted.values,by=list(pa.model$data$year),FUN=sd)
#SEd<-aggregate(cpue.model$fitted.values,by=list(cpue.model$data$year),FUN=sd)
#SEPn<-aggregate(pa.model$fitted.values,by=list(pa.model$data$year),FUN=length)
#SEdn<-aggregate(cpue.model$fitted.values,by=list(cpue.model$data$year),FUN=length)
#SEP<-SEP$x/sqrt(SEPn$x)
#SEd<-SEd$x/sqrt(SEdn$x)
#Covdp<-rho*(SEd*SEP) # From Lo 1992
#Covdp<-cov(pa.predict$fit,cpue.predict$fit) #Using cov function

###COV term using dependent variables
Covdp<-cov(pa.predict$fit^2,cpue.predict$fit^2)+mean(pa.predict$fit^2)*mean(cpue.predict$fit^2)-
  (cov(pa.predict$fit,cpue.predict$fit)+mean(pa.predict$fit)*mean(cpue.predict$fit))^2
#Covdp<-cov(d$x^2,P$x^2)+P$x^2*d$x^2-(cov(P$x,d$x)+(P$x*d$x))^2

VardP2<-Vard*P$x^2
VarPd2<-VarP*d$x^2
#Vari<-VardP2+VarPd2+2*d$x*P$x*Covdp
Vari<-VardP2+VarPd2+Covdp
  glm.index.sd<-sqrt(Vari)  

glm.index<-exp(glm.raw.index) 
upper_CI<-exp(glm.raw.index+1.96*glm.index.sd)
lower_CI<-exp(glm.raw.index-1.96*glm.index.sd)

  return(data.frame(years=yrs,glm.index,lower_CI=lower_CI,upper_CI=upper_CI))
}

#' A function to predict CPUE-index from a delta-lognormal GLM with bootstrapped errors
#'
#' This function uses a binomial GLM and a log-normal GLM to compute a delta-lognormal
#' index of abundance with associated 95% confidence intervals estimated by the boot-
#' strapping the data. A dummy variable for year must be included to get the annual abundance
#' index. The function back-transforms the index and the CI's.
#' @param pa.model binomial GLM predicting presence or absence
#' @param cpue.model log-normal GLM predicting abundance where catch is positive
#' @param boot_reps number of bootstrap replications (default = 500)
#' @keywords delta-lognormal glm, survey abundance index, bootstrap
#' @export
#' @examples
#' glm.pa.xvars<-c("inverts","slope","btemp","bdepth")
#' glm.pa.yvar<-ifelse(PA.data[species.name]>0,1,0)
#' glm.pa.form <- as.formula(paste("glm.pa.yvar ~", paste(glm.pa.xvars,collapse="+"),"+as.factor(year)",sep=""))
#' pa.glm <- glm(glm.pa.form, family = binomial, data = PA.data)
#' 
#' CPUE.data<-subset(PA.data,PA.data[species.name]>0)
#' glm.cpue.yvar<-unlist(log(CPUE.data[species.name]))
#' glm.cpue.xvars<-c("inverts","slope","btemp","bdepth")
#' glm.cpue.form <- as.formula(paste("glm.cpue.yvar ~", paste(glm.cpue.xvars,collapse="+"),"+as.factor(year)",sep=""))
#' cpue.glm <- glm(glm.cpue.form, family = gaussian, data = CPUE.data)
#' 
#' predict.glm.bindex(pa.glm,cpue.glm)
#' 
#######bootstrapping method

predict.glm.bindex<-function(pa.model,cpue.model,boot_reps=500){
  pa.data<-pa.model$data
  pa.meds<-t(data.frame(apply(pa.data,MARGIN=2,FUN=median)))
  pa.meds.names<-colnames(pa.meds)
  yrs<-unique(pa.data$year)

  pa.meds<-data.frame(rep.row(pa.meds,length(yrs)))
  colnames(pa.meds)<-pa.meds.names
  pa.meds$year<-yrs
  
  pa.predict<-predict(pa.model,newdata=pa.meds,type="response",se.fit=TRUE)
  cpue.predict<-predict(cpue.model,newdata=pa.meds,type="response",se.fit=TRUE)
  glm.raw.index<-pa.predict$fit*cpue.predict$fit
  glm.index<-exp(glm.raw.index)

index_ests<-data.frame(array(0,dim=c(0,2)))
colnames(index_ests)<-c("yrs","temp.index")
for(i in 1:boot_reps){
  bootdata1<- sample(1:length(pa.model$data[,1]), replace=TRUE)
  bootdata<-data.frame(pa.model$data[bootdata1,])
  glm.pa.yvar<-pa.model$y[bootdata1]
  boot.pa<-glm(pa.model$formula,family=binomial,data=bootdata)
  
  bootdata1<-sample(1:length(cpue.model$data[,1]), replace=TRUE)
  bootdata<-data.frame(cpue.model$data[bootdata1,])
  glm.cpue.yvar<-cpue.model$y[bootdata1]
  boot.cpue<-glm(cpue.glm$formula,family=gaussian,data=bootdata)
  
  #temp1<-predict.glm.bindex(boot.pa,boot.cpue)
  pa.pboot<-predict(boot.pa,newdata=pa.meds,type="response",se.fit=TRUE)
  cpue.pboot<-predict(boot.cpue,newdata=pa.meds,type="response",se.fit=TRUE)
  temp.raw.index<-pa.pboot$fit*cpue.pboot$fit
  temp.index<-exp(temp.raw.index)  
  temp.index<-cbind(yrs,temp.index)
  
  index_ests<-rbind(index_ests,temp.index)
   }

index_est_sd<-aggregate(index_ests$temp.index,by=list(index_ests$yrs),FUN=sd)
upper_bootCI<-glm.index+1.96*index_est_sd$x#/sqrt(boot_reps)
lower_bootCI<-glm.index-1.96*index_est_sd$x#/sqrt(boot_reps)

return(data.frame(years=yrs,glm.index,lower_bootCI=lower_bootCI,upper_bootCI=upper_bootCI))
}
############################################################################
########repeat row-column functions
#' Function to repeat the columns of a matrix
#'
#' This function repeats the columns of a matrix, x for n times
#' @param x matrix or data frame
#' @param n number of times to repeat the columns of the matrix x
#' @keywords repeat rows
#' @export
#' @examples
#' 
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#' Function to repeat the rows of a matrix
#'
#' This function repeats the rows of a matrix, x for n times
#' @param x matrix or data frame
#' @param n number of times to repeat the rows of the matrix x
#' @keywords repeat rows
#' @export
#' @examples
#' 
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#######################################################################################
#######################################################################################
#######################################################################################


#' A function to predict CPUE-index from a delta-lognormal GAM with delta method errors
#'
#' This function uses a binomial GAM and a log-normal GAM to compute a delta-lognormal
#' index of abundance with associated 95% confidence intervals estimated by the delta
#' method. A dummy variable for year must be included to get the annual abundance
#' index. The function back-transforms the index and the CI's.
#' @param pa.model binomial GAM predicting presence or absence
#' @param cpue.model log-normal GAM predicting abundance where catch is positive
#' @keywords delta-lognormal gam, survey abundance index, delta method
#' @export
#' @examples
#' gam.pa.xvars<-c("inverts","slope","btemp","bdepth")
#' gam.pa.yvar<-ifelse(PA.data[species.name]>0,1,0)
#' gam.pa.form <- as.formula(paste("gam.pa.yvar ~", paste(gam.pa.xvars,collapse="+"),"+as.factor(year)",sep=""))
#' pa.gam <- gam(gam.pa.form, family = binomial, data = PA.data)
#' 
#' CPUE.data<-subset(PA.data,PA.data[species.name]>0)
#' gam.cpue.yvar<-unlist(log(CPUE.data[species.name]))
#' gam.cpue.xvars<-c("inverts","slope","btemp","bdepth")
#' gam.cpue.form <- as.formula(paste("gam.cpue.yvar ~", paste(gam.cpue.xvars,collapse="+"),"+as.factor(year)",sep=""))
#' cpue.gam <- gam(gam.cpue.form, family = gaussian, data = CPUE.data)
#' 
#' predict.gam.index(pa.gam,cpue.gam)
#' 


###Delta-lognormal GAM with Delta method errors  
predict.gam.index<-function(pa.model,cpue.model){
  pa.data<-pa.model$data
  cpue.data<-cpue.model$data
  pa.meds<-t(data.frame(apply(pa.data,MARGIN=2,FUN=median)))
  pa.meds.names<-colnames(pa.meds)
  yrs<-unique(pa.data$year)
  
  pa.meds<-data.frame(rep.row(pa.meds,length(yrs)))
  colnames(pa.meds)<-pa.meds.names
  pa.meds$year<-yrs
  
  pa.predict<-predict(pa.model,newdata=pa.meds,type="response",se.fit=TRUE)
  cpue.predict<-predict(cpue.model,newdata=pa.meds,type="response",se.fit=TRUE)
  gam.raw.index<-pa.predict$fit*cpue.predict$fit
  
  d<-aggregate(cpue.model$fitted.values,by=list(cpue.model$data$year),FUN=mean)
  P<-aggregate(pa.model$fitted.values,by=list(pa.model$data$year),FUN=mean)
  Vard<-cpue.predict$se.fit^2
  VarP<-pa.predict$se.fit^2
  
  pa.obsyr<-aggregate(pa.model$y,by=list(pa.data$year),FUN=mean)
  cpue.obsyr<-aggregate(cpue.model$y,by=list(cpue.data$year),FUN=mean)
  
  ###COV term using observations
  #rho<-cor(pa.obsyr$x,cpue.obsyr$x)
  #SEP<-aggregate(pa.model$y,by=list(pa.data$year),FUN=sd)
  #SEd<-aggregate(cpue.model$y,by=list(cpue.data$year),FUN=sd)
  #SEPn<-aggregate(pa.model$y,by=list(pa.data$year),FUN=length)
  #SEdn<-aggregate(cpue.model$y,by=list(cpue.data$year),FUN=length)
  #SEP<-SEP$x/sqrt(SEPn$x)
  #SEd<-SEd$x/sqrt(SEdn$x)
  #Covdp<-cov(pa.obsyr$x,cpue.obsyr$x)
  
  ###COV term using predictions
  #rho<-cor(pa.predict$fit,cpue.predict$fit)
  #SEP<-aggregate(pa.model$fitted.values,by=list(pa.model$data$year),FUN=sd)
  #SEd<-aggregate(cpue.model$fitted.values,by=list(cpue.model$data$year),FUN=sd)
  #SEPn<-aggregate(pa.model$fitted.values,by=list(pa.model$data$year),FUN=length)
  #SEdn<-aggregate(cpue.model$fitted.values,by=list(cpue.model$data$year),FUN=length)
  #SEP<-SEP$x/sqrt(SEPn$x)
  #SEd<-SEd$x/sqrt(SEdn$x)
  #Covdp<-rho*(SEd*SEP) # From Lo 1992
  #Covdp<-cov(pa.predict$fit,cpue.predict$fit) #Using cov function
  
  ###COV term using dependent variables
  Covdp<-cov(pa.predict$fit^2,cpue.predict$fit^2)+mean(pa.predict$fit^2)*mean(cpue.predict$fit^2)-
    (cov(pa.predict$fit,cpue.predict$fit)+mean(pa.predict$fit)*mean(cpue.predict$fit))^2
  #Covdp<-cov(d$x^2,P$x^2)+P$x^2*d$x^2-(cov(P$x,d$x)+(P$x*d$x))^2
  
  VardP2<-Vard*P$x^2
  VarPd2<-VarP*d$x^2
  #Vari<-VardP2+VarPd2+2*d$x*P$x*Covdp
  Vari<-VardP2+VarPd2+Covdp
  gam.index.sd<-sqrt(Vari)  
  
  gam.index<-exp(gam.raw.index) 
  upper_CI<-exp(gam.raw.index+1.96*gam.index.sd)
  lower_CI<-exp(gam.raw.index-1.96*gam.index.sd)
  
  return(data.frame(years=yrs,gam.index,lower_CI=lower_CI,upper_CI=upper_CI))
}

#' A function to predict CPUE-index from a delta-lognormal GAM with bootstrapped errors
#'
#' This function uses a binomial GAM and a log-normal GAM to compute a delta-lognormal
#' index of abundance with associated 95% confidence intervals estimated by the boot-
#' strapping the data. A dummy variable for year must be included to get the annual abundance
#' index. The function back-transforms the index and the CI's.
#' @param pa.model binomial GAM predicting presence or absence
#' @param cpue.model log-normal GAM predicting abundance where catch is positive
#' @param boot_reps number of bootstrap replications (default = 500)
#' @keywords delta-lognormal gam, survey abundance index, bootstrap
#' @export
#' @examples
#' gam.pa.xvars<-c("inverts","slope","btemp","bdepth")
#' gam.pa.yvar<-ifelse(PA.data[species.name]>0,1,0)
#' gam.pa.form <- as.formula(paste("gam.pa.yvar ~", paste(gam.pa.xvars,collapse="+"),"+as.factor(year)",sep=""))
#' pa.gam <- gam(gam.pa.form, family = binomial, data = PA.data)
#' 
#' CPUE.data<-subset(PA.data,PA.data[species.name]>0)
#' gam.cpue.yvar<-unlist(log(CPUE.data[species.name]))
#' gam.cpue.xvars<-c("inverts","slope","btemp","bdepth")
#' gam.cpue.form <- as.formula(paste("gam.cpue.yvar ~", paste(gam.cpue.xvars,collapse="+"),"+as.factor(year)",sep=""))
#' cpue.gam <- gam(gam.cpue.form, family = gaussian, data = CPUE.data)
#' 
#' predict.gam.bindex(pa.gam,cpue.gam)
#' 
#######bootstrapping method

predict.gam.bindex<-function(pa.model,cpue.model,boot_reps=500){
  pa.data<-pa.model$data
  pa.meds<-t(data.frame(apply(pa.data,MARGIN=2,FUN=median)))
  pa.meds.names<-colnames(pa.meds)
  yrs<-unique(pa.data$year)
  
  pa.meds<-data.frame(rep.row(pa.meds,length(yrs)))
  colnames(pa.meds)<-pa.meds.names
  pa.meds$year<-yrs
  
  pa.predict<-predict(pa.model,newdata=pa.meds,type="response",se.fit=TRUE)
  cpue.predict<-predict(cpue.model,newdata=pa.meds,type="response",se.fit=TRUE)
  gam.raw.index<-pa.predict$fit*cpue.predict$fit
  gam.index<-exp(gam.raw.index)
  
  index_ests<-data.frame(array(0,dim=c(0,2)))
  colnames(index_ests)<-c("yrs","temp.index")
  for(i in 1:boot_reps){
    bootdata1<- sample(1:length(pa.model$data[,1]), replace=TRUE)
    bootdata<-data.frame(pa.model$data[bootdata1,])
    gam.pa.yvar<-pa.model$y[bootdata1]
    boot.pa<-gam(pa.model$formula,family=binomial,data=bootdata)
    
    bootdata1<-sample(1:length(cpue.model$data[,1]), replace=TRUE)
    bootdata<-data.frame(cpue.model$data[bootdata1,])
    gam.cpue.yvar<-cpue.model$y[bootdata1]
    boot.cpue<-gam(cpue.gam$formula,family=gaussian,data=bootdata)
    
    #temp1<-predict.gam.bindex(boot.pa,boot.cpue)
    pa.pboot<-predict(boot.pa,newdata=pa.meds,type="response",se.fit=TRUE)
    cpue.pboot<-predict(boot.cpue,newdata=pa.meds,type="response",se.fit=TRUE)
    temp.raw.index<-pa.pboot$fit*cpue.pboot$fit
    temp.index<-exp(temp.raw.index)  
    temp.index<-cbind(yrs,temp.index)
    
    index_ests<-rbind(index_ests,temp.index)
  }
  
  index_est_sd<-aggregate(index_ests$temp.index,by=list(index_ests$yrs),FUN=sd)
  upper_bootCI<-gam.index+1.96*index_est_sd$x#/sqrt(boot_reps)
  lower_bootCI<-gam.index-1.96*index_est_sd$x#/sqrt(boot_reps)
  
  return(data.frame(years=yrs,gam.index,lower_bootCI=lower_bootCI,upper_bootCI=upper_bootCI))
}