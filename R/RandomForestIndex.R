#' A function to predict CPUE-index from a random forest model with bootstrapped errors
#'
#' This function uses a random forest model to produce and index of abundance with 
#' associated 95% confidence intervals estimated by the boot-strapping the data. 
#' A dummy variable for year must be included to get the annual abundance
#' index. The function back-transforms the index and the CI's.
#' @param data data used to formulate the random forest model
#' @param rf.yvar y variables (CPUE) in the random forest model
#' @param rf.xvars x variables in the random forest model
#' @param rf.form formula for the random forest model
#' @param rf.model random forest model to bootstrap
#' @param boot_reps number of bootstrap replications (default = 500)
#' @keywords random forest, survey abundance index, bootstrap
#' @export
#' @examples
#' rf.cpue <- randomForest(rf.form,data=PA.data, mtry=3, ntree=1000, importance=TRUE, do.trace=250, keep.forest=TRUE) 
#' print(rf.cpue)	
#' rf.index<-predict.rf.index(PA.data,rf.yvar,rf.xvars,rf.form,rf.cpue)
#' 
#######bootstrapping method

predict.rf.index<-function(data,rf.yvar,rf.xvars,rf.form,rf.model,boot_reps=500){
#data<-PA.data
#rf.model<-rf.cpue
  pa.meds<-as.numeric(t(apply(cbind(data[rf.xvars]),MARGIN=2,FUN="median",na.rm=TRUE)))
  pa.meds<-t(matrix(pa.meds,byrow=FALSE))
#  pa.meds.names<-colnames(data)
  year<-unique(data$year)
  
  pa.meds<-rep.row(pa.meds,length(year))
  colnames(pa.meds)<-rf.xvars
  #colnames(pa.meds)[which(colnames(pa.meds)=="year")]<-"y1"
  pa.meds<-data.frame(pa.meds,year)
  pa.meds<-pa.meds[order(pa.meds$year),]
  
  rf.raw.index<-predict(rf.model,newdata=pa.meds,type="response")
  rf.index<-exp(rf.raw.index)-.5*min(subset(rf.yvar,rf.yvar>0))

  index_ests<-data.frame(array(0,dim=c(0,2)))
  colnames(index_ests)<-c("yrs","temp.index")
  for(i in 1:boot_reps){
    bootdata1<- sample(1:length(data[,1]), replace=TRUE)
    bootdata<-data.frame(data[bootdata1,])
    rf.yvar<-rf.yvar[bootdata1]
    boot.pa<-randomForest(rf.form,data=data, mtry=3, ntree=1000, importance=TRUE ) 
    #print(mean(boot.pa$rsq))

    #temp1<-predict.glm.bindex(boot.pa,boot.cpue)
    temp.raw.index<-predict(boot.pa,newdata=pa.meds,type="response")
    rf.index<-exp(rf.raw.index)
    temp.index<-exp(temp.raw.index)-.5*min(subset(rf.yvar,rf.yvar>0))  
    temp.index<-cbind(yrs=pa.meds$year,temp.index)
    
    index_ests<-rbind(index_ests,temp.index)
  }
  
  index_est_sd<-aggregate(index_ests$temp.index,by=list(index_ests$yrs),FUN=sd)
  upper_bootCI<-rf.index+1.96*index_est_sd$x#/sqrt(boot_reps)
  lower_bootCI<-rf.index-1.96*index_est_sd$x#/sqrt(boot_reps)
  
  return(data.frame(years=pa.meds$year,rf.index,lower_bootCI=lower_bootCI,upper_bootCI=upper_bootCI))
}