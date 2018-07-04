
#' Function to estimate Biomass or CPUE from stratified-random survey design
#'
#' This function estimates a design-based survey estimator for a stratified-random survey 
#' design. Either a population biomass or CPUE can be output. At the moment the function is
#' set up to incorporate the Gulf of Alaska or the Aleutian Islands bottom trawl survey.
#' It requires CPUE data, Strata designations and total area (in km^2). It can also be used 
#' to calculate an estimate for only the trawlable portion of the strata.
#' 
#' @param Catch Catch-per-unit-of-effort for each bottom trawl haul
#' @param Year Survey year
#' @param Strata Strata designation for each haul 
#' @param Strata_area Area of the strata in km2 where each haul occurs
#' @param Region Either "GOA" (the default) or "AI"
#' @param Proportion_trawlable The proportion of the strata that can be sampled by the survey
#' trawl (default is 1)
#' @param method Either "Biomass" or "CPUE"
#' @keywords stratified random survey, survey abundance index
#' @export
#' @examples
#' design.index<-Stratified_CPUE(Design.data$juvenile_POP_CPUE,Design.data$year,Design.data$STRATUM,Design.data$AREA_KM2,"GOA",1, "CPUE")
#' design.index<-Stratified_CPUE(Design.data$juvenile_POP_CPUE,Design.data$year,Design.data$STRATUM,Design.data$AREA_KM2,"AI",Design.data$Trawlable, "Biomass")

####FUNCTION FOR STRATIFIED BIOMASS ESTIMATE#######################
Stratified_CPUE<-function(Catch,Year,Strata,Strata_area,Region="GOA",Proportion_trawlable=1, method="Biomass"){
#Strata<-Design.data$STRATUM
#Catch<-Design.data$juvenile_POP_CPUE
#Year<-Design.data$year
#Strata_area<-Design.data$AREA_KM2
if(length(Proportion_trawlable)==1){Proportion_trawlable<-rep(Proportion_trawlable,length(Strata_area))}

stratas<-unique(Strata)
years<-unique(Year)
Strat_index1<-array(dim=c(length(years),6))
if(Region=="GOA"){Strat_index2<-array(dim=c(59,8))}
if(Region=="AI"){Strat_index2<-array(dim=c(45,8))}

if(method=="Biomass"){
for(i in 1:length(years)){
for(j in 1:length(stratas)){
CPUE<-subset(Catch,Strata==stratas[j]&Year==years[i])
if(length(CPUE)>1){
CPUEh<-mean(CPUE)
varCPUEh<-var(CPUE)

area<-mean(subset(Strata_area,Strata==stratas[j]&Year==years[i]),na.rm=TRUE)
prop<-mean(subset(Proportion_trawlable,Strata==stratas[j]&Year==years[i]),na.rm=TRUE)
area<-area*prop
Biomassh<-CPUEh*area
varBiomassh<-area^2*varCPUEh/length(CPUE)
tval<-qt(1-.05/2,df=(length(CPUE)-1))
Strat_index2[j,1]<-stratas[j]
Strat_index2[j,2]<-area
Strat_index2[j,3]<-CPUEh
Strat_index2[j,4]<-varCPUEh
Strat_index2[j,5]<-Biomassh
Strat_index2[j,6]<-varBiomassh
Strat_index2[j,7]<-Biomassh-tval*sqrt(varBiomassh)
Strat_index2[j,8]<-Biomassh+tval*sqrt(varBiomassh)
}
if(length(CPUE)<1){Strat_index2[j,1]<-stratas[j]}
}
Biomassy<-sum(Strat_index2[,5],na.rm=TRUE)
varBiomassy<-sum(Strat_index2[,6],na.rm=TRUE)
tval<-qt(1-.05/2,df=(length(stratas)-1))
#tval<-1.96
Strat_index1[i,1]<-years[i]
Strat_index1[i,2]<-Biomassy
Strat_index1[i,3]<-varBiomassy
Strat_index1[i,4]<-sqrt(varBiomassy)
Strat_index1[i,5]<-Biomassy-tval*sqrt(varBiomassy)
Strat_index1[i,6]<-Biomassy+tval*sqrt(varBiomassy)
}
colnames(Strat_index1)<-c("Year","Biomass","Var","SE","Lower_CI","Upper_CI")
Strat_index1<-Strat_index1[order(Strat_index1[,1]),]
return(data.frame(Strat_index1))}

#############################################################################

if(method=="CPUE"){
  for(i in 1:length(years)){
    for(j in 1:length(stratas)){
      CPUE<-subset(Catch,Strata==stratas[j]&Year==years[i])
      if(length(CPUE)>1){
        CPUEh<-mean(CPUE)
        varCPUEh<-var(CPUE)
        
        area<-mean(subset(Strata_area,Strata==stratas[j]&Year==years[i]),na.rm=TRUE)
        prop<-mean(subset(Proportion_trawlable,Strata==stratas[j]&Year==years[i]),na.rm=TRUE)
        area<-area*prop
        Biomassh<-CPUEh*area
        varBiomassh<-area^2*varCPUEh/length(CPUE)
        tval<-qt(1-.05/2,df=(length(CPUE)-1))
        Strat_index2[j,1]<-stratas[j]
        Strat_index2[j,2]<-area
        Strat_index2[j,3]<-CPUEh
        Strat_index2[j,4]<-varCPUEh
        Strat_index2[j,5]<-length(CPUE)
        Strat_index2[j,6]<-varBiomassh
        Strat_index2[j,7]<-Biomassh-tval*sqrt(varBiomassh)
        Strat_index2[j,8]<-Biomassh+tval*sqrt(varBiomassh)
      }
      if(length(CPUE)<1){Strat_index2[j,1]<-stratas[j]}
    }
    CPUEy<-sum(Strat_index2[,3]*Strat_index2[,2],na.rm=TRUE)/sum(Strat_index2[,2],na.rm=TRUE)
    varCPUEy<-sum((Strat_index2[,2]/sum(Strat_index2[,2],na.rm=TRUE))^2*(Strat_index2[,4]/Strat_index2[,5]),na.rm=TRUE)
    tval<-qt(1-.05/2,df=(length(stratas)-1))
    #tval<-1.96
    Strat_index1[i,1]<-years[i]
    Strat_index1[i,2]<-CPUEy
    Strat_index1[i,3]<-varCPUEy
    Strat_index1[i,4]<-sqrt(varCPUEy)
    Strat_index1[i,5]<-CPUEy-tval*sqrt(varCPUEy)
    Strat_index1[i,6]<-CPUEy+tval*sqrt(varCPUEy)
  }
  colnames(Strat_index1)<-c("Year","CPUE","Var","SE","Lower_CI","Upper_CI")
  Strat_index1<-Strat_index1[order(Strat_index1[,1]),]
  return(data.frame(Strat_index1))}
}


#' Function to get the strata area for a stratum
#'
#' This function reads in the GOA and AI stratum tables so that a stratum area can be designated 
#' for each bottom trawl haul where the stratum is known.
#' 
#' @param data data set where the stratum area is needed
#' @param strata_column Strata designation for each haul
#' @param region Either "GOA" or "AI"
#' @keywords stratified random survey
#' @export
#' @examples
#' Design.data<-get_strata_area(Juvenile_POP_Data,"STRATUM","GOA")


get_strata_area<-function(data,strata_column, region){

  if(region=="GOA"){
    data(goa.strata.area)
    strata_area<-merge(data,goa.strata.area,by.x=strata_column,by.y="STRATUM",all.x=TRUE)
    return(strata_area)
  } 
  if(region=="AI"){
    data(goa.strata.area)
    strata_area<-merge(data,ai.strata.area,by.x=strata_column,by.y="STRATUM",all.x=TRUE)
    return(strata_area)
  }
  
}

#' Function to calculate root-mean-squared-error for a model
#'
#' This function calculates the rmse for a model given the observed and predicted
#' observations. It was stolen from Brian Stock.
#' 
#' @param obs CPUE observations
#' @param pred Predicted CPUE from a model
#' @keywords model evaluation
#' @export
#' @examples
#' calc_RMSE(glm.cpue$fitted.values,Juvenile_POP_Data$juvenile_POP_CPUE)
calc_RMSE <- function(pred, obs){
  RMSE <- round(sqrt(mean((pred-obs)^2)),3)
  return(RMSE)
}






