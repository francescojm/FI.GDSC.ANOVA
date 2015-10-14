#' @title ANOVA preamble
#' 
#' @description
#' \code{gdscANOVA_Preamble(ANOVA_setting_file)} returns a list of setting variables 
#' initilized with the values contained in the ANOVA_settig_file [see an example of this file
#' in the HTML_templates_and_elements/FI.GDSC.ANOVA.Settings_file_template.csv]
#' 
gdscANOVA_Preamble<-function(ANOVA_setting_file='HTML_templates_and_elements/FI.GDSC.ANOVA.Settings_file_template.csv'){
  
  settings<-read.csv(ANOVA_setting_file,stringsAsFactors=FALSE,header=FALSE)
  settings<-settings[setdiff(1:nrow(settings),c(1:4,13,15)),]
  
  GLOBAL_SETTINGS<-settings$V2
  names(GLOBAL_SETTINGS)<-settings$V3
  
  GLOBAL_SETTINGS<-as.list(GLOBAL_SETTINGS)
  
  GLOBAL_SETTINGS$gdscANOVA.settings.includeMSI_Factor<-as.logical(GLOBAL_SETTINGS$gdscANOVA.settings.includeMSI_Factor)
  GLOBAL_SETTINGS$gdsdANOVA.settings.additionalFeaturesOnly<-as.logical(GLOBAL_SETTINGS$gdsdANOVA.settings.additionalFeaturesOnly)
  GLOBAL_SETTINGS$gdscANOVA.settings.excludeHyperMetData<-as.logical(GLOBAL_SETTINGS$gdscANOVA.settings.excludeHyperMetData)
  GLOBAL_SETTINGS$gdscANOVA.settings.featFactorPopulationTh<-as.numeric(GLOBAL_SETTINGS$gdscANOVA.settings.featFactorPopulationTh)
  GLOBAL_SETTINGS$gdscANOVA.settings.MSIfactorPopulationTh<-as.numeric(GLOBAL_SETTINGS$gdscANOVA.settings.MSIfactorPopulationTh)
  GLOBAL_SETTINGS$gdscANOVA.settings.pval_TH<-as.numeric(GLOBAL_SETTINGS$gdscANOVA.settings.pval_TH)
  GLOBAL_SETTINGS$gdscANOVA.settings.FDR_TH<-as.numeric(GLOBAL_SETTINGS$gdscANOVA.settings.FDR_TH)
  
  if (!as.logical(GLOBAL_SETTINGS$gdscANOVA.settings.additionalFeatures)){GLOBAL_SETTINGS$gdscANOVA.settings.additionalFeatures<-NULL}
  
  return(GLOBAL_SETTINGS)
}


#' @title ANOVA preamble loading
#' 
#' @description
#' \code{gdscANOVA_Preamble_loading(global_setting_list)} loads and export in the workspace, as global variables 
#' objects and list, whose path is contained in the global_setting_list
#' 
#' 
gdscANOVA_Preamble_loading<-function(global_setting_list){
  load(global_setting_list$gdscANOVA.screening.fn)
  IC50s<<-IC50s
  
  load(paste(global_setting_list$gdscANOVA.cellLineFeatures.dir,
             global_setting_list$gdscANOVA.settings.CELL_LINES,
             global_setting_list$gdscANOVA.cellLineFeatures.fn_suffix,sep='')) 
  TOTALBEM<<-MoBEM
   
  load(global_setting_list$gdscANOVA.drugProperties.fn)
  DRUG_PROPS<<-DRUG_PROPS
  
  load(global_setting_list$gdscANOVA.maxTestedConc.fn)
 
  maxConcTested<<-maxConcTested
 # load(gdscANOVA.screening.DRcurves.fn)
  load(global_setting_list$annotations.master_list.fn)
  MASTER_LIST<<-MASTER_LIST
 
#   load(annotations.cosmic2sampleMap.fn)
#   
   if(length(global_setting_list$gdscANOVA.settings.additionalFeatures)>0){
     print('Loading Additional Features...')
     load(global_setting_list$gdscANOVA.additionalFeatures.fn)
     gdscANOVA.additionalFeatures<<-BEM
     print('+ Done!')
   }else{
     gdscANOVA.additionalFeatures<<-NULL
   }
}



