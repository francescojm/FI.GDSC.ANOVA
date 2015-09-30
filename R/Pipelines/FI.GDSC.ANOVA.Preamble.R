gdscANOVA_Preamble<-function(ANOVA_setting_file='HTML_templates_and_elements/FI.GDSC.ANOVA.Settings_file_template.csv'){
  
  settings<-read.csv(ANOVA_setting_file,stringsAsFactors=FALSE)
  settings<-settings[setdiff(1:nrow(settings),c(1:3,11,13))]
  
  GLOBAL_SETTINGS<-settings$Analysis.Description
  
  
  load(gdscANOVA.screening.fn)
  load(paste(gdscANOVA.cellLineFeatures.dir,gdscANOVA.settings.CELL_LINES,gdscANOVA.cellLineFeatures.fn_suffix,sep=''))
  
  TOTALBEM <-MoBEM
  
  load(gdscANOVA.drugProperties.fn)
  load(gdscANOVA.drugOwnership.fn)
  load(gdscANOVA.maxTestedConc.fn)
  load(gdscANOVA.screening.DRcurves.fn)
  
  load(annotations.master_list.fn)
  load(annotations.cosmic2sampleMap.fn)
  
  if(gdscANOVA.settings.additionalFeatures){
    print('Loading Additional Features...')
    load(gdscANOVA.additionalFeatures.fn)
    gdscANOVA.additionalFeatures=BEM
    print('+ Done!')
  }else{
    gdscANOVA.additionalFeatures=NULL
  }
  
  
  
  
}

