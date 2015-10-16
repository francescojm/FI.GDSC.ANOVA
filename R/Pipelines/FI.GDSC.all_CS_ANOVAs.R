

DRUG_DOMAIN<-'All available in v18'
GDSC_SETTINGS_FILE_FOLDER<-'../../Data/ANOVA_setting_files/all_Analyses_Oct2015/'

createIndividual_GDSCSettingFile<-function(ctype,useMsi){
  gdscSettingsTemplate<-read.csv('HTML_templates_and_elements/FI.GDSC.ANOVA.Settings_file_template.csv',stringsAsFactors=FALSE,header=FALSE)
  
  gdscSettingsTemplate[1,2]<-paste(ctype,'analysis GDSC Settings file')
  gdscSettingsTemplate[2,2]<-Sys.info()['user']  
  gdscSettingsTemplate[3,2]<-date()
  gdscSettingsTemplate[17,2]<-DRUG_DOMAIN
  if(ctype!='PANCAN'){
    gdscSettingsTemplate[18,2]<-'CS'  
  }else{
    gdscSettingsTemplate[18,2]<-'PANCAN'
  }
  gdscSettingsTemplate[19,2]<-useMsi
  if(ctype!='PANCAN'){
    gdscSettingsTemplate[20,2]<-ctype  
  }else{
    gdscSettingsTemplate[20,2]<-'PANCAN'
  }
  
  gdscSettingsTemplate[29,2]<-35
  
  
  
  fn<-paste(ctype,'_',DRUG_DOMAIN,'_',Sys.Date(),'.csv',sep='')
  
  write.table(gdscSettingsTemplate,quote=FALSE,row.names=FALSE,file=paste(GDSC_SETTINGS_FILE_FOLDER,fn,sep=''),sep=',',col.names = FALSE)
  return(fn)
}



CancerTypes<-dir('../../Data/Cell_lines/MultiOmicBEMs/')
CancerTypes<-unlist(str_split(CancerTypes,'_simple_MOBEM.rdata'))[seq(1,length(CancerTypes)*2,2)]
CancerTypes<-setdiff(CancerTypes,c('UCEC','PRAD'))
load('../../Data/Cell_lines/Annotations/IncludeMSIfactor.rdata')

nCancerTypes<-length(CancerTypes)


for (i in 1:nCancerTypes){
  
  
  print(paste('- Executing ',CancerTypes[i],'-specific ANOVA',sep=''))
  print('          Assembling GDSC-Settings Files')
  fn<-createIndividual_GDSCSettingFile(ctype = CancerTypes[i],useMsi = as.logical(IncludeMSIfactor[CancerTypes[i]]))
      
  ANOVA_setting_file<-paste(GDSC_SETTINGS_FILE_FOLDER,fn,sep='')
  source('R/Pipelines/FI.GDSC.ANOVA.Pipeline.R')
  source('R/Pipelines/FI.GDSC.ANOVA.Graphic_Pipeline.R')
  
  rm(list=setdiff(ls(),c('i',
                         'nCancerTypes',
                         'CancerTypes',
                         'IncludeMSIfactor',
                         'DRUG_DOMAIN',
                         'GDSC_SETTINGS_FILE_FOLDER',
                         'createIndividual_GDSCSettingFile')))
  

}
