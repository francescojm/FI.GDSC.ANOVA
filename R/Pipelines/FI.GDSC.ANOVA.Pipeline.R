ANOVA_setting_file<-'../../Data/ANOVA_setting_files/FI.GDSC.ANOVA.Settings_file_template.csv'

DATE<-Sys.Date()
TIME<-Sys.time()
SYSINFO<-Sys.info()

cat(paste('\n\n**** gdscANOVA analysis - ',TIME,' ****',sep=''))
cat(paste('\nExecuted by ',SYSINFO[7],' on ',SYSINFO[4],sep=''))
cat('\n*********************************************************')

cat('\n\n- Loading gdscANOVA packages...')
library(beeswarm)
library(stringr)
library(qvalue)
library(djvMixedIC50)
library(scales)
library(roxygen2)
cat('\n+ Done!\n\n')

cat('- Running gdscANOVA preamble (loading R objects, setting paths and file names)')
source('R/FI.GDSC.ANOVA.Preamble_Library.R')
GDSCANOVA_SETTINGS<-gdscANOVA_Preamble(ANOVA_setting_file = ANOVA_setting_file)
gdscANOVA_Preamble_loading(GDSCANOVA_SETTINGS)
cat('\n+ Done!\n\n')

cat('\n\n- Loading gdscANOVA functions for Data Manipulations, Statistics and User interaction...')
source('R/FI.GDSC.ANOVA.DataManipulation_library.R')
source('R/FI.GDSC.ANOVA.Statistics_library.R')
source('R/FI.GDSC.ANOVA.UserInteraction.R')
cat('\n+ Done!\n\n')
 
TIME<-str_replace_all(TIME,'[:]','-')
TIME<-str_replace_all(TIME,'[ ]','_')
 
MAIN_DIR<-paste(GDSCANOVA_SETTINGS$gdscANOVA.results.dir,'/',
                 GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES,'_',
                 GDSCANOVA_SETTINGS$gdscANOVA.settings.DRUG_domain,'_',
                 TIME,'/',sep='')

INPUT_DIR<-paste(MAIN_DIR,'INPUT/',sep='')
OUTPUT_DIR<-paste(MAIN_DIR,'OUTPUT/',sep='')

if(!exists(MAIN_DIR)){
  dir.create(MAIN_DIR)
}

if(!exists(INPUT_DIR)){
  dir.create(INPUT_DIR)
}

if(!exists(OUTPUT_DIR)){
  dir.create(OUTPUT_DIR)
}

cat(paste('- Retrieving data for the selected drug domain (',GDSCANOVA_SETTINGS$gdscANOVA.settings.DRUG_domain,')...',sep=''))
cat(paste('\n\t\t'))
 
IC50s<-gdscANOVA_createDrugDataInput(GDSCANOVA_SETTINGS$gdscANOVA.settings.DRUG_domain)
cat(paste('\n\t\t(',GDSCANOVA_SETTINGS$gdscANOVA.settings.SCREENING_VERSION,')',sep=''))
cat('\n+ Done!\n\n')
  
cat('- Assembling input features')
cat(paste('\n\t\t'))
InputFeatures<-gdscANOVA_createInputFeatures(additional_features=GDSCANOVA_SETTINGS$gdscANOVA.settings.additionalFeatures,
                                               oneGeneOnly = GDSCANOVA_SETTINGS$gescANOVA.settings.oneGeneOnly,
                                               excludeHyperMetData = GDSCANOVA_SETTINGS$gdscANOVA.settings.excludeHyperMetData,
                                               additional_features_only = GDSCANOVA_SETTINGS$gdsdANOVA.settings.additionalFeaturesOnly)
cat('\n+ Done!\n\n')
  
save(InputFeatures,file=paste(INPUT_DIR,'InputFeatures.rdata',sep=''))
save(IC50s,file=paste(INPUT_DIR,'IC50s.rdata',sep=''))
  
DIAGNOSTICS<-gdscANOVA_diagnostics()
  
gdscANOVA_create_systemInfos()
 
TOTRES<-gdscANOVA_totalANOVA(fn=paste(OUTPUT_DIR,'ANOVA_results',sep=''))
 
write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(OUTPUT_DIR,'ANOVA_results','.txt',sep=''))
 
save(TOTRES,file=paste(OUTPUT_DIR,'ANOVA_results','.rdata',sep=''))
