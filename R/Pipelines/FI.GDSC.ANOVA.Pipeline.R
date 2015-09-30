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
cat('\n+ Done!\n\n')

#cat('- Running gdscANOVA preamble (loading R objects, setting paths and file names)')
#source('R/FI.GDSC.ANOVA.Preamble.R')
#cat('\n+ Done!\n\n')

cat('\n\n- Loading gdscANOVA functions for Data Manipulations and Statistics...')
source('R/FI.GDSC.ANOVA.DataManipulation_library.R')
source('R/FI.GDSC.ANOVA.Statistics_library.R')
source('R/FI.GDSC.ANOVA.UserInteraction.R')
cat('\n+ Done!\n\n')

TIME<-str_replace_all(TIME,'[:]','-')
TIME<-str_replace_all(TIME,'[ ]','_')

MAIN_DIR<-paste(gdscANOVA.results.dir,'/',
                gdscANOVA.settings.CELL_LINES,'_',
                gdscANOVA.settings.DRUG_domain,'_',
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

cat(paste('- Retrieving data for the selected drug domain (',gdscANOVA.settings.DRUG_domain,')...',sep=''))
cat(paste('\n\t\t'))

IC50s<-gdscANOVA_createDrugDataInput(gdscANOVA.settings.DRUG_domain)
cat(paste('\n\t\t(',gdscANOVA.settings.SCREENING_VERSION,')',sep=''))
cat('\n+ Done!\n\n')
 
cat('- Assembling input features')
cat(paste('\n\t\t'))
InputFeatures<-gdscANOVA_createInputFeatures(additional_features=gdscANOVA.additionalFeatures,
                                              oneGeneOnly = gescANOVA.settings.oneGeneOnly,
                                              excludeHyperMetData = gdscANOVA.settings.excludeHyperMetData,
                                              additional_features_only = gdsdANOVA.settings.additionalFeaturesOnly)
cat('\n+ Done!\n\n')
 
save(InputFeatures,file=paste(INPUT_DIR,'InputFeatures.rdata',sep=''))
save(IC50s,file=paste(INPUT_DIR,'IC50s.rdata',sep=''))
 
DIAGNOSTICS<-gdscANOVA_diagnostics()
 
gdscANOVA_create_systemInfos()

TOTRES<-gdscANOVA_totalANOVA(fn=paste(OUTPUT_DIR,'ANOVA_results',sep=''))

write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(OUTPUT_DIR,'ANOVA_results','.txt',sep=''))

save(TOTRES,file=paste(OUTPUT_DIR,'ANOVA_results','.rdata',sep=''))


