## ANOVA paths and filenames

annotations.master_list.fn<-'../../Data/Cell_lines/Annotations/MASTER_LIST_20150115.rdata'
annotations.cosmic2sampleMap.fn<-'../../Data/Cell_lines/Annotations/DanielCosmic2sampleNames.rdata'

gdscANOVA.screening.fn<-'../../Data/DrugScreening/simplyIC50s-Final-RMSE-Filtered-0.3.rdata'
gdscANOVA.screening.DRcurves.fn<-'../../Data/DrugScreening/DoseResponseCurves.rdata'
#gdscANOVA.additionalFeatures.fn<-'../../Data/DR.Lancet/Cluster_BEMs/20150925//COREAD.rdata'
gdscANOVA.results.dir<-'../../Results/DR.Lancet/ANOVA/'
gdscANOVA.cellLineFeatures.dir<-'../../Data/Cell_lines/MultiOmicBEMs/'

gdscANOVA.maxTestedConc.fn<-'../../Data/DrugScreening/maximalTestedConcentrations.rdata'
gdscANOVA.drugProperties.fn<-'../../Data/DrugScreening/v17_DRUG_PROPS_01072014.rdata'
gdscANOVA.drugOwnership.fn<-'../../Data/DrugScreening/DRUG_BY_COMPANIES_20150126.rdata'


## ANOVA file name suffixes
gdscANOVA.cellLineFeatures.fn_suffix<-'_simple_MOBEM.rdata'

## ANOVA settings
gdscANOVA.settings.SCREENING_VERSION<-paste('v17 - up to',Sys.Date())
gdscANOVA.settings.DRUG_domain<-"GDSC1000_paper_set"
#gdscANOVA.settings.CELL_LINES<-'COREAD'
#gdscANOVA.settings.includeMSI_Factor<-TRUE
gdscANOVA.settings.analysisType<-'CS'
gdscANOVA.settings.additionalFeatures<-TRUE
gdsdANOVA.settings.additionalFeaturesOnly<-TRUE
gdscANOVA.settings.excludeHyperMetData<-TRUE
gescANOVA.settings.oneGeneOnly<-NULL
gdscANOVA.settings.resPackageHeader<-paste('DR.Lancet Project ANOVA',gdscANOVA.settings.CELL_LINES,'specific result package')
gdscANOVA.settings.resPackageIncludeSangerLogo<-FALSE
gdscANOVA.settings.resPackageIncludeNKILogo<-FALSE
gdscANOVA.settings.resPackageIncludeEBIlogo<-TRUE
gdscANOVA.settings.featFactorPopulationTh<-3
gdscANOVA.settings.MSIfactorPopulationTh<-2

gdscANOVA.settings.pval_correction_method<-'fdr'

if (gdscANOVA.settings.CELL_LINES == 'PANCAN'){
  gdscANOVA.settings.pval_TH<-Inf  
}else{
  gdscANOVA.settings.pval_TH<-Inf
}

gdscANOVA.settings.FDR_TH<-25

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


