
source('R/FI.GDSC.ANOVA.Graphic_Library.R')
current_dir<-MAIN_DIR

load(paste(current_dir,'OUTPUT/ANOVA_results.rdata',sep=''))
load(paste(current_dir,'INPUT/InputFeatures.rdata',sep=''))

OUT_GRAPH_DIR<-paste(current_dir,'OUTPUT/GRAPHICS/',sep='')

if(!file.exists(OUT_GRAPH_DIR)){
  dir.create(OUT_GRAPH_DIR)
}

association_plot_dir<-paste(OUT_GRAPH_DIR,'association_plots/',sep='')

if(!file.exists(association_plot_dir)){
  dir.create(association_plot_dir)
}

range<-which(as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<gdscANOVA.settings.FDR_TH & as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<gdscANOVA.settings.pval_TH)


if (length(range)>0){
  gdscANOVA_createSCATTERS(range=range,PATH=association_plot_dir)  
}

if(length(range)>5){
  FeatureSummary<-gdscANOVA_computeFeatureStats(redTOTRES=TOTRES,fdrTH = gdscANOVA.settings.FDR_TH,pvalTH = gdscANOVA.settings.pval_TH)
  DrugSummary<-gdscANOVA_computeDrugStats(redTOTRES=TOTRES,fdrTH = gdscANOVA.settings.FDR_TH,pvalTH = gdscANOVA.settings.pval_TH)
}

delta<-sign(as.numeric(TOTRES[,"FEATURE_deltaMEAN_IC50"]))*as.numeric(TOTRES[,"FEATURE_IC50_effect_size"])
pval<-as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])
fdr<-as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])
N<-as.numeric(TOTRES[,"N_FEATURE_pos"])
labels<-paste(TOTRES[,'Drug name'],' [',TOTRES[,'Drug Target'],']\n',TOTRES[,'FEATURE'])

png(paste(current_dir,'OUTPUT/GRAPHICS/comprehensive_volcanoPlot.png',sep=''),768,1024)
gdscANOVA_volcanoPlot_T(delta=delta,pval=pval,qvals=fdr,N=N,effth=1,fdrth = gdscANOVA.settings.FDR_TH)
dev.off()

vdPATH<-paste(paste(current_dir,'OUTPUT/GRAPHICS/DRUG_volcanos/',sep=''))
if(!file.exists(vdPATH)){
  dir.create(vdPATH) 
}

gdscANOVA_allDRUGS_volcanoPlots(TOTRES,PATH=vdPATH)

vfPATH<-paste(paste(current_dir,'OUTPUT/GRAPHICS/FEATURE_volcanos/',sep=''))
if(!file.exists(vfPATH)){
  dir.create(vfPATH)
}

gdscANOVA_allFEATURES_volcanoPlots(TOTRES,PATH=vfPATH)
