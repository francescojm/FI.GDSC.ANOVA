
CancerTypes<-dir('../../Data//DR.Lancet//Cluster_BEMs//20150925/')
CancerTypes<-unlist(str_split(CancerTypes,'.rdata'))[seq(1,length(CancerTypes)*2,2)]
CancerTypes<-setdiff(CancerTypes,c('CESC','ESCA','MESO','GBM'))


load('../../Data/Cell_lines/Annotations/IncludeMSIfactor.rdata')

nCancerTypes<-length(CancerTypes)

for (i in 1:nCancerTypes){
  print(paste('- Executing ',CancerTypes[i],'-specific ANOVA',sep=''))
  
  
  gdscANOVA.settings.CELL_LINES<-CancerTypes[i]
  gdscANOVA.settings.includeMSI_Factor<-(IncludeMSIfactor[i]==1)
  gdscANOVA.additionalFeatures.fn<-paste('../../Data/DR.Lancet/Cluster_BEMs/20150925//',CancerTypes[i],'.rdata',sep='')
  
  cat('- Running gdscANOVA preamble (loading R objects, setting paths and file names)')
  source('R/FI.GDSC.ANOVA.Preamble.R')
  cat('\n+ Done!\n\n')
    
  source('R/FI.GDSC.ANOVA.Pipeline.R')
  source('R/FI.GDSC.ANOVA.Graphic_Pipeline.R')
  source('R/FI.GDSC.ANOVA_ind_pack_creation.R')
  
  vars<-ls()
  vars<-setdiff(vars,c('CancerTypes','IncludeMSIfactor','nCancerTypes'))
  rm(list=vars)
  
  print('+ DONE!')
}
