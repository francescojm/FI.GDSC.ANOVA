
DRUG_DOMAINS<-c("All available in v18")


for (DRUG_DOMAIN in DRUG_DOMAINS){
  source('~/Desktop/GitHubRepositories/CODE/FI.GDSC.ANOVA/R/Pipelines/FI.GDSC.ANOVA_ind_pack_creation.R')
  rm(list=setdiff(ls(),c('DRUG_DOMAIN','DRUG_DOMAINS')))
}