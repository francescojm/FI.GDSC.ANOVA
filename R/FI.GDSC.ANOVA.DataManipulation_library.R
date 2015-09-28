###########################################################################
#                                                                         #
#   Project: GDSC1000_ANOVA                                               #
#                                                                         #  
#   File: DataManipualation_library.R                                     #  
#                                                                         #
###########################################################################

## data manipulation
gdscANOVA_diagnostics<-function(){
  featFactorPopulationTh=gdscANOVA.settings.featFactorPopulationTh
  MSIfactorPopulationTh=gdscANOVA.settings.MSIfactorPopulationTh
  print('- Calculating the number of feasible tests...')
  
  drugs<-colnames(IC50s)
  features<-rownames(InputFeatures$BEM)
  
  ndrugs<-length(drugs)
  nfeatures<-length(features)
  
  nDrugFeat_combos<-ndrugs*nfeatures
  
  feasibleTests<-0
  
  pb<-my.Pbar.create(max = ndrugs*nfeatures)
  flag<-0
  
  for (d in 1:ndrugs){
    
    commonC<-intersect(rownames(IC50s),colnames(InputFeatures$BEM))
    commonC<-commonC[!is.na(IC50s[commonC,drugs[d]])]
    
    IC50pattern<-IC50s[commonC,drugs[d]]
    
    
    commonC<-commonC[!is.na(IC50s[commonC,drugs[d]])]
    TISSUE_FACTOR<-InputFeatures$TISSUES[commonC]
    MSI_FACTOR<-InputFeatures$MSI_VARIABLE[commonC]
    
    for (f in 1:nfeatures){
      flag<-flag+1
      my.Pbar.progres(pb = pb,i = flag)
      FEATpattern<-InputFeatures$BEM[features[f],commonC]
      
      FEATpattern<-FEATpattern[!is.na(IC50pattern)]
      IC50pattern<-IC50pattern[!is.na(IC50pattern)]
      
      FEATpattern[which(FEATpattern==0)]<-'neg'
      FEATpattern[which(FEATpattern=='1')]<-'pos'
      
      TISSUEpattern<-as.factor(TISSUE_FACTOR[names(IC50pattern)])
      MSIpattern<-as.factor(MSI_FACTOR[names(IC50pattern)])
      
      if (gdscANOVA.settings.includeMSI_Factor & length(which(FEATpattern=='pos'))>=featFactorPopulationTh &
            length(which(FEATpattern=='neg'))>=featFactorPopulationTh &
            length(which(MSIpattern==0))>=MSIfactorPopulationTh &
            length(which(MSIpattern==1))>=MSIfactorPopulationTh){
        feasibleTests<-feasibleTests+1
      } 
      if (!gdscANOVA.settings.includeMSI_Factor & length(which(FEATpattern=='pos'))>=featFactorPopulationTh &
            length(which(FEATpattern=='neg'))>=featFactorPopulationTh){
        feasibleTests<-feasibleTests+1
      } 
      
    }
  }
  
  my.Pbar.close(pb)
  print('+ Done!')
  
  print(paste('n. of feasible tests =',feasibleTests))
  print(paste(format(100*feasibleTests/nDrugFeat_combos,digits = 3),'% of all the drug/feature combinations',sep=''))
  
  return(list(N_DRUGS=ndrugs,N_FEATURES=nfeatures,N_COMBOS=nDrugFeat_combos,N_FEASIBLE_TESTS=feasibleTests,
              PERC_FEAS_TESTS=100*feasibleTests/nDrugFeat_combos))
}
gdscANOVA_create_systemInfos<-function(){
  
  ANALYSIS_SYSTEMS_INFOS<-list()
  
  ANALYSIS_SYSTEMS_INFOS$DATE<-DATE
  ANALYSIS_SYSTEMS_INFOS$TIME<-TIME
  ANALYSIS_SYSTEMS_INFOS$USER<-SYSINFO[7]
  ANALYSIS_SYSTEMS_INFOS$MACHINE<-SYSINFO[4]
  ANALYSIS_SYSTEMS_INFOS$SCREENING_VERSION<-gdscANOVA.settings.SCREENING_VERSION
  ANALYSIS_SYSTEMS_INFOS$DRUG_DOMAIN<-gdscANOVA.settings.DRUG_domain
  if(gdscANOVA.settings.CELL_LINES!='PANCAN'){
    analysis_type<-paste(gdscANOVA.settings.CELL_LINES,'specific')
  }else{
    analysis_type<-gdscANOVA.settings.CELL_LINES
  }
  
  ANALYSIS_SYSTEMS_INFOS$CELL_LINES_DOMAIN<-analysis_type
  ANALYSIS_SYSTEMS_INFOS$TISSUE_FACTOR<-gdscANOVA.settings.CELL_LINES=='PANCAN'
  ANALYSIS_SYSTEMS_INFOS$MSI_FACTOR<-gdscANOVA.settings.includeMSI_Factor
  
  ANALYSIS_SYSTEMS_INFOS$featFactorPopulationTh<-gdscANOVA.settings.featFactorPopulationTh
  if(gdscANOVA.settings.includeMSI_Factor){
    ANALYSIS_SYSTEMS_INFOS$MSI_FACTOR<-gdscANOVA.settings.MSIfactorPopulationTh
  }else{
    ANALYSIS_SYSTEMS_INFOS$MSI_FACTOR<-NA
  }
  
  ANALYSIS_SYSTEMS_INFOS$N_DRUGS<-DIAGNOSTICS$N_DRUGS
  ANALYSIS_SYSTEMS_INFOS$N_FEATURES<-DIAGNOSTICS$N_FEATURES
  ANALYSIS_SYSTEMS_INFOS$N_COMBOS<-DIAGNOSTICS$N_COMBOS
  ANALYSIS_SYSTEMS_INFOS$N_FEASIBLE_TESTS<-DIAGNOSTICS$N_FEASIBLE_TESTS
  ANALYSIS_SYSTEMS_INFOS$PERC_FEAS_TESTS<-DIAGNOSTICS$PERC_FEAS_TESTS
  
  save(ANALYSIS_SYSTEMS_INFOS,file=paste(MAIN_DIR,'SYSTEM_INFOS.rdata',sep=''))
}
gdscANOVA_createDrugDataInput<-function(DRUG_DOMAIN,diagnostic=TRUE){
  
  drugIdxs<-which(DRUG_BY_COMPANIES[,gdscANOVA.settings.DRUG_domain]==1)
  
  ndDomain<-length(drugIdxs)
  
  cat(paste(ndDomain,'drugs in the selected domain'))
  
  drugIdxs<-drugIdxs[which(DRUG_BY_COMPANIES$dataAvailableInV17[drugIdxs]==1)]
  ndDomain<-length(drugIdxs)
  cat(paste(' of which ',ndDomain,' available in the current version of the screening',sep=''))
  
  DRUGIDS<-as.character(DRUG_BY_COMPANIES$DRUG_ID[drugIdxs])
  IC50s<-IC50s[,DRUGIDS]
  
  if (diagnostic){
    
    totalSc<-nrow(IC50s)
    drugperCell<-rowSums(sign(abs(IC50s)),na.rm=TRUE)
    totalScAV<-length(which(drugperCell>0))
    cat(paste('\n\t\t\t*d* ',totalSc,' cell lines listed in the screening (of which ',totalScAV,' tested againts at least 1 drug)',sep=''))
    totalSc<-ncol(IC50s)
    cellperDrug<-colSums(sign(abs(IC50s)),na.rm=TRUE)
    totalScAV<-length(which(cellperDrug>0))
    cat(paste('\n\t\t\t*d* ',totalSc,' drugs listed in the screening (of which ',totalScAV,' tested on at least 1 cell line)',sep=''))
    cat(paste('\n\t\t\t*d* ',' median number of drugs screened on each cell line = ',median(drugperCell),' (min = ',min(drugperCell),', max = ',max(drugperCell),')',sep=''))
    cat(paste('\n\t\t\t*d* ',' median number of cell lines screened against each drug = ',median(cellperDrug),' (min = ',min(cellperDrug),', max = ',max(cellperDrug),')',sep=''))
  }
  
  return(IC50s)
}
gdscANOVA_createInputFeatures<-function(additional_features=NULL,additional_features_only=FALSE,
                                        excludeHyperMetData=TRUE,oneGeneOnly=NULL){
  
  if(additional_features_only){
    BEM<-additional_features
  }else{
    additional_features<-additional_features[,colnames(TOTALBEM)]
    BEM<-rbind(TOTALBEM,additional_features)  
  }
  
  rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ':',replacement = '_')
  rownames(BEM)<-str_replace_all(rownames(BEM),pattern = ' ',replacement = '_')
  
  if(excludeHyperMetData){
    idxs<-which(!str_detect(rownames(BEM),'HypMET'))
    BEM<-BEM[idxs,]  
  }
  
  if(length(oneGeneOnly)>0){
    
    for (kk in 1:length(oneGeneOnly)){
      if (kk == 1){
        idxs<-which(str_detect(rownames(BEM),oneGeneOnly[kk]))
      }else{
        idxs<-union(idxs,which(str_detect(rownames(BEM),oneGeneOnly[kk]))) 
      }
      
    }
    
    if (length(idxs)==1){
      BEM<-matrix(BEM,1,ncol(BEM),dimnames = list(rownames(BEM)[idxs],colnames(BEM)))
    }else{
      BEM<-BEM[idxs,]  
    }
  }
  
  
  
  idxs<-which(rowSums(abs(BEM))>2)
  
  if (length(idxs)==1){
    BEM<-matrix(BEM,1,ncol(BEM),dimnames = list(rownames(BEM)[idxs],colnames(BEM)))
  }else{
    BEM<-BEM[idxs,]  
  }
  
  if(length(idxs)>1){
    cat(paste('\n\t','Features with n > 2 positive samples = ',nrow(BEM)))
    cat('\n\tMerging features with identical positive patterns:')
    BEM<-my.compress_identical_patterns(BEM)
  }
  BEM<-abs(BEM)
  cat(paste('\n\t','Features with n > 2 positive samples = ',nrow(BEM)))
  
  TISSUE_VARIABLE<-gdscANOVA_createTissueVariable(colnames(BEM))
  
  cat(paste('\t*d* ','Tissue variable with',length(unique(TISSUE_VARIABLE)),' factors'))
  
  if (nrow(BEM)>1){
    BEM<-BEM[,names(TISSUE_VARIABLE)]  
  }else{
    BEM<-matrix(BEM,nrow = 1,ncol = ncol(BEM),dimnames = list(rownames(BEM),colnames(BEM)))
  }
  
  
  MSI_VARIABLE<-as.character(MASTER_LIST[as.character(colnames(BEM)),'MMR'])
  
  MSI_VARIABLE[which(MSI_VARIABLE!='MSI-H')]<-0
  MSI_VARIABLE[which(MSI_VARIABLE!='0')]<-1
  MSI_VARIABLE<-as.numeric(MSI_VARIABLE)
  names(MSI_VARIABLE)<-colnames(BEM)
  
  return(list(BEM=BEM,TISSUES=TISSUE_VARIABLE,MSI_VARIABLE=MSI_VARIABLE))
}

gdscANOVA_createTissueVariable<-function(CID){
  
  TISSUE_VARIABLE<-MASTER_LIST$GDSC.description_1
  names(TISSUE_VARIABLE)<-rownames(MASTER_LIST)
  
  TISSUE_VARIABLE[which(TISSUE_VARIABLE=='digestive_system')]<-MASTER_LIST$GDSC.description_2[which(TISSUE_VARIABLE=='digestive_system')]
  TISSUE_VARIABLE[which(TISSUE_VARIABLE=='urogenital_system')]<-MASTER_LIST$GDSC.description_2[which(TISSUE_VARIABLE=='urogenital_system')]
  
  TISSUE_VARIABLE<-TISSUE_VARIABLE[CID]
  
  S<-summary(as.factor(TISSUE_VARIABLE))
  
  TISSUE_VARIABLE<-TISSUE_VARIABLE[is.element(TISSUE_VARIABLE,names(which(S>2)))]
  
  
  return(TISSUE_VARIABLE)
}

my.compress_identical_patterns<-function(DATA){
  
  flag<-1
  
  while(flag){
    JD<-as.matrix(dist(DATA,method='binary'))
    JD[lower.tri(JD,diag=TRUE)]<-Inf
    
    toMix<-which(JD==0,arr.ind=TRUE)
    
    if (length(toMix)>0){
      toMix<-toMix[1,]
      tm1<-row.names(JD)[toMix[1]]
      tm2<-row.names(JD)[toMix[2]]
      cat(paste('\n\t\t\t\t\t\t\tMerging',tm1,'and',tm2))  
      rownames(DATA)[toMix[2]]<-paste(tm1,tm2,sep=', ')
      DATA<-DATA[-toMix[1],]
    }else{
      flag<-0
    }
    
  }
  
  print('Done!')
  return(DATA)
}