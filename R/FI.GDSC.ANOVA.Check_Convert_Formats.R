gdscANOVA_convertScreeningFormat<-function(screeningFile,savingDir){
  load(screeningFile)
  
  dataM<-ic50_matrix[,2:ncol(ic50_matrix)]
  rownames(dataM)<-ic50_matrix[,1]
  
  did<-colnames(dataM)
  did<-unlist(str_split(did,'_'))[seq(1,ncol(dataM)*2,2)]
  udid<-unique(did)
  
  totalTREATED<-rep(NA,length(udid))
  names(totalTREATED)<-udid
  multipleConc<-rep(NA,length(udid))
  names(multipleConc)<-udid
  
  CONCMAT<-matrix(0,nrow = length(udid),ncol = 3,dimnames = list(udid,c('c1','c2','c3')))
  
  for (i in 1:length(udid)){
    IDX<-which(did==udid[i])
    if(length(IDX)>1){
      cr<-colSums(!is.na(dataM[,IDX]))
      
      if (length(cr)==2){cr<-c(cr,0)}
      totalTREATED[i]<-length(which(rowSums(!is.na(dataM[,IDX]))>0))
      
      ss<-which(rowSums(!is.na(dataM[,IDX]))>1)
      multipleConc[i]<-length(which(rowSums(!is.na(dataM[,IDX]))>1))
      
      CONCMAT[i,]<-cr
    }
    else{
      totalTREATED[i]<-sum(!is.na(dataM[,IDX]))
      CONCMAT[i,]<-c(totalTREATED[i],0,0)
      multipleConc[i]<-0
    }
  }
  
  dd<-cbind(CONCMAT,totalTREATED,multipleConc,multipleConc/totalTREATED)
  dd<-dd[order(multipleConc/totalTREATED,decreasing=TRUE),]
  
  id<-which(dd[,6]>0)
  dd<-dd[id,]
  
  colnames(dd)<-c('max_conc_1','max_conc_2','max_conc_3','total number of IC50s','IC50 from multiple max concentrations','perc')
  #   write.table(dd,file='../../Data/DrugScreening/September2015/MultipleMaxConcDrugs.txt',quote=FALSE,sep='\t')
  
  toMerge<-rownames(dd)[which(dd[,5]<100)]
  
  temp_maxConcTested<-matrix(NA,nrow(dataM),ncol(dataM),dimnames = list(rownames(dataM),colnames(dataM)))
  
  for (i in 1:ncol(temp_maxConcTested)){
    idx<-which(!is.na(dataM[,i]))
    temp_maxConcTested[idx,i]<-as.character(as.numeric(unlist(str_split(colnames(dataM)[i],'_'))[2]))
  }
  
  
  
  IC50s<-dataM[,which(!is.element(did,toMerge))]
  maxConcTested<-matrix(as.character(temp_maxConcTested[,which(!is.element(did,toMerge))]),
                        nrow = nrow(temp_maxConcTested),ncol=length(which(!is.element(did,toMerge))),
                        dimnames = list(rownames(temp_maxConcTested),colnames(temp_maxConcTested)[which(!is.element(did,toMerge))]))
  
  for (i in 1:length(toMerge)){
    idxs<-which(did==toMerge[i])
    IC50s<-cbind(IC50s,rowMeans(dataM[,idxs],na.rm = TRUE))
    
    tmp<-temp_maxConcTested[,idxs]
    tmp[which(is.na(tmp))]<-''
    if (length(idxs)==2){toAdd<-str_trim(paste(tmp[,1],tmp[,2],sep=' '))}
    else{toAdd<-str_trim(paste(tmp[,1],tmp[,2],tmp[,3],sep=' '))}
    maxConcTested<-cbind(maxConcTested,toAdd)
    
    colnames(maxConcTested)[ncol(maxConcTested)]<-paste(toMerge[i],paste(setdiff(unique(c(tmp)),''),collapse='/'),sep='_')
    colnames(IC50s)[ncol(IC50s)]<-toMerge[i]<-paste(toMerge[i],paste(setdiff(unique(c(tmp)),''),collapse='/'),sep='_')
    
  }
  
  save(IC50s,file=paste(savingDir,'IC50s.rdata',sep=''))
  save(maxConcTested,file=paste(savingDir,'maximalTestedConcentrations.rdata',sep=''))
  
}
gdscANOVA_convertDrugOwnershipFormat<-function(drug_prop_file,savingDir){
  
  load('../../Data/DrugScreening/GDSC1000 paper freeze/DRUG_BY_COMPANIES_20150126.rdata')
  
  inGDSC100Paper<-DRUG_BY_COMPANIES$DRUG_ID[which(DRUG_BY_COMPANIES$GDSC1000_paper_set==1)]
  
  tmp<-readRDS(drug_prop_file)
  
  
  load(GDSCANOVA_SETTINGS$gdscANOVA.screening.fn)
  drugDomains<-c('All available in v18','GDSC1000 paper','Web Released', 'Sanger Internal',sort(setdiff(unique(tmp$OWNED_BY),NA)))

  DRUG_BY_COMPANIES<-matrix(0,nrow = nrow(tmp),ncol=length(drugDomains)+1,dimnames = list(tmp$DRUG_ID,c('DRUG_ID',drugDomains)))
  
  drugIds<-unlist(str_split(colnames(IC50s),'_'))[seq(1,ncol(IC50s)*2,2)]
  
  DRUG_BY_COMPANIES[,1]<-as.numeric(rownames(DRUG_BY_COMPANIES))
  
  DRUG_BY_COMPANIES[drugIds,'All available in v18']<-1
  DRUG_BY_COMPANIES[as.character(inGDSC100Paper),'GDSC1000 paper']<-1
  DRUG_BY_COMPANIES[which(tmp$WEBRELEASE=='Y'),'Web Released']<-1
  DRUG_BY_COMPANIES[is.na(tmp$OWNED_BY),'Sanger Internal']<-1
  
  companies<-sort(setdiff(unique(tmp$OWNED_BY),NA))
  
  for (i in 1:length(companies)){
    did<-tmp$DRUG_ID[which(tmp$OWNED_BY==companies[i])]
    DRUG_BY_COMPANIES[as.character(did),companies[i]]<-1
  }
  
  return(DRUG_BY_COMPANIES)
}

