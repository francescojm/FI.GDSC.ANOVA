#############################################################################
#                                                                           #
#   Project: FI.GDSC.ANOVA                                                  #
#   Version: 1.0                                                            #
#   Description: Scripts and R functions used to perform a systematic       #
#                analysis of variance to identify statistical interactions  #
#                                                                           #
#   File: FI.GDSC.ANOVA.Statistics_library.R                                #
#                                                                           #
#   Copyright (c) 2014 - 2019, EMBL - European Bioinformatics Institute     #
#                                                                           #
#       Author: Francesco Iorio (iorio@ebi.ac.uk)                           #
#       Distributed under the Apache 2.0 License.                           #
#                                                                           #
#############################################################################


#' @title Effect Size of a statistical interaction through Cohen's D computation
#' @author Francesco Iorio - \email{iorio@@ebi.ac.uk}
#' @description
#' \code{gdscANOVA_cohens_d(x,y)} \cr\cr returns the effect size of the statistical interaction between a variable and an underlying binary factor.
#' It takes in input two vectors containing a partition of the variable observations based on the corresponding factor value (positive or negative).
#' @usage gdscANOVA_cohens_d(x,y) 
#' @param x Numerical vector containing the variable observations corresponding to the positive values of the underlying factor.
#' @param y Numerical vector containing the variable observations corresponding to the negative values of the underlying factor.   
#' @return The Cohen's d, defined as the difference between the means of the two vector, divided by their pooled standard deviation.       
#' @seealso gdscANOVA_glass_Ds
gdscANOVA_cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  cd  <- md/csd                        ## cohen's d
}

#' @title Effect Size of a statistical interaction through individual Glass Delta
#' @author Francesco Iorio - \email{iorio@@ebi.ac.uk}
#' @description
#' \code{gdscANOVA_glass_Ds(x,y)} \cr\cr returns two numerical scores quantifying the effect sizes of the statistical interaction between a variable and an
#' underlying binary factor. It takes in input two vectors containing a partition of the variable observations based on the corresponding factor value
#' (positive or negative).
#' @usage gdscANOVA_glass_Ds(x,y) 
#' @param x Numerical vector containing the variable observations corresponding to the positive values of the underlying factor.
#' @param y Numerical vector containing the variable observations corresponding to the negative values of the underlying factor.   
#' @return A list of two Glass Deltas, defined as the difference between the means of the two vector, divided by the standard deviations of the two groups
#' respectively       
#' @seealso gdscANOVA_glass_Ds
gdscANOVA_glass_Ds<-function(x,y){
  md<-abs(mean(x)-mean(y))
  g1<-md/sd(x)
  g2<-md/sd(y)
  return(list(g1=g1,g2=g2))
}

#' @title Individual ANOVA test
#' @author Francesco Iorio - \email{iorio@@ebi.ac.uk}
#' @description
#' \code{gdscANOVA_individualANOVA(DRUG_ID, FEATURE, display=TRUE, printTOfig=FALSE, PATH='', FN='', FDR=NA, OUTPUT_PATH='')} \cr\cr
#' executes an individual ANOVA test to evaluate the statistical interaction between the IC50 values of a given drug, across cell lines, and an underlying
#' binary factor. \cr\cr
#' Calling this function requires the existance of the following global variables (described in  and created by the functions in
#' \code{FI.GDSC.ANOVA.Preamble_Library.R}) \cr
#' \itemize{
#'   \item IC50s = Matrix with IC50 values
#'   \item InputFeatures = List containing the input features to be correlated with drug response
#'   \item GDSCANOVA_SETTINGS = List of settings variables
#'}
#' @usage A
#' @param A
#' @param A
#' @return A        
#' @seealso A
#'            
gdscANOVA_individualANOVA<-function(DRUG_ID,FEATURE,display=TRUE,printTOfig=FALSE,PATH='',FN='',FDR=NA,OUTPUT_PATH=''){
  
  if (printTOfig){
    #png(paste(PATH,FN,'_00_scatterSet.png',sep=''),width=793.92,height=1122.24)
    png(paste(PATH,FN,'_00_scatterSet.png',sep=''),width=1200,height=800)
  }
  
  commonC<-intersect(rownames(IC50s),colnames(InputFeatures$BEM))
  
  commonC<-commonC[!is.na(IC50s[commonC,DRUG_ID])]
  
  
  TISSUE_FACTOR<-InputFeatures$TISSUES[commonC]
  MSI_FACTOR<-InputFeatures$MSI_VARIABLE[commonC]
  
  IC50pattern<-IC50s[commonC,DRUG_ID]
  names(IC50pattern)<-commonC
  
  
  FEATpattern<-InputFeatures$BEM[FEATURE,commonC]
  
  FEATpattern<-FEATpattern[!is.na(IC50pattern)]
  IC50pattern<-IC50pattern[!is.na(IC50pattern)]
  
  FEATpattern[which(FEATpattern==0)]<-'neg'
  FEATpattern[which(FEATpattern=='1')]<-'pos'
  FEATpattern<-as.factor(FEATpattern)
  
  TISSUEpattern<-as.factor(TISSUE_FACTOR[names(IC50pattern)])
  MSIpattern<-as.factor(MSI_FACTOR[names(IC50pattern)])
  
  Y <- IC50pattern
    
  featFactorPopulationTh<-GDSCANOVA_SETTINGS$gdscANOVA.settings.featFactorPopulationTh
  MSIfactorPopulationTh<-GDSCANOVA_SETTINGS$gdscANOVA.settings.MSIfactorPopulationTh
  
  A<-(GDSCANOVA_SETTINGS$gdscANOVA.settings.includeMSI_Factor & length(which(FEATpattern=='pos'))>=featFactorPopulationTh &
        length(which(FEATpattern=='neg'))>=featFactorPopulationTh &
        length(which(MSIpattern==0))>=MSIfactorPopulationTh &
        length(which(MSIpattern==1))>=MSIfactorPopulationTh)
  
  B<-(!GDSCANOVA_SETTINGS$gdscANOVA.settings.includeMSI_Factor & length(which(FEATpattern=='pos'))>=featFactorPopulationTh &
        length(which(FEATpattern=='neg'))>=featFactorPopulationTh)
  
  if (A | B){
    
    if(GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType=='PANCAN'){
      fit <- aov(Y ~ TISSUEpattern+MSIpattern+FEATpattern)  
    }else{
      if(GDSCANOVA_SETTINGS$gdscANOVA.settings.includeMSI_Factor){
        fit <- aov(Y ~ MSIpattern+FEATpattern)
      }else{
        fit <- aov(Y ~ FEATpattern)
      }
    }
    
    tfit<-t.test(IC50pattern~FEATpattern,var.equal = TRUE)
    
    if (display){
      gdscANOVA_scatterSets(IC50pattern,FEATpattern,DRUG_ID=DRUG_ID,FEATURE=FEATURE)
    }     
    
    Y<-anova(fit)
    
    if(GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType=='PANCAN'){
      FEATURE_PVAL<-Y[3,5]
      MSI_PVAL<-Y[2,5]
      tissue_PVAL<-Y[1,5]
    }else{
      if(GDSCANOVA_SETTINGS$gdscANOVA.settings.includeMSI_Factor){
        FEATURE_PVAL<-Y[2,5]
        MSI_PVAL<-Y[1,5]  
      }else{
        FEATURE_PVAL<-Y[1,5]
        MSI_PVAL<-NA
      }
      tissue_PVAL<-NA
    }
    
    FEATURE_IC50_WTT_pvalue<-tfit$p.value
    
    pos_IC50_MEAN<-mean(IC50pattern[FEATpattern=='pos'])
    neg_IC50_MEAN<-mean(IC50pattern[FEATpattern=='neg'])
    
    deltaMEAN_IC50<-pos_IC50_MEAN-neg_IC50_MEAN
    
    pos_IC50_sd<-sd(IC50pattern[FEATpattern=='pos'])
    neg_IC50_sd<-sd(IC50pattern[FEATpattern=='neg'])
    
    Npos<-length(which(FEATpattern=='pos'))
    Nneg<-length(which(FEATpattern=='neg'))
    
    EFFECTSIZE_IC50<-gdscANOVA_cohens_d(IC50pattern[FEATpattern=='pos'],IC50pattern[FEATpattern=='neg'])
    GLASS_d<-gdscANOVA_glass_Ds(IC50pattern[FEATpattern=='pos'],IC50pattern[FEATpattern=='neg'])
    
    maxC<-log(as.numeric(unique(maxConcTested[,DRUG_ID])))
    
    maxC<-maxC[which(!is.na(maxC))]
    
    if (length(maxC)==1){
      maxC<-c(maxC,NA)
    }
    
    RES<-matrix(c(FEATURE,
                  DRUG_ID,
                  DRUG_PROPS[str_split(DRUG_ID,'_')[[1]][1],"OWNED_BY"],    
                  as.character(DRUG_PROPS[str_split(DRUG_ID,'_')[[1]][1],"DRUG_NAME"]),
                  as.character(DRUG_PROPS[str_split(DRUG_ID,'_')[[1]][1],"PUTATIVE_TARGET"]),
                  Npos,Nneg,
                  maxC[1],maxC[2],
                  pos_IC50_MEAN,
                  neg_IC50_MEAN,
                  deltaMEAN_IC50,
                  pos_IC50_sd,
                  neg_IC50_sd,
                  EFFECTSIZE_IC50,
                  GLASS_d$g1,
                  GLASS_d$g2,
                  FEATURE_PVAL,tissue_PVAL,MSI_PVAL,
                  FEATURE_IC50_WTT_pvalue),1,21,
                dimnames=list('1',c('FEATURE',
                                    'Drug id',
                                    'Owned_by',
                                    'Drug name',
                                    'Drug Target',
                                    'N_FEATURE_pos','N_FEATURE_neg',
                                    'log max.Conc.tested','log max.Conc.tested2',                                 
                                    'FEATUREpos_logIC50_MEAN',
                                    'FEATUREneg_logIC50_MEAN',
                                    'FEATURE_deltaMEAN_IC50',                                                           
                                    'FEATUREpos_IC50_sd',
                                    'FEATUREneg_IC50_sd',
                                    'FEATURE_IC50_effect_size',
                                    'FEATUREpos_Glass_delta',
                                    'FEATUREneg_Glass_delta',
                                    'FEATURE_ANOVA_pval','Tissue_ANOVA_pval','MSI_ANOVA_pval',                                 
                                    'FEATURE_IC50_T_pval')))
  } else {
    Npos<-length(which(FEATpattern=='pos'))
    Nneg<-length(which(FEATpattern=='neg'))
    RES<-matrix(c(FEATURE,
                  DRUG_ID,
                  DRUG_PROPS[str_split(DRUG_ID,'_')[[1]][1],"OWNED_BY"],
                  as.character(DRUG_PROPS[DRUG_ID,"DRUG_NAME"]),
                  as.character(DRUG_PROPS[DRUG_ID,"PUTATIVE_TARGET"]),Npos,Nneg,rep(NA,14)),
                1,21,dimnames=list('1',c('FEATURE',
                                         'Drug id',
                                         'Owned by',
                                         'Drug name',
                                         'Drug Target',
                                         'N_FEATURE_pos','N_FEATURE_neg',
                                         'log max.Conc.tested','log max.Conc.tested2',                                 
                                         'FEATUREpos_logIC50_MEAN',
                                         'FEATUREneg_logIC50_MEAN',
                                         'FEATURE_deltaMEAN_IC50',                                                           
                                         'FEATUREpos_IC50_sd',
                                         'FEATUREneg_IC50_sd',
                                         'FEATURE_IC50_effect_size',
                                         'FEATUREpos_Glass_delta',
                                         'FEATUREneg_Glass_delta',
                                         'FEATURE_ANOVA_pval','Tissue_ANOVA_pval','MSI_ANOVA_pval',                                 
                                         'FEATURE_IC50_T_pval')))
  }
  
  
  
  if (display){
    drug_id<-str_split(DRUG_ID,'_')[[1]][1]
    mtext(paste(DRUG_ID,' ',DRUG_PROPS[drug_id,"DRUG_NAME"],' [',DRUG_PROPS[drug_id,"PUTATIVE_TARGET"],'] \n ',FEATURE,sep=''),
          side = 3, line = -4, outer = TRUE,cex=1.5)
    
    plot(0,0,col=NA,xaxt='n',yaxt='n',frame.plot = FALSE,xlab='',ylab='')
    legend('top',c('sample with negative feature','sample with positive feature'),pch=16,col=c(GRAY,BLUE),pt.cex = 2.5,cex=2,inset = c(0,0.2),bty = 'n')
    legend('center',c('mean','mean +/- sd/2','mean +/- sd','median'),lty = c(1,1,2,1),col=c('red','purple','purple','black'),lwd=c(4,4,4,6),cex=2,bty = 'n')
    legend('bottom',c('- Bottom and Top of the boxes indicates 25 and 75 percentiles','- Whiskers indicates min and max values',
                      '[exc. outliers, i.e. more (resp. less) 2/3 upper (resp. lower) 25 percentile]'),xjust = 0,inset = c(0,0.25),bty = 'n')
    
    if (printTOfig){
      dev.off()
    }
    
    if(GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES=='PANCAN'){
      if (printTOfig){
        png(paste(PATH,FN,'_01_TissueWisk.png',sep=''),width=793.92,height=1122.24)
      }
      
      gdscANOVA_whiskerPlots(IC50pattern,FEATpattern,TISSUEpattern,DRUG_ID=DRUG_ID,FEATURE=FEATURE,diciture='Cancer-Type')
      
      if (printTOfig){
        dev.off()
      }
    }
    if(GDSCANOVA_SETTINGS$gdscANOVA.settings.includeMSI_Factor){
      
      MSIpatternLit<-rep('MSI-stable',length(MSIpattern))
      names(MSIpatternLit)<-names(MSIpattern)
      
      MSIpatternLit[which(MSIpattern==1)]<-'MSI-instable'
      MSIpatternLit<-as.factor(MSIpatternLit)
      
      if (printTOfig){
        png(paste(PATH,FN,'_02_MSIWisk.png',sep=''),width=793.92,height=1122.24)
      }
      
      gdscANOVA_whiskerPlots(IC50pattern,FEATpattern,MSIpatternLit,DRUG_ID=DRUG_ID,FEATURE=FEATURE,diciture='MS-instability')
      
      if (printTOfig){
        dev.off()
      }
    }
  }
  
  return(RES)  
}

gdscANOVA_singleDrugANOVA<-function(DRUG_ID,verbose=TRUE){
  
  FEATURES<-names(which(rowSums(InputFeatures$BEM)>=3))
  nFEATURES<-length(FEATURES)
  
  if(verbose){
    drug_id<-str_split(DRUG_ID,'_')[[1]][1]
    print(paste('Running Single Drug ANOVA:',DRUG_ID,' - ',
                DRUG_PROPS[drug_id,"DRUG_NAME"],' [',DRUG_PROPS[drug_id,"PUTATIVE_TARGET"],']',sep=''))
    pb <- txtProgressBar(min=1,max=nFEATURES,style=3)
  }
  
  for (i in 1:nFEATURES){
    
    if(verbose){
      setTxtProgressBar(pb, i)
    }
    if (i ==1){
      TOTRES<-gdscANOVA_individualANOVA(DRUG_ID=DRUG_ID,FEATURE=FEATURES[i],display=FALSE)
    } else{
      TOTRES<-rbind(TOTRES,gdscANOVA_individualANOVA(DRUG_ID=DRUG_ID,FEATURE=FEATURES[i],display=FALSE))
    }
  }
  
  idxs<-which(!is.na(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])))
  
  TOTRES<-TOTRES[idxs,]
  
  if (length(idxs)==1){
    TOTRES<-matrix(TOTRES,1,length(TOTRES),dimnames=list(1,names(TOTRES)))
  }
  
  
  if(verbose){
    Sys.sleep(1)
    close(pb)
  }
  
  
  
  return(TOTRES)
}
gdscANOVA_totalANOVA<-function(fn){
  
  DRUGS<-names(which(colSums(!is.na(IC50s))>=6))
  
  nDRUGS<-length(DRUGS)
  
  print('+ Running ANOVA')
  pb <- txtProgressBar(min=1,max=nDRUGS,style=3)
  
  flag<-1
  
  for (i in 1:nDRUGS){
    
    setTxtProgressBar(pb, i)
    currentRES<-gdscANOVA_singleDrugANOVA(DRUGS[i],verbose=FALSE)
    
    if (flag == 1 & ncol(currentRES)>0){
      TOTRES<-currentRES
      flag<-0
    } else{
      if (ncol(currentRES)>0){
        TOTRES<-rbind(TOTRES,currentRES)
      }
    }
    if (i%%50==0){
      write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(fn,'.txt',sep=''))
    }
  }
  
  TOTRES<-TOTRES[order(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])),]
  
  if(GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_correction_method!='qvalue'){
    FDR<-p.adjust(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"]),method='fdr')
  }else{
    Q<-qvalue(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"]))
    FDR<-Q$qvalue
  }
  
  TOTRES<-cbind(TOTRES,FDR*100)
  colnames(TOTRES)[ncol(TOTRES)]<-'ANOVA FEATURE FDR %'
  
  assoc_id<-paste('a',1:nrow(TOTRES),sep='')
  
  TOTRES<-cbind(assoc_id,TOTRES)
  
  colnames(TOTRES)[1]<-'assoc_id'
  
  write.table(TOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(fn,'.txt',sep=''))
  
  save(TOTRES,file=paste(fn,'.rdata',sep=''))
  
  if(length(range)>1){
    Sys.sleep(1)
    close(pb)
  }
  
  return(TOTRES)
}



