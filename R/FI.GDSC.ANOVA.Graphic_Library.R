## graphics
BLUE<-rgb(0,0,255,100,maxColorValue=255)
GRAY<-rgb(180,180,180,180,maxColorValue=255)



gdscANOVA_createSCATTERS<-function(range,PATH){
  print('Creating idividual association plots')
  
  if (length(range)==1){M<-2}
  else{M<-max(range)}
  pb <- txtProgressBar(min=1,max=M,style=3)
  for (i in range){
    setTxtProgressBar(pb,i)
    gdscANOVA_individualANOVA(DRUG_ID=TOTRES[i,'Drug id'],FEATURE=TOTRES[i,'FEATURE'],PATH=PATH,
                              printTOfig=TRUE,FN=TOTRES[i,'assoc_id'],FDR=as.numeric(TOTRES[i,"ANOVA FEATURE FDR %"]))
    
  }
  
  Sys.sleep(1)
  close(pb)
  
  print('done!')
}
gdscANOVA_scatterSets<-function(IC50pattern,MUTpattern,DRUG_ID,FEATURE){
  
  par(mfrow=c(1,2))
  
  cols<-rep(GRAY,length(IC50pattern))
  cols[which(MUTpattern=='pos')]<-BLUE      
  
  orDRUG_ID<-DRUG_ID
  
  DRUG_ID<-str_split(DRUG_ID,'_')[[1]][1]
  par(mar=c(4,6,5,4))
  
  boxplot(IC50pattern~MUTpattern, col=NA,frame.plot=FALSE,boxwex=0.5,outline=FALSE,lwd=2,cex.axis=2,cex.lab=2,
          ylim=c(min(IC50pattern),max(IC50pattern)),ylab=paste(DRUG_PROPS[DRUG_ID,'DRUG_NAME'],'log(IC50)'))
  par(new=TRUE)
  beeswarm(IC50pattern~MUTpattern,ylim=c(min(IC50pattern),max(IC50pattern)),col=c(GRAY,BLUE),pch=16,xaxt='n',
           yaxt='n',xlab='',ylab='',corral='wrap',cex=3)
  par(new=TRUE)
  boxplot(IC50pattern~MUTpattern, col=NA,frame.plot=FALSE,boxwex=0.5,outline=FALSE,lwd=2,names=c('',''),yaxt='n',xaxt='n',
          ylim=c(min(IC50pattern),max(IC50pattern)))
  
  #DRC<-gdscANOVA_retrieveDScurve(names(IC50pattern),drug_id=DRUG_ID)
  
  
  
  umax<-as.numeric(unique(maxConcTested[,orDRUG_ID]))
  umax<-log(umax[!is.na(umax)])
  

  for (i in 1:length(umax)){
    abline(h=umax[i],lty=2,lwd=1)  
  }
  
  Mpos<-mean(IC50pattern[MUTpattern=='pos'])
  Mneg<-mean(IC50pattern[MUTpattern=='neg'])
  
  SDpos<-sd(IC50pattern[MUTpattern=='pos'])
  SDneg<-sd(IC50pattern[MUTpattern=='neg'])
  
  lines(x=c(1.85,2.15),y=c(Mpos,Mpos),col='red',lwd=4)
  lines(x=c(0.85,1.15),y=c(Mneg,Mneg),col='red',lwd=4)
  
  lines(x=c(1.80,2.20),y=c(Mpos+SDpos/2,Mpos+SDpos/2),col='purple',lwd=4)
  lines(x=c(1.80,2.20),y=c(Mpos-SDpos/2,Mpos-SDpos/2),col='purple',lwd=4)
  
  lines(x=c(1.80,2.20),y=c(Mpos+SDpos,Mpos+SDpos),col='purple',lwd=2,lty=2)
  lines(x=c(1.80,2.20),y=c(Mpos-SDpos,Mpos-SDpos),col='purple',lwd=2,lty=2)
  
  lines(x=c(0.80,1.20),y=c(Mneg+SDneg/2,Mneg+SDneg/2),col='purple',lwd=4)
  lines(x=c(0.80,1.20),y=c(Mneg-SDneg/2,Mneg-SDneg/2),col='purple',lwd=4)
  
  lines(x=c(0.80,1.20),y=c(Mneg+SDneg,Mneg+SDneg),col='purple',lwd=2,lty=2)
  lines(x=c(0.80,1.20),y=c(Mneg-SDneg,Mneg-SDneg),col='purple',lwd=2,lty=2)
  
}
gdscANOVA_whiskerPlots<-function(IC50pattern,MUTpattern,TISSUEpattern,DRUG_ID,FEATURE,diciture){
  
  
  M<-tapply(IC50pattern,list(MUTpattern,TISSUEpattern),'mean')
  N<-tapply(IC50pattern,list(MUTpattern,TISSUEpattern),'length')
  
  id<-which(colSums(N>1)==2)
  
  if (length(id)>1){
    M<-M[,id]
    N<-N[,id]
  } else{
    M<-matrix(M[,id],nrow=2,ncol=1,dimnames=list(c('pos','neg'),colnames(M)[id]))
    N<-matrix(N[,id],nrow=2,ncol=1,dimnames=list(c('pos','neg'),colnames(N)[id]))
  }
  
  
  delta<-M[1,]-M[2,]
  
  if (length(id)>1){
    M<-M[,order(delta)]
    N<-N[,order(delta)]
  }
  
  rem_tissues<-colnames(N)
  nrt<-length(rem_tissues)
  
  TP<-rep(NA,nrt)
  names(TP)<-rem_tissues
  
  if (length(TP)>0){
    for (i in 1:length(TP)){
      
      cid<-which(TISSUEpattern==rem_tissues[i])
      cIC50p<-IC50pattern[cid]
      cMUTp<-MUTpattern[cid]
      
      
      
      res<-t.test(cIC50p~cMUTp)
      
      TP[i]<-res$p.value
    }
    
    
    l1<-rownames(M)
    l2<-colnames(M)
    
    factors<-rep(NA,length(l1)*length(l2))
    flag<-0
    
    for (j in 1:length(l2)){
      for (i in 1:length(l1)){
        factors[MUTpattern==l1[i] & TISSUEpattern==l2[j]]<-flag  
        flag<-flag+1
      }
    }
    
    par(mar=c(5,8,4,2))
    par(las=2)
    par(fig=c(0,1,0,1))
    
    cols<-rep(BLUE,length(l2)*2)
    cols[seq(1,length(l2)*2-1,2)]<-GRAY
    
    if(length(cols)==2){
      cols<-cols[2:1]
    }
    
    par(mar=c(5,15,4,4))
    
    #   l2[seq(1,13,2)]<-paste('mut',l2[seq(1,13,2)],sep=' / ')
    #   l2[seq(2,14,2)]<-paste('wt',l2[seq(2,14,2)],sep=' / ')
    
    significance<-rep('',length(TP))
    
    names(significance)<-names(TP)
    
    significance[which(TP<0.05)]<-'*'
    significance[which(TP<0.01)]<-'**'
    significance[which(TP<0.001)]<-'***'
    
    NAMES<-c(rbind(paste(l2,l1[1]),paste(l2,l1[2])))
    
    if (length(NAMES)<=2){
      NAMES[seq(1,length(NAMES),2)]<-paste(significance,NAMES[seq(1,length(NAMES),2)])
    } else{
      NAMES[seq(2,length(NAMES),2)]<-paste(significance,NAMES[seq(2,length(NAMES),2)])
    }
    IC50pattern<-IC50pattern[which(!is.na(factors))]
    factors<-factors[which(!is.na(factors))]
    
    
    boxplot(IC50pattern~factors,horizontal=TRUE,names=NAMES,col=cols,cex=2,cex.axis=1.2,cex.main=2,cex.lab=1.2,
            cex.axis=0.6,xlab=paste(DRUG_ID,DRUG_PROPS[DRUG_ID,"DRUG_NAME"],'logIC50'),
            ylim=c(min(IC50pattern),max(IC50pattern)),xaxt='n',main=paste('FEATURE/',diciture,' interactions',sep=''))
    
    
    axis(1)
    
    par(new=TRUE)
    beeswarm(IC50pattern~factors,horizontal=TRUE,col=cols,corral='gutter',pch=16,cex=2,
             cex.axis=0.6,xlab='',xlim=c(min(IC50pattern),max(IC50pattern)),
             xaxt='n',yaxt='n',ylab='')
    
    count<-tapply(IC50pattern,factors,'length')
    
    axis(4,at=1:length(unique(factors)),count)
    
    for (i in 1:length(unique(factors))){
      abline(h=i,col='gray',lty=2)
    }
    
    
    
    legend('topleft',inset=c(0.25,0),bty='n',fill=c(GRAY,BLUE),y.intersp=1.5,
           legend=c(paste(FEATURE,'negative'),
                    paste(FEATURE,'positive')))
    
    legend('topleft',bty='n',y.intersp=1.5,
           legend=c('* p < 0.05', '** p < 0.01','*** p < 0.001'))
    
  }
}
gdscANOVA_volcanoPlot_T<-function(delta,pval,qvals,N,
                                  minN=NA,maxN=NA,
                                  fdrth=GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH,effth=1,LIM=Inf,pointLabels=NULL,fdr=NULL,main='',fontsize=1){
  
  
  par(mar=c(5,7,4,4))
  
  nullcol<-rgb(0,0,0,25,maxColorValue=255)
  redcol<-rgb(255,0,0,100,maxColorValue=255)
  greencol<-rgb(0,255,0,100,maxColorValue=255)
  
  COL<-rep(nullcol,length(pval))
  
  if(GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES=='PANCAN'){
    COL[which(qvals<=fdrth & delta<=-effth)]<-greencol
    COL[which(qvals<=fdrth & delta>=effth)]<-redcol  
  }else{
    COL[which(qvals<=fdrth & pval<=GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH & delta>0)]<-redcol
    COL[which(qvals<=fdrth & pval<=GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH & delta<0)]<-greencol
  }
  
  
  if(is.na(minN) & is.na(maxN)){
    pwsizes<-5*(N-min(N))/(max(N)-min(N))+1
  }else{
    pwsizes<-5*(N-minN)/(maxN-minN)+1
  }
  
  
  ii<-which(-log10(pval)<LIM)
  
  XL<-max(as.numeric(abs(delta))+1)
  YL<-max(-log10(pval[ii]))
  
  plot(delta[ii],-log10(pval[ii]),xlim=c(-XL-0.5,XL+0.5),ylim=c(0,YL+5),yaxt='n',yaxs='i',pch=16,frame.plot=FALSE,
       col=COL[ii],xlab='signed effect size',ylab='-log10 (p-value)',cex=pwsizes[ii],main=main,cex.axis=2,cex.lab=2)
  abline(v=0,col='black')
  
  if(length(fdr)==0){
    fdrLim20<-max(pval[which(qvals<GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)]) 
    fdrLim1<-max(pval[which(qvals<10)]) 
    fdrLim001<-max(pval[which(qvals<1)]) 
  } else{
    fdrLim20<-fdr[1]
    fdrLim1<-fdr[2]
    fdrLim001<-fdr[3]
  }
  
  abline(h=-log10(0.001),col='darkgray',lty=5,lwd=1.5)
  
  abline(h=-log10(fdrLim20),col='darkgray',lty=2,lwd=1.5)
  
  abline(h=-log10(fdrLim1),col='darkgray',lty=3,lwd=1.5)
  abline(h=-log10(fdrLim001),col='darkgray',lty=4,lwd=1.5)
  
  #   text(-XL-0.5,-log10(0.05)+0.05,'p=0.05',cex=0.7,col='darkgray')  
  #   text(-XL-0.5,-log10(fdrLim20)+0.05,'FDR 20%',cex=0.7,col='gray')
  #   text(-XL-0.5,-log10(fdrLim1)+0.05,'FDR 1%',cex=0.7,col='gray')
  #   text(-XL-0.5,-log10(fdrLim001)+0.05,'FDR 0.01%',cex=0.7,col='gray')
  
  if (GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES=='PANCAN'){
    legend('topleft',legend=(c('FDR 1%','FDR 10%',paste('FDR ',GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH,'%',sep=''),'p 0.001')),cex=1,lty=c(4,3,2,5),lwd=1.5,y.intersp=1.5,
           col='gray',bty='n')  
  }else{
    legend('topleft',legend=(c('FDR 1%','FDR 10%',paste('FDR ',GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH,'%',sep=''),'p 0.001')),cex=1,lty=c(4,3,2,5),lwd=1.5,y.intersp=1.5,
           col='gray',bty='n')
  }
  
  
  if(effth>0){
    legend('topleft',inset=c(0.7,0),legend=(paste('| effect size | =',effth)),cex=1,lty=1,lwd=2,col='gray',bty='n',xjust=0)
  }
  
  abline(v=-effth,col='lightgray',lty=1,lwd=2)
  abline(v=effth,col='lightgray',lty=1,lwd=2)
  
  axis(2,las=1,col='gray')
  axis(2, at=c(0,max(-log10(pval))+10), labels=c("",""), lwd.ticks=0)
  axis(2, at=seq(0 , max(-log10(pval)+10), by=5), lwd=0, lwd.ticks=1,las=1,col='gray')
  
  
  if (length(pointLabels)>0){
    id<-which(abs(delta)>effth & pval < fdrLim20) 
    if (length(id)>0){
      text(delta[id],-log10(pval[id]),pointLabels[id],pos=3,cex=fontsize)
    }
    
  } 
}
gdscANOVA_allDRUGS_volcanoPlots<-function(TOTRES,PATH=''){
  
  print('- Generating Individual Drug Volcano plots')
  DRUG_IDS<-unique(TOTRES[,'Drug id'])
  
  nDRUGS<-length(DRUG_IDS)
  pb <- txtProgressBar(min=1,max=nDRUGS,style=3)
  for (i in 1:nDRUGS){
    setTxtProgressBar(pb,i)
    gdscANOVA_DRUG_volcanoPlot(TOTRES,DRUG_IDS[i],print=TRUE,PATH=PATH)
  }
  
  Sys.sleep(1)
  close(pb)
  print('Done!')
}
gdscANOVA_DRUG_volcanoPlot<-function(TOTRES,drug_id,print=FALSE,PATH=''){
  id<-which(TOTRES[,"Drug id"]==drug_id)
  
  minN<-min(as.numeric(TOTRES[,"N_FEATURE_pos"]))
  maxN<-max(as.numeric(TOTRES[,"N_FEATURE_pos"]))
  
  qvals<-as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])
  pvals<-as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])
  
  fdr<-c(max(pvals[which(qvals<GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)]),max(pvals[which(qvals<1)]),max(pvals[which(qvals<0.01)]))
  
  if(length(id)>1){
    TOTRES<-TOTRES[id,]
  }else{
    TOTRES<-matrix(TOTRES[id,],1,ncol(TOTRES),dimnames=list(1,colnames(TOTRES)))
  }
  labels<-TOTRES[,'FEATURE']
  
  labels[which(as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])>=GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)]<-''
  MAIN<-paste(drug_id,' - ',DRUG_PROPS[str_split(drug_id,'_')[[1]][1],"DRUG_NAME"],' [',DRUG_PROPS[str_split(drug_id,'_')[[1]][1],"PUTATIVE_TARGET"],']',sep='')
  
  if(print){
    FN<-paste(PATH,str_replace(drug_id,pattern = '/',replacement = '_'),'.png',sep='')
    png(FN,768,1024)
  }
  
  delta<-sign(as.numeric(TOTRES[,"FEATURE_deltaMEAN_IC50"]))*as.numeric(TOTRES[,"FEATURE_IC50_effect_size"])
  pval<-as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])
  qval<-as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])
  N<-as.numeric(TOTRES[,"N_FEATURE_pos"])
  
  gdscANOVA_volcanoPlot_T(delta=delta,fontsize=1,
                          pval=pval,
                          qvals=qval,
                          N=N,minN=minN,maxN=maxN,effth=0,fdr=fdr,main=MAIN,
                          pointLabels=labels)
  
  if(print){
    dev.off()
  }
  
  
}
gdscANOVA_allFEATURES_volcanoPlots<-function(TOTRES,PATH=''){
  print('- Generating Individual Feature Volcano plots')
  FEATURES<-unique(TOTRES[,"FEATURE"])
  nFEATURES<-length(FEATURES)
  
  if (nFEATURES>1){
    pb <- txtProgressBar(min=1,max=nFEATURES,style=3)
  }
  for (i in 1:nFEATURES){
    if (nFEATURES>1){setTxtProgressBar(pb,i)}
    gdscANOVA_FEATURE_volcanoPlot(FEATURES[i],cTOTRES=TOTRES,print=TRUE,PATH=PATH)
  }
  Sys.sleep(1)
  if (nFEATURES>1){close(pb)}
}
gdscANOVA_FEATURE_volcanoPlot<-function(FEATURE,print=FALSE,cTOTRES,PATH=''){
  
  id<-which(cTOTRES[,"FEATURE"]==FEATURE)
  
  
  minN<-min(as.numeric(cTOTRES[,"N_FEATURE_pos"]))
  maxN<-max(as.numeric(cTOTRES[,"N_FEATURE_pos"]))
  
  qvals<-as.numeric(cTOTRES[,"ANOVA FEATURE FDR %"])
  pvals<-as.numeric(cTOTRES[,"FEATURE_ANOVA_pval"])
  
  fdr<-c(max(pvals[which(qvals<GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)]),max(pvals[which(qvals<1)]),fdrLim001<-max(pvals[which(qvals<0.01)]))
  
  if (length(id)==1){cTOTRES<-t(as.matrix(cTOTRES[id,],1,ncol(cTOTRES)))}
  else{cTOTRES<-cTOTRES[id,]}
  
  MAIN<-FEATURE
  if(print){
    FN<-paste(PATH,FEATURE,'.png',sep='')
    png(FN,768,1024)
  }
  
  labels<-paste(cTOTRES[,"Drug id"],'-',cTOTRES[,"Drug name"])
  
  labels[which(as.numeric(cTOTRES[,"ANOVA FEATURE FDR %"])>=GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)]<-''
  MAIN<-FEATURE
  
  
  delta<-sign(as.numeric(cTOTRES[,"FEATURE_deltaMEAN_IC50"]))*as.numeric(cTOTRES[,"FEATURE_IC50_effect_size"])
  pval<-as.numeric(cTOTRES[,"FEATURE_ANOVA_pval"])
  qval<-as.numeric(cTOTRES[,"ANOVA FEATURE FDR %"])
  N<-as.numeric(cTOTRES[,"N_FEATURE_pos"])
  
  gdscANOVA_volcanoPlot_T(delta=delta,
                          fontsize=1,
                          fdrth=GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH,
                          pval=pval,
                          qvals=qval,
                          N=N,
                          minN=minN,
                          maxN=maxN,
                          effth=0,
                          fdr=fdr,
                          main=MAIN,
                          pointLabels=labels)
  
  if(print){
    dev.off()
  }
  
}

## result summaries
gdscANOVA_computeFeatureStats<-function(redTOTRES,fdrTH=30,pvalTH=Inf,save=TRUE,display=TRUE,printTOfig=TRUE,PATH=NULL){
  
  print('- Computing feature-wise associations summary')
  
  idSENS<-which(as.numeric(redTOTRES[,'ANOVA FEATURE FDR %'])<fdrTH & as.numeric(redTOTRES[,"FEATURE_ANOVA_pval"])<pvalTH & as.numeric(redTOTRES[,'FEATURE_deltaMEAN_IC50'])<0)
  idRES<-which(as.numeric(redTOTRES[,'ANOVA FEATURE FDR %'])<fdrTH & as.numeric(redTOTRES[,"FEATURE_ANOVA_pval"])<pvalTH & as.numeric(redTOTRES[,'FEATURE_deltaMEAN_IC50'])>0)
  
  if(length(idSENS)>1){
    sensredTOTRES<-redTOTRES[idSENS,]
  }else{
    sensredTOTRES<-matrix(redTOTRES[idSENS,],1,ncol(redTOTRES),dimnames=list(1,colnames(redTOTRES)))
  }
  
  if(length(idRES)>1){
    resredTOTRES<-redTOTRES[idRES,]
  }else{
    resredTOTRES<-matrix(redTOTRES[idRES,],1,ncol(redTOTRES),dimnames=list(1,colnames(redTOTRES)))
  }
  
  tmp<-as.factor(sensredTOTRES[,'FEATURE'])
  if (!is.na(tmp)){
    sensitizations<-summary(tmp,maxsum=length(levels(as.factor(sensredTOTRES[,'FEATURE']))))
  }else{
    sensitizations<-NULL
  }
  tmp<-as.factor(resredTOTRES[,'FEATURE'])
  if (!is.na(tmp)){
    resistizations<-summary(tmp,maxsum=length(levels(as.factor(resredTOTRES[,'FEATURE']))))
  }else{
    resistizations<-NULL
  }
  
  if(length(idRES)>0 & length(idSENS)>0){
    gdomain<-union(names(sensitizations),names(resistizations))
  } else{
    if (length(idRES)>0){
      gdomain<-names(resistizations)
    }else{
      if (length(idSENS)>0){
        gdomain<-names(sensitizations)
      } else{gdomain<-NULL}
    }
  }
  
  
  totalgdomain<-unique(redTOTRES[,'FEATURE'])
  
  Gperc<-length(gdomain)/length(totalgdomain)
  
  S<-rep(0,length(gdomain))
  names(S)<-sort(gdomain)
  
  if (length(idSENS)>0){
    S[names(sensitizations)]<-sensitizations
  }
  
  R<-rep(0,length(gdomain))
  names(R)<-sort(gdomain)
  
  if (length(idRES)>0){
    R[names(resistizations)]<-resistizations
  }
  
  TOTAL<-cbind(S,R,S+R)
  
  if(length(gdomain)>1){
    TOTAL<-TOTAL[order(TOTAL[,3],decreasing=TRUE),]
  }
  
  colnames(TOTAL)<-c('sens assoc','res assoc','tot')
  
  gnames<-unique(rownames(TOTAL))
  
  
  if (nrow(TOTAL)>1){
    n_altered_samples<-rowSums(InputFeatures$BEM[rownames(TOTAL),])
  }else{
    n_altered_samples<-sum(InputFeatures$BEM[rownames(TOTAL),])
  }
  
  TOTAL<-cbind(n_altered_samples,TOTAL)
  
  if(nrow(TOTAL)>1){
    TOTAL<-TOTAL[,c(1,4,2,3)]
  } else{
    TOTAL<-matrix(TOTAL[,c(1,4,2,3)],1,ncol(TOTAL),dimnames=list(rownames(TOTAL),colnames(TOTAL))) 
  }
  
  
  if (display){
    if(printTOfig){
      if(length(PATH)==0){
        png(file=paste(current_dir,'OUTPUT','/GRAPHICS/Top50FeatureSummary.png',sep=''),width=1300,height=2000)
      }else{
        png(file=paste(PATH,'.png',sep=''),width=1300,height=2000)
      }
    }
    par(mar=c(9,52,8,4))
    
    n<-nrow(TOTAL)
    
    if (n > 50){
      n<-50
    }
    if (n>1){
      subTOTAL<-TOTAL[1:n,]  
    }else{
      subTOTAL<-TOTAL
    }
    
    par(xpd=TRUE)
    barplot(t(subTOTAL[seq(nrow(subTOTAL),1,-1),c("sens assoc","res assoc")]),horiz=TRUE,las=2,col=c('purple','orange'),cex.names=2,cex.main=3,
            cex.axis=2,xlab='',cex.lab=2,
            border=FALSE,xlim=c(0,max(as.numeric(subTOTAL[,"sens assoc"])+as.numeric(subTOTAL[,"res assoc"]))+1),
            main=paste('Top-50 features most frequently\nassociated with drug response'))
    legend('bottomright',c('sensitivity','resistance'),fill=c('purple','orange'),cex=3,border=FALSE,bty='n')
    
    text((max(as.numeric(subTOTAL[,"sens assoc"])+as.numeric(subTOTAL[,"res assoc"]))+10)/2,-5.5,'n. significant associations',cex=3)
    if(printTOfig){
      dev.off()
    }
  }
  
  if (save){
    TOTAL<-cbind(rownames(TOTAL),TOTAL)
    colnames(TOTAL)[1]<-'Feature'
    colnames(TOTAL)[5]<-'res assoc'
    
    if(length(PATH)==0){
      write.table(TOTAL,quote=FALSE,row.names=F,sep='\t',file=paste(current_dir,'OUTPUT','/FeatureSummary.txt',sep=''))
    }else{
      write.table(TOTAL,quote=FALSE,row.names=F,sep='\t',file=paste(PATH,'.txt',sep=''))
    }
    FeatureSummary<-TOTAL
    if(length(PATH)==0){
      save(FeatureSummary,file=paste(current_dir,'OUTPUT','/FeatureSummary.rdata',sep=''))
    } else{
      save(FeatureSummary,file=paste(PATH,'.rdata',sep=''))
    }
  }
  print('+ done!')
  return(TOTAL)
}
gdscANOVA_computeDrugStats<-function(redTOTRES,fdrTH=30,pvalTH=Inf,save=TRUE,display=TRUE,printTOfig=TRUE,PATH=NULL){
  print('- Computing drug-wise associations summary')
  
  
  idSENS<-which(as.numeric(redTOTRES[,'ANOVA FEATURE FDR %'])<fdrTH & as.numeric(redTOTRES[,"FEATURE_ANOVA_pval"])<pvalTH & as.numeric(redTOTRES[,'FEATURE_deltaMEAN_IC50'])<0)
  idRES<-which(as.numeric(redTOTRES[,'ANOVA FEATURE FDR %'])<fdrTH & as.numeric(redTOTRES[,"FEATURE_ANOVA_pval"])<pvalTH & as.numeric(redTOTRES[,'FEATURE_deltaMEAN_IC50'])>0)
  
  
  if(length(idSENS)>1){
    sensredTOTRES<-redTOTRES[idSENS,]
  }else{
    if(length(idSENS)==1){
      sensredTOTRES<-matrix(redTOTRES[idSENS,],1,ncol(redTOTRES),dimnames=list(1,colnames(redTOTRES)))
    }else{
      sensredTOTRES<-NULL
    }
  }
  
  if(length(idRES)>1){
    resredTOTRES<-redTOTRES[idRES,]
  }else{
    if(length(idRES)==1){
      resredTOTRES<-matrix(redTOTRES[idRES,],1,ncol(redTOTRES),dimnames=list(1,colnames(redTOTRES)))  
    }else{
      resredTOTRES<-NULL
    }
  }
  
  
  if(length(sensredTOTRES)>0){
    sensitizations<-summary(as.factor(sensredTOTRES[,"Drug id"]),
                            maxsum=length(levels(as.factor(sensredTOTRES[,"Drug id"]))))  
  }else{
    sensitizations<-NULL
  }
  
  if(length(resredTOTRES)>0){
    resistizations<-summary(as.factor(resredTOTRES[,"Drug id"]),
                            maxsum=length(levels(as.factor(resredTOTRES[,"Drug id"]))))  
  }else{
    resistizations<-NULL
  }
  
  DrugDomain<-union(names(sensitizations),names(resistizations))
  
  totalddomain<-unique(redTOTRES[,"Drug id"])
  
  S<-rep(0,length(DrugDomain))
  names(S)<-sort(DrugDomain)
  
  S[names(sensitizations)]<-sensitizations
  
  R<-rep(0,length(DrugDomain))
  names(R)<-sort(DrugDomain)
  
  R[names(resistizations)]<-resistizations
  
  TOTAL<-cbind(S,R,S+R)
  
  if(nrow(TOTAL)>1){
    TOTAL<-TOTAL[order(TOTAL[,3],decreasing=TRUE),]
  }
  
  colnames(TOTAL)<-c('sens assoc','res assoc','tot')
  
  DRUGS<-paste(as.character(rownames(TOTAL)),' - ',DRUG_PROPS[unlist(str_split(as.character(rownames(TOTAL)),'_'))[seq(1,2*nrow(TOTAL),2)],'DRUG_NAME'],
               ' [',DRUG_PROPS[unlist(str_split(as.character(rownames(TOTAL)),'_'))[seq(1,2*nrow(TOTAL),2)],'PUTATIVE_TARGET'],']',sep='')
  
  TOTAL<-cbind(DRUGS,TOTAL)
  
  if (display){
    if(printTOfig){
      if(length(PATH)==0){
        png(file=paste(current_dir,'OUTPUT','/GRAPHICS/Top50DrugSummary.png',sep=''),width=2000,height=2000)
      } else{
        png(file=paste(PATH,'.png',sep=''),width=2000,height=2000)
      }
    }
    par(mar=c(9,62,8,4))
    
    n<-nrow(TOTAL)
    
    if (n > 50){
      n<-50
    }
    subTOTAL<-TOTAL[seq(n,1,-1),]
    
    par(xpd=TRUE)
    
    barplot(t(subTOTAL[,c("sens assoc","res assoc")]),horiz=TRUE,las=2,col=c('purple','orange'),cex.names=2,cex.main=3,
            cex.axis=2,xlab='',cex.lab=2,names.arg=subTOTAL[1:n,"DRUGS"],
            border=FALSE,xlim=c(0,max(as.numeric(subTOTAL[,"sens assoc"])+as.numeric(subTOTAL[,"res assoc"]))),
            main=paste('Top-50 drugs whose response\n is frequently associated with a feature'))
    legend('bottomright',c('sensitivity','resistance'),fill=c('purple','orange'),cex=3,border=FALSE,bty='n')
    
    text((max(as.numeric(subTOTAL[,"sens assoc"])+as.numeric(subTOTAL[,"res assoc"]))+3)/2,-5.5,'n. significant associations',cex=3)
    if(printTOfig){
      dev.off()
    }
  }
  
  if (save){
    
    if(length(PATH)==0){
      write.table(TOTAL,quote=FALSE,row.names=F,sep='\t',file=paste(current_dir,'OUTPUT','/DrugSummary.txt',sep=''))
    }else{
      write.table(TOTAL,quote=FALSE,row.names=F,sep='\t',file=paste(PATH,'.txt',sep=''))
    }
    
    DrugSummary<-TOTAL
    if(length(PATH)==0){
      save(DrugSummary,file=paste(current_dir,'OUTPUT','/DrugSummary.rdata',sep=''))
    }else{
      save(DrugSummary,file=paste(PATH,'.rdata',sep=''))
    }
  }
  
  print('+ done!')
  return(TOTAL)
}
gdscANOVA_howManySignHits<-function(TOTRES,YLIM=2000,MAIN=''){
  
  FDRranges<-seq(5,50,5)
  
  nths<-length(FDRranges)
  
  significantHits<-rep(0,nths)
  significantMeaningfullHits<-rep(0,nths)
  significantMeaningfullStrongHits<-rep(0,nths)
  significantMeaningfullSuperStrongHits<-rep(0,nths)
  
  names(significantHits)<-FDRranges
  names(significantMeaningfullHits)<-FDRranges
  names(significantMeaningfullStrongHits)<-FDRranges
  names(significantMeaningfullSuperStrongHits)<-FDRranges
  
  for (i in 1:nths){
    significantHits[i]<-length(which(as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<FDRranges[i]))
    
    MC1<-as.numeric(TOTRES[,"log max.Conc.tested"])
    MC2<-as.numeric(TOTRES[,"log max.Conc.tested2"])
    MC2[is.na(MC2)]<- -Inf
    
    idxs<-which(as.numeric(TOTRES[,"ANOVA FEATURE FDR %"]) < FDRranges[i] & (as.numeric(TOTRES[,"FEATUREpos_logIC50_MEAN"])< MC1 |
                                                                               (as.numeric(TOTRES[,"FEATUREpos_logIC50_MEAN"])< MC2) |
                                                                               (as.numeric(TOTRES[,"FEATUREneg_logIC50_MEAN"])< MC1 |
                                                                                  (as.numeric(TOTRES[,"FEATUREneg_logIC50_MEAN"])< MC2))))
    
    significantMeaningfullHits[i]<-length(idxs)
    
    sidxs<-which(as.numeric(TOTRES[idxs,"FEATUREpos_Glass_delta"])>=1 | 
                   as.numeric(TOTRES[idxs,"FEATUREneg_Glass_delta"])>=1)
    
    significantMeaningfullStrongHits[i]<-length(sidxs)
    
    sidxs<-idxs[which(as.numeric(TOTRES[idxs,"FEATUREpos_Glass_delta"])>=1 & 
                        as.numeric(TOTRES[idxs,"FEATUREneg_Glass_delta"])>=1)]
    
    
    if (i==5){
      
      MIPS<-sidxs
    }
    significantMeaningfullSuperStrongHits[i]<-length(sidxs)
    
    
  }
  
  toPlot<-t(cbind(significantHits,
                  significantMeaningfullHits,
                  significantMeaningfullStrongHits,
                  significantMeaningfullSuperStrongHits))
  
  bplt<-barplot(toPlot,ylim=c(0,YLIM),col=c('aquamarine4','cyan2','cornflowerblue','aquamarine'),
                xlab='FDR threshold',ylab='n.associations',las=1,beside = TRUE,border=FALSE,main=MAIN)
  
  text(y= toPlot+75, x= bplt+0.25, labels=as.character(toPlot), xpd=TRUE,cex = 0.7,las=2, srt=90)
  
  legend('topleft',legend = c('1) significant','2) 1 + meaningful','3) 2 + strong','4) 2 + very strong'),
         fill=c('aquamarine4','cyan2','cornflowerblue','aquamarine'),border=FALSE)
  
  return(MIPS)
  
}

gdscANOVA_familyBasedDrugSummary<-function(TOTRES){
  
  d_families<-colnames(DRUG_ANALYSIS_SET)[c(8:11,19,7,18)]
  
  #d_families<-colnames(DRUG_ANALYSIS_SET)[7]
  nd_families<-length(d_families)
  
  for (j in 1:4){
    events<-miniPathwayEvents[[j]]
    
    for (i in 1:nd_families){
      d_family<-d_families[i]
      dids<-annotation.drugs.retrieve.family(d_family)
      tmp<-gdscANOVA_summarySubMatrix(TOTRES,events,dids)
      if (i == 1){
        currentMiniPathMat<-tmp$subM
        currentMiniPathFDR<-tmp$FDRs
      }else{
        currentMiniPathMat<-cbind(currentMiniPathMat,tmp$subM)
        currentMiniPathFDR<-cbind(currentMiniPathFDR,tmp$FDRs)
      }
    }
    
    if (j == 1){
      TOTALPMAT<-currentMiniPathMat
      TOTALPFDR<-currentMiniPathFDR
    }else{
      TOTALPMAT<-rbind(TOTALPMAT,currentMiniPathMat)
      TOTALPFDR<-rbind(TOTALPFDR,currentMiniPathFDR)
    }
  }
  
  TOTALPMAT[TOTALPMAT>4]<-4
  TOTALPMAT[TOTALPMAT< -4]<- -4
  
  DF<-DrugFamilyVector[colnames(TOTALPMAT)]
  
  ANNOTATIONS<- data.frame(Var1 = factor(DF))
  ANNOTATIONS$Var1<-factor(ANNOTATIONS$Var1,levels=unique(DF))
  
  print(DRUG_ANALYSIS_SET[colnames(TOTALPMAT),"PUTATIVE_TARGET"])
  rownames(ANNOTATIONS)<-colnames(TOTALPMAT)
  
  TOTALPMAT<-TOTALPMAT[rownames(TOTALPMAT)!='loss_cnaPANCAN194_(FGFR3)',]
  TOTALPMAT<-TOTALPMAT[rownames(TOTALPMAT)!='gain_cnaPANCAN395_(AKT1,HSP90AA1,PPP2R5C)',]
  TOTALPMAT<-TOTALPMAT[rownames(TOTALPMAT)!='loss_cnaPANCAN2_(STK11)',]
  TOTALPMAT<-TOTALPMAT[rownames(TOTALPMAT)!='CDKN2A_mut',]
  TOTALPMAT<-TOTALPMAT[rownames(TOTALPMAT)!='gain_cnaPANCAN1_(CCNE1)',]
  TOTALPMAT<-TOTALPMAT[rownames(TOTALPMAT)!='loss_cnaPANCAN415_(B2M,BUB1B,MGA,TP53BP1)',]
  
  TOTALPFDR<-TOTALPFDR[rownames(TOTALPFDR)!='loss_cnaPANCAN194_(FGFR3)',]
  TOTALPFDR<-TOTALPFDR[rownames(TOTALPFDR)!='gain_cnaPANCAN395_(AKT1,HSP90AA1,PPP2R5C)',]
  TOTALPFDR<-TOTALPFDR[rownames(TOTALPFDR)!='loss_cnaPANCAN2_(STK11)',]
  TOTALPFDR<-TOTALPFDR[rownames(TOTALPFDR)!='CDKN2A_mut',]
  TOTALPFDR<-TOTALPFDR[rownames(TOTALPFDR)!='gain_cnaPANCAN1_(CCNE1)',]
  TOTALPFDR<-TOTALPFDR[rownames(TOTALPFDR)!='loss_cnaPANCAN415_(B2M,BUB1B,MGA,TP53BP1)',]
  
  rownames(TOTALPMAT)<-NULL
  rownames(TOTALPFDR)<-NULL
  pheatmap(TOTALPMAT,annotation = ANNOTATIONS,file='Fi_gdsc1000_GRAPHICS/widgets/tmp.pdf',width = 15,border_color = NA,
           cluster_rows = FALSE,cluster_cols = FALSE,col=colorRampPalette(c('purple','white','orange'))(100))  
  
  TOTALPFDR[which(TOTALPFDR>=57.80)]<-NA
  TOTALPFDR[which(TOTALPFDR>=20)]<- -1
  TOTALPFDR[which(TOTALPFDR>=10 & TOTALPFDR<20)]<- -2
  TOTALPFDR[which(TOTALPFDR>=0 & TOTALPFDR<10)]<- -3
  
  
  print(unique(c(TOTALPFDR)))
  pheatmap(abs(TOTALPFDR),annotation = ANNOTATIONS,legend_breaks = c(0,1,2,3,4),legend_labels = c(0,1,2,3,4),
           file='Fi_gdsc1000_GRAPHICS/widgets/tmp2.pdf',width = 15,border_color = NA,
           cluster_rows = FALSE,cluster_cols = FALSE,col=c('bisque3','darkgray','black'))  
}
gdscANOVA_summarySubMatrix<-function(TOTRES,events,dids){
  
  events<-str_replace_all(events,':','_')
  events<-str_replace_all(events,' ','_')
  
  redTOTRES<-TOTRES[which(is.element(TOTRES[,"Drug id"],dids) & is.element(TOTRES[,"FEATURE"],events)),]
  
  pvalues<- -log10(as.numeric(redTOTRES[,"FEATURE_ANOVA_pval"]))
  deltas<- as.numeric(redTOTRES[,"FEATURE_deltaMEAN_IC50"])
  
  drscores<-deltas*pvalues
  
  subM<-matrix(0,nrow = length(events),ncol=length(dids),dimnames = list(events,dids))
  FDRs<-matrix(NA,nrow = length(events),ncol=length(dids),dimnames = list(events,dids))
  
  for (i in 1:nrow(redTOTRES)){
    subM[redTOTRES[i,"FEATURE"],redTOTRES[i,"Drug id"]]<-drscores[i] 
    FDRs[redTOTRES[i,"FEATURE"],redTOTRES[i,"Drug id"]]<-as.numeric(redTOTRES[i,"ANOVA FEATURE FDR %"]) 
    
  }
  
  
  
  return(list(subM=subM,FDRs=FDRs))
}