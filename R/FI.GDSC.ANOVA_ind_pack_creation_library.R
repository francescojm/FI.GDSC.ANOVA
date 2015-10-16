
load(GDSCANOVA_SETTINGS$gdscANOVA.drugProperties.fn)
load(GDSCANOVA_SETTINGS$gdscANOVA.screening.fn)
load(GDSCANOVA_SETTINGS$gdscANOVA.maxTestedConc.fn)
load(GDSCANOVA_SETTINGS$gdscANOVA.drugOwnership.fn)



gdscANOVA_RP_copy_html_elements<-function(PATH){
  
  fl<-dir('HTML_templates_and_elements/images/')
  
  for (i in 1:length(fl)){
    file.copy(paste('HTML_templates_and_elements/images/',fl[i],sep=''),
              paste(PATH,fl[i],sep=''))  
  }
  
  
}

gdscANOVA_RP_create_individual_associations_html<-function(superSet){
    
  drug_ids_panp<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                                 DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
  
  
  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  idxs<-which(is.element(did_in_totres,drug_ids_panp))
  
  propTOTRES<-TOTRES[idxs,]
  
  range<-which(as.numeric(propTOTRES[,"ANOVA FEATURE FDR %"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)
  
  assoc_id<-propTOTRES[range,"assoc_id"]
  
  
  
  n_assoc<-length(assoc_id)
  fromD<-paste(unparsed_results_dir,superSet,'/OUTPUT/GRAPHICS/association_plots/',sep='')    
  fromD<-paste(fromD,assoc_id,'_00_scatterSet.png',sep='')
  file.copy(from=fromD,to=packages_DD_DATA_HTMLEL_ASSOC_dir)  
  fromD<-paste(unparsed_results_dir,superSet,'/OUTPUT/GRAPHICS/association_plots/',sep='')
  fromD<-paste(fromD,assoc_id,'_01_TissueWisk.png',sep='')
  file.copy(from=fromD,to=packages_DD_DATA_HTMLEL_ASSOC_dir)
  fromD<-paste(unparsed_results_dir,superSet,'/OUTPUT/GRAPHICS/association_plots/',sep='')
  fromD<-paste(fromD,assoc_id,'_02_MSIWisk.png',sep='')
  file.copy(from=fromD,to=packages_DD_DATA_HTMLEL_ASSOC_dir)
  
  if (n_assoc>1){
    pb <- txtProgressBar(min=1,max=n_assoc,style=3)
  }
  
  if (n_assoc>=1){  
    for (i in 1:n_assoc){
        if (n_assoc>1){setTxtProgressBar(pb,i)}
        gdscANOVA_RP_create_single_assoc_html(assoc_id[i])
    }
  }
  
  if (n_assoc>1){
    Sys.sleep(1)
    close(pb)
  }
}

gdscANOVA_RP_create_single_assoc_html<-function(id){

  logo1<-paste('<img src=',paste('../IMAGES/sanger-logo.png',sep=''),' title=sanger-logo align="center"/>')
  logo2<-paste('<img src=',paste('../IMAGES/logo-nki.png',sep=''),' title=sanger-logo align="center"/>')
  logo3<-paste('<img src=',paste('../IMAGES/EBI_logo.png',sep=''),' title=sanger-logo align="center"/>')
  
  logos<-''
  
  logos<-paste(logos,logo1,'\n') 
  logos<-paste(logos,logo3,'\n') 
  
  
  superheader <- paste('<br><br><center><font size=+2 face="Arial">','Result package from the GDSC1000 project','</font><br></center>\n')
  
  author<-paste('<br><center><font font size=-1 face="Arial" color="gray"><i>Francesco Iorio, EMBL - European Bioinformatics Institute','</i></font><br></center>\n')
  
  timestamp<-paste('<center><font size=-2 face="Arial" color="gray"> result parsing started on ',Sys.time(),'</font></center>',sep='')
  
  contact<-paste('<center><font font size=-1 face="Arial"><a href="mailto:iorio@ebi.ac.uk?Subject=from%20GDSC%20Data-Package%20user" target="_top">iorio@ebi.ac.uk</a></font></center>')
  
  
  if (GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType=='CS'){
    ad<-paste(  GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES,'specific')
  }else{
    ad<-'panCancer'
  }
  
  owner<-paste('<br><br><center><font size=+1 face="Arial">Drug Domain: <b>',DrugDomain,'</b>, Analysis Domain: <b>',ad,'</b></font><br></center><br><br>',sep='')
  
  body<-paste(logos,superheader,author,timestamp,contact,owner)
  
  testo <- paste('<font size=+1 face="Arial">Individual association analysis:',id,'</font><br>\n')
  
  body<-paste(body,testo)
  
  resLINE<-TOTRES[which(TOTRES[,"assoc_id"]==id),]
  resLINE<-matrix(resLINE,nrow = 1,ncol=length(resLINE),dimnames = list(NULL,colnames(TOTRES)))
  
  AT<-gdscANOVA_RP_all_assoc_summary_body(resLINE)
    
  NSIG<-AT$n
  AT<-AT$body
  
  body<-paste(body,AT)
  
  body<-paste(body,'<br>')

  scatters<-paste('<center><img src=',paste('./',id,'_00_scatterSet.png',sep=''),' width=1200 height=800 title=scatter-set align="center"/></center>')
  
  body<-paste(body,scatters)
  
  body<-paste(body,'<br>')
  if(GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType!='CS'){
    scatters<-paste('<img src=',paste('./',id,'_01_TissueWisk.png',sep=''),' width=500 height=707 title=tissue-wisk align="center"/>')
    body<-paste(body,scatters)
    
    body<-paste(body,'<br>')
  }
  
  if(GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType!='CS' & GDSCANOVA_SETTINGS$gdscANOVA.settings.includeMSI_Factor){
    scatters<-paste('<img src=',paste('./',id,'_02_MSIWisk.png',sep=''),' width=505 height=707 title=msi-wisk align="center"/>')
    body<-paste(body,scatters)
  }
  
  write(body, file = paste(packages_DD_DATA_HTMLEL_ASSOC_dir,id,".html",sep=''),append=FALSE)
   
  
}
gdscANOVA_RP_all_assoc_summary_body<-function(redTOTRES){
  
  if (is.vector(redTOTRES)){
    redTOTRES<-matrix(redTOTRES,nrow = 1,ncol=length(redTOTRES),dimnames = list(NULL,names(redTOTRES)))
  }
  
  if(  GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES=='PANCAN'){
    range<-which(as.numeric(redTOTRES[,"ANOVA FEATURE FDR %"])<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)  
  }else{
    range<-which(as.numeric(redTOTRES[,"ANOVA FEATURE FDR %"])<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH &
                   as.numeric(redTOTRES[,"FEATURE_ANOVA_pval"])< GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH)
  }
  
  
  
  if (length(range)>0){
    DELTAS<-as.numeric(redTOTRES[,"FEATURE_deltaMEAN_IC50"])
    limits<- c(-max(abs(DELTAS)),max(abs(DELTAS)))
    DELTAMcolors<-gdscANOVA_RP_map2color(DELTAS,colorRampPalette(c('green','white','red'))(200),limits=limits)
  }
  
  header <- paste('<center><table border=1 width=100% cellspacing=1 cellpadding=2 cols=7 style="font-family: Arial; font-size: 10px">\n
                  <tr><td><b><center>association id</center></b></td>
                  <td><b><center>FEATURE</center></b></td>
                  <td><b><center>Drug id</center></b></td>
                  <td><b><center>Propr</center></b></td>
                  <td><b><center>Public</center></b></td>
                  <td><b><center>Owned by</center></b></td>
                  <td><b><center>Drug Name</center></b></td>
                  <td><b><center>Drug Target</center></b></td>
                  <td><b><center>N FEATURE pos</center></b></td>
                  <td><b><center>N FEATURE neg</center></b></td>
                  <td><b><center>log max Conc tested</center></b></td>
                  <td><b><center>FEATURE neg logIC50 MEAN</center></b></td>
                  <td><b><center>FEATURE pos logIC50 MEAN</center></b></td>
                  <td><b><center>FEATURE deltaMEAN IC50</center></b></td>
                  <td><b><center>FEATURE IC50 effectSize</center></b></td>
                  <td><b><center>FEATURE neg Glass delta</center></b></td>
                  <td><b><center>FEATURE pos Glass delta</center></b></td>
                  <td><b><center>ANOVA FEATURE pval</center></b></td>
                  <td><b><center>Tissue ANOVA pval</center></b></td>
                  <td><b><center>MSI ANOVA pval</center></b></td>
                  <td><b><center>ANOVA FEATURE FDR %</center></b></td></tr>\n')
  
  tail<-paste('</tr>\n</table></center>')
  
  DID<-unlist(str_split(redTOTRES[range,'Drug id'],'_'))[seq(1,nrow(redTOTRES)*2,2)]
  
  posit<-match(DID,DRUG_BY_COMPANIES[,'DRUG_ID'])
  
  owners<-DRUG_PROPS[as.character(DID),'OWNED_BY']
  
  publicDRUG<-DRUG_BY_COMPANIES[posit,"Web Released"]
  proprDRUG<-DRUG_BY_COMPANIES[posit,DrugDomain]
  
  publicDRUGCOLOR<-rep('white',length(publicDRUG))
  publicDRUGCOLOR[which(publicDRUG==1)]<-'blue'
  
  proprDRUGCOLOR<-rep('white',length(publicDRUG))
  proprDRUGCOLOR[which(proprDRUG==1)]<-'blue'
  
  
  body<-'\n'
  
  if (length(range)>0){
    
    for (i in 1:length(range)){
      
      currentLine<-paste('<td><center><a href="../associations/',redTOTRES[i,"assoc_id"],'.html" target="_blank">',redTOTRES[i,1],'</a></center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center><a href="../FEATURES/',as.character(redTOTRES[i,"FEATURE"]),'.html" target="_blank">',as.character(redTOTRES[i,2]),'</a></center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center><a href="../DRUGS/',as.character(redTOTRES[i,"Drug id"]),'.html" target="_blank">',as.character(redTOTRES[i,3]),'</a></center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      currentLine<-paste('<td bgcolor=',proprDRUGCOLOR[i],'><center>',as.character(proprDRUG[i]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td bgcolor=',publicDRUGCOLOR[i],'><center>',as.character(publicDRUG[i]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      currentLine<-paste('<td><center>',owners[i],'</center></td>\n',sep='')
      body<-paste(body,currentLine)
                  
      currentLine<-paste('<td><center>',as.character(redTOTRES[i,"Drug name"]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center>',as.character(redTOTRES[i,"Drug Target"]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center>',as.character(redTOTRES[i,"N_FEATURE_pos"]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center>',as.character(redTOTRES[i,"N_FEATURE_neg"]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      maxConc<-as.numeric(redTOTRES[i,c("log max.Conc.tested","log max.Conc.tested2")])
      maxConc<-maxConc[!is.na(maxConc)]
      maxConc<-format(maxConc,digits=3)
      maxConc<-paste(maxConc,collapse=', ')
      currentLine<-paste('<td><center>',maxConc,'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      for (k in c("FEATUREneg_logIC50_MEAN","FEATUREpos_logIC50_MEAN")){
        currentLine<-paste('<td><center>',format(as.numeric(redTOTRES[i,k]),digits=3),'</center></td>\n',sep='')
        body<-paste(body,currentLine)
      }
      
      
      currentLine<-paste('<td bgcolor=',DELTAMcolors[i],'><center>',format(as.numeric(redTOTRES[i,'FEATURE_deltaMEAN_IC50']),digits=3),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      effectS<-as.numeric(redTOTRES[i,"FEATURE_IC50_effect_size"])
      
      if(effectS>=2){effectS<-paste(format(effectS,digits=2),'**')
                     COLOR='#0099FF'}
      else{if(effectS>=1){effectS<-paste(format(effectS,digits=2),'*')
                          COLOR='#00CCFF'}
           else{if(effectS>=0.5){effectS<-paste(format(effectS,digits=2),'.')
                                 COLOR='#99FFFF'}
                else{
                  effectS<-format(effectS,digits=2)
                  COLOR='white'
                }
           }
      }
      
      currentLine<-paste('<td bgcolor=',COLOR,'><center>',effectS,'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      
      effectS<-as.numeric(redTOTRES[i,"FEATUREneg_Glass_delta"])
      
      if(effectS>=2){effectS<-paste(format(effectS,digits=2),'**')
                     COLOR='#0099FF'}
      else{if(effectS>=1){effectS<-paste(format(effectS,digits=2),'*')
                          COLOR='#00CCFF'}
           else{if(effectS>=0.5){effectS<-paste(format(effectS,digits=2),'.')
                                 COLOR='#99FFFF'}
                else{
                  effectS<-format(effectS,digits=2)
                  COLOR='white'
                }
           }
      }
      
      currentLine<-paste('<td bgcolor=',COLOR,'><center>',effectS,'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      
      
      effectS<-as.numeric(redTOTRES[i,"FEATUREpos_Glass_delta"])
      
      if(effectS>=2){effectS<-paste(format(effectS,digits=2),'**')
                     COLOR='#0099FF'}
      else{if(effectS>=1){effectS<-paste(format(effectS,digits=2),'*')
                          COLOR='#00CCFF'}
           else{if(effectS>=0.5){effectS<-paste(format(effectS,digits=2),'.')
                                 COLOR='#99FFFF'}
                else{
                  effectS<-format(effectS,digits=2)
                  COLOR='white'
                }
           }
      }
      
      currentLine<-paste('<td bgcolor=',COLOR,'><center>',effectS,'</center></td>\n',sep='')
      body<-paste(body,currentLine)
    
      
      MP<-format(as.numeric(redTOTRES[i,"FEATURE_ANOVA_pval"]),digits=3,scientific=TRUE)          
      currentLine<-paste('<td><center>',MP,'</center></td>\n',sep='')           
      body<-paste(body,currentLine)
      MP<-format(as.numeric(redTOTRES[i,"Tissue_ANOVA_pval"]),digits=3,scientific=TRUE)          
      currentLine<-paste('<td><center>',MP,'</center></td>\n',sep='')           
      body<-paste(body,currentLine)
      MP<-format(as.numeric(redTOTRES[i,"MSI_ANOVA_pval"]),digits=3,scientific=TRUE)          
      currentLine<-paste('<td><center>',MP,'</center></td>\n',sep='')           
      body<-paste(body,currentLine)
      MP<-format(as.numeric(redTOTRES[i,"ANOVA FEATURE FDR %"]),digits=2)          
      currentLine<-paste('<td><center>',MP,'</center></td>\n',sep='')           
      body<-paste(body,currentLine)
      
      
      body<-paste(body,'</tr>\n')
    }
    
    
    body<-paste(body,'</table>\n')
    body<-paste(header,body,tail,sep='')
  
  }
  
  return(list(body=body,n=length(range)))
}
gdscANOVA_RP_all_assoc_summary<-function(superSet,redTOTRES,PATH){
  
  element_folder<-paste('./DATA/HTML_elements/')
  

    
  logo1<-paste('<img src=','./IMAGES/sanger-logo.png',' title=sanger-logo align="center"/>',sep='')
  logo2<-paste('<img src=','./IMAGES/logo-nki.png',' title=sanger-logo align="center"/>',sep='')
  logo3<-paste('<img src=','./IMAGES/EBI_logo.png',' title=sanger-logo align="center"/>',sep='')
  
  logos<-paste(logo1,logo3,'\n')
  
  superheader <- paste('<br><br><center><font size=+2 face="Arial">','Result package from the GDSC1000 project','</font><br></center>\n')
  
  author<-paste('<br><center><font font size=-1 face="Arial" color="gray"><i>Francesco Iorio, EMBL - European Bioinformatics Institute','</i></font><br></center>\n')
  
  timestamp<-paste('<center><font size=-2 face="Arial" color="gray"> result parsing started on ',Sys.time(),'</font></center>',sep='')
  
  contact<-paste('<center><font font size=-1 face="Arial"><a href="mailto:iorio@ebi.ac.uk?Subject=from%20GDSC%20Data-Package%20user" target="_top">iorio@ebi.ac.uk</a></font></center>')
  
  if (GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType=='CS'){
    ad<-paste(  GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES,'specific')
  }else{
    ad<-'panCancer'
  }
  
  owner<-paste('<br><br><center><font size=+1 face="Arial">Drug Domain: <b>',DrugDomain,'</b>, Analysis Domain: <b>',ad,'</b></font><br></center><br><br>',sep='')
  
  
  print('creating all associations summary')
  range<-which(as.numeric(redTOTRES[,"ANOVA FEATURE FDR %"])<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH &
                 as.numeric(redTOTRES[,"FEATURE_ANOVA_pval"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH)
  
  
  nn<-length(range)
  
  if (nn<=1){
    nn<-2
  }
  
  pb <- txtProgressBar(min=1,max=nn,style=3)
  if (length(range)>0){
    DELTAS<-as.numeric(redTOTRES[,"FEATURE_deltaMEAN_IC50"])
    limits<- c(-max(abs(DELTAS)),max(abs(DELTAS)))
    DELTAMcolors<-gdscANOVA_RP_map2color(DELTAS,colorRampPalette(c('green','white','red'))(200),limits=limits)
    }
  
  header <- paste('<center><font size=+2 face="Arial">MANOVA Results Summary</font><br><center>\n
                  <center><table border=1 width=100% cellspacing=1 cellpadding=2 cols=7 style="font-family: Arial; font-size: 10px">\n
                  <tr><td><b><center>association id</center></b></td>
                  <td><b><center>FEATURE</center></b></td>
                  <td><b><center>Drug id</center></b></td>
                  <td><b><center>Propr</center></b></td>
                  <td><b><center>Public</center></b></td>
                  <td><b><center>Owned by</center></b></td>
                  <td><b><center>Drug Name</center></b></td>
                  <td><b><center>Drug Target</center></b></td>
                  <td><b><center>N FEATURE pos</center></b></td>
                  <td><b><center>N FEATURE neg</center></b></td>
                  <td><b><center>log max Conc tested</center></b></td>
                  <td><b><center>FEATURE pos logIC50 MEAN</center></b></td>
                  <td><b><center>FEATURE neg logIC50 MEAN</center></b></td>
                  <td><b><center>FEATURE deltaMEAN IC50</center></b></td>
                  <td><b><center>FEATURE IC50 effectSize</center></b></td>
                  <td><b><center>FEATURE neg Glass delta</center></b></td>
                  <td><b><center>FEATURE pos Glass delta</center></b></td>
                  <td><b><center>ANOVA FEATURE pval</center></b></td>
                  <td><b><center>ANOVA Tissue pval</center></b></td>
                  <td><b><center>ANOVA MSI pval</center></b></td>
                  <td><b><center>ANOVA FEATURE %FDR</center></b></td></tr>\n')
  
  tail<-paste('</tr>\n</table></center>')
  
  body<-'\n'
  
  if (length(range)>0){
    
  
  DID<-unlist(str_split(redTOTRES[range,'Drug id'],'_'))[seq(1,length(range)*2,2)]
  
  posit<-match(DID,DRUG_BY_COMPANIES[,'DRUG_ID'])
  
  publicDRUG<-DRUG_BY_COMPANIES[posit,'Web Released']
  proprDRUG<-DRUG_BY_COMPANIES[posit,DRUG_DOMAIN]
  ownerDRUG<-DRUG_PROPS[posit,"OWNED_BY"]
  
  ownerDRUG[is.na(ownerDRUG)]<''
  
  publicDRUGCOLOR<-rep('white',length(publicDRUG))
  publicDRUGCOLOR[which(publicDRUG==1)]<-'blue'
  
  proprDRUGCOLOR<-rep('white',length(publicDRUG))
  proprDRUGCOLOR[which(proprDRUG==1)]<-'blue'
    
  body<-'\n'
  

    for (i in 1:length(range)){
      setTxtProgressBar(pb,i)
      
      currentLine<-paste('<td><center><a href="./associations/',redTOTRES[i,"assoc_id"],'.html" target="_blank">',redTOTRES[i,1],'</a></center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center><a href="./FEATURES/',as.character(redTOTRES[i,"FEATURE"]),'.html" target="_blank">',as.character(redTOTRES[i,2]),'</a></center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center><a href="./DRUGS/',as.character(redTOTRES[i,"Drug id"]),'.html" target="_blank">',as.character(redTOTRES[i,3]),'</a></center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      currentLine<-paste('<td bgcolor=',proprDRUGCOLOR[i],'><center>',as.character(proprDRUG[i]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td bgcolor=',publicDRUGCOLOR[i],'><center>',as.character(publicDRUG[i]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      currentLine<-paste('<td><center>',ownerDRUG[i],'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      currentLine<-paste('<td><center>',as.character(redTOTRES[i,"Drug name"]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center>',as.character(redTOTRES[i,"Drug Target"]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center>',as.character(redTOTRES[i,"N_FEATURE_pos"]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      currentLine<-paste('<td><center>',as.character(redTOTRES[i,"N_FEATURE_neg"]),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      maxConc<-as.numeric(redTOTRES[i,c("log max.Conc.tested","log max.Conc.tested2")])
      maxConc<-maxConc[!is.na(maxConc)]
      maxConc<-format(maxConc,digits=3)
      maxConc<-paste(maxConc,collapse=', ')
      currentLine<-paste('<td><center>',maxConc,'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      for (k in c("FEATUREpos_logIC50_MEAN","FEATUREneg_logIC50_MEAN")){
        currentLine<-paste('<td><center>',format(as.numeric(redTOTRES[i,k]),digits=3),'</center></td>\n',sep='')
        body<-paste(body,currentLine)
      }
      
      
      currentLine<-paste('<td bgcolor=',DELTAMcolors[i],'><center>',format(as.numeric(redTOTRES[i,'FEATURE_deltaMEAN_IC50']),digits=3),'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      effectS<-as.numeric(redTOTRES[i,"FEATURE_IC50_effect_size"])
      
      if(effectS>=2){effectS<-paste(format(effectS,digits=2),'**')
                     COLOR='#0099FF'}
      else{if(effectS>=1){effectS<-paste(format(effectS,digits=2),'*')
                          COLOR='#00CCFF'}
           else{if(effectS>=0.5){effectS<-paste(format(effectS,digits=2),'.')
                                 COLOR='#99FFFF'}
                else{
                  effectS<-format(effectS,digits=2)
                  COLOR='white'
                }
           }
      }
      
      
      
      currentLine<-paste('<td bgcolor=',COLOR,'><center>',effectS,'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      effectS<-as.numeric(redTOTRES[i,"FEATUREneg_Glass_delta"])
      
      if(effectS>=2){effectS<-paste(format(effectS,digits=2),'**')
                     COLOR='#0099FF'}
      else{if(effectS>=1){effectS<-paste(format(effectS,digits=2),'*')
                          COLOR='#00CCFF'}
           else{if(effectS>=0.5){effectS<-paste(format(effectS,digits=2),'.')
                                 COLOR='#99FFFF'}
                else{
                  effectS<-format(effectS,digits=2)
                  COLOR='white'
                }
           }
      }
      
      currentLine<-paste('<td bgcolor=',COLOR,'><center>',effectS,'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      
      
      effectS<-as.numeric(redTOTRES[i,"FEATUREpos_Glass_delta"])
      
      if(effectS>=2){effectS<-paste(format(effectS,digits=2),'**')
                     COLOR='#0099FF'}
      else{if(effectS>=1){effectS<-paste(format(effectS,digits=2),'*')
                          COLOR='#00CCFF'}
           else{if(effectS>=0.5){effectS<-paste(format(effectS,digits=2),'.')
                                 COLOR='#99FFFF'}
                else{
                  effectS<-format(effectS,digits=2)
                  COLOR='white'
                }
           }
      }
      
      currentLine<-paste('<td bgcolor=',COLOR,'><center>',effectS,'</center></td>\n',sep='')
      body<-paste(body,currentLine)
      
      
      MP<-format(as.numeric(redTOTRES[i,"FEATURE_ANOVA_pval"]),digits=3,scientific=TRUE)          
      currentLine<-paste('<td><center>',MP,'</center></td>\n',sep='')           
      body<-paste(body,currentLine)
      MP<-format(as.numeric(redTOTRES[i,"Tissue_ANOVA_pval"]),digits=3,scientific=TRUE)          
      currentLine<-paste('<td><center>',MP,'</center></td>\n',sep='')           
      body<-paste(body,currentLine)
      MP<-format(as.numeric(redTOTRES[i,"MSI_ANOVA_pval"]),digits=3,scientific=TRUE)          
      currentLine<-paste('<td><center>',MP,'</center></td>\n',sep='')           
      body<-paste(body,currentLine)
      MP<-format(as.numeric(redTOTRES[i,"ANOVA FEATURE FDR %"]),digits=2)          
      currentLine<-paste('<td><center>',MP,'</center></td>\n',sep='')           
      body<-paste(body,currentLine)
      
      
      body<-paste(body,'</tr>\n')
    }
    
    
    body<-paste(body,'</table>\n')
    body<-paste(logos,superheader,author,timestamp,contact,owner,header,body,tail,sep='')
    
    Sys.sleep(1)
    close(pb)
  }
  

  write(body, file = paste(PATH,"000_Significant_Hits.html",sep=''),append=FALSE)
}


gdscANOVA_RP_map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
gdscANOVA_RP_create_individual_drug_html<-function(superSet){
  
  
  drug_ids_panp<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                        DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
  
  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  idxs<-which(is.element(did_in_totres,drug_ids_panp))
  
  DIDS<-unique(TOTRES[idxs,"Drug id"])
  
  
  n_drugs<-length(DIDS)
  
  fromD<-paste(unparsed_results_dir,superSet,'/OUTPUT/GRAPHICS/DRUG_volcanos/',sep='')    
  fromD<-paste(fromD,str_replace(DIDS,'/','_'),'.png',sep='')
  file.copy(from=fromD,to=packages_DD_DATA_HTMLEL_DRUG_dir)  
  
  if (length(idxs)>0){
    
    pb <- txtProgressBar(min=1,max=n_drugs,style=3)
  
    for (i in 1:n_drugs){
      setTxtProgressBar(pb,i)
      gdscANOVA_RP_create_single_drug_html(DIDS[i])
    }
  
    Sys.sleep(1)
    close(pb)
  }
}
gdscANOVA_RP_create_single_drug_html<-function(id){
  
  logo1<-paste('<img src=',paste('../IMAGES/sanger-logo.png',sep=''),' title=sanger-logo align="center"/>')
  logo2<-paste('<img src=',paste('../IMAGES/logo-nki.png',sep=''),' title=sanger-logo align="center"/>')
  logo3<-paste('<img src=',paste('../IMAGES/EBI_logo.png',sep=''),' title=sanger-logo align="center"/>')
  
  logos<-''
  
  logos<-paste(logos,logo1,'\n') 
  logos<-paste(logos,logo3,'\n') 
  
  superheader <- paste('<br><br><center><font size=+2 face="Arial">','Result package from the GDSC1000 project','</font><br></center>\n')
  
  author<-paste('<br><center><font font size=-1 face="Arial" color="gray"><i>Francesco Iorio, EMBL - European Bioinformatics Institute','</i></font><br></center>\n')
  
  timestamp<-paste('<center><font size=-2 face="Arial" color="gray"> result parsing started on ',Sys.time(),'</font></center>',sep='')
  
  contact<-paste('<center><font font size=-1 face="Arial"><a href="mailto:iorio@ebi.ac.uk?Subject=from%20GDSC%20Data-Package%20user" target="_top">iorio@ebi.ac.uk</a></font></center>')
  
  if (GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType=='CS'){
    ad<-paste(  GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES,'specific')
  }else{
    ad<-'panCancer'
  }
  
  owner<-paste('<br><br><center><font size=+1 face="Arial">Drug Domain: <b>',DrugDomain,'</b>, Analysis Domain: <b>',ad,'</b></font><br></center><br><br>',sep='')
    
  body<-paste(logos,superheader,author,timestamp,contact,owner)
  
  testo <- paste('<font size=+1 face="Arial">Individual drug analysis:',id,'</font><br>\n')
  
  body<-paste(body,testo)
  
  drugName<-DRUG_PROPS[str_split(id,'_')[[1]][1],"DRUG_NAME"]
  drugSynonyms<-DRUG_PROPS[str_split(id,'_')[[1]][1],"SYNONYMS"]
  drugBrandName<-DRUG_PROPS[str_split(id,'_')[[1]][1],"OWNED_BY"]
  drugTarget<-DRUG_PROPS[str_split(id,'_')[[1]][1],"PUTATIVE_TARGET"]
  
  concRange<-unique(as.numeric(maxConcTested[,id]))
  
  concRange<-concRange[!is.na(concRange)]
  
  minC<-concRange/(4^4)
  
  concRange<-paste(minC,'to',concRange,'uM')
  
  if (length(concRange)>1){
    concRange<-paste(concRange,collapse=' and ')
  }
    
  ncellScreened<-length(which(!is.na(IC50s[,id])))
  
  drug_details <- paste('<br><br><font size=+1 face="Arial"><b>Drug details and screening numbers</b></font><br>\n')
  
  body<-paste(body,drug_details)
  
  bunch<-paste('<br><font face="Arial">Drug Name:',drugName,'</font><br>\n')
  body<-paste(body,bunch)
  bunch<-paste('<font face="Arial">Drug ID:',id,'</font><br>\n')
  body<-paste(body,bunch)
  bunch<-paste('<font face="Arial">Synonyms:',drugSynonyms,'</font><br>\n')
  body<-paste(body,bunch)
  bunch<-paste('<font face="Arial">Owned By:',drugBrandName,'</font><br>\n')
  body<-paste(body,bunch)
  bunch<-paste('<font face="Arial">Target:',drugTarget,'</font><br>\n')
  body<-paste(body,bunch)
  bunch<-paste('<br><font face="Arial">Number of cell lines screened:',ncellScreened,'</font><br>\n')
  body<-paste(body,bunch)
  
  bunch<-paste('<font face="Arial">Screening concentration range(s):',concRange,'</font><br>\n')
  body<-paste(body,bunch)
  
  idxs<-which(TOTRES[,'Drug id']==id)
  
  propTOTRES<-TOTRES[idxs,]
  
  AT<-gdscANOVA_RP_all_assoc_summary_body(redTOTRES = propTOTRES)
  NSIG<-AT$n
  AT<-AT$body
  bunch <- paste('<br><br><font size=+1 face="Arial"><b>Statistically significant associations identified =',NSIG,'</b></font><br>\n')
  body <- paste(body,bunch)
  
  body<-paste(body,AT)
  
  bunch <- paste('<br><br><font size=+1 face="Arial"><b>All-tests volcano plot</b></font><br>\n')
  body <- paste(body,bunch)
  
  vp<-paste('<center><img src=',paste('./',id,'.png',sep=''),' title=single-d-volcano align="center"/></center>')
  
  body <- paste(body,vp)

  write(body, file = paste(packages_DD_DATA_HTMLEL_DRUG_dir,str_replace(id,'/','_'),".html",sep=''),append=FALSE)
}

gdscANOVA_RP_create_individual_feature_html<-function(superSet){
  
  drug_ids_panp<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                        DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
  
  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  idxs<-which(is.element(did_in_totres,drug_ids_panp))
  
  propTOTRES<-TOTRES[idxs,]
  
  feat<-unique(propTOTRES[,'FEATURE'])
  
  n_feat<-length(feat)
  
  nn<-n_feat
  
  if (length(idxs)>0){
    
  
    if (nn<=1){nn<-2}
    pb <- txtProgressBar(min=1,max=nn,style=3)
  
      for (i in 1:n_feat){
        setTxtProgressBar(pb,i)
        gdscANOVA_RP_create_single_feat_html(propTOTRES=propTOTRES,ff = feat[i],superSet = superSet)
      }
  
      Sys.sleep(1)
      close(pb)
    }
}

gdscANOVA_RP_create_single_feat_html<-function(propTOTRES,ff,superSet){
  ORF<-ff
  gdscANOVA_FEATURE_volcanoPlot(FEATURE=ff,print=TRUE,cTOTRES=propTOTRES,PATH=packages_DD_DATA_HTMLEL_FEATURE_dir)
  
  logo1<-paste('<img src=',paste('../IMAGES/sanger-logo.png',sep=''),' title=sanger-logo align="center"/>')
  logo2<-paste('<img src=',paste('../IMAGES/logo-nki.png',sep=''),' title=sanger-logo align="center"/>')
  logo3<-paste('<img src=',paste('../IMAGES/EBI_logo.png',sep=''),' title=sanger-logo align="center"/>')
  
  logos<-''
  
  logos<-paste(logos,logo1,'\n') 
  logos<-paste(logos,logo3,'\n') 
  
  superheader <- paste('<br><br><center><font size=+2 face="Arial">','Result package from the GDSC1000 project','</font><br></center>\n')
  
 
  author<-paste('<br><center><font font size=-1 face="Arial" color="gray"><i>Francesco Iorio, EMBL - European Bioinformatics Institute','</i></font><br></center>\n')
   
  timestamp<-paste('<center><font size=-2 face="Arial" color="gray"> result parsing started on ',Sys.time(),'</font></center>',sep='')
   
  contact<-paste('<center><font font size=-1 face="Arial"><a href="mailto:iorio@ebi.ac.uk?Subject=from%20GDSC%20Data-Package%20user" target="_top">iorio@ebi.ac.uk</a></font></center>')
   
  if (GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType=='CS'){
    ad<-paste(GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES,'specific')
  }else{
    ad<-'panCancer'
  }
  
  owner<-paste('<br><br><center><font size=+1 face="Arial">Drug Domain: <b>',DrugDomain,'</b>, Analysis Domain: <b>',ad,'</b></font><br></center><br><br>',sep='')
  
  body<-paste(logos,superheader,author,timestamp,contact,owner)
  
  testo <- paste('<font size=+1 face="Arial">Individual Feature analysis:',ff,'</font><br>\n')
  
  body<-paste(body,testo)
   
  if (str_detect(ff,'mut')){
     
     id<-which(propTOTRES[,'FEATURE']==ff)[1]
     featDescription<-paste('Binary feature equal to 1 for samples harboring mutations in',
                            str_sub(ff,1,str_length(ff)-4))
  } else{
     mods<-str_sub(ff,1,4)
     ttt<-str_sub(ff,6,str_length(ff))
     if (mods=='gain'){mods<-'amplified'}
     else{mods<-'deleted'}
     featDescription<-paste('Binary feature equal to 1 for samples having the genomic region(s)',ttt,mods,'(see genomic region decoding file for further infos on this region)')
  }
    
  feature_details <- paste('<br><br><font size=+1 face="Arial"><b>Feature details and cell panel numbers</b></font><br>\n')
  body<-paste(body,feature_details)
  
  
  bunch<-paste('<br><font face="Arial">',featDescription,'</font><br>\n')
  body<-paste(body,bunch)
  
  npositive_samples<-sum(InputFeatures$BEM[ORF,])
  bunch<-paste('<br><font face="Arial">Number of cell lines positive for this feature:',npositive_samples,'</font><br>\n')
  
  body<-paste(body,bunch)
  
  idxs<-which(propTOTRES[,'FEATURE']==ORF)
  
  propTOTRES<-propTOTRES[idxs,]
  
  AT<-gdscANOVA_RP_all_assoc_summary_body(propTOTRES)
  NSIG<-AT$n
  AT<-AT$body
  bunch <- paste('<br><br><font size=+1 face="Arial"><b>Statistically significant associations identified =',NSIG,'</b></font><br>\n')
  body <- paste(body,bunch)
  
  body<-paste(body,AT)
  bunch <- paste('<br><br><font size=+1 face="Arial"><b>All-tests volcano plot</b></font><br>\n')
  body <- paste(body,bunch)
  
  vp<-paste('<center><img src=',paste('./',ORF,'.png',sep=''),' title=single-f-volcano align="center"/></center>')
  
  body <- paste(body,vp)
  
  
  write(body, file = paste(packages_DD_DATA_HTMLEL_FEATURE_dir,ORF,".html",sep=''),append=FALSE)
}
gdscANOVA_RP_create_MANOVA_INPUT_OUTPUT_files<-function(superSet){
  
  load(paste(GDSCANOVA_SETTINGS$gdscANOVA.results.dir,superSet,'/INPUT/','InputFeatures.rdata',sep=''))
  load(GDSCANOVA_SETTINGS$annotations.master_list.fn)
  
  drug_ids<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                     DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
  
  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  
  idxs<-which(is.element(did_in_totres,drug_ids))
  
  proprietaryTOTRES<-TOTRES[idxs,]
  
  save(proprietaryTOTRES,file=paste(packages_DD_DATA_OUTPUT_dir,'ANOVA_results.rdata',sep=''))
  write.table(proprietaryTOTRES,quote=FALSE,row.names=F,sep='\t',file=paste(packages_DD_DATA_OUTPUT_dir,'ANOVA_results.txt',sep=''))
  
  genomicDATA<-InputFeatures$BEM
  
  drug_ids<-intersect(drug_ids,colnames(IC50s))
  IC50DATA<-IC50s[,drug_ids]
  
  MSIDATA<-InputFeatures$MSI_VARIABLE
  
  TISSUEDATA<-InputFeatures$TISSUES
  
  commonCL<-intersect(row.names(IC50DATA),colnames(genomicDATA))
  
  genomicDATA<-genomicDATA[,commonCL]
  
  IC50DATA<-IC50DATA[commonCL,]
  
  MSIDATA<-MSIDATA[commonCL]
  
  TISSUEDATA<-TISSUEDATA[commonCL]
  
  samplesNames<-MASTER_LIST[commonCL,"Analysis.Set.Name"]
  
    
    mutData<-genomicDATA[which(!str_detect(rownames(genomicDATA),pattern = 'gain') & !str_detect(rownames(genomicDATA),pattern = 'loss') &
                               !str_detect(rownames(genomicDATA),pattern = 'HypMUT')),]
    
    cnvData<-genomicDATA[!which(str_detect(rownames(genomicDATA),pattern = 'mut') & !str_detect(rownames(genomicDATA),pattern = 'HypMUT')),]
    
    HmetData<-genomicDATA[which(str_detect(rownames(genomicDATA),pattern = 'HypMUT')),]
        
    manova_INPUT<-t(rbind(mutData,cnvData,HmetData))
    
    cl<-intersect(rownames(manova_INPUT),rownames(IC50s))
    
    dids<-unlist(str_split(colnames(IC50s),'_'))[seq(1,ncol(IC50s)*2,2)]
    idxs<-which(DRUG_BY_COMPANIES[dids,'Web Released']>0 | DRUG_BY_COMPANIES[dids,DRUG_DOMAIN]>0)  
  
    manova_INPUT<-cbind(cl,MASTER_LIST[cl,"Analysis.Set.Name"],InputFeatures$TISSUES[cl],
                  InputFeatures$MSI_VARIABLE[cl],t(InputFeatures$BEM[,cl]),IC50s[cl,idxs])
        
    colnames(manova_INPUT)[1:4]<-c('COSMIC ID','Sample Name','Tissue Factor Value','MS-instability Factor Value')
    colnames(manova_INPUT)[(5+nrow(InputFeatures$BEM)):ncol(manova_INPUT)]<-paste('Drug_',colnames(manova_INPUT)[(5+nrow(InputFeatures$BEM)):ncol(manova_INPUT)],'_IC50',sep='')
    write.table(manova_INPUT,file=paste(packages_DD_DATA_INPUT_dir,'ANOVA_input.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)

}
gdscANOVA_RP_create_DRUG_DECODE_file<-function(){
  
  drug_ids<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                     DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
  
  
  DRUG_DECODE<-DRUG_PROPS[drug_ids,]
  
  write.table(DRUG_DECODE,file=paste(packages_DD_DATA_INPUT_dir,'DRUG_DECODE.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
}
gdscANOVA_RP_create_GenomicRegions_DECODE_file<-function(){
  
  load('HTML_templates_and_elements/cna_region_annotations.rdata')
  
  idxs<-which(str_detect(rownames(cna_region_annotations),GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES))
  
  GENOMIC_REGION_DECODE<-cna_region_annotations[idxs,]
  GENOMIC_REGION_DECODE<-cbind(rownames(GENOMIC_REGION_DECODE),GENOMIC_REGION_DECODE)
  colnames(GENOMIC_REGION_DECODE)[1]<-'Region Id'
  
  write.table(GENOMIC_REGION_DECODE[,c(1,2,4,5,6,7,8,13,14,15,16)],
              file=paste(packages_DD_DATA_INPUT_dir,'GENOMIC_REGIONS_DECODE.txt',sep=''),sep='\t',quote=FALSE,row.names=FALSE)
  
}
gdscANOVA_RP_create_comprehensive_vp<-function(){
  drug_ids_panp<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                        DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
  
 
  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  idxs<-which(is.element(did_in_totres,drug_ids_panp))
  
  panpTOTRES<-TOTRES[idxs,]
  
  delta<-sign(as.numeric(panpTOTRES[,"FEATURE_deltaMEAN_IC50"]))*as.numeric(panpTOTRES[,"FEATURE_IC50_effect_size"])
  pval<-as.numeric(panpTOTRES[,"FEATURE_ANOVA_pval"])
  fdr<-as.numeric(panpTOTRES[,"ANOVA FEATURE FDR %"])
  N<-as.numeric(panpTOTRES[,"N_FEATURE_pos"])
  labels<-paste(panpTOTRES[,'Drug name'],panpTOTRES[,'Drug Target'],panpTOTRES[,'FEATURE'])
  
  png(paste(packages_DD_DATA_OUTPUT_dir,'comprehensive_volcanoPlot.png',sep=''),768,1024)
  gdscANOVA_volcanoPlot_T(delta=delta,pval=pval,qvals=fdr,N=N,effth=1)
  dev.off()
}

gdscANOVA_RP_create_summaries<-function(){
  drug_ids_panp<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                        DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
  
  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  idxs<-which(is.element(did_in_totres,drug_ids_panp))
  
  panpTOTRES<-TOTRES[idxs,]
  
  gdscANOVA_computeFeatureStats(redTOTRES=panpTOTRES,fdrTH =   GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH,pvalTH = GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH,
                                save=TRUE,display=TRUE,printTOfig=TRUE,PATH=paste(packages_DD_DATA_OUTPUT_dir,
                                                                                             'FeatureSummary',sep=''))
  gdscANOVA_computeDrugStats(redTOTRES=panpTOTRES,fdrTH =   GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH,pvalTH = GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH,
                             save=TRUE,display=TRUE,printTOfig=TRUE,PATH=paste(packages_DD_DATA_OUTPUT_dir,
                                                                                          'DrugSummary',sep=''))
}

gdscANOVA_RP_start_page_creation<-function(superSet,PATH,packageName){
  
  
  element_folder<-paste('./DATA/HTML_elements/')
  
  URLsuperSet<-str_replace_all(superSet,' ','%20')
  
  logo1<-paste('<img src=',paste('./',URLsuperSet,'/DATA/HTML_elements/IMAGES/sanger-logo.png',sep=''),' title=sanger-logo align="center"/>')
  logo2<-paste('<img src=',paste('./',URLsuperSet,'/DATA/HTML_elements/IMAGES/logo-nki.png',sep=''),' title=sanger-logo align="center"/>')
  logo3<-paste('<img src=',paste('./',URLsuperSet,'/DATA/HTML_elements/IMAGES/EBI_logo.png',sep=''),' title=sanger-logo align="center"/>')
  
  logos<-''
  
  logos<-paste(logos,logo1,'\n') 
  logos<-paste(logos,logo3,'\n') 
  
  superheader <- paste('<br><br><center><font size=+2 face="Arial">','Result package from the GDSC1000 project','</font><br></center>\n')
  
  author<-paste('<br><center><font font size=-1 face="Arial" color="gray"><i>Francesco Iorio, EMBL - European Bioinformatics Institute','</i></font><br></center>\n')
  
  timestamp<-paste('<center><font size=-2 face="Arial" color="gray"> result parsing started on ',Sys.time(),'</font></center>',sep='')
  
  contact<-paste('<center><font font size=-1 face="Arial"><a href="mailto:iorio@ebi.ac.uk?Subject=from%20GDSC%20Data-Package%20user" target="_top">iorio@ebi.ac.uk</a></font></center>')
  
  
  if (GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType=='CS'){
    ad<-paste(  GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES,'specific')
  }else{
    ad<-'panCancer'
  }
  
  owner<-paste('<br><br><center><font size=+1 face="Arial">Drug Domain: <b>',DrugDomain,'</b>, Analysis Domain: <b>',ad,'</b></font><br></center><br><br>',sep='')
  
  
  body<-paste(logos,superheader,author,timestamp,contact,owner)
  
  
  
  analysis_info<-gdscANOVA_RP_retrieveAnalysisInfo(superSet,paste(unparsed_results_dir,superSet,'/SYSTEM_INFOS.rdata',sep=''))
  body<-paste(body,analysis_info)
  
  title<-paste('<br><font font size=+1 face="Arial"><b><a href=" ./',URLsuperSet,'/DATA/HTML_elements/000_Significant_Hits.html" target="_blank">
  Explore all the significant associations</a></b></font><br>\n',sep='')
  
  body<-paste(body,title)
  
  link1<-paste('<br><font font face="Arial"><a href="',' ./',URLsuperSet,'/DATA/INPUT/ANOVA_input.txt',
               '">ANOVA input data as tab delimited txt file (genomic features and drug response) </a></font><br>\n',sep='')
  
  body<-paste(body,link1)
  
  link1<-paste('<br><font font face="Arial"><a href="',paste(' ./',URLsuperSet,'/DATA/INPUT/DRUG_DECODE.txt',sep=''),
               '">Drug annotations as tab delimited txt file </a></font><br>\n')
  
  body<-paste(body,link1)
  
  link1<-paste('<br><font font face="Arial"><a href="',paste(' ./',URLsuperSet,'/DATA/INPUT/GENOMIC_REGIONS_DECODE.txt',sep=''),
               '">Copy number altered chromosomal region annotations as tab delimited txt file </a></font><br><hr>\n')
  
  body<-paste(body,link1)
  
  
  idxs<-which(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH &
                as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)
  
  if(length(idxs)>2){
    results_summary<-gdscANOVA_RP_analysis_summary(superSet=superSet,DrugDomain=DrugDomain)
    body<-paste(body,results_summary)
    
  }
  
  results_exploration<-gdscANOVA_RP_results_browsing(superSet=superSet,DrugDomain=DrugDomain)
  body<-paste(body,results_exploration)
  
  write(body, file = paste(PATH,'Home.html',sep=''),append=FALSE)
}


gdscANOVA_RP_retrieveAnalysisInfo<-function(superSet,ANALYSIS_SYSTEMS_INFOS_FN){
  
  load(ANALYSIS_SYSTEMS_INFOS_FN)
  load(paste(GDSCANOVA_SETTINGS$gdscANOVA.results.dir,superSet,'/INPUT/','InputFeatures.rdata',sep=''))
  
  title<-paste('<br><br><font size=+1 face="Arial"><b>Super Analysis Details</b></font><br>\n')
  
  screeningVersion<-paste('<br><font font size=-2 face="Arial">Screening Version:',ANALYSIS_SYSTEMS_INFOS$SCREENING_VERSION,'</font><br>\n')
  drugDomain<-paste('<font size=-2 face="Arial">Drug Domain:',ANALYSIS_SYSTEMS_INFOS$DRUG_DOMAIN,'</font><br>\n')
  time<-paste('<font size=-2 face="Arial">Time:',ANALYSIS_SYSTEMS_INFOS$TIME,'</font><br>\n')
  machine<-paste('<font size=-2 face="Arial">Machine:',ANALYSIS_SYSTEMS_INFOS$MACHINE,'</font><br>\n')
  user<-paste('<font size=-2 face="Arial">User:',ANALYSIS_SYSTEMS_INFOS$USER,'</font><br>\n')
  
  
  total_number_of_performed_tests<-ANALYSIS_SYSTEMS_INFOS$N_FEASIBLE_TESTS
  perc_of_performed_tests<-ANALYSIS_SYSTEMS_INFOS$PERC_FEAS_TESTS
  
  
  ntest<-paste('<br><font size=-1 face="Arial">Total Number of ANOVA tests performed: ',total_number_of_performed_tests,
               ' (',format(perc_of_performed_tests,digits=2),'% of total number of drug/feature combos)','</font><br>\n',sep='')
    
  total_number_drugs<-length(unique(TOTRES[,'Drug id']))
  ndrug<-paste('<font size=-1 face="Arial">Total number of tested drugs:',total_number_drugs,'</font><br>\n')
  
  total_number_features<-length(unique(TOTRES[,'FEATURE']))
  nfeat<-paste('<font size=-1 face="Arial">Total number of tested genomic features (mutated driver genes and copy number altered genomic regions):',total_number_features,'</font><br>\n')
  
  total_number_cell_lines<-length(intersect(colnames(InputFeatures$BEM),rownames(IC50s)))
  ncell<-paste('<font size=-1 face="Arial">Total number of screened cell lines:',total_number_cell_lines,'</font><br><br>')
  
  MSIinclusion<-paste('<font size=-1 face="Arial">MicroSatellite instability included as factor = ',
                      as.logical(GDSCANOVA_SETTINGS$gdscANOVA.settings.includeMSI_Factor),'</font><br><br>\n<hr>')
  
    
  total_number_of_significant_ass<-length(which(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH &
                                                  as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH))
  
  total_number_of_significant_sens_ass<-length(which(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH &
                                                  as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH &
                                                 as.numeric(TOTRES[,"FEATURE_deltaMEAN_IC50"])<0))

  total_number_of_significant_res_ass<-length(which(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH &
                                                     as.numeric(TOTRES[,"ANOVA FEATURE FDR %"]<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH) &
                                               as.numeric(TOTRES[,"FEATURE_deltaMEAN_IC50"])>0))


  
  nsigHits<-paste('<font size=-1 face="Arial">Total number of significant associations: ',total_number_of_significant_ass,' (',total_number_of_significant_sens_ass,
' for sensitivity and ', total_number_of_significant_res_ass,' for resistance)</font><br><br>\n',sep='')
  
  thresholds<-paste('<font size=-1 face="Arial">p-value significance threshold: ',GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH,
                    ', % FDR significance threshold: ',  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH,'</font><br><br>\n',sep='')

  pp<-as.numeric(TOTRES[which(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH &
                     as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH),"FEATURE_ANOVA_pval"])
  
  minPval<-format(min(pp),digits=3,scientific = TRUE)
  maxPval<-format(max(pp),digits=3,scientific = TRUE)
  pvalranges<-paste('<font size=-1 face="Arial">range of significant p-values: [',minPval,', ',maxPval,']','</font><br><br>\n',sep='')
  
  pp<-as.numeric(TOTRES[which(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH &
                              as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<  GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH),"ANOVA FEATURE FDR %"])

  minPval<-format(min(pp),digits=3,scientific = FALSE)
  maxPval<-format(max(pp),digits=3,scientific = FALSE)
  
  fdrranges<-paste('<font size=-1 face="Arial">range of significant % FDRs: [',minPval,', ',maxPval,']','</font><br><br>\n<hr>',sep='')

  body<-paste(title,screeningVersion,drugDomain,time,machine,user,ntest,ndrug,nfeat,ncell,MSIinclusion,nsigHits,thresholds,pvalranges,fdrranges)
    
  return(body)
  
}

gdscANOVA_RP_analysis_summary<-function(superSet,DrugDomain){
  
  URLsuperSet<-str_replace_all(superSet,' ','%20')
  
  title<-paste('<br><br><font size=+1 face="Arial"><b>Results summaries</b></font><br>\n')
  
  vpall<-paste('<img src= ./',URLsuperSet,'/DATA/OUTPUT/comprehensive_volcanoPlot.png width=510 height=680 title=comp_vp align="center"/>',sep='')
  
  volcanoPlots<-paste(vpall,'\n')
  
  bodytext<-paste('<br><br><font font size=+1 face="Arial"><b>- Comprehensive volcano plots</b></font><br>\n')
  
  body<-paste(title,bodytext,volcanoPlots)
  
  
  vpall<-paste('<img src= ./',URLsuperSet,'/DATA/OUTPUT/FeatureSummary.png width=429 height=660 title=comp_vp align="center"/>',sep='')
  featureSummaries<-paste(vpall,'\n')
  
  bodytext<-paste('<br><br><br><font font size=+1 face="Arial"><b>- Features most frequently associated with drug response:</b></font><br><br>\n')
  
  body<-paste(body,bodytext,featureSummaries)
  
  link1<-paste('<br><br><font font size=-1 face="Arial"><a href="',paste(' ./',URLsuperSet,'/DATA/OUTPUT/FeatureSummary.txt',sep=''),
               '">right click on this link to save summary infos for all the features as tab delimited txt file </a></font><br>\n')
  
  body<-paste(body,link1)
  
  
  vpall<-paste('<img src= ./',URLsuperSet,'/DATA/OUTPUT/DrugSummary.png width=500 height=500 title=comp_vp align="center"/>',sep='')
  
  featureSummaries<-paste(vpall,'\n')
  
  bodytext<-paste('<br><br><br><font font size=+1 face="Arial"><b>- Drugs whose response is frequently associated with a feature:
                  </b></font><br><br>\n')
  
  body<-paste(body,bodytext,featureSummaries)
  
  link1<-paste('<br><br><font font size=-1 face="Arial"><a href="',paste('./',URLsuperSet,'/DATA/OUTPUT/DrugSummary.txt',sep=''),
               '">right click on this link to save summary infos for all the features as tab delimited txt file </a></font><br>\n')
  
  body<-paste(body,link1)
  
  
  return(body)
}


gdscANOVA_RP_results_browsing<-function(superSet=superSet,DrugDomain=DrugDomain){
  
   URLsuperSet<-str_replace_all(superSet,' ','%20')
   body<-''
   
   drug_ids_panp<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                         DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
   
   did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
   idxs<-which(is.element(did_in_totres,drug_ids_panp))
   
   redTOTRES<-TOTRES[idxs,]
   
   gdscANOVA_RP_all_assoc_summary(superSet=superSet,redTOTRES=redTOTRES,PATH=(packages_DD_DATA_HTMLEL_dir))
   
   body<-paste(body,'<br><hr>')
   
   bunch<-'<br><font font size=+1 face="Arial"><b>Drug-wise associations browser</b></font><br>\n'
   
   body<-paste(body,bunch)
   
   drug_ids_panp<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                         DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])
  
  didintotres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(2,nrow(TOTRES)*2,2)]

  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  idxs<-which(is.element(did_in_totres,drug_ids_panp))
  
  DIDS<-sort(unique(TOTRES[idxs,"Drug id"]))
  
  
  crunk<-''
  for (i in 1:length(DIDS)){
    
    voce<-paste(DIDS[i],' - ',DRUG_PROPS[str_split(DIDS[i],'_')[[1]][1],'DRUG_NAME'],sep='')
    bunch<-paste('<center><font face="Arial" font size=-1><a href="./',URLsuperSet,'/DATA/HTML_elements/DRUGS/',DIDS[i],'.html" target="_blank">',voce,'</a></font></center>\n',sep='')
    crunk<-paste(crunk,bunch)
  }
  
  body<-paste(body,crunk)
  
  
  body<-paste(body,'<br><hr>')
  
  
  bunch<-'<br><font font size=+1 face="Arial"><b>Feature-wise associations browser</b></font><br>\n'
  
  body<-paste(body,bunch)
  
  drug_ids_panp<-as.character(DRUG_BY_COMPANIES[which(DRUG_BY_COMPANIES[,DrugDomain]==1 |
                                                        DRUG_BY_COMPANIES[,"Web Released"]==1),'DRUG_ID'])

  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  
  did_in_totres<-unlist(str_split(TOTRES[,'Drug id'],'_'))[seq(1,nrow(TOTRES)*2,2)]
  idxs<-which(is.element(did_in_totres,drug_ids_panp))
  
  FF<-sort(unique(TOTRES[idxs,"FEATURE"]))
  
  
  crunk<-''
  for (i in 1:length(FF)){
    
    voce<-FF[i]
    bunch<-paste('<center><font face="Arial" font size=-1><a href="./',URLsuperSet,'/DATA/HTML_elements/FEATURES/',FF[i],'.html" target="_blank">',voce,'</a></font></center>\n',sep='')
    crunk<-paste(crunk,bunch)
  }
  
  body<-paste(body,crunk)
  
  
  body<-paste(body,'<br><hr>')
  
  
 return(body)
  
}











