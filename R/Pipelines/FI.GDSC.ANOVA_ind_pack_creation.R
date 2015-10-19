

library(stringr)

source('R/FI.GDSC.ANOVA.Preamble_Library.R')
source('R/FI.GDSC.ANOVA.Graphic_Library.R')

#source('R/FI.GDSC.ANOVA_ind_pack_creation_library.R')

destination_dir<-'../../RESULTS/ANOVA/PARSED_PACKAGES/'

unparsed_results_dir<-'../../RESULTS/ANOVA/'

settingFile_dir<-'../../Data/ANOVA_setting_files/all_Analyses_Oct2015/'

performed_analyses<-setdiff(dir(unparsed_results_dir),'PARSED_PACKAGES')

if (!file.exists(destination_dir)){
  dir.create(destination_dir)
}

finaldestination_dir<-paste(destination_dir,DRUG_DOMAIN,sep='')

if (!file.exists(finaldestination_dir)){
  dir.create(finaldestination_dir)
}


DrugDomain<-DRUG_DOMAIN
 
print(paste('Creating result package for drug domain:',DrugDomain))
  


for (currentAnalysis in performed_analyses){
  
  print(paste('Generating analysis Specific Package:',currentAnalysis))
  
  GDSCANOVA_SETTINGS<-gdscANOVA_Preamble(ANOVA_setting_file = paste(settingFile_dir,currentAnalysis,'.csv',sep=''))
  
  load(paste(GDSCANOVA_SETTINGS$gdscANOVA.results.dir,currentAnalysis,'/INPUT/','InputFeatures.rdata',sep=''))
  source('R/FI.GDSC.ANOVA_ind_pack_creation_library.R')
  
  load(paste(unparsed_results_dir,currentAnalysis,'/OUTPUT/ANOVA_results.rdata',sep=''))
  
  
  current_dir<-currentAnalysis
  
  print('creating folder structure')
  
  packages_DD_dir<-paste(destination_dir,DrugDomain,'/',currentAnalysis,sep='')
  
  if (!file.exists(packages_DD_dir)){
    dir.create(packages_DD_dir)  
  }
  
  packages_DD_DATA_dir<-paste(packages_DD_dir,'/DATA/',sep='')
  
  if (!file.exists(packages_DD_DATA_dir)){
    dir.create(packages_DD_DATA_dir)
  }
  
  packages_DD_DATA_INPUT_dir<-paste(packages_DD_DATA_dir,'INPUT/',sep='')
  
  if (!file.exists(packages_DD_DATA_INPUT_dir)){
    dir.create(packages_DD_DATA_INPUT_dir)
  }
  
  packages_DD_DATA_OUTPUT_dir<-paste(packages_DD_DATA_dir,'OUTPUT/',sep='')
  
  if (!file.exists(packages_DD_DATA_OUTPUT_dir)){
    dir.create(packages_DD_DATA_OUTPUT_dir)
  }
  
  packages_DD_DATA_HTMLEL_dir<-paste(packages_DD_DATA_dir,'HTML_elements/',sep='')
  
  if (!file.exists(packages_DD_DATA_HTMLEL_dir)){
    dir.create(packages_DD_DATA_HTMLEL_dir)
  }
  
  packages_DD_DATA_HTMLEL_ASSOC_dir<-paste(packages_DD_DATA_HTMLEL_dir,'associations/',sep='')
  
  if (!file.exists(packages_DD_DATA_HTMLEL_ASSOC_dir)){
    dir.create(packages_DD_DATA_HTMLEL_ASSOC_dir)
  }
  
  packages_DD_DATA_HTMLEL_DRUG_dir<-paste(packages_DD_DATA_HTMLEL_dir,'DRUGS/',sep='')
  
  if (!file.exists(packages_DD_DATA_HTMLEL_DRUG_dir)){
    dir.create(packages_DD_DATA_HTMLEL_DRUG_dir)
  }
  
  packages_DD_DATA_HTMLEL_FEATURE_dir<-paste(packages_DD_DATA_HTMLEL_dir,'FEATURES/',sep='')
  
  if (!file.exists(packages_DD_DATA_HTMLEL_FEATURE_dir)){
    dir.create(packages_DD_DATA_HTMLEL_FEATURE_dir)
  }
  
  packages_DD_DATA_HTMLEL_images_dir<-paste(packages_DD_DATA_HTMLEL_dir,'IMAGES/',sep='')
  
  if (!file.exists(packages_DD_DATA_HTMLEL_images_dir)){
    dir.create(packages_DD_DATA_HTMLEL_images_dir)
  }
  
  
  print('copying HTML elements')
  gdscANOVA_RP_copy_html_elements(packages_DD_DATA_HTMLEL_images_dir)
  
  print('creating individual associations htmls')
  gdscANOVA_RP_create_individual_associations_html(superSet = current_dir)
  
  print('creating individual drug htmls')
  gdscANOVA_RP_create_individual_drug_html(superSet = current_dir)
  
  print('creating individual features htmls')
  gdscANOVA_RP_create_individual_feature_html(superSet = current_dir)
  
  print('creating MANOVA input/output files') 
  gdscANOVA_RP_create_MANOVA_INPUT_OUTPUT_files(current_dir)
  
  ###################################################################
  
  print('creating DRUG DECODE files') 
  gdscANOVA_RP_create_DRUG_DECODE_file()
  
   
  print('creating GENOMIC REGIONS DECODE files') 
  gdscANOVA_RP_create_GenomicRegions_DECODE_file()
  
  print('creating comprehensive volcano plots') 
  gdscANOVA_RP_create_comprehensive_vp()
  
  idxs<-which(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.pval_TH &
                as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<GDSCANOVA_SETTINGS$gdscANOVA.settings.FDR_TH)
  
  if (length(idxs)>10){
    print('creating drug/feature summaries') 
    gdscANOVA_RP_create_summaries()
  }
  
  print('creating HTML index page')
  gdscANOVA_RP_start_page_creation(superSet=current_dir,PATH = packages_DD_dir,packageName = DrugDomain)
  
  if (currentAnalysis==performed_analyses[1]){
    finalDiag<-gdscANOVA_RP_finalDiagnosis(TOTRES)  
  }else{
    finalDiag<-rbind(finalDiag,gdscANOVA_RP_finalDiagnosis(TOTRES))
  }
}

rownames(finalDiag)<-performed_analyses

URLdir<-str_replace_all(current_dir,' ','%20')
URLDRUG_DOMAIN<-str_replace_all(DRUG_DOMAIN,' ','%20')

logo1<-paste('<img src=','./',URLDRUG_DOMAIN,'/',URLdir,'/DATA/HTML_elements/IMAGES/sanger-logo.png',' title=sanger-logo align="center"/>',sep='')
logo3<-paste('<img src=','./',URLDRUG_DOMAIN,'/',URLdir,'/DATA/HTML_elements/IMAGES/EBI_logo.png',' title=sanger-logo align="center"/>',sep='')

logos<-paste(logo1,logo3,'\n')

superheader <- paste('<br><br><center><font size=+2 face="Arial">','Result package from the GDSC1000 project - Collaborator: ',DRUG_DOMAIN,'</font><br></center>\n')

author<-paste('<br><center><font font size=-1 face="Arial" color="gray"><i>Francesco Iorio, EMBL - European Bioinformatics Institute','</i></font><br></center>\n')

timestamp<-paste('<center><font size=-2 face="Arial" color="gray"> result parsing started on ',Sys.time(),'</font></center>',sep='')

contact<-paste('<center><font font size=-1 face="Arial"><a href="mailto:iorio@ebi.ac.uk?Subject=from%20GDSC%20Data-Package%20user" target="_top">iorio@ebi.ac.uk</a></font></center>')

if (GDSCANOVA_SETTINGS$gdscANOVA.settings.analysisType=='CS'){
  ad<-paste(  GDSCANOVA_SETTINGS$gdscANOVA.settings.CELL_LINES,'specific')
}else{
  ad<-'panCancer'
}

owner<-paste('<br><br><center><font size=+1 face="Arial">Drug Domain: <b>',DrugDomain,'</b> proprietary compounds and public compounds</font><br></center><br><br>',sep='')

header <- paste('<center><font size=+2 face="Arial">Individual ANOVAs Results Summary</font><br><center>\n
                  <center><table border=1 width=100% cellspacing=1 cellpadding=2 cols=7 style="font-family: Arial; font-size: 13px">\n
                  <tr><td><b><center>Analysis name:</center></b></td>
                  <td><b><center>Number of hits:</center></b></td>
                  <td><b><center>n. Involved Proprietary compounds:</center></b></td>
                  <td><b><center>out of:</center></b></td>
                  <td><b><center>n. Involved Public compounds:</center></b></td>
                  <td><b><center>out of:</center></b></td></tr>\n')

tail<-paste('</tr>\n</table></center>')

body<-'\n'

body<-paste(body,logos,superheader,author,timestamp,contact,owner,header)


for (i in 1:nrow(finalDiag)){
  currentLine<-paste('<td><center><a href="./',DRUG_DOMAIN,'/',paste(rownames(finalDiag)[i],'Home.html',sep=''),'" target="_blank">',rownames(finalDiag)[i],'</a></center></td>\n',sep='')
  body<-paste(body,currentLine)
  currentLine<-paste('<td><center>',finalDiag[i,1],'<center></td>\n',sep='')
  body<-paste(body,currentLine)
  currentLine<-paste('<td><center>',finalDiag[i,2],'<center></td>\n',sep='')
  body<-paste(body,currentLine)
  currentLine<-paste('<td><center>',finalDiag[i,3],'<center></td>\n',sep='')
  body<-paste(body,currentLine)
  currentLine<-paste('<td><center>',finalDiag[i,4],'<center></td>\n',sep='')
  body<-paste(body,currentLine)
  currentLine<-paste('<td><center>',finalDiag[i,5],'<center></td></tr>\n',sep='')
  body<-paste(body,currentLine)
}

body<-paste(body,tail)
write(body, file = paste(destination_dir,'/',DRUG_DOMAIN,'_home.html',sep=''),append=FALSE)



