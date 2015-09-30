source('R/FI.GDSC.ANOVA_ind_pack_creation_library.R')

DRUG_DOMAINS<-gdscANOVA.settings.DRUG_domain

packages_dir<-paste(current_dir,'PARSED_RESULTS/',sep='')

if (!file.exists(packages_dir)){
  dir.create(packages_dir)
}

DrugDomain<-DRUG_DOMAINS
 
print(paste('Creating result package for drug domain:',DrugDomain))
  
print('creating folder structure')
  
packages_DD_dir<-paste(packages_dir,DrugDomain,'/',sep='')
  
if (!file.exists(packages_DD_dir)){
  dir.create(packages_DD_dir)  
}

packages_DD_DATA_dir<-paste(packages_DD_dir,'DATA/',sep='')

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
gdscANOVA_RP_create_individual_associations_html(current_dir)

print('creating individual drug htmls')
gdscANOVA_RP_create_individual_drug_html(superSet = current_dir)

print('creating individual features htmls')
gdscANOVA_RP_create_individual_feature_html(current_dir)
 
print('creating MANOVA input/output files') 
gdscANOVA_RP_create_MANOVA_INPUT_OUTPUT_files(current_dir)

print('creating DRUG DECODE files') 
gdscANOVA_RP_create_DRUG_DECODE_file()

# 
#print('creating GENOMIC REGIONS DECODE files') 
#create_GenomicRegions_DECODE_file()
# #   

print('creating comprehensive volcano plots') 
gdscANOVA_RP_create_comprehensive_vp()
 
idxs<-which(as.numeric(TOTRES[,"FEATURE_ANOVA_pval"])<gdscANOVA.settings.pval_TH &
              as.numeric(TOTRES[,"ANOVA FEATURE FDR %"])<gdscANOVA.settings.FDR_TH)

if (length(idxs)>2){
  print('creating drug/feature summaries') 
  gdscANOVA_RP_create_summaries()
}
 
print('creating HTML index page')
gdscANOVA_RP_start_page_creation(superSet=current_dir,packages_DD_dir,DrugDomain)
