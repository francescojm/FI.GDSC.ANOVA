
gdscANOVA.settings.CELL_LINES<-'PANCAN'
gdscANOVA.settings.additionalFeatures<-FALSE

source('R/FI.GDSC.ANOVA.Preamble.R')


logo3<-paste('<img src=',paste('../BLCA_GDSC1000_paper_set_2015-09-28_21-20-38//PARSED_RESULTS//GDSC1000_paper_set//DATA//HTML_elements//IMAGES//EBI_logo.png',sep=''),' title=sanger-logo align="center"/>')

logos<-''

# if (gdscANOVA.settings.resPackageIncludeSangerLogo){
#   logos<-paste(logos,logo1,'\n') 
# }
#  if (gdscANOVA.settings.resPackageIncludeNKILogo){
#    logos<-paste(logos,logo2,'\n') 
# }
  logos<-paste(logos,logo3,'\n') 



# superheader <- paste('<br><br><center><font size=+2 face="Arial">',
#                      gdscANOVA.settings.resPackageHeader,
#                      '</font><br></center>\n',sep='')

superheader <- paste('<br><br><center><font size=+2 face="Arial">',
                     'DR.Lancet all ANOVA results index',
                     '</font><br></center>\n',sep='')

author<-paste('<br><center><font font size=-1 face="Arial" color="gray"><i>Francesco Iorio, EMBL - European Bioinformatics Institute','</i></font><br></center>\n')

timestamp<-paste('<center><font size=-2 face="Arial" color="gray"> result assembling started on ',Sys.time(),'</font></center>',sep='')

contact<-paste('<center><font font size=-1 face="Arial"><a href="mailto:iorio@ebi.ac.uk?Subject=from%20GDSC%20Data-Package%20user" target="_top">iorio@ebi.ac.uk</a></font></center>')

body<-paste(logos,superheader,author,timestamp,contact)

fn<-dir('../../RESULTS/DR.Lancet//ANOVA/')

fn<-setdiff(fn,"all_CS_packages")


bunch<-paste('<center><font face="Arial" font size=+2>Cancer Specific Analyses:</a></font></center>\n',sep='')

body<-paste(body,'<br><br>',bunch)
crunk<-''

for (i in 1:length(fn)){   
  bunch<-paste('<center><font face="Arial" font size=+1><a href="../',fn[i],'/PARSED_RESULTS/GDSC1000_paper_set//Home.html" target="_blank">',fn[i],' (',fn[i],')</a></font></center>\n',sep='')
  crunk<-paste(crunk,'<br>',bunch)
}

body<-paste(body,'<br>',crunk)

write(body, file = paste('../../RESULTS/DR.Lancet//ANOVA//all_CS_packages/','Home.html',sep=''),append=FALSE)




