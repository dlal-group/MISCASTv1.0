library(shinydashboard)
library(shiny)
library(shinyjs)
library(ggplot2)
library(readr)
library(DT)
library(ggrepel)
library(tm)
library(stringr)
library(shinythemes)
library(dplyr)
library(plyr)
library(tidyr)
library(scales)
library(ggpubr)
library(shinycssloaders)
library(plotly)
library(ggExtra)
library(ggalt)

jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

jscode <- "
shinyjs.runMolart = function(params) { 
runMolart(params);
}

shinyjs.runMsa = function(params) { 
runMsa(params);
}
"
#a global variable
reload_var <- FALSE

molartContainerId <- "molartContainer"
#protaelContainerResearchId <- "protaelContainerAAwFResearch"
msaContainerId <- "msaContainer"


## data load -- start ##
#table for genes with protein subclass name and whether it is there as controlled missense or patient missense mutation

pmissOrgmissGene <- read.csv("csvFiles/Gene_Info.csv")
protein_class_def <- read.delim("csvFiles/protein_class_def_genes_counts.txt", header = T, sep = "\t", stringsAsFactors = F, colClasses = "character")

#one category table for all
psuperclass_wise_proteins <- read_delim("data/psuperclass_wise_proteins_prots_genes_with_name.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
psuperclass_wise_proteins$category <- "By Protein Major Class"
func_wise_propteins <- rbind(psuperclass_wise_proteins)

# AA - SS - ASA - PTM property file
gmiss_dssp_ptm <- read.delim("data/gmiss_prp_164915.txt", header = T, sep = "\t", stringsAsFactors = F, colClasses = "character")
pmiss_dssp_ptm <- read.delim("data/pmiss_prp_32924.txt", header = T, sep = "\t", stringsAsFactors = F, colClasses = "character")

gmiss_prp_all <- read_delim("data/gmiss_4907_gene_492755_sav_prop.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
pmiss_dssp_prp <- read_delim("data/pmiss_1330_gene_32923_sav_prop.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

path_genes <- data.frame(unique(pmiss_dssp_prp$geneName))
colnames(path_genes) <- "genename"
gmiss_dssp_prp <- subset(gmiss_prp_all, geneName %in% path_genes$genename)

gmiss_dssp_uniprot <- read_delim("data/gmiss_1330_gene_164915_sav_uniprot_comb.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
pmiss_dssp_uniprot <- read_delim("data/pmiss_1330_gene_32923_sav_uniprot_comb.txt", "\t", escape_double = FALSE, trim_ws = TRUE)


# gene- protein- transcript
GPTLPC_list <- read.delim("csvFiles/GPTLPC_1330.txt", header = T, sep = "\t", stringsAsFactors = F, colClasses = "character")
gene_1330_cv <- as.character(GPTLPC_list[,1])
feature_index_df <- read.delim("csvFiles/feature_index_to_name.txt", header = T, sep = "\t", stringsAsFactors = F, colClasses = "character")

aa_info <- read.delim("csvFiles/aa_info.txt", header = T, sep = "\t", stringsAsFactors = F)

PC24_PF40_pathogenic <- read_delim("data/PC24_PF40_pathogenic.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(PC24_PF40_pathogenic)[1] <- "pcindex"
colnames(PC24_PF40_pathogenic)[2] <- "pcname"
colnames(PC24_PF40_pathogenic)[3] <- "pcnamelower"
PC24_PF40_population <- read_delim("data/PC24_PF40_population.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(PC24_PF40_population)[1] <- "pcindex"
colnames(PC24_PF40_population)[2] <- "pcname"
colnames(PC24_PF40_population)[3] <- "pcnamelower"

bond_type_gmiss <- read.csv("csvFiles/bond_type_gmiss.csv")
bond_type_pmiss <- read.csv("csvFiles/bond_type_pmiss.csv")
dist_gmiss <- read.csv("csvFiles/dist_gmiss_new.csv")
dist_pmiss <- read.csv("csvFiles/dist_pmiss_new.csv")

## data load -- end ##
gnomAD_link = paste(c("http://gnomad.broadinstitute.org/", ""), collapse = "")
ClinVar_link = paste(c("https://www.ncbi.nlm.nih.gov/clinvar/", ""), collapse = "")
HGMD_link = paste(c("http://www.hgmd.cf.ac.uk/ac/index.php", ""), collapse = "")
SC_link = paste(c("https://www.broadinstitute.org/stanley", ""), collapse = "")
Broad_link = paste(c("https://www.broadinstitute.org/", ""), collapse = "")
PDB_link = paste(c("https://www.rcsb.org/", ""), collapse = "")
DSSP_link = paste(c("https://swift.cmbi.umcn.nl/gv/dssp/", ""), collapse = "")
PDBsum_link = paste(c("http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=index.html", ""), collapse = "")
SIFTS_link = paste(c("https://www.ebi.ac.uk/pdbe/docs/sifts/", ""), collapse = "")
PANTHER_link = paste(c("http://www.pantherdb.org/panther/ontologies.jsp?", ""), collapse = "")
PPS_link = paste(c("https://www.phosphosite.org/homeAction", ""), collapse = "")
Uniprot_link = paste(c("https://www.uniprot.org/", ""), collapse = "")


############################################################################
############################################################################
############################################################################
# Define server logic required to draw a histogram

shinyServer(function(input, output,session) {
  
  observeEvent(input$upperlayer,{
    if(input$upperlayer == "page1"){
      if(reload_var == TRUE){
        reload_var <<- FALSE     #the change will be visible from other sessions
        js$reset() # this should restart the app
      }
    }
    else{
      reload_var <<- TRUE
    }
  })
  
  hideTab(inputId = "mainTabset", target = "Protein class wise characterization of variants")
  
  
  observeEvent(input$submit,{
    if(input$homeSideBarTabSetPanel == 'Select a Gene')
    {
      if(input$geneSelected == '' ){
        showModal(modalDialog(
          title = "Warning Message",
          div(tags$p("Not a Valid Option!", style = "text-align: center;")),
          size = 's',
          easyClose = TRUE
          
        ))
      }
    }
    else if(input$homeSideBarTabSetPanel == 'Protein Class')
    {
      if(input$pclassNameselected == '' ){
        showModal(modalDialog(
          title = "Warning Message",
          div(tags$p("Not a Valid Option!", style = "text-align: center;")),
          size = 's',
          easyClose = TRUE
          
        ))
      }
    }
  })
  
  observeEvent(input$reportTsubmit,{
    #if(input$homeSideBarTabSetPanel == 'Select a Gene' || input$homeSideBarTabSetPanel == 'Protein Class')
    #{
      if(input$reportTgeneSelected == ''){
        showModal(modalDialog(
          title = "Warning Message",
          div(tags$p("Not a Valid Option!", style = "text-align: center;")),
          size = 's',
          easyClose = TRUE
          
        ))
      }
    #}
  })
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel != 'Protein Class')
    {
      hideTab(inputId = "mainTabset", target = "Protein class wise characterization of variants")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Protein Class')
    {
      showTab(inputId = "mainTabset", target = "Protein class wise characterization of variants")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Protein Class')
    {
      hideTab(inputId = "mainTabset", target = "2D visualization")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Protein Class')
    {
      hideTab(inputId = "mainTabset", target = "1D visualization")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Protein Class')
    {
      hideTab(inputId = "mainTabset", target = "Feature table")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Protein Class')
    {
      hideTab(inputId = "mainTabset", target = "3D visualization")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Protein Class')
    {
      hideTab(inputId = "mainTabset", target = "Index table")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene')
    {
      showTab(inputId = "mainTabset", target = "2D visualization")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene')
    {
      showTab(inputId = "mainTabset", target = "1D visualization")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene')
    {
      showTab(inputId = "mainTabset", target = "Feature table")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene')
    {
      showTab(inputId = "mainTabset", target = "3D visualization")
    }
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene')
    {
      showTab(inputId = "mainTabset", target = "Index table")
    }
  })
  
  observeEvent(input$submit, {
    updateNavbarPage(session, "mainTabset",
                     selected = "Information")
  })
  
  output$sidebarHeader = renderUI(
    {
      tags$b("Selection")
    }
  )
  
  # Report Track / Research Track#
  observeEvent(input$researchB, {
    updateNavbarPage(session, "upperlayer",
                     selected = "researchTrack")
  })
  
  
  observeEvent(input$reportB, {
    updateNavbarPage(session, "upperlayer",
                     selected = "reportTrack")
  })
  # Report Track / Research Track#
  
  # Report Track information show up#
  observeEvent(
    input$reportTsubmit,{
      if(input$reportTgeneSelected != ""){
        shinyjs::show("report_INFOwellpanel")
      }
      else{
        shinyjs::hide("report_INFOwellpanel")
        
      }
    }
  )
  
  observeEvent(
    input$aa_Selected,{
      
      if(input$aa_Selected!= ""){
        shinyjs::show("below_report_INFOwellpanel")
      }
    }
  )
  
  observeEvent(
    input$aa_Selected_aaFeat,{
      
      if(input$aa_Selected_aaFeat!= ""){
        shinyjs::show("below_report_aaFeat_wellPanel")
      }
    }
  )
  observeEvent(input$trackButton,{
    if(input$trackInput != ''){
      shinyjs::hide(id  = "queryHide")
      if(input$trackInput == "Variant Summary Report"){
        updateNavbarPage(session, "upperlayer", selected = "singleTrack")
        shinyjs::show(id = "reportTrackHide")
        shinyjs::hide(id = "researchTrackHide")
        
      }
      else if(input$trackInput == "Variant Analysis Suite"){
        updateNavbarPage(session, "upperlayer", selected = "singleTrack")
        shinyjs::show(id = "researchTrackHide")
        shinyjs::hide(id = "reportTrackHide")
      }
    }
  })
  observeEvent(input$submit,{
    if(input$homeSideBarTabSetPanel == 'Select a Gene')
    {
      if(input$geneSelected != '' ){
        shinyjs::show(id  = "queryHide")
        shinyjs::hide(id = "researchTrackHide")
      }
    }
    else if(input$homeSideBarTabSetPanel == 'Protein Class')
    {
      if(input$pclassNameselected != '' ){
        shinyjs::show(id  = "queryHide")
        shinyjs::hide(id = "researchTrackHide")
      }
    }
    
  })
  observeEvent(input$researchB,{
    updateNavbarPage(session, "upperlayer", selected = "singleTrack")
    shinyjs::show(id = "researchTrackHide")
    shinyjs::hide(id = "reportTrackHide")
    
  })
  observeEvent(input$reportB,{
    updateNavbarPage(session, "upperlayer", selected = "singleTrack")
    shinyjs::show(id = "reportTrackHide")
    shinyjs::hide(id = "researchTrackHide")
  })
  
  observeEvent(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ptm6_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      ptm6_loli_gene_info <- read_delim(ptm6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    
    
    #ptm6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/sf3d_final/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    
    this_gene_path_pop_ptm6 <- ptm6_loli_gene_info[((ptm6_loli_gene_info$`Pathogenic Mutation`!="none") | (ptm6_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_aa_ptm6 <- as.data.frame(this_gene_path_pop_ptm6$`Amino Acid Index`)
    this_gene_acetyl <- as.data.frame(this_gene_path_pop_ptm6$`Acetylation Detail`)
    this_gene_methyl <- as.data.frame(this_gene_path_pop_ptm6$`Methylation Detail`)
    this_gene_gcln <- as.data.frame(this_gene_path_pop_ptm6$`O.GclNAc Detail`)
    this_gene_phos <- as.data.frame(this_gene_path_pop_ptm6$`Phosphorylation Detail`)
    this_gene_sumoy <- as.data.frame(this_gene_path_pop_ptm6$`Sumoylation Detail`)
    this_gene_ubiq <- as.data.frame(this_gene_path_pop_ptm6$`Ubiquitination Detail`)
    this_gene_mut_ptm6 <- cbind(this_gene_aa_ptm6, this_gene_acetyl, this_gene_methyl, this_gene_gcln, this_gene_phos, this_gene_sumoy, this_gene_ubiq)
    colnames(this_gene_mut_ptm6) <- c("aaindex", "acetyl", "methyl", "gcln", "phos", "sumoy", "ubiq")
    mut_count_ptm6 <- nrow(this_gene_mut_ptm6)
    
    for(i in 1:mut_count_ptm6){
      if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ptm6$var_type[i] <- "pop"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] == "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "path"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "both"
      }
    }
    
    #this_gene_mut_ptm6 <- this_gene_mut_ptm6[this_gene_mut_ptm6$bonddetail != "-", ]
    
    this_gene_ptm6_ext <- data.frame(ptmdist = as.numeric(), aaind = as.numeric(), vartype = as.character(), ptmtype = as.character())
    for(i in 1:mut_count_ptm6){
      #print(i)
      this_aaind = this_gene_mut_ptm6$aaindex[i]
      this_vartype = this_gene_mut_ptm6$var_type[i]
      
      this_acetyldetail = as.character(this_gene_mut_ptm6$acetyl[i])
      this_methyldetail = as.character(this_gene_mut_ptm6$methyl[i])
      this_gclndetail = as.character(this_gene_mut_ptm6$gcln[i])
      this_phosdetail = as.character(this_gene_mut_ptm6$phos[i])
      this_sumoydetail = as.character(this_gene_mut_ptm6$sumoy[i])
      this_ubiqdetail = as.character(this_gene_mut_ptm6$ubiq[i])
      
      if((this_acetyldetail == "-") && (this_methyldetail == "-") && (this_gclndetail == "-") && (this_phosdetail == "-") && (this_sumoydetail == "-") && (this_ubiqdetail == "-")){
        this_gene_ptm6_ext_entry = data.frame(-25, this_aaind, this_vartype, "No Annotation")
        colnames(this_gene_ptm6_ext_entry) = c("ptmdist", "aaind", "vartype", "ptmtype")
        this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
      }else{
        if(this_acetyldetail != "-"){
          this_acetyllist = unlist(strsplit(this_acetyldetail, split = "[;/]"))
          this_acetyllist_num = suppressWarnings(as.numeric(this_acetyllist))
          this_acetyllist_num = this_acetyllist_num[!is.na(this_acetyllist_num)]
          this_acetyllist_num_ord = sort(this_acetyllist_num)
          if(length(this_acetyllist_num_ord) > 10){
            this_acetyllist_num_ord = this_acetyllist_num_ord[1:10]
          }
          this_acetyllist_num_df = as.data.frame(this_acetyllist_num_ord)
          colnames(this_acetyllist_num_df) <- c("ptmdist")
          this_acetyllist_num_df$aaind <- this_aaind
          this_acetyllist_num_df$vartype <- this_vartype
          this_acetyllist_num_df$ptmtype <- "Acetylation"
          
          this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_acetyllist_num_df)
        }
        if(this_methyldetail != "-"){
          this_methyllist = unlist(strsplit(this_methyldetail, split = "[;/]"))
          this_methyllist_num = suppressWarnings(as.numeric(this_methyllist))
          this_methyllist_num = this_methyllist_num[!is.na(this_methyllist_num)]
          this_methyllist_num_ord = sort(this_methyllist_num)
          if(length(this_methyllist_num_ord) > 10){
            this_methyllist_num_ord = this_methyllist_num_ord[1:10]
          }
          this_methyllist_num_df = as.data.frame(this_methyllist_num_ord)
          colnames(this_methyllist_num_df) <- c("ptmdist")
          this_methyllist_num_df$aaind <- this_aaind
          this_methyllist_num_df$vartype <- this_vartype
          this_methyllist_num_df$ptmtype <- "Methylation"
          
          this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_methyllist_num_df)
        }
        if(this_gclndetail != "-"){
          this_gclnlist = unlist(strsplit(this_gclndetail, split = "[;/]"))
          this_gclnlist_num = suppressWarnings(as.numeric(this_gclnlist))
          this_gclnlist_num = this_gclnlist_num[!is.na(this_gclnlist_num)]
          this_gclnlist_num_ord = sort(this_gclnlist_num)
          if(length(this_gclnlist_num_ord) > 30){
            this_gclnlist_num_ord = this_gclnlist_num_ord[1:30]
          }
          this_gclnlist_num_df = as.data.frame(this_gclnlist_num_ord)
          colnames(this_gclnlist_num_df) <- c("ptmdist")
          this_gclnlist_num_df$aaind <- this_aaind
          this_gclnlist_num_df$vartype <- this_vartype
          this_gclnlist_num_df$ptmtype <- "O.GclNAc"
          
          this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gclnlist_num_df)
        }
        if(this_phosdetail != "-"){
          this_phoslist = unlist(strsplit(this_phosdetail, split = "[;/]"))
          this_phoslist_num = suppressWarnings(as.numeric(this_phoslist))
          this_phoslist_num = this_phoslist_num[!is.na(this_phoslist_num)]
          this_phoslist_num_ord = sort(this_phoslist_num)
          if(length(this_phoslist_num_ord) > 10){
            this_phoslist_num_ord = this_phoslist_num_ord[1:10]
          }
          this_phoslist_num_df = as.data.frame(this_phoslist_num_ord)
          colnames(this_phoslist_num_df) <- c("ptmdist")
          this_phoslist_num_df$aaind <- this_aaind
          this_phoslist_num_df$vartype <- this_vartype
          this_phoslist_num_df$ptmtype <- "Phosphorylation"
          
          this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_phoslist_num_df)
        }
        if(this_sumoydetail != "-"){
          this_sumoylist = unlist(strsplit(this_sumoydetail, split = "[;/]"))
          this_sumoylist_num = suppressWarnings(as.numeric(this_sumoylist))
          this_sumoylist_num = this_sumoylist_num[!is.na(this_sumoylist_num)]
          this_sumoylist_num_ord = sort(this_sumoylist_num)
          if(length(this_sumoylist_num_ord) > 10){
            this_sumoylist_num_ord = this_sumoylist_num_ord[1:10]
          }
          this_sumoylist_num_df = as.data.frame(this_sumoylist_num_ord)
          colnames(this_sumoylist_num_df) <- c("ptmdist")
          this_sumoylist_num_df$aaind <- this_aaind
          this_sumoylist_num_df$vartype <- this_vartype
          this_sumoylist_num_df$ptmtype <- "Sumoylation"
          
          this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_sumoylist_num_df)
        }
        if(this_ubiqdetail != "-"){
          this_ubiqlist = unlist(strsplit(this_ubiqdetail, split = "[;/]"))
          this_ubiqlist_num = suppressWarnings(as.numeric(this_ubiqlist))
          this_ubiqlist_num = this_ubiqlist_num[!is.na(this_ubiqlist_num)]
          this_ubiqlist_num_ord = sort(this_ubiqlist_num)
          if(length(this_ubiqlist_num_ord) > 10){
            this_ubiqlist_num_ord = this_ubiqlist_num_ord[1:10]
          }
          this_ubiqlist_num_df = as.data.frame(this_ubiqlist_num_ord)
          colnames(this_ubiqlist_num_df) <- c("ptmdist")
          this_ubiqlist_num_df$aaind <- this_aaind
          this_ubiqlist_num_df$vartype <- this_vartype
          this_ubiqlist_num_df$ptmtype <- "Ubiquitination"
          
          this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_ubiqlist_num_df)
        }
      }
    }
    
    this_gene_ptm6_ext_acet = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Acetylation", ]
    this_gene_ptm6_ext_met = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Methylation", ]
    this_gene_ptm6_ext_ogcl = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "O.GclNAc", ]
    this_gene_ptm6_ext_phos = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Phosphorylation", ]
    this_gene_ptm6_ext_sumo = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Sumoylation", ]
    this_gene_ptm6_ext_ubiq = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Ubiquitination", ]
    
    if(nrow(this_gene_ptm6_ext_acet) == 0 
       && nrow(this_gene_ptm6_ext_met) == 0 
       && nrow(this_gene_ptm6_ext_ogcl) == 0 
       && nrow(this_gene_ptm6_ext_phos) == 0 
       && nrow(this_gene_ptm6_ext_sumo) == 0 
       && nrow(this_gene_ptm6_ext_ubiq) == 0)
    {
      shinyjs::hide("PTMtext")
      shinyjs::hide("ptmcaption")
    }
    else
    {
      shinyjs::show("PTMtext")
      shinyjs::show("ptmcaption")
    }
    
    if(nrow(this_gene_ptm6_ext_acet) != 0){
      shinyjs::show("Acet")
    }
    else{
      shinyjs::hide("Acet")
    }
    
    if(nrow(this_gene_ptm6_ext_met) != 0){
      shinyjs::show("Met")
    }
    else{
      shinyjs::hide("Met")
    }
    
    if(nrow(this_gene_ptm6_ext_ogcl) != 0){
      shinyjs::show("Ogcl")
    }
    else{
      shinyjs::hide("Ogcl")
    }
    
    if(nrow(this_gene_ptm6_ext_phos) != 0){
      shinyjs::show("Phos")
    }
    else{
      shinyjs::hide("Phos")
    }
    
    if(nrow(this_gene_ptm6_ext_sumo) != 0){
      shinyjs::show("Sumo")
    }
    else{
      shinyjs::hide("Sumo")
    }
    
    if(nrow(this_gene_ptm6_ext_ubiq) != 0){
      shinyjs::show("Ubiq")
    }
    else{
      shinyjs::hide("Ubiq")
    }}
    
  })
  
  observeEvent(input$submit, {
    
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      up6_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      up6_loli_gene_info <- read_delim(up6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_up6 <- nrow(up6_loli_gene_info)
    
    
    #up6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/miscastapp/miscast01/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_up6 <- nrow(up6_loli_gene_info)
    
    up6_loli_gene_info <- up6_loli_gene_info[((up6_loli_gene_info$`Pathogenic Mutation`!="none") | (up6_loli_gene_info$`Population Mutation`!="none")), ]
    
    this_gene_aa_up6_all <- data.frame(aaindex = as.numeric(), subtype = as.character(), uptype = as.character(), vartype = as.character())
    
    #functional site
    available_funcsite <- length(unique(up6_loli_gene_info$`Functional Site`))
    if(available_funcsite > 1){
      up6_loli_funcsite_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Functional Site` != "-"))
      this_gene_aa_funcsite <- as.data.frame(up6_loli_funcsite_info$`Amino Acid Index`)
      this_gene_funcsite <- as.data.frame(up6_loli_funcsite_info$`Functional Site`)
      this_gene_aa_funcsite <- cbind(this_gene_aa_funcsite, this_gene_funcsite)
      colnames(this_gene_aa_funcsite) <- c("aaindex", "subtype")
      this_gene_aa_funcsite$uptype <- "Functional Site"
      this_gene_aa_funcsite_len <- nrow(this_gene_aa_funcsite)
      
      for(i in 1:this_gene_aa_funcsite_len){
        if((up6_loli_funcsite_info$`Population Mutation`[i] != "none") && (up6_loli_funcsite_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_funcsite$vartype[i] <- "pop"
        }else if((up6_loli_funcsite_info$`Population Mutation`[i] == "none") && (up6_loli_funcsite_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_funcsite$vartype[i] <- "path"
        }else if((up6_loli_funcsite_info$`Population Mutation`[i] != "none") && (up6_loli_funcsite_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_funcsite$vartype[i] <- "both"
        }else{
          this_gene_aa_funcsite$vartype[i] <- "none"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_funcsite)
      
    }
    
    #functional/binding region
    available_bindreg <- length(unique(up6_loli_gene_info$`Functional Binding Region`))
    if(available_bindreg > 1){
      up6_loli_bindreg_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Functional Binding Region` != "-"))
      this_gene_aa_bindreg <- as.data.frame(up6_loli_bindreg_info$`Amino Acid Index`)
      this_gene_bindreg <- as.data.frame(up6_loli_bindreg_info$`Functional Binding Region`)
      this_gene_aa_bindreg <- cbind(this_gene_aa_bindreg, this_gene_bindreg)
      colnames(this_gene_aa_bindreg) <- c("aaindex", "subtype")
      this_gene_aa_bindreg$uptype <- "Functional/Binding Region"
      this_gene_aa_bindreg_len <- nrow(this_gene_aa_bindreg)
      
      for(i in 1:this_gene_aa_bindreg_len){
        if((up6_loli_bindreg_info$`Population Mutation`[i] != "none") && (up6_loli_bindreg_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_bindreg$vartype[i] <- "pop"
        }else if((up6_loli_bindreg_info$`Population Mutation`[i] == "none") && (up6_loli_bindreg_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_bindreg$vartype[i] <- "path"
        }else if((up6_loli_bindreg_info$`Population Mutation`[i] != "none") && (up6_loli_bindreg_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_bindreg$vartype[i] <- "both"
        }else{
          this_gene_aa_bindreg$vartype[i] <- "none"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_bindreg)
      
    }
    
    #sequence mmotif/region
    available_seqmotif <- length(unique(up6_loli_gene_info$`Sequence Motif Region`))
    if(available_seqmotif > 1){
      up6_loli_seqmotif_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Sequence Motif Region` != "-"))
      this_gene_aa_seqmotif <- as.data.frame(up6_loli_seqmotif_info$`Amino Acid Index`)
      this_gene_seqmotif <- as.data.frame(up6_loli_seqmotif_info$`Sequence Motif Region`)
      this_gene_aa_seqmotif <- cbind(this_gene_aa_seqmotif, this_gene_seqmotif)
      colnames(this_gene_aa_seqmotif) <- c("aaindex", "subtype")
      this_gene_aa_seqmotif$uptype <- "Sequence Motif/Region"
      this_gene_aa_seqmotif_len <- nrow(this_gene_aa_seqmotif)
      
      for(i in 1:this_gene_aa_seqmotif_len){
        if((up6_loli_seqmotif_info$`Population Mutation`[i] != "none") && (up6_loli_seqmotif_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_seqmotif$vartype[i] <- "pop"
        }else if((up6_loli_seqmotif_info$`Population Mutation`[i] == "none") && (up6_loli_seqmotif_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_seqmotif$vartype[i] <- "path"
        }else if((up6_loli_seqmotif_info$`Population Mutation`[i] != "none") && (up6_loli_seqmotif_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_seqmotif$vartype[i] <- "both"
        }else{
          this_gene_aa_seqmotif$vartype[i] <- "none"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_seqmotif)
      
    }
    
    #modular domain
    available_moddom <- length(unique(up6_loli_gene_info$`Modular Domain`))
    if(available_moddom > 1){
      up6_loli_moddom_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Modular Domain` != "-"))
      this_gene_aa_moddom <- as.data.frame(up6_loli_moddom_info$`Amino Acid Index`)
      this_gene_moddom <- as.data.frame(up6_loli_moddom_info$`Modular Domain`)
      this_gene_aa_moddom <- cbind(this_gene_aa_moddom, this_gene_moddom)
      colnames(this_gene_aa_moddom) <- c("aaindex", "subtype")
      this_gene_aa_moddom$uptype <- "Modular Domain"
      this_gene_aa_moddom_len <- nrow(this_gene_aa_moddom)
      
      for(i in 1:this_gene_aa_moddom_len){
        if((up6_loli_moddom_info$`Population Mutation`[i] != "none") && (up6_loli_moddom_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_moddom$vartype[i] <- "pop"
        }else if((up6_loli_moddom_info$`Population Mutation`[i] == "none") && (up6_loli_moddom_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_moddom$vartype[i] <- "path"
        }else if((up6_loli_moddom_info$`Population Mutation`[i] != "none") && (up6_loli_moddom_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_moddom$vartype[i] <- "both"
        }else{
          this_gene_aa_moddom$vartype[i] <- "none"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_moddom)
    }
    
    #modified residues
    available_modres <- length(unique(up6_loli_gene_info$`Modified Residues`))
    if(available_modres > 1){
      up6_loli_modres_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Modified Residues` != "-"))
      this_gene_aa_modres <- as.data.frame(up6_loli_modres_info$`Amino Acid Index`)
      this_gene_modres <- as.data.frame(up6_loli_modres_info$`Modified Residues`)
      this_gene_aa_modres <- cbind(this_gene_aa_modres, this_gene_modres)
      colnames(this_gene_aa_modres) <- c("aaindex", "subtype")
      this_gene_aa_modres$uptype <- "Modified Residues"
      this_gene_aa_modres_len <- nrow(this_gene_aa_modres)
      
      for(i in 1:this_gene_aa_modres_len){
        if((up6_loli_modres_info$`Population Mutation`[i] != "none") && (up6_loli_modres_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_modres$vartype[i] <- "pop"
        }else if((up6_loli_modres_info$`Population Mutation`[i] == "none") && (up6_loli_modres_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_modres$vartype[i] <- "path"
        }else if((up6_loli_modres_info$`Population Mutation`[i] != "none") && (up6_loli_modres_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_modres$vartype[i] <- "both"
        }else{
          this_gene_aa_modres$vartype[i] <- "none"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_modres)
      
    }
    
    #molecular processing
    available_molproc <- length(unique(up6_loli_gene_info$`Molecular Processing`))
    if(available_molproc > 1){
      up6_loli_molproc_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Molecular Processing` != "-"))
      this_gene_aa_molproc <- as.data.frame(up6_loli_molproc_info$`Amino Acid Index`)
      this_gene_molproc <- as.data.frame(up6_loli_molproc_info$`Molecular Processing`)
      this_gene_aa_molproc <- cbind(this_gene_aa_molproc, this_gene_molproc)
      colnames(this_gene_aa_molproc) <- c("aaindex", "subtype")
      this_gene_aa_molproc$uptype <- "Molecular Processing Associated Region"
      this_gene_aa_molproc_len <- nrow(this_gene_aa_molproc)
      
      for(i in 1:this_gene_aa_molproc_len){
        if((up6_loli_molproc_info$`Population Mutation`[i] != "none") && (up6_loli_molproc_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_molproc$vartype[i] <- "pop"
        }else if((up6_loli_molproc_info$`Population Mutation`[i] == "none") && (up6_loli_molproc_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_molproc$vartype[i] <- "path"
        }else if((up6_loli_molproc_info$`Population Mutation`[i] != "none") && (up6_loli_molproc_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_molproc$vartype[i] <- "both"
        }else{
          this_gene_aa_molproc$vartype[i] <- "none"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_molproc)
      
    }
    
    droplevels(this_gene_aa_up6_all$subtype)
    
    this_gene_aa_up6_all_fs = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Functional Site", ]
    this_gene_aa_up6_all_fbr = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Functional/Binding Region", ]
    this_gene_aa_up6_all_smr = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Sequence Motif/Region", ]
    this_gene_aa_up6_all_md = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Modular Domain", ]
    this_gene_aa_up6_all_mr = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Modified Residue", ]
    this_gene_aa_up6_all_molp = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Molecular Processing Associated Region", ]
    
    if(nrow(this_gene_aa_up6_all_fs) == 0 
       && nrow(this_gene_aa_up6_all_fbr) == 0 
       && nrow(this_gene_aa_up6_all_smr) == 0 
       && nrow(this_gene_aa_up6_all_md) == 0 
       && nrow(this_gene_aa_up6_all_mr) == 0 
       && nrow(this_gene_aa_up6_all_molp) == 0)
    {
      shinyjs::hide("uniprottext")
    }
    else
    {
      shinyjs::show("uniprottext")
    }
    
    if(nrow(this_gene_aa_up6_all_fs) != 0)
    {
      shinyjs::show("func_site")
    }
    else
    {
      shinyjs::hide("func_site")
    }
    
    if(nrow(this_gene_aa_up6_all_molp) != 0)
    {
      shinyjs::show("mol_proc")
    }
    else
    {
      shinyjs::hide("mol_proc")
    }
    
    if(nrow(this_gene_aa_up6_all_fbr) != 0)
    {
      shinyjs::show("func_bind")
    }
    else
    {
      shinyjs::hide("func_bind")
    }
    
    if(nrow(this_gene_aa_up6_all_smr) != 0)
    {
      shinyjs::show("seq_mot")
    }
    else
    {
      shinyjs::hide("seq_mot")
    }
    
    if(nrow(this_gene_aa_up6_all_md) != 0)
    {
      shinyjs::show("mod_dom")
    }
    else
    {
      shinyjs::hide("mod_dom")
    }
    
    if(nrow(this_gene_aa_up6_all_mr) != 0)
    {
      shinyjs::show("mod_res")
    }
    else
    {
      shinyjs::hide("mod_res")
    }}
    
    
  })
  
  report_info <- eventReactive(input$reportTsubmit, {
    if(input$reportTgeneSelected != ""){
      print(input$reportTgeneSelected)
      gene_protein_trans_class <- subset(GPTLPC_list, geneName == input$reportTgeneSelected)
      uniprotIDreport = gene_protein_trans_class$uniprotAc
      transcriptIDreport = gene_protein_trans_class$transcript
      pclassIDreport = gene_protein_trans_class$pclass
      lenIDreport = gene_protein_trans_class$length
      
      pclassIds = unlist(strsplit(pclassIDreport, ", "))
      pclassIds = c("c0", pclassIds)
      
      pclassNames_andF <- subset(PC24_PF40_pathogenic, pcindex %in% pclassIds)
      pclassNAMEreport = paste0(pclassNames_andF$pcname, collapse=", ")
      
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      report_gene_name = paste("gene_wise_info/",input$reportTgeneSelected,".txt",sep='')
      report_gene_wise_info <- read_delim(report_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      #write.table(report_gene_wise_info, file,sep = sep,row.names = FALSE)
      gene_length <- nrow(report_gene_wise_info)
      updated_choice <- c(1:gene_length) 
      
      updateSelectizeInput(session, "aa_Selected",
                           label = "Amino Acid Index",
                           choices = updated_choice,
                           #selected = NULL,
                           options = list(
                             placeholder = 'Amino Acid Index',
                             onInitialize = I('function() { this.setValue(""); }')
                           )
                           
      )
      
      gene_card_link = c("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", input$reportTgeneSelected)
      gene_card_link_str = paste(gene_card_link, collapse = "")
      uniprotac_link = c("https://www.uniprot.org/uniprot/", uniprotIDreport)
      uniprotac_link_str = paste(uniprotac_link, collapse = "")
      
      #cat(gene_card_link_str)
      HTML(
        paste(
          "<b>","<font size = +1 color = Black face = Arial>","<strong>Gene:</strong> ","</font>","</b><em><a href=", gene_card_link_str, " target=_blank>", input$reportTgeneSelected,"</a></em><br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>UniProt-AC:</strong> ","</font>","</b><a href=", uniprotac_link_str, " target=_blank>",uniprotIDreport,"</a><br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>Protein length:</strong> ","</font>","</b>",lenIDreport," amino acids<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>Transcript:</strong> ","</font>","</b>",transcriptIDreport,"<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>Annotated Protein class(es) for ","<i>", input$reportTgeneSelected, "</i>", ": </strong> ","</font>","</b>",pclassNAMEreport,"<br>", "<br></br>"
        )
      )
    }
  })
  
  ## aa wise feature page
  observeEvent(
    input$submit,{
      if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
        shinyjs::show("aaFeature_wellPanel")
        #shinyjs::show("report_FEATwellpanel")
      }
    }
  )
  
  uniprot_source_link = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/"
  
  #documentation
  output$doc = renderText(
    {
      HTML(
        paste(
          "<p><span style='font-size: 20pt;'>Tutorials and Resources</span></p>", "<hr>",
          "<p><span style='font-size: 16pt; font-family: arial, helvetica, sans-serif;color: windowtext;'>Check out tutorials and application showcases of MISCAST in <a href='https://www.youtube.com/playlist?list=PL0U7anUEN5ksUGn8sWfipgRz3FIFdMcBL' target='_blank'><strong><span style='font-size: 20pt; color: #ff0000;'>YouTube</span></strong></a></span></p>",
          "<br><br>",
          "<p><span style='font-size: 20pt;'>Protein Feature Mining and Annotation</span></p>", "<hr>",
          #"<p><span style='font-size: 16pt; font-family: arial, helvetica, sans-serif;'><strong>Protein Structural, Physicochemical, and Functional Feature Mining and Annotation:</strong></span></p>",
          "<p><span style='font-size: 16pt; font-family: arial, helvetica, sans-serif;color: windowtext;'>We annotated amino acid residues with&nbsp;forty&nbsp;different&nbsp;structural, physicochemical, functional features, grouped&nbsp;into seven main actegories:</span></p>",
          "<ol>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>3-class secondary structure (feature count: 3)</span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>8-class secondary structure (feature count: 8)</span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>Residue exposure level (feature count: 5)</span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>Physicochemical properties of amino acids (feature count: 8)</span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>Protein-protein interaction types (feature count: 4)</span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>Post-translational modification types (feature count: 6)</span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>UniProt-based functional features (feature count: 6)</span></li>",
          "</ol>","<br>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 14pt;'>3-class secondary structure (feature count: 3)</span></strong></span></p>",
          "<p style='text-align: justify;'><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='color: windowtext;'>Protein secondary structure is the three-dimensional (3D) form of local segments of proteins that spontaneously form as an intermediate element before the protein folds into its three-dimensional tertiary structure. For proteins with multiple solved structures available in the PDB, we came up with a consensus annotation representing the maximum agreement of structure type assignments from all the available structures.</span></span></p>", 
          #Secondary structure is defined by the pattern of hydrogen bonds between the side chain atoms of amino acids in the peptide backbone.
          "<p style='text-align: justify;'><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif; color: windowtext;'><strong>Data source:</strong> <a href='https://swift.cmbi.umcn.nl/gv/dssp/' target='_blank'>DSSP</a>, a database of secondary structure assignments for all entries in the Protein Data Bank (PDB). Based on the hydrogen bond energy, DSSP asigns eight different types of secondary structures to the amino acid residues in protein structures, which can be broadly grouped into three categories.</span></p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 12pt; color: windowtext;'>Features:</span></strong></span></p>",
          "<ol>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>&beta;-strand/sheet</strong>: &beta;-strand or &beta;-sheet</span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Helices</strong>: 3<sub>10</sub>-helix, &alpha;-helix, &pi;-helix</span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Coils</strong>: random loop, bend or turn</span></li>",
          "</ol>",
          "<p>&nbsp;</p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 14pt;'>8-class secondary structure (feature count: 8)</span></strong></span></p>",
          #"<p style='text-align: justify;'><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='color: windowtext;'>Protein secondary structure is the three-dimensional (3D) form of local segments of proteins that spontaneously form as an intermediate element before the protein folds into its three-dimensional tertiary structure. Secondary structure is defined by the pattern of hydrogen bonds between the side chain atoms of amino acids in the peptide backbone.</span></span></p>",
          "<p><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif; color: windowtext;'><strong>Data source:</strong> <a href='https://swift.cmbi.umcn.nl/gv/dssp/' target='_blank'>DSSP</a> (see above).</span></p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;color: windowtext;'><strong><span style='font-size: 12pt; color: windowtext;'>Features:</span></strong></span></p>",
          "<ol>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>&beta;-strand</strong> <span style='color: windowtext;'>&ndash; residue in isolated </span><span style='color: windowtext;'>&beta;-bridge</span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>&beta;-sheet</strong> <span style='color: windowtext;'>&ndash; residue in the extended strand that participates in &beta;-</span><span style='color: windowtext;'>ladder</span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>3<sub>10</sub>-helix</strong> <span style='color: windowtext;'>&ndash; hydrogen bond between <em>i</em><sup>th</sup> and (<em>i</em>+3)<sup>th</sup> residue to build each helical turn</span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>&alpha;-helix</strong> <span style='color: windowtext;'>&ndash; hydrogen bond between <em>i</em><sup>th</sup> and (<em>i</em>+4)<sup>th</sup> residue to build each helical turn</span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>&pi;-helix</strong> <span style='color: windowtext;'>&ndash; hydrogen bond between <em>i</em><sup>th</sup> and (<em>i</em>+5)<sup>th</sup> residue to build each helical turn</span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Bend</strong> <span style='color: windowtext;'>&ndash; residues of high curvature where the angle between C<sub>i</sub>C<sub>i+2</sub> and C<sub>i-2</sub>C<sub>i</sub> is at least 70</span><span style='color: windowtext;'>&deg;</span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Turn</strong>&nbsp;<span style='color: windowtext;'>&ndash; hydrogen bonded turn</span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Loop</strong> <span style='color: windowtext;'>&ndash; random loop where no other rule applies</span></span></li>",
          "</ol>",
          "<p>&nbsp;</p>",    
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 14pt;'>Residue exposure level (feature count: 5)</span></strong></span></p>",
          "<p style='text-align: justify;'><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='color: windowtext;'>The solvent-exposure level of an amino acid residue is defined by its relative accessible surface area (RSA) in the protein's tertiary structure. The accessible surface area (ASA) of a biomolecule&nbsp;is its surface area that is accessible to a solvent which is measured in the unit of Angstroms squared (Å<sup>2</sup>). Here based on the RSA, each amino acid is assigned with one out of five different exposure levels.</span></span></p>",
          "<p style='text-align: justify;'><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;color: windowtext;'><strong>Data source:</strong> <a href='https://swift.cmbi.umcn.nl/gv/dssp/' target='_blank'>DSSP</a>. Given the 3D atomic coordinates of a protein structure, the DSSP program calculates that residue&rsquo;s water-exposed surface area, referred to as the accessible surface area (ASA). <span style='font-size: 12pt; color: windowtext;'>The Relative ASA or RSA of a protein residue is calculated by normalizing the ASA of that residue by the surface area of the same type of residue in a reference state. We used the ASA normalizing values derived by </span><span style='font-size: 12pt; color: windowtext;'>[<a href='https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0080635' target='_blank'>Tien et.al, 2013</a>]</span><span style='font-size: 12pt; color: windowtext;'> using the Gly-X-Gly tripeptide as the reference state for a given residue X. For proteins with multiple solved structures available in the PDB, we came up with a consensus RSA per amino acid which represents the median RSA out of all available structures.</span></span></p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 12pt;color: windowtext;'>Features:</span></strong></span></p>",
          "<ol>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Core</strong>: RSA &lt; 5%</span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Buried</strong>: 5% &lt;= RSA &lt; 25%</span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Medium-buried</strong>: 2<span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>5% &lt;= RSA &lt; 50%</span></span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Medium-exposed</strong>: <span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>50% &lt;= RSA &lt; 75%</span></span></span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Exposed</strong>: <span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>RSA &gt; 75%</span></span></span></li>",
          "</ol>",
          "<p>&nbsp;</p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 14pt;'>Physicochemical properties of amino acid (feature count: 8)</span></strong></span></p>",
          "<p style='text-align: justify;'><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='color: windowtext;'><span style='color: windowtext;'><span style='color: windowtext; background: white;'>Amino acids are the building blocks of peptides and proteins. While they all have common elements of an amine group, a carboxyl group and a side chain, different functional groups that comprise the side chain give each amino acid distinct physicochemical properties that influence protein conformation and function.</span></span></span></span></p>",
          "<p style='text-align: justify;'><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif; color: windowtext;'><strong>Data source:</strong> <span style='font-size: 12pt;'><span style='color: windowtext;'><span style='font-size: 12pt; color: windowtext; background: white;'>&nbsp;Based on the property of the side chain, the 20 natural amino acids are grouped into eight categories</span></span></span><span style='font-size: 12pt; color: windowtext;'>.&nbsp;</span></span></p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 12pt; color: windowtext;'>Features:</span></strong></span></p>",
          "<ol>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Aliphatic</strong>: <span style='color: windowtext; background: white;'>Alanine (Ala/A), Isoleucine (Ile/I), Leucine (Leu/L), Methionine (Met/M), Valine (Val/V)</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><strong>Aromatic</strong>: <span style='color: windowtext; background: white;'>Phenylalanine (Phe/F), Tryptophan (Trp/W), Tyrosine (Tyr/Y)</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><strong>Hydrophobic</strong>:&nbsp;aliphatic and aromatic amino acids</span></li>",
          "<li><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif;'><strong>Positively-charged</strong>: <span style='color: windowtext; background: white;'>Arginine (Arg/R), Histidine (His/H), Lysine (Lys/K)</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='font-family: arial, helvetica, sans-serif;'><strong>Negatively-charged</strong>: <span style='color: windowtext; background: white;'>Aspartic acid (Asp/D), Glutamic acid (Glu/E)</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='font-family: arial, helvetica, sans-serif;'><strong>Neutral</strong>: <span style='color: windowtext; background: white;'>Asparagine (Asn/N), Glutamine (Gln/Q), Serine (Ser/S), Threonine (Thr/T)</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='font-family: arial, helvetica, sans-serif;'><strong>Polar</strong>: positively-charged, negatively-charged, and neutral amino acids</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='font-family: arial, helvetica, sans-serif;'><strong>Special</strong>: </span></span>",
          "<ul>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='font-family: arial, helvetica, sans-serif;'><span style='color: windowtext; background: white;'>Proline (Pro/P): doesn't have hydrogen for backbone hydrogen bond formation</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='font-family: arial, helvetica, sans-serif;'><span style='color: windowtext; background: white;'>Glycine (Gly/G): the smallest amino acid, which is very flexible</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='font-family: arial, helvetica, sans-serif;'><span style='color: windowtext; background: white;'>Cystine (Cys/C): has a very reactive sulfhydryl group</span></span></span></li>",
          "</ul>",
          "</li>",
          "</ol>",
          "<p>&nbsp;</p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 14pt;'>Protein-protein interactions (feature count: 4)</span></strong></span></p>",
          "<p style='text-align: justify;'><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='color: windowtext;'><span style='color: windowtext;'><span style='color: windowtext; background: white;'>The amino acid residues are annotated with interaction types between amino acid pairs at distance &le; 4&nbsp;Å on the protein-protein interface.&nbsp;</span></span></span></span></p>",
          "<p style='text-align: justify;'><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif; color: windowtext;'><strong>Data source:</strong> <span style='color: windowtext;'><span style='color: windowtext; background: white;'> <span style='color: windowtext; background: white;'>T<span style='font-family: arial, helvetica, sans-serif;'>he&nbsp;<a href='http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=index.html' target='_blank'>PDBsum</a> database. <span style='color: windowtext; background: white;'></span><span style='color: windowtext;'>The &ldquo;prot<em>-</em>prot&rdquo; entries of the PDBsum database </span><span style='color: windowtext; background: white;'>was parsed to collect the residue-residue interaction annotations. If there are multiple solved structures available of a gene in multiple complexes, one amino acid can possibly be assigned to have multiple bonds and, we collected records of all types of bonds in which one residue was involved in different conformations.</span></span></span></span></span></span></p>",
          "<p><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><strong><span style='color: windowtext;'>Features:</span></strong></span></p>",
          "<ol>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='font-family: arial, helvetica, sans-serif;'><strong>Disulfide bond</strong>:<span style='color: windowtext; background: white;'> a covalent bond between cysteine side chains, yielding a disulfide bridge, which is an key determinant of protein 3D structure</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><strong>Salt-bridge ionic interaction</strong>: an&nbsp;<span style='color: windowtext; background: white;'>an ionic bond between oppositely charged residues. </span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><strong>Hydrogen bond</strong>: <span style='font-size: 12pt; color: windowtext; background: white;'>electrostatic interaction between two atoms bearing partial negative charges, that share a partially positively charged hydrogen.</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><strong>Nonbonded Van der&nbsp;Waals interaction</strong>: <span style='font-size: 12pt; color: windowtext; background: white;'>a bond formed by the transient and weak Van der Waals attraction between two close atoms. </span></span></li>",
          "</ol>",
          "<p>&nbsp;</p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 14pt;'>Post-translational modifications (feature count: 6)</span></strong></span></p>",
          "<p style='text-align: justify;'><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='color: windowtext;'><span style='color: windowtext; background: white;'>Post-translational modification (PTM) refers to the covalent and generally enzymatically-mediated modification of proteins to form the mature protein.&nbsp; </span></span></span></p>",
          "<p style='text-align: justify;'><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif; color: windowtext;'><strong>Data source:</strong> <span style='font-size: 12pt;'><span style='color: windowtext;'><span style='font-size: 12pt; color: windowtext; background: white;'> The <a href='https://www.phosphosite.org/homeAction.action' target='_blank'>PhosphoSitePlus</a> database.<span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 12pt; color: windowtext;'> We collected six diffrenet types of PTM annotations from the database and mapped them onto the protein structures. We then computed the spatial distances between amino acids on the structure and different PTM sites. We annotated each amino acid with the nearest distance to six different PTM sites. </span></span></span></span></span></span></span></p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 12pt; color: windowtext;'>Features:</span></strong></span></p>",
          "<ol>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Acetylation</strong>: introduces an acetyl group.</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Methylation</strong>: addition of a methyl group.&nbsp;</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><strong>O.GlcNAc</strong>: also known as &beta;-linked N-acetylglucosamine, is a form of protein glycosylation.&nbsp;</span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Phosphorylation</strong>: attachment of a phosphoryl group.&nbsp;</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Sumoylation</strong>: addition of SUMOs (small ubiquitin-like modifiers).</span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Ubiquitination</strong>: attachment of a ubiquitin.&nbsp;</span></span></li>",
          "</ol>",
          "<p>&nbsp;</p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 14pt;'>UniProt-based functional features (feature count: 6)</span></strong></span></p>",
          "<p style='text-align: justify;'><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'><span style='color: windowtext;'><span style='color: windowtext; background: white;'>The amino acids are annotated with 25 different function<span style='font-family: arial, helvetica, sans-serif;'>al features, collected from <a href='ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_beds/' target='_blank'>UniProt genome annotation tracks</a>. Guided by the <a href='https://www.ebi.ac.uk/training/online/course/UniProt-exploring-protein-sequence-and-functional-information/sequence-features' target='_blank'>EMBL-EBI training material</a> and similarity, these features are grouped into six categories. <span style='font-size: 12pt; color: windowtext;'>If an amino acid residue has any of the sub-feature of a category, that residue is annotated by the corresponding feature</span></span></span></span></p>",
          "<p style='text-align: justify;'><span style='font-size: 12pt; font-family: arial, helvetica, sans-serif; color: windowtext;'><strong>Data source:</strong><a href='https://www.uniprot.org/' target='_blank'>UniProt</a> database</p>",
          "<p><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 12pt; color: windowtext;'>Features:</span></strong></span></p>",
          "<ol>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Functional site</strong></span></span>",
          "<ul>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Active site <span style='font-size: 12pt; color: windowtext;'>&ndash; <span style='font-family: arial, helvetica, sans-serif;'>amino acids that are directly involved in the activity of an enzyme</span></span></span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Metal binding site <span style='font-size: 12pt; color: windowtext;'>&ndash; binding site for a metal ion</span></span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Binding site <span style='font-size: 12pt; color: windowtext;'>&ndash; binding site for any chemical groups (co-enzyme, prosthetic group, etc.)</span> </span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Site <span style='font-size: 12pt; color: windowtext;'>&ndash; any interesting amino single amino acid site</span></span></span></span></li>",
          "</ul>",
          "</li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Functional/binding region</strong></span></span>",
          "<ul>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Zinc finger <span style='font-size: 12pt; color: windowtext;'>&ndash; position of zinc fingers within the protein</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>DNA binding region <span style='font-size: 12pt; color: windowtext;'>&ndash; position of DNA binding region</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Nucleotide binding region <span style='font-size: 12pt; color: windowtext;'>&ndash; position of nucleotide phosphate binding region</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Calcium-binding region <span style='font-size: 12pt; color: windowtext;'>&ndash; position of calcium-binding region</span></span></span></li>",
          "</ul>",
          "</li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 12pt;'>Sequence motif/region</span></strong></span>",
          "<ul>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>Region<span style='font-size: 16px;'> <span style='font-size: 12pt; color: windowtext;'>&ndash; region of interest in the sequence</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>Repeat<span style='font-size: 16px;'> <span style='font-size: 12pt; color: windowtext;'>&ndash; repeated sequence motifs</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>Coiled-coil<span style='font-size: 16px;'> <span style='font-size: 12pt; color: windowtext;'>&ndash; regions of coiled-coil within a protein</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif; font-size: 12pt;'>Motif<span style='font-size: 16px;'> <span style='font-size: 12pt; color: windowtext;'>&ndash; short (around 20 amino acid) sequence motif of biological interest</span></span></span></li>",
          "</ul>",
          "</li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><strong><span style='font-size: 16px;'>Modular domain</span></strong></span>",
          "<ul>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Domain <span style='font-size: 12pt; color: windowtext;'>&ndash; </span><span style='font-size: 12pt; color: windowtext;'>positions &nbsp;of modular domains in protein</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Topological domain <span style='font-size: 12pt; color: windowtext;'>&ndash; location of non-membrane regions of membrane-spanning proteins</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Transmembrane <span style='font-size: 12pt; color: windowtext;'>&ndash; a membrane-spanning region</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Intramembrane <span style='font-size: 12pt; color: windowtext;'>&ndash; region located in a membrane without crossing it</span></span></span></li>",
          "</ul>",
          "</li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Molecular-processing-associated region</strong>:</span></span>",
          "<ul>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Peptide <span style='font-size: 12pt; color: windowtext;'>&ndash; extent of an active peptide in the mature protein</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Transit peptide <span style='font-size: 12pt; color: windowtext;'>&ndash; extent of a transit peptide for organelle targeting</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Signal peptide <span style='font-size: 12pt; color: windowtext;'>&ndash; sequence targeting proteins to the secretory pathway or periplasmic space</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>Propeptide <span style='font-size: 12pt; color: windowtext;'>&ndash; part of a protein that is cleaved during maturation or activation</span></span></span></li>",
          "</ul>",
          "</li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'><strong>Modified residues</strong>:</span></span>",
          "<ul>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>modified residue <span style='font-size: 12.0pt; font-family: 'Times New Roman',serif; color: windowtext;'>&ndash; modified residues by posttranslational modification excluding lipids, glycans and protein crosslinks</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>lipidation <span style='font-size: 12.0pt; font-family: 'Times New Roman',serif; color: windowtext;'>&ndash; covalently attached lipid groups</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>disulfide bond <span style='font-size: 12.0pt; font-family: 'Times New Roman',serif; color: windowtext;'>&ndash; cysteine residues participating in disulfide bonds</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>cross-link <span style='font-size: 12.0pt; font-family: 'Times New Roman',serif; color: windowtext;'>&ndash; residues participating in covalent linkages between proteins</span></span></span></li>",
          "<li><span style='font-family: arial, helvetica, sans-serif;'><span style='font-size: 16px;'>glycosylation <span style='font-size: 12.0pt; font-family: 'Times New Roman',serif; color: windowtext;'>&ndash; covalently attached glycan group.</span></span></span></li>",
          "</ul>",
          "</li>",
          "</ol>"
        )
      )
    }
  )
  
  #documentation
  output$about_page = renderText(
    {
      gnomAD_link = paste(c("http://gnomad.broadinstitute.org/", ""), collapse = "")
      ClinVar_link = paste(c("https://www.ncbi.nlm.nih.gov/clinvar/", ""), collapse = "")
      HGMD_link = paste(c("http://www.hgmd.cf.ac.uk/ac/index.php", ""), collapse = "")
      SC_link = paste(c("https://www.broadinstitute.org/stanley", ""), collapse = "")
      Broad_link = paste(c("https://www.broadinstitute.org/", ""), collapse = "")
      PDB_link = paste(c("https://www.rcsb.org/", ""), collapse = "")
      DSSP_link = paste(c("https://swift.cmbi.umcn.nl/gv/dssp/", ""), collapse = "")
      PDBsum_link = paste(c("http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode=index.html", ""), collapse = "")
      SIFTS_link = paste(c("https://www.ebi.ac.uk/pdbe/docs/sifts/", ""), collapse = "")
      PANTHER_link = paste(c("http://www.pantherdb.org/panther/ontologies.jsp?", ""), collapse = "")
      PPS_link = paste(c("https://www.phosphosite.org/homeAction", ""), collapse = "")
      Uniprot_link = paste(c("https://www.uniprot.org/", ""), collapse = "")
      
      HTML(
        paste(
          "<p><span style='font-size: 24pt;'>About MISCAST</span></p>", "<hr>",
          "<color = Black face = Arial align = justify>", "<div align='justify'>",
          "<b>MISCAST</b> is developed at the ",
          "<a href=", Broad_link, " target=_blank>","Broad Institute","</a>","of ",
          "MIT and Harvard, by a combined effort from the Genetics and Therapeutics groups of",
          "<a href=", SC_link, " target=_blank>","Stanley Center</a>.",
          "The goal of MISCAST is to visualize and analyze single amino-acid-altering missense variants on protein sequence and 3-dimensional structure, and thereby forecast their biological impact.",
          "The project combines structural information from the largest database of protein 3D structures, the ",
          "<a href=", PDB_link, " target=_blank>","the protein data bank","</a>"," (PDB), ",
          "with genetic variation data from freely accessible large-scale databases such as the ",
          "<a href=", gnomAD_link, " target=_blank>","genome aggregation database","</a>"," (gnomAD), ",
          "<a href=", ClinVar_link, " target=_blank>","ClinVar","</a>","and the ",
          "<a href=", HGMD_link, " target=_blank>","human gene mutation database","</a>"," (HGMD).",
          "<br><br>",
          "The dataset in this web server spans 1,330 human genes, with 406,449 population variants from gnomAD (release 2.1.1) and ",
          "54,137 variants from the ClinVar (February, 2019 release, pathogenic and likely-pathogenic) and HGMD (Professional release 2018.4, disease mutations) databases, ", "and >17k protein 3D structures solved in human. ",
          "Using EMBL-EBI supported ",
          "<a href=", SIFTS_link, " target=_blank>","SIFTS","</a>",
          " resource, the missense variants were mapped onto protein 3D structures using an automated pipeline. ",
          "The amino acid residues of 1,330 genes were annotated with seven structural, physicochemical, and functional",
          "features comprised of 40 subtypes, collected from multiple resources spanning ",
          "<a href=", DSSP_link, "target=_blank>","DSSP</a>,",
          "<a href=", PDBsum_link, "target=_blank>","PDBsum</a>,",
          "<a href=", PPS_link, " target=_blank>","PhosphoSitePlus</a>,",
          "<a href=", PANTHER_link, " target=_blank>","PANTHER</a>,",
          "<a href=", Uniprot_link, " target=_blank>","UniProt</a>.",
          "For the details about feature set annotation and mining, we refer the user to read the documentation page. ",
          "<br><br>",
          "The website provides two separate tracks: <strong>Variant Summary Report</strong> and <strong>Variant Analysis Suite</strong>.",
          "<br><br>",
          "The <strong>Variant Summary Report</strong> track provides an amino acid-wise report that summarizes its features, and compares with",
          "precomputed pathogenic  and population variant-associated features (Ref:<a href='https://www.biorxiv.org/content/10.1101/693259v1' target=_blank>","Paper</a>). ",
          #"<br><br>",
          "Here, the impact of an amino acid substitution is estimated by the difference of the number of pathogenic  and population variant-associated features of the reference amino acid, termed as Pathogenic 3D Feature index (P3DFi). The amino acid residues located in vulnerable 3D sites will have a higher number of pathogenic variant-associated features (P3DFi > 0) whereas amino acids substituted in benign variants are expected to have a greater number of population variant-associated features (P3DFi < 0). The P3DFi score per amino acid of a protein can be computed taking the <i>in general</i> characteristic 3D features associated to pathogenic and population variants of all proteins as well the features that are found specifically critical for the function of the protein (as defined by the protein class). ",
          "<br><br>",
          "The <strong>Variant Analysis Suite</strong> track provides a comprehensive visualization and analysis platform to inspect missense variants in protein sequence and structure space alongside a rich collection of protein features (described in the Documentation page). Besides interactive visualization, all the gene-wise annotation tracks are downloadable as text file, and the annotated protein structures are exportable in PyMOL for downstream analysis.",
          "<br><br><br><br>",
          "<p><span style='font-size: 20pt;'>Contact and Citation information</span></p>", "<hr>",
          "For bug reports, comments, questions or suggestions, please contact", "<a href='mailto:sumaiya@broadinstitute.org?cc=dlal@broadinstitute.org,arthurc@broadinstitute.org'>here</a>.<br><br>",
          "We request that any use of data obtained from the MISCAST cite our bioRxiv submission:<br></br> Sumaiya Iqbal, Jakob B. Jespersen, Patrick May, Eduardo Perez-Palma, David Hoksza, Henrike O. Heyne, Shehab S. Ahmed, Zaara T. Rifat, M. Sohel Rahman, Kasper Lage, Aarno Palotie, Jeffrey R. Cottrell, Florence F. Wagner, Mark J. Daly, Arthur J. Campbell, Dennis Lal. <i>Insights into Protein Structural, Physicochemical, and Functional Consequences of Missense Variants in 1,330 Disease-associated Human Genes</i>. <strong>bioRxiv</strong> 693259 (2019); DOI: <a href=https://www.biorxiv.org/content/10.1101/693259v1 target=_blank> https://doi.org/10.1101/69325 </a>",
          #"There is no need to include us as authors on your manuscript, unless we contributed specific advice or analysis for your work.",
          "<br><br><br><br>",
          "<p><span style='font-size: 20pt;'>Data usuage and Disclaimer</span></p>", "<hr>",
          "For the benefit of the wider biomedical community, all data here are open. Users can freely download and search the data, and we encourage the use and publication of results generated from these data. ", 
          "There are absolutely no restrictions or embargoes on the publication of results derived from MISCAST data. However, MISCAST aggregates data from many different other resources, thus, there might be inconsistencies around versions. ",
          "Please contact us with any concerns, we will do our best to address them. ",
          "<br><br>",
          "<strong><i>Disclaimer</i></strong>: Please note that all content on this website is provided for academic and research purposes only. It has not undergone clinical validation, and should not be used for medical advice, diagnosis, or treatment.",
          "<br><br><br><br>",
          "<p><span style='font-size: 20pt;'>Development Team</span></p>", "<hr>",
          "<table style='width: 100%; border-collapse: collapse; border-style: none;' border='0'>",
          "<tbody>",
          "<tr>",
          "<td style='width: 32.9673%; vertical-align: top;'>",
          "<p><span style='font-size: 14pt; color: #000000;'><strong>Principal Investigators</strong></span></p>",
          "<p><span style='color: #000000;'>Dennis Lal</span></p>",
          "<p><span style='color: #000000;'>Mark Daly</span></p>",
          "<p><span style='color: #000000;'>Florence F. Wagner</span></p>",
          "<p><span style='color: #000000;'>Arthur J. Campbell</span></p>",
          "<p><span style='color: #000000;'>Jeffrey R. Cottrell</span></p>",
          "</td>",
          "<td style='width: 32.9673%; vertical-align: top;'>",
          "<p><span style='font-size: 14pt; color: #000000;'><strong>Analysis Team</strong></span></p>",
          "<p><span style='color: #000000;'>Sumaiya Iqbal</span></p>",
          "<p><span style='color: #000000;'>Patrick May</span></p>",
          "<p><span style='color: #000000;'>Jakob Berg Jespersen</span></p>",
          "<p><span style='color: #000000;'>Eduardo Perez-Palma</span></p>",
          "<p><span style='color: #000000;'>Henrike Heyne</span></p>",
          "<p>&nbsp;</p>",
          "</td>",
          "<td style='width: 32.9673%; vertical-align: top;'>",
          "<p><span style='font-size: 14pt; color: #000000;'><strong>Website Team</strong></span></p>",
          "<p><span style='color: #000000;'>Sumaiya Iqbal</span></p>",
          "<p><span style='color: #000000;'>David Hoksza</span></p>",
          "<p><span style='color: #000000;'>Shehab Sarar Ahmed</span></p>",
          "<p><span style='color: #000000;'>Zaara Tasnim Rifat</span></p>",
          "<p><span style='color: #000000;'>Eduardo Perez-Palma</span></p>",
          "</td>",
          "</tr>",
          "</tbody>",
          "</table>",
          "<br>",
          "<hr>",
          "</div>"
        )
      )
    }
  )
  
  report_info_aa_wise <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      gene_protein_trans_class <- subset(GPTLPC_list, geneName == input$geneSelected)
      uniprotIDreport = gene_protein_trans_class$uniprotAc
      transcriptIDreport = gene_protein_trans_class$transcript
      pclassIDreport = gene_protein_trans_class$pclass
      lenIDreport = gene_protein_trans_class$length
      
      pclassIds = unlist(strsplit(pclassIDreport, ", "))
      pclassIds = c("c0", pclassIds)
      pclassNames_andF <- subset(PC24_PF40_pathogenic, pcindex %in% pclassIds)
      pclassNAMEreport = paste0(pclassNames_andF$pcname, collapse=", ")
      
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      report_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      report_gene_wise_info <- read_delim(report_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      #write.table(report_gene_wise_info, file,sep = sep,row.names = FALSE)
      gene_length <- nrow(report_gene_wise_info)
      updated_choice <- c(1:gene_length) 
      
      updateSelectizeInput(session, "aa_Selected_aaFeat",
                           label = "Amino Acid Index",
                           choices = updated_choice,
                           #selected = NULL,
                           options = list(
                             placeholder = 'Amino Acid Index',
                             onInitialize = I('function() { this.setValue(""); }')
                           )
                           
      )
      
      gene_card_link = c("https://www.genecards.org/cgi-bin/carddisp.pl?gene=", input$geneSelected)
      gene_card_link_str = paste(gene_card_link, collapse = "")
      uniprotac_link = c("https://www.uniprot.org/uniprot/", uniprotIDreport)
      uniprotac_link_str = paste(uniprotac_link, collapse = "")
      
      #cat(gene_card_link_str)
      HTML(
        paste(
          #"<b>","<font size = +1 color = Black face = Arial>","<strong>Gene:</strong> ","</font>","</b><em><a href=", gene_card_link_str, " target=_blank>", input$reportTgeneSelected,"</a></em><br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>UniProt-AC:</strong> ","</font>","</b><a href=", uniprotac_link_str, " target=_blank>",uniprotIDreport,"</a><br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>Protein length:</strong> ","</font>","</b>",lenIDreport," amino acids<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>Transcript:</strong> ","</font>","</b>",transcriptIDreport,"<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>Protein class:</strong> ","</font>","</b>",pclassNAMEreport,"<br>", "<br></br>"
        )
      )
    }
  })
  
  output$reportTinformation_aaFeat = renderText(
    {
      report_info_aa_wise()
    }
  )
  
  
  ## AA wise feature
  
  report_aaFeature <- eventReactive({
    input$aa_Selected_aaFeat
    input$submit
  },{
    #reportGeneName <- input$reportTgeneSelected
    if(input$aa_Selected_aaFeat != ""){
      
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      report_gene_nameF = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      report_gene_wise_infoF <- read_delim(report_gene_nameF, "\t", escape_double = FALSE, trim_ws = TRUE)
      
      aa_wise_feature <- subset(report_gene_wise_infoF, `Amino Acid Index` == input$aa_Selected_aaFeat)
      amino_acid_reportF <- aa_wise_feature$`Amino Acid`
      
      thisaa1 <- amino_acid_reportF
      thisaaindex <- input$aa_Selected_aaFeat
      this_aa_info <- subset(aa_info, `one` == thisaa1)
      thisaa3 <- this_aa_info$three
      thisaafull <- this_aa_info$full
      thisaaprp <- this_aa_info$prp
      thisaacomm <- this_aa_info$comm
      
      aaFull_reportF <- c(thisaafull, " (", thisaa3, "/", thisaa1, ") at position ", thisaaindex)
      aaFull_reportF_str <- paste(aaFull_reportF, collapse = "")
      
      aaShorter <- c(thisaa3, "/", thisaa1, thisaaindex)
      aaShorter_str <- paste(aaShorter, collapse = "")
      
      # SS3 and SS8
      structure = 0;
      ss3_reportF <- aa_wise_feature$`3-class DSSP Secondary Structure Properties`
      ss8_reportF <- aa_wise_feature$`8-class DSSP Secondary Structure Properties`
      
      if(ss8_reportF == "-"){
        local_conformSS3 <- "annotation not available"
        local_conformSS3Table <- "none"
        local_conformSS8 <- "annotation not available"
        local_conformSS8Table <- "none"
        
        structure = 0
        aa_SS3_flag <- 0
        aa_SS8_flag <- 0
        }else{
          structure = 1
          
          aa_SS3_flag <- 0
          if(grepl("beta", ss3_reportF)){
            local_conformSS3 <- "&beta;-strand/sheet"
            aa_SS3_flag <- 1
          }else if(grepl("Heli", ss3_reportF)){
            local_conformSS3 <- "helix"
            aa_SS3_flag <- 2
          }else if(grepl("coil", ss3_reportF)){
            local_conformSS3 <- "coil"
            aa_SS3_flag <- 3
          }
          
          aa_SS8_flag <- 0
          if(grepl("strand", ss8_reportF)){
            local_conformSS8 <- "&beta;-strand"
            aa_SS8_flag <- 1
          }else if(grepl("sheet", ss8_reportF)){
            local_conformSS8 <- "&beta;-sheet"
            aa_SS8_flag <- 2
          }else if(grepl("310", ss8_reportF)){
            local_conformSS8 <- "3<sub>10</sub>-helix"
            aa_SS8_flag <- 3
          }else if(grepl("alpha", ss8_reportF)){
            local_conformSS8 <- "&alpha;-helix"
            aa_SS8_flag <- 4
          }else if(grepl("pi", ss8_reportF)){
            local_conformSS8 <- "&pi;-helix"
            aa_SS8_flag <- 5
          }else if(grepl("loop", ss8_reportF)){
            local_conformSS8 <- "coil-loop"
            aa_SS8_flag <- 8
          }else if(grepl("bend", ss8_reportF)){
            local_conformSS8 <- "coil-bend"
            aa_SS8_flag <- 6
          }else if(grepl("turn", ss8_reportF)){
            local_conformSS8 <- "coil-turn"
            aa_SS8_flag <- 7
          }
      }
      
      exposure_reportF <- aa_wise_feature$`Exposure`
      acc_reportF <- aa_wise_feature$`Accessible Surface Area`
      rsa_reportF <- aa_wise_feature$`Relative Accessible Surface Area`
      if(exposure_reportF == "-"){
        solv_exposure <- "annotation not available"
        solv_exposureTable <- "none"
        aa_ASA_flag <- 0
      }else{
        rsa_reportF_num <- as.numeric(rsa_reportF)
        aa_ASA_flag <- 0
        if(rsa_reportF_num < 0.05){
          rsaReport <- "(<5% relative solvent accessibility)"
          aa_ASA_flag <- 1
        }else if((rsa_reportF_num >= 0.05) && (rsa_reportF_num < 0.25)){
          rsaReport <- "(5 - 25% relative solvent accessibility)"
          aa_ASA_flag <- 2
        }else if((rsa_reportF_num >= 0.25) && (rsa_reportF_num < 0.50)){
          rsaReport <- "(25 - 50% relative accessible surface area)"
          aa_ASA_flag <- 3
        }else if((rsa_reportF_num >= 0.50) && (rsa_reportF_num < 0.75)){
          rsaReport <- "(50 - 75% relative solvent accessibility)"
          aa_ASA_flag <- 4
        }else if(rsa_reportF_num > 0.75){
          rsaReport <- "(>75% relative solvent accessibility)"
          aa_ASA_flag <- 5
        }
        solv_exposure <- paste(exposure_reportF, rsaReport, sep = " ")
        solv_exposure <- paste(solv_exposure, acc_reportF, sep = ", solvent accessible surface area: ")
        solv_exposure <- paste(solv_exposure, "&#8491;<sup>2</sup>", sep = " ")
      }
      
      
      chemp_reportF <- aa_wise_feature$`Physiochemical Property`
      if(thisaacomm != "na"){
        chem_reportF_list <- c(thisaaprp, " ", thisaacomm)
      }else{
        chem_reportF_list <- c(thisaaprp)
      }
      
      aa_CHEM_flag <- c()
      if(thisaaprp == 'hydrophobic & aliphatic'){
        aa_CHEM_flag <- c(1, 2)  
      }else if(thisaaprp == 'hydrophobic & aromatic'){
        aa_CHEM_flag <- c(1, 3)  
      }else if(thisaaprp == 'polar & positively-charged'){
        aa_CHEM_flag <- c(4, 5)  
      }else if(thisaaprp == 'polar & negatively-charged'){
        aa_CHEM_flag <- c(4, 6)  
      }else if(thisaaprp == 'polar & neutral'){
        aa_CHEM_flag <- c(4, 7)  
      }else if(thisaaprp == 'special'){
        aa_CHEM_flag <- c(8)  
      }
      
      
      chem_reportF_str <- paste(chem_reportF_list, collapse="")
      chem_reportF_strTable <- paste(chem_reportF_list, collapse="")
      
      aa_BOND_flag <- c()
      bond_reportF <- aa_wise_feature$`Bond Detail`
      if(bond_reportF == "-"){
        bondinfo_reportF <- "annotation not available"
        bondinfo_reportFTable <- "none"
        aa_BOND_flag <- c(aa_BOND_flag, 0)
      }else{
        boncollectTable <- c("")
        bondcollect <- c("<ul>")
        allbonds <- unlist(strsplit(bond_reportF, split="[;]"))
        for(i in 1:length(allbonds)){
          if(grepl("-D/", as.character(bond_reportF))){
            aa_BOND_flag <- c(aa_BOND_flag, 1)
            boncollectTable <- c(boncollectTable, "disulfie bond")
            thisbondinfo <- unlist(strsplit(allbonds[i], split="[-]"))
            partnerinfo <- unlist(strsplit(thisbondinfo[3], split=":"))
            thisbondinfoTypeLen <- unlist(strsplit(thisbondinfo[2], split="[/]"))
            bondcollect <- c(bondcollect, "<li>Disulfide bond with ", partnerinfo[3], " at structure position ", partnerinfo[2], ", bond length = ", thisbondinfoTypeLen[2], " &#8491;</li>")
          }
          if(grepl("-N/", as.character(bond_reportF))){
            aa_BOND_flag <- c(aa_BOND_flag, 4)
            boncollectTable <- c(boncollectTable, "nonbonded Van der Waal interaction")
            thisbondinfo <- unlist(strsplit(allbonds[i], split="[-]"))
            partnerinfo <- unlist(strsplit(thisbondinfo[3], split=":"))
            thisbondinfoTypeLen <- unlist(strsplit(thisbondinfo[2], split="[/]"))
            bondcollect <- c(bondcollect, "<li>Nonbonded Van der Waal interaction with ", partnerinfo[3], " at structure position ", partnerinfo[2], ", bond length = ", thisbondinfoTypeLen[2], " &#8491;</li>")
          }
          if(grepl("-H/", as.character(bond_reportF))){
            aa_BOND_flag <- c(aa_BOND_flag, 3)
            boncollectTable <- c(boncollectTable, "hydrogen bond")
            thisbondinfo <- unlist(strsplit(allbonds[i], split="[-]"))
            partnerinfo <- unlist(strsplit(thisbondinfo[3], split=":"))
            thisbondinfoTypeLen <- unlist(strsplit(thisbondinfo[2], split="[/]"))
            bondcollect <- c(bondcollect, "<li>Hydrogen bond with ", partnerinfo[3], " at structure position ", partnerinfo[2], ", bond length = ", thisbondinfoTypeLen[2], " &#8491;</li>")
          }
          if(grepl("-S/", as.character(bond_reportF))){
            aa_BOND_flag <- c(aa_BOND_flag, 2)
            boncollectTable <- c(boncollectTable, "salt-bridge interaction")
            thisbondinfo <- unlist(strsplit(allbonds[i], split="[-]"))
            partnerinfo <- unlist(strsplit(thisbondinfo[3], split=":"))
            thisbondinfoTypeLen <- unlist(strsplit(thisbondinfo[2], split="[/]"))
            bondcollect <- c(bondcollect, "<li>Salt-bridge ionic interaction with ", partnerinfo[3], " at structure position ", partnerinfo[2], ", bond length = ", thisbondinfoTypeLen[2], " &#8491;</li>")
          }
        }
        bondcollect <- c(bondcollect, "</ul>")
        bondinfo_reportF <- paste(bondcollect, collapse="")
        boncollectTable <- unique(boncollectTable)
        boncollectTable_count = length(boncollectTable)

        if(boncollectTable_count > 2){
          bondinfo_reportFTable <- paste(boncollectTable, collapse=",")
          bondinfo_reportFTable <- substring(bondinfo_reportFTable, 2)
        }else{
          bondinfo_reportFTable <- paste(boncollectTable, collapse="")
        }
        
      }
      
      #PTM6
      ptmcollect = c("<ul>")
      ptmcollectTable = c("")
      ptmcollectTable_count = 0;
      aa_PTM_flag <- c()
      acetyl_f <- aa_wise_feature$`Acetylation Detail`
      if(acetyl_f != "-"){
        ptmcollect <- c(ptmcollect, "<li>Acetylation: ")
        acetyl_info_all <- unlist(strsplit(acetyl_f, "[;]"))
        acetylcollect = c("<ul>")
        for(i in 1:length(acetyl_info_all)){
          acetyl_info <- unlist(strsplit(acetyl_info_all[i], "[/]"))
          acetyl_struc <- unlist(strsplit(acetyl_info[1], "[.]"))
          
          acetyl_pdb <- acetyl_struc[1]
          acetyl_chain <- acetyl_struc[2]
          acetyl_struc_pos <- acetyl_struc[3]
          acetyl_dist <- acetyl_info[2]
          
          if(as.numeric(acetyl_dist) == 0.0){
            acetyldetails <- c("<li>", acetyl_pdb, ":", acetyl_chain, ", structure position: ", acetyl_struc_pos, ", ", thisaa1, thisaaindex, " is an acetylation site (", "distance = ", acetyl_dist, "&#8491;)", "</li>")
          }else{
            acetyldetails <- c("<li>", acetyl_pdb, ":", acetyl_chain, ", structure position: ", acetyl_struc_pos, ", distance to ", thisaa1, thisaaindex, ": ", acetyl_dist, "&#8491;", "</li>")
          }
          
          acetyldetails_str <- paste(acetyldetails, collapse = "")
          acetylcollect <- c(acetylcollect, acetyldetails_str)
          
          if(as.numeric(acetyl_info[2]) < 10.0){
            ptmcollectTable <- c(ptmcollectTable, "acetylation")
            aa_PTM_flag <- c(aa_PTM_flag, 1)
          }
        }
        acetylcollect <- c(acetylcollect, "</ul>")
        acetylcollect_all<- paste(acetylcollect, collapse = "")
        ptmcollect <- c(ptmcollect, acetylcollect_all, "</li>")
      }else{
        ptmcollect <- c(ptmcollect, "")
      }
      
      
      methyl_f <- aa_wise_feature$`Methylation Detail`
      if(methyl_f != "-"){
        ptmcollect <- c(ptmcollect, "<li>Methylation: ")
        methyl_info_all <- unlist(strsplit(methyl_f, "[;]"))
        methylcollect = c("<ul>")
        for(i in 1:length(methyl_info_all)){
          methyl_info <- unlist(strsplit(methyl_info_all[i], "[/]"))
          methyl_struc <- unlist(strsplit(methyl_info[1], "[.]"))
          
          methyl_pdb <- methyl_struc[1]
          methyl_chain <- methyl_struc[2]
          methyl_struc_pos <- methyl_struc[3]
          methyl_dist <- methyl_info[2]
          
          if(as.numeric(methyl_dist) == 0.0){
            methyldetails <- c("<li>", methyl_pdb, ":", methyl_chain, ", structure position: ", methyl_struc_pos, ", ", thisaa1, thisaaindex, " is a methylation site (", "distance = ", methyl_dist, "&#8491;)", "</li>")
          }else{
            methyldetails <- c("<li>", methyl_pdb, ":", methyl_chain, ", structure position: ", methyl_struc_pos, ", distance to ", thisaa1, thisaaindex, ": ", methyl_dist, "&#8491;", "</li>")
          }
          
          methyldetails_str <- paste(methyldetails, collapse = "")
          methylcollect <- c(methylcollect, methyldetails_str)
          
          if(as.numeric(methyl_info[2]) < 10.0){
            ptmcollectTable <- c(ptmcollectTable, "methylation")
            aa_PTM_flag <- c(aa_PTM_flag, 2)
          }
        }
        methylcollect <- c(methylcollect, "</ul>")
        methylcollect_all<- paste(methylcollect, collapse = "")
        ptmcollect <- c(ptmcollect, methylcollect_all, "</li>")
      }else{
        ptmcollect <- c(ptmcollect, "")
      }
      
      gcln_f <- aa_wise_feature$`O.GclNAc Detail`
      if(gcln_f != "-"){
        ptmcollect <- c(ptmcollect, "<li>O.GclNAc: ")
        gcln_info_all <- unlist(strsplit(gcln_f, "[;]"))
        gclncollect = c("<ul>")
        for(i in 1:length(gcln_info_all)){
          gcln_info <- unlist(strsplit(gcln_info_all[i], "[/]"))
          gcln_struc <- unlist(strsplit(gcln_info[1], "[.]"))
          
          gcln_pdb <- gcln_struc[1]
          gcln_chain <- gcln_struc[2]
          gcln_struc_pos <- gcln_struc[3]
          gcln_dist <- gcln_info[2]
          
          if(as.numeric(gcln_dist) == 0.0){
            gclndetails <- c("<li>", gcln_pdb, ":", gcln_chain, ", structure position: ", gcln_struc_pos, ", ", thisaa1, thisaaindex, " is a O.GclNAc site (", "distance = ", gcln_dist, "&#8491;)", "</li>")
          }else{
            gclndetails <- c("<li>", gcln_pdb, ":", gcln_chain, ", structure position: ", gcln_struc_pos, ", distance to ", thisaa1, thisaaindex, ": ", gcln_dist, "&#8491;", "</li>")
          }
          
          gclndetails_str <- paste(gclndetails, collapse = "")
          gclncollect <- c(gclncollect, gclndetails_str)
          
          if(as.numeric(gcln_info[2]) < 10.0){
            ptmcollectTable <- c(ptmcollectTable, "O.GclNAc")
            aa_PTM_flag <- c(aa_PTM_flag, 3)
          }
        }
        gclncollect <- c(gclncollect, "</ul>")
        gclncollect_all<- paste(gclncollect, collapse = "")
        ptmcollect <- c(ptmcollect, gclncollect_all, "</li>")
      }else{
        ptmcollect <- c(ptmcollect, "")
      }
      
      phos_f <- aa_wise_feature$`Phosphorylation Detail`
      if(phos_f != "-"){
        ptmcollect <- c(ptmcollect, "<li>Phosphorylation: ")
        phos_info_all <- unlist(strsplit(phos_f, "[;]"))
        phoscollect = c("<ul>")
        for(i in 1:length(phos_info_all)){
          phos_info <- unlist(strsplit(phos_info_all[i], "[/]"))
          phos_struc <- unlist(strsplit(phos_info[1], "[.]"))
          
          phos_pdb <- phos_struc[1]
          phos_chain <- phos_struc[2]
          phos_struc_pos <- phos_struc[3]
          phos_dist <- phos_info[2]
          
          if(as.numeric(phos_dist) == 0.0){
            phosdetails <- c("<li>", phos_pdb, ":", phos_chain, ", structure position: ", phos_struc_pos, ", ", thisaa1, thisaaindex, " is a phosphorylation site (", "distance = ", phos_dist, "&#8491;)", "</li>")
          }else{
            phosdetails <- c("<li>", phos_pdb, ":", phos_chain, ", structure position: ", phos_struc_pos, ", distance to ", thisaa1, thisaaindex, ": ", phos_dist, "&#8491;", "</li>")
          }
          
          phosdetails_str <- paste(phosdetails, collapse = "")
          phoscollect <- c(phoscollect, phosdetails_str)
          
          if(as.numeric(phos_info[2]) < 10.0){
            ptmcollectTable <- c(ptmcollectTable, "phosphorylation")
            aa_PTM_flag <- c(aa_PTM_flag, 4)
          }
        }
        phoscollect <- c(phoscollect, "</ul>")
        phoscollect_all<- paste(phoscollect, collapse = "")
        ptmcollect <- c(ptmcollect, phoscollect_all, "</li>")
      }else{
        ptmcollect <- c(ptmcollect, "")
      }
      
      sumoy_f <- aa_wise_feature$`Sumoylation Detail`
      if(sumoy_f != "-"){
        ptmcollect <- c(ptmcollect, "<li>Sumoylation: ")
        sumoy_info_all <- unlist(strsplit(sumoy_f, "[;]"))
        sumoycollect = c("<ul>")
        for(i in 1:length(sumoy_info_all)){
          sumoy_info <- unlist(strsplit(sumoy_info_all[i], "[/]"))
          sumoy_struc <- unlist(strsplit(sumoy_info[1], "[.]"))
          
          sumoy_pdb <- sumoy_struc[1]
          sumoy_chain <- sumoy_struc[2]
          sumoy_struc_pos <- sumoy_struc[3]
          sumoy_dist <- sumoy_info[2]
          
          if(as.numeric(sumoy_dist) == 0.0){
            sumoydetails <- c("<li>", sumoy_pdb, ":", sumoy_chain, ", structure position: ", sumoy_struc_pos, ", ", thisaa1, thisaaindex, " is a sumoylation site (", "distance = ", sumoy_dist, "&#8491;)", "</li>")
          }else{
            sumoydetails <- c("<li>", sumoy_pdb, ":", sumoy_chain, ", structure position: ", sumoy_struc_pos, ", distance to ", thisaa1, thisaaindex, ": ", sumoy_dist, "&#8491;", "</li>")
          }
          
          sumoydetails_str <- paste(sumoydetails, collapse = "")
          sumoycollect <- c(sumoycollect, sumoydetails_str)
          
          if(as.numeric(sumoy_info[2]) < 10.0){
            ptmcollectTable <- c(ptmcollectTable, "sumoylation")
            aa_PTM_flag <- c(aa_PTM_flag, 5)
          }
        }
        sumoycollect <- c(sumoycollect, "</ul>")
        sumoycollect_all<- paste(sumoycollect, collapse = "")
        ptmcollect <- c(ptmcollect, sumoycollect_all, "</li>")
      }else{
        ptmcollect <- c(ptmcollect, "")
      }
      
      ubiq_f <- aa_wise_feature$`Ubiquitination Detail`
      if(ubiq_f != "-"){
        ptmcollect <- c(ptmcollect, "<li>Ubiquitinition: ")
        ubiq_info_all <- unlist(strsplit(ubiq_f, "[;]"))
        ubiqcollect = c("<ul>")
        for(i in 1:length(ubiq_info_all)){
          ubiq_info <- unlist(strsplit(ubiq_info_all[i], "[/]"))
          ubiq_struc <- unlist(strsplit(ubiq_info[1], "[.]"))
          
          ubiq_pdb <- ubiq_struc[1]
          ubiq_chain <- ubiq_struc[2]
          ubiq_struc_pos <- ubiq_struc[3]
          ubiq_dist <- ubiq_info[2]
          
          if(as.numeric(ubiq_dist) == 0.0){
            ubiqdetails <- c("<li>", ubiq_pdb, ":", ubiq_chain, ", structure position: ", ubiq_struc_pos, ", ", thisaa1, thisaaindex, " is a ubiquitinition site (", "distance = ", ubiq_dist, "&#8491;)", "</li>")
          }else{
            ubiqdetails <- c("<li>", ubiq_pdb, ":", ubiq_chain, ", structure position: ", ubiq_struc_pos, ", distance to ", thisaa1, thisaaindex, ": ", ubiq_dist, "&#8491;", "</li>")
          }
          
          ubiqdetails_str <- paste(ubiqdetails, collapse = "")
          ubiqcollect <- c(ubiqcollect, ubiqdetails_str)
          
          if(as.numeric(ubiq_info[2]) < 10.0){
            ptmcollectTable <- c(ptmcollectTable, "ubiquitination")
            aa_PTM_flag <- c(aa_PTM_flag, 6)
          }
        }
        ubiqcollect <- c(ubiqcollect, "</ul>")
        ubiqcollect_all<- paste(ubiqcollect, collapse = "")
        ptmcollect <- c(ptmcollect, ubiqcollect_all, "</li>")
      }else{
        ptmcollect <- c(ptmcollect, "")
      }
      
      ptmcollect <- c(ptmcollect, "</ul>")
      ptmcollect_str <- paste(ptmcollect, collapse="")
      if(ptmcollect_str == "<ul></ul>"){
        ptminfo_reportFTable <- "none"
        if(structure == 0){
          ptminfo_reportF <- "annotation not available"
        }else{
          ptminfo_reportF <- "annotation not available"
        }
      }else{
        #ptmcollect <- c("<strong>: </strong>", ptmcollect)
        ptminfo_reportF <- paste(ptmcollect, collapse="")
        
        ptmcollectTable <- unique(ptmcollectTable)
        ptminfo_reportFtest <- paste(ptmcollectTable, collapse="")
        if(ptminfo_reportFtest == ""){
          ptminfo_reportFTable <- "none"
        }else{
          ptminfo_reportFTable <- paste(ptmcollectTable, collapse=", ") 
          ptminfo_reportFTable <- substring(ptminfo_reportFTable, 3)
        }
      }
      
      #UP6
      upcollect = c("<ul>")
      upcollectTable <- c("")
      aa_UP_flag <- c()
      func_site <- aa_wise_feature$`Functional Site`
      if(func_site != "-"){
        aa_UP_flag <- c(aa_UP_flag, 1)
        upcollect <- c(upcollect, "<li>Functional site: ")
        upcollectTable <- c(upcollectTable, "functional site")
        active_site <- aa_wise_feature$`Active Site`
        metal_site <- aa_wise_feature$`Metal Binding Site`
        bind_site <- aa_wise_feature$`binding Site`
        site_site <- aa_wise_feature$`Site`
        
        func_site_all = c("<ul>")
        if(active_site != "-"){
          active_site_info <- unlist(strsplit(active_site, ":"))
          func_site_all <- c(func_site_all, "<li>Active site: ", active_site_info[1], " (", trimws(active_site_info[2]), ")</li>")
        }
        if(metal_site != "-"){
          metal_site_info <- unlist(strsplit(metal_site, ":"))
          func_site_all <- c(func_site_all, "<li>Metal binding site: ", metal_site_info[1], " (", trimws(metal_site_info[2]), ")</li>")
        }
        if(bind_site != "-"){
          bind_site_info <- unlist(strsplit(bind_site, ":"))
          func_site_all <- c(func_site_all, "<li>Binding site: ", bind_site_info[1], " (", trimws(bind_site_info[2]), ")</li>")
        }
        if(site_site != "-"){
          site_site_info <- unlist(strsplit(site_site, ":"))
          func_site_all <- c(func_site_all, "<li>Site: ", site_site_info[1], " (", trimws(site_site_info[2]), ")</li>")
        }
        func_site_all <- c(func_site_all, "</ul>")
        upcollect <- c(upcollect, func_site_all, "</li>")
      }else{
        upcollect <- c(upcollect, "")
      }
      
      bind_reg <- aa_wise_feature$`Functional Binding Region`
      if(bind_reg != "-"){
        aa_UP_flag <- c(aa_UP_flag, 2)
        upcollect <- c(upcollect, "<li>Functional/binding region: ")
        upcollectTable <- c(upcollectTable, "functional/binding region")
        
        zinc_fing <- aa_wise_feature$`Zinc Finger`
        dna_bind <- aa_wise_feature$`Dna Binding`
        np_bind <- aa_wise_feature$`Nucleotide Phosphate Binding`
        ca_bind <- aa_wise_feature$`Calcium Binding`
        
        bind_reg_all = c("<ul>")
        if(zinc_fing != "-"){
          zinc_fing_info <- unlist(strsplit(zinc_fing, ":"))
          bind_reg_all <- c(bind_reg_all, "<li>Zinc finger: ", zinc_fing_info[1], " (", trimws(zinc_fing_info[2]), ")</li>")
        }
        if(dna_bind != "-"){
          dna_bind_info <- unlist(strsplit(dna_bind, ":"))
          bind_reg_all <- c(bind_reg_all, "<li>DNA binding region: ", dna_bind_info[1], " (", trimws(dna_bind_info[2]), ")</li>")
        }
        if(np_bind != "-"){
          np_bind_info <- unlist(strsplit(np_bind, ":"))
          bind_reg_all <- c(bind_reg_all, "<li>Nucleotide phosphate binding region: ", np_bind_info[1], " (", trimws(np_bind_info[2]), ")</li>")
        }
        if(ca_bind != "-"){
          ca_bind_info <- unlist(strsplit(ca_bind, ":"))
          bind_reg_all <- c(bind_reg_all, "<li>Ca binding region: ", ca_bind_info[1], " (", trimws(ca_bind_info[2]), ")</li>")
        }
        bind_reg_all <- c(bind_reg_all, "</ul>")
        upcollect <- c(upcollect, bind_reg_all, "</li>")
      }else{
        upcollect <- c(upcollect, "")
      }
      
      seq_motif <- aa_wise_feature$`Sequence Motif Region`
      if(seq_motif != "-"){
        aa_UP_flag <- c(aa_UP_flag, 3)
        upcollect <- c(upcollect, "<li>Sequence motif/region: ")
        upcollectTable <- c(upcollectTable, "sequence motif/region")
        
        region <- aa_wise_feature$`Region`
        repeatr <- aa_wise_feature$`Repeat`
        coiledc <- aa_wise_feature$`Coiled Coil`
        motif <- aa_wise_feature$`Motif`
        
        seq_motif_all = c("<ul>")
        if(region != "-"){
          region_info <- unlist(strsplit(region, ";"))
          seq_motif_all <- c(seq_motif_all, "<li>Region: ")
          for(n_region_info in 1:length(region_info)){
            each_region_info <- region_info[n_region_info]
            each_region_info_parts <- unlist(strsplit(each_region_info, ":"))
            if(n_region_info == length(region_info)){
              seq_motif_all <- c(seq_motif_all, each_region_info_parts[1], " (", trimws(each_region_info_parts[2]), ")")
            }else{
              seq_motif_all <- c(seq_motif_all, each_region_info_parts[1], " (", trimws(each_region_info_parts[2]), "), ")
            }
          }
          seq_motif_all <- c(seq_motif_all,"</li>")
        }
        if(repeatr != "-"){
          repeatr_info <- unlist(strsplit(repeatr, ":"))
          seq_motif_all <- c(seq_motif_all, "<li>Repeat: ", repeatr_info[1], " (", trimws(repeatr_info[2]), ")</li>")
        }
        if(coiledc != "-"){
          coiledc_info <- unlist(strsplit(coiledc, ":"))
          seq_motif_all <- c(seq_motif_all, "<li>Coiled coil: ", coiledc_info[1], " (", trimws(coiledc_info[2]), ")</li>")
        }
        if(motif != "-"){
          motif_info <- unlist(strsplit(motif, ":"))
          seq_motif_all <- c(seq_motif_all, "<li>Motif: ", motif_info[1], " (", trimws(motif_info[2]), ")</li>")
        }
        seq_motif_all <- c(seq_motif_all, "</ul>")
        upcollect <- c(upcollect, seq_motif_all, "</li>")
      }else{
        upcollect <- c(upcollect, "")
      }
      
      mod_domain <- aa_wise_feature$`Modular Domain`
      if(mod_domain != "-"){
        aa_UP_flag <- c(aa_UP_flag, 4)
        upcollect <- c(upcollect, "<li>Modular domain: ")
        upcollectTable <- c(upcollectTable, "modular domain")
        
        domain <- aa_wise_feature$`Domain`
        top_domain <- aa_wise_feature$`Topological Domain`
        trans_domain <- aa_wise_feature$`Transmembrane`
        intra_domain <- aa_wise_feature$`Intramembrane`
        
        mod_domain_all = c("<ul>")
        if(domain != "-"){
          domain_info <- unlist(strsplit(domain, ":"))
          mod_domain_all <- c(mod_domain_all, "<li>Domain: ", domain_info[1], " (", trimws(domain_info[2]), ")</li>")
        }
        if(top_domain != "-"){
          top_domain_info <- unlist(strsplit(top_domain, ":"))
          mod_domain_all <- c(mod_domain_all, "<li>Topological domain: ", top_domain_info[1], " (", trimws(top_domain_info[2]), ")</li>")
        }
        if(trans_domain != "-"){
          trans_domain_info <- unlist(strsplit(trans_domain, ":"))
          mod_domain_all <- c(mod_domain_all, "<li>Transmembrane: ", trans_domain_info[1], " (", trimws(trans_domain_info[2]), ")</li>")
        }
        if(intra_domain != "-"){
          intra_domain_info <- unlist(strsplit(intra_domain, ":"))
          mod_domain_all <- c(mod_domain_all, "<li>Intramembrane: ", intra_domain_info[1], " (", trimws(intra_domain_info[2]), ")</li>")
        }
        mod_domain_all <- c(mod_domain_all, "</ul>")
        upcollect <- c(upcollect, mod_domain_all, "</li>")
      }else{
        upcollect <- c(upcollect, "")
      }
      
      mol_proc <- aa_wise_feature$`Molecular Processing`
      if(mol_proc != "-"){
        aa_UP_flag <- c(aa_UP_flag, 5)
        upcollect <- c(upcollect, "<li>Molecular processing associated region: ")
        upcollectTable <- c(upcollectTable, "molecular processing associated region")
        
        peptide <- aa_wise_feature$`Peptide`
        propep <- aa_wise_feature$`Propetide`
        trans_pep <- aa_wise_feature$`Transit Peptide`
        sig_pep <- aa_wise_feature$`Signal Peptide`
        
        mol_proc_all = c("<ul>")
        if(peptide != "-"){
          peptide_info <- unlist(strsplit(peptide, ":"))
          mol_proc_all <- c(mol_proc_all, "<li>Peptide: ", peptide_info[1], " (", trimws(peptide_info[2]), ")</li>")
        }
        if(propep != "-"){
          propep_info <- unlist(strsplit(propep, ":"))
          mol_proc_all <- c(mol_proc_all, "<li>Pro-peptide: ", propep_info[1], " (", trimws(propep_info[2]), ")</li>")
        }
        if(trans_pep != "-"){
          trans_pep_info <- unlist(strsplit(trans_pep, ":"))
          mol_proc_all <- c(mol_proc_all, "<li>Transit peptide: ", trans_pep_info[1], " (", trimws(trans_pep_info[2]), ")</li>")
        }
        if(sig_pep != "-"){
          sig_pep_info <- unlist(strsplit(sig_pep, ":"))
          mol_proc_all <- c(mol_proc_all, "<li>Signal peptide: ", sig_pep_info[1], " (", trimws(sig_pep_info[2]), ")</li>")
        }
        mol_proc_all <- c(mol_proc_all, "</ul>")
        upcollect <- c(upcollect, mol_proc_all, "</li>")
      }else{
        upcollect <- c(upcollect, "")
      }
      
      mod_res <- aa_wise_feature$`Modified Residues`
      if(mod_res != "-"){
        aa_UP_flag <- c(aa_UP_flag, 6)
        upcollect <- c(upcollect, "<li>Modified residues: ")
        upcollectTable <- c(upcollectTable, "modified residues")
        
        modified <- aa_wise_feature$`Modified Residue`
        lipid <- aa_wise_feature$`Lipidation`
        glycos <- aa_wise_feature$`Glycosylation`
        clinks <- aa_wise_feature$`Cross Links`
        dsulf <- aa_wise_feature$`Disulfide Bond`
        
        mod_res_all = c("<ul>")
        if(modified != "-"){
          peptide_info <- unlist(strsplit(modified, ":"))
          mod_res_all <- c(mod_res_all, "<li>PTM-modified residues: ", peptide_info[1], " (", trimws(peptide_info[2]), ")</li>")
        }
        if(lipid != "-"){
          lipid_info <- unlist(strsplit(lipid, ":"))
          mod_res_all <- c(mod_res_all, "<li>Lipidation: ", lipid_info[1], " (", trimws(lipid_info[2]), ")</li>")
        }
        if(glycos != "-"){
          glycos_info <- unlist(strsplit(glycos, ":"))
          mod_res_all <- c(mod_res_all, "<li>Glycosylation: ", glycos_info[1], " (", trimws(glycos_info[2]), ")</li>")
        }
        if(clinks != "-"){
          clinks_info <- unlist(strsplit(clinks, ":"))
          mod_res_all <- c(mod_res_all, "<li>Cross-links: ", clinks_info[1], " (", trimws(clinks_info[2]), ")</li>")
        }
        if(dsulf != "-"){
          dsulf_info <- unlist(strsplit(dsulf, ":"))
          mod_res_all <- c(mod_res_all, "<li>Disulfide bond: ", dsulf_info[1], " (", trimws(dsulf_info[2]), ")</li>")
        }
        mod_res_all <- c(mod_res_all, "</ul>")
        upcollect <- c(upcollect, mod_res_all, "</li>")
      }else{
        upcollect <- c(upcollect, "")
      }
      
      upcollect <- c(upcollect, "</ul>")
      upcollect_str <- paste(upcollect, collapse="")
      
      if(upcollect_str == "<ul></ul>"){
        upinfo_reportF <- "annotation not available"
        upinfo_reportFTable <- "none"
      }else{
        #upcollect <- c("[amino acid and its position in sequence (feature description)]<strong>: </strong>", upcollect)
        upinfo_reportF <- paste(upcollect, collapse="")
        upinfo_reportFTable <- paste(upcollectTable, collapse = ", ")
        upinfo_reportFTable <- substring(upinfo_reportFTable, 3)
      }
      
      #pclass info/table - start
      
      # summary report - pclass#
      gene_protein_trans_class <- subset(GPTLPC_list, geneName == input$geneSelected)
      pclassIDreport = gene_protein_trans_class$pclass
      
      pclassIDreport = gene_protein_trans_class$pclass #pclass
      pclassIds = unlist(strsplit(pclassIDreport, ", "))
      pclassIds = c("c0", pclassIds)
      
      PC24_PF40_path <- subset(PC24_PF40_pathogenic, pcindex %in% pclassIds)
      PC24_PF40_pop <- subset(PC24_PF40_population, pcindex %in% pclassIds)
      
      npclasses = nrow(PC24_PF40_path)
      #print(npclasses)
      table_code <- c("<strong>III) Comparison of the features of ", aaShorter_str, " with Protein class-specific pathogenic and population variant associated features: </strong><br>","[Features of the selected amino acid is shown in","<span style='color: #ff0000;'>"," red</span> (","<span style='color: #0000ff;'>","blue","</span>) if matched with ", "<span style='color: #ff0000;'>", "pathogenic</span> (","<span style='color: #0000ff;'>","population</span>) variant-associated features]", "<br></br>")
      for(i in 1:npclasses){
        pcname <- PC24_PF40_path$pcname[i]
        pcnamelower <- PC24_PF40_path$pcnamelower[i]
        
        local_conformSS3Table <- local_conformSS3
        local_conformSS3Table <- paste("<span style='color: #000000;'>",local_conformSS3Table,"</span>", collapse = "");
        
        local_conformSS8Table <- local_conformSS8
        local_conformSS8Table <- paste("<span style='color: #000000;'>",local_conformSS8Table,"</span>", collapse = "");
        
        solv_exposureTable <- exposure_reportF
        solv_exposureTable <- paste("<span style='color: #000000;'>",solv_exposureTable,"</span>", collapse = "");
        
        #pathogenic feature for pclass -- start#
        pss3path_df <- PC24_PF40_path[i,4:6]
        pss3path_list <- c()
        for(fn in 1:ncol(pss3path_df))
        {
          pfindex = as.character(names(pss3path_df[1,fn]))
          fflag = as.integer(pss3path_df[fn])

          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss3path_list <- c(pss3path_list, fdes)
            if(aa_SS3_flag == fn){
              local_conformSS3Table <- local_conformSS3
              local_conformSS3Table <- paste("<span style='color: #ff0000;'>",local_conformSS3Table,"</span>", collapse = "");
            }
          }
        }
        if(length(pss3path_list) > 0){
          pss3path_str <- paste(pss3path_list, collapse = ", ")
        }else{
          pss3path_str <- "none"
        }
        
        pss8path_df <- PC24_PF40_path[i,7:14]
        pss8path_list <- c()
        for(fn in 1:ncol(pss8path_df))
        {
          pfindex = as.character(names(pss8path_df[1,fn]))
          fflag = as.integer(pss8path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss8path_list <- c(pss8path_list, fdes)
            
            if(aa_SS8_flag == fn){
              local_conformSS8Table <- local_conformSS8
              local_conformSS8Table <- paste("<span style='color: #ff0000;'>",local_conformSS8Table,"</span>", collapse = "");
            }
          }
        }
        if(length(pss8path_list) > 0){
          pss8path_str <- paste(pss8path_list, collapse = ", ")
        }else{
          pss8path_str <- "none"
        }
        
        pasapath_df <- PC24_PF40_path[i,15:19]
        pasapath_list <- c()
        for(fn in 1:ncol(pasapath_df))
        {
          pfindex = as.character(names(pasapath_df[1,fn]))
          fflag = as.integer(pasapath_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pasapath_list <- c(pasapath_list, fdes)
            if(aa_ASA_flag == fn){
              solv_exposureTable <- exposure_reportF
              solv_exposureTable <- paste("<span style='color: #ff0000;'>",solv_exposureTable,"</span>", collapse = "");
            }
          }
        }
        if(length(pasapath_list) > 0){
          pasapath_str <- paste(pasapath_list, collapse = ", ")
        }else{
          pasapath_str <- "none"
        }
        
        pcp8path_df <- PC24_PF40_path[i,20:27]
        pcp8path_list <- c()
        for(fn in 1:ncol(pcp8path_df))
        {
          pfindex = as.character(names(pcp8path_df[1,fn]))
          fflag = as.integer(pcp8path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pcp8path_list <- c(pcp8path_list, fdes)
          }
        }
        if(length(pcp8path_list) > 0){
          pcp8path_str <- paste(pcp8path_list, collapse = ", ")
          pcp8path_str <- paste(pcp8path_str, "amino acids", sep = " ")
        }else{
          pcp8path_str <- "none"
        }
        
        pib4path_df <- PC24_PF40_path[i,28:31]
        pib4path_list <- c()
        for(fn in 1:ncol(pib4path_df))
        {
          pfindex = as.character(names(pib4path_df[1,fn]))
          fflag = as.integer(pib4path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pib4path_list <- c(pib4path_list, fdes)
          }
        }
        if(length(pib4path_list) > 0){
          pib4path_str <- paste(pib4path_list, collapse = ", ")
          pib4path_str <- paste(pib4path_str, "interactions", sep = " ")
        }else{
          pib4path_str <- "none"
        }
        
        pptm6path_df <- PC24_PF40_path[i,32:37]
        pptm6path_list <- c()
        for(fn in 1:ncol(pptm6path_df))
        {
          pfindex = as.character(names(pptm6path_df[1,fn]))
          fflag = as.integer(pptm6path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pptm6path_list <- c(pptm6path_list, fdes)
          }
        }
        if(length(pptm6path_list) > 0){
          pptm6path_str <- paste(pptm6path_list, collapse = ", ")
        }else{
          pptm6path_str <- "none"
        }
        
        pup6path_df <- PC24_PF40_path[i,38:43]
        pup6path_list <- c()
        for(fn in 1:ncol(pup6path_df))
        {
          pfindex = as.character(names(pup6path_df[1,fn]))
          fflag = as.integer(pup6path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pup6path_list <- c(pup6path_list, fdes)
          }
        }
        if(length(pup6path_list) > 0){
          pup6path_str <- paste(pup6path_list, collapse = ", ")
        }else{
          pup6path_str <- "none"
        }
        #pathogenic feature for pclass -- end#
        
        #popultaion feature for pclass -- start#
        pss3pop_df <- PC24_PF40_pop[i,4:6]
        pss3pop_list <- c()
        for(fn in 1:ncol(pss3pop_df))
        {
          pfindex = as.character(names(pss3pop_df[1,fn]))
          fflag = as.integer(pss3pop_df[fn])

          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss3pop_list <- c(pss3pop_list, fdes)
            if(aa_SS3_flag == fn){
              local_conformSS3Table <- local_conformSS3
              local_conformSS3Table <- paste("<span style='color: #0000ff;'>",local_conformSS3Table,"</span>", collapse = "");
            }
          }
        }
        if(length(pss3pop_list) > 0){
          pss3pop_str <- paste(pss3pop_list, collapse = ", ")
        }else{
          pss3pop_str <- "none"
        }
        
        pss8pop_df <- PC24_PF40_pop[i,7:14]
        pss8pop_list <- c()
        for(fn in 1:ncol(pss8pop_df))
        {
          pfindex = as.character(names(pss8pop_df[1,fn]))
          fflag = as.integer(pss8pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss8pop_list <- c(pss8pop_list, fdes)
            if(aa_SS8_flag == fn){
              local_conformSS8Table <- local_conformSS8
              local_conformSS8Table <- paste("<span style='color: #0000ff;'>",local_conformSS8Table,"</span>", collapse = "");
            }
          }
        }
        if(length(pss8pop_list) > 0){
          pss8pop_str <- paste(pss8pop_list, collapse = ", ")
        }else{
          pss8pop_str <- "none"
        }
        
        pasapop_df <- PC24_PF40_pop[i,15:19]
        pasapop_list <- c()
        for(fn in 1:ncol(pasapop_df))
        {
          pfindex = as.character(names(pasapop_df[1,fn]))
          fflag = as.integer(pasapop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pasapop_list <- c(pasapop_list, fdes)
            if(aa_ASA_flag == fn){
              solv_exposureTable <- exposure_reportF
              solv_exposureTable <- paste("<span style='color: #0000ff;'>",solv_exposureTable,"</span>", collapse = "");
            }
          }
        }
        if(length(pasapop_list) > 0){
          pasapop_str <- paste(pasapop_list, collapse = ", ")
        }else{
          pasapop_str <- "none"
        }
        
        pcp8pop_df <- PC24_PF40_pop[i,20:27]
        pcp8pop_list <- c()
        for(fn in 1:ncol(pcp8pop_df))
        {
          pfindex = as.character(names(pcp8pop_df[1,fn]))
          fflag = as.integer(pcp8pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pcp8pop_list <- c(pcp8pop_list, fdes)
          }
        }
        if(length(pcp8pop_list) > 0){
          pcp8pop_str <- paste(pcp8pop_list, collapse = ", ")
          pcp8pop_str <- paste(pcp8pop_str, "amino acids", sep = " ")
        }else{
          pcp8pop_str <- "none"
        }
        
        pib4pop_df <- PC24_PF40_pop[i,28:31]
        pib4pop_list <- c()
        for(fn in 1:ncol(pib4pop_df))
        {
          pfindex = as.character(names(pib4pop_df[1,fn]))
          fflag = as.integer(pib4pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pib4pop_list <- c(pib4pop_list, fdes)
          }
        }
        if(length(pib4pop_list) > 0){
          pib4pop_str <- paste(pib4pop_list, collapse = ", ")
          pib4pop_str <- paste(pib4pop_str, "interactions", sep = " ")
        }else{
          pib4pop_str <- "none"
        }
        
        
        pptm6pop_df <- PC24_PF40_pop[i,32:37]
        pptm6pop_list <- c()
        for(fn in 1:ncol(pptm6pop_df))
        {
          pfindex = as.character(names(pptm6pop_df[1,fn]))
          fflag = as.integer(pptm6pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pptm6pop_list <- c(pptm6pop_list, fdes)
          }
        }
        if(length(pptm6pop_list) > 0){
          pptm6pop_str <- paste(pptm6pop_list, collapse = ", ")
        }else{
          pptm6pop_str <- "none"
        }
        
        pup6pop_df <- PC24_PF40_pop[i,38:43]
        pup6pop_list <- c()
        for(fn in 1:ncol(pup6pop_df))
        {
          pfindex = as.character(names(pup6pop_df[1,fn]))
          fflag = as.integer(pup6pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pup6pop_list <- c(pup6pop_list, fdes)
          }
        }
        if(length(pup6pop_list) > 0){
          pup6pop_str <- paste(pup6pop_list, collapse = ", ")
        }else{
          pup6pop_str <- "none"
        }
        #population feature for pclass -- end#
        if(i == 1){
          pcname <- "All protein"
        }
        
        # color bond info
        aa_BOND_flag <- aa_BOND_flag[aa_BOND_flag != 0]
        aa_BOND_flag <- unique(aa_BOND_flag)
        #print(aa_BOND_flag)
        if(length(aa_BOND_flag) > 0){
          bondinfo_reportFTable_strlist <- c()
          for(aa_feat in 1:length(aa_BOND_flag)){
            aa_feature_flag <- as.integer(aa_BOND_flag[aa_feat])
            flag_path <- as.integer(pib4path_df[aa_feature_flag])
            flag_pop <- as.integer(pib4pop_df[aa_feature_flag])
            
            if((aa_feature_flag == 1) && (flag_path == 1)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #ff0000;'>disulfide bond</span>");
            }else if((aa_feature_flag == 1) && (flag_pop == 1)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #0000ff;'>disulfide bond</span>");
            }else if((aa_feature_flag == 1) && (flag_pop == 0) && (flag_path == 0)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #000000;'>disulfide bond</span>");
            }
            
            if((aa_feature_flag == 2) && (flag_path == 1)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #ff0000;'>salt-bridge interaction</span>");
            }else if((aa_feature_flag == 2) && (flag_pop == 1)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #0000ff;'>salt-bridge interaction</span>");
            }else if((aa_feature_flag == 2) && (flag_pop == 0) && (flag_path == 0)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #000000;'>salt-bridge interaction</span>");
            }
            
            if((aa_feature_flag == 3) && (flag_path == 1)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #ff0000;'>hydrogen bond</span>");
            }else if((aa_feature_flag == 3) && (flag_pop == 1)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #0000ff;'>hydrogen bond</span>");
            }else if((aa_feature_flag == 3) && (flag_pop == 0) && (flag_path == 0)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #000000;'>hydrogen bond</span>");
            }
            
            if((aa_feature_flag == 4) && (flag_path == 1)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #ff0000;'>nonbonded Van der Waals interaction</span>");
            }else if((aa_feature_flag == 4) && (flag_pop == 1)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #0000ff;'>nonbonded Van der Waals interaction</span>");
            }else if((aa_feature_flag == 4) && (flag_pop == 0) && (flag_path == 0)){
              bondinfo_reportFTable_strlist <- c(bondinfo_reportFTable_strlist, "<span style='color: #000000;'>nonbonded Van der Waals interaction</span>");
            }
          }
          
          bondinfo_reportFTable <- paste(bondinfo_reportFTable_strlist, collapse = ", ")
          #print(bondinfo_reportFTable)
        }
        
        # color CHEM info
        aa_CHEM_flag <- aa_CHEM_flag[aa_CHEM_flag != 0]
        aa_CHEM_flag <- unique(aa_CHEM_flag)
        #print(aa_CHEM_flag)
        if(length(aa_CHEM_flag) > 0){
          CHEMinfo_reportFTable_strlist <- c()
          for(aa_feat in 1:length(aa_CHEM_flag)){
            aa_feature_flag <- as.integer(aa_CHEM_flag[aa_feat])
            flag_path <- as.integer(pcp8path_df[aa_feature_flag])
            flag_pop <- as.integer(pcp8pop_df[aa_feature_flag])
            
            if((aa_feature_flag == 1) && (flag_path == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>hydrophobic</span>");
            }else if((aa_feature_flag == 1) && (flag_pop == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>hydrophobic</span>");
            }else if((aa_feature_flag == 1) && (flag_pop == 0) && (flag_path == 0)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #000000;'>hydrophobic</span>");
            }
            
            if((aa_feature_flag == 2) && (flag_path == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>aliphatic</span>");
            }else if((aa_feature_flag == 2) && (flag_pop == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>aliphatic</span>");
            }else if((aa_feature_flag == 2) && (flag_pop == 0) && (flag_path == 0)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #000000;'>aliphatic</span>");
            }
            
            if((aa_feature_flag == 3) && (flag_path == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>aromatic</span>");
            }else if((aa_feature_flag == 3) && (flag_pop == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>aromatic</span>");
            }else if((aa_feature_flag == 3) && (flag_pop == 0) && (flag_path == 0)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #000000;'>aromatic</span>");
            }
            
            if((aa_feature_flag == 4) && (flag_path == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>polar</span>");
            }else if((aa_feature_flag == 4) && (flag_pop == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>polar</span>");
            }else if((aa_feature_flag == 4) && (flag_pop == 0) && (flag_path == 0)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #000000;'>polar</span>");
            }
            
            if((aa_feature_flag == 5) && (flag_path == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>positively-charged</span>");
            }else if((aa_feature_flag == 5) && (flag_pop == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>positively-charged</span>");
            }else if((aa_feature_flag == 5) && (flag_pop == 0) && (flag_path == 0)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #000000;'>positively-charged</span>");
            }
            
            if((aa_feature_flag == 6) && (flag_path == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>negatively-charged</span>");
            }else if((aa_feature_flag == 6) && (flag_pop == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>negatively-charged</span>");
            }else if((aa_feature_flag == 6) && (flag_pop == 0) && (flag_path == 0)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #000000;'>negatively-charged</span>");
            }
            
            if((aa_feature_flag == 7) && (flag_path == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>neutral</span>");
            }else if((aa_feature_flag == 7) && (flag_pop == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>neutral</span>");
            }else if((aa_feature_flag == 7) && (flag_pop == 0) && (flag_path == 0)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #000000;'>neutral</span>");
            }
            
            if((aa_feature_flag == 8) && (flag_path == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>special</span>");
            }else if((aa_feature_flag == 8) && (flag_pop == 1)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>special</span>");
            }else if((aa_feature_flag == 8) && (flag_pop == 0) && (flag_path == 0)){
              CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "<span style='color: #000000;'>special</span>");
            }
          }
          
          #CHEMinfo_reportFTable_strlist <- c(CHEMinfo_reportFTable_strlist, "amino acid")
          chem_reportF_strTable <- paste(CHEMinfo_reportFTable_strlist, collapse = ", ")
          #print(CHEMinfo_reportFTable)
        }
        # color PTM info
        aa_PTM_flag <- aa_PTM_flag[aa_PTM_flag != 0]
        aa_PTM_flag <- unique(aa_PTM_flag)
        #print(aa_PTM_flag)
        if(length(aa_PTM_flag) > 0){
          PTMinfo_reportFTable_strlist <- c()
          for(aa_feat in 1:length(aa_PTM_flag)){
            aa_feature_flag <- as.integer(aa_PTM_flag[aa_feat])
            flag_path <- as.integer(pptm6path_df[aa_feature_flag])
            flag_pop <- as.integer(pptm6pop_df[aa_feature_flag])
            
            if((aa_feature_flag == 1) && (flag_path == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>acetylation</span>");
            }else if((aa_feature_flag == 1) && (flag_pop == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>acetylation</span>");
            }else if((aa_feature_flag == 1) && (flag_pop == 0) && (flag_path == 0)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #000000;'>acetylation</span>");
            }
            
            if((aa_feature_flag == 2) && (flag_path == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>methylation</span>");
            }else if((aa_feature_flag == 2) && (flag_pop == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>methylation</span>");
            }else if((aa_feature_flag == 2) && (flag_pop == 0) && (flag_path == 0)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #000000;'>methylation</span>");
            }
            
            if((aa_feature_flag == 3) && (flag_path == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>O.GclNAc</span>");
            }else if((aa_feature_flag == 3) && (flag_pop == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>O.GclNAc</span>");
            }else if((aa_feature_flag == 3) && (flag_pop == 0) && (flag_path == 0)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #000000;'>O.GclNAc</span>");
            }
            
            if((aa_feature_flag == 4) && (flag_path == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>phosphorylation</span>");
            }else if((aa_feature_flag == 4) && (flag_pop == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>phosphorylation</span>");
            }else if((aa_feature_flag == 4) && (flag_pop == 0) && (flag_path == 0)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #000000;'>phosphorylation</span>");
            }
            
            if((aa_feature_flag == 5) && (flag_path == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>sumoylation</span>");
            }else if((aa_feature_flag == 5) && (flag_pop == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>sumoylation</span>");
            }else if((aa_feature_flag == 5) && (flag_pop == 0) && (flag_path == 0)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #000000;'>sumoylation</span>");
            }
            
            if((aa_feature_flag == 6) && (flag_path == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #ff0000;'>ubiquitination</span>");
            }else if((aa_feature_flag == 6) && (flag_pop == 1)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #0000ff;'>ubiquitination</span>");
            }else if((aa_feature_flag == 6) && (flag_pop == 0) && (flag_path == 0)){
              PTMinfo_reportFTable_strlist <- c(PTMinfo_reportFTable_strlist, "<span style='color: #000000;'>ubiquitination</span>");
            }
          }
          
          ptminfo_reportFTable <- paste(PTMinfo_reportFTable_strlist, collapse = ", ")
          #print(PTMinfo_reportFTable)
        }
        
        # color UP info
        aa_UP_flag <- aa_UP_flag[aa_UP_flag != 0]
        aa_UP_flag <- unique(aa_UP_flag)
        #print(aa_UP_flag)
        if(length(aa_UP_flag) > 0){
          UPinfo_reportFTable_strlist <- c()
          for(aa_feat in 1:length(aa_UP_flag)){
            aa_feature_flag <- as.integer(aa_UP_flag[aa_feat])
            flag_path <- as.integer(pup6path_df[aa_feature_flag])
            flag_pop <- as.integer(pup6pop_df[aa_feature_flag])
            
            if((aa_feature_flag == 1) && (flag_path == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #ff0000;'>functional sites</span>");
            }else if((aa_feature_flag == 1) && (flag_pop == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #0000ff;'>functional sites</span>");
            }else if((aa_feature_flag == 1) && (flag_pop == 0) && (flag_path == 0)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #000000;'>functional sites</span>");
            }
            
            if((aa_feature_flag == 2) && (flag_path == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #ff0000;'>functional/binding regions</span>");
            }else if((aa_feature_flag == 2) && (flag_pop == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #0000ff;'>functional/binding regions</span>");
            }else if((aa_feature_flag == 2) && (flag_pop == 0) && (flag_path == 0)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #000000;'>functional/binding regions</span>");
            }
            
            if((aa_feature_flag == 3) && (flag_path == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #ff0000;'>sequence motifs/regions</span>");
            }else if((aa_feature_flag == 3) && (flag_pop == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #0000ff;'>sequence motifs/regions</span>");
            }else if((aa_feature_flag == 3) && (flag_pop == 0) && (flag_path == 0)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #000000;'>sequence motifs/regions</span>");
            }
            
            if((aa_feature_flag == 4) && (flag_path == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #ff0000;'>modular domains</span>");
            }else if((aa_feature_flag == 4) && (flag_pop == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #0000ff;'>modular domains</span>");
            }else if((aa_feature_flag == 4) && (flag_pop == 0) && (flag_path == 0)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #000000;'>modular domains</span>");
            }
            
            if((aa_feature_flag == 5) && (flag_path == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #ff0000;'>molecular-processing-associated regions</span>");
            }else if((aa_feature_flag == 5) && (flag_pop == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #0000ff;'>molecular-processing-associated regions</span>");
            }else if((aa_feature_flag == 5) && (flag_pop == 0) && (flag_path == 0)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #000000;'>molecular-processing-associated regions</span>");
            }
            
            if((aa_feature_flag == 6) && (flag_path == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #ff0000;'>modified residues</span>");
            }else if((aa_feature_flag == 6) && (flag_pop == 1)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #0000ff;'>modified residues</span>");
            }else if((aa_feature_flag == 6) && (flag_pop == 0) && (flag_path == 0)){
              UPinfo_reportFTable_strlist <- c(UPinfo_reportFTable_strlist, "<span style='color: #000000;'>modified residues</span>");
            }
          }
          
          upinfo_reportFTable <- paste(UPinfo_reportFTable_strlist, collapse = ", ")
          #print(UPinfo_reportFTable)
        }
        table_code <- c(table_code,
                        "<strong>",LETTERS[i], ") ", pcname, " class", "</strong>", "<br>",                
                        "<table style='border-collapse: collapse; width: 100%; height: 144px;' border='1'>",
                        "<tbody>",
                        "<tr style='height: 18px;'>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;Feature category</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;Feature of this amino acid:", aaShorter_str, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;<span style='color: #ff0000;'>Pathogenic</span>-variant-associated features in <strong>", pcnamelower, "</strong></td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;<span style='color: #0000ff;'>Population</span>-variant-associated features in <strong>", pcnamelower, "</strong></td>",
                        "</tr>",
                        "<tr style='height: 18px;'>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;1. 3-class secondary structure:</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;",local_conformSS3Table,"</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pss3path_str, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pss3pop_str, "</td>",
                        "</tr>",
                        "<tr style='height: 18px;'>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;2. 8-class secondary structure:</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;",local_conformSS8Table,"</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pss8path_str, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pss8pop_str, "</td>",
                        "</tr>",
                        "<tr style='height: 18px;'>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;3. Residue exposure level:</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;",solv_exposureTable,"</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pasapath_str, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pasapop_str, "</td>",
                        "</tr>",
                        "<tr style='height: 18px;'>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;4. Physicochemical property:</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;",chem_reportF_strTable,"</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pcp8path_str, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pcp8pop_str, "</td>",
                        "</tr>",
                        "<tr style='height: 18px;'>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;5. Protein-protein interaction/bond type:</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", bondinfo_reportFTable, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pib4path_str, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pib4pop_str, "</td>",
                        "</tr>",
                        "<tr style='height: 18px;'>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;6. Post-translational modification (<10 &#8491; distance):</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", ptminfo_reportFTable, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pptm6path_str, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pptm6pop_str, "</td>",
                        "</tr>",
                        "<tr style='height: 18px;'>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;7. UniProt-based functional feature:</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", upinfo_reportFTable, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pup6path_str, "</td>",
                        "<td style='width: 25%; height: 18px;'>&nbsp;", pup6pop_str, "</td>",
                        "</tr>",
                        "</tbody>",
                        "</table>", "<br></br>")
      }
      # #pclass info/table - end
      
      table_code_str <- paste(table_code, collapse = "")
      
      aa_pop_var_details <- aa_wise_feature$`Population Mutation (details)`
      aa_pop_var_details <- gsub("AC", "allele count", aa_pop_var_details)
      aa_pop_var_details <- gsub("AN", "allele number", aa_pop_var_details)
      aa_pop_var_details <- gsub("AF", "allele frequency", aa_pop_var_details)
      aa_path_var_details <- aa_wise_feature$`Pathogenic Mutation (details)`
      

      HTML(
        paste(
          "<b>","<font size = +1 color = Black face = Arial>","<strong>I) Known missense variations of amino acid residue,","</font>","</b>",aaFull_reportF_str,"</strong><br></br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>gnomAD database (population variants): </strong>","</font>","</b>",aa_pop_var_details, "<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>ClinVar/HGMD database (likely-/pathogenic variants): </strong>","</font>","</b>",aa_path_var_details, "<br></br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>II) Features of amino acid residue,","</font>","</b>",aaFull_reportF_str,"</strong><br></br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>1. 3-class secondary structure: </strong>","</font>","</b>",local_conformSS3, "<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>2. 8-class secondary structure: </strong>","</font>","</b>",local_conformSS8, "<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>3. Exposure to solvent: </strong>","</font>","</b>",solv_exposure,"<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>4. Physicochemical property: </strong>","</font>","</b>",chem_reportF_str, "amino acid", "<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>5. Protein-protein interactions: </strong>","</font>","</b>",bondinfo_reportF,"<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>6. Distance to post-translational modification sites on 3D structure: </strong>","</font>","</b>",ptminfo_reportF,"<br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>7. UniProt-based functional features: </strong>","</font>","</b>",upinfo_reportF,
          "<br></br>",
          table_code_str
        )
      )
    }
  })
  
  report_summReport <- eventReactive({
    input$aa_Selected
    input$reportTsubmit
  },{
    #reportGeneName <- input$reportTgeneSelected  
    if(input$aa_Selected != ""){  
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      report_gene_nameF = paste("gene_wise_info/",input$reportTgeneSelected,".txt",sep='')
      report_gene_wise_infoF <- read_delim(report_gene_nameF, "\t", escape_double = FALSE, trim_ws = TRUE)
      geneNsummary <- input$reportTgeneSelected
      
      aa_wise_feature <- subset(report_gene_wise_infoF, `Amino Acid Index` == input$aa_Selected)
      #AA
      amino_acid_reportF <- aa_wise_feature$`Amino Acid`
      
      thisaa1 <- amino_acid_reportF
      thisaaindex <- input$aa_Selected
      this_aa_info <- subset(aa_info, `one` == thisaa1)
      thisaa3 <- this_aa_info$three
      thisaafull <- this_aa_info$full
      thisaaprp <- this_aa_info$prp
      thisaacomm <- this_aa_info$comm
      
      aa_mention <- c(thisaafull, " ", thisaaindex, " of <em>", geneNsummary, "</em>")
      aa_mention_str <- paste(aa_mention, collapse = "")
      
      line_aa <- c("<strong>Amino Acid Residue:</strong> ", thisaafull, " (", thisaa3, "/", thisaa1, "), position: ", thisaaindex)
      line_aa_str <- paste(line_aa, collapse = "")
      
      path_var <- aa_wise_feature$`Pathogenic Mutation`
      pop_var <- aa_wise_feature$`Population Mutation`
      
      path_var_detail <- aa_wise_feature$`Pathogenic Mutation (details)`
      pop_var_detail <- aa_wise_feature$`Population Mutation (details)`
      pop_var_detail <- gsub("AC", "allele count", pop_var_detail)
      pop_var_detail <- gsub("AN", "allele number", pop_var_detail)
      pop_var_detail <- gsub("AF", "allele frequency", pop_var_detail)
      
      if(pop_var == "none"){
        line_aa_pop <- c("<strong>Known variants in gnomAD (population variants):</strong> none") 
      }else{
        line_aa_pop <- c("<strong>Known variants in gnomAD (population variants):</strong> ",pop_var_detail)
      }
      if(path_var == "none"){
        line_aa_path <- c("<strong>Known variants in ClinVar/HGMD (pathogenic variants):</strong> none") 
      }else{
        line_aa_path <- c("<strong>Known variants in ClinVar/HGMD (pathogenic variants):</strong> ",path_var_detail)
      }
      line_aa_path_str <- paste(line_aa_path, collapse = "")
      line_aa_pop_str <- paste(line_aa_pop, collapse = "")
      
      #SS3
      ss3_reportF <- aa_wise_feature$`3-class DSSP Secondary Structure Properties`
      #SS8
      ss8_reportF <- aa_wise_feature$`8-class DSSP Secondary Structure Properties`
      if(ss8_reportF == "-"){
        line1_1 <- c(aa_mention_str, " doesn't have coordinates in protein structure. ")
        ss3_value = 0;
        ss8_value = 0
      }else{
        if(grepl("beta", ss3_reportF)){
          ss3_value = 1
          ss3_reportFstr <- "&beta;-strand/sheet type"
        }else if(grepl("Heli", ss3_reportF)){
          ss3_value = 2
          ss3_reportFstr <- "helix type"
        }else if(grepl("coil", ss3_reportF)){
          ss3_value = 3
          ss3_reportFstr <- "coil type"
        }
        
        if(grepl("strand", ss8_reportF)){
          ss8_value = 1
          ss8_reportFstr <- "a &beta;-strand"
        }else if(grepl("sheet", ss8_reportF)){
          ss8_value = 2
          ss8_reportFstr <- "a &beta;-sheet"
        }else if(grepl("310", ss8_reportF)){
          ss8_value = 3
          ss8_reportFstr <- "a 3<sub>10</sub>-helix"
        }else if(grepl("alpha", ss8_reportF)){
          ss8_value = 4
          ss8_reportFstr <- "an &alpha;-helix"
        }else if(grepl("pi", ss8_reportF)){
          ss8_value = 5
          ss8_reportFstr <- "a &pi;-helix"
        }else if(grepl("loop", ss8_reportF)){
          ss8_value = 6
          ss8_reportFstr <- "a coil-loop"
        }else if(grepl("bend", ss8_reportF)){
          ss8_value = 7
          ss8_reportFstr <- "a coil-bend"
        }else if(grepl("turn", ss8_reportF)){
          ss8_value = 8
          ss8_reportFstr <- "a coil-turn"
        }
        
        line1_1 <- c(aa_mention_str, " forms ", ss3_reportFstr, " secondary structure, specifically, ", ss8_reportFstr, " in 8-class secondary structure classification. ")
      }
      #ASA
      exposure_reportF <- aa_wise_feature$`Exposure`
      acc_reportF <- aa_wise_feature$`Accessible Surface Area`
      rsa_reportF <- aa_wise_feature$`Relative Accessible Surface Area`
      if(exposure_reportF == "-"){
        asa_value = -1
        line2_1 <- c(aa_mention_str, " doesn't have coordinates in protein structure. ")
      }else{
        rsa_reportF_num <- as.numeric(rsa_reportF)
        if(rsa_reportF_num < 0.05){
          exposure_reportFstr <- "core"
          asa_value = 1
          rsaReport <- "(<5% relative solvent accessibility)"
        }else if((rsa_reportF_num >= 0.05) && (rsa_reportF_num < 0.25)){
          exposure_reportFstr <- "buried"
          asa_value = 2
          rsaReport <- "(5 - 25% relative solvent accessibility)"
        }else if((rsa_reportF_num >= 0.25) && (rsa_reportF_num < 0.50)){
          exposure_reportFstr <- "medium-buried"
          asa_value = 3
          rsaReport <- "(25 - 50% relative accessible surface area)"
        }else if((rsa_reportF_num >= 0.50) && (rsa_reportF_num < 0.75)){
          exposure_reportFstr <- "medium-exposed"
          asa_value = 4
          rsaReport <- "(50 - 75% relative solvent accessibility)"
        }else if(rsa_reportF_num > 0.75){
          exposure_reportFstr <- "exposed"
          asa_value = 5
          rsaReport <- "(>75% relative solvent accessibility)"
        }
        line2_1 <- c(aa_mention_str, " is a ", exposure_reportFstr, " residue ", rsaReport, " in 3D, and ")
      }
      #CP8
      chemp_reportF <- aa_wise_feature$`Physiochemical Property`
      if(thisaacomm != "na"){
        line3_1 <- c(" ", thisaafull, " is a ", thisaaprp, " amino acid ", thisaacomm, ". ")
      }else{
        line3_1 <- c(" ", thisaafull, " is a ", thisaaprp, " amino acid. ")
      }
      thisaaprp_list <- unlist(strsplit(thisaaprp, "[&]"))
      
      #IB4
      bond_reportF <- aa_wise_feature$`Bond Detail`
      bondcollect <- c("")
      if(bond_reportF == "-"){
        bondcollect <- c(bondcollect, "none")
        bondstr <- paste(bondcollect, collapse=", ")
        bondstr <- substring(bondstr, 3)
        thisaaIB_list <- unlist(strsplit(bondstr, "[,]"))
        line4_1 <- c(aa_mention_str, " is not annotated with any protein-protein interaction type.")
      }else{
        if(grepl("-D/", as.character(bond_reportF))){
          bondcollect <- c(bondcollect, "disulfide bond")
        }
        if(grepl("-N/", as.character(bond_reportF))){
          bondcollect <- c(bondcollect, "nonbonded van der Waals")
        }
        if(grepl("-H/", as.character(bond_reportF))){
          bondcollect <- c(bondcollect, "hydrogen bond")
        }
        if(grepl("-S/", as.character(bond_reportF))){
          bondcollect <- c(bondcollect, "salt-bridge ionic")
        }
        
        bondstr <- paste(bondcollect, collapse=", ")
        bondstr <- substring(bondstr, 3)
        line4_1 <- c(aa_mention_str, " is a known ", bondstr, " interaction site(s) in solved protein complex structure.")
        thisaaIB_list <- unlist(strsplit(bondstr, "[,]"))
      }
      
      #PTM6
      ptmcollect = c("")
      acetyl_f <- aa_wise_feature$`Acetylation Detail` 
      if(acetyl_f != "-"){
        #paste("distance to acetylation site(s)", acetyl_f,sep = ": ")
        acetyl_info_all <- unlist(strsplit(acetyl_f, "[;]"))
        for(i in 1:length(acetyl_info_all)){
          acetyl_info <- unlist(strsplit(acetyl_info_all[i], "[/]"))
          if(as.numeric(acetyl_info[2]) < 10.0){
            ptmcollect <- c(ptmcollect, "acetylation")
            break;
          }
        }
      }
      methyl_f <- aa_wise_feature$`Methylation Detail` 
      if(methyl_f != "-"){
        #paste("distance to methylation site(s)", methyl_f,sep = ": ")
        methyl_info_all <- unlist(strsplit(methyl_f, "[;]"))
        for(i in 1:length(methyl_info_all)){
          methyl_info <- unlist(strsplit(methyl_info_all[i], "[/]"))
          if(as.numeric(methyl_info[2]) < 10.0){
            ptmcollect <- c(ptmcollect, "methylation")
            break;
          }
        }
      }
      gcln_f <- aa_wise_feature$`O.GclNAc Detail` 
      if(gcln_f != "-"){
        #paste("distance to O.GclNAc site(s)", gcln_f,sep = ": ")
        gcln_info_all <- unlist(strsplit(gcln_f, "[;]"))
        for(i in 1:length(gcln_info_all)){
          gcln_info <- unlist(strsplit(gcln_info_all[i], "[/]"))
          if(as.numeric(gcln_info[2]) < 10.0){
            ptmcollect <- c(ptmcollect, "O.GclNAc")
            break;
          }
        }
      }
      phos_f <- aa_wise_feature$`Phosphorylation Detail` 
      if(phos_f != "-"){
        #paste("distance to phosphorylation site(s)", phos_f,sep = ": ")
        phos_info_all <- unlist(strsplit(phos_f, "[;]"))
        for(i in 1:length(phos_info_all)){
          phos_info <- unlist(strsplit(phos_info_all[i], "[/]"))
          if(as.numeric(phos_info[2]) < 10.0){
            ptmcollect <- c(ptmcollect, "phosphorylation")
            break;
          }
        }
      }
      sumoy_f <- aa_wise_feature$`Sumoylation Detail` 
      if(sumoy_f != "-"){
        #paste("distance to sumoylation sites", sumoy_f,sep = ": ")
        sumoy_info_all <- unlist(strsplit(sumoy_f, "[;]"))
        for(i in 1:length(sumoy_info_all)){
          sumoy_info <- unlist(strsplit(sumoy_info_all[i], "[/]"))
          if(as.numeric(sumoy_info[2]) < 10.0){
            ptmcollect <- c(ptmcollect, "sumoylation")
            break;
          }
        }
      }
      ubiq_f <- aa_wise_feature$`Ubiquitination Detail` 
      if(ubiq_f != "-"){
        #paste("distance to ubiquitination sites", ubiq_f,sep = ": ")
        ubiq_info_all <- unlist(strsplit(ubiq_f, "[;]"))
        for(i in 1:length(ubiq_info_all)){
          ubiq_info <- unlist(strsplit(ubiq_info_all[i], "[/]"))
          if(as.numeric(ubiq_info[2]) < 10.0){
            ptmcollect <- c(ptmcollect, "ubiquitination")
            break;
          }
        }
      }
      
      if(length(ptmcollect) == 1){
        ptmcollect <- c(ptmcollect, "none")
        ptmcollectstr <- paste(ptmcollect, collapse=", ")
        ptmcollectstr <- substring(ptmcollectstr, 3)
        thisaaPTM_list <- unlist(strsplit(ptmcollectstr, "[,]"))
        line5_1 <- c(aa_mention_str, " doesn't stay at close proximity to any post-translational modification site on protein structure.")
      }else{
        ptmcollectstr <- paste(ptmcollect, collapse=", ")
        ptmcollectstr <- substring(ptmcollectstr, 3)
        thisaaPTM_list <- unlist(strsplit(ptmcollectstr, "[,]"))
        line5_1 <- c(aa_mention_str, " stays at close proximity (less than 10 &#8491; distance in 3D) to ", ptmcollectstr, " site(s) on protein structure.")
      }
      
      #functional site
      upcollect = c("")
      upcollect_nick = c("")
      upcollect_lab = c()
      func_site <- aa_wise_feature$`Functional Site`
      if(func_site != "-"){
        upcollect_this <- paste(c("a known ", as.character(func_site), " (functional site)"), collapse = "")
        upcollect <- c(upcollect, upcollect_this)
        upcollect_nick <- c(upcollect_nick, "sites")
        upcollect_lab <- c(upcollect_lab, 1)
      }else{
        upcollect <- c(upcollect, "")
      }
      bind_region <- aa_wise_feature$`Functional Binding Region`
      if(bind_region != "-"){
        upcollect_nick <- c(upcollect_nick, "binding")
        upcollect_this <- paste(c("part of a ", as.character(bind_region), " region (functional/binding region)"), collapse = "")
        upcollect <- c(upcollect, upcollect_this)
        upcollect_lab <- c(upcollect_lab, 2)
      }else{
        upcollect <- c(upcollect, "")
      }
      seq_motif <- aa_wise_feature$`Sequence Motif Region`
      if(seq_motif != "-"){
        upcollect_nick <- c(upcollect_nick, "motifs")
        upcollect_this <- paste(c("part of a ", as.character(seq_motif)), collapse = "")
        upcollect <- c(upcollect, upcollect_this)
        upcollect_lab <- c(upcollect_lab, 3)
      }else{
        upcollect <- c(upcollect, "")
      }
      domain <- aa_wise_feature$`Modular Domain`
      if(domain != "-"){
        upcollect_nick <- c(upcollect_nick, "modular")
        upcollect_this <- paste(c("part of a ", as.character(domain), " (modular domain)"), collapse = "")
        upcollect <- c(upcollect, upcollect_this)
        upcollect_lab <- c(upcollect_lab, 4)
      }else{
        upcollect <- c(upcollect, "")
      }
      mol_process <- aa_wise_feature$`Molecular Processing`
      if(mol_process != "-"){
        upcollect_nick <- c(upcollect_nick, "molecular")
        upcollect_this <- paste(c("part of ", as.character(mol_process), " region (molecular processing associated region)"), collapse = "")
        upcollect <- c(upcollect, upcollect_this)
        upcollect_lab <- c(upcollect_lab, 5)
      }else{
        upcollect <- c(upcollect, "")
      }
      mod_res <- aa_wise_feature$`Modified Residues`
      if(mod_res != "-"){
        upcollect_nick <- c(upcollect_nick, "modified")
        mod_res_split <- unlist(strsplit(mod_res, split="_"))
        mod_res_split_tog <- paste(mod_res_split, collapse = " ")
        upcollect_this <- paste(c("annotated as a ", as.character(mod_res_split_tog), " site (modified residues)"), collapse = "")
        upcollect <- c(upcollect, upcollect_this)
        upcollect_lab <- c(upcollect_lab, 6)
      }else{
        upcollect <- c(upcollect, "")
      }
      
      upcollect <- upcollect[upcollect != ""]
      #print(upcollect)
      if(length(upcollect_lab) > 1){
        upcollect_multiple <- paste(upcollect, collapse = ", ")
        upcollect_multiple <- sub(",([^,]*)$", " and\\1", upcollect_multiple)
        upcollect <- c(aa_mention_str, " is ", upcollect_multiple, ".")
        line6_1 <- paste(upcollect, collapse="")
        
      }else if((length(upcollect_lab) <= 1) && (length(upcollect_lab) > 0)){
        upcollect <- c(aa_mention_str, " is ", upcollect, ".")
        line6_1 <- paste(upcollect, collapse="")
      }else{
        line6_1 = ""
      }
      
      
      if(line6_1 == ""){
        upcollect_nick <- c(upcollect_nick, "none")
        upcollect_nickstr <- paste(upcollect_nick, collapse=", ")
        upcollect_nickstr <- substring(upcollect_nickstr, 3)
        thisaaUP_list <- unlist(strsplit(upcollect_nickstr, "[,]"))
        line6_1 <- c(" This amino acid is not annotated with any UniProt-based functional feature. ")
      }else{
        upcollect_nickstr <- paste(upcollect_nick, collapse=", ")
        upcollect_nickstr <- substring(upcollect_nickstr, 3)
        thisaaUP_list <- unlist(strsplit(upcollect_nickstr, "[,]"))
        line6_1 <- upcollect
      }
      
      
      # summary report - pclass#
      ## new
      gene_protein_trans_class <- subset(GPTLPC_list, geneName == geneNsummary)
      pclassIDreport = gene_protein_trans_class$pclass
      
      pclassIDreport = gene_protein_trans_class$pclass #pclass
      pclassIds = unlist(strsplit(pclassIDreport, ", "))
      pclassIds = c("c0", pclassIds)
      
      PC24_PF40_path <- subset(PC24_PF40_pathogenic, pcindex %in% pclassIds)
      PC24_PF40_pop <- subset(PC24_PF40_population, pcindex %in% pclassIds)
      
      npclasses = nrow(PC24_PF40_path)
      ## new
      
      # pclass wise report
      
      summ = c("")
      line1 = c("")
      for(i in 2:npclasses){
        pcname <- PC24_PF40_path$pcname[i]
        pcnamelower <- PC24_PF40_path$pcnamelower[i]
        
        summ0 <- c("<br><strong>Detailed Report for ", aa_mention_str, " as ", pcname, ": </strong>")
        summ0str <- paste(summ0, collapse="")
        
        #pathogenic feature for pclass -- start#
        pss3path_df <- PC24_PF40_path[i,4:6]
        pss3path_list <- c()
        for(fn in 1:ncol(pss3path_df))
        {
          pfindex = as.character(names(pss3path_df[1,fn]))
          fflag = as.integer(pss3path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss3path_list <- c(pss3path_list, fdes)
          }
        }
        if(length(pss3path_list) > 0){
          pss3path_str <- paste(pss3path_list, collapse = ", ")
          pss3path <- sub(",([^,]*)$", " and\\1", pss3path_str)
        }else{
          pss3path <- "none"
        }
        
        pss8path_df <- PC24_PF40_path[i,7:14]
        pss8path_list <- c()
        for(fn in 1:ncol(pss8path_df))
        {
          pfindex = as.character(names(pss8path_df[1,fn]))
          fflag = as.integer(pss8path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss8path_list <- c(pss8path_list, fdes)
          }
        }
        if(length(pss8path_list) > 0){
          pss8path_str <- paste(pss8path_list, collapse = ", ")
          pss8path <- sub(",([^,]*)$", " and\\1", pss8path_str)
        }else{
          pss8path <- "none"
        }
        
        pasapath_df <- PC24_PF40_path[i,15:19]
        pasapath_list <- c()
        for(fn in 1:ncol(pasapath_df))
        {
          pfindex = as.character(names(pasapath_df[1,fn]))
          fflag = as.integer(pasapath_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pasapath_list <- c(pasapath_list, fdes)
          }
        }
        if(length(pasapath_list) > 0){
          pasapath_str <- paste(pasapath_list, collapse = ", ")
          pasapath <- sub(",([^,]*)$", " and\\1", pasapath_str)
        }else{
          pasapath <- "none"
        }
        
        pcp8path_df <- PC24_PF40_path[i,20:27]
        pcp8path_list <- c()
        for(fn in 1:ncol(pcp8path_df))
        {
          pfindex = as.character(names(pcp8path_df[1,fn]))
          fflag = as.integer(pcp8path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pcp8path_list <- c(pcp8path_list, fdes)
          }
        }
        if(length(pcp8path_list) > 0){
          pcp8path_str <- paste(pcp8path_list, collapse = ", ")
          pcp8path <- sub(",([^,]*)$", " and\\1", pcp8path_str)
        }else{
          pcp8path <- "none"
        }
        
        #IB - path
        pib4path_df <- PC24_PF40_path[i,28:31]
        pib4path_list <- c()
        for(fn in 1:ncol(pib4path_df))
        {
          pfindex = as.character(names(pib4path_df[1,fn]))
          fflag = as.integer(pib4path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pib4path_list <- c(pib4path_list, fdes)
          }
        }
        if(length(pib4path_list) > 0){
          pib4path_str <- paste(pib4path_list, collapse = ", ")
          pib4path <- sub(",([^,]*)$", " and\\1", pib4path_str)
        }else{
          pib4path <- "none"
        }
        
        #PTM - pop
        pptm6path_df <- PC24_PF40_path[i,32:37]
        pptm6path_list <- c()
        for(fn in 1:ncol(pptm6path_df))
        {
          pfindex = as.character(names(pptm6path_df[1,fn]))
          fflag = as.integer(pptm6path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pptm6path_list <- c(pptm6path_list, fdes)
          }
        }
        if(length(pptm6path_list) > 0){
          pptm6path_str <- paste(pptm6path_list, collapse = ", ")
          pptm6path <- sub(",([^,]*)$", " and\\1", pptm6path_str)
        }else{
          pptm6path <- "none"
        }
        
        #UP - path
        pup6path_df <- PC24_PF40_path[i,38:43]
        pup6path_list <- c()
        for(fn in 1:ncol(pup6path_df))
        {
          pfindex = as.character(names(pup6path_df[1,fn]))
          fflag = as.integer(pup6path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pup6path_list <- c(pup6path_list, fdes)
          }
        }
        if(length(pup6path_list) > 0){
          pup6path_str <- paste(pup6path_list, collapse = ", ")
          pup6path <- sub(",([^,]*)$", " and\\1", pup6path_str)
        }else{
          pup6path <- "none"
        }
        #pathogenic feature for pclass -- end#
        
        #popultaion feature for pclass -- start#
        pss3pop_df <- PC24_PF40_pop[i,4:6]
        pss3pop_list <- c()
        for(fn in 1:ncol(pss3pop_df))
        {
          pfindex = as.character(names(pss3pop_df[1,fn]))
          fflag = as.integer(pss3pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss3pop_list <- c(pss3pop_list, fdes)
          }
        }
        if(length(pss3pop_list) > 0){
          pss3pop_str <- paste(pss3pop_list, collapse = ", ")
          pss3pop <- sub(",([^,]*)$", " and\\1", pss3pop_str)
        }else{
          pss3pop <- "none"
        }
        
        pss8pop_df <- PC24_PF40_pop[i,7:14]
        pss8pop_list <- c()
        for(fn in 1:ncol(pss8pop_df))
        {
          pfindex = as.character(names(pss8pop_df[1,fn]))
          fflag = as.integer(pss8pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss8pop_list <- c(pss8pop_list, fdes)
          }
        }
        if(length(pss8pop_list) > 0){
          pss8pop_str <- paste(pss8pop_list, collapse = ", ")
          pss8pop <- sub(",([^,]*)$", " and\\1", pss8pop_str)
        }else{
          pss8pop <- "none"
        }
        
        pasapop_df <- PC24_PF40_pop[i,15:19]
        pasapop_list <- c()
        for(fn in 1:ncol(pasapop_df))
        {
          pfindex = as.character(names(pasapop_df[1,fn]))
          fflag = as.integer(pasapop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pasapop_list <- c(pasapop_list, fdes)
          }
        }
        if(length(pasapop_list) > 0){
          pasapop_str <- paste(pasapop_list, collapse = ", ")
          pasapop <- sub(",([^,]*)$", " and\\1", pasapop_str)
        }else{
          pasapop <- "none"
        }
        
        pcp8pop_df <- PC24_PF40_pop[i,20:27]
        pcp8pop_list <- c()
        for(fn in 1:ncol(pcp8pop_df))
        {
          pfindex = as.character(names(pcp8pop_df[1,fn]))
          fflag = as.integer(pcp8pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pcp8pop_list <- c(pcp8pop_list, fdes)
          }
        }
        if(length(pcp8pop_list) > 0){
          pcp8pop_str <- paste(pcp8pop_list, collapse = ", ")
          pcp8pop <- sub(",([^,]*)$", " and\\1", pcp8pop_str)
        }else{
          pcp8pop <- "none"
        }
        
        #IB - pop
        pib4pop_df <- PC24_PF40_pop[i,28:31]
        pib4pop_list <- c()
        for(fn in 1:ncol(pib4pop_df))
        {
          pfindex = as.character(names(pib4pop_df[1,fn]))
          fflag = as.integer(pib4pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pib4pop_list <- c(pib4pop_list, fdes)
          }
        }
        if(length(pib4pop_list) > 0){
          pib4pop_str <- paste(pib4pop_list, collapse = ", ")
          #pib4pop_str <- paste(pib4pop_str, "interactions", sep = " ")
          pib4pop <- sub(",([^,]*)$", " and\\1", pib4pop_str)
        }else{
          pib4pop <- "none"
        }
        
        #PTM - pop
        pptm6pop_df <- PC24_PF40_pop[i,32:37]
        pptm6pop_list <- c()
        for(fn in 1:ncol(pptm6pop_df))
        {
          pfindex = as.character(names(pptm6pop_df[1,fn]))
          fflag = as.integer(pptm6pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pptm6pop_list <- c(pptm6pop_list, fdes)
          }
        }
        if(length(pptm6pop_list) > 0){
          pptm6pop_str <- paste(pptm6pop_list, collapse = ", ")
          pptm6pop <- sub(",([^,]*)$", " and\\1", pptm6pop_str)
        }else{
          pptm6pop <- "none"
        }
        
        #UP - pop
        pup6pop_df <- PC24_PF40_pop[i,38:43]
        pup6pop_list <- c()
        for(fn in 1:ncol(pup6pop_df))
        {
          pfindex = as.character(names(pup6pop_df[1,fn]))
          fflag = as.integer(pup6pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pup6pop_list <- c(pup6pop_list, fdes)
          }
        }
        if(length(pup6pop_list) > 0){
          pup6pop_str <- paste(pup6pop_list, collapse = ", ")
          pup6pop <- sub(",([^,]*)$", " and\\1", pup6pop_str)
        }else{
          pup6pop <- "none"
        }
        

        #population feature for pclass -- end#
        
        #pathogenic feature for pclass -- start - sentence construction#
        if(pss3path != "none"){
          if(pss8path != "none"){
            summ1_2 <- c(pcname, "s in general have significant enrichment of pathogenic variants in ", trimws(pss3path), " (out of three classes), and in ", trimws(pss8path), " (out of eight classes). ")
          }else{
            summ1_2 <- c(pcname, "s in general have significant enrichment of pathogenic variants in ", trimws(pss3path), " (out of three classes). ")
          }
        }else{
          if(pss8path != "none"){
            summ1_2 <- c(pcname, "s in general have significant enrichment of pathogenic variants in ", trimws(pss8path), " (out of eight classes). ")
          }else{
            summ1_2 <- c(pcname, "s in general have no significant enrichment of pathogenic variants in protein secondary structure. ")
          }
        }
        #pathogenic feature for pclass -- end - sentence construction#
        
        #population feature for pclass -- start - sentence construction#
        if(pss3pop != "none"){
          if(pss8pop != "none"){
            summ1_3 <- c("And, population variants in ", pcnamelower, "s are significantly enriched on ", trimws(pss3pop), ", and on ", trimws(pss8pop), " in terms of 3-class and 8-class secondary structure classifications, respectively. ")
          }else{
            summ1_3 <- c("And, population variants in ", pcnamelower, "s are significantly enriched on ", trimws(pss3pop), " (out of three classes). ")
          }
        }else{
          if(pss8pop != "none"){
            summ1_3 <- c("And, population variants in ", pcnamelower, "s are significantly enriched on ", trimws(pss8pop), " (out of eight classes). ")
          }else{
            summ1_3 <- c(pcname, "s have no significant enrichment of population variants in protein secondary structure. ")
          }
        }
        #population feature for pclass -- end - sentence construction#
        
        #pathogenic feature for pclass -- start - sentence construction#
        if(pasapath != "none"){
          summ2_2 <- c(pcnamelower, "s carry a significant enrichment of pathogenic variants in ", trimws(pasapath), " residues. ")
        }else{
          summ2_2 <- c(" solvent accessibility is not a significantly associated features of pathogenic variants in ", pcname, "s. ")
        }
        #pathogenic feature for pclass -- end - sentence construction#
        
        #population feature for pclass -- start - sentence construction#
        if(pasapop != "none"){
          summ2_3 <- c("The population variants in ", pcnamelower, "s mostly affect the ", trimws(pasapop), " residues. ")
        }else{
          summ2_3 <- c(" The solvent accessibility is not a significantly associated features of population variants in ", pcname, "s. ")
        }
        #population feature for pclass -- end - sentence construction#
        
        #pathogenic feature for pclass -- start - sentence construction#
        if(pcp8path != "none"){
          summ3_2 <- c("The pathogenic missense variants in ", pcnamelower, "s mutate the ", trimws(pcp8path), " residues significantly more than other types of residues. ")
        }else{
          summ3_2 <- c("Physicochemical properties of amino acid is not a significantly associated features of pathogenic variants in ", pcname, ". ")
        }
        #pathogenic feature for pclass -- end - sentence construction#
        
        #population feature for pclass -- start - sentence construction#
        if(pcp8pop != "none"){
          summ3_3 <- c("In contrast, the population missense variants in ", pcnamelower, "s mutate the ", trimws(pcp8pop), " residues significantly more than other residue types. ")
        }else{
          summ3_3 <- c("Physicochemical properties of amino acid is not a significantly associated features of population variants in ", pcname, ". ")
        }
        #population feature for pclass -- end - sentence construction#
        
        #pathogenic feature for pclass -- start - sentence construction#
        if((pib4path != "none") && (pib4pop != "none")){
          summ4_2 <- c(" Pathogenic missense variants in ", pcnamelower,  "s are significantly more frequent on ", trimws(pib4path), " interaction sites. ")
          summ4_3 <- c(" And, population variants in ", pcnamelower,  "s are significantly more frequent on ", trimws(pib4pop), " interaction sites. ")
        }else if((pib4path != "none") && (pib4pop == "none")){
          summ4_2 <- c(" Pathogenic variants in ", pcnamelower,  "s are significantly more frequent on ", trimws(pib4path), " interaction sites. ")
          summ4_3 <- c(" And, population variants in ", pcnamelower,  "s are not significantly more frequent on any interaction sites. ")
        }else if((pib4path == "none") && (pib4pop != "none")){
          summ4_2 <- c(" Population variants in ", pcnamelower,  "s are significantly more frequent on ", trimws(pib4pop), " interaction sites. ")
          summ4_3 <- c(" And, pathogenic variants in ", pcnamelower,  "s are not significantly more frequent on any interaction sites. ")
        }else if((pib4path == "none") && (pib4pop == "none")){
          summ4_2 <- c(" There is no significant enrichment of pathogenic ")
          summ4_3 <- c(" or population variants on protein-protein interaction sites. ")
        }
        #population feature for pclass -- end - sentence construction#
        
        #pathogenic and population feature for pclass -- start - sentence construction#
        if((pptm6path != "none") && (pptm6pop != "none")){
          summ5_2 <- c(" Pathogenic missense variants in ", pcnamelower,  "s are significantly enriched on/around 10 &#8491; distance to ", trimws(pptm6path), " sites. ")
          summ5_3 <- c(" Population missense variants in ", pcnamelower,  "s are significantly enriched on/around 10 &#8491; distance to ", trimws(pptm6pop), " sites. ")
        }else if((pptm6path != "none") && (pptm6pop == "none")){
          summ5_2 <- c(" Pathogenic missense variants in ", pcnamelower,  "s are significantly enriched on/around 10 &#8491; distance to ", trimws(pptm6path), " sites. ")
          summ5_3 <- c(" In comparison, population missense variants in ", pcnamelower,  "s are not significantly enriched within 10 &#8491; spatial distance to any post-translational modification sites. ")
        }else if((pptm6path == "none") && (pptm6pop != "none")){
          summ5_2 <- c(" Population missense variants in ", pcnamelower,  "s are significantly enriched on/around 10 &#8491; distance to ", trimws(pptm6pop), " sites. ")
          summ5_3 <- c(" However, pathogenic missense variants in ", pcnamelower,  "s are not significantly enriched within 10 &#8491; spatial distance to any post-translational modification sites. ")
        }else if((pptm6path == "none") && (pptm6pop == "none")){
          summ5_2 <- c(" There is no significant enrichment of pathogenic ")
          summ5_3 <- c(" or population variants within 10 &#8491; spatial distance to any post-translational modification sites in ", pcnamelower,  "s. ")
        }
        #pahogenic and population feature for pclass -- end - sentence construction#
        
        #pathogenic feature for pclass -- end - sentence construction#
        if(pup6path != "none"){
          summ6_2 <- c(" ", pcname, "s in general show significant enrichment of pathogenic variants in ", trimws(pup6path), ".")
        }else{
          summ6_2 <- c(" ", pcname, "s in general show no significant enrichment of pathogenic variants in any functional features from UniProt.")
        }
        #pathogenic feature for pclass -- end - sentence construction#
        
        #population feature for pclass -- end - sentence construction#
        if(pup6pop != "none"){
          summ6_3 <- c(" On the other hand, ", pcnamelower, "s have significant enrichment of population variants in ", trimws(pup6pop), ". ")
        }else{
          summ6_3 <- c("and, have no significant enrichment of population variants in any functional features from UniProt.")
        }
        #population feature for pclass -- end - sentence construction#
        
        summ <- c(summ, summ0, line1_1, summ1_2, summ1_3, line2_1, summ2_2, summ2_3, line3_1, summ3_2, summ3_3, line4_1, summ4_2, summ4_3, line5_1, summ5_2, summ5_3, line6_1, summ6_2, summ6_3, "<br>")
      }
      #pathogenic feature for pclass -- end - sentence construction
      
      summstr <- paste(summ, collapse="")
      
      # pclass wise report
      
      HTML(
        paste(
          "<br>",
          "<b>","<font size = +1 color = Black face = Arial>","&nbsp;&nbsp;","</font>","</b>",summstr, "<br>"
        )
      )
    }
  })
  
  output$reportTsummreportPlot <- renderPlot({
    if((input$aa_Selected != "" && input$reportTgeneSelected != "")){  
      tgene <- input$reportTgeneSelected
      taa <- input$aa_Selected
      
      tgene_info <- GPTLPC_list[GPTLPC_list$geneName == tgene, ]
      tgene_pcs <- trimws(c("c0",unlist(strsplit(tgene_info$pclass, split = ","))))
      
      
      tgene_pcs_pathpop_stat <- data.frame(gene=as.character(), unipos=as.character(), unires=as.character(), match_count=as.character(), var_type=as.character(), cindex=as.character())
      aa_match_path <- c("")  
      aa_match_pop <- c("")
      for(i in 1:length(tgene_pcs)){
        pcindex <- tgene_pcs[i]
        pcname_row <- PC24_PF40_population[PC24_PF40_population$pcindex == pcindex, ]
        class_name_mod <- paste(c(LETTERS[i], ") ", pcname_row$pcname), collapse="")
        
        pc_pathstat_filename = paste(c("25/", pcindex,"_path.txt"),collapse ='')
        pc_pathstat <- read_delim(pc_pathstat_filename, "\t", escape_double = FALSE, trim_ws = TRUE)
        pc_popstat_filename = paste(c("25/", pcindex,"_pop.txt"),collapse ='')
        pc_popstat <- read_delim(pc_popstat_filename, "\t", escape_double = FALSE, trim_ws = TRUE)
        
        tgene_pc_pathpop_stat_filename <- paste(c("1330_pc/", tgene, "_fmatch_", pcindex, ".txt"),collapse ='')
        tgene_pc_pathpop_stat <- read_delim(tgene_pc_pathpop_stat_filename, "\t", escape_double = FALSE, trim_ws = TRUE)
        
        aa_pc_match_row <- tgene_pc_pathpop_stat[tgene_pc_pathpop_stat$aapos == taa, ]
        aa_pc_match_path <- as.numeric(aa_pc_match_row$match_path)
        aa_match_path <- c(aa_match_path, aa_pc_match_path)
        aa_pc_match_pop <- as.numeric(aa_pc_match_row$match_pop)
        aa_match_pop <- c(aa_match_pop, aa_pc_match_pop)
        
        pc_pathstat$var <- "Pathogenic"
        pc_popstat$var <- "Population"
        
        colnames(pc_pathstat) <- c("gene", "unipos", "unires", "match_count", "var_type")
        colnames(pc_popstat) <- c("gene", "unipos", "unires", "match_count", "var_type")
        
        pc_pathpop_stat <- rbind(pc_pathstat, pc_popstat)
        pc_pathpop_stat$cindex <- pcindex
        pc_pathpop_stat$cname <- class_name_mod
        
        tgene_pcs_pathpop_stat <- rbind(tgene_pcs_pathpop_stat, pc_pathpop_stat)
      }
      
      
      facet_header_text <- as.vector(unique(tgene_pcs_pathpop_stat$cname))
      
      path_ann <- data.frame(
        cname   = facet_header_text,
        x     = c(rep(1, length(tgene_pcs))),
        y     = as.integer(aa_match_path[2:(length(tgene_pcs)+1)])
      )
      pop_ann <- data.frame(
        cname   = facet_header_text,
        x     = c(rep(2, length(tgene_pcs))),
        y     = as.integer(aa_match_pop[2:(length(tgene_pcs)+1)])
      )
      dodge <- position_dodge(width = T)
      
      ggplot(data = tgene_pcs_pathpop_stat, aes(tgene_pcs_pathpop_stat$var_type, tgene_pcs_pathpop_stat$match_count)) +
        geom_boxplot(aes(fill=factor(tgene_pcs_pathpop_stat$var_type)),width=.3, position = dodge, alpha = 0.2) +
        geom_point(data = path_ann, mapping = aes(x = path_ann$x, y = path_ann$y), colour = "red", size = 3) +
        geom_point(data = pop_ann, mapping = aes(x = pop_ann$x, y = pop_ann$y), colour = "blue", size = 3) +
        theme_bw(base_size = 16) +
        theme(legend.position = "top", legend.direction = "horizontal") +
        labs(title="Distribution of number of characteristic 3D features of amino acids affacted by pathogenic and population \n variants (for all proteins together and individually for protein classes)", 
             subtitle = "The red and blue circles show the number of pathogenic and population variant-associated 3D features of the selected amino acid, respectively",
             x="Variant types",
             y="",
             fill="Variant Type",
             caption = "") +
        scale_fill_manual(name="", 
                          values = c("Population"="blue", "Pathogenic"="red")) +
        facet_wrap(. ~cname, scales = "free_y") 
    }
  })
  
  report_summReportHeadlines <- eventReactive({
    input$aa_Selected
    input$reportTsubmit
  },{
    print(input$reportTgeneSelected)
    #reportGeneName <- input$reportTgeneSelected  
    if(input$aa_Selected != ""){  
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      report_gene_nameF = paste("gene_wise_info/",input$reportTgeneSelected,".txt",sep='')
      report_gene_wise_infoF <- read_delim(report_gene_nameF, "\t", escape_double = FALSE, trim_ws = TRUE)
      geneNsummary <- input$reportTgeneSelected
      
      aa_wise_feature <- subset(report_gene_wise_infoF, `Amino Acid Index` == input$aa_Selected)
      #AA
      amino_acid_reportF <- aa_wise_feature$`Amino Acid`
      
      thisaa1 <- amino_acid_reportF
      thisaaindex <- input$aa_Selected
      this_aa_info <- subset(aa_info, `one` == thisaa1)
      thisaa3 <- this_aa_info$three
      thisaafull <- this_aa_info$full
      thisaaprp <- this_aa_info$prp
      thisaacomm <- this_aa_info$comm
      
      aa_mention <- c(thisaafull, " ", thisaaindex, " of <em>", geneNsummary, "</em>")
      aa_mention_str <- paste(aa_mention, collapse = "")
      
      line_aa <- c("<strong>Amino Acid Residue:</strong> ", thisaafull, " (", thisaa3, "/", thisaa1, "), position: ", thisaaindex)
      line_aa_str <- paste(line_aa, collapse = "")
      
      path_var <- aa_wise_feature$`Pathogenic Mutation`
      pop_var <- aa_wise_feature$`Population Mutation`
      
      path_var_detail <- aa_wise_feature$`Pathogenic Mutation (details)`
      pop_var_detail <- aa_wise_feature$`Population Mutation (details)`
      pop_var_detail <- gsub("AC", "allele count", pop_var_detail)
      pop_var_detail <- gsub("AN", "allele number", pop_var_detail)
      pop_var_detail <- gsub("AF", "allele frequency", pop_var_detail)
      
      if(pop_var == "none"){
        line_aa_pop <- c("<strong>Known variants in gnomAD (population variants):</strong> none") 
      }else{
        line_aa_pop <- c("<strong>Known variants in gnomAD (population variants):</strong> ", pop_var_detail)
      }
      if(path_var == "none"){
        line_aa_path <- c("<strong>Known variants in ClinVar/HGMD (pathogenic variants):</strong> none") 
      }else{
        line_aa_path <- c("<strong>Known variants in ClinVar/HGMD (pathogenic variants):</strong> ", path_var_detail)
      }
      line_aa_path_str <- paste(line_aa_path, collapse = "")
      line_aa_pop_str <- paste(line_aa_pop, collapse = "")
      
      # summary report - pclass#
      
      tgene <- input$reportTgeneSelected
      taa <- input$aa_Selected
      
      tgene_info <- GPTLPC_list[GPTLPC_list$geneName == tgene, ]
      tgene_pcs <- trimws(c("c0",unlist(strsplit(tgene_info$pclass, split = ","))))

      
      # pclass wise report
      line0 <- c("<em>", geneNsummary, "</em> is grouped into")
      
      npclasses = length(tgene_pcs)
      
      summ = c("")
      line1 = c("")
      for(i in 1:npclasses){
        pcindex <- tgene_pcs[i]
        pcname_row <- PC24_PF40_population[PC24_PF40_population$pcindex == pcindex, ]
        class_name_mod <- paste(c(LETTERS[i], ") ", pcname_row$pcname), collapse="")
      
        pcname <- pcname_row$pcname
        pcnamelower <- pcname_row$pcnamelower
        
        tgene_pc_pathpop_stat_filename <- paste(c("1330_pc/", tgene, "_fmatch_", pcindex, ".txt"),collapse ='')
        tgene_pc_pathpop_stat <- read_delim(tgene_pc_pathpop_stat_filename, "\t", escape_double = FALSE, trim_ws = TRUE)
        
        aa_pc_match_row <- tgene_pc_pathpop_stat[tgene_pc_pathpop_stat$aapos == taa, ]
        
        aa_pc_match_path <- as.numeric(aa_pc_match_row$match_path)
        pclasspathfeature <- as.numeric(aa_pc_match_row$total_path)
        pclasspathfeature_mean_match <- as.numeric(aa_pc_match_row$mean_match_path)
        
        aa_pc_match_pop <- as.numeric(aa_pc_match_row$match_pop)
        pclasspopfeature <- as.numeric(aa_pc_match_row$total_pop)
        pclasspopfeature_mean_match <- as.numeric(aa_pc_match_row$mean_match_pop)
        
        this_amino_acid_path = aa_pc_match_path;
        this_amino_acid_pop = aa_pc_match_pop;
        
        if(i == 1){
          line0 <- c(line0, " ", "<strong>", LETTERS[i], ") ", pcnamelower, "</strong>")
        }else if(i == npclasses){
          line0 <- c(line0, " and ", "<strong>", LETTERS[i], ") ", pcnamelower, "</strong>")
        }else{
          line0 <- c(line0, ", ", "<strong>", LETTERS[i], ") ", pcnamelower, "</strong>")
        }
        
        
        if(i == 1){
          line1 <- c(line1, "<br></br>", "<strong>", LETTERS[i], ") ", pcname, "</strong>", "<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", "Collectively for ", pcnamelower, ", we identified ", "<strong>", pclasspathfeature, "</strong>", " features that are significantly associated with pathogenic variants and ", "<strong>", pclasspopfeature, "</strong>", " features that are significantly associated with population variants. ")
          if(this_amino_acid_path > pclasspathfeature_mean_match){
            line1 <- c(line1, aa_mention_str, " possesses ", "<strong>", this_amino_acid_path, " out of ", pclasspathfeature, "</strong>", " pathogenic variant-associated features (> average match obtained by the known pathogenic variants with 3D-coordinates in ", pcnamelower, ", see the plots below) and ")
          }else{
            line1 <- c(line1, aa_mention_str, " possesses ", "<strong>", this_amino_acid_path, " out of ", pclasspathfeature, "</strong>", " pathogenic variant-associated features (< average match obtained by the known pathogenic variants with 3D-coordinates in ", pcnamelower, ", see the plots below) and ")
          }
          if(this_amino_acid_pop > pclasspopfeature_mean_match){
            line1 <- c(line1, "<strong>", this_amino_acid_pop, " out of ", pclasspopfeature, "</strong>", " population variant-associated features (> average match obtained by the known population variants with 3D-coordinates in ", pcnamelower, ".")
          }else{
            line1 <- c(line1, "<strong>", this_amino_acid_pop, " out of ", pclasspopfeature, "</strong>", " population variant-associated features (< average match obtained by the known population variants with 3D-coordinates in ", pcnamelower, ".")
          }
        }else{
          line1 <- c(line1, "<br></br>", "<strong>", LETTERS[i], ") ", pcname, "</strong>", "<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", "For ", pcnamelower, "s in general, we identified ", "<strong>", pclasspathfeature, "</strong>", " features that are significantly associated with pathogenic variants and ", "<strong>", pclasspopfeature, "</strong>", " features that are significantly associated with population variants. ")
          if(this_amino_acid_path > pclasspathfeature_mean_match){
            line1 <- c(line1, aa_mention_str, " possesses ", "<strong>", this_amino_acid_path, " out of ", pclasspathfeature, "</strong>", " pathogenic variant-associated features (> average match obtained by the known pathogenic variants with 3D-coordinates in ", pcnamelower, "s, see the plots below) and ")
          }else{
            line1 <- c(line1, aa_mention_str, " possesses ", "<strong>", this_amino_acid_path, " out of ", pclasspathfeature, "</strong>", " pathogenic variant-associated features (< average match obtained by the known pathogenic variants with 3D-coordinates in ", pcnamelower, "s, see the plots below) and ")
          }
          if(this_amino_acid_pop > pclasspopfeature_mean_match){
            line1 <- c(line1, "<strong>", this_amino_acid_pop, " out of ", pclasspopfeature, "</strong>", " population variant-associated features (> average match obtained by the known population variants with 3D-coordinates in ", pcnamelower, ".")
          }else{
            line1 <- c(line1, "<strong>", this_amino_acid_pop, " out of ", pclasspopfeature, "</strong>", " population variant-associated features (< average match obtained by the known population variants with 3D-coordinates in ", pcnamelower, ".")
          }
        }
        line1 <- c(line1, " Pathogenic 3D Feature index, ", "<strong><em>", "P3DFi", "</em></strong>", " = pathogenic variant-associated features -  population variant-associated features = ", this_amino_acid_path, " - ", this_amino_acid_pop, " = <strong>", this_amino_acid_path - this_amino_acid_pop, "</strong>.")
        
      }
      #pathogenic feature for pclass -- end - sentence construction
      
      line0 <- c(line0, ".", " We evaluated the burden of pathogenic and population variants in a total of 40 protein structural, physicochemical, and functional features (see Documentation for feature definition and ascertainment). ")
      line0str <- paste(line0, collapse="")
      
      
      line1str <- paste(line1, collapse="")
      
      # pclass wise report
      
      HTML(
        paste(
          "<b>","<font size = +1 color = Black face = Arial>","<strong>I) Amino acid information: ","</strong>","</font>","</b>", "<br>",
          "<b>","<font size = +1 color = Black face = Arial>","</font>","</b>",line_aa_str, "<br>", line_aa_pop_str, "<br>", line_aa_path_str,"<br></br>",
          "<b>","<font size = +1 color = Black face = Arial>","<strong>II) Protein class-specific summary Report for ",aa_mention_str,":</strong>","</font>","</b>", "<br>", line0str,
          "<b>","<font size = +1 color = Black face = Arial>","</font>","</b>","&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", line1str,"<br>"
        )
      )
    }
  })
  
  output$reportTinformation = renderText(
    {
      report_info()
    }
  )
  
  output$reportTaafeature = renderText(
    {
      if(input$aa_Selected_aaFeat != "")
      {
        report_aaFeature()
      }
    }
  )
  
  output$reportTsummreport = renderText(
    {
      if((input$reportTgeneSelected != "" && input$aa_Selected != "" ))
      {
        report_summReport()
      }
    }
  )
  
  output$reportTsummreportHeadlines = renderText(
    {
      if((input$reportTgeneSelected != "" && input$aa_Selected != ""))
      {
        report_summReportHeadlines()
      }
    }
  )
  
  details_information2 <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      
      gene_name = input$geneSelected #input gene
      gene_protein_trans_class <- subset(GPTLPC_list, geneName == gene_name) #crossref info for that gene
      uniprotIDresearch = as.character(gene_protein_trans_class$uniprotAc) #uniprot
      uniprotac_link = c("https://www.uniprot.org/uniprot/", uniprotIDresearch)
      uniprotac_link_str = paste(uniprotac_link, collapse = "")
      protein_name = as.character(gene_protein_trans_class$Protein_Name) #protein name
      
      pclassIDreport = gene_protein_trans_class$pclass #pclass
      pclassIds = unlist(strsplit(pclassIDreport, ", "))
      pclassIds = c("c0", pclassIds)
      
      PC24_PF40_path <- subset(PC24_PF40_pathogenic, pcindex %in% pclassIds)
      PC24_PF40_pop <- subset(PC24_PF40_population, pcindex %in% pclassIds)
      
      class_list_ele <- length(pclassIds)
      pClassName_df <- PC24_PF40_path$pcname
      pClassName_str = unlist(strsplit(pClassName_df, ", "))
      proteinClass_list_str <- paste(LETTERS[1:class_list_ele], pClassName_str, sep=") ")
      pclass_all_str <- paste(proteinClass_list_str, collapse=", ")
      
      this_pclass <- c("")
      summary_gene <- c("")
      for(i in 1:class_list_ele){
        pcname <- PC24_PF40_path$pcname[i]
        pcnamelower <- PC24_PF40_path$pcnamelower[i]
        
        if(i == 1){
          this_pclass <- c(this_pclass, "<strong>", "<font size = +1 color = Black face = Arial>", LETTERS[i], ") ", "<i>", as.character(gene_name), "</i>", " as part of ", pcname, ":</font></strong> <br>")
        }else{
          this_pclass <- c(this_pclass, "<strong>", "<font size = +1 color = Black face = Arial>", LETTERS[i], ") ", "<i>", as.character(gene_name), "</i>", " as ", pcname, ":</font></strong> <br>")
        }
        
        #pathogenic feature for pclass -- start#
        pss3path_df <- PC24_PF40_path[i,4:6]
        pss3path_list <- c()
        for(fn in 1:ncol(pss3path_df))
        {
          pfindex = as.character(names(pss3path_df[1,fn]))
          fflag = as.integer(pss3path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
            pss3path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
            pss3path_list <- c(pss3path_list, pss3path_for)
          }
        }
        if(length(pss3path_list) > 0){
          pss3path_str <- paste(pss3path_list, collapse = ", ")
        }else{
          pss3path_str <- "none"
        }
        
        
        pss8path_df <- PC24_PF40_path[i,7:14]
        pss8path_list <- c()
        for(fn in 1:ncol(pss8path_df))
        {
          pfindex = as.character(names(pss8path_df[1,fn]))
          fflag = as.integer(pss8path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss8path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>"), collapse = "")
            pss8path_list <- c(pss8path_list, pss8path_for)
          }
        }
        if(length(pss8path_list) > 0){
          pss8path_str <- paste(pss8path_list, collapse = ", ")
        }else{
          pss8path_str <- "none"
        }
        
        pasapath_df <- PC24_PF40_path[i,15:19]
        pasapath_list <- c()
        for(fn in 1:ncol(pasapath_df))
        {
          pfindex = as.character(names(pasapath_df[1,fn]))
          fflag = as.integer(pasapath_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
            pasapath_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
            pasapath_list <- c(pasapath_list, pasapath_for)
          }
        }
        if(length(pasapath_list) > 0){
          pasapath_str <- paste(pasapath_list, collapse = ", ")
        }else{
          pasapath_str <- "none"
        }

        pcp8path_df <- PC24_PF40_path[i,20:27]
        pcp8path_list <- c()
        for(fn in 1:ncol(pcp8path_df))
        {
          pfindex = as.character(names(pcp8path_df[1,fn]))
          fflag = as.integer(pcp8path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
            pcp8path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
            pcp8path_list <- c(pcp8path_list, pcp8path_for)
          }
        }
        if(length(pcp8path_list) > 0){
          pcp8path_str <- paste(pcp8path_list, collapse = ", ")
          pcp8path_str <- paste(pcp8path_str, "amino acids", sep = " ")
        }else{
          pcp8path_str <- "none"
        }
        
        pib4path_df <- PC24_PF40_path[i,28:31]
        pib4path_list <- c()
        for(fn in 1:ncol(pib4path_df))
        {
          pfindex = as.character(names(pib4path_df[1,fn]))
          fflag = as.integer(pib4path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pib4path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>"), collapse = "")
            pib4path_list <- c(pib4path_list, pib4path_for)
          }
        }
        if(length(pib4path_list) > 0){
          pib4path_str <- paste(pib4path_list, collapse = ", ")
          pib4path_str <- paste(pib4path_str, "interactions", sep = " ")
        }else{
          pib4path_str <- "none"
        }
        
        pptm6path_df <- PC24_PF40_path[i,32:37]
        pptm6path_list <- c()
        for(fn in 1:ncol(pptm6path_df))
        {
          pfindex = as.character(names(pptm6path_df[1,fn]))
          fflag = as.integer(pptm6path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pptm6path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>"), collapse = "")
            pptm6path_list <- c(pptm6path_list, pptm6path_for)
          }
        }
        if(length(pptm6path_list) > 0){
          pptm6path_str <- paste(pptm6path_list, collapse = ", ")
        }else{
          pptm6path_str <- "none"
        }
        
        pup6path_df <- PC24_PF40_path[i,38:43]
        pup6path_list <- c()
        for(fn in 1:ncol(pup6path_df))
        {
          pfindex = as.character(names(pup6path_df[1,fn]))
          fflag = as.integer(pup6path_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
            pup6path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
            pup6path_list <- c(pup6path_list, pup6path_for)
          }
        }
        if(length(pup6path_list) > 0){
          pup6path_str <- paste(pup6path_list, collapse = ", ")
        }else{
          pup6path_str <- "none"
        }
        #pathogenic feature for pclass -- end#
        
        #popultaion feature for pclass -- start#
        pss3pop_df <- PC24_PF40_pop[i,4:6]
        pss3pop_list <- c()
        for(fn in 1:ncol(pss3pop_df))
        {
          pfindex = as.character(names(pss3pop_df[1,fn]))
          fflag = as.integer(pss3pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
            pss3pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
            pss3pop_list <- c(pss3pop_list, pss3pop_for)
          }
        }
        if(length(pss3pop_list) > 0){
          pss3pop_str <- paste(pss3pop_list, collapse = ", ")
        }else{
          pss3pop_str <- "none"
        }
        
        pss8pop_df <- PC24_PF40_pop[i,7:14]
        pss8pop_list <- c()
        for(fn in 1:ncol(pss8pop_df))
        {
          pfindex = as.character(names(pss8pop_df[1,fn]))
          fflag = as.integer(pss8pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pss8pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>"), collapse = "")
            pss8pop_list <- c(pss8pop_list, pss8pop_for)
          }
        }
        if(length(pss8pop_list) > 0){
          pss8pop_str <- paste(pss8pop_list, collapse = ", ")
        }else{
          pss8pop_str <- "none"
        }
        
        pasapop_df <- PC24_PF40_pop[i,15:19]
        pasapop_list <- c()
        for(fn in 1:ncol(pasapop_df))
        {
          pfindex = as.character(names(pasapop_df[1,fn]))
          fflag = as.integer(pasapop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
            pasapop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
            pasapop_list <- c(pasapop_list, pasapop_for)
          }
        }
        if(length(pasapop_list) > 0){
          pasapop_str <- paste(pasapop_list, collapse = ", ")
        }else{
          pasapop_str <- "none"
        }
        
        pcp8pop_df <- PC24_PF40_pop[i,20:27]
        pcp8pop_list <- c()
        for(fn in 1:ncol(pcp8pop_df))
        {
          pfindex = as.character(names(pcp8pop_df[1,fn]))
          fflag = as.integer(pcp8pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
            pcp8pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
            pcp8pop_list <- c(pcp8pop_list, pcp8pop_for)
          }
        }
        if(length(pcp8pop_list) > 0){
          pcp8pop_str <- paste(pcp8pop_list, collapse = ", ")
          pcp8pop_str <- paste(pcp8pop_str, "amino acids", sep = " ")
        }else{
          pcp8pop_str <- "none"
        }
        
        pib4pop_df <- PC24_PF40_pop[i,28:31]
        pib4pop_list <- c()
        for(fn in 1:ncol(pib4pop_df))
        {
          pfindex = as.character(names(pib4pop_df[1,fn]))
          fflag = as.integer(pib4pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pib4pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>"), collapse = "")
            pib4pop_list <- c(pib4pop_list, pib4pop_for)
          }
        }
        if(length(pib4pop_list) > 0){
          pib4pop_str <- paste(pib4pop_list, collapse = ", ")
          pib4pop_str <- paste(pib4pop_str, "interactions", sep = " ")
        }else{
          pib4pop_str <- "none"
        }
        
        pptm6pop_df <- PC24_PF40_pop[i,32:37]
        pptm6pop_list <- c()
        for(fn in 1:ncol(pptm6pop_df))
        {
          pfindex = as.character(names(pptm6pop_df[1,fn]))
          fflag = as.integer(pptm6pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            pptm6pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>"), collapse = "")
            pptm6pop_list <- c(pptm6pop_list, pptm6pop_for)
          }
        }
        if(length(pptm6pop_list) > 0){
          pptm6pop_str <- paste(pptm6pop_list, collapse = ", ")
        }else{
          pptm6pop_str <- "none"
        }
        
        pup6pop_df <- PC24_PF40_pop[i,38:43]
        pup6pop_list <- c()
        for(fn in 1:ncol(pup6pop_df))
        {
          pfindex = as.character(names(pup6pop_df[1,fn]))
          fflag = as.integer(pup6pop_df[fn])
          if(fflag == 1){
            fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
            fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
            pup6pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
            pup6pop_list <- c(pup6pop_list, pup6pop_for)
          }
        }
        if(length(pup6pop_list) > 0){
          pup6pop_str <- paste(pup6pop_list, collapse = ", ")
        }else{
          pup6pop_str <- "none"
        }
        #population feature for pclass -- end#
        
        # see Documentation for feature definition and ascertainment
        
        this_pclass <- c(this_pclass, "<font size = -1 color = Black face = Arial>", "<strong>Features with significant enrichment of <span style='color: #993300';>pathogenic variants</span> in ", pcnamelower , "s: </strong>")
        this_pclass <- c(this_pclass, "<ol>")
        this_pclass <- c(this_pclass, "<li>", "3-class secondary structure: ", pss3path_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "8-class secondary structure: ", pss8path_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "Residue exposure level: ", pasapath_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "Physicochemical&nbsp;property: ", pcp8path_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "Protein-protein interaction/bond type: ", pib4path_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "Post-translational modification (within &lt;10 &#8491 spatial distance): ", pptm6path_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "UniProt-based function features: ", pup6path_str, "</li>")
        this_pclass <- c(this_pclass, "</ol>")
        this_pclass <- c(this_pclass, "<strong>Features with significant enrichment of <span style='color: #000080';>population variants</span> in ", pcnamelower , "s: </strong>")
        this_pclass <- c(this_pclass, "<ol>")
        this_pclass <- c(this_pclass, "<li>", "3-class secondary structure: ", pss3pop_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "8-class secondary structure: ", pss8pop_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "Residue exposure level: ", pasapop_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "Physicochemical&nbsp;property: ", pcp8pop_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "Protein-protein interaction/bond type: ", pib4pop_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "Post-translational modification (within &lt;10 &#8491 spatial distance): ", pptm6pop_str, "</li>")
        this_pclass <- c(this_pclass, "<li>", "UniProt-based function features: ", pup6pop_str, "</li>")
        this_pclass <- c(this_pclass, "</ol></font><br>")
      }
      
      summary_gene <- paste(this_pclass, collapse = "")
      
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      gene_info_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      gene_info_gene_name_info <- read_delim(gene_info_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_info <- nrow(gene_info_gene_name_info)
      gene_aa_path_var <- nrow(gene_info_gene_name_info[gene_info_gene_name_info$`Pathogenic Mutation` != "none", ])
      gene_aa_path_var_wstruc <- nrow(gene_info_gene_name_info[(gene_info_gene_name_info$`Pathogenic Mutation` != "none") & (gene_info_gene_name_info$`8-class DSSP Secondary Structure Properties` != "-"), ])
      gene_aa_pop_var <- nrow(gene_info_gene_name_info[gene_info_gene_name_info$`Population Mutation` != "none", ])
      gene_aa_pop_var_wstruc <- nrow(gene_info_gene_name_info[(gene_info_gene_name_info$`Population Mutation` != "none") & (gene_info_gene_name_info$`8-class DSSP Secondary Structure Properties` != "-"), ])
      bioRxiv_link = "https://doi.org/10.1101/693259"
        
      HTML(
        
        paste(
          "<i>","<font size = +3 color = Black face = Arial>",gene_name,"</i>","</font>","<br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","UniProt-AC: ","</font>","</b>",
          "<font size = -1 color = Black face = Arial>","<a href=", uniprotac_link_str, " target=_blank>", uniprotIDresearch,"</a><br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Protein Name: ","</font>","</b>",
          "<font size = -1 color = Black face = Arial>",protein_name,"<br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Protein length: ","</font>","</b>",
          "<font size = -1 color = Black face = Arial>",gene_length_info, " amino acid residues", "<br>","<br>",
          "<font size = +1 color = Black face = Arial>","Amino acid residues with pathogenic (ClinVar/HGMD) missense variants: ", gene_aa_path_var, "</font>", "<br>",
          "<font size = +1 color = Black face = Arial>","Residues with pathogenic variants and mapped on the 3D structure: ", gene_aa_path_var_wstruc, "</font>", "<br>", "<br>",
          "<font size = +1 color = Black face = Arial>","Amino acid residues with population (gnomAD) missense variants: ", gene_aa_pop_var, "</font>", "<br>",
          "<font size = +1 color = Black face = Arial>","Residues with population variants and mapped on the 3D structure: ", gene_aa_pop_var_wstruc, "</font>", "<br>", "<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Annotated Protein Class(es) for", "<i>", as.character(gene_name), "</i>", ": ","</font>","</b>",
          "<font size = +1 color = Black face = Arial>",pclass_all_str,"</font><br>","<br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Protein class-specific feature associated with missense variants", " (Ref: <a href=", bioRxiv_link, " target=_blank>paper link</a>)","</font>","</b>", "<br>",
          "<font size = -1 color = Black face = Arial>", "(see Documentation for feature definition and ascertainment)" ,"</font><br></br>",
          "<font color = Black face = Arial>",summary_gene
        )
      )
      
    }
    else if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      pclassSpecific_def <- subset(protein_class_def,superclassName == input$pclassNameselected)
      pclass_name = pclassSpecific_def$superclassName
      pclass_def = pclassSpecific_def$class_definition
      pclass_def_id = pclassSpecific_def$superclassID
      
      pclass_gnomadvariant = pclassSpecific_def$gnomAD_variant_mapped_in_this_study
      pclass_patientvariant = pclassSpecific_def$patient_variants_mapped_n_this_study

      pclass_gene_names_list = sort(unlist(strsplit(as.character(pclassSpecific_def$genes), split = ', ')))
      pclass_gene_count = length(pclass_gene_names_list)
      pclass_gene_names <- paste(pclass_gene_names_list, collapse = ', ')
      
      pclassIds <- pclass_def_id
      PC24_PF40_path <- subset(PC24_PF40_pathogenic, pcindex %in% pclassIds)
      PC24_PF40_pop <- subset(PC24_PF40_population, pcindex %in% pclassIds)
      
      pcname <- PC24_PF40_path$pcname
      pcnamelower <- PC24_PF40_path$pcnamelower
      
      i = 1
      #pathogenic feature for pclass -- start#
      pss3path_df <- PC24_PF40_path[i,4:6]
      pss3path_list <- c()
      for(fn in 1:ncol(pss3path_df))
      {
        pfindex = as.character(names(pss3path_df[1,fn]))
        fflag = as.integer(pss3path_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
          pss3path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
          pss3path_list <- c(pss3path_list, pss3path_for)
        }
      }
      if(length(pss3path_list) > 0){
        pss3path_str <- paste(pss3path_list, collapse = ", ")
      }else{
        pss3path_str <- "none"
      }
      
      
      pss8path_df <- PC24_PF40_path[i,7:14]
      pss8path_list <- c()
      for(fn in 1:ncol(pss8path_df))
      {
        pfindex = as.character(names(pss8path_df[1,fn]))
        fflag = as.integer(pss8path_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          pss8path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>"), collapse = "")
          pss8path_list <- c(pss8path_list, pss8path_for)
        }
      }
      if(length(pss8path_list) > 0){
        pss8path_str <- paste(pss8path_list, collapse = ", ")
      }else{
        pss8path_str <- "none"
      }
      
      pasapath_df <- PC24_PF40_path[i,15:19]
      pasapath_list <- c()
      for(fn in 1:ncol(pasapath_df))
      {
        pfindex = as.character(names(pasapath_df[1,fn]))
        fflag = as.integer(pasapath_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
          pasapath_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
          pasapath_list <- c(pasapath_list, pasapath_for)
        }
      }
      if(length(pasapath_list) > 0){
        pasapath_str <- paste(pasapath_list, collapse = ", ")
      }else{
        pasapath_str <- "none"
      }
      
      pcp8path_df <- PC24_PF40_path[i,20:27]
      pcp8path_list <- c()
      for(fn in 1:ncol(pcp8path_df))
      {
        pfindex = as.character(names(pcp8path_df[1,fn]))
        fflag = as.integer(pcp8path_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
          pcp8path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
          pcp8path_list <- c(pcp8path_list, pcp8path_for)
        }
      }
      if(length(pcp8path_list) > 0){
        pcp8path_str <- paste(pcp8path_list, collapse = ", ")
        pcp8path_str <- paste(pcp8path_str, "amino acids", sep = " ")
      }else{
        pcp8path_str <- "none"
      }
      
      pib4path_df <- PC24_PF40_path[i,28:31]
      pib4path_list <- c()
      for(fn in 1:ncol(pib4path_df))
      {
        pfindex = as.character(names(pib4path_df[1,fn]))
        fflag = as.integer(pib4path_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          pib4path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>"), collapse = "")
          pib4path_list <- c(pib4path_list, pib4path_for)
        }
      }
      if(length(pib4path_list) > 0){
        pib4path_str <- paste(pib4path_list, collapse = ", ")
        pib4path_str <- paste(pib4path_str, "interactions", sep = " ")
      }else{
        pib4path_str <- "none"
      }
      
      pptm6path_df <- PC24_PF40_path[i,32:37]
      pptm6path_list <- c()
      for(fn in 1:ncol(pptm6path_df))
      {
        pfindex = as.character(names(pptm6path_df[1,fn]))
        fflag = as.integer(pptm6path_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          pptm6path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>"), collapse = "")
          pptm6path_list <- c(pptm6path_list, pptm6path_for)
        }
      }
      if(length(pptm6path_list) > 0){
        pptm6path_str <- paste(pptm6path_list, collapse = ", ")
      }else{
        pptm6path_str <- "none"
      }
      
      pup6path_df <- PC24_PF40_path[i,38:43]
      pup6path_list <- c()
      for(fn in 1:ncol(pup6path_df))
      {
        pfindex = as.character(names(pup6path_df[1,fn]))
        fflag = as.integer(pup6path_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
          pup6path_for <- paste(c("<span style='color: #993300;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
          pup6path_list <- c(pup6path_list, pup6path_for)
        }
      }
      if(length(pup6path_list) > 0){
        pup6path_str <- paste(pup6path_list, collapse = ", ")
      }else{
        pup6path_str <- "none"
      }
      #pathogenic feature for pclass -- end#
      
      #popultaion feature for pclass -- start#
      pss3pop_df <- PC24_PF40_pop[i,4:6]
      pss3pop_list <- c()
      for(fn in 1:ncol(pss3pop_df))
      {
        pfindex = as.character(names(pss3pop_df[1,fn]))
        fflag = as.integer(pss3pop_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
          pss3pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
          pss3pop_list <- c(pss3pop_list, pss3pop_for)
        }
      }
      if(length(pss3pop_list) > 0){
        pss3pop_str <- paste(pss3pop_list, collapse = ", ")
      }else{
        pss3pop_str <- "none"
      }
      
      pss8pop_df <- PC24_PF40_pop[i,7:14]
      pss8pop_list <- c()
      for(fn in 1:ncol(pss8pop_df))
      {
        pfindex = as.character(names(pss8pop_df[1,fn]))
        fflag = as.integer(pss8pop_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          pss8pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>"), collapse = "")
          pss8pop_list <- c(pss8pop_list, pss8pop_for)
        }
      }
      if(length(pss8pop_list) > 0){
        pss8pop_str <- paste(pss8pop_list, collapse = ", ")
      }else{
        pss8pop_str <- "none"
      }
      
      pasapop_df <- PC24_PF40_pop[i,15:19]
      pasapop_list <- c()
      for(fn in 1:ncol(pasapop_df))
      {
        pfindex = as.character(names(pasapop_df[1,fn]))
        fflag = as.integer(pasapop_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
          pasapop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
          pasapop_list <- c(pasapop_list, pasapop_for)
        }
      }
      if(length(pasapop_list) > 0){
        pasapop_str <- paste(pasapop_list, collapse = ", ")
      }else{
        pasapop_str <- "none"
      }
      
      pcp8pop_df <- PC24_PF40_pop[i,20:27]
      pcp8pop_list <- c()
      for(fn in 1:ncol(pcp8pop_df))
      {
        pfindex = as.character(names(pcp8pop_df[1,fn]))
        fflag = as.integer(pcp8pop_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
          pcp8pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
          pcp8pop_list <- c(pcp8pop_list, pcp8pop_for)
        }
      }
      if(length(pcp8pop_list) > 0){
        pcp8pop_str <- paste(pcp8pop_list, collapse = ", ")
        pcp8pop_str <- paste(pcp8pop_str, "amino acids", sep = " ")
      }else{
        pcp8pop_str <- "none"
      }
      
      pib4pop_df <- PC24_PF40_pop[i,28:31]
      pib4pop_list <- c()
      for(fn in 1:ncol(pib4pop_df))
      {
        pfindex = as.character(names(pib4pop_df[1,fn]))
        fflag = as.integer(pib4pop_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          pib4pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>"), collapse = "")
          pib4pop_list <- c(pib4pop_list, pib4pop_for)
        }
      }
      if(length(pib4pop_list) > 0){
        pib4pop_str <- paste(pib4pop_list, collapse = ", ")
        pib4pop_str <- paste(pib4pop_str, "interactions", sep = " ")
      }else{
        pib4pop_str <- "none"
      }
      
      pptm6pop_df <- PC24_PF40_pop[i,32:37]
      pptm6pop_list <- c()
      for(fn in 1:ncol(pptm6pop_df))
      {
        pfindex = as.character(names(pptm6pop_df[1,fn]))
        fflag = as.integer(pptm6pop_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          pptm6pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>"), collapse = "")
          pptm6pop_list <- c(pptm6pop_list, pptm6pop_for)
        }
      }
      if(length(pptm6pop_list) > 0){
        pptm6pop_str <- paste(pptm6pop_list, collapse = ", ")
      }else{
        pptm6pop_str <- "none"
      }
      
      pup6pop_df <- PC24_PF40_pop[i,38:43]
      pup6pop_list <- c()
      for(fn in 1:ncol(pup6pop_df))
      {
        pfindex = as.character(names(pup6pop_df[1,fn]))
        fflag = as.integer(pup6pop_df[fn])
        if(fflag == 1){
          fdes <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdes'])
          fdes_xntd <- as.character(feature_index_df[feature_index_df$findex == pfindex, 'fdesextnd'])
          pup6pop_for <- paste(c("<span style='color: #000080;'><strong>", fdes, "</span></strong>", " (", fdes_xntd, ")"), collapse = "")
          pup6pop_list <- c(pup6pop_list, pup6pop_for)
        }
      }
      if(length(pup6pop_list) > 0){
        pup6pop_str <- paste(pup6pop_list, collapse = ", ")
      }else{
        pup6pop_str <- "none"
      }
      #population feature for pclass -- end#
      
      # see Documentation for feature definition and ascertainment
      
      this_pclass <- c("")
      this_pclass <- c(this_pclass, "<font size = -1 color = Black face = Arial>", "<strong>Features with significant enrichment of <span style='color: #993300';>pathogenic variants</span> in ", pcnamelower , "s: </strong>")
      this_pclass <- c(this_pclass, "<ol>")
      this_pclass <- c(this_pclass, "<li>", "3-class secondary structure: ", pss3path_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "8-class secondary structure: ", pss8path_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "Residue exposure level: ", pasapath_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "Physicochemical&nbsp;property: ", pcp8path_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "Protein-protein interaction/bond type: ", pib4path_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "Post-translational modification (within &lt;10 &#8491 spatial distance): ", pptm6path_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "UniProt-based function features: ", pup6path_str, "</li>")
      this_pclass <- c(this_pclass, "</ol>")
      this_pclass <- c(this_pclass, "<strong>Features with significant enrichment of <span style='color: #000080';>population variants</span> in ", pcnamelower , "s: </strong>")
      this_pclass <- c(this_pclass, "<ol>")
      this_pclass <- c(this_pclass, "<li>", "3-class secondary structure: ", pss3pop_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "8-class secondary structure: ", pss8pop_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "Residue exposure level: ", pasapop_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "Physicochemical&nbsp;property: ", pcp8pop_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "Protein-protein interaction/bond type: ", pib4pop_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "Post-translational modification (within &lt;10 &#8491 spatial distance): ", pptm6pop_str, "</li>")
      this_pclass <- c(this_pclass, "<li>", "UniProt-based function features: ", pup6pop_str, "</li>")
      this_pclass <- c(this_pclass, "</ol></font><br>")
      
      summary_gene <- paste(this_pclass, collapse = "")
      
      protein_class_def_url = "http://www.pantherdb.org/panther/ontologies.jsp"
      bioRxiv_link = "https://doi.org/10.1101/693259"
      
      HTML(
        paste(
          "<font size = +3 color = Black face = Arial>",pclass_name,"</font>","<br>","<br>",
          #"<font size = -1 color = Black face = Arial>",pclass_name,"<br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Protein Class Definition","</font>","</b>",
          "<font size = -1 color = Black face = Arial>", "<a href=", protein_class_def_url, " target=_blank>","(collected)","</a>","<br>",pclass_def,"<br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Annotated Genes in", pclass_name, "Class (count: ",pclass_gene_count,")","</font>","</b>","<br>",
          "<i>","<font size = -1 color = Black face = Arial>",pclass_gene_names,"</font>","</i>","<br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Number of population (gnomAD) missense Variants in", pclass_name, "class mapped on 3D structure:","</b>",pclass_gnomadvariant,"</font>", "<br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Number of pathogenic (ClinVar/HGMD) missense variants in", pclass_name, "class mapped on 3D structure: ","</b>",pclass_patientvariant,"</font>","<br>","<br>",
          "<b>","<font size = +1 color = Black face = Arial>","Protein class-specific feature associated with missense variants", " (Ref: <a href=", bioRxiv_link, " target=_blank>paper link</a>)","</font>","</b>", "<br>",
          "<font size = -1 color = Black face = Arial>",summary_gene
        )
      )
    }
  })
  
  twod_vis_gene_name <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      gene_name = geneSpecificpmissOrgmiss$Gene_Name
      proteinClass = as.character(geneSpecificpmissOrgmiss$Protein_Class_Name)
      
      HTML(
        paste(
          #"<font size = +1 color = Black face = Arial>","Gene","</font>", "<br>",
          "<br><br>",
          "<i>","<font size = +3 color = Black face = Arial>",gene_name,"</i>","</font>","<br>","<br>"
        )
      )
    }
  })
  
  threed_vis_gene_name <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      gene_name = geneSpecificpmissOrgmiss$Gene_Name
      proteinClass = as.character(geneSpecificpmissOrgmiss$Protein_Class_Name)
      
      HTML(
        paste(
          #"<font size = +1 color = Black face = Arial>","Gene","</font>", "<br>",
          "<br><br>",
          "<i>","<font size = +3 color = Black face = Arial>",gene_name,"</i>","</font>","<br>","<br>"
        )
      )
    }
  })
  
  aawisefeaturetable_vis_gene_name <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      gene_name = geneSpecificpmissOrgmiss$Gene_Name
      proteinClass = as.character(geneSpecificpmissOrgmiss$Protein_Class_Name)
      
      HTML(
        paste(
          #"<font size = +1 color = Black face = Arial>","Gene","</font>", "<br>",
          "<br><br>",
          "<i>","<font size = +3 color = Black face = Arial>",gene_name,"</i>","</font>","<br>","<br>"
        )
      )
    }
  })
  
  tabularview_vis_gene_name <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      gene_name = geneSpecificpmissOrgmiss$Gene_Name
      proteinClass = as.character(geneSpecificpmissOrgmiss$Protein_Class_Name)
      
      HTML(
        paste(
          #"<font size = +1 color = Black face = Arial>","Gene","</font>", "<br>",
          "<br><br>",
          "<i>","<font size = +3 color = Black face = Arial>",gene_name,"</i>","</font>","<br>","<br>"
        )
      )
    }
  })
  
  indextabularview_vis_gene_name <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      gene_name = geneSpecificpmissOrgmiss$Gene_Name
      proteinClass = as.character(geneSpecificpmissOrgmiss$Protein_Class_Name)
      
      HTML(
        paste(
          #"<font size = +1 color = Black face = Arial>","Gene","</font>", "<br>",
          "<br><br>",
          "<i>","<font size = +3 color = Black face = Arial>",gene_name,"</i>","</font>","<br>","<br>"
        )
      )
    }
  })
  
  pclasswise_vis_gene_name <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      pclassSpecific_def <- subset(protein_class_def,superclassName == input$pclassNameselected)
      pclass_name = pclassSpecific_def$superclassName
      
      HTML(
        paste(
          "<br><br>",
          "<font size = +3 color = Black face = Arial>",pclass_name,"<br>","<br>"
        )
      )
    }
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { 
      if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
        paste(input$geneSelected, "txt", sep=".")
      }
      
      #else if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      #  paste(input$pclassNameselected, input$filetype, sep=".")
      #}
    },
    content = function(file) {
      
      
      if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
        sep <- "\t"
        gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
        gene_wise_info <- read_delim(gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
        write.table(gene_wise_info, file,sep = sep,row.names = FALSE)
      }
    })
  
  output$indexTabledownloadData <- downloadHandler(
    filename = function() { 
      if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
        paste(input$geneSelected, "txt", sep=".")
      }
      
      #else if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      #  paste(input$pclassNameselected, input$filetype, sep=".")
      #}
    },
    content = function(file) {
      
      
      if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
        sep <- "\t"
        gene_name = paste("index_table_files/",input$geneSelected,".txt",sep='')
        gene_wise_info <- read_delim(gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
        write.table(gene_wise_info, file,sep = sep,row.names = FALSE)
      }
    })
  
  
  observeEvent(input$prot_feat_to_doc,{
    #updateTabItems(session, "mainTabset", "Documentation" )
    updateNavbarPage(session, "upperlayer", "Documentation")
  })
  output$Information2 = renderText(
    {
      details_information2()
    }
  )
  
  output$pClassVisGeneName = renderText(
    {
      if(input$homeSideBarTabSetPanel == 'Protein Class')
      {
        pclasswise_vis_gene_name()
      }
    }
  )
  
  output$TWODVisGeneName = renderText(
    {
      if(input$homeSideBarTabSetPanel == 'Select a Gene')
      {
        twod_vis_gene_name()
      }
    }
  )
  
  output$THREEDVisGeneName = renderText(
    {
      if(input$homeSideBarTabSetPanel == 'Select a Gene')
      {
        twod_vis_gene_name()
      }
    }
  )
  
  output$aaWiseFeatureTableGeneName = renderText(
    {
      if(input$homeSideBarTabSetPanel == 'Select a Gene')
      {
        aawisefeaturetable_vis_gene_name()
      }
    }
  )
  
  output$TablularViewGeneName = renderText(
    {
      if(input$homeSideBarTabSetPanel == 'Select a Gene')
      {
        tabularview_vis_gene_name()
      }
    }
  )
  
  output$IndexTablularViewGeneName = renderText(
    {
      if(input$homeSideBarTabSetPanel == 'Select a Gene')
      {
        indextabularview_vis_gene_name()
      }
    }
  )
  
  #output$aminoAcidSelection = renderText(
  #  {
  #HTML(
  #  paste(
  #     "<p align = right>", "<font size = \"4\">","Select Rows: " , "</font>","</p>"
  #  )
  #)
  #    tags$p(tags$font( size = "4"),"Select Rows: ")
  #  }
  #)
  
  output$structureaminoAcidSelection = renderText(
    {
      HTML(
        paste(
          "<p align = center>","Select Rows: ","</p>"
        )
      )
    }
  )
  
  output$IndexTableSummary = renderText(
    {
      HTML(
        paste(#\"#EFEFEF\"
          "<font size = -1 color = \"#5F5F5F\">", "Genomic Index of the Amino Acid (transcript>chromosome:strand:position1,position2,position3)" , "</font>" , "<br>",
          "<font size = \"-1\" color = \"#5F5F5F\">", "Structure Index of the Amino Acid (PDB id:chain:structure position)" , "</font>","<br>",
          "<font size = -1 color = \"#5F5F5F\">", "<strong>(.)</strong> in the structure index column indicates 3D-coordinate is not available." , "</font>"
        )
      )
    }
  )
  
  output$asaplotSummary = renderText(
    {
      HTML(
        paste(#\"#EFEFEF\"
          "<font size = -1 color = \"#5F5F5F\">", "-25 is a dummy value for residues with no structure" , "</font>" , "<br>"
        )
      )
    }
  )
  
  output$ptmplotsummary = renderText(
    {
      HTML(
        paste(#\"#EFEFEF\"
          "<font size = -1 color = \"#5F5F5F\">", "Only three smallest distances are visualized in 2D plots. For full annotation, download text file from Feature table tab" , "</font>" , "<br>"
        )
      )
    }
  )
  
  output$FeatureTableSummary = renderText(
    {
      HTML(
        paste(#\"#EFEFEF\"
          "<font size = -1 color = \"#5F5F5F\">", "<strong>Population mutation details format,</strong> AC = allele count, AN = allele number, and AF  = allele frequency" , "</font>" , "<br>",
          "<font size = -1 color = \"#5F5F5F\">", "<strong>Bond details format,</strong> interacting chain 1:amino acid position 1:amino acid 1-BOND/distance-interacting chain 2:amino acid position 2:amino acid 2" , "</font>", "<br>",
          "<font size = -1 color = \"#5F5F5F\">", "<strong>BOND annotation,</strong> D = disulfide bond, S = salt-bridge ionic bond, H = hydrogen bond, N = nonbonded Var der waals interaction" , "</font>", "<br>",
          "<font size = -1 color = \"#5F5F5F\">", "<strong>PTM details format,</strong> PDB id.chain id.PTM position on structure/spatial distance from this amino acid to the PTM position" , "</font>","<br>",
          "<font size = -1 color = \"#5F5F5F\">", "<strong>(*)</strong>: Maximum of two annotations are shown for these features. To see full annotation, please click <strong>Download</strong> on the top-left corner of the table." , "</font>","<br>",
          "<font size = -1 color = \"#5F5F5F\">", "<strong>(-)</strong> in the feature columns indicate annotation is not available." , "</font>"
        )
      )
    }
  )
  
  
  show_gcount <- eventReactive(input$submit, {
    category_wise_gmiss <- subset(func_wise_propteins, pclassName == input$pclassNameselected, select = c('gmiss_genes')) 
    categorywise_gmiss_genes <- data.frame(unique(data.frame(separate_rows(func_wise_gmiss_genes(), 'gmiss_genes'))))
    #hnsd
    ignore_null_gmiss_count <- subset(categorywise_gmiss_genes,categorywise_gmiss_genes$gmiss_genes != "null")
    no_of_gmiss <- nrow(ignore_null_gmiss_count)
    
    HTML(
      paste(
        "<b>","<font size=+1 color = Black face = Arial>","Genes with Population Mutation (",no_of_gmiss,")","</font>","</b>"
      )
    )
  })
  output$show_gmiss_count = renderText(
    {
      show_gcount()
    }
  )
  
  show_category_wise_gmiss_genes <- eventReactive(input$submit, {
    category_wise_gmiss <- subset(func_wise_propteins, pclassName == input$pclassNameselected, select = c('gmiss_genes')) 
    categorywise_gmiss_genes <- data.frame(unique(data.frame(separate_rows(func_wise_gmiss_genes(), 'gmiss_genes'))))
    #hnsd
    
    HTML(
      paste(
        "<i>","<font size=-1 color = Black face = Arial>",categorywise_gmiss_genes[,1],", ","</font>","</i>"
      )
    )
    
    
    
  })
  
  output$gmiss_geneInfo = renderText(
    {
      show_category_wise_gmiss_genes()
    }
  )
  
  show_pcount <- eventReactive(input$submit, {
    category_wise_pmiss <- subset(func_wise_propteins, pclassName == input$pclassNameselected, select = c('pmiss_genes')) 
    categorywise_pmiss_genes <- data.frame(unique(data.frame(separate_rows(func_wise_pmiss_genes(), 'pmiss_genes'))))
    #hnsd
    ignore_null_pmiss_count <- subset(categorywise_pmiss_genes,categorywise_pmiss_genes$pmiss_genes != "null")
    no_of_pmiss <- nrow(ignore_null_pmiss_count)
    
    HTML(
      paste(
        "<b>","<font size=+1 color = Black face = Arial>","Genes with Pathogenic Mutation (",no_of_pmiss,")","</font>","</b>"
      )
    )
  })
  output$show_pmiss_count = renderText(
    {
      show_pcount()
    }
  )
  
  show_category_wise_pmiss_genes <- eventReactive(input$submit, {
    category_wise_pmiss <- subset(func_wise_propteins, pclassName == input$pclassNameselected, select = c('pmiss_genes')) 
    categorywise_pmiss_genes <- data.frame(unique(data.frame(separate_rows(func_wise_pmiss_genes(), 'pmiss_genes'))))
    
    HTML(
      paste(
        
        "<i>","<font size=-1color = Black face = Arial>",categorywise_pmiss_genes[,1],", ","</font>","</i>"
      )
    )
    
    
    
  })
  
  output$pmiss_geneInfo = renderText(
    {
      show_category_wise_pmiss_genes()
    }
  )
  
  
  func_wise_gmiss_genes <- eventReactive(input$submit, {
    subset(func_wise_propteins, pclassName == input$pclassNameselected, select = c('gmiss_genes')) 
  })
  func_wise_gmiss_genes_separate <- eventReactive(input$submit, {
    data.frame(unique(data.frame(separate_rows(func_wise_gmiss_genes(), 'gmiss_genes')))) 
  })
  
  ### category-wise gene/protein list with patient missense variants mappable to soleved protein structure##
  func_wise_pmiss_genes <- eventReactive(input$submit, {
    subset(func_wise_propteins, pclassName == input$pclassNameselected, select = c('pmiss_genes')) 
  })
  func_wise_pmiss_genes_separate <- eventReactive(input$submit, {
    data.frame(unique(data.frame(separate_rows(func_wise_pmiss_genes(), 'pmiss_genes')))) 
  })
  
  ## Protein class sepecific data load
  #new
  class_specific_info <- eventReactive(input$submit, {
    protein_class_def[protein_class_def$superclassName == input$pclassNameselected,] 
  })
  
  
  class_specific_gene <- eventReactive(input$submit, {
    data.frame(unique(data.frame(separate_rows(class_specific_info(), 'genes')))) 
  })
  
  gmiss_dssp_ptm_category <- eventReactive(input$submit, {
    subset(gmiss_dssp_prp, geneName %in% class_specific_gene()$genes)
  })
  
  pmiss_dssp_ptm_category <- eventReactive(input$submit, {
    subset(pmiss_dssp_prp, geneName %in% class_specific_gene()$genes)
  })
  
  total_gmiss <- eventReactive(input$submit, { nrow(gmiss_dssp_ptm_category()) })
  total_pmiss <- eventReactive(input$submit, { nrow(pmiss_dssp_ptm_category()) })
  
  
  ## PLOTS -- SS8 --- START ##
  #SS8bar
  output$SS8bar <- renderPlot({
    #SS8 bar plot
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    
    gmiss_dsspSS_count_val <- data_frame(as.numeric(c(nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'B', ]), nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'E', ]),nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'G', ]), nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'H', ]), nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'I', ]),nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == '.', ]), nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'S', ]), nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'T', ]))))
    gmiss_dsspSS_lab <- data_frame(as.character(c("B","E","G","H","I",".","S","T")))
    gmiss_dsspSS_count_raw <- data.frame(cbind(gmiss_dsspSS_lab,gmiss_dsspSS_count_val), stringsAsFactors = FALSE)
    gmiss_dsspSS_count <- setNames(gmiss_dsspSS_count_raw, c("dsspSS", "sscount"))
    
    gmiss_dsspSS_count_prop <- gmiss_dsspSS_count$sscount / gmiss_class_total 
    gmiss_dsspSS_count_label <- data.frame(rep("Population Variant",nrow(gmiss_dsspSS_count))) 
    gmiss_dsspSS_count_freq_rest <- sum(gmiss_dsspSS_count$sscount) - gmiss_dsspSS_count$sscount 
    gmiss_dsspSS_count_comp <- cbind(gmiss_dsspSS_count,gmiss_dsspSS_count_prop,gmiss_dsspSS_count_label,gmiss_dsspSS_count_freq_rest) 
    gmiss_dsspSS_count_comp_formatted <- setNames(gmiss_dsspSS_count_comp, c("dsspSS", "freq", "prop", "variant","freq_rest")) 
    
    pmiss_dsspSS_count_val <- data_frame(as.numeric(c(nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'B', ]), nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'E', ]),nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'G', ]), nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'H', ]), nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'I', ]),nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == '.', ]), nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'S', ]), nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'T', ]))))
    pmiss_dsspSS_lab <- data_frame(as.character(c("B","E","G","H","I",".","S","T")))
    pmiss_dsspSS_count_raw <- data.frame(cbind(pmiss_dsspSS_lab,pmiss_dsspSS_count_val), stringsAsFactors = FALSE)
    pmiss_dsspSS_count <- setNames(pmiss_dsspSS_count_raw, c("dsspSS", "sscount"))
    
    pmiss_dsspSS_count_prop <- pmiss_dsspSS_count$sscount / pmiss_class_total 
    pmiss_dsspSS_count_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_dsspSS_count))) 
    pmiss_dsspSS_count_freq_rest <- sum(pmiss_dsspSS_count$sscount) - pmiss_dsspSS_count$sscount 
    pmiss_dsspSS_count_comp <- cbind(pmiss_dsspSS_count,pmiss_dsspSS_count_prop,pmiss_dsspSS_count_label,pmiss_dsspSS_count_freq_rest) 
    pmiss_dsspSS_count_comp_formatted <- setNames(pmiss_dsspSS_count_comp, c("dsspSS", "freq", "prop", "variant","freq_rest")) 
    
    dsspSS_proportion <- rbind(gmiss_dsspSS_count_comp_formatted,pmiss_dsspSS_count_comp_formatted)
    
    ggplot(dsspSS_proportion, aes(dsspSS_proportion$dsspSS, dsspSS_proportion$prop)) +   
      geom_bar(aes(fill = dsspSS_proportion$variant), position = "dodge", stat="identity", colour="black", width = 0.7) +
      ggtitle("Distribution of Variants in Structure Types") +
      xlab("8-class secondary structure types") + 
      ylab("Proportion of variants") +
      scale_y_continuous(labels = percent) +
      labs(fill="Variant Type") +
      theme_bw(base_size=14) +
      theme(legend.position = "top", legend.direction = "horizontal") +
      scale_fill_manual(name="", 
                        values = c("Population Variant"="blue", "Pathogenic Variant"="red"))+
      scale_x_discrete(labels=c("."="Loop","B"=expression(paste(beta,"-strand")),"E"=expression(paste(beta,"-sheet")),"G"=expression(paste("3"["10"],"-helix")),"H"=expression(paste(alpha,"-helix")),"I"=expression(paste(pi,"-helix")),"S"="Bend","T"="Turn"))
  })
  
  #SS8forest
  output$SS8forest <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    #SS8forest
    #pmiss counts for fischer test
    class_strand_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "B", ])
    class_sheet_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "E", ])
    class_alpha_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "H", ])
    class_3_10_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "G", ])
    class_pi_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "I", ])
    class_loop_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == ".", ])
    class_bend_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "S", ])
    class_turn_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "T", ])
    
    #gmiss counts for fischer's test
    class_strand_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "B", ])
    class_sheet_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "E", ])
    class_alpha_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "H", ])
    class_3_10_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "G", ])
    class_pi_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "I", ])
    class_loop_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == ".", ])
    class_bend_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "S", ])
    class_turn_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "T", ])
    
    #strand
    fisher_marix <- matrix(c(class_strand_pmiss_count, pmiss_class_total-class_strand_pmiss_count,class_strand_gmiss_count,gmiss_class_total-class_strand_gmiss_count), nrow=2)
    fisher_strand <- fisher.test(fisher_marix)
    
    #sheet
    fisher_marix <- matrix(c(class_sheet_pmiss_count, pmiss_class_total-class_sheet_pmiss_count,class_sheet_gmiss_count,gmiss_class_total-class_sheet_gmiss_count), nrow=2)
    fisher_sheet <- fisher.test(fisher_marix)
    
    #3_10
    fisher_marix <- matrix(c(class_3_10_pmiss_count, pmiss_class_total-class_3_10_pmiss_count,class_3_10_gmiss_count,gmiss_class_total-class_3_10_gmiss_count), nrow=2)
    fisher_3_10 <- fisher.test(fisher_marix)
    
    #alpha
    fisher_marix <- matrix(c(class_alpha_pmiss_count, pmiss_class_total-class_alpha_pmiss_count,class_alpha_gmiss_count,gmiss_class_total-class_alpha_gmiss_count), nrow=2)
    fisher_alpha <- fisher.test(fisher_marix)
    
    #pi
    fisher_marix <- matrix(c(class_pi_pmiss_count, pmiss_class_total-class_pi_pmiss_count,class_pi_gmiss_count,gmiss_class_total-class_pi_gmiss_count), nrow=2)
    fisher_pi <- fisher.test(fisher_marix)
    
    #loop
    fisher_marix <- matrix(c(class_loop_pmiss_count, pmiss_class_total-class_loop_pmiss_count,class_loop_gmiss_count,gmiss_class_total-class_loop_gmiss_count), nrow=2)
    fisher_loop <- fisher.test(fisher_marix)
    
    #bend
    fisher_marix <- matrix(c(class_bend_pmiss_count, pmiss_class_total-class_bend_pmiss_count,class_bend_gmiss_count,gmiss_class_total-class_bend_gmiss_count), nrow=2)
    fisher_bend <- fisher.test(fisher_marix)
    
    #turn
    fisher_marix <- matrix(c(class_turn_pmiss_count, pmiss_class_total-class_turn_pmiss_count,class_turn_gmiss_count,gmiss_class_total-class_turn_gmiss_count), nrow=2)
    fisher_turn <- fisher.test(fisher_marix)
    
    #Forest Plot -- comm
    #Forest Plot -- comm
    odds_ss8  <- c(fisher_bend$estimate, fisher_turn$estimate, fisher_loop$estimate,fisher_3_10$estimate, fisher_alpha$estimate, fisher_pi$estimate, fisher_strand$estimate,fisher_sheet$estimate) 
    lower_ss8 <- c(fisher_bend$conf.int[1], fisher_turn$conf.int[1], fisher_loop$conf.int[1], fisher_3_10$conf.int[1], fisher_alpha$conf.int[1], fisher_pi$conf.int[1], fisher_strand$conf.int[1],fisher_sheet$conf.int[1]) 
    upper_ss8 <- c(fisher_bend$conf.int[2], fisher_turn$conf.int[2], fisher_loop$conf.int[2], fisher_3_10$conf.int[2], fisher_alpha$conf.int[2], fisher_pi$conf.int[2], fisher_strand$conf.int[2],fisher_sheet$conf.int[2]) 
    p_value_ss8 <- c(fisher_bend$p.value, fisher_turn$p.value, fisher_loop$p.value, fisher_3_10$p.value, fisher_alpha$p.value, fisher_pi$p.value, fisher_strand$p.value,fisher_sheet$p.value) 
    label_ss8 <- c('1.1.Bend', '1.2.Turn', '1.3.Loop', '2.1.3/10-Helix', '2.2.a-Helix', '2.3.pi-Helix', '3.1.b-Strand', '3.2.b-Sheet') 
    for_dsspSS_raw <- data.frame(label_ss8, odds_ss8, lower_ss8, upper_ss8, p_value_ss8) 
    for_dsspSS <- setNames(for_dsspSS_raw, c("structure_type", "odds", "lower", "upper", "p_value"))  
    
    pcut = 5.0e-05;
    for_dsspSS$odds <- ifelse(for_dsspSS$odds == 0.0, 0.05, for_dsspSS$odds)
    for_dsspSS$lower <- ifelse(for_dsspSS$odds == 0.0, 0.05, for_dsspSS$lower)
    for_dsspSS$lower <- ifelse(for_dsspSS$lower < 0.05, 0.05, for_dsspSS$lower)
    for_dsspSS$odds <- ifelse(for_dsspSS$odds > 20.0, 20.0, for_dsspSS$odds)
    for_dsspSS$upper <- ifelse(for_dsspSS$odds > 20.0, 20.0, for_dsspSS$upper)
    for_dsspSS$upper <- ifelse(for_dsspSS$upper > 20.0, 20.0, for_dsspSS$upper)
    for_dsspSS$p_value <- ifelse(for_dsspSS$p_value < 1.0e-300, 1.0e-300, for_dsspSS$p_value)
    for_dsspSS$p_value_label <- ifelse(for_dsspSS$p_value < pcut, paste(sprintf("%0.2e", for_dsspSS$p_value), "*)", sep = "("), sprintf("%0.2e", for_dsspSS$p_value))
    
    ggplot(for_dsspSS, aes(x=for_dsspSS$structure_type, y=for_dsspSS$odds, ymin=for_dsspSS$lower, ymax=for_dsspSS$upper))+
      geom_errorbar(aes(ymin = for_dsspSS$lower, ymax = for_dsspSS$upper), width = 0.15, size=1,colour=ifelse(for_dsspSS$p_value < pcut, "black", "darkgrey")) +
      geom_point(size=4, colour=ifelse(for_dsspSS$odds < 1.0, "blue", "red")) +
      geom_text(data = for_dsspSS,aes(x=rep(1.4:8.4,1), y=for_dsspSS$odds), label = for_dsspSS$p_value_label, fontface = "bold", size = 4,hjust="inward") +
      geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
      #scale_shape_identity() +
      coord_flip() +  # flip coordinates (puts labels on y axis)
      #scale_y_log10() +
      theme_bw(base_size = 14) +  # use a white background
      scale_y_log10(breaks=c(0.05,0.1,0.3,1.0,3.0,10.0,20.0), limits=c(0.05,20.0))+
      theme(legend.position = "top", legend.direction = "horizontal") +
      labs(title="Pathogenic vs. Population Variant Enrichment Analysis", 
           subtitle="Odds>1.0 (<1.0) indicates enrichment in pathogenic (populations) variants\np-value cut-off = 5.0e-05",
           x="8-class secondary strcuture types",
           y=expression(paste("OR (95% CI) "["log10"],""))) +
      scale_x_discrete(labels=c("1.3.Loop"="loop","3.1.b-Strand"=expression(paste(beta,"-strand")),"3.2.b-Sheet"=expression(paste(beta,"-sheet")),"2.1.3/10-Helix"=expression(paste("3"["10"],"-helix")),"2.2.a-Helix"=expression(paste(alpha,"-helix")),"2.3.pi-Helix"=expression(paste(pi,"-helix")),"1.1.Bend"="bend","1.2.Turn"="turn"))
    
    
  })
  
  ## PLOTS -- SS8 --- END ##
  
  
  ## PLOTS -- SS3 --- START ##
  # count -- all for bar plot #
  output$SS3bar <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    
    gmiss_conSS_count_val <- data_frame(as.numeric(c(nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'B', ])+nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'E', ]),nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'G', ])+nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'H', ])+nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'I', ]),nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == '.', ])+nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'S', ])+nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == 'T', ]))))
    gmiss_conSS_lab <- data_frame(as.character(c("E","H","C")))
    gmiss_conSS_count_raw <- data.frame(cbind(gmiss_conSS_lab,gmiss_conSS_count_val), stringsAsFactors = FALSE)
    gmiss_conSS_count <- setNames(gmiss_conSS_count_raw, c("conSS", "sscount"))
    
    gmiss_conSS_count_prop <- gmiss_conSS_count$sscount / gmiss_class_total 
    gmiss_conSS_count_label <- data.frame(rep("Population Variant",nrow(gmiss_conSS_count))) 
    gmiss_conSS_count_freq_rest <- sum(gmiss_conSS_count$sscount) - gmiss_conSS_count$sscount 
    gmiss_conSS_count_comp <- cbind(gmiss_conSS_count,gmiss_conSS_count_prop,gmiss_conSS_count_label,gmiss_conSS_count_freq_rest) 
    gmiss_conSS_count_comp_formatted <- setNames(gmiss_conSS_count_comp, c("conSS", "freq", "prop", "variant","freq_rest")) 
    
    pmiss_conSS_count_val <- data_frame(as.numeric(c(nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'B', ])+nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'E', ]),nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'G', ])+nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'H', ])+nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'I', ]),nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == '.', ])+nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'S', ])+nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == 'T', ]))))
    pmiss_conSS_lab <- data_frame(as.character(c("E","H","C")))
    pmiss_conSS_count_raw <- data.frame(cbind(pmiss_conSS_lab,pmiss_conSS_count_val), stringsAsFactors = FALSE)
    pmiss_conSS_count <- setNames(pmiss_conSS_count_raw, c("conSS", "sscount"))
    
    pmiss_conSS_count_prop <- pmiss_conSS_count$sscount / pmiss_class_total 
    pmiss_conSS_count_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_conSS_count))) 
    pmiss_conSS_count_freq_rest <- sum(pmiss_conSS_count$sscount) - pmiss_conSS_count$sscount 
    pmiss_conSS_count_comp <- cbind(pmiss_conSS_count,pmiss_conSS_count_prop,pmiss_conSS_count_label,pmiss_conSS_count_freq_rest) 
    pmiss_conSS_count_comp_formatted <- setNames(pmiss_conSS_count_comp, c("conSS", "freq", "prop", "variant","freq_rest")) 
    
    conSS_proportion <- rbind(gmiss_conSS_count_comp_formatted,pmiss_conSS_count_comp_formatted)
    
    
    ggplot(conSS_proportion, aes(conSS_proportion$conSS, conSS_proportion$prop)) +   
      geom_bar(aes(fill = conSS_proportion$variant), position = "dodge", stat="identity", colour="black", width = 0.7) +
      ggtitle("Distribution of Variants in Structure Types") +
      xlab("3-class secondary structure types") + 
      ylab("Proportion of variants") +
      scale_y_continuous(labels = percent) +
      labs(fill="Variant type") +
      theme_bw(base_size=14) +
      theme(legend.position = "top", legend.direction = "horizontal") +
      scale_fill_manual(name="", 
                        values = c("Population Variant"="blue", "Pathogenic Variant"="red"))+
      scale_x_discrete(labels=c("C"="coils","E"=expression(paste(beta,"-strand/sheet")),"H"="helices"))
    
  })
  
  #SS3 - forest plot prepareation#
  output$SS3forest <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    
    #pmiss counts for fischer test
    class_strand_dssp_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "B", ])
    class_sheet_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "E", ])
    class_alpha_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "H", ])
    class_3_10_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "G", ])
    class_pi_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "I", ])
    class_loop_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == ".", ])
    class_bend_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "S", ])
    class_turn_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dsspSS == "T", ])
    
    class_strand_pmiss_count <- class_strand_dssp_pmiss_count + class_sheet_pmiss_count
    class_helix_pmiss_count <- class_alpha_pmiss_count + class_3_10_pmiss_count + class_pi_pmiss_count
    class_coil_pmiss_count <- class_loop_pmiss_count + class_bend_pmiss_count + class_turn_pmiss_count
    
    #gmiss counts for fischer's test
    class_strand_dssp_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "B", ])
    class_sheet_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "E", ])
    class_alpha_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "H", ])
    class_3_10_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "G", ])
    class_pi_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "I", ])
    class_loop_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == ".", ])
    class_bend_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "S", ])
    class_turn_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dsspSS == "T", ])
    
    class_strand_gmiss_count <- class_strand_dssp_gmiss_count + class_sheet_gmiss_count
    class_helix_gmiss_count <- class_alpha_gmiss_count + class_3_10_gmiss_count + class_pi_gmiss_count
    class_coil_gmiss_count <- class_loop_gmiss_count + class_bend_gmiss_count + class_turn_gmiss_count
    
    fisher_marix <- matrix(c(class_strand_pmiss_count, pmiss_class_total-class_strand_pmiss_count,class_strand_gmiss_count,gmiss_class_total-class_strand_gmiss_count), nrow=2)
    fisher_strand <- fisher.test(fisher_marix)
    
    fisher_marix <- matrix(c(class_helix_pmiss_count, pmiss_class_total-class_helix_pmiss_count,class_helix_gmiss_count,gmiss_class_total-class_helix_gmiss_count), nrow=2)
    fisher_helix <- fisher.test(fisher_marix)
    
    fisher_marix <- matrix(c(class_coil_pmiss_count, pmiss_class_total-class_coil_pmiss_count,class_coil_gmiss_count,gmiss_class_total-class_coil_gmiss_count), nrow=2)
    fisher_coil <- fisher.test(fisher_marix)
    
    odds <- c(fisher_coil$estimate,fisher_helix$estimate,fisher_strand$estimate) 
    lower <- c(fisher_coil$conf.int[1],fisher_helix$conf.int[1],fisher_strand$conf.int[1]) 
    upper <- c(fisher_coil$conf.int[2],fisher_helix$conf.int[2],fisher_strand$conf.int[2]) 
    p_value <- c(fisher_coil$p.value,fisher_helix$p.value,fisher_strand$p.value) 
    label <- c('Coil', 'Helix', 'Beta') 
    for_conSS_raw <- data.frame(label, odds, lower, upper, p_value)
    for_conSS <- setNames(for_conSS_raw, c("structure_type", "odds", "lower", "upper", "p_value"))  
    
    pcut = 5.0e-05;
    for_conSS$odds <- ifelse(for_conSS$odds == 0.0, 0.05, for_conSS$odds)
    for_conSS$lower <- ifelse(for_conSS$odds == 0.0, 0.05, for_conSS$lower)
    for_conSS$lower <- ifelse(for_conSS$lower < 0.05, 0.05, for_conSS$lower)
    for_conSS$odds <- ifelse(for_conSS$odds > 20.0, 20.0, for_conSS$odds)
    for_conSS$upper <- ifelse(for_conSS$odds > 20.0, 20.0, for_conSS$upper)
    for_conSS$upper <- ifelse(for_conSS$upper > 20.0, 20.0, for_conSS$upper)
    for_conSS$p_value <- ifelse(for_conSS$p_value < 1.0e-300, 1.0e-300, for_conSS$p_value)
    for_conSS$p_value_label <- ifelse(for_conSS$p_value < pcut, paste(sprintf("%0.2e", for_conSS$p_value), "*)", sep = "("), sprintf("%0.2e", for_conSS$p_value))
    
    ggplot(for_conSS, aes(x=for_conSS$structure_type, y=for_conSS$odds, ymin=for_conSS$lower, ymax=for_conSS$upper))+
      geom_errorbar(aes(ymin = for_conSS$lower, ymax = for_conSS$upper), width = 0.15, size=1,colour=ifelse(for_conSS$p_value < pcut, "black", "darkgrey")) +
      geom_point(size=4, colour=ifelse(for_conSS$odds < 1.0, "blue", "red")) +
      geom_text(data = for_conSS,aes(x=c(2.3,3.3,1.3), y=for_conSS$odds), label = for_conSS$p_value_label, fontface = "bold", size = 4) +
      geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
      coord_flip() +  # flip coordinates (puts labels on y axis)
      #scale_y_log10() +
      theme_bw(base_size = 14) +  # use a white background
      scale_y_log10(breaks=c(0.05,0.1,0.3,1.0,3.0,10.0,20.0), limits=c(0.2,5.0))+
      labs(title="Pathogenic vs. Population Variant Enrichment Analysis", 
           subtitle="Odds>1.0 (<1.0) indicates enrichment in pathogenic (populations) variants\np-value cut-off = 5.0e-05",
           #caption="Source: mpg",
           x="3-class secondary strcuture types",
           y=expression(paste("OR (95% CI) "["log10"],"")))+
      theme(legend.position = "top", legend.direction = "horizontal") +
      scale_x_discrete(labels=c("Coil"="coils","Helix"="helices","Beta"=expression(paste(beta,"-strand/sheet"))))
    
  })
  ## PLOTS -- SS3 --- END ##
  
  ## PLOTS -- ASA --- START ##
  ## ASA: _all ##
  output$ASAbox <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    
    gmiss_asa_raw <- data.frame(gmiss_prp_class$asa) 
    gmiss_asa_raw_label <- data.frame(rep("Population Variant",nrow(gmiss_asa_raw))) 
    gmiss_asa_raw_comb <- cbind(gmiss_asa_raw, gmiss_asa_raw_label) 
    gmiss_asa_formatted <- setNames(gmiss_asa_raw_comb, c("asa", "type")) 
    
    pmiss_asa_raw <- data.frame(pmiss_prp_class$asa) 
    pmiss_asa_raw_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_asa_raw))) 
    pmiss_asa_raw_comb <- cbind(pmiss_asa_raw, pmiss_asa_raw_label) 
    pmiss_asa_formatted <- setNames(pmiss_asa_raw_comb, c("asa", "type"))
    
    asa_combined <- rbind(gmiss_asa_formatted,pmiss_asa_formatted) 
    
    dodge <- position_dodge(width = T) 
    
    my_comparisons <- list( c("Population Variant", "Pathogenic Variant")) 
    pmiss_max_asa = max(pmiss_prp_class$asa) 
    gmiss_max_asa = max(gmiss_prp_class$asa)
    pmiss_rasa_med = median(pmiss_prp_class$asa)
    gmiss_rasa_med = median(gmiss_prp_class$asa)
    
    ggplot(asa_combined, aes(asa_combined$type, asa_combined$asa)) +
      geom_violin(position = dodge) +
      geom_boxplot(aes(fill=factor(asa_combined$type)),width=.1, position = dodge) +
      stat_compare_means(comparisons = my_comparisons) +
      ggplot2::annotate("text", x = 0.8, y = gmiss_max_asa, label = sprintf("median: %.2f", gmiss_rasa_med)) +
      ggplot2::annotate("text", x = 2.1, y = gmiss_max_asa, label = sprintf("median: %.2f", pmiss_rasa_med)) +
      theme_bw(base_size = 14) +
      theme(legend.position = "top", legend.direction = "horizontal") +
      labs(title="Distribution of Accessible Surface Area", 
           x="Variant types",
           y=expression(paste("Accessible surface area (",ring(A)^"2", ")")),
           fill="Variant Type") +
      scale_fill_manual(name="", 
                        values = c("Population Variant"="blue", "Pathogenic Variant"="red"))
    
  })
  
  output$RSAforest <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    
    #ExposureForest
    #pmiss counts for fischer test
    class_core_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$rsaCorrected < 0.05, ])
    class_buried_pmiss_count <- nrow(pmiss_prp_class[((pmiss_prp_class$rsaCorrected >= 0.05) & (pmiss_prp_class$rsaCorrected < 0.25)), ])
    class_medburied_pmiss_count <- nrow(pmiss_prp_class[((pmiss_prp_class$rsaCorrected >= 0.25) & (pmiss_prp_class$rsaCorrected < 0.50)), ])
    class_medexposed_pmiss_count <- nrow(pmiss_prp_class[((pmiss_prp_class$rsaCorrected >= 0.50) & (pmiss_prp_class$rsaCorrected < 0.75)), ])
    class_exposed_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$rsaCorrected >= 0.75, ])
    
    #gmiss counts for fischer's test
    class_core_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$rsaCorrected < 0.05, ])
    class_buried_gmiss_count <- nrow(gmiss_prp_class[((gmiss_prp_class$rsaCorrected >= 0.05) & (gmiss_prp_class$rsaCorrected < 0.25)), ])
    class_medburied_gmiss_count <- nrow(gmiss_prp_class[((gmiss_prp_class$rsaCorrected >= 0.25) & (gmiss_prp_class$rsaCorrected < 0.50)), ])
    class_medexposed_gmiss_count <- nrow(gmiss_prp_class[((gmiss_prp_class$rsaCorrected >= 0.50) & (gmiss_prp_class$rsaCorrected < 0.75)), ])
    class_exposed_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$rsaCorrected >= 0.75, ])
    
    #core
    fisher_marix <- matrix(c(class_core_pmiss_count, pmiss_class_total-class_core_pmiss_count,class_core_gmiss_count,gmiss_class_total-class_core_gmiss_count), nrow=2)
    fisher_core <- fisher.test(fisher_marix)
    
    #buried
    fisher_marix <- matrix(c(class_buried_pmiss_count, pmiss_class_total-class_buried_pmiss_count,class_buried_gmiss_count,gmiss_class_total-class_buried_gmiss_count), nrow=2)
    fisher_buried <- fisher.test(fisher_marix)
    
    #medburied
    fisher_marix <- matrix(c(class_medburied_pmiss_count, pmiss_class_total-class_medburied_pmiss_count,class_medburied_gmiss_count,gmiss_class_total-class_medburied_gmiss_count), nrow=2)
    fisher_medburied <- fisher.test(fisher_marix)
    
    #medexposed
    fisher_marix <- matrix(c(class_medexposed_pmiss_count, pmiss_class_total-class_medexposed_pmiss_count,class_medexposed_gmiss_count,gmiss_class_total-class_medexposed_gmiss_count), nrow=2)
    fisher_medexposed <- fisher.test(fisher_marix)
    
    #exposed
    fisher_marix <- matrix(c(class_exposed_pmiss_count, pmiss_class_total-class_exposed_pmiss_count,class_exposed_gmiss_count,gmiss_class_total-class_exposed_gmiss_count), nrow=2)
    fisher_exposed <- fisher.test(fisher_marix)
    
    #Forest Plot -- comm
    odds_rsa5  <- c(fisher_exposed$estimate, fisher_medexposed$estimate, fisher_medburied$estimate, fisher_buried$estimate,fisher_core$estimate) 
    lower_rsa5 <- c(fisher_exposed$conf.int[1], fisher_medexposed$conf.int[1], fisher_medburied$conf.int[1], fisher_buried$conf.int[1],fisher_core$conf.int[1]) 
    upper_rsa5 <- c(fisher_exposed$conf.int[2], fisher_medexposed$conf.int[2], fisher_medburied$conf.int[2], fisher_buried$conf.int[2],fisher_core$conf.int[2]) 
    p_value_rsa5 <- c(fisher_exposed$p.value, fisher_medexposed$p.value, fisher_medburied$p.value, fisher_buried$p.value,fisher_core$p.value) 
    label_rsa5 <- c('1.1.Exposed', '1.2.Medium-exposed', '1.3.Medium-buried', '2.1.Buried', '2.2.Core') 
    for_dsspASA_raw <- data.frame(label_rsa5, odds_rsa5, lower_rsa5, upper_rsa5, p_value_rsa5) 
    for_dsspASA <- setNames(for_dsspASA_raw, c("exposure_type", "odds", "lower", "upper", "p_value"))  
    
    pcut = 5.0e-05;
    for_dsspASA$odds <- ifelse(for_dsspASA$odds == 0.0, 0.05, for_dsspASA$odds)
    for_dsspASA$lower <- ifelse(for_dsspASA$odds == 0.0, 0.05, for_dsspASA$lower)
    for_dsspASA$lower <- ifelse(for_dsspASA$lower < 0.05, 0.05, for_dsspASA$lower)
    for_dsspASA$odds <- ifelse(for_dsspASA$odds > 20.0, 20.0, for_dsspASA$odds)
    for_dsspASA$upper <- ifelse(for_dsspASA$odds > 20.0, 20.0, for_dsspASA$upper)
    for_dsspASA$upper <- ifelse(for_dsspASA$upper > 20.0, 20.0, for_dsspASA$upper)
    for_dsspASA$p_value <- ifelse(for_dsspASA$p_value < 1.0e-300, 1.0e-300, for_dsspASA$p_value)
    for_dsspASA$p_value_label <- ifelse(for_dsspASA$p_value < pcut, paste(sprintf("%0.2e", for_dsspASA$p_value), "*)", sep = "("), sprintf("%0.2e", for_dsspASA$p_value))
    
    ggplot(for_dsspASA, aes(x=for_dsspASA$exposure_type, y=for_dsspASA$odds, ymin=for_dsspASA$lower, ymax=for_dsspASA$upper))+
      geom_errorbar(aes(ymin = for_dsspASA$lower, ymax = for_dsspASA$upper), width = 0.15, size=1,colour=ifelse(for_dsspASA$p_value < pcut, "black", "darkgrey")) +
      geom_point(size=5, colour=ifelse(for_dsspASA$odds < 1.0, "blue", "red")) +
      geom_text(data = for_dsspASA,aes(x=rep(1.4:5.4,1), y=for_dsspASA$odds), label = for_dsspASA$p_value_label, fontface = "bold", size = 4,hjust="inward") +
      geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
      coord_flip() +  # flip coordinates (puts labels on y axis)
      theme_bw(base_size = 14) +  # use a white background
      scale_y_log10(breaks=c(0.05,0.1,0.3,1.0,3.0,10.0,20.0), limits=c(0.05,20.0))+
      theme(legend.position = "top", legend.direction = "horizontal") +
      labs(title="Pathogenic vs. Population Variant Enrichment Analysis", 
           subtitle="Odds>1.0 (<1.0) indicates enrichment in pathogenic (populations) variants\np-value cut-off = 5.0e-05",
           x="Exposure levels (% accessible surface area)",
           y=expression(paste("Odds (95% CI) "["log10"],""))) +
      scale_x_discrete(labels=c("1.1.Exposed"="exposed\n(>75%)","1.2.Medium-exposed"="medium-exposed\n(25-50)%","1.3.Medium-buried"="medium-buried\n(5-25)%","2.1.Buried"="buried\n(5-25)%","2.2.Core"="core\n(< 5%)"))
    
  })
  
  ## PLOTS -- chem --- START ##chemProp
  output$chemPropbar <- renderPlot({
    #Chemi8 bar plot
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    
    class_aliphatic_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "A", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "I", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "L", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "M", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "V", ]))
    class_aromatic_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "F", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "W", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "Y", ]))
    class_hydrophobic_gmiss_count <- sum(class_aliphatic_gmiss_count,class_aromatic_gmiss_count)
    class_positively_charged_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "H", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "K", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "R", ]))
    class_negatively_charged_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "D", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "E", ]))
    class_neutral_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "N", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "Q", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "S", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "T", ]))
    class_polar_gmiss_count <- sum(class_positively_charged_gmiss_count,class_negatively_charged_gmiss_count,class_neutral_gmiss_count)
    class_special_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "C", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "G", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "P", ]))
    
    gmiss_chem_count_val <- data_frame(c(class_aliphatic_gmiss_count, class_aromatic_gmiss_count,class_hydrophobic_gmiss_count, class_positively_charged_gmiss_count, class_negatively_charged_gmiss_count,class_neutral_gmiss_count, class_polar_gmiss_count, class_special_gmiss_count))
    gmiss_chem_lab <- data_frame(as.character(c("aliphatic","aromatic","hydrophobic","pos_charged","neg_charged","neutral","polar","special")))
    gmiss_chem_count_raw <- data.frame(cbind(gmiss_chem_lab,gmiss_chem_count_val), stringsAsFactors = FALSE)
    gmiss_chem_count <- setNames(gmiss_chem_count_raw, c("chem", "aacount"))
    
    gmiss_chem_count_prop <- gmiss_chem_count$aacount / gmiss_class_total 
    gmiss_chem_count_label <- data.frame(rep("Population Variant",nrow(gmiss_chem_count))) 
    gmiss_chem_count_freq_rest <- sum(gmiss_chem_count$aacount) - gmiss_chem_count$aacount 
    gmiss_chem_count_comp <- cbind(gmiss_chem_count,gmiss_chem_count_prop,gmiss_chem_count_label,gmiss_chem_count_freq_rest) 
    gmiss_chem_count_comp_formatted <- setNames(gmiss_chem_count_comp, c("chem", "freq", "prop", "variant","freq_rest")) 
    
    class_aliphatic_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "A", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "I", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "L", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "M", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "V", ]))
    class_aromatic_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "F", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "W", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "Y", ]))
    class_hydrophobic_pmiss_count <- sum(class_aliphatic_pmiss_count,class_aromatic_pmiss_count)
    class_positively_charged_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "H", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "K", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "R", ]))
    class_negatively_charged_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "D", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "E", ]))
    class_neutral_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "N", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "Q", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "S", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "T", ]))
    class_polar_pmiss_count <- sum(class_positively_charged_pmiss_count,class_negatively_charged_pmiss_count,class_neutral_pmiss_count)
    class_special_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "C", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "G", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "P", ]))
    
    pmiss_chem_count_val <- data_frame(c(class_aliphatic_pmiss_count, class_aromatic_pmiss_count,class_hydrophobic_pmiss_count, class_positively_charged_pmiss_count, class_negatively_charged_pmiss_count,class_neutral_pmiss_count, class_polar_pmiss_count, class_special_pmiss_count))
    pmiss_chem_lab <- data_frame(as.character(c("aliphatic","aromatic","hydrophobic","pos_charged","neg_charged","neutral","polar","special")))
    pmiss_chem_count_raw <- data.frame(cbind(pmiss_chem_lab,pmiss_chem_count_val), stringsAsFactors = FALSE)
    pmiss_chem_count <- setNames(pmiss_chem_count_raw, c("chem", "aacount"))
    
    pmiss_chem_count_prop <- pmiss_chem_count$aacount / pmiss_class_total 
    pmiss_chem_count_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_chem_count))) 
    pmiss_chem_count_freq_rest <- sum(pmiss_chem_count$aacount) - pmiss_chem_count$aacount 
    pmiss_chem_count_comp <- cbind(pmiss_chem_count,pmiss_chem_count_prop,pmiss_chem_count_label,pmiss_chem_count_freq_rest) 
    pmiss_chem_count_comp_formatted <- setNames(pmiss_chem_count_comp, c("chem", "freq", "prop", "variant","freq_rest")) 
    
    chem_proportion <- rbind(gmiss_chem_count_comp_formatted,pmiss_chem_count_comp_formatted)
    
    ggplot(chem_proportion, aes(chem_proportion$chem, chem_proportion$prop)) +   
      geom_bar(aes(fill = chem_proportion$variant), position = "dodge", stat="identity", colour="black", width = 0.7) +
      ggtitle("Distribution of Variants in Amino Acids of Different Properties") +
      xlab("Pysicochemical properties") + 
      ylab("Proportion of variants") +
      scale_y_continuous(labels = percent) +
      labs(fill="Variant Type") +
      theme_bw(base_size=14) +
      theme(legend.position = "top", legend.direction = "horizontal") +
      scale_fill_manual(name="", 
                        values = c("Population Variant"="blue", "Pathogenic Variant"="red"))+
      scale_x_discrete(labels=c("aliphatic"="Aliphatic","aromatic"="Aromatic","hydrophobic"="Hydro-\nphobic","pos_charged"="Positively\ncharged","neg_charged"="Negatively\ncharged","neutral"="Neutral","polar"="Polar","special"="Special")) +
      theme(axis.text.x = element_text(angle=0, vjust=0.6))
    
    
  })
  
  output$chemPropforest <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    #chemforest
    #pmiss counts for fischer test
    class_aliphatic_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "A", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "I", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "L", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "M", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "V", ]))
    class_aromatic_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "F", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "W", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "Y", ]))
    class_hydrophobic_pmiss_count <- sum(class_aliphatic_pmiss_count,class_aromatic_pmiss_count)
    class_positively_charged_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "H", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "K", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "R", ]))
    class_negatively_charged_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "D", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "E", ]))
    class_neutral_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "N", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "Q", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "S", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "T", ]))
    class_polar_pmiss_count <- sum(class_positively_charged_pmiss_count,class_negatively_charged_pmiss_count,class_neutral_pmiss_count)
    class_special_pmiss_count <- sum(nrow(pmiss_prp_class[pmiss_prp_class$refAA == "C", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "G", ]), nrow(pmiss_prp_class[pmiss_prp_class$refAA == "P", ]))
    
    #gmiss counts for fischer's test
    class_aliphatic_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "A", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "I", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "L", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "M", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "V", ]))
    class_aromatic_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "F", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "W", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "Y", ]))
    class_hydrophobic_gmiss_count <- sum(class_aliphatic_gmiss_count,class_aromatic_gmiss_count)
    class_positively_charged_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "H", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "K", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "R", ]))
    class_negatively_charged_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "D", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "E", ]))
    class_neutral_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "N", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "Q", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "S", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "T", ]))
    class_polar_gmiss_count <- sum(class_positively_charged_gmiss_count,class_negatively_charged_gmiss_count,class_neutral_gmiss_count)
    class_special_gmiss_count <- sum(nrow(gmiss_prp_class[gmiss_prp_class$refAA == "C", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "G", ]), nrow(gmiss_prp_class[gmiss_prp_class$refAA == "P", ]))
    
    #aliphatic
    fisher_marix <- matrix(c(class_aliphatic_pmiss_count, pmiss_class_total-class_aliphatic_pmiss_count,class_aliphatic_gmiss_count,gmiss_class_total-class_aliphatic_gmiss_count), nrow=2)
    fisher_aliphatic <- fisher.test(fisher_marix)
    
    #aromatic
    fisher_marix <- matrix(c(class_aromatic_pmiss_count, pmiss_class_total-class_aromatic_pmiss_count,class_aromatic_gmiss_count,gmiss_class_total-class_aromatic_gmiss_count), nrow=2)
    fisher_aromatic <- fisher.test(fisher_marix)
    
    #hydrophobic
    fisher_marix <- matrix(c(class_hydrophobic_pmiss_count, pmiss_class_total-class_hydrophobic_pmiss_count,class_hydrophobic_gmiss_count,gmiss_class_total-class_hydrophobic_gmiss_count), nrow=2)
    fisher_hydrophobic <- fisher.test(fisher_marix)
    
    #positively_charged
    fisher_marix <- matrix(c(class_positively_charged_pmiss_count, pmiss_class_total-class_positively_charged_pmiss_count,class_positively_charged_gmiss_count,gmiss_class_total-class_positively_charged_gmiss_count), nrow=2)
    fisher_positively_charged <- fisher.test(fisher_marix)
    
    #negatively_charged
    fisher_marix <- matrix(c(class_negatively_charged_pmiss_count, pmiss_class_total-class_negatively_charged_pmiss_count,class_negatively_charged_gmiss_count,gmiss_class_total-class_negatively_charged_gmiss_count), nrow=2)
    fisher_negatively_charged <- fisher.test(fisher_marix)
    
    #neutral
    fisher_marix <- matrix(c(class_neutral_pmiss_count, pmiss_class_total-class_neutral_pmiss_count,class_neutral_gmiss_count,gmiss_class_total-class_neutral_gmiss_count), nrow=2)
    fisher_neutral <- fisher.test(fisher_marix)
    
    #polar
    fisher_marix <- matrix(c(class_polar_pmiss_count, pmiss_class_total-class_polar_pmiss_count,class_polar_gmiss_count,gmiss_class_total-class_polar_gmiss_count), nrow=2)
    fisher_polar <- fisher.test(fisher_marix)
    
    #special
    fisher_marix <- matrix(c(class_special_pmiss_count, pmiss_class_total-class_special_pmiss_count,class_special_gmiss_count,gmiss_class_total-class_special_gmiss_count), nrow=2)
    fisher_special <- fisher.test(fisher_marix)
    
    #Forest Plot -- comm
    odds_cp8  <- c(fisher_special$estimate, fisher_polar$estimate, fisher_neutral$estimate,fisher_negatively_charged$estimate, fisher_positively_charged$estimate, fisher_hydrophobic$estimate, fisher_aromatic$estimate,fisher_aliphatic$estimate) 
    lower_cp8 <- c(fisher_special$conf.int[1], fisher_polar$conf.int[1], fisher_neutral$conf.int[1], fisher_negatively_charged$conf.int[1], fisher_positively_charged$conf.int[1], fisher_hydrophobic$conf.int[1], fisher_aromatic$conf.int[1],fisher_aliphatic$conf.int[1]) 
    upper_cp8 <- c(fisher_special$conf.int[2], fisher_polar$conf.int[2], fisher_neutral$conf.int[2], fisher_negatively_charged$conf.int[2], fisher_positively_charged$conf.int[2], fisher_hydrophobic$conf.int[2], fisher_aromatic$conf.int[2],fisher_aliphatic$conf.int[2]) 
    p_value_cp8 <- c(fisher_special$p.value, fisher_polar$p.value, fisher_neutral$p.value, fisher_negatively_charged$p.value, fisher_positively_charged$p.value, fisher_hydrophobic$p.value, fisher_aromatic$p.value,fisher_aliphatic$p.value) 
    label_cp8 <- c('1.1.Special', '1.2.Polar', '1.3.Neutral', '2.1.Neg_charged', '2.2.Pos_charged', '2.3.hydrophobic', '3.1.aromatic', '3.2.aliphatic') 
    for_chem_raw <- data.frame(label_cp8, odds_cp8, lower_cp8, upper_cp8, p_value_cp8) 
    for_chem <- setNames(for_chem_raw, c("chem_type", "odds", "lower", "upper", "p_value"))  
    
    pcut = 5.0e-05;
    for_chem$odds <- ifelse(for_chem$odds == 0.0, 0.05, for_chem$odds)
    for_chem$lower <- ifelse(for_chem$odds == 0.0, 0.05, for_chem$lower)
    for_chem$lower <- ifelse(for_chem$lower < 0.05, 0.05, for_chem$lower)
    for_chem$odds <- ifelse(for_chem$odds > 20.0, 20.0, for_chem$odds)
    for_chem$upper <- ifelse(for_chem$odds > 20.0, 20.0, for_chem$upper)
    for_chem$upper <- ifelse(for_chem$upper > 20.0, 20.0, for_chem$upper)
    for_chem$p_value <- ifelse(for_chem$p_value < 1.0e-300, 1.0e-300, for_chem$p_value)
    for_chem$p_value_label <- ifelse(for_chem$p_value < pcut, paste(sprintf("%0.2e", for_chem$p_value), "*)", sep = "("), sprintf("%0.2e", for_chem$p_value))
    
    ggplot(for_chem, aes(x=for_chem$chem_type, y=for_chem$odds, ymin=for_chem$lower, ymax=for_chem$upper))+
      geom_errorbar(aes(ymin = for_chem$lower, ymax = for_chem$upper), width = 0.15, size=1,colour=ifelse(for_chem$p_value < pcut, "black", "darkgrey")) +
      geom_point(size=5, colour=ifelse(for_chem$odds < 1.0, "blue", "red")) +
      geom_text(data = for_chem,aes(x=rep(1.4:8.4,1), y=for_chem$odds), label = for_chem$p_value_label, fontface = "bold", size = 4,hjust="inward") +
      geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
      #scale_shape_identity() +
      coord_flip() +  # flip coordinates (puts labels on y axis)
      #scale_y_log10() +
      theme_bw(base_size = 14) +  # use a white background
      scale_y_log10(breaks=c(0.05,0.1,0.3,1.0,3.0,10.0,20.0), limits=c(0.05,20.0))+
      theme(legend.position = "top", legend.direction = "horizontal") +
      labs(title="Pathogenic vs. Population Variant Enrichment Analysis", 
           subtitle="Odds>1.0 (<1.0) indicates enrichment in pathogenic (populations) variants\np-value cut-off = 5.0e-05",
           x="Physiochemmical properties of amino acids",
           y=expression(paste("Odds (95% CI) "["log10"],""))) +
      scale_x_discrete(labels=c("1.1.Special"="special\n(C,P,G)","1.2.Polar"="polar\n(H,K,R,D,E,N,Q,S,T)","1.3.Neutral"="Neutral\n(N,Q,S,T)","2.1.Neg_charged"="negatively-charged\n(D,E)","2.2.Pos_charged"="positively-charged\n(H,K,R)","2.3.hydrophobic"="hydrophobic\n(F,W,Y,A,I,L,M,V)","3.1.aromatic"="aromatic\n(F,W,Y)","3.2.aliphatic"="aliphatic\n(A,I,L,M,V)"))
    
  })
  
  ## PLOTS -- CHEM --- END ##
  
  
  # Bar Plot of Bonds
  output$BONDbar <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
    }
    
    gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
    gmiss_prp_class_bond <- subset(gmiss_prp_class,gmiss_prp_class$bond_type != 'na')
    gmiss_prp_class_bond_all <- data.frame(separate_rows(gmiss_prp_class_bond,'bond_type'))
    total_gmiss_bond <- nrow(gmiss_prp_class_bond_all)
    
    pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
    pmiss_prp_class_bond <- subset(pmiss_prp_class,pmiss_prp_class$bond_type != 'na')
    pmiss_prp_class_bond_all <- data.frame(separate_rows(pmiss_prp_class_bond,'bond_type'))
    total_pmiss_bond <- nrow(pmiss_prp_class_bond_all)
    
    gmiss_bond_count_val <- data_frame(as.numeric(c(nrow(gmiss_prp_class_bond_all[gmiss_prp_class_bond_all$bond_type == 'D', ]),nrow(gmiss_prp_class_bond_all[gmiss_prp_class_bond_all$bond_type == 'H', ]),nrow(gmiss_prp_class_bond_all[gmiss_prp_class_bond_all$bond_type == 'N', ]),nrow(gmiss_prp_class_bond_all[gmiss_prp_class_bond_all$bond_type == 'S', ]))))
    gmiss_bond_lab <- data_frame(as.character(c("D","H","N","S"))) 
    gmiss_bond_count_raw <- data.frame(cbind(gmiss_bond_lab,gmiss_bond_count_val), stringsAsFactors = FALSE) 
    gmiss_bond_count <- setNames(gmiss_bond_count_raw, c("bond_type", "bondCount")) 
    
    gmiss_bond_count_prop <- (gmiss_bond_count$bondCount / total_gmiss_bond) 
    gmiss_bond_count_label <- data.frame(rep("Population Variant",nrow(gmiss_bond_count))) 
    gmiss_bond_count_freq_rest <- sum(gmiss_bond_count$bondCount) - gmiss_bond_count$bondCount 
    gmiss_bond_count_comp <- cbind(gmiss_bond_count,gmiss_bond_count_prop,gmiss_bond_count_label,gmiss_bond_count_freq_rest) 
    gmiss_bond_count_comp_formatted <- setNames(gmiss_bond_count_comp, c("bond_type", "freq", "prop", "variant","freq_rest")) 
    
    pmiss_bond_count_val <- data_frame(as.numeric(c(nrow(pmiss_prp_class_bond_all[pmiss_prp_class_bond_all$bond_type == 'D', ]),nrow(pmiss_prp_class_bond_all[pmiss_prp_class_bond_all$bond_type == 'H', ]),nrow(pmiss_prp_class_bond_all[pmiss_prp_class_bond_all$bond_type == 'N', ]),nrow(pmiss_prp_class_bond_all[pmiss_prp_class_bond_all$bond_type == 'S', ]))))
    pmiss_bond_lab <- data_frame(as.character(c("D","H","N","S"))) 
    pmiss_bond_count_raw <- data.frame(cbind(pmiss_bond_lab,pmiss_bond_count_val), stringsAsFactors = FALSE) 
    pmiss_bond_count <- setNames(pmiss_bond_count_raw, c("bond_type", "bondCount")) 
    
    pmiss_bond_count_prop <- (pmiss_bond_count$bondCount / total_pmiss_bond) 
    pmiss_bond_count_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_bond_count))) 
    pmiss_bond_count_freq_rest <- sum(pmiss_bond_count$bondCount) - pmiss_bond_count$bondCount 
    pmiss_bond_count_comp <- cbind(pmiss_bond_count,pmiss_bond_count_prop,pmiss_bond_count_label,pmiss_bond_count_freq_rest) 
    pmiss_bond_count_comp_formatted <- setNames(pmiss_bond_count_comp, c("bond_type", "freq", "prop", "variant","freq_rest")) 
    
    bond_proportion <- rbind(gmiss_bond_count_comp_formatted,pmiss_bond_count_comp_formatted) 
    
    ggplot(bond_proportion, aes(bond_proportion$bond_type, bond_proportion$prop)) +   
      geom_bar(aes(fill = bond_proportion$variant), position = "dodge", stat="identity", colour="black", width = 0.7) +
      ggtitle("Distribution of Variants in Amino Acids of Interaction Types") +
      xlab("Protein-protein interaction types") + 
      ylab("Proportion of variants") +
      scale_y_continuous(labels = percent) +
      labs(fill="Variant Type") +
      theme_bw(base_size=14) +
      theme(legend.position = "top", legend.direction = "horizontal") +
      scale_fill_manual(name="", 
                        values = c("Population Variant"="blue", "Pathogenic Variant"="red"))+
      scale_x_discrete(labels=c("D"="Disulfide\nbond","H"="Hydrogen\nbond","N"="Nonbonded Van der\nWaals interaction","S"="Salt-bridge ionic\ninteraction")) +
      theme(axis.text.x = element_text(angle=0, vjust=0.6))
    
  })
  
  #forest plot of bond
  output$BONDforest <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
    }
    
    gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
    gmiss_prp_class_bond <- subset(gmiss_prp_class,gmiss_prp_class$bond_type != 'na')
    gmiss_prp_class_bond_all <- data.frame(separate_rows(gmiss_prp_class_bond,'bond_type'))
    total_gmiss_bond <- nrow(gmiss_prp_class_bond_all)
    
    pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
    pmiss_prp_class_bond <- subset(pmiss_prp_class,pmiss_prp_class$bond_type != 'na')
    pmiss_prp_class_bond_all <- data.frame(separate_rows(pmiss_prp_class_bond,'bond_type'))
    total_pmiss_bond <- nrow(pmiss_prp_class_bond_all)
    
    #pmiss counts for fischer test
    class_disulfide_pmiss_count <- nrow(pmiss_prp_class_bond_all[pmiss_prp_class_bond_all$bond_type == "D", ])
    class_saltbridge_pmiss_count <- nrow(pmiss_prp_class_bond_all[pmiss_prp_class_bond_all$bond_type == "S", ])
    class_hydrogen_pmiss_count <- nrow(pmiss_prp_class_bond_all[pmiss_prp_class_bond_all$bond_type == "H", ])
    class_nonbonded_pmiss_count <- nrow(pmiss_prp_class_bond_all[pmiss_prp_class_bond_all$bond_type == "N", ])
    
    #gmiss counts for fischer's test
    class_disulfide_gmiss_count <- nrow(gmiss_prp_class_bond_all[gmiss_prp_class_bond_all$bond_type == "D", ])
    class_saltbridge_gmiss_count <- nrow(gmiss_prp_class_bond_all[gmiss_prp_class_bond_all$bond_type == "S", ])
    class_hydrogen_gmiss_count <- nrow(gmiss_prp_class_bond_all[gmiss_prp_class_bond_all$bond_type == "H", ])
    class_nonbonded_gmiss_count <- nrow(gmiss_prp_class_bond_all[gmiss_prp_class_bond_all$bond_type == "N", ])
    
    #disulfide
    fisher_marix <- matrix(c(class_disulfide_pmiss_count, total_pmiss_bond-class_disulfide_pmiss_count,class_disulfide_gmiss_count,total_gmiss_bond-class_disulfide_gmiss_count), nrow=2)
    fisher_disulfide <- fisher.test(fisher_marix)
    
    #saltbridge
    fisher_marix <- matrix(c(class_saltbridge_pmiss_count, total_pmiss_bond-class_saltbridge_pmiss_count,class_saltbridge_gmiss_count,total_gmiss_bond-class_saltbridge_gmiss_count), nrow=2)
    fisher_saltbridge <- fisher.test(fisher_marix)
    
    #hydrogen
    fisher_marix <- matrix(c(class_hydrogen_pmiss_count, total_pmiss_bond-class_hydrogen_pmiss_count,class_hydrogen_gmiss_count,total_gmiss_bond-class_hydrogen_gmiss_count), nrow=2)
    fisher_hydrogen <- fisher.test(fisher_marix)
    
    #nonbonded
    fisher_marix <- matrix(c(class_nonbonded_pmiss_count, total_pmiss_bond-class_nonbonded_pmiss_count,class_nonbonded_gmiss_count,total_gmiss_bond-class_nonbonded_gmiss_count), nrow=2)
    fisher_nonbonded <- fisher.test(fisher_marix)
    
    #Forest Plot -- comm
    odds_ib4  <- c(fisher_nonbonded$estimate, fisher_hydrogen$estimate, fisher_saltbridge$estimate,fisher_disulfide$estimate) 
    lower_ib4 <- c(fisher_nonbonded$conf.int[1], fisher_hydrogen$conf.int[1], fisher_saltbridge$conf.int[1],fisher_disulfide$conf.int[1]) 
    upper_ib4 <- c(fisher_nonbonded$conf.int[2], fisher_hydrogen$conf.int[2], fisher_saltbridge$conf.int[2],fisher_disulfide$conf.int[2]) 
    p_value_ib4 <- c(fisher_nonbonded$p.value, fisher_hydrogen$p.value, fisher_saltbridge$p.value,fisher_disulfide$p.value) 
    label_ib4 <- c('1.1.Nonbonded', '1.2.Hydrogen', '1.3.Saltbridge', '2.1.Disulfide') 
    for_bond_raw <- data.frame(label_ib4, odds_ib4, lower_ib4, upper_ib4, p_value_ib4) 
    for_bond <- setNames(for_bond_raw, c("bond_type", "odds", "lower", "upper", "p_value"))  
    
    pcut = 5.0e-05;
    for_bond$odds <- ifelse(for_bond$odds == 0.0, 0.05, for_bond$odds)
    for_bond$lower <- ifelse(for_bond$odds == 0.0, 0.05, for_bond$lower)
    for_bond$lower <- ifelse(for_bond$lower < 0.05, 0.05, for_bond$lower)
    for_bond$odds <- ifelse(for_bond$odds > 20.0, 20.0, for_bond$odds)
    for_bond$upper <- ifelse(for_bond$odds > 20.0, 20.0, for_bond$upper)
    for_bond$upper <- ifelse(for_bond$upper > 20.0, 20.0, for_bond$upper)
    for_bond$p_value <- ifelse(for_bond$p_value < 1.0e-300, 1.0e-300, for_bond$p_value)
    for_bond$p_value_label <- ifelse(for_bond$p_value < pcut, paste(sprintf("%0.2e", for_bond$p_value), "*)", sep = "("), sprintf("%0.2e", for_bond$p_value))
    
    ggplot(for_bond, aes(x=for_bond$bond_type, y=for_bond$odds, ymin=for_bond$lower, ymax=for_bond$upper))+
      geom_errorbar(aes(ymin = for_bond$lower, ymax = for_bond$upper), width = 0.15, size=1,colour=ifelse(for_bond$p_value < pcut, "black", "darkgrey")) +
      geom_point(size=5, colour=ifelse(for_bond$odds < 1.0, "blue", "red")) +
      geom_text(data = for_bond,aes(x=rep(1.4:4.4,1), y=for_bond$odds), label = for_bond$p_value_label, fontface = "bold", size = 4,hjust="inward") +
      geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
      coord_flip() +  # flip coordinates (puts labels on y axis)
      theme_bw(base_size = 14) +  # use a white background
      scale_y_log10(breaks=c(0.05,0.1,0.3,1.0,3.0,10.0,20.0), limits=c(0.05,20.0))+
      theme(legend.position = "top", legend.direction = "horizontal") +
      labs(title="Pathogenic vs. Population Variant Enrichment Analysis", 
           subtitle="Odds>1.0 (<1.0) indicates enrichment in pathogenic (populations) variants\np-value cut-off = 5.0e-05",
           x="Protein-protein interaction types",
           y=expression(paste("Odds (95% CI) "["log10"],""))) +
      scale_x_discrete(labels=c("1.1.Nonbonded"="nonbonded Van der\nWaals interaction","1.2.Hydrogen"="Hydrogen\nbond","1.3.Saltbridge"="Salt-bridge ionic\ninteraction","2.1.Disulfide"="Disulfide\nbond"))
    
  })
  
  #PTMforest
  output$PTMforest <- renderPlot({
    #methylation
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    
    class_acetyl_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dist_acet < 10.0, ])
    class_methyl_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dist_meth < 10.0, ])
    class_gclnac_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dist_glcn < 10.0, ])
    class_phosphoryl_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dist_phos < 10.0, ])
    class_sumoyl_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dist_sumo < 10.0, ])
    class_ubiquitin_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$dist_ubiq < 10.0, ])
    
    class_acetyl_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dist_acet < 10.0, ])
    class_methyl_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dist_meth < 10.0, ])
    class_gclnac_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dist_glcn < 10.0, ])
    class_phosphoryl_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dist_phos < 10.0, ])
    class_sumoyl_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dist_sumo < 10.0, ])
    class_ubiquitin_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$dist_ubiq < 10.0, ])
    
    #acetyl
    fisher_marix <- matrix(c(class_acetyl_pmiss_count, pmiss_class_total-class_acetyl_pmiss_count,class_acetyl_gmiss_count,gmiss_class_total-class_acetyl_gmiss_count), nrow=2)
    fisher_acetyl <- fisher.test(fisher_marix)
    
    #methyl
    fisher_marix <- matrix(c(class_methyl_pmiss_count, pmiss_class_total-class_methyl_pmiss_count,class_methyl_gmiss_count,gmiss_class_total-class_methyl_gmiss_count), nrow=2)
    fisher_methyl <- fisher.test(fisher_marix)
    
    #gclnac
    fisher_marix <- matrix(c(class_gclnac_pmiss_count, pmiss_class_total-class_gclnac_pmiss_count,class_gclnac_gmiss_count,gmiss_class_total-class_gclnac_gmiss_count), nrow=2)
    fisher_gclnac <- fisher.test(fisher_marix)
    
    #phosphoryl
    fisher_marix <- matrix(c(class_phosphoryl_pmiss_count, pmiss_class_total-class_phosphoryl_pmiss_count,class_phosphoryl_gmiss_count,gmiss_class_total-class_phosphoryl_gmiss_count), nrow=2)
    fisher_phosphoryl <- fisher.test(fisher_marix)
    
    #sumoyl
    fisher_marix <- matrix(c(class_sumoyl_pmiss_count, pmiss_class_total-class_sumoyl_pmiss_count,class_sumoyl_gmiss_count,gmiss_class_total-class_sumoyl_gmiss_count), nrow=2)
    fisher_sumoyl <- fisher.test(fisher_marix)
    
    #ubiquitin
    fisher_marix <- matrix(c(class_ubiquitin_pmiss_count, pmiss_class_total-class_ubiquitin_pmiss_count,class_ubiquitin_gmiss_count,gmiss_class_total-class_ubiquitin_gmiss_count), nrow=2)
    fisher_ubiquitin <- fisher.test(fisher_marix)
    
    
    #Forest Plot -- comm
    odds_ptm6  <- c(fisher_ubiquitin$estimate,fisher_gclnac$estimate, fisher_phosphoryl$estimate, fisher_sumoyl$estimate, fisher_acetyl$estimate,fisher_methyl$estimate)
    lower_ptm6 <- c(fisher_ubiquitin$conf.int[1], fisher_gclnac$conf.int[1], fisher_phosphoryl$conf.int[1], fisher_sumoyl$conf.int[1], fisher_acetyl$conf.int[1],fisher_methyl$conf.int[1])
    upper_ptm6 <- c(fisher_ubiquitin$conf.int[2], fisher_gclnac$conf.int[2], fisher_phosphoryl$conf.int[2], fisher_sumoyl$conf.int[2], fisher_acetyl$conf.int[2],fisher_methyl$conf.int[2])
    p_value_ptm6 <- c(fisher_ubiquitin$p.value, fisher_gclnac$p.value, fisher_phosphoryl$p.value, fisher_sumoyl$p.value, fisher_acetyl$p.value,fisher_methyl$p.value)
    label_ptm6 <- c('1.1.ubiquitin', '1.2.gclnac', '1.3.phosphoryl', '2.1.sumoyl', '2.2.acetyl', '2.3.methyl')
    for_ptm_raw <- data.frame(label_ptm6, odds_ptm6, lower_ptm6, upper_ptm6, p_value_ptm6)
    for_ptm <- setNames(for_ptm_raw, c("ptm_type", "odds", "lower", "upper", "p_value"))
    
    pcut = 5.0e-05;
    for_ptm$odds <- ifelse(for_ptm$odds == 0.0, 0.05, for_ptm$odds)
    for_ptm$lower <- ifelse(for_ptm$odds == 0.0, 0.05, for_ptm$lower)
    for_ptm$lower <- ifelse(for_ptm$lower < 0.05, 0.05, for_ptm$lower)
    for_ptm$odds <- ifelse(for_ptm$odds > 20.0, 20.0, for_ptm$odds)
    for_ptm$upper <- ifelse(for_ptm$odds > 20.0, 20.0, for_ptm$upper)
    for_ptm$upper <- ifelse(for_ptm$upper > 20.0, 20.0, for_ptm$upper)
    for_ptm$p_value <- ifelse(for_ptm$p_value < 1.0e-300, 1.0e-300, for_ptm$p_value)
    for_ptm$p_value_label <- ifelse(for_ptm$p_value < pcut, paste(sprintf("%0.2e", for_ptm$p_value), "*)", sep = "("), sprintf("%0.2e", for_ptm$p_value))
    
    ggplot(for_ptm, aes(x=for_ptm$ptm_type, y=for_ptm$odds, ymin=for_ptm$lower, ymax=for_ptm$upper))+
      geom_errorbar(aes(ymin = for_ptm$lower, ymax = for_ptm$upper), width = 0.15, size=1,colour=ifelse(for_ptm$p_value < pcut, "black", "darkgrey")) +
      geom_point(size=5, colour=ifelse(for_ptm$odds < 1.0, "blue", "red")) +
      geom_text(data = for_ptm,aes(x=rep(1.4:6.4,1), y=for_ptm$odds), label = for_ptm$p_value_label, fontface = "bold", size = 4,hjust="inward") +
      geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
      #scale_shape_identity() +
      coord_flip() +  # flip coordinates (puts labels on y axis)
      #scale_y_log10() +
      theme_bw(base_size = 14) +  # use a white background
      scale_y_log10(breaks=c(0.05,0.1,0.3,1.0,3.0,10.0,20.0), limits=c(0.05,20.0))+
      theme(legend.position = "top", legend.direction = "horizontal") +
      labs(title=expression(paste("Pathogenic vs. Population Variant ( near < 10 ",ring(A), " to PTMs) Enrichment")),
           subtitle="Odds>1.0 (<1.0) indicates enrichment in pathogenic (populations) variants\np-value cut-off = 5.0e-05",
           x="Post-translational modification types",
           y=expression(paste("Odds (95% CI) "["log10"],""))) +
      scale_x_discrete(labels=c("1.1.ubiquitin"="ubiquitination","1.2.gclnac"="O.GclNAc","1.3.phosphoryl"="phosphorylation","2.1.sumoyl"="sumoylation","2.2.acetyl"="acetylation","2.3.methyl"="methylation"))
    
  })
  
  #PTMbox
  output$PTMbox <- renderPlot({
    #methylation
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
      
      gmiss_prp_class <- subset(gmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      gmiss_class_total <- nrow(gmiss_prp_class)
      pmiss_prp_class <- subset(pmiss_dssp_prp, geneName %in% class_specific_gene$genes)
      pmiss_class_total <- nrow(pmiss_prp_class)
    }
    
    ptm_all <- data.frame(dist = as.numeric(), vartype = as.character(), ptmtype = as.character())
    
    #acetylation
    gmiss_prp_class_acetyl <- subset(gmiss_prp_class,gmiss_prp_class$dist_acet < 100)
    gmiss_acetyl_raw <- data.frame(gmiss_prp_class_acetyl$dist_acet)
    gmiss_acetyl_raw_label <- data.frame(rep("Population Variant",nrow(gmiss_acetyl_raw)))
    gmiss_acetyl_raw_comb <- cbind(gmiss_acetyl_raw, gmiss_acetyl_raw_label)
    gmiss_acetyl_formatted <- setNames(gmiss_acetyl_raw_comb, c("dist", "vartype"))
    gmiss_acetyl_count <- nrow(gmiss_acetyl_formatted)
    
    pmiss_prp_class_acetyl <- subset(pmiss_prp_class,pmiss_prp_class$dist_acet < 100)
    pmiss_acetyl_raw <- data.frame(pmiss_prp_class_acetyl$dist_acet)
    pmiss_acetyl_raw_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_acetyl_raw)))
    pmiss_acetyl_raw_comb <- cbind(pmiss_acetyl_raw, pmiss_acetyl_raw_label)
    pmiss_acetyl_formatted <- setNames(pmiss_acetyl_raw_comb, c("dist", "vartype"))
    pmiss_acetyl_count <- nrow(pmiss_acetyl_formatted)
    
    if ((gmiss_acetyl_count > 0) | (pmiss_acetyl_count > 0)){
      acetyl_combined <- rbind(gmiss_acetyl_formatted,pmiss_acetyl_formatted)
      acetyl_combined$ptmtype <- "Acetylation"
      
      ptm_all <- rbind(ptm_all, acetyl_combined)
    }
    
    #methylation
    gmiss_prp_class_methyl <- subset(gmiss_prp_class,gmiss_prp_class$dist_meth < 100)
    gmiss_methyl_raw <- data.frame(gmiss_prp_class_methyl$dist_meth)
    gmiss_methyl_raw_label <- data.frame(rep("Population Variant",nrow(gmiss_methyl_raw)))
    gmiss_methyl_raw_comb <- cbind(gmiss_methyl_raw, gmiss_methyl_raw_label)
    gmiss_methyl_formatted <- setNames(gmiss_methyl_raw_comb, c("dist", "vartype"))
    gmiss_methyl_count <- nrow(gmiss_methyl_formatted)
    
    pmiss_prp_class_methyl <- subset(pmiss_prp_class,pmiss_prp_class$dist_meth < 100)
    pmiss_methyl_raw <- data.frame(pmiss_prp_class_methyl$dist_meth)
    pmiss_methyl_raw_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_methyl_raw)))
    pmiss_methyl_raw_comb <- cbind(pmiss_methyl_raw, pmiss_methyl_raw_label)
    pmiss_methyl_formatted <- setNames(pmiss_methyl_raw_comb, c("dist", "vartype"))
    pmiss_methyl_count <- nrow(pmiss_methyl_formatted)
    
    if ((gmiss_methyl_count > 0) | (pmiss_methyl_count > 0)){
      methyl_combined <- rbind(gmiss_methyl_formatted,pmiss_methyl_formatted)
      methyl_combined$ptmtype <- "Methylation"
      
      ptm_all <- rbind(ptm_all, methyl_combined)
    }
    
    #O.GclNAc
    gmiss_prp_class_glcn <- subset(gmiss_prp_class,gmiss_prp_class$dist_glcn < 100)
    gmiss_glcn_raw <- data.frame(gmiss_prp_class_glcn$dist_glcn)
    gmiss_glcn_raw_label <- data.frame(rep("Population Variant",nrow(gmiss_glcn_raw)))
    gmiss_glcn_raw_comb <- cbind(gmiss_glcn_raw, gmiss_glcn_raw_label)
    gmiss_glcn_formatted <- setNames(gmiss_glcn_raw_comb, c("dist", "vartype"))
    gmiss_glcn_count <- nrow(gmiss_glcn_formatted)
    
    pmiss_prp_class_glcn <- subset(pmiss_prp_class,pmiss_prp_class$dist_glcn < 100)
    pmiss_glcn_raw <- data.frame(pmiss_prp_class_glcn$dist_glcn)
    pmiss_glcn_raw_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_glcn_raw)))
    pmiss_glcn_raw_comb <- cbind(pmiss_glcn_raw, pmiss_glcn_raw_label)
    pmiss_glcn_formatted <- setNames(pmiss_glcn_raw_comb, c("dist", "vartype"))
    pmiss_glcn_count <- nrow(pmiss_glcn_formatted)
    
    if ((gmiss_glcn_count > 0) | (pmiss_glcn_count > 0)){
      glcn_combined <- rbind(gmiss_glcn_formatted,pmiss_glcn_formatted)
      glcn_combined$ptmtype <- "O.GlcNAc"
      
      ptm_all <- rbind(ptm_all, glcn_combined)
    }
    
    #phosphorylation
    gmiss_prp_class_phos <- subset(gmiss_prp_class,gmiss_prp_class$dist_phos < 100)
    gmiss_phos_raw <- data.frame(gmiss_prp_class_phos$dist_phos)
    gmiss_phos_raw_label <- data.frame(rep("Population Variant",nrow(gmiss_phos_raw)))
    gmiss_phos_raw_comb <- cbind(gmiss_phos_raw, gmiss_phos_raw_label)
    gmiss_phos_formatted <- setNames(gmiss_phos_raw_comb, c("dist", "vartype"))
    gmiss_phos_count <- nrow(gmiss_phos_formatted)
    
    pmiss_prp_class_phos <- subset(pmiss_prp_class,pmiss_prp_class$dist_phos < 100)
    pmiss_phos_raw <- data.frame(pmiss_prp_class_phos$dist_phos)
    pmiss_phos_raw_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_phos_raw)))
    pmiss_phos_raw_comb <- cbind(pmiss_phos_raw, pmiss_phos_raw_label)
    pmiss_phos_formatted <- setNames(pmiss_phos_raw_comb, c("dist", "vartype"))
    pmiss_phos_count <- nrow(pmiss_phos_formatted)
    
    if ((gmiss_phos_count > 0) | (pmiss_phos_count > 0)){
      phos_combined <- rbind(gmiss_phos_formatted,pmiss_phos_formatted)
      phos_combined$ptmtype <- "Phosphorylation"
      
      ptm_all <- rbind(ptm_all, phos_combined)
    }
    
    #sumoylation
    gmiss_prp_class_sumo <- subset(gmiss_prp_class,gmiss_prp_class$dist_sumo < 100)
    gmiss_sumo_raw <- data.frame(gmiss_prp_class_sumo$dist_sumo)
    gmiss_sumo_raw_label <- data.frame(rep("Population Variant",nrow(gmiss_sumo_raw)))
    gmiss_sumo_raw_comb <- cbind(gmiss_sumo_raw, gmiss_sumo_raw_label)
    gmiss_sumo_formatted <- setNames(gmiss_sumo_raw_comb, c("dist", "vartype"))
    gmiss_sumo_count <- nrow(gmiss_sumo_formatted)
    
    pmiss_prp_class_sumo <- subset(pmiss_prp_class,pmiss_prp_class$dist_sumo < 100)
    pmiss_sumo_raw <- data.frame(pmiss_prp_class_sumo$dist_sumo)
    pmiss_sumo_raw_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_sumo_raw)))
    pmiss_sumo_raw_comb <- cbind(pmiss_sumo_raw, pmiss_sumo_raw_label)
    pmiss_sumo_formatted <- setNames(pmiss_sumo_raw_comb, c("dist", "vartype"))
    pmiss_sumo_count <- nrow(pmiss_sumo_formatted)
    
    if ((gmiss_sumo_count > 0) | (pmiss_sumo_count > 0)){
      sumo_combined <- rbind(gmiss_sumo_formatted,pmiss_sumo_formatted)
      sumo_combined$ptmtype <- "Sumoylation"
      
      ptm_all <- rbind(ptm_all, sumo_combined)
    }
    
    #ubiquitination
    gmiss_prp_class_ubiq <- subset(gmiss_prp_class,gmiss_prp_class$dist_ubiq < 100)
    gmiss_ubiq_raw <- data.frame(gmiss_prp_class_ubiq$dist_ubiq)
    gmiss_ubiq_raw_label <- data.frame(rep("Population Variant",nrow(gmiss_ubiq_raw)))
    gmiss_ubiq_raw_comb <- cbind(gmiss_ubiq_raw, gmiss_ubiq_raw_label)
    gmiss_ubiq_formatted <- setNames(gmiss_ubiq_raw_comb, c("dist", "vartype"))
    gmiss_ubiq_count <- nrow(gmiss_ubiq_formatted)
    
    pmiss_prp_class_ubiq <- subset(pmiss_prp_class,pmiss_prp_class$dist_ubiq < 100)
    pmiss_ubiq_raw <- data.frame(pmiss_prp_class_ubiq$dist_ubiq)
    pmiss_ubiq_raw_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_ubiq_raw)))
    pmiss_ubiq_raw_comb <- cbind(pmiss_ubiq_raw, pmiss_ubiq_raw_label)
    pmiss_ubiq_formatted <- setNames(pmiss_ubiq_raw_comb, c("dist", "vartype"))
    pmiss_ubiq_count <- nrow(pmiss_ubiq_formatted)
    
    if ((gmiss_ubiq_count > 0) | (pmiss_ubiq_count > 0)){
      ubiq_combined <- rbind(gmiss_ubiq_formatted,pmiss_ubiq_formatted)
      ubiq_combined$ptmtype <- "Ubiquitination"
      
      ptm_all <- rbind(ptm_all, ubiq_combined)
    }
    
    dodge <- position_dodge(width = T)
    my_comparisons <- list( c("Population Variant", "Pathogenic Variant"))
    
    ggplot(ptm_all, aes(ptm_all$vartype, ptm_all$dist)) +
      geom_violin(position = dodge) +
      geom_boxplot(aes(fill=factor(ptm_all$vartype)),width=.1, position = dodge) +
      theme_bw(base_size = 14) +
      theme(legend.position = "top", legend.direction = "horizontal") +
      labs(title=expression(paste("Distribution of Distances ( < 100 ",ring(A), ") to PTM sites from Variant Location in 3D")),
           x="Variant types",
           y=expression(paste("Spatial distance to PTM sites (",ring(A), ")")),
           fill="Variant Type") +
      scale_fill_manual(name="",
                        values = c("Population Variant"="blue", "Pathogenic Variant"="red")) +
      facet_wrap( ~ ptm_all$ptmtype, ncol=2)
  })
  
  
  
  # Bar Plot of Bonds
  output$upFuncBar <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
    }
    
    gmiss_prp_class <- subset(gmiss_dssp_uniprot, geneName %in% class_specific_gene$genes)
    gmiss_class_total <- nrow(gmiss_prp_class)
    pmiss_prp_class <- subset(pmiss_dssp_uniprot, geneName %in% class_specific_gene$genes)
    pmiss_class_total <- nrow(pmiss_prp_class)
    
    gmiss_uniprot_count_val <- data_frame(as.numeric(c(nrow(gmiss_prp_class[gmiss_prp_class$functional_site == 1, ]), nrow(gmiss_prp_class[gmiss_prp_class$functional_region == 1, ]),nrow(gmiss_prp_class[gmiss_prp_class$sequence_motif_region == 1, ]), nrow(gmiss_prp_class[gmiss_prp_class$domain == 1, ]), nrow(gmiss_prp_class[gmiss_prp_class$molecular_processing == 1, ]),nrow(gmiss_prp_class[gmiss_prp_class$modified_aa == 1, ]))))
    gmiss_uniprot_lab <- data_frame(as.character(c("funcsite","bindreg","seqmotif","moddom","molproc","modres")))
    gmiss_uniprot_count_raw <- data.frame(cbind(gmiss_uniprot_lab,gmiss_uniprot_count_val), stringsAsFactors = FALSE)
    gmiss_uniprot_count <- setNames(gmiss_uniprot_count_raw, c("uniprot", "aacount"))
    
    gmiss_uniprot_count_prop <- gmiss_uniprot_count$aacount / gmiss_class_total 
    gmiss_uniprot_count_label <- data.frame(rep("Population Variant",nrow(gmiss_uniprot_count))) 
    gmiss_uniprot_count_freq_rest <- sum(gmiss_uniprot_count$aacount) - gmiss_uniprot_count$aacount 
    gmiss_uniprot_count_comp <- cbind(gmiss_uniprot_count,gmiss_uniprot_count_prop,gmiss_uniprot_count_label,gmiss_uniprot_count_freq_rest) 
    gmiss_uniprot_count_comp_formatted <- setNames(gmiss_uniprot_count_comp, c("uniprot", "freq", "prop", "variant","freq_rest")) 
    
    pmiss_uniprot_count_val <- data_frame(as.numeric(c(nrow(pmiss_prp_class[pmiss_prp_class$functional_site == 1, ]), nrow(pmiss_prp_class[pmiss_prp_class$functional_region == 1, ]),nrow(pmiss_prp_class[pmiss_prp_class$sequence_motif_region == 1, ]), nrow(pmiss_prp_class[pmiss_prp_class$domain == 1, ]), nrow(pmiss_prp_class[pmiss_prp_class$molecular_processing == 1, ]),nrow(pmiss_prp_class[pmiss_prp_class$modified_aa == 1, ]))))
    pmiss_uniprot_lab <- data_frame(as.character(c("funcsite","bindreg","seqmotif","moddom","molproc","modres")))
    pmiss_uniprot_count_raw <- data.frame(cbind(pmiss_uniprot_lab,pmiss_uniprot_count_val), stringsAsFactors = FALSE)
    pmiss_uniprot_count <- setNames(pmiss_uniprot_count_raw, c("uniprot", "aacount"))
    
    pmiss_uniprot_count_prop <- pmiss_uniprot_count$aacount / pmiss_class_total 
    pmiss_uniprot_count_label <- data.frame(rep("Pathogenic Variant",nrow(pmiss_uniprot_count))) 
    pmiss_uniprot_count_freq_rest <- sum(pmiss_uniprot_count$aacount) - pmiss_uniprot_count$aacount 
    pmiss_uniprot_count_comp <- cbind(pmiss_uniprot_count,pmiss_uniprot_count_prop,pmiss_uniprot_count_label,pmiss_uniprot_count_freq_rest) 
    pmiss_uniprot_count_comp_formatted <- setNames(pmiss_uniprot_count_comp, c("uniprot", "freq", "prop", "variant","freq_rest")) 
    
    uniprot_proportion <- rbind(gmiss_uniprot_count_comp_formatted,pmiss_uniprot_count_comp_formatted)
    
    ggplot(uniprot_proportion, aes(uniprot_proportion$uniprot, uniprot_proportion$prop)) +   
      geom_bar(aes(fill = uniprot_proportion$variant), position = "dodge", stat="identity", colour="black", width = 0.7) +
      ggtitle("Distribution of Variants in Functional Features") +
      xlab("Functional feature types") + 
      ylab("Proportion of variants") +
      scale_y_continuous(labels = percent) +
      labs(fill="Variant Type") +
      theme_bw(base_size=14) +
      theme(legend.position = "top", legend.direction = "horizontal") +
      scale_fill_manual(name="", 
                        values = c("Population Variant"="blue", "Pathogenic Variant"="red"))+
      scale_x_discrete(labels=c("funcsite"="Functional\nsite","bindreg"="Functional/binding\nregion","seqmotif"="Sequence\nmotif/region","moddom"="Modular\ndomain","molproc"="Molecular processing\nassociated region","modres"="Modified\nresidue"))
    
  })
  
  #forest plot of bond
  output$upFuncForest <- renderPlot({
    if(input$homeSideBarTabSetPanel == 'Protein Class' && input$pclassNameselected != ''){
      class_specific_info <- protein_class_def[protein_class_def$superclassName == input$pclassNameselected,]
      class_specific_gene <- data.frame(unique(data.frame(separate_rows(class_specific_info, 'genes'))))
    }
    
    gmiss_prp_class <- subset(gmiss_dssp_uniprot, geneName %in% class_specific_gene$genes)
    gmiss_class_total <- nrow(gmiss_prp_class)
    pmiss_prp_class <- subset(pmiss_dssp_uniprot, geneName %in% class_specific_gene$genes)
    pmiss_class_total <- nrow(pmiss_prp_class)
    
    #pmiss counts for fischer test
    class_funcsite_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$functional_site == 1, ])
    class_bindreg_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$functional_region == 1, ])
    class_seqmotif_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$sequence_motif_region == 1, ])
    class_moddom_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$domain == 1, ])
    class_molproc_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$molecular_processing == 1, ])
    class_modres_pmiss_count <- nrow(pmiss_prp_class[pmiss_prp_class$modified_aa == 1, ])
    
    #gmiss counts for fischer's test
    class_funcsite_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$functional_site == 1, ])
    class_bindreg_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$functional_region == 1, ])
    class_seqmotif_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$sequence_motif_region == 1, ])
    class_moddom_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$domain == 1, ])
    class_molproc_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$molecular_processing == 1, ])
    class_modres_gmiss_count <- nrow(gmiss_prp_class[gmiss_prp_class$modified_aa == 1, ])
    
    #funcsite
    fisher_marix <- matrix(c(class_funcsite_pmiss_count, pmiss_class_total-class_funcsite_pmiss_count,class_funcsite_gmiss_count,gmiss_class_total-class_funcsite_gmiss_count), nrow=2)
    fisher_funcsite <- fisher.test(fisher_marix)
    
    #bindreg
    fisher_marix <- matrix(c(class_bindreg_pmiss_count, pmiss_class_total-class_bindreg_pmiss_count,class_bindreg_gmiss_count,gmiss_class_total-class_bindreg_gmiss_count), nrow=2)
    fisher_bindreg <- fisher.test(fisher_marix)
    
    #seqmotif
    fisher_marix <- matrix(c(class_seqmotif_pmiss_count, pmiss_class_total-class_seqmotif_pmiss_count,class_seqmotif_gmiss_count,gmiss_class_total-class_seqmotif_gmiss_count), nrow=2)
    fisher_seqmotif <- fisher.test(fisher_marix)
    
    #moddom
    fisher_marix <- matrix(c(class_moddom_pmiss_count, pmiss_class_total-class_moddom_pmiss_count,class_moddom_gmiss_count,gmiss_class_total-class_moddom_gmiss_count), nrow=2)
    fisher_moddom <- fisher.test(fisher_marix)
    
    #molproc
    fisher_marix <- matrix(c(class_molproc_pmiss_count, pmiss_class_total-class_molproc_pmiss_count,class_molproc_gmiss_count,gmiss_class_total-class_molproc_gmiss_count), nrow=2)
    fisher_molproc <- fisher.test(fisher_marix)
    
    #modres
    fisher_marix <- matrix(c(class_modres_pmiss_count, pmiss_class_total-class_modres_pmiss_count,class_modres_gmiss_count,gmiss_class_total-class_modres_gmiss_count), nrow=2)
    fisher_modres <- fisher.test(fisher_marix)
    
    
    #Forest Plot -- comm
    odds_up6  <- c(fisher_modres$estimate,fisher_seqmotif$estimate, fisher_moddom$estimate, fisher_molproc$estimate, fisher_funcsite$estimate,fisher_bindreg$estimate) 
    lower_up6 <- c(fisher_modres$conf.int[1], fisher_seqmotif$conf.int[1], fisher_moddom$conf.int[1], fisher_molproc$conf.int[1], fisher_funcsite$conf.int[1],fisher_bindreg$conf.int[1]) 
    upper_up6 <- c(fisher_modres$conf.int[2], fisher_seqmotif$conf.int[2], fisher_moddom$conf.int[2], fisher_molproc$conf.int[2], fisher_funcsite$conf.int[2],fisher_bindreg$conf.int[2]) 
    p_value_up6 <- c(fisher_modres$p.value, fisher_seqmotif$p.value, fisher_moddom$p.value, fisher_molproc$p.value, fisher_funcsite$p.value,fisher_bindreg$p.value) 
    label_up6 <- c('1.1.modres', '1.2.seqmotif', '1.3.moddom', '2.1.molproc', '2.2.funcsite', '2.3.bindreg') 
    for_uniprot_raw <- data.frame(label_up6, odds_up6, lower_up6, upper_up6, p_value_up6) 
    for_uniprot <- setNames(for_uniprot_raw, c("uniprot_type", "odds", "lower", "upper", "p_value"))  
    
    pcut = 5.0e-05;
    for_uniprot$odds <- ifelse(for_uniprot$odds == 0.0, 0.05, for_uniprot$odds)
    for_uniprot$lower <- ifelse(for_uniprot$odds == 0.0, 0.05, for_uniprot$lower)
    for_uniprot$lower <- ifelse(for_uniprot$lower < 0.05, 0.05, for_uniprot$lower)
    for_uniprot$odds <- ifelse(for_uniprot$odds > 20.0, 20.0, for_uniprot$odds)
    for_uniprot$upper <- ifelse(for_uniprot$odds > 20.0, 20.0, for_uniprot$upper)
    for_uniprot$upper <- ifelse(for_uniprot$upper > 20.0, 20.0, for_uniprot$upper)
    for_uniprot$p_value <- ifelse(for_uniprot$p_value < 1.0e-300, 1.0e-300, for_uniprot$p_value)
    for_uniprot$p_value_label <- ifelse(for_uniprot$p_value < pcut, paste(sprintf("%0.2e", for_uniprot$p_value), "*)", sep = "("), sprintf("%0.2e", for_uniprot$p_value))
    
    ggplot(for_uniprot, aes(x=for_uniprot$uniprot_type, y=for_uniprot$odds, ymin=for_uniprot$lower, ymax=for_uniprot$upper))+
      geom_errorbar(aes(ymin = for_uniprot$lower, ymax = for_uniprot$upper), width = 0.15, size=1,colour=ifelse(for_uniprot$p_value < pcut, "black", "darkgrey")) +
      geom_point(size=5, colour=ifelse(for_uniprot$odds < 1.0, "blue", "red")) +
      geom_text(data = for_uniprot,aes(x=rep(1.4:6.4,1), y=for_uniprot$odds), label = for_uniprot$p_value_label, fontface = "bold", size = 4,hjust="inward") +
      geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
      coord_flip() +  # flip coordinates (puts labels on y axis)
      theme_bw(base_size = 14) +  # use a white background
      scale_y_log10(breaks=c(0.05,0.1,0.3,1.0,3.0,10.0,20.0), limits=c(0.05,20.0))+
      theme(legend.position = "top", legend.direction = "horizontal") +
      labs(title="Pathogenic vs. Population Variant Enrichment Analysis", 
           subtitle="Odds>1.0 (<1.0) indicates enrichment in pathogenic (populations) variants\np-value cut-off = 5.0e-05",
           x="UniProt-based functional feature types",
           y=expression(paste("Odds (95% CI) "["log10"],""))) +
      scale_x_discrete(labels=c("1.1.modres"="modified\nresidue","1.2.seqmotif"="sequence\nmotif/region","1.3.moddom"="modular\ndomain","2.1.molproc"="molecular processing\nassociated region","2.2.funcsite"="functional\nsite","2.3.bindreg"="functional/binding\nregion"))
    
    
  })
  
  output$molart <- renderUI({
    geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
    uniprotId <- geneSpecificpmissOrgmiss$uniprotAC
    
    gene_wise_info=read_delim(paste("gene_wise_info/",input$geneSelected,".txt",sep=''), "\t", escape_double = FALSE, trim_ws = TRUE)
    shinyjs::js$runMolart(uniprotId=uniprotId, gene_wise_info=gene_wise_info, containerId=molartContainerId)
  })
  
  ##### 3 class lolipop plot
  gmiss_conss_ptm_specific_protein <- eventReactive(input$submit, {
    #subset(gmiss_dssp_ptm, uniprotAc == input$proteinSelected)
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      subset(gmiss_dssp_ptm, uniprotAc == geneSpecificpmissOrgmiss$uniprotAC)
    }
    else if(input$homeSideBarTabSetPanel == 'Protein' && input$proteinSelected != ''){
      subset(gmiss_dssp_ptm, uniprotAc == input$proteinSelected)
    }
  })
  
  gmiss_conss_ptm_specific_protein_formatted <- eventReactive(input$submit, {
    gmiss_conss_ptm_specific_protein()[, c("unipos","conSS")]
  })
  
  gmiss_loli_label <- eventReactive(input$submit, { data.frame(rep("Population Variant",nrow(gmiss_conss_ptm_specific_protein_formatted()))) })
  gmiss_conss_ptm_specific_protein_comb <- eventReactive(input$submit, { cbind(gmiss_conss_ptm_specific_protein_formatted(), gmiss_loli_label()) })
  gmiss_SS3_loli_final <- eventReactive(input$submit, { setNames(gmiss_conss_ptm_specific_protein_comb(), c("unipos", "conSS","type")) })
  
  pmiss_conss_ptm_specific_protein <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      subset(pmiss_dssp_ptm, uniprotAc == geneSpecificpmissOrgmiss$uniprotAC)
    }
    else if(input$homeSideBarTabSetPanel == 'Protein' && input$proteinSelected != ''){
      subset(pmiss_dssp_ptm, uniprotAc == input$proteinSelected)
    }
  })
  protein_len_loli <- eventReactive(input$submit, {
    GPTLPC_list[GPTLPC_list$geneName == input$geneSelected,4]
  })
  
  pmiss_conss_ptm_specific_protein_formatted <- eventReactive(input$submit, {
    pmiss_conss_ptm_specific_protein()[, c("unipos","conSS")]
  })
  
  pmiss_loli_label <- eventReactive(input$submit, { data.frame(rep("Pathogenic Variant",nrow(pmiss_conss_ptm_specific_protein_formatted()))) })
  pmiss_conss_ptm_specific_protein_comb <- eventReactive(input$submit, { cbind(pmiss_conss_ptm_specific_protein_formatted(), pmiss_loli_label()) })
  pmiss_SS3_loli_final <- eventReactive(input$submit, { setNames(pmiss_conss_ptm_specific_protein_comb(), c("unipos", "conSS","type")) })
  
  
  conss_ptm_specific_protein_combined <- eventReactive(input$submit,{
    rbind(gmiss_SS3_loli_final(),pmiss_SS3_loli_final())
  })
  
  output$SS3loli = renderPlotly(
    {
      if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
        #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
        ss3_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
        ss3_loli_gene_info <- read_delim(ss3_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
        gene_length <- nrow(ss3_loli_gene_info)
      }
      
      #this_gene_path_pop <- ss3_loli_gene_info[((ss3_loli_gene_info$`Pathogenic Mutation`!="none") | (ss3_loli_gene_info$`Population Mutation`!="none")), ]
      this_gene_path_pop <- ss3_loli_gene_info
      this_gene_aa <- as.data.frame(this_gene_path_pop$`Amino Acid Index`)
      this_gene_ss3 <- as.data.frame(this_gene_path_pop$`3-class DSSP Secondary Structure Properties`)
      this_gene_mut <- cbind(this_gene_aa, this_gene_ss3)
      colnames(this_gene_mut) <- c("aaindex", "ss3type")
      
      mut_count <- nrow(this_gene_mut)
      
      for(i in 1:mut_count){
        if((this_gene_path_pop$`Population Mutation`[i] != "none") && (this_gene_path_pop$`Pathogenic Mutation`[i] =="none")){
          this_gene_mut$var_type[i] <- "Population variant"
        }else if((this_gene_path_pop$`Population Mutation`[i] == "none") && (this_gene_path_pop$`Pathogenic Mutation`[i] !="none")){
          this_gene_mut$var_type[i] <- "Pathogenic variant"
        }else if((this_gene_path_pop$`Population Mutation`[i] != "none") && (this_gene_path_pop$`Pathogenic Mutation`[i] !="none")){
          this_gene_mut$var_type[i] <- "Population and Pathogenic variant"
        }
        else
        {
          this_gene_mut$var_type[i] <- "No variant: no color"
        }
      }
      temp <- list()
      for(i in c(1:mut_count)){
        shape <- list()
        shape[['type']] <- "line"
        shape[['xref']] <- "x"
        shape[['yref']] <- "y"
        shape[['x0']] <- this_gene_mut$aaindex[i]
        shape[['y0']] <- -1
        shape[['x1']] <- this_gene_mut$aaindex[i]
        shape[['y1']] <- this_gene_mut$ss3type[i]
        shape[['line']] <- list(
          width = .6,
          color = "grey"
        )
        temp[[i]] = shape
        
      }
      #print(temp[1])
      pal <- c("blue", "red", "orange","white")
      pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
      
      p <- plot_ly(data = this_gene_mut, x = ~aaindex, y = ~ss3type, color = ~var_type, colors = pal, type = 'scatter',mode = 'markers',
                   text = ~paste("Amino acid: ", aaindex,", ",ss3_loli_gene_info$`Amino Acid`[aaindex]),
                   hoverinfo = 'text'
      ) %>%
        layout(
          xaxis = list(
            zeroline = FALSE,
            showline = TRUE,
            mirror = "ticks",
            title = "Amino acid index",
            tickangle = -65,
            range = c(-9,gene_length+10)
          ),
          yaxis = list(
            zeroline = FALSE,
            showline = TRUE,
            tickmode = "array",
            mirror = "ticks",
            categoryorder = "array",
            categoryarray = c("-","coils","beta-strand/sheet","Helices"),
            title = "Secondary structure types",
            tickvals = c("Helices","beta-strand/sheet","coils","-"),
            ticktext = c('helices',"\u03B2-strand/sheet",'coils','no structure')
          ),
          legend = list(
            x=.01,
            y=1.27,
            traceorder='normal',
            font = list(size = 11)
          ),
          annotations = list(
            x=.03,
            y=1.31,
            font = list(size = 12),
            xref='paper',
            yref='paper',
            text='Variant Type',
            showarrow=FALSE
          ),
          shapes = temp
        ) %>%
        config(displayModeBar = T)
      
      
    }
  )
  
  
  ## SS8 lollipop chart
  
  output$SS8loli <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ss8_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      ss8_loli_gene_info <- read_delim(ss8_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ss8 <- nrow(ss8_loli_gene_info)
    }
    
    #this_gene_path_pop_ss8 <- ss8_loli_gene_info[((ss8_loli_gene_info$`Pathogenic Mutation`!="none") | (ss8_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_ss8 <- ss8_loli_gene_info
    this_gene_aa_ss8 <- as.data.frame(this_gene_path_pop_ss8$`Amino Acid Index`)
    this_gene_ss8 <- as.data.frame(this_gene_path_pop_ss8$`8-class DSSP Secondary Structure Properties`)
    this_gene_mut_ss8 <- cbind(this_gene_aa_ss8, this_gene_ss8)
    colnames(this_gene_mut_ss8) <- c("aaindex", "ss8type")
    
    mut_count_ss8 <- nrow(this_gene_mut_ss8)
    
    for(i in 1:mut_count_ss8){
      if((this_gene_path_pop_ss8$`Population Mutation`[i] != "none") && (this_gene_path_pop_ss8$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ss8$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_ss8$`Population Mutation`[i] == "none") && (this_gene_path_pop_ss8$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ss8$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_ss8$`Population Mutation`[i] != "none") && (this_gene_path_pop_ss8$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ss8$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_ss8$var_type[i] <- "No variant: no color"
      }
    }
    
    temp <- list()
    for(i in c(1:mut_count_ss8)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_mut_ss8$aaindex[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_mut_ss8$aaindex[i]
      shape[['y1']] <- this_gene_mut_ss8$ss8type[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_mut_ss8, x = ~aaindex, y = ~ss8type, color = ~var_type, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",ss8_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_ss8+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          categoryorder = "array",
          categoryarray = c("-","310-helix","alpha-helix","pi-helix","bend","beta-sheet","beta-strand","loop","turn"),
          title = "Secondary structure types",
          tickvals = c("-","310-helix","alpha-helix","pi-helix","bend","beta-sheet","beta-strand","loop","turn"),
          ticktext = c('no structure',"3<sub>10</sub>-helix","\u03B1-helix","\u213C-helix","coil-bend","\u03B2-sheet","\u03B2-strand",'coil-loop','coil-turn')
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)
    
      
  })
  
  ## Chem Structure Type lollipop chart
  
  
  output$chemloli <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      chem_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      chem_loli_gene_info <- read_delim(chem_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_chem <- nrow(chem_loli_gene_info)
    }
    
    #this_gene_path_pop_chem <- chem_loli_gene_info[((chem_loli_gene_info$`Pathogenic Mutation`!="none") | (chem_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_chem <- chem_loli_gene_info
    this_gene_aa_chem <- as.data.frame(this_gene_path_pop_chem$`Amino Acid Index`)
    this_gene_chem <- as.data.frame(this_gene_path_pop_chem$`Amino Acid`)
    this_gene_mut_chem <- cbind(this_gene_aa_chem, this_gene_chem)
    colnames(this_gene_mut_chem) <- c("aaindex", "aatype")
    mut_count_chem <- nrow(this_gene_mut_chem)
    
    for(i in 1:mut_count_chem){
      if((this_gene_path_pop_chem$`Population Mutation`[i] != "none") && (this_gene_path_pop_chem$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_chem$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_chem$`Population Mutation`[i] == "none") && (this_gene_path_pop_chem$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_chem$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_chem$`Population Mutation`[i] != "none") && (this_gene_path_pop_chem$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_chem$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_chem$var_type[i] <- "No variant: no color"
      }
    }
    temp <- list()
    for(i in c(1:mut_count_chem)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_mut_chem$aaindex[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_mut_chem$aaindex[i]
      shape[['y1']] <- this_gene_mut_chem$aatype[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_mut_chem, x = ~aaindex, y = ~aatype, color = ~var_type, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",chem_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_chem+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          categoryorder = "array",
          categoryarray = c("A","I","L","M","V","F","W","Y","H","K","R","D","E","N","Q","S","T","C","G","P"),
          title = "Twenty amino acids and their physicochemical properties",
          tickvals = c("A","I","L","M","V","F","W","Y","H","K","R","D","E","N","Q","S","T","C","G","P"),
          ticktext = c( "Aliphatic (A)","Aliphatic (I)","Aliphatic (L)","Aliphatic (M)","Aliphatic (V)","Aromatic (F)","Aromatic (F)","Aromatic (Y)","Positively-charged (H)","Positively-charged (K)","Positively-charged (R)","Negatively-charged (D)","Negatively-charged (E)","Neutral (N)","Neutral(Q)","Neutral (S)","Neutral (T)","Special (C)","Special (G)","Special (P)")
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
          #xanchor = "center",
          #yanchor = "top"
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)  
    
  })
  
  #accessible surface area
  
  output$asaloli <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      asa_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      asa_loli_gene_info <- read_delim(asa_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_asa <- nrow(asa_loli_gene_info)
    }
    
    #this_gene_path_pop_asa <- asa_loli_gene_info[((asa_loli_gene_info$`Pathogenic Mutation`!="none") | (asa_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_asa <- asa_loli_gene_info
    this_gene_aa_asa <- as.data.frame(this_gene_path_pop_asa$`Amino Acid Index`)
    this_gene_asa <- as.data.frame(this_gene_path_pop_asa$`Accessible Surface Area`)
    this_gene_mut_asa <- cbind(this_gene_aa_asa, this_gene_asa)
    colnames(this_gene_mut_asa) <- c("aaindex", "asa")
    this_gene_mut_asa[this_gene_mut_asa == -1] <- -25
    mut_count_asa <- nrow(this_gene_mut_asa)
    
    for(i in 1:mut_count_asa){
      if((this_gene_path_pop_asa$`Population Mutation`[i] != "none") && (this_gene_path_pop_asa$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_asa$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_asa$`Population Mutation`[i] == "none") && (this_gene_path_pop_asa$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_asa$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_asa$`Population Mutation`[i] != "none") && (this_gene_path_pop_asa$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_asa$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_asa$var_type[i] <- "No variant: no color"
      }
    }
    
   
    
    temp <- list()
    for(i in c(1:mut_count_asa)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_mut_asa$aaindex[i]
      shape[['y0']] <- 0
      shape[['x1']] <- this_gene_mut_asa$aaindex[i]
      shape[['y1']] <- this_gene_mut_asa$asa[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_mut_asa, x = ~aaindex, y = ~asa, color = ~var_type, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",asa_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_asa+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "linear",
          mirror = "ticks",
          tick0 = -25,
          dtick = 25,
          title = "Accessible surface area (\u212B<sup>2</sup>)"
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
          #xanchor = "center",
          #yanchor = "top"
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)
    
  })
  
  ##### bond type lolipop plot
  gmiss_conss_ptm_specific_protein_bond <- eventReactive(input$submit, {
    #subset(gmiss_dssp_ptm, uniprotAc == input$proteinSelected)
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      subset(bond_type_gmiss, UniprotAC == as.character(geneSpecificpmissOrgmiss$uniprotAC))#geneSpecificpmissOrgmiss$uniprotAC
    }
    else if(input$homeSideBarTabSetPanel == 'Protein' && input$proteinSelected != ''){
      subset(bond_type_gmiss, UniprotAC == input$proteinSelected)
    }
  })
  
  gmiss_conss_ptm_specific_protein_formatted_bond <- eventReactive(input$submit, {
    gmiss_conss_ptm_specific_protein_bond()[, c("unipos","bond_type","count")]
  })
  
  gmiss_loli_label_bond <- eventReactive(input$submit, { data.frame(rep("Population Variant",nrow(gmiss_conss_ptm_specific_protein_formatted_bond()))) })
  gmiss_conss_ptm_specific_protein_comb_bond <- eventReactive(input$submit, { cbind(gmiss_conss_ptm_specific_protein_formatted_bond(), gmiss_loli_label_bond()) })
  gmiss_bond_loli_final <- eventReactive(input$submit, { setNames(gmiss_conss_ptm_specific_protein_comb_bond(), c("unipos", "bond_type","count","type")) })
  
  pmiss_conss_ptm_specific_protein_bond <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      subset(bond_type_pmiss, UniprotAC == as.character(geneSpecificpmissOrgmiss$uniprotAC))#geneSpecificpmissOrgmiss$uniprotAC
    }
    else if(input$homeSideBarTabSetPanel == 'Protein' && input$proteinSelected != ''){
      subset(bond_type_pmiss, UniprotAC == input$proteinSelected)
    }
  })
  
  pmiss_conss_ptm_specific_protein_formatted_bond <- eventReactive(input$submit, {
    pmiss_conss_ptm_specific_protein_bond()[, c("unipos","bond_type","count")]
  })
  
  pmiss_loli_label_bond <- eventReactive(input$submit, { data.frame(rep("Pathogenic Variant",nrow(pmiss_conss_ptm_specific_protein_formatted_bond()))) })
  pmiss_conss_ptm_specific_protein_comb_bond <- eventReactive(input$submit, { cbind(pmiss_conss_ptm_specific_protein_formatted_bond(), pmiss_loli_label_bond()) })
  pmiss_bond_loli_final <- eventReactive(input$submit, { setNames(pmiss_conss_ptm_specific_protein_comb_bond(), c("unipos", "bond_type","count","type")) })
  
  
  conss_ptm_specific_protein_combined_bond <- eventReactive(input$submit,{
    rbind(gmiss_bond_loli_final(),pmiss_bond_loli_final())
    #ranges$x <- NULL
    #ranges$y <- NULL
    
  })
  
  bond_names <- list(
    'H'="Hydrogen",
    'D'="Disulphide",
    'N'="Nonbonded",
    'S'="Salt-bridge"
  )
  
  bond_labeller <- function(variable,value){
    return(bond_names[value])
  }
  
  output$bondloli <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ib4_loli_gene_name = paste("gene_wise_info_2D/",input$geneSelected,".txt",sep='')
      ib4_loli_gene_info <- read_delim(ib4_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ib4 <- nrow(ib4_loli_gene_info)
    }
  
    
    for(i in 1:gene_length_ib4){
      if(ib4_loli_gene_info$`8-class DSSP Secondary Structure Properties`[i] != "-"){
        ib4_loli_gene_info$`Bond Detail`[i] <- ib4_loli_gene_info$`Bond Detail`[i]
      }else if(ib4_loli_gene_info$`8-class DSSP Secondary Structure Properties`[i] == "-"){
        ib4_loli_gene_info$`Bond Detail`[i] <- "nos"
      }
    }
    
    #this_gene_path_pop_ib4 <- ib4_loli_gene_info[((ib4_loli_gene_info$`Pathogenic Mutation`!="none") | (ib4_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_ib4 <- ib4_loli_gene_info
    this_gene_aa_ib4 <- as.data.frame(this_gene_path_pop_ib4$`Amino Acid Index`)
    this_gene_ib4 <- as.data.frame(this_gene_path_pop_ib4$`Bond Detail`)
    this_gene_mut_ib4 <- cbind(this_gene_aa_ib4, this_gene_ib4)
    colnames(this_gene_mut_ib4) <- c("aaindex", "bonddetail")
    mut_count_ib4 <- nrow(this_gene_mut_ib4)
    
    for(i in 1:mut_count_ib4){
      if((this_gene_path_pop_ib4$`Population Mutation`[i] != "none") && (this_gene_path_pop_ib4$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ib4$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_ib4$`Population Mutation`[i] == "none") && (this_gene_path_pop_ib4$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ib4$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_ib4$`Population Mutation`[i] != "none") && (this_gene_path_pop_ib4$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ib4$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_ib4$var_type[i] <- "No variant: no color"
      }
    }

    
    this_gene_ib4_ext <- data.frame(aaind = as.numeric(), bondtype = as.character(), vartype = as.character())
    for(i in 1:mut_count_ib4){
      this_aaind = this_gene_mut_ib4$aaindex[i]
      this_vartype = this_gene_mut_ib4$var_type[i]
      this_bonddetail = this_gene_mut_ib4$bonddetail[i]
      
      if("-" == this_bonddetail){
        this_gene_ib4_ext_entry = data.frame(this_aaind, "nob", this_vartype)
        colnames(this_gene_ib4_ext_entry) = c("aaind", "bondtype", "vartype")
        this_gene_ib4_ext <- rbind(this_gene_ib4_ext, this_gene_ib4_ext_entry)
      }
      else if("nos" == this_bonddetail){
        this_gene_ib4_ext_entry = data.frame(this_aaind, "nos", this_vartype)
        colnames(this_gene_ib4_ext_entry) = c("aaind", "bondtype", "vartype")
        this_gene_ib4_ext <- rbind(this_gene_ib4_ext, this_gene_ib4_ext_entry)
      }
      else if(grepl("N", as.character(this_bonddetail))){
        this_gene_ib4_ext_entry = data.frame(this_aaind, "nonbonded", this_vartype)
        colnames(this_gene_ib4_ext_entry) = c("aaind", "bondtype", "vartype")
        this_gene_ib4_ext <- rbind(this_gene_ib4_ext, this_gene_ib4_ext_entry)
      }
      else if(grepl("D", as.character(this_bonddetail))){
        this_gene_ib4_ext_entry = data.frame(this_aaind, "disulfide", this_vartype)
        colnames(this_gene_ib4_ext_entry) = c("aaind", "bondtype", "vartype")
        this_gene_ib4_ext <- rbind(this_gene_ib4_ext, this_gene_ib4_ext_entry)
      }
      else if(grepl("H", as.character(this_bonddetail))){
        this_gene_ib4_ext_entry = data.frame(this_aaind, "hydrogen", this_vartype)
        colnames(this_gene_ib4_ext_entry) = c("aaind", "bondtype", "vartype")
        this_gene_ib4_ext <- rbind(this_gene_ib4_ext, this_gene_ib4_ext_entry)
      }
      else if(grepl("S", as.character(this_bonddetail))){
        this_gene_ib4_ext_entry = data.frame(this_aaind, "saltbridge", this_vartype)
        colnames(this_gene_ib4_ext_entry) = c("aaind", "bondtype", "vartype")
        this_gene_ib4_ext <- rbind(this_gene_ib4_ext, this_gene_ib4_ext_entry)
      }
    }
    
    temp <- list()
    for(i in c(1:mut_count_ib4)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_ib4_ext$aaind[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_ib4_ext$aaind[i]
      shape[['y1']] <- this_gene_ib4_ext$bondtype[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_ib4_ext, x = ~aaind, y = ~bondtype, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaind,", ",ib4_loli_gene_info$`Amino Acid`[aaind]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_ib4+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          categoryorder = "array",
          categoryarray = c("nos", "nob", "nonbonded","disulfide","saltbridge","hydrogen"),
          title = "interaction types (if available)",
          tickvals = c( "nonbonded","disulfide","saltbridge","hydrogen", "nos", "nob"),
          ticktext = c( "nonbonded","disulfide","saltbridge","hydrogen","no structure","no bond")
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
          #xanchor = "center",
          #yanchor = "top"
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)
    
    
  })
  
  
  ##### dist type lolipop plot
  
  gmiss_conss_ptm_specific_protein_dist <- eventReactive(input$submit, {
    #subset(gmiss_dssp_ptm, uniprotAc == input$proteinSelected)
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      subset(dist_gmiss, UniprotAC == as.character(geneSpecificpmissOrgmiss$uniprotAC))#geneSpecificpmissOrgmiss$uniprotAC
    }
    else if(input$homeSideBarTabSetPanel == 'Protein' && input$proteinSelected != ''){
      subset(dist_gmiss, UniprotAC == input$proteinSelected)
    }
  })
  
  gmiss_conss_ptm_specific_protein_formatted_dist <- eventReactive(input$submit, {
    gmiss_conss_ptm_specific_protein_dist()[, c("unipos","dist_type","value")]
  })
  
  gmiss_loli_label_dist <- eventReactive(input$submit, { data.frame(rep("Population Variant",nrow(gmiss_conss_ptm_specific_protein_formatted_dist()))) })
  gmiss_conss_ptm_specific_protein_comb_dist <- eventReactive(input$submit, { cbind(gmiss_conss_ptm_specific_protein_formatted_dist(), gmiss_loli_label_dist()) })
  gmiss_dist_loli_final <- eventReactive(input$submit, { setNames(gmiss_conss_ptm_specific_protein_comb_dist(), c("unipos", "dist_type","value","type")) })
  
  pmiss_conss_ptm_specific_protein_dist <- eventReactive(input$submit, {
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      geneSpecificpmissOrgmiss <- subset(pmissOrgmissGene,Gene_Name == input$geneSelected)
      subset(dist_pmiss, UniprotAC == as.character(geneSpecificpmissOrgmiss$uniprotAC))#geneSpecificpmissOrgmiss$uniprotAC
    }
    else if(input$homeSideBarTabSetPanel == 'Protein' && input$proteinSelected != ''){
      subset(dist_pmiss, UniprotAC == input$proteinSelected)
    }
  })
  
  pmiss_conss_ptm_specific_protein_formatted_dist <- eventReactive(input$submit, {
    pmiss_conss_ptm_specific_protein_dist()[, c("unipos","dist_type","value")]
  })
  
  pmiss_loli_label_dist <- eventReactive(input$submit, { data.frame(rep("Pathogenic Variant",nrow(pmiss_conss_ptm_specific_protein_formatted_dist()))) })
  pmiss_conss_ptm_specific_protein_comb_dist <- eventReactive(input$submit, { cbind(pmiss_conss_ptm_specific_protein_formatted_dist(), pmiss_loli_label_dist()) })
  pmiss_dist_loli_final <- eventReactive(input$submit, { setNames(pmiss_conss_ptm_specific_protein_comb_dist(), c("unipos", "dist_type","value","type")) })
  
  ranges_dist_phos <- reactiveValues(x = NULL, y = NULL)
  
  conss_ptm_specific_protein_combined_dist <- eventReactive(input$submit,{
    rbind(gmiss_dist_loli_final(),pmiss_dist_loli_final())
    #ranges$x <- NULL
    #ranges$y <- NULL
    
  })
  
  dist_names <- list(
    'dist_phos'="Phosphorylation",
    'dist_sumo'="Sumoylation",
    'dist_acet'="Acetylation",
    'dist_meth'="Methylation",
    'dist_ubiq'="Ubiquitination",
    'dist_gcln'="GclNAc",
    'dist_galn'="GlaNAc"
  )
  
  dist_labeller <- function(variable,value){
    return(dist_names[value])
  }
  
  output$distloliphos <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ptm6_loli_gene_name = paste("gene_wise_info_2D/",input$geneSelected,".txt",sep='')
      ptm6_loli_gene_info <- read_delim(ptm6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    }
    
    #ptm6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/sf3d_final/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    
    #this_gene_path_pop_ptm6 <- ptm6_loli_gene_info[((ptm6_loli_gene_info$`Pathogenic Mutation`!="none") | (ptm6_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_ptm6 <- ptm6_loli_gene_info
    this_gene_aa_ptm6 <- as.data.frame(this_gene_path_pop_ptm6$`Amino Acid Index`)
    this_gene_acetyl <- as.data.frame(this_gene_path_pop_ptm6$`Acetylation Detail`)
    this_gene_methyl <- as.data.frame(this_gene_path_pop_ptm6$`Methylation Detail`)
    this_gene_gcln <- as.data.frame(this_gene_path_pop_ptm6$`O.GclNAc Detail`)
    this_gene_phos <- as.data.frame(this_gene_path_pop_ptm6$`Phosphorylation Detail`)
    this_gene_sumoy <- as.data.frame(this_gene_path_pop_ptm6$`Sumoylation Detail`)
    this_gene_ubiq <- as.data.frame(this_gene_path_pop_ptm6$`Ubiquitination Detail`)
    this_gene_mut_ptm6 <- cbind(this_gene_aa_ptm6, this_gene_acetyl, this_gene_methyl, this_gene_gcln, this_gene_phos, this_gene_sumoy, this_gene_ubiq)
    colnames(this_gene_mut_ptm6) <- c("aaindex", "acetyl", "methyl", "gcln", "phos", "sumoy", "ubiq")
    mut_count_ptm6 <- nrow(this_gene_mut_ptm6)
    
    for(i in 1:mut_count_ptm6){
      if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] == "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_ptm6$var_type[i] <- "No variant: no color"
      }
    }
    
    
    this_gene_ptm6_ext <- data.frame(aaind = as.numeric(), ptmtype = as.character(), ptmdist = as.numeric(), vartype = as.character())
    for(i in 1:mut_count_ptm6){
      this_aaind = this_gene_mut_ptm6$aaindex[i]
      this_vartype = this_gene_mut_ptm6$var_type[i]
      
      this_acetyldetail = as.character(this_gene_mut_ptm6$acetyl[i])
      this_methyldetail = as.character(this_gene_mut_ptm6$methyl[i])
      this_gclndetail = as.character(this_gene_mut_ptm6$gcln[i])
      this_phosdetail = as.character(this_gene_mut_ptm6$phos[i])
      this_sumoydetail = as.character(this_gene_mut_ptm6$sumoy[i])
      this_ubiqdetail = as.character(this_gene_mut_ptm6$ubiq[i])
      
      if((this_acetyldetail == "-") && (this_methyldetail == "-") && (this_gclndetail == "-") && (this_phosdetail == "-") && (this_sumoydetail == "-") && (this_ubiqdetail == "-")){
        this_gene_ptm6_ext_entry = data.frame(this_aaind, "No Annotation", -25, this_vartype)
        colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
        this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
      }else{
        if(this_phosdetail != "-"){
          this_phoslist = unlist(strsplit(this_phosdetail, split = "[;]"))
          for(j in 1:length(this_phoslist)){

            this_phosdist = as.numeric(this_phoslist[j])
              
            if(this_phosdist < 10.0){
              this_gene_ptm6_ext_entry = data.frame(this_aaind, "Phosphorylation", this_phosdist, this_vartype)
              colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
              this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
            }
          }
        }
      }
    }
    
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype != "No Annotation", ]
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Phosphorylation", ]
      
      mut_count_ptm6 = nrow(this_gene_ptm6_ext)
      
      if(nrow(this_gene_ptm6_ext) != 0){
      
        temp <- list()
        for(i in c(1:mut_count_ptm6)){
          shape <- list()
          shape[['type']] <- "line"
          shape[['xref']] <- "x"
          shape[['yref']] <- "y"
          shape[['x0']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y0']] <- -1
          shape[['x1']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y1']] <- this_gene_ptm6_ext$ptmdist[i]
          shape[['line']] <- list(
            width = .6,
            color = "grey"
          )
          temp[[i]] = shape
          
        }
        #print(temp[1])
        pal <- c("blue", "red", "orange","white")
        pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
        
        p <- plot_ly(data = this_gene_ptm6_ext, x = ~aaind, y = ~ptmdist, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                     text = ~paste("Amino acid: ", aaind,", ",ptm6_loli_gene_info$`Amino Acid`[aaind]),
                     hoverinfo = 'text'
        ) %>%
          layout(
            showlegend = T,
            xaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              mirror = "ticks",
              title = "Amino acid index",
              tickangle = -65,
              range = c(-9,gene_length_ptm6+10)
            ),
            yaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              tickmode = "linear",
              mirror = "ticks",
              #categoryorder = "array",
              #categoryarray = c("nos", "nob", "nonbonded","disulfide","saltbridge","hydrogen"),
              title = "Distance (< 10 ,\u212B) from variant to Phosphorylation",
              tick0 = -20,
              dtick = 2.5
              #tickvals = c( "nonbonded","disulfide","saltbridge","hydrogen", "nos", "nob"),
              #ticktext = c( "nonbonded","disulfide","saltbridge","hydrogen","no structure","no bond")
            ),
            legend = list(
              x=.01,
              y=1.27,
              traceorder='normal',
              font = list(size = 11)
              #xanchor = "center",
              #yanchor = "top"
            ),
            annotations = list(
              x=.03,
              y=1.31,
              font = list(size = 12),
              xref='paper',
              yref='paper',
              text='Variant Type',
              showarrow=FALSE
            ),
            shapes = temp
          ) %>%
          config(displayModeBar = T)
    }
  })
  
  output$distloliacet <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ptm6_loli_gene_name = paste("gene_wise_info_2D/",input$geneSelected,".txt",sep='')
      ptm6_loli_gene_info <- read_delim(ptm6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    }
    
    #ptm6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/sf3d_final/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    
    #this_gene_path_pop_ptm6 <- ptm6_loli_gene_info[((ptm6_loli_gene_info$`Pathogenic Mutation`!="none") | (ptm6_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_ptm6 <- ptm6_loli_gene_info
    this_gene_aa_ptm6 <- as.data.frame(this_gene_path_pop_ptm6$`Amino Acid Index`)
    this_gene_acetyl <- as.data.frame(this_gene_path_pop_ptm6$`Acetylation Detail`)
    this_gene_methyl <- as.data.frame(this_gene_path_pop_ptm6$`Methylation Detail`)
    this_gene_gcln <- as.data.frame(this_gene_path_pop_ptm6$`O.GclNAc Detail`)
    this_gene_phos <- as.data.frame(this_gene_path_pop_ptm6$`Phosphorylation Detail`)
    this_gene_sumoy <- as.data.frame(this_gene_path_pop_ptm6$`Sumoylation Detail`)
    this_gene_ubiq <- as.data.frame(this_gene_path_pop_ptm6$`Ubiquitination Detail`)
    this_gene_mut_ptm6 <- cbind(this_gene_aa_ptm6, this_gene_acetyl, this_gene_methyl, this_gene_gcln, this_gene_phos, this_gene_sumoy, this_gene_ubiq)
    colnames(this_gene_mut_ptm6) <- c("aaindex", "acetyl", "methyl", "gcln", "phos", "sumoy", "ubiq")
    mut_count_ptm6 <- nrow(this_gene_mut_ptm6)
    
    for(i in 1:mut_count_ptm6){
      if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] == "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_ptm6$var_type[i] <- "No variant: no color"
      }
    }
    
    #this_gene_mut_ptm6 <- this_gene_mut_ptm6[this_gene_mut_ptm6$bonddetail != "-", ]
    
    this_gene_ptm6_ext <- data.frame(aaind = as.numeric(), ptmtype = as.character(), ptmdist = as.numeric(), vartype = as.character())
    for(i in 1:mut_count_ptm6){
      this_aaind = this_gene_mut_ptm6$aaindex[i]
      this_vartype = this_gene_mut_ptm6$var_type[i]
      
      this_acetyldetail = as.character(this_gene_mut_ptm6$acetyl[i])
      this_methyldetail = as.character(this_gene_mut_ptm6$methyl[i])
      this_gclndetail = as.character(this_gene_mut_ptm6$gcln[i])
      this_phosdetail = as.character(this_gene_mut_ptm6$phos[i])
      this_sumoydetail = as.character(this_gene_mut_ptm6$sumoy[i])
      this_ubiqdetail = as.character(this_gene_mut_ptm6$ubiq[i])
      
      if((this_acetyldetail == "-") && (this_methyldetail == "-") && (this_gclndetail == "-") && (this_phosdetail == "-") && (this_sumoydetail == "-") && (this_ubiqdetail == "-")){
        this_gene_ptm6_ext_entry = data.frame(this_aaind, "No Annotation", -25, this_vartype)
        colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
        this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
      }else{
        if(this_acetyldetail != "-"){
          this_acetyllist = unlist(strsplit(this_acetyldetail, split = "[;]"))
          for(j in 1:length(this_acetyllist)){
            this_acetyldist = as.numeric(this_acetyllist[j])
            if(this_acetyldist < 10.0){
              this_gene_ptm6_ext_entry = data.frame(this_aaind, "Acetylation", this_acetyldist, this_vartype)
              colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
              this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
            }
          }
        }
      }
    }
    
    this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype != "No Annotation", ]
    this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Acetylation", ]
    
    mut_count_ptm6 = nrow(this_gene_ptm6_ext)
    
    
    if(nrow(this_gene_ptm6_ext) != 0){
      
      temp <- list()
      for(i in c(1:mut_count_ptm6)){
        shape <- list()
        shape[['type']] <- "line"
        shape[['xref']] <- "x"
        shape[['yref']] <- "y"
        shape[['x0']] <- this_gene_ptm6_ext$aaind[i]
        shape[['y0']] <- -1
        shape[['x1']] <- this_gene_ptm6_ext$aaind[i]
        shape[['y1']] <- this_gene_ptm6_ext$ptmdist[i]
        shape[['line']] <- list(
          width = .6,
          color = "grey"
        )
        temp[[i]] = shape
        
      }
      #print(temp[1])
      pal <- c("blue", "red", "orange","white")
      pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
      
      p <- plot_ly(data = this_gene_ptm6_ext, x = ~aaind, y = ~ptmdist, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                   text = ~paste("Amino acid: ", aaind,", ",ptm6_loli_gene_info$`Amino Acid`[aaind]),
                   hoverinfo = 'text'
      ) %>%
        layout(
          showlegend = T,
          xaxis = list(
            zeroline = FALSE,
            showline = TRUE,
            mirror = "ticks",
            title = "Amino acid index",
            tickangle = -65,
            range = c(-9,gene_length_ptm6+10)
            
          ),
          yaxis = list(
            zeroline = FALSE,
            showline = TRUE,
            tickmode = "linear",
            mirror = "ticks",
            #categoryorder = "array",
            #categoryarray = c("nos", "nob", "nonbonded","disulfide","saltbridge","hydrogen"),
            title = "Distance (< 10 ,\u212B) from variant to Acetylation",
            tick0 = -20,
            dtick = 2.5
            #tickvals = c( "nonbonded","disulfide","saltbridge","hydrogen", "nos", "nob"),
            #ticktext = c( "nonbonded","disulfide","saltbridge","hydrogen","no structure","no bond")
          ),
          legend = list(
            x=.01,
            y=1.27,
            traceorder='normal',
            font = list(size = 11)
            #xanchor = "center",
            #yanchor = "top"
          ),
          annotations = list(
            x=.03,
            y=1.31,
            font = list(size = 12),
            xref='paper',
            yref='paper',
            text='Variant Type',
            showarrow=FALSE
          ),
          shapes = temp
        ) %>%
        config(displayModeBar = T)
        
      }
  })
  
  output$distlolimet <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ptm6_loli_gene_name = paste("gene_wise_info_2D/",input$geneSelected,".txt",sep='')
      ptm6_loli_gene_info <- read_delim(ptm6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    }
    
    #ptm6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/sf3d_final/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    
    #this_gene_path_pop_ptm6 <- ptm6_loli_gene_info[((ptm6_loli_gene_info$`Pathogenic Mutation`!="none") | (ptm6_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_ptm6 <- ptm6_loli_gene_info
    this_gene_aa_ptm6 <- as.data.frame(this_gene_path_pop_ptm6$`Amino Acid Index`)
    this_gene_acetyl <- as.data.frame(this_gene_path_pop_ptm6$`Acetylation Detail`)
    this_gene_methyl <- as.data.frame(this_gene_path_pop_ptm6$`Methylation Detail`)
    this_gene_gcln <- as.data.frame(this_gene_path_pop_ptm6$`O.GclNAc Detail`)
    this_gene_phos <- as.data.frame(this_gene_path_pop_ptm6$`Phosphorylation Detail`)
    this_gene_sumoy <- as.data.frame(this_gene_path_pop_ptm6$`Sumoylation Detail`)
    this_gene_ubiq <- as.data.frame(this_gene_path_pop_ptm6$`Ubiquitination Detail`)
    this_gene_mut_ptm6 <- cbind(this_gene_aa_ptm6, this_gene_acetyl, this_gene_methyl, this_gene_gcln, this_gene_phos, this_gene_sumoy, this_gene_ubiq)
    colnames(this_gene_mut_ptm6) <- c("aaindex", "acetyl", "methyl", "gcln", "phos", "sumoy", "ubiq")
    mut_count_ptm6 <- nrow(this_gene_mut_ptm6)
    
    for(i in 1:mut_count_ptm6){
      if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] == "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_ptm6$var_type[i] <- "No variant: no color"
      }
    }
    
    #this_gene_mut_ptm6 <- this_gene_mut_ptm6[this_gene_mut_ptm6$bonddetail != "-", ]
    
    this_gene_ptm6_ext <- data.frame(aaind = as.numeric(), ptmtype = as.character(), ptmdist = as.numeric(), vartype = as.character())
    for(i in 1:mut_count_ptm6){
      this_aaind = this_gene_mut_ptm6$aaindex[i]
      this_vartype = this_gene_mut_ptm6$var_type[i]
      
      this_acetyldetail = as.character(this_gene_mut_ptm6$acetyl[i])
      this_methyldetail = as.character(this_gene_mut_ptm6$methyl[i])
      this_gclndetail = as.character(this_gene_mut_ptm6$gcln[i])
      this_phosdetail = as.character(this_gene_mut_ptm6$phos[i])
      this_sumoydetail = as.character(this_gene_mut_ptm6$sumoy[i])
      this_ubiqdetail = as.character(this_gene_mut_ptm6$ubiq[i])
      
      if((this_acetyldetail == "-") && (this_methyldetail == "-") && (this_gclndetail == "-") && (this_phosdetail == "-") && (this_sumoydetail == "-") && (this_ubiqdetail == "-")){
        this_gene_ptm6_ext_entry = data.frame(this_aaind, "No Annotation", -25, this_vartype)
        colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
        this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
      }else{
        if(this_methyldetail != "-"){
          this_methyllist = unlist(strsplit(this_methyldetail, split = "[;]"))
          for(j in 1:length(this_methyllist)){
            #this_methylall = unlist(strsplit(this_methyllist[j], split = "[/]"))
            #this_methyldist = as.numeric(this_methylall[2])
            this_methyldist = as.numeric(this_methyllist[j])
            
            if(this_methyldist < 10.0){
              this_gene_ptm6_ext_entry = data.frame(this_aaind, "Methylation", this_methyldist, this_vartype)
              colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
              this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
            }
          }
        }
      }
    }
    
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype != "No Annotation", ]
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Methylation", ]
      
      
      mut_count_ptm6 = nrow(this_gene_ptm6_ext)
      if(nrow(this_gene_ptm6_ext) != 0){
        
        temp <- list()
        for(i in c(1:mut_count_ptm6)){
          shape <- list()
          shape[['type']] <- "line"
          shape[['xref']] <- "x"
          shape[['yref']] <- "y"
          shape[['x0']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y0']] <- -1
          shape[['x1']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y1']] <- this_gene_ptm6_ext$ptmdist[i]
          shape[['line']] <- list(
            width = .6,
            color = "grey"
          )
          temp[[i]] = shape
          
        }
        #print(temp[1])
        pal <- c("blue", "red", "orange","white")
        pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
        
        p <- plot_ly(data = this_gene_ptm6_ext, x = ~aaind, y = ~ptmdist, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                     text = ~paste("Amino acid: ", aaind,", ",ptm6_loli_gene_info$`Amino Acid`[aaind]),
                     hoverinfo = 'text'
        ) %>%
          layout(
            showlegend = T,
            xaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              mirror = "ticks",
              title = "Amino acid index",
              tickangle = -65,
              range = c(-9,gene_length_ptm6+10)
            ),
            yaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              tickmode = "linear",
              mirror = "ticks",
              #categoryorder = "array",
              #categoryarray = c("nos", "nob", "nonbonded","disulfide","saltbridge","hydrogen"),
              title = "Distance (< 10 ,\u212B) from variant to Methylation",
              tick0 = -20,
              dtick = 2.5
              #tickvals = c( "nonbonded","disulfide","saltbridge","hydrogen", "nos", "nob"),
              #ticktext = c( "nonbonded","disulfide","saltbridge","hydrogen","no structure","no bond")
            ),
            legend = list(
              x=.01,
              y=1.27,
              traceorder='normal',
              font = list(size = 11)
              #xanchor = "center",
              #yanchor = "top"
            ),
            annotations = list(
              x=.03,
              y=1.31,
              font = list(size = 12),
              xref='paper',
              yref='paper',
              text='Variant Type',
              showarrow=FALSE
            ),
            shapes = temp
          ) %>%
          config(displayModeBar = T)
      }
  })
  
  output$distloliogcl <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ptm6_loli_gene_name = paste("gene_wise_info_2D/",input$geneSelected,".txt",sep='')
      ptm6_loli_gene_info <- read_delim(ptm6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    }
    
    #ptm6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/sf3d_final/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    
    #this_gene_path_pop_ptm6 <- ptm6_loli_gene_info[((ptm6_loli_gene_info$`Pathogenic Mutation`!="none") | (ptm6_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_ptm6 <- ptm6_loli_gene_info
    this_gene_aa_ptm6 <- as.data.frame(this_gene_path_pop_ptm6$`Amino Acid Index`)
    this_gene_acetyl <- as.data.frame(this_gene_path_pop_ptm6$`Acetylation Detail`)
    this_gene_methyl <- as.data.frame(this_gene_path_pop_ptm6$`Methylation Detail`)
    this_gene_gcln <- as.data.frame(this_gene_path_pop_ptm6$`O.GclNAc Detail`)
    this_gene_phos <- as.data.frame(this_gene_path_pop_ptm6$`Phosphorylation Detail`)
    this_gene_sumoy <- as.data.frame(this_gene_path_pop_ptm6$`Sumoylation Detail`)
    this_gene_ubiq <- as.data.frame(this_gene_path_pop_ptm6$`Ubiquitination Detail`)
    this_gene_mut_ptm6 <- cbind(this_gene_aa_ptm6, this_gene_acetyl, this_gene_methyl, this_gene_gcln, this_gene_phos, this_gene_sumoy, this_gene_ubiq)
    colnames(this_gene_mut_ptm6) <- c("aaindex", "acetyl", "methyl", "gcln", "phos", "sumoy", "ubiq")
    mut_count_ptm6 <- nrow(this_gene_mut_ptm6)
    
    for(i in 1:mut_count_ptm6){
      if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] == "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_ptm6$var_type[i] <- "No variant: no color"
      }
    }
    
    #this_gene_mut_ptm6 <- this_gene_mut_ptm6[this_gene_mut_ptm6$bonddetail != "-", ]
    
    this_gene_ptm6_ext <- data.frame(aaind = as.numeric(), ptmtype = as.character(), ptmdist = as.numeric(), vartype = as.character())
    for(i in 1:mut_count_ptm6){
      this_aaind = this_gene_mut_ptm6$aaindex[i]
      this_vartype = this_gene_mut_ptm6$var_type[i]
      
      this_acetyldetail = as.character(this_gene_mut_ptm6$acetyl[i])
      this_methyldetail = as.character(this_gene_mut_ptm6$methyl[i])
      this_gclndetail = as.character(this_gene_mut_ptm6$gcln[i])
      this_phosdetail = as.character(this_gene_mut_ptm6$phos[i])
      this_sumoydetail = as.character(this_gene_mut_ptm6$sumoy[i])
      this_ubiqdetail = as.character(this_gene_mut_ptm6$ubiq[i])
      
      if((this_acetyldetail == "-") && (this_methyldetail == "-") && (this_gclndetail == "-") && (this_phosdetail == "-") && (this_sumoydetail == "-") && (this_ubiqdetail == "-")){
        this_gene_ptm6_ext_entry = data.frame(this_aaind, "No Annotation", -25, this_vartype)
        colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
        this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
      }else{
        
        if(this_gclndetail != "-"){
          this_gclnlist = unlist(strsplit(this_gclndetail, split = "[;]"))
          for(j in 1:length(this_gclnlist)){
            #this_gclnall = unlist(strsplit(this_gclnlist[j], split = "[/]"))
            #this_gclndist = as.numeric(this_gclnall[2])
            this_gclndist = as.numeric(this_gclnlist[j])
            
            if(this_gclndist < 10.0){
              this_gene_ptm6_ext_entry = data.frame(this_aaind, "O.GclNAc", this_gclndist, this_vartype)
              colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
              this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
            }
          }
        }
      }
    }
    
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype != "No Annotation", ]
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "O.GclNAc", ]
      mut_count_ptm6 = nrow(this_gene_ptm6_ext)
      if(nrow(this_gene_ptm6_ext) != 0){
        
        temp <- list()
        for(i in c(1:mut_count_ptm6)){
          shape <- list()
          shape[['type']] <- "line"
          shape[['xref']] <- "x"
          shape[['yref']] <- "y"
          shape[['x0']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y0']] <- -1
          shape[['x1']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y1']] <- this_gene_ptm6_ext$ptmdist[i]
          shape[['line']] <- list(
            width = .6,
            color = "grey"
          )
          temp[[i]] = shape
          
        }
        #print(temp[1])
        pal <- c("blue", "red", "orange","white")
        pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
        
        p <- plot_ly(data = this_gene_ptm6_ext, x = ~aaind, y = ~ptmdist, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                     text = ~paste("Amino acid: ", aaind,", ",ptm6_loli_gene_info$`Amino Acid`[aaind]),
                     hoverinfo = 'text'
        ) %>%
          layout(
            showlegend = T,
            xaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              mirror = "ticks",
              title = "Amino acid index",
              tickangle = -65,
              range = c(-9,gene_length_ptm6+10)
              
            ),
            yaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              tickmode = "linear",
              mirror = "ticks",
              #categoryorder = "array",
              #categoryarray = c("nos", "nob", "nonbonded","disulfide","saltbridge","hydrogen"),
              title = "Distance (< 10 ,\u212B) from variant to O.GclNAc",
              tick0 = -20,
              dtick = 2.5
              #tickvals = c( "nonbonded","disulfide","saltbridge","hydrogen", "nos", "nob"),
              #ticktext = c( "nonbonded","disulfide","saltbridge","hydrogen","no structure","no bond")
            ),
            legend = list(
              x=.01,
              y=1.27,
              traceorder='normal',
              font = list(size = 11)
              #xanchor = "center",
              #yanchor = "top"
            ),
            annotations = list(
              x=.03,
              y=1.31,
              font = list(size = 12),
              xref='paper',
              yref='paper',
              text='Variant Type',
              showarrow=FALSE
            ),
            shapes = temp
          ) %>%
          config(displayModeBar = T)
      }
  })
  
  output$distlolisumo <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ptm6_loli_gene_name = paste("gene_wise_info_2D/",input$geneSelected,".txt",sep='')
      ptm6_loli_gene_info <- read_delim(ptm6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    }
    
    #ptm6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/sf3d_final/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    
    #this_gene_path_pop_ptm6 <- ptm6_loli_gene_info[((ptm6_loli_gene_info$`Pathogenic Mutation`!="none") | (ptm6_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_ptm6 <- ptm6_loli_gene_info
    this_gene_aa_ptm6 <- as.data.frame(this_gene_path_pop_ptm6$`Amino Acid Index`)
    this_gene_acetyl <- as.data.frame(this_gene_path_pop_ptm6$`Acetylation Detail`)
    this_gene_methyl <- as.data.frame(this_gene_path_pop_ptm6$`Methylation Detail`)
    this_gene_gcln <- as.data.frame(this_gene_path_pop_ptm6$`O.GclNAc Detail`)
    this_gene_phos <- as.data.frame(this_gene_path_pop_ptm6$`Phosphorylation Detail`)
    this_gene_sumoy <- as.data.frame(this_gene_path_pop_ptm6$`Sumoylation Detail`)
    this_gene_ubiq <- as.data.frame(this_gene_path_pop_ptm6$`Ubiquitination Detail`)
    this_gene_mut_ptm6 <- cbind(this_gene_aa_ptm6, this_gene_acetyl, this_gene_methyl, this_gene_gcln, this_gene_phos, this_gene_sumoy, this_gene_ubiq)
    colnames(this_gene_mut_ptm6) <- c("aaindex", "acetyl", "methyl", "gcln", "phos", "sumoy", "ubiq")
    mut_count_ptm6 <- nrow(this_gene_mut_ptm6)
    
    for(i in 1:mut_count_ptm6){
      if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] == "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_ptm6$var_type[i] <- "No variant: no color"
      }
    }
    
    #this_gene_mut_ptm6 <- this_gene_mut_ptm6[this_gene_mut_ptm6$bonddetail != "-", ]
    
    this_gene_ptm6_ext <- data.frame(aaind = as.numeric(), ptmtype = as.character(), ptmdist = as.numeric(), vartype = as.character())
    for(i in 1:mut_count_ptm6){
      this_aaind = this_gene_mut_ptm6$aaindex[i]
      this_vartype = this_gene_mut_ptm6$var_type[i]
      
      this_acetyldetail = as.character(this_gene_mut_ptm6$acetyl[i])
      this_methyldetail = as.character(this_gene_mut_ptm6$methyl[i])
      this_gclndetail = as.character(this_gene_mut_ptm6$gcln[i])
      this_phosdetail = as.character(this_gene_mut_ptm6$phos[i])
      this_sumoydetail = as.character(this_gene_mut_ptm6$sumoy[i])
      this_ubiqdetail = as.character(this_gene_mut_ptm6$ubiq[i])
      
      if((this_acetyldetail == "-") && (this_methyldetail == "-") && (this_gclndetail == "-") && (this_phosdetail == "-") && (this_sumoydetail == "-") && (this_ubiqdetail == "-")){
        this_gene_ptm6_ext_entry = data.frame(this_aaind, "No Annotation", -25, this_vartype)
        colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
        this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
      }else{
        
        if(this_sumoydetail != "-"){
          this_sumoylist = unlist(strsplit(this_sumoydetail, split = "[;]"))
          for(j in 1:length(this_sumoylist)){
            #this_sumoyall = unlist(strsplit(this_sumoylist[j], split = "[/]"))
            #this_sumoydist = as.numeric(this_sumoyall[2])
            this_sumoydist = as.numeric(this_sumoylist[j])
            
            if(this_sumoydist < 10.0){
              this_gene_ptm6_ext_entry = data.frame(this_aaind, "Sumoylation", this_sumoydist, this_vartype)
              colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
              this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
            }
          }
        }
      }
    }
    
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype != "No Annotation", ]
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Sumoylation", ]
      
      mut_count_ptm6 = nrow(this_gene_ptm6_ext)
      
      if(nrow(this_gene_ptm6_ext) != 0){
        
        temp <- list()
        for(i in c(1:mut_count_ptm6)){
          shape <- list()
          shape[['type']] <- "line"
          shape[['xref']] <- "x"
          shape[['yref']] <- "y"
          shape[['x0']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y0']] <- -1
          shape[['x1']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y1']] <- this_gene_ptm6_ext$ptmdist[i]
          shape[['line']] <- list(
            width = .6,
            color = "grey"
          )
          temp[[i]] = shape
          
        }
        #print(temp[1])
        pal <- c("blue", "red", "orange","white")
        pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
        
        p <- plot_ly(data = this_gene_ptm6_ext, x = ~aaind, y = ~ptmdist, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                     text = ~paste("Amino acid: ", aaind,", ",ptm6_loli_gene_info$`Amino Acid`[aaind]),
                     hoverinfo = 'text'
        ) %>%
          layout(
            showlegend = T,
            xaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              mirror = "ticks",
              title = "Amino acid index",
              tickangle = -65,
              range = c(-9,gene_length_ptm6+10)
              
            ),
            yaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              tickmode = "linear",
              mirror = "ticks",
              #categoryorder = "array",
              #categoryarray = c("nos", "nob", "nonbonded","disulfide","saltbridge","hydrogen"),
              title = "Distance (< 10 ,\u212B) from variant to Sumoylation",
              tick0 = -20,
              dtick = 2.5
              #tickvals = c( "nonbonded","disulfide","saltbridge","hydrogen", "nos", "nob"),
              #ticktext = c( "nonbonded","disulfide","saltbridge","hydrogen","no structure","no bond")
            ),
            legend = list(
              x=.01,
              y=1.27,
              traceorder='normal',
              font = list(size = 11)
              #xanchor = "center",
              #yanchor = "top"
            ),
            annotations = list(
              x=.03,
              y=1.31,
              font = list(size = 12),
              xref='paper',
              yref='paper',
              text='Variant Type',
              showarrow=FALSE
            ),
            shapes = temp
          ) %>%
          config(displayModeBar = T)
      }
  })
  
  output$distloliubiq <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      ptm6_loli_gene_name = paste("gene_wise_info_2D/",input$geneSelected,".txt",sep='')
      ptm6_loli_gene_info <- read_delim(ptm6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    }
    
    #ptm6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/sf3d_final/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_ptm6 <- nrow(ptm6_loli_gene_info)
    
    #this_gene_path_pop_ptm6 <- ptm6_loli_gene_info[((ptm6_loli_gene_info$`Pathogenic Mutation`!="none") | (ptm6_loli_gene_info$`Population Mutation`!="none")), ]
    this_gene_path_pop_ptm6 <- ptm6_loli_gene_info
    this_gene_aa_ptm6 <- as.data.frame(this_gene_path_pop_ptm6$`Amino Acid Index`)
    this_gene_acetyl <- as.data.frame(this_gene_path_pop_ptm6$`Acetylation Detail`)
    this_gene_methyl <- as.data.frame(this_gene_path_pop_ptm6$`Methylation Detail`)
    this_gene_gcln <- as.data.frame(this_gene_path_pop_ptm6$`O.GclNAc Detail`)
    this_gene_phos <- as.data.frame(this_gene_path_pop_ptm6$`Phosphorylation Detail`)
    this_gene_sumoy <- as.data.frame(this_gene_path_pop_ptm6$`Sumoylation Detail`)
    this_gene_ubiq <- as.data.frame(this_gene_path_pop_ptm6$`Ubiquitination Detail`)
    this_gene_mut_ptm6 <- cbind(this_gene_aa_ptm6, this_gene_acetyl, this_gene_methyl, this_gene_gcln, this_gene_phos, this_gene_sumoy, this_gene_ubiq)
    colnames(this_gene_mut_ptm6) <- c("aaindex", "acetyl", "methyl", "gcln", "phos", "sumoy", "ubiq")
    mut_count_ptm6 <- nrow(this_gene_mut_ptm6)
    
    for(i in 1:mut_count_ptm6){
      if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] =="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] == "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Pathogenic variant"
      }else if((this_gene_path_pop_ptm6$`Population Mutation`[i] != "none") && (this_gene_path_pop_ptm6$`Pathogenic Mutation`[i] !="none")){
        this_gene_mut_ptm6$var_type[i] <- "Population and Pathogenic variant"
      }
      else{
        this_gene_mut_ptm6$var_type[i] <- "No variant: no color"
      }
    }
    
    #this_gene_mut_ptm6 <- this_gene_mut_ptm6[this_gene_mut_ptm6$bonddetail != "-", ]
    
    this_gene_ptm6_ext <- data.frame(aaind = as.numeric(), ptmtype = as.character(), ptmdist = as.numeric(), vartype = as.character())
    for(i in 1:mut_count_ptm6){
      this_aaind = this_gene_mut_ptm6$aaindex[i]
      this_vartype = this_gene_mut_ptm6$var_type[i]
      
      this_acetyldetail = as.character(this_gene_mut_ptm6$acetyl[i])
      this_methyldetail = as.character(this_gene_mut_ptm6$methyl[i])
      this_gclndetail = as.character(this_gene_mut_ptm6$gcln[i])
      this_phosdetail = as.character(this_gene_mut_ptm6$phos[i])
      this_sumoydetail = as.character(this_gene_mut_ptm6$sumoy[i])
      this_ubiqdetail = as.character(this_gene_mut_ptm6$ubiq[i])
      
      if((this_acetyldetail == "-") && (this_methyldetail == "-") && (this_gclndetail == "-") && (this_phosdetail == "-") && (this_sumoydetail == "-") && (this_ubiqdetail == "-")){
        this_gene_ptm6_ext_entry = data.frame(this_aaind, "No Annotation", -25, this_vartype)
        colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
        this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
      }else{
        if(this_ubiqdetail != "-"){
          this_ubiqlist = unlist(strsplit(this_ubiqdetail, split = "[;]"))
          for(j in 1:length(this_ubiqlist)){
            #this_ubiqall = unlist(strsplit(this_ubiqlist[j], split = "[/]"))
            #this_ubiqdist = as.numeric(this_ubiqall[2])
            this_ubiqdist = as.numeric(this_ubiqlist[j])
            
            if(this_ubiqdist < 10.0){
              this_gene_ptm6_ext_entry = data.frame(this_aaind, "Ubiquitination", this_ubiqdist, this_vartype)
              colnames(this_gene_ptm6_ext_entry) = c("aaind", "ptmtype", "ptmdist", "vartype")
              this_gene_ptm6_ext <- rbind(this_gene_ptm6_ext, this_gene_ptm6_ext_entry)
            }
          }
        }
      }
    }
    
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype != "No Annotation", ]
      this_gene_ptm6_ext = this_gene_ptm6_ext[this_gene_ptm6_ext$ptmtype == "Ubiquitination", ]
      
      mut_count_ptm6 = nrow(this_gene_ptm6_ext)
      
      if(nrow(this_gene_ptm6_ext) != 0){
        
        temp <- list()
        for(i in c(1:mut_count_ptm6)){
          shape <- list()
          shape[['type']] <- "line"
          shape[['xref']] <- "x"
          shape[['yref']] <- "y"
          shape[['x0']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y0']] <- -1
          shape[['x1']] <- this_gene_ptm6_ext$aaind[i]
          shape[['y1']] <- this_gene_ptm6_ext$ptmdist[i]
          shape[['line']] <- list(
            width = .6,
            color = "grey"
          )
          temp[[i]] = shape
          
        }
        #print(temp[1])
        pal <- c("blue", "red", "orange","white")
        pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
        
        p <- plot_ly(data = this_gene_ptm6_ext, x = ~aaind, y = ~ptmdist, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                     text = ~paste("Amino acid: ", aaind,", ",ptm6_loli_gene_info$`Amino Acid`[aaind]),
                     hoverinfo = 'text'
        ) %>%
          layout(
            showlegend = T,
            xaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              mirror = "ticks",
              title = "Amino acid index",
              tickangle = -65,
              range = c(-9,gene_length_ptm6+10)
              
            ),
            yaxis = list(
              zeroline = FALSE,
              showline = TRUE,
              tickmode = "linear",
              mirror = "ticks",
              #categoryorder = "array",
              #categoryarray = c("nos", "nob", "nonbonded","disulfide","saltbridge","hydrogen"),
              title = "Distance (< 10 ,\u212B) from variant to Ubiquitination",
              tick0 = -20,
              dtick = 2.5
              #tickvals = c( "nonbonded","disulfide","saltbridge","hydrogen", "nos", "nob"),
              #ticktext = c( "nonbonded","disulfide","saltbridge","hydrogen","no structure","no bond")
            ),
            legend = list(
              x=.01,
              y=1.27,
              traceorder='normal',
              font = list(size = 11)
              #xanchor = "center",
              #yanchor = "top"
            ),
            annotations = list(
              x=.03,
              y=1.31,
              font = list(size = 12),
              xref='paper',
              yref='paper',
              text='Variant Type',
              showarrow=FALSE
            ),
            shapes = temp
          ) %>%
          config(displayModeBar = T)
      }
  })
  
  output$funclolifs <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      up6_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      up6_loli_gene_info <- read_delim(up6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_up6 <- nrow(up6_loli_gene_info)
    }
    
    #up6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/miscastapp/miscast01/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_up6 <- nrow(up6_loli_gene_info)
    
    #up6_loli_gene_info <- up6_loli_gene_info[((up6_loli_gene_info$`Pathogenic Mutation`!="none") | (up6_loli_gene_info$`Population Mutation`!="none")), ]
    
    this_gene_aa_up6_all <- data.frame(aaindex = as.numeric(), subtype = as.character(), uptype = as.character(), vartype = as.character())
    
    #functional site
    available_funcsite <- length(unique(up6_loli_gene_info$`Functional Site`))
    if(available_funcsite > 1){
      up6_loli_funcsite_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Functional Site` != "-"))
      this_gene_aa_funcsite <- as.data.frame(up6_loli_funcsite_info$`Amino Acid Index`)
      this_gene_funcsite <- as.data.frame(up6_loli_funcsite_info$`Functional Site`)
      this_gene_aa_funcsite <- cbind(this_gene_aa_funcsite, this_gene_funcsite)
      colnames(this_gene_aa_funcsite) <- c("aaindex", "subtype")
      this_gene_aa_funcsite$uptype <- "Functional Site"
      this_gene_aa_funcsite_len <- nrow(this_gene_aa_funcsite)
      
      for(i in 1:this_gene_aa_funcsite_len){
        if((up6_loli_funcsite_info$`Population Mutation`[i] != "none") && (up6_loli_funcsite_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_funcsite$vartype[i] <- "Population variant"
        }else if((up6_loli_funcsite_info$`Population Mutation`[i] == "none") && (up6_loli_funcsite_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_funcsite$vartype[i] <- "Pathogenic variant"
        }else if((up6_loli_funcsite_info$`Population Mutation`[i] != "none") && (up6_loli_funcsite_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_funcsite$vartype[i] <- "Population and Pathogenic variant"
        }else{
          this_gene_aa_funcsite$vartype[i] <- "No variant: no color"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_funcsite)
      
    }
    
    droplevels(this_gene_aa_up6_all$subtype)
    this_gene_aa_up6_all = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Functional Site", ]
    mut_count_func <- nrow(this_gene_aa_up6_all)
    temp <- list()
    for(i in c(1:mut_count_func)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y1']] <- this_gene_aa_up6_all$subtype[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_aa_up6_all, x = ~aaindex, y = ~subtype, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",up6_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        showlegend = T,
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_up6+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          #categoryorder = "array",
          #categoryarray = c("active site"),
          title = "Functional Site"
          #tickvals = c("active site"),
          #ticktext = c("active site")
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
          #xanchor = "center",
          #yanchor = "top"
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)
    
    
    
  })
  
  output$funclolimolp <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      up6_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      up6_loli_gene_info <- read_delim(up6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_up6 <- nrow(up6_loli_gene_info)
    }
    
    #up6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/miscastapp/miscast01/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_up6 <- nrow(up6_loli_gene_info)
    
    #up6_loli_gene_info <- up6_loli_gene_info[((up6_loli_gene_info$`Pathogenic Mutation`!="none") | (up6_loli_gene_info$`Population Mutation`!="none")), ]
    
    this_gene_aa_up6_all <- data.frame(aaindex = as.numeric(), subtype = as.character(), uptype = as.character(), vartype = as.character())
    
    #molecular processing
    available_molproc <- length(unique(up6_loli_gene_info$`Molecular Processing`))
    if(available_molproc > 1){
      up6_loli_molproc_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Molecular Processing` != "-"))
      this_gene_aa_molproc <- as.data.frame(up6_loli_molproc_info$`Amino Acid Index`)
      this_gene_molproc <- as.data.frame(up6_loli_molproc_info$`Molecular Processing`)
      this_gene_aa_molproc <- cbind(this_gene_aa_molproc, this_gene_molproc)
      colnames(this_gene_aa_molproc) <- c("aaindex", "subtype")
      this_gene_aa_molproc$uptype <- "Molecular Processing Associated Region"
      this_gene_aa_molproc_len <- nrow(this_gene_aa_molproc)
      
      for(i in 1:this_gene_aa_molproc_len){
        if((up6_loli_molproc_info$`Population Mutation`[i] != "none") && (up6_loli_molproc_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_molproc$vartype[i] <- "Population variant"
        }else if((up6_loli_molproc_info$`Population Mutation`[i] == "none") && (up6_loli_molproc_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_molproc$vartype[i] <- "Pathogenic variant"
        }else if((up6_loli_molproc_info$`Population Mutation`[i] != "none") && (up6_loli_molproc_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_molproc$vartype[i] <- "Population and Pathogenic variant"
        }else{
          this_gene_aa_molproc$vartype[i] <- "No variant: no color"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_molproc)
      
    }
    
    droplevels(this_gene_aa_up6_all$subtype)
    this_gene_aa_up6_all = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Molecular Processing Associated Region", ]
    
    mut_count_func <- nrow(this_gene_aa_up6_all)
    temp <- list()
    for(i in c(1:mut_count_func)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y1']] <- this_gene_aa_up6_all$subtype[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_aa_up6_all, x = ~aaindex, y = ~subtype, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",up6_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        showlegend = T,
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_up6+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          #categoryorder = "array",
          #categoryarray = c("signal peptide"),
          title = "Molecular Processing Associated Region"
          #tickvals = c("signal peptide"),
          #ticktext = c("signal peptide")
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
          #xanchor = "center",
          #yanchor = "top"
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)  
    
    
    
  })
  
  output$funclolifbr <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      up6_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      up6_loli_gene_info <- read_delim(up6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_up6 <- nrow(up6_loli_gene_info)
    }
    
    #up6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/miscastapp/miscast01/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_up6 <- nrow(up6_loli_gene_info)
    
    #up6_loli_gene_info <- up6_loli_gene_info[((up6_loli_gene_info$`Pathogenic Mutation`!="none") | (up6_loli_gene_info$`Population Mutation`!="none")), ]
    
    this_gene_aa_up6_all <- data.frame(aaindex = as.numeric(), subtype = as.character(), uptype = as.character(), vartype = as.character())
    
    #functional/binding region
    available_bindreg <- length(unique(up6_loli_gene_info$`Functional Binding Region`))
    if(available_bindreg > 1){
      up6_loli_bindreg_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Functional Binding Region` != "-"))
      this_gene_aa_bindreg <- as.data.frame(up6_loli_bindreg_info$`Amino Acid Index`)
      this_gene_bindreg <- as.data.frame(up6_loli_bindreg_info$`Functional Binding Region`)
      this_gene_aa_bindreg <- cbind(this_gene_aa_bindreg, this_gene_bindreg)
      colnames(this_gene_aa_bindreg) <- c("aaindex", "subtype")
      this_gene_aa_bindreg$uptype <- "Functional/Binding Region"
      this_gene_aa_bindreg_len <- nrow(this_gene_aa_bindreg)
      
      for(i in 1:this_gene_aa_bindreg_len){
        if((up6_loli_bindreg_info$`Population Mutation`[i] != "none") && (up6_loli_bindreg_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_bindreg$vartype[i] <- "Population variant"
        }else if((up6_loli_bindreg_info$`Population Mutation`[i] == "none") && (up6_loli_bindreg_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_bindreg$vartype[i] <- "Pathogenic variant"
        }else if((up6_loli_bindreg_info$`Population Mutation`[i] != "none") && (up6_loli_bindreg_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_bindreg$vartype[i] <- "Population and Pathogenic variant"
        }else{
          this_gene_aa_bindreg$vartype[i] <- "No variant: no color"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_bindreg)
      
    }
    
    droplevels(this_gene_aa_up6_all$subtype)
    this_gene_aa_up6_all = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Functional/Binding Region", ]
    mut_count_func <- nrow(this_gene_aa_up6_all)
    temp <- list()
    for(i in c(1:mut_count_func)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y1']] <- this_gene_aa_up6_all$subtype[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_aa_up6_all, x = ~aaindex, y = ~subtype, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",up6_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        showlegend = T,
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_up6+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          #categoryorder = "array",
          #categoryarray = c("active site"),
          title = "Functional/Binding Region"
          #tickvals = c("active site"),
          #ticktext = c("active site")
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
          #xanchor = "center",
          #yanchor = "top"
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T) 
    
    
    
  })
  
  output$funclolismr <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      up6_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      up6_loli_gene_info <- read_delim(up6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_up6 <- nrow(up6_loli_gene_info)
    }
    
    #up6_loli_gene_info <- read_delim("/Users/sumaiya/Broad_Work/VCF_HAIL_VEP/App/miscastapp/miscast01/gene_wise_info/CDKL5.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
    #gene_length_up6 <- nrow(up6_loli_gene_info)
    
    #up6_loli_gene_info <- up6_loli_gene_info[((up6_loli_gene_info$`Pathogenic Mutation`!="none") | (up6_loli_gene_info$`Population Mutation`!="none")), ]
    
    this_gene_aa_up6_all <- data.frame(aaindex = as.numeric(), subtype = as.character(), uptype = as.character(), vartype = as.character())
    
    #sequence mmotif/region
    available_seqmotif <- length(unique(up6_loli_gene_info$`Sequence Motif Region`))
    if(available_seqmotif > 1){
      up6_loli_seqmotif_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Sequence Motif Region` != "-"))
      this_gene_aa_seqmotif <- as.data.frame(up6_loli_seqmotif_info$`Amino Acid Index`)
      this_gene_seqmotif <- as.data.frame(up6_loli_seqmotif_info$`Sequence Motif Region`)
      this_gene_aa_seqmotif <- cbind(this_gene_aa_seqmotif, this_gene_seqmotif)
      colnames(this_gene_aa_seqmotif) <- c("aaindex", "subtype")
      this_gene_aa_seqmotif$uptype <- "Sequence Motif/Region"
      this_gene_aa_seqmotif_len <- nrow(this_gene_aa_seqmotif)
      
      for(i in 1:this_gene_aa_seqmotif_len){
        if((up6_loli_seqmotif_info$`Population Mutation`[i] != "none") && (up6_loli_seqmotif_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_seqmotif$vartype[i] <- "Population variant"
        }else if((up6_loli_seqmotif_info$`Population Mutation`[i] == "none") && (up6_loli_seqmotif_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_seqmotif$vartype[i] <- "Pathogenic variant"
        }else if((up6_loli_seqmotif_info$`Population Mutation`[i] != "none") && (up6_loli_seqmotif_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_seqmotif$vartype[i] <- "Population and Pathogenic variant"
        }else{
          this_gene_aa_seqmotif$vartype[i] <- "No variant: no color"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_seqmotif)
      
    }
    
    droplevels(this_gene_aa_up6_all$subtype)
    this_gene_aa_up6_all = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Sequence Motif/Region", ]
    mut_count_func <- nrow(this_gene_aa_up6_all)
    temp <- list()
    for(i in c(1:mut_count_func)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y1']] <- this_gene_aa_up6_all$subtype[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_aa_up6_all, x = ~aaindex, y = ~subtype, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",up6_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        showlegend = T,
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_up6+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          #categoryorder = "array",
          #categoryarray = c("active site"),
          title = "Sequence Motif/Region"
          #tickvals = c("active site"),
          #ticktext = c("active site")
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
          #xanchor = "center",
          #yanchor = "top"
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)
    
    
    
  })
  
  output$funclolimr <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      up6_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      up6_loli_gene_info <- read_delim(up6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_up6 <- nrow(up6_loli_gene_info)
    }

    this_gene_aa_up6_all <- data.frame(aaindex = as.numeric(), subtype = as.character(), uptype = as.character(), vartype = as.character())
    
    #modified residues
    available_modres <- length(unique(up6_loli_gene_info$`Modified Residues`))
    if(available_modres > 1){
      up6_loli_modres_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Modified Residues` != "-"))
      this_gene_aa_modres <- as.data.frame(up6_loli_modres_info$`Amino Acid Index`)
      this_gene_modres <- as.data.frame(up6_loli_modres_info$`Modified Residues`)
      this_gene_aa_modres <- cbind(this_gene_aa_modres, this_gene_modres)
      colnames(this_gene_aa_modres) <- c("aaindex", "subtype")
      this_gene_aa_modres$uptype <- "Modified Residues"
      this_gene_aa_modres_len <- nrow(this_gene_aa_modres)
      
      for(i in 1:this_gene_aa_modres_len){
        if((up6_loli_modres_info$`Population Mutation`[i] != "none") && (up6_loli_modres_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_modres$vartype[i] <- "Population variant"
        }else if((up6_loli_modres_info$`Population Mutation`[i] == "none") && (up6_loli_modres_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_modres$vartype[i] <- "Pathogenic variant"
        }else if((up6_loli_modres_info$`Population Mutation`[i] != "none") && (up6_loli_modres_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_modres$vartype[i] <- "Population and Pathogenic variant"
        }else{
          this_gene_aa_modres$vartype[i] <- "No variant: no color"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_modres)
      
    }
    
    droplevels(this_gene_aa_up6_all$subtype)
    this_gene_aa_up6_all = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Modified Residues", ]
    mut_count_func <- nrow(this_gene_aa_up6_all)
    temp <- list()
    for(i in c(1:mut_count_func)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y1']] <- this_gene_aa_up6_all$subtype[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_aa_up6_all, x = ~aaindex, y = ~subtype, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",up6_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        showlegend = T,
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_up6+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          title = "Modified Residues"
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)
    
    
    
  })
  
  output$funclolimd <- renderPlotly({
    if(input$homeSideBarTabSetPanel == 'Select a Gene' && input$geneSelected != ''){
      #sep <- switch(input$filetype, "csv" = ",", "txt" = "\t")
      up6_loli_gene_name = paste("gene_wise_info/",input$geneSelected,".txt",sep='')
      up6_loli_gene_info <- read_delim(up6_loli_gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_length_up6 <- nrow(up6_loli_gene_info)
    }
    
    this_gene_aa_up6_all <- data.frame(aaindex = as.numeric(), subtype = as.character(), uptype = as.character(), vartype = as.character())
    
    #modular domain
    available_moddom <- length(unique(up6_loli_gene_info$`Modular Domain`))
    if(available_moddom > 1){
      up6_loli_moddom_info <- as.data.frame(subset(up6_loli_gene_info, up6_loli_gene_info$`Modular Domain` != "-"))
      this_gene_aa_moddom <- as.data.frame(up6_loli_moddom_info$`Amino Acid Index`)
      this_gene_moddom <- as.data.frame(up6_loli_moddom_info$`Modular Domain`)
      this_gene_aa_moddom <- cbind(this_gene_aa_moddom, this_gene_moddom)
      colnames(this_gene_aa_moddom) <- c("aaindex", "subtype")
      this_gene_aa_moddom$uptype <- "Modular Domain"
      this_gene_aa_moddom_len <- nrow(this_gene_aa_moddom)
      
      for(i in 1:this_gene_aa_moddom_len){
        if((up6_loli_moddom_info$`Population Mutation`[i] != "none") && (up6_loli_moddom_info$`Pathogenic Mutation`[i] =="none")){
          this_gene_aa_moddom$vartype[i] <- "Population variant"
        }else if((up6_loli_moddom_info$`Population Mutation`[i] == "none") && (up6_loli_moddom_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_moddom$vartype[i] <- "Pathogenic variant"
        }else if((up6_loli_moddom_info$`Population Mutation`[i] != "none") && (up6_loli_moddom_info$`Pathogenic Mutation`[i] !="none")){
          this_gene_aa_moddom$vartype[i] <- "Population and Pathogenic variant"
        }else{
          this_gene_aa_moddom$vartype[i] <- "No variant: no color"
        }
      }
      
      this_gene_aa_up6_all <- rbind(this_gene_aa_up6_all, this_gene_aa_moddom)
    }
    
    droplevels(this_gene_aa_up6_all$subtype)
    this_gene_aa_up6_all = this_gene_aa_up6_all[this_gene_aa_up6_all$uptype == "Modular Domain", ]
    mut_count_func <- nrow(this_gene_aa_up6_all)
    temp <- list()
    for(i in c(1:mut_count_func)){
      shape <- list()
      shape[['type']] <- "line"
      shape[['xref']] <- "x"
      shape[['yref']] <- "y"
      shape[['x0']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y0']] <- -1
      shape[['x1']] <- this_gene_aa_up6_all$aaindex[i]
      shape[['y1']] <- this_gene_aa_up6_all$subtype[i]
      shape[['line']] <- list(
        width = .6,
        color = "grey"
      )
      temp[[i]] = shape
      
    }
    #print(temp[1])
    pal <- c("blue", "red", "orange","white")
    pal <- setNames(pal, c("Population variant", "Pathogenic variant", "Population and Pathogenic variant","No variant: no color"))
    
    p <- plot_ly(data = this_gene_aa_up6_all, x = ~aaindex, y = ~subtype, color = ~vartype, colors = pal, type = 'scatter',mode = 'markers',
                 text = ~paste("Amino acid: ", aaindex,", ",up6_loli_gene_info$`Amino Acid`[aaindex]),
                 hoverinfo = 'text'
    ) %>%
      layout(
        showlegend = T,
        xaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          mirror = "ticks",
          title = "Amino acid index",
          tickangle = -65,
          range = c(-9,gene_length_up6+10)
          
        ),
        yaxis = list(
          zeroline = FALSE,
          showline = TRUE,
          tickmode = "array",
          mirror = "ticks",
          #categoryorder = "array",
          #categoryarray = c("active site"),
          title = "Modular Domain"
          #tickvals = c("active site"),
          #ticktext = c("active site")
        ),
        legend = list(
          x=.01,
          y=1.27,
          traceorder='normal',
          font = list(size = 11)
          #xanchor = "center",
          #yanchor = "top"
        ),
        annotations = list(
          x=.03,
          y=1.31,
          font = list(size = 12),
          xref='paper',
          yref='paper',
          text='Variant Type',
          showarrow=FALSE
        ),
        shapes = temp
      ) %>%
      config(displayModeBar = T)  
    
    
    
  })
  
  output$full_gene_wise_info = DT::renderDataTable({
    if(input$geneSelected != '' && input$tableFilterSelected=="All amino acids"){
      gene_name = paste("gene_wise_info_table/",input$geneSelected,".txt",sep='')
      gene_wise_inf<- read_delim(gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      gene_wise_info_filtered <- gene_wise_inf
      
      colnames(gene_wise_info_filtered)[14] <- "Acetylation Detail(*)"
      colnames(gene_wise_info_filtered)[15] <- "Methylation Detail(*)"
      colnames(gene_wise_info_filtered)[17] <- "O.GclNAc Detail(*)"
      colnames(gene_wise_info_filtered)[18] <- "Phosphorylation Detail(*)"
      colnames(gene_wise_info_filtered)[19] <- "Sumoylation Detail(*)"
      colnames(gene_wise_info_filtered)[20] <- "Ubiquitination Detail(*)"
      colnames(gene_wise_info_filtered)[7] <- "Bond Detail(*)"
      
      if(length(input$show_columns_of_Feature_Table) == 0){
        DT::datatable(gene_wise_inf,
                      class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                      options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = TRUE, #fixedHeader = TRUE,paging = FALSE,
                                     fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                     columnDefs = list(list(width = '20px', targets = '_all'))))
      }

      
      else{
        DT::datatable(gene_wise_info_filtered[,input$show_columns_of_Feature_Table,drop = FALSE],
                    class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                    options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = TRUE, #fixedHeader = TRUE,paging = FALSE,
                                   fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                   columnDefs = list(list(width = '20px', targets = '_all'))))
      }
    }
    
  }
  )
  
  output$gnomad_gene_wise_info = DT::renderDataTable({
    if(input$geneSelected != '' && input$tableFilterSelected=="Amino acids with population variants"){
      gene_name = paste("gene_wise_info_table/",input$geneSelected,".txt",sep='')
      gene_wise_inf<- read_delim(gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      #gene_wise_info<-gene_wise_inf[-grep('Disease', colnames(gene_wise_inf))]
      #gene_wise_info_filtered <- gene_wise_info[-grep('DBSNP ID', colnames(gene_wise_info))]
      gene_wise_info_filtered <- gene_wise_inf
      
      colnames(gene_wise_info_filtered)[14] <- "Acetylation Detail(*)"
      colnames(gene_wise_info_filtered)[15] <- "Methylation Detail(*)"
      colnames(gene_wise_info_filtered)[17] <- "O.GclNAc Detail(*)"
      colnames(gene_wise_info_filtered)[18] <- "Phosphorylation Detail(*)"
      colnames(gene_wise_info_filtered)[19] <- "Sumoylation Detail(*)"
      colnames(gene_wise_info_filtered)[20] <- "Ubiquitination Detail(*)"
      colnames(gene_wise_info_filtered)[7] <- "Bond Detail(*)"
      
      gnomad_info <- subset(gene_wise_info_filtered,gene_wise_info_filtered$'Population Mutation' != "none")
      if(length(input$show_columns_of_Feature_Table) == 0){
        DT::datatable(gnomad_info,
                      class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                      options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = TRUE, #fixedHeader = TRUE,paging = FALSE,
                                     fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                     columnDefs = list(list(width = '20px', targets = '_all'))))
      }
      else{
        DT::datatable(gnomad_info[,input$show_columns_of_Feature_Table,drop = FALSE],
                      class = 'cell-border stripe',rownames = FALSE, extensions = c("FixedColumns"),
                      options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = FALSE,# fixedHeader = TRUE,paging = FALSE,
                                     fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                     columnDefs = list(list(width = '20px', targets = '_all'))))
      }
    }
    
  })
  output$patient_gene_wise_info = DT::renderDataTable({
    if(input$geneSelected != '' && input$tableFilterSelected=="Amino acids with pathogenic variants"){
      gene_name = paste("gene_wise_info_table/",input$geneSelected,".txt",sep='')
      gene_wise_inf<- read_delim(gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      #gene_wise_info<-gene_wise_inf[-grep('Disease', colnames(gene_wise_inf))]
      #gene_wise_info_filtered <- gene_wise_info[-grep('DBSNP ID', colnames(gene_wise_info))]
      gene_wise_info_filtered <- gene_wise_inf
      
      colnames(gene_wise_info_filtered)[14] <- "Acetylation Detail(*)"
      colnames(gene_wise_info_filtered)[15] <- "Methylation Detail(*)"
      colnames(gene_wise_info_filtered)[17] <- "O.GclNAc Detail(*)"
      colnames(gene_wise_info_filtered)[18] <- "Phosphorylation Detail(*)"
      colnames(gene_wise_info_filtered)[19] <- "Sumoylation Detail(*)"
      colnames(gene_wise_info_filtered)[20] <- "Ubiquitination Detail(*)"
      colnames(gene_wise_info_filtered)[7] <- "Bond Detail(*)"
      
      patient_info <- subset(gene_wise_info_filtered,gene_wise_info_filtered$'Pathogenic Mutation' != "none")
      if(length(input$show_columns_of_Feature_Table) == 0){
        DT::datatable(patient_info,
                      class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                      options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = TRUE, #fixedHeader = TRUE,paging = FALSE,
                                     fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                     columnDefs = list(list(width = '20px', targets = '_all'))))
      }
      else{
        DT::datatable(patient_info[,input$show_columns_of_Feature_Table,drop = FALSE],
                      class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                      options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = FALSE,
                                     fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                     columnDefs = list(list(width = '20px', targets = '_all'))))
      }
      
    }
    
  })
  
  output$structure_all_gene_wise_info = DT::renderDataTable({
    if(input$geneSelected != '' && input$indextableFilterSelected=="All amino acids"){
      gene_name = paste("index_table_files/",input$geneSelected,".txt",sep='')
      gene_wise_inf<- read_delim(gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      
      if(length(input$show_columns_of_Index_Table) == 0){
        DT::datatable(gene_wise_inf,
                      class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                      options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = TRUE, #fixedHeader = TRUE,paging = FALSE,
                                     fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                     columnDefs = list(list(width = '20px', targets = '_all'))))
      }

      else{
        DT::datatable(gene_wise_inf[,input$show_columns_of_Index_Table,drop = FALSE],
                    class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                    options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = TRUE, #fixedHeader = TRUE,paging = FALSE,
                                   fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                   columnDefs = list(list(width = '20px', targets = '_all'))))
      }
    }
    
  })
  
  output$structure_gene_wise_info = DT::renderDataTable({
    if(input$geneSelected != '' && input$indextableFilterSelected=="Amino acids with structure"){
      gene_name = paste("index_table_files/",input$geneSelected,".txt",sep='')
      gene_wise_inf<- read_delim(gene_name, "\t", escape_double = FALSE, trim_ws = TRUE)
      structure_info <- subset(gene_wise_inf,gene_wise_inf$'Structure Index of the Amino Acid' != ".")
      
      if(length(input$show_columns_of_Index_Table) == 0){
        DT::datatable(gene_wise_inf,
                      class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                      options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = TRUE, #fixedHeader = TRUE,paging = FALSE,
                                     fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                     columnDefs = list(list(width = '20px', targets = '_all'))))
      }

      else{
        DT::datatable(structure_info[,input$show_columns_of_Index_Table,drop = FALSE],
                    class = 'cell-border stripe',rownames = FALSE,extensions = c("FixedColumns"),
                    options = list(pageLength = 10, autoWidth = FALSE, scrollX = TRUE, scrollY = TRUE, #fixedHeader = TRUE,paging = FALSE,
                                   fixedColumns = list(leftColumns = 2,rightColumns = 0),
                                   columnDefs = list(list(width = '20px', targets = '_all'))))
      }
    }
    
  })
  

  
  output$msaResearch <- renderUI({
    
    gene_wise_info=read_delim(paste("gene_wise_info/",input$geneSelected,".txt",sep=''), "\t", escape_double = FALSE, trim_ws = TRUE)
    shinyjs::js$runMsa(containerId=msaContainerId, aaInputId="aa_Selected_aaFeat", 
                       sequence=paste(gene_wise_info$`Amino Acid`, collapse=""),
                       pathogenic=gene_wise_info$`Pathogenic Mutation`,
                       population=gene_wise_info$`Population Mutation`)
  })
  
})
