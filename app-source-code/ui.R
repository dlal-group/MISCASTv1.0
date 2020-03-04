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

## java script -- start ##
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

jscode <- "
shinyjs.runMolart = function(params) { 
runMolart(params);
}

shinyjs.runMsa = function(params) { 
runMsa(params);
}

shinyjs.disableTab = function(name) {
  var tab = $('.nav li a[data-value=' + name + ']');
tab.bind('click.tab', function(e) {
e.preventDefault();
return false;
});
tab.addClass('disabled');
}
"

#a global variable
reload_var <- FALSE

molartContainerId <- "molartContainer"
msaContainerId <- "msaContainer"
## java script -- end ##

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


navbarPageWithInputs <- function(..., inputs) {
  navbar <- navbarPage(...)
  form <- tags$form(class = "navbar-form", inputs)
  navbar[[3]][[1]]$children[[1]] <- htmltools::tagAppendChild(
    navbar[[3]][[1]]$children[[1]], form)
  navbar
}
# Define UI for application that draws a histogram
shinyUI(fluidPage(title = "MISCAST",
                        theme = "styles.css",
                        tags$head(tags$style("#UBIQstat{background-color: White; text-align: center;}div.box-header {
                                             text-align: center;
                                             }")),
  tags$head(tags$style("#ASAstat{background-color: White;text-align: center;}")),
  tags$head(tags$style("#RSAstat{background-color: White;text-align: center;}")),
  tags$head(tags$style("#ACETstat{background-color: White;text-align: center;}")),
  tags$head(tags$style("#GLCNstat{background-color: White;text-align: center;}")),
  tags$head(tags$style("#GALNstat{background-color: White;text-align: center;}")),
  tags$head(tags$style("#PHOSstat{background-color: White;text-align: center;}")),
  tags$head(tags$style("#SUMOstat{background-color: White;text-align: center;}")),
  tags$head(tags$style("#METHstat{background-color: White;text-align: center;}")),
  
  
  #Molart
  tags$head(tags$script(src="js/molart.js", type="text/javascript")),
  tags$head(tags$script(src="js/molart_controller.js", type="text/javascript")),
  tags$head(tags$script(src="js/color.min.js", type="text/javascript")),
  
  #msa.biojs.net
  tags$head(tags$script(src="js/msa.min.gz.js", type="text/javascript")),
  tags$head(tags$script(src="js/msa_controller.js", type="text/javascript")),
  #tags$head(tags$script(type="text/javascript", src="js/custom.js")),

  
  useShinyjs(),
  shinyjs::extendShinyjs(text = jsResetCode),
  shinyjs::extendShinyjs(text = jscode), #script parameter is not working (can't find the js file)
  
  #I have changed navbarPage with navbarPageWithInputs
  navbarPageWithInputs(id = "upperlayer", tags$b("MISCAST"),
             tabPanel(value = "page1",
                      tags$p("Home"),
                      fluidRow(
                        column(12, align = 'center',
                               HTML(paste("<h1><img src = 'miscast.png', height = 100, width = 600></h1>"))
                               #tags$img(src = 'miscast.png', height = 100, width = 500, align = 'center')
                        )
                      ),
                      fluidRow(
                        HTML(paste("<h1>Welcome to MIssense variant to protein StruCture Analysis web SuiTe</h1>"))
                      ),
                      fluidRow(
                        column(2),
                        column(8, align = 'center',
                          HTML(
                            paste(
                              "<br><br>",
                              "<color = Black face = Arial align = justify>", "<div align='justify'>",
                              "<b>MISCAST</b> is developed at the ",
                              "<a href=", Broad_link, " target=_blank>","Broad Institute","</a>","of ",
                              "MIT and Harvard, by a combined effort from the Genetics and Therapeutics group of",
                              "<a href=", SC_link, " target=_blank>","Stanley Center</a>.",
                              "The goal of MISCAST is to visualize and analyze single amino-acid-altering missense variants on protein sequence and 3-dimensional structure, and thereby forecast their biological impact.",
                              "<br><br>",
                              #"The dataset provided on this website spans 1,330 human genes, including 406,449 population variants from <a href = ",gnomAD_link," target=_blank>gnomAD</a> (release 2.1.1), ",
                              #"54,137 pathogenic variants and disease mutations from the <a href= ",ClinVar_link," target=_blank>ClinVar</a> (February, 2019 release) and <a href= ",HGMD_link," target=_blank>HGMD</a>, (Professional release 2018.4) databases, respectively,", "and >14k molecularly-solved human protein 3D structures from the Protein Data Bank (<a href=",PDB_link," target=_blank>PDB</a>). ",
                              "The website provides two primary tracks. The <strong>Variant Summary Report</strong> track provides an amino-acid-wise report, summarizing its pathogenic and population variant-associated protein features.",
                              "The <strong>Variant Analysis Suite</strong> provides a platform for detailed exploration, analysis, and visualization of protein features and missense variants on 1D, 2D and 3D protein structure space.",
                              "<br><br>",
                              "<strong>Disclaimer</strong>: <i>Please note that all content on this website is provided for academic and research purposes only. It has not undergone clinical validation, and should not be used for medical advice, diagnosis, or treatment. For best experience, please open Variant Analysis Suite on a laptop/desktop instead of any smaller screen.</i>",
                              "</div>"
                            )
                          )
                        ),
                        column(2)
                      ),
                      br(),
                      fluidRow(
                        column(2),
                        column(3, align = 'center',
                               actionButton(inputId = "reportB", 
                                            label = HTML(paste("<font size = +2><strong>Variant Summary Report</strong></font>")),style = "color: white; 
                                            background-color: #0059b3; width: 320px;")
                               )
                        ,
                        column(3, offset = 2, align = 'center',
                               actionButton(inputId = "researchB", 
                                            label = HTML(paste("<font size = +2><strong>Variant Analysis Suite</strong></font>")),style = "color: white; 
                                            background-color: #0059b3; width: 320px;")
                        ),
                        column(1)
                      ),
                      br(),
                      fluidRow(
                        column(3),
                        column(6, align = 'center',
                               tags$a(href = 'https://www.broadinstitute.org/', target='_blank',
                                      tags$img(src = 'BroadInstLogoforDigitalRGB.png', height = 80, width = 300, align = 'center') 
                                     ),
                          HTML(paste("<p><span style='font-size: 10pt; color: #000000;'>Copyright &copy; Broad Institute. All rights reserved.</span></p>"))
                        ),
                        column(3)
                      )
             ),
             tabPanel(value = "Documentation",tags$p("Documentation"),
                      mainPanel(style = "background-color: White;",
                                htmlOutput("doc")    
                      )
                      
             ),
             tabPanel(value = "About",tags$p("About"), 
                      mainPanel(style = "background-color: White;",
                                htmlOutput("about_page")    
                      )
             ),
             
             tabPanel(value = "singleTrack",
                      tags$p("Tracks"),
                      " ",
                      div(id = "reportTrackHide", fluidPage(
                        id = "reportTrack",
                        mainPanel(width = 12,
                                  fluidRow(
                                    column(3),
                                    column(width = 2, align = 'center', offset = 0, style='padding:6px;', h3("Select Gene:")),
                                    column(width = 3,  offset = 0, style='padding:0px;',
                                           selectizeInput(inputId = "reportTgeneSelected", label = "", choices = gene_1330_cv,
                                                          options = list(placeholder = "Gene Name",
                                                                         maxOptions = 3000)
                                           )),
                                    column(1, offset = 0, align = 'center', style='padding:23px;',
                                           actionButton(inputId = "reportTsubmit", 
                                                        label = "Submit",style = "color: white; 
                                                        background-color: #0059b3;font-size : 16px;")
                                           )
                                  ),
                                  br(),
                                  fluidRow(
                                    #column(10, align = 'center', offset = 0, style='padding:2px;',
                                    div(id = "report_INFOwellpanel",
                                        wellPanel(style = "background-color: White;",
                                                  htmlOutput("reportTinformation"),
                                                  selectizeInput(inputId = "aa_Selected", label = 'Select Amino Acid Index',
                                                                 choices = "", width = '300px',
                                                                 options = list(
                                                                   placeholder = 'Please select an option below',
                                                                   onInitialize = I('function() { this.setValue(""); }')
                                                                 )),
                                                  #checkboxInput("checkbox", label = HTML(paste("<strong><i>", "<span style='font-size: 14pt; color: #ff0000;'>", "Please check: All content on this website is provided for academic and research purposes only. It has not undergone clinical validation, and should not be used for medical advice, diagnosis, or treatment.", "</span>", "</i></strong>")), value = FALSE),
                                                  div(id = "below_report_INFOwellpanel",
                                                      wellPanel(style = "background-color: White;",
                                                                #fluidRow(HTML(paste("<strong>","&nbsp;&nbsp;", "I) Amino acid information: ", "</strong>"))),
                                                                htmlOutput("reportTsummreportHeadlines"),
                                                                br(),
                                                                br(),
                                                                plotOutput("reportTsummreportPlot", width = "80%", height = "550px"),
                                                                htmlOutput("reportTsummreport")
                                                                
                                                                #h3("here!")
                                                      ),
                                                      fluidRow(HTML(paste("&nbsp;&nbsp;", "<strong>Disclaimer: </strong>", "<i>", "Please note that all content on this website is provided for academic and research purposes only. It has not undergone clinical validation, and should not be used for medical advice, diagnosis, or treatment.", "</i>")))
                                                  ) %>% shinyjs::hidden()
                                                  
                                        )
                                    ) %>% shinyjs::hidden()
                                  )
                        )
                      ))%>% shinyjs::hidden(),#fluidpage1
                      div(id = "researchTrackHide", fluidPage(
                        id = "researchTrack",
                        br(),
                        fluidPage(
                                  fluidRow(
                                    HTML(paste("<h1>Explore Protein Features of Missense Variants by Gene</h1>"))
                                  ),
                                  br(),
                                  br(),
                                  fluidRow(
                                    column(3),
                                    column(6, align = 'center',
                                           wellPanel(
                                             style = "background-color: White;",
                                             fluidRow(
                                               column(12, align = 'center',
                                                      tabsetPanel(
                                                        tabPanel("Select a Gene",
                                                                 selectizeInput(inputId = "geneSelected", label = "", choices = c("", gene_1330_cv),
                                                                                options = list(placeholder = "Gene Name",
                                                                                               maxOptions = 3000)
                                                                 )
                                                        )
                                                        # ,
                                                        # tabPanel("Protein Class",
                                                        #          selectizeInput(inputId = "pclassNameselected" ,
                                                        #                         label = "",
                                                        #                         choices = c(
                                                        #                           "","Receptor",
                                                        #                           "Signaling Molecule",
                                                        #                           "Kinase",
                                                        #                           "Phosphatase",
                                                        #                           "Protease",
                                                        #                           "Enzyme Modulator",
                                                        #                           "Calcium-Binding Protein",
                                                        #                           "Transcription Factor",
                                                        #                           "Nucleic Acid Binding",
                                                        #                           "Transporter",
                                                        #                           "Transfer/Carrier Protein",
                                                        #                           "Cell Adhesion Molecule",
                                                        #                           "Cytoskeletal Protein",
                                                        #                           "Extracellular Matrix Protein",
                                                        #                           "Cell Junction Protein",
                                                        #                           "Oxidoreductase",
                                                        #                           "Transferase",
                                                        #                           "Hydrolase",
                                                        #                           "Lyase",
                                                        #                           "Isomerase",
                                                        #                           "Ligase",
                                                        #                           "Defense/Immunity protein",
                                                        #                           "Membrane Traffic Protein",
                                                        #                           "Chaperone"
                                                        #                           
                                                        #                         ),
                                                        #                         options = list(placeholder = "Protein Class")
                                                        #          )
                                                        #          
                                                        # )
                                                        ,
                                                        id = "homeSideBarTabSetPanel"
                                                      )
                                               )
                                             )
                                             
                                             ,
                                             actionButton(inputId = "submit", 
                                                          label = "Submit",style = "color: white; 
                                                          background-color: #0059b3; "),
                                             style = "width: 500px;"
                                             )
                                    ),
                                    column(3)
                                  )
                        )
                        
                      ))%>% shinyjs::hidden(),#fluidpage2
                  
                      div(id = "queryHide", fluidPage(
                        tabsetPanel(id = 'mainTabset',
                                    tabPanel( value = "Information","Information",
                                              br(),
                                              br(),
                                              fluidRow(
                                                htmlOutput("Information2")
                                              )
                                              
                                              
                                    ),
                                    tabPanel(
                                      value = "1D visualization", "1D visualization",
                                      mainPanel(style = "background-color: White;",
                                                htmlOutput("aaWiseFeatureTableGeneName")
                                      ),
                                      fluidRow(
                                        div(id = "aaFeature_wellPanel",
                                            wellPanel(style = "background-color: White;",
                                                      htmlOutput("reportTinformation_aaFeat"),
                                                      fluidRow(
                                                        column(width = 12,
                                                               div(id=msaContainerId, class="msa"),
                                                               uiOutput("msaResearch")
                                                        ),
                                                        br()
                                                      ),
                                                      shinyjs::hidden(selectizeInput(inputId = "aa_Selected_aaFeat", label = 'Amino Acid Index',
                                                                                     choices = "", width = '300px',
                                                                                     options = list(
                                                                                       placeholder = 'Please select an option below',
                                                                                       onInitialize = I('function() { this.setValue(""); }')
                                                                                     ))),
                                                      div(id = "below_report_aaFeat_wellPanel",
                                                          wellPanel(style = "background-color: White;",
                                                                    htmlOutput("reportTaafeature")
                                                                    #h3("here!")
                                                          )
                                                      ) %>% shinyjs::hidden()
                                                      
                                            )
                                        ) %>% shinyjs::hidden()
                                      )
                                    ),
                                    tabPanel(value = "2D visualization","2D visualization",
                                             mainPanel(style = "background-color: White;",
                                                       htmlOutput("TWODVisGeneName")
                                             ),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "3-class Secondary Structure", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(width = 8,offset = 2, 
                                                      withSpinner(plotlyOutput("SS3loli", height = 500, width = 900
                                                      )
                                                      )
                                               )                                            
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "8-class Secondary Structure", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("SS8loli", height = 500, width = 900))
                                                      )
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "Accessible Surface Area", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             br(),
                                             br(),
                                              
                                             fluidRow(
                                               column
                                               (width = 8, offset = 2, 
                                                      withSpinner(plotlyOutput("asaloli", height = 500, width = 900))
                                               )
                                             ),
                                             br(),
                                             fluidRow(
                                               column(5),
                                               column(6,htmlOutput("asaplotSummary"))
                                             ),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "Twenty Amino Acids and Physiochemical Property", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(width = 8, offset = 2, 
                                                      withSpinner(plotlyOutput("chemloli", height = 500, width = 900
                                                      ))
                                               )
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "Protein-protein Interaction types", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             br(),
                                             fluidRow(
                                               column(width = 8, offset = 2, 
                                                      withSpinner(plotlyOutput("bondloli", height = 500, width = 900)) 
                                               )
                                             ),
                                             div(id = "PTMtext",br(),br(),br(),fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "Distance to Post-translational Modification (PTM) Types", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ))%>% shinyjs::hidden(),
                                             
                                             div(id="Acet",br(),br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("distloliacet", height = 500, width = 900))
                                               )
                                             )) %>% shinyjs::hidden(),
                                             
                                             div(id="Met",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("distlolimet", height = 500, width = 900))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             
                                             div(id="Ogcl",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("distloliogcl", height = 500, width = 900))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             
                                             div(id="Phos",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("distloliphos", height = 500, width = 900))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             
                                             div(id="Sumo",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("distlolisumo", height = 500, width = 900))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             
                                             div(id="Ubiq",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("distloliubiq", height = 500, width = 900))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             div(
                                               id = "ptmcaption",br(),
                                               fluidRow(
                                                 column(3),
                                                 column(8,
                                                        htmlOutput("ptmplotsummary")
                                                 )
                                               )
                                             )%>% shinyjs::hidden(),
                                             
                                             div(id="uniprottext",br(),fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "UniProt-based Functional Features", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ))%>% shinyjs::hidden(),
                                             div(id="func_site",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("funclolifs", height = 500, width = 900
                                                      ))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             div(id="mol_proc",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("funclolimolp", height = 500, width = 900
                                                      ))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             div(id="func_bind",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("funclolifbr", height = 500, width = 900
                                                      ))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             div(id="seq_mot",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("funclolismr", height = 500, width = 900
                                                      ))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             div(id="mod_dom",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("funclolimd", height = 500, width = 900
                                                      ))
                                               )
                                             ))%>% shinyjs::hidden(),
                                             div(id="mod_res",br(),fluidRow(
                                               column(width = 8, offset = 2,
                                                      withSpinner(plotlyOutput("funclolimr", height = 500, width = 900
                                                      ))
                                               )
                                             ))%>% shinyjs::hidden()
                                    ),
                                    
                                    tabPanel(value = "3D visualization","3D visualization",
                                             mainPanel(style = "background-color: White;",
                                                       htmlOutput("THREEDVisGeneName")
                                             ),
                                             fluidRow(
                                               column(width = 12,
                                                      div(id=molartContainerId, class="molart"),
                                                      uiOutput("molart")
                                               )
                                             )
                                    ),
                                    tabPanel(value = "Protein class wise characterization of variants","Protein class wise characterization of variants",
                                             mainPanel(style = "background-color: White;",
                                                       htmlOutput("pClassVisGeneName")
                                             ),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "3-class Secondary Structure", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             fluidRow(
                                               br(),
                                               column(6, withSpinner(plotOutput("SS3bar", height = 450, width = 550))),
                                               column(6, withSpinner(plotOutput("SS3forest", height = 500, width = 600))),
                                               br()
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "8-class Secondary Structure", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             fluidRow(
                                               br(),
                                               column(6, withSpinner(plotOutput("SS8bar", height = 500, width = 600))),
                                               column(6, withSpinner(plotOutput("SS8forest", height = 500, width = 600))),
                                               br()
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "Acessible Surface Area/Residue Exposure", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             fluidRow(
                                               br(),
                                               column(6, withSpinner(plotOutput("ASAbox", height = 500, width = 600))),
                                               column(6, withSpinner(plotOutput("RSAforest", height = 500, width = 600))),
                                               br()
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "Amino Acids and Physiochemical Property", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             fluidRow(
                                               br(),
                                               column(6, withSpinner(plotOutput("chemPropbar", height = 500, width = 600))),
                                               column(6, withSpinner(plotOutput("chemPropforest", height = 500, width = 600))),
                                               br()
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "Protein-protein Interactions", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             fluidRow(
                                               br(),
                                               column(6, withSpinner(plotOutput("BONDbar", height = 500, width = 600))),
                                               column(6, withSpinner(plotOutput("BONDforest", height = 500, width = 600))),
                                               br()
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "Post-translational Modifications (PTMs)", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             fluidRow(
                                               br(),
                                               column(6, withSpinner(plotOutput("PTMbox", height = 730, width = 630))),
                                               column(6, withSpinner(plotOutput("PTMforest", height = 500, width = 630))),
                                               br()
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             br(),
                                             fluidRow(
                                               column(7, HTML(
                                                 paste(
                                                   "<font size = +2 color = #0059b3 face = Arial align = left>",
                                                   "&nbsp;", "UniProt-based Functional Features", "</font>"
                                                 )
                                               )),
                                               column(5)
                                             ),
                                             fluidRow(
                                               br(),
                                               column(6, withSpinner(plotOutput("upFuncBar", height = 500, width = 600))),
                                               column(6, withSpinner(plotOutput("upFuncForest", height = 500, width = 600))),
                                               br()
                                             ),
                                             br(),
                                             br(),
                                             br(),
                                             br()
                                    ),
                                    tabPanel(
                                      value = "Feature table","Feature table",
                                      
                                      fluidRow(
                                                htmlOutput("TablularViewGeneName")
                                                
                                      ),
                                      fluidRow(column(4),div(style="display:inline-block;vertical-align:top;width:100px;margin-left:10px;",tags$br(tags$h6(align = "left","Select Rows: "))),
                                               div(style="display:inline-block;vertical-align:top;",selectizeInput(inputId = "tableFilterSelected",label = "",choices = c("All amino acids","Amino acids with population variants","Amino acids with pathogenic variants")))
                                               ),
                                      fluidRow(column(4),
                                        div(style="display:inline-block;vertical-align:top;width:100px;margin-left:10px;",downloadButton("downloadData", "Download"))
                                      ),
                                      br(),
                                      fluidRow(
                                        sidebarPanel(
                                          checkboxGroupInput("show_columns_of_Feature_Table","Select Columns:",
                                                             choices = c("Amino Acid Index",
                                                                         "Amino Acid",
                                                                         "Population Mutation",
                                                                         "Population Mutation (details)",
                                                                         "Pathogenic Mutation",
                                                                         "Pathogenic Mutation (details)",
                                                                         "Bond Detail(*)",
                                                                         "Accessible Surface Area",
                                                                         "Relative Accessible Surface Area",
                                                                         "8-class DSSP Secondary Structure Properties",
                                                                         "3-class DSSP Secondary Structure Properties",
                                                                         "Exposure",
                                                                         "Physiochemical Property",
                                                                         "Acetylation Detail(*)",
                                                                         "Methylation Detail(*)",
                                                                         "O.GclNAc Detail(*)",
                                                                         "Phosphorylation Detail(*)",
                                                                         "Sumoylation Detail(*)",
                                                                         "Ubiquitination Detail(*)",
                                                                         "Active Site",
                                                                         "Metal Binding Site",
                                                                         "binding Site",
                                                                         "Site",
                                                                         "Zinc Finger",
                                                                         "Dna Binding",
                                                                         "Nucleotide Phosphate Binding",
                                                                         "Calcium Binding",
                                                                         "Region",
                                                                         "Repeat",
                                                                         "Coiled Coil",
                                                                         "Topological Domain",
                                                                         "Peptide",
                                                                         "Modified Residue",
                                                                         "Motif",
                                                                         "Domain",
                                                                         "Intramembrane",
                                                                         "Propetide",
                                                                         "Lipidation",
                                                                         "Glycosylation",
                                                                         "Cross Links",
                                                                         "Disulfide Bond",
                                                                         "Transit Peptide",
                                                                         "Transmembrane",
                                                                         "Signal Peptide",
                                                                         "Functional Site",
                                                                         "Functional Binding Region",
                                                                         "Sequence Motif Region",
                                                                         "Modular Domain",
                                                                         "Molecular Processing",
                                                                         "Modified Residues"
                                                             ),
                                                             selected = c("Amino Acid Index",
                                                                          "Amino Acid",
                                                                          "Population Mutation",
                                                                          "Pathogenic Mutation")
                                          )
                                        ),
                                        column(8,div(DT::dataTableOutput("full_gene_wise_info"),style = "font-size:70%")),
                                        column(8,div(DT::dataTableOutput("gnomad_gene_wise_info"),style = "font-size:70%")),
                                        column(8,div(DT::dataTableOutput("patient_gene_wise_info"),style = "font-size:70%"))
                                        ,
                                        br(),
                                        fluidRow(
                                          column(4),
                                          column(7,htmlOutput("FeatureTableSummary"))
                                        )
                                      )
                                    ),
                                    tabPanel(
                                      value = "Index table","Index table",
                                      
                                      #mainPanel(style = "background-color: White;",
                                      #          htmlOutput("IndexTablularViewGeneName")
                                                
                                      #),
                                      fluidRow(
                                        htmlOutput("IndexTablularViewGeneName")
                                      ),
                                      fluidRow(column(4),div(style="display:inline-block;vertical-align:top;width:100px;margin-left:10px;",tags$br(tags$h6(align = "left","Select Rows: "))),
                                               div(style="display:inline-block;vertical-align:top;",selectizeInput(inputId = "indextableFilterSelected",label = "",choices = c("All amino acids","Amino acids with structure")))
                                               #column(1,radioButtons("indexTablefiletype", "File type:",choices = c("csv", "txt")))
                                               ),
                                      fluidRow(column(4),
                                        div(style="display:inline-block;vertical-align:top;width:100px;margin-left:10px;",downloadButton("indexTabledownloadData", "Download"))
                                      ),
                                      fluidRow(
                                        sidebarPanel(
                                          checkboxGroupInput("show_columns_of_Index_Table","Select Columns:",choices = c(
                                            "Amino Acid Index in Protein Sequence",
                                            "Amino Acid",
                                            "Genomic Index of the Amino Acid",
                                            "Structure Index of the Amino Acid"),
                                            selected = c(
                                              "Amino Acid Index in Protein Sequence",
                                              "Amino Acid",
                                              "Genomic Index of the Amino Acid",
                                              "Structure Index of the Amino Acid"
                                            )
                                          )
                                        )
                                        ,
                                        column(8,div(DT::dataTableOutput("structure_gene_wise_info"),style = "font-size:70%")),
                                        column(8,div(DT::dataTableOutput("structure_all_gene_wise_info"),style = "font-size:70%"))
                                        
                                      ),
                                      br(),
                                      fluidRow(
                                        column(4),
                                        column(7,htmlOutput("IndexTableSummary"))
                                      )
                                    )
                        )
                      ))%>% shinyjs::hidden()#query page
                    )
             ,
             
             inputs = fluidRow(
                               div(style="display: inline-block;vertical-align:top; width: 250px;",selectInput(inputId = "trackInput",label = NULL, choices = c("","Variant Summary Report","Variant Analysis Suite"))),
                               div(style="display: inline-block;vertical-align:top; width: 100px;",actionButton(inputId = 'trackButton',label = "Select Track", padding='12px', style = "color: white; background-color: #0059b3;font-size : 14px;"))
             )
  )
        
            
             
  )
  
  
)