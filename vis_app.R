library(shiny)
library(shinyjs)
library(shinythemes)
library(plotly)
library(stringr)
library(ggplot2);  theme_set(theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=11, color="black"),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black", size=0.75),
                                   axis.ticks = element_line(color = "black", size=0.75),
                                   axis.text = element_text(size=6, color="black"), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), 
                                   axis.text.y = element_text(size = 12) , legend.title = element_text(size = 11),
                                   legend.text = element_text(size = 12)))
options(stringsAsFactors=FALSE)


# SINGLE CELL TUMOR DATA #################################

#print(getwd())
print("Initial load: single cell data...")

# colors
subcloneCols <- c("green4","red4","#527cde","orange2")
names(subcloneCols) <- c("Normal","C1","C2","C3")
varColors <- viridisLite::inferno(50)
exprColors <- viridisLite::cividis(50)

# single cell cohort annotation
snRnaSeqAnn<- read.delim("data/G34_MYC_N_PRDM6_samples_ann.detailed.txt")
rownames(snRnaSeqAnn) <- snRnaSeqAnn$Sample
snRnaSeqAnn <- snRnaSeqAnn[ snRnaSeqAnn$Status == "Primary",]
#print(rownames(snRnaSeqAnn))

subcloneRna <- read.delim("data/RNA_clones_info.txt")
rownames(subcloneRna) <- subcloneRna$Sample
subcloneRna$Sample <- NULL

print("Load snRNA-seq UMAP results...")
snRnaSeqInput <- list()
rnaInput <- list.files("data/snRNAseq")
for(fName in rnaInput ) {
    print(fName)
    sId = str_split(fName, "\\.")[[1]][1]
    snRnaSeqInput[[sId]] = paste0("data/snRNAseq/",fName)
}


print("Load snATAC-seq UMAP results...")
snAtacSeqInput <- list()
atacInput <- list.files("data/snATACseq")
for(fName in atacInput ) {
  print(fName)
  sId = str_split(fName, "\\.")[[1]][1]
  snAtacSeqInput[[sId]] = paste0("data/snATACseq/",fName)
}


# RNAseq UMAP annotation
rnaAnnTypes <- c("Subclone","Group34_A","Group34_B","Group34_C",
              "MYC","MYCN","SNCAIP","PRDM6")
# atacSeqseq UMAP annotation
atacAnnTypes <- c("Subclone","MYC","MYCN","SNCAIP","PRDM6")


# WGS TUMOR DATA #################################

message("Initial load: whole genome sequencing data...")

# colors
time.colors <- rep(c("#4FB12B", "#176A02"), 2)
names(time.colors) <- c("MRCA", "ECA", "non-amplified", "amplified")

group.colors <- c(`MB, G3` = "goldenrod1", `MB, G4` = "darkgreen", `MB, SHH CHL AD` = "firebrick",
                  `MB, WNT` = "blue", `MB, SHH INF` = "red")
g34.subgroup.colors <- c(MB_G34_I = "#cd89d9", MB_G34_II = "#A34D7C", MB_G34_III = "#C4A908", 
                         MB_G34_IV = "#FFFF00",
                         MB_G34_V = "#ADCA02", MB_G34_VI = "#89B395", MB_G34_VII = "#9BC2E6", MB_G34_VIII = "#127134")

cbp1 <- rev(c("#999999", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

sig.colors <- c(SBS1 = "#FF0000", SBS3 = "#00A08A", SBS5 = "#F2AD00",
                SBS8 = "#F98400", SBS10a = "#5BBCD6", SBS10b = "#ECCBAE", Clock = "#046C9A")

time.shapes <- c("< ECA" = 0, "ECA" = 15, "> ECA, < MRCA" = 4, "MRCA" = 16, "> MRCA" = 1)
  
# whole genome sequencing cohort annotation
wgsAnn <- read.delim("data/WGS_sample_information.txt")
rownames(wgsAnn) <- wgsAnn$MRCA_ID
tumors.wgs <- sort(rownames(wgsAnn))

## CNV results
message("Load CNVs called from WGS data")
WGSCNVs <- list()
wgsInput <- list.files("data/WGS/CNV/", pattern = "CNV_plot", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  WGSCNVs[[sId]] = fName
}

## Variant allele frequencies 
message("Load variant allele frequencies")
WGSVafs <- list()
wgsInput <- list.files("data/WGS/VAF/", pattern = "VAFs", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  WGSVafs[[sId]] = fName
}
WGSVafFits <- list()
wgsInput <- list.files("data/WGS/VAF/", pattern = "Fit_VAF_distr", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  WGSVafFits[[sId]] = fName
}

## Signature results
message("Load mutational signature results")
WGSSigs <- list()
wgsInput <- list.files("data/WGS/MutationalSignatures/", pattern = "Mutational_signature", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.txt")[[1]][1]
  sId = gsub(".*_", "", sId)
  WGSSigs[[sId]] = fName
}

## Mutational densities
message("Load mutational densities")
WGSMutdensities <- list()
wgsInput <- list.files("data/WGS/MRCA_ECA/", pattern = "Mutational_densities", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  WGSMutdensities[[sId]] = fName
}
EventTimelines <- list()
wgsInput <- list.files("data/WGS/MRCA_ECA/", pattern = "Timeline", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  EventTimelines[[sId]] = fName
}
MRCAdensities <- read.delim("data/MRCA_densities.txt")
ECAdensities <- read.delim("data/ECA_densities.txt")

## Comparison to MutationTimeR
message("Load MutationTimeR")
WGSMutTimeR <- list()
wgsInput <- list.files("data/WGS/MutationTimeR/", pattern = "Lachesis_", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  WGSMutTimeR[[sId]] = fName
}
WGSEvoMutTimeR <- list()
wgsInput <- list.files("data/WGS/MutationTimeR/", pattern = "Timeline_", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  WGSEvoMutTimeR[[sId]] = fName
}

## Mobster result
message("Load Mobster results")
Mobster <- list()
wgsInput <- list.files("data/WGS/Mobster/", pattern = "Mobster_", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  Mobster[[sId]] = fName
}

## Neutral evolution fits
message("Load neutral evolution fits")
NeutralFits <- list()
wgsInput <- list.files("data/WGS/Neutral_fits/", pattern = "Neutral_fit", full.names = T)
for(fName in wgsInput ) {
  sId = str_split(fName, "\\.")[[1]][1]
  sId = gsub(".*_", "", sId)
  NeutralFits[[sId]] = fName
}

# default start #################################

curId <- "MB272"
print(snRnaSeqInput[[curId]])
dmRes <- read.delim(snRnaSeqInput[[curId]])
print(colnames(dmRes))
atacRes <- read.delim(snAtacSeqInput[[curId]])
print("Loaded default info")

######### COOKIES CONTROL #################


# addResourcePath("js", "www")
# 
# jsCode <- '
#   shinyjs.getcookie = function(params) {
#     var cookie = Cookies.get("terms_accepted");
#     if (typeof cookie !== "undefined") {
#       Shiny.onInputChange("jscookie", cookie);
#     } else {
#       var cookie = "";
#       Shiny.onInputChange("jscookie", cookie);
#     }
#   }
#   shinyjs.setcookie = function(params) {
#     Cookies.set("terms_accepted", escape(params), { expires: 10 });  
#     Shiny.onInputChange("jscookie", params);
#   }
#   shinyjs.rmcookie = function(params) {
#     Cookies.remove("terms_accepted");
#     Shiny.onInputChange("jscookie", "");
#   }
# '
# 


################ GUI ###########################


ui <- fluidPage(  
  # tags$head(
  #   tags$script(src = "js/js.cookie.js")
  # ),
  #useShinyjs(),
  #extendShinyjs(text = jsCode,functions=c("getcookie","setcookie","rmcookie")),
  theme=shinytheme("spacelab"),
  navbarPage(title = "MB G34 onco aberrations",
       tabPanel("Single nuclei RNA/ATAC",
             fluidRow(
                 column(width = 5,
                      fluidRow (
                        column( 6,
                           selectInput("inputScData", "Sample:",
                                       choices=rownames(snRnaSeqAnn),selected = curId),
                           selectInput("data_type", "Data type:", choices =NULL),
                           h5(htmlOutput("snTumGroup")),
                           h5(htmlOutput("snTumSubgroup")),
                           h5(textOutput("status")),
                           h5(textOutput("age")),
                           h5(textOutput("gender")),
                           h5(uiOutput("somatic")),
                           h5(uiOutput("subclones")),
                           selectInput("selGroup","UMAP annotation:",
                                       rnaAnnTypes,selected = "Subclone")
                     ),
                        column( 3,
                                h5("Phylogeny:"),
                                tags$div(
                                  uiOutput("image_cnv_tree",align = "center")
                                )
                        )
                     )
                  ), 
                 column(width = 7,
                        div(
                          h4("CNV bulk (methylation)", align = "left"),
                          tags$div (
                            style = "margin-left: 50px", offset = 0, # move to the left
                            uiOutput("image_cnv_bulk")
                          )
                        )
                 )
             ),
             fluidRow(
                 column(width = 5,
                     plotlyOutput("dimplot", height = "400px")
                 ),
                 column(width = 7,
                        div(
                          h4("CNV single nuclei ", align = "left"),
                          uiOutput("image_cnv_snRNA")
                        )
                 )
             )
       ),
       tabPanel("Whole genome sequencing", 
                fluidRow(
                  column(width = 4,
                         selectInput(inputId = "inputWGSData", label = "Sample:",
                                     choices=tumors.wgs, selected = curId),
                         h5(htmlOutput("mnp11")),
                         h5(htmlOutput("mnp12")),
                         h5(textOutput("ageWGS")),
                         h5(textOutput("genderWGS")),
                         radioButtons('EvoAna', 'Evolutionary analysis:', choices = "initiation"),
                         selectInput(inputId = "selCN","Copy number",
                                     choices=NULL)
                  ),  
                  column(width = 8,
                         div(
                           h5("Copy number profile", align = "left"),
                           tags$div (
                            # plotOutput("cnvplot", height = "300px")
                             uiOutput("image_cnvplot", height = "300px")
                           )
                         )
                  )
                ),
                fluidRow(
                  column(width = 4,
                         div(
                           h5("Mutational signatures", align = "left"),
                           tags$div (
                             plotOutput("mutsigplot", height = "300px")
                           )
                         )
                  ),
                  column(width = 5,
                         div(
                           h5("Single nucleotide variants", align = "left"),
                           tags$div (
                             plotOutput("vafplot", height = "250px")
                           )
                         )
                  ),
                  column(width = 3,
                         tags$div (
                           style = "margin-top: 50px", offset = 0, # move to the left
                           uiOutput("EvoChoices")
                         )
                         )
                  ),
                uiOutput("EvoPlots")
       ),
       tabPanel("Help & About",
                fluidRow(
                  column(8,
                         h3("About the tool"),
                         p("This ShinyApp focuses on the interactive inspection of somatic changes of medulloblastoma 
                            Group 3/4 brain tumors. For more details please refer to the original manuscript:"),
                            a("Okonechnikov, Joshi, Körber et al \"Medulloblastoma oncogene aberrations are not 
                              involved in tumor initiation, but essential for disease progression and therapy resistance\""  , href = "https://www.biorxiv.org/content/10.1101/2024.02.09.579690", target = "_blank"),
                         p(""),
                         p("The application consists of two main parts."), 
                         h4("Single nuclei RNA/ATAC"),
                         p("The first part is the visualization of copy number variance in single cell resolution 
                            from Group 3/4 cases with specific somatic changes: MYC/MYCN ampflications and SNCAIP/PRDM6 
                           translocations."),
                         p("Target sample can be selected from the box \"Sample\", while the single nuclei data type  
                           from \"Data type\". Details of a target tumor tissue sample are futher provided below."),
                         p("The visualization plots include UMAP (bottom-left), control bulk CNV profile  obtained from  
                         methylation array (top-right) and single cell CNV profiles computed with InferCNV (bottom-right)."),
                         p("UMAP interactive plot has specific vizualization options"),
                         tags$ul(
                              tags$li("Subclone annotation across cells (same colors are used to mark subclones in single cell CNV profile)"), 
                              tags$li("Only in snRNA-seq: enrichment of medulloblastoma group 3/4 proliferation (G34_A), progenitor-like (G34_B) and differentation (G34_C) signals"), 
                              tags$li("Gene expression (snRNA-seq) and gene body signal enrichment (snATAC-seq) for main targets")
    
                          ),
                         p(""),
                         h4("Whole genome sequencing"),
                         p("The second part visualizes results on medulloblastoma evolution that were obtained from deep whole genome sequencing data.
                           For each sample, the following information is shown:"),
                         tags$ul(
                           tags$li("meta data (medulloblastoma group, group3/4 subgroup, age at diagnosis, gender)"),
                           tags$li("the allele-specific copy number profile"),
                           tags$li("the distribution of somatic single nucleotide variants (SNVs) explained by particular mutational signatures"),
                           tags$li("the variant allele frequency (VAF) distribution of SNVs, stratified by copy number, and the binomial density distributions fitted to the clonal VAF peaks"),
                           tags$li("densities of non-amplified and amplified SNVs per genomic segment; moreover, SNV densities at an early common ancestor (ECA) and at the most recent common ancestor (MRCA)"),
                           tags$li("a comparison of mutation densities to MutationTimeR (Gerstung et al., Nature, 2020); here, SNV densities at chromosomal gains relative to the SNV densities at MRCA are shown."),
                           tags$li("if subclonal evolution was analyzed for this sample:"),
                           tags$ul(
                             tags$li("subclonal deconvolution of cancer cell fractions with Mobster (Caravagna et al., Nature Genetics, 2020)"),
                             tags$li("copy numbers-stratified fits of a population-genetics model of variant accumulation in expanding tumors to the cumulative VAF distribution"),
                           )
                         ),
                         p("Target sample can be selected from the box \"Sample\". The box \"Evolutionary analysis\" allows to display plots focusing on tumor initiation (timeline of acquiring clonal chromosomal gains) 
                         or growth (plots showing tumor evolution as inferred from the VAF distribution of subclonal variants). If \"tumor initiation\" is chosen, the box \"Mutation density method\" allows to switch between results obtained
                         with the methodology described in this paper and with MutationTimeR (Gerstung et al., Nature, 2020). Finally, the box \"Copy number\" allows to select a copy number for which the VAF distribution and evolutionary fits are shown.
                           "),
                         p("")
                         
                  
                  )
                )
       ),
  ),
  # image on right side of header
  tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><a href=\"https://www.kitz-heidelberg.de/\"><img src=\"logo.png\" alt=\"alt\" style=\"float:right;width:94px;height:45px;padding-top:10px;\"> </a></div>');
    console.log(header)")
              
  ),
  tags$footer(tags$em("For questions & help:",tags$a(href="mailto:k.okonechnikov@kitz-heidelberg.de", "contact e-mail")), align = "center"),
)

server <- function(input, output,session) {
  
  ### INITIAL MESSAGE
  
  # query_modal <- modalDialog(
  #   title = "G34 Onco Aberrations Disclaimer & Terms of Use",
  #   p("Disclaimer"),
  #   p(" G34 Onco Aberrations ShinyApp focuses on  on interactive visualization 
  #   of somatic properties of medulloblastoma G34 tumors obtained from sequencing data from
  #   published bulk and single cell profiles. More details about datasets are provided 
  #   in additional application page “Data description”. ", style = "font-family: 'times'; font-si14pt"),
  #   p("Terms of Use"),
  #   tags$ul(
  #     tags$li("Users of this shiny-app are not able to upload files."),
  #     tags$li("This shiny-app is only used for visualization purposes, demonstration of already computed 
  #             results with additional visualization controls."),
  #     tags$li("Since this shiny-app is constantly being further developed, there may also be changes in the use itself.")
  #     , style = "font-family: 'times'; font-si14pt"),
  #       easyClose = F,
  #   footer = tagList(
  #     actionButton("close_dialog", "Accept Disclaimer & Terms of Use"),
  #     actionButton("decline","Decline")
  #   )
  # )
  # 
  # js$getcookie()
  # observeEvent(input$jscookie, {
  #   if (!is.null(input$jscookie) && (input$jscookie != "") ) {
  #     print("cookie found")
  #     print(input$jscookie)
  #     # clean
  #     # js$rmcookie()
  #   }
  #   else {
  #     print("cookie not found")
  #     showModal(query_modal)
  #   }
  # })

  # observeEvent(input$close_dialog, {
  #   js$setcookie("terms_accepted")
  #   removeModal()
  # })
  # 
  # observeEvent(input$decline, {
  #   stopApp()
  # })
 
  ### SINGLE CELL FOCUS
  
  getScTum <- eventReactive(input$inputScData, {
      #print("Update!")
      curId <<- input$inputScData
      dmRes <<- read.delim(snRnaSeqInput[[curId]])
      print(paste("Change dataset to",curId))
      #if (!is.null(selected_sample) && selected_sample %in% names(sample_data)) {
      dataTypes = c("snRNA-seq")
      #if (snRnaSeqAnn[curId,"snATAC"] == 1){
      if (curId %in% names(snAtacSeqInput)){
          atacRes <<- read.delim(snAtacSeqInput[[curId]])
          dataTypes  <- c(dataTypes,"snATAC-seq")
      } 
      updateSelectInput(session, "data_type", choices = dataTypes)
      
      curId
  })
  
  getAnnGroup <- reactive({
      input$selGroup
  })
  
  getDataType <- reactive({
    if (input$data_type == "snRNA-seq") {
        annTypes <- rnaAnnTypes      
    } else {
       annTypes <- atacAnnTypes
    }
    updateSelectInput(session, "selGroup", choices = annTypes,selected = "Subclone")
    input$data_type
  })
  
  output$snTumGroup <- renderUI({ 
      tId <- getScTum()
      HTML(paste0("Group: <b>",snRnaSeqAnn[tId,"Group"],"</b>"))
  })
 
  output$snTumSubgroup <- renderUI({ 
    tId <- getScTum()
    HTML(paste0("Subgroup: <b>",snRnaSeqAnn[tId,"Subgroup"],"</b>")) 
  })
  
  
  output$status <- renderText({ 
    tId <- getScTum()
    paste("Status:",snRnaSeqAnn[tId,"Status"]) 
  })
  
  output$age <- renderText({ 
    tId <- getScTum()
    paste0("Age (years): ",snRnaSeqAnn[tId,"Age"]) 
  })
  
  output$gender <- renderText({
    tId <- getScTum()
    paste0("Gender: ",snRnaSeqAnn[tId,"Gender"])
  })
  
  output$somatic <- renderText({
    tId <- getScTum()
    sIds <- colnames(snRnaSeqAnn)[9:11] 
    sIds[ snRnaSeqAnn[tId,9:11] == 1 ]
    paste0("Aberration: ",paste( sIds[ snRnaSeqAnn[tId,9:11] == 1 ],collapse = ","))
  })
  
  output$subclones <- renderText({
    tId <- getScTum()
    dataType <- getDataType()
    if (dataType == "snRNA-seq") {
      tRes <- dmRes[dmRes$Subclone != "Normal",]
    } else {
      tRes <- atacRes[ atacRes$Subclone != "Normal",]
    }
    combProps <- summary(as.factor(tRes$Subclone)) *100 / nrow(tRes)
    
    subcloneInfo <- c()
    for (i in 1:(length(combProps))) {
        targProps = combProps[i]
        sc <- names(targProps)
        cType <- subcloneRna[tId,sc]
        sc <- ifelse(cType == "", sc,paste0(sc,"-",cType))
        subcloneInfo = c(subcloneInfo, sprintf("%s~%.0f%%" , sc,targProps))
    }
    #print(subcloneInfo)
    paste("Subclones:",paste(subcloneInfo, collapse = ", "))
  })
  
  
  output$image_cnv_bulk <- renderUI({
    selected_id <-  getScTum()
    
    image_src <- paste0("cnvBulk/",selected_id,"_meth.png?version=1")
    if (!is.null(image_src)) {
      img_tag <- tags$img(src = image_src, width = "600px")
    } else {
      img_tag <- tags$p("No image selected.")
    }
    
    # Return the generated img tag
    img_tag
  })
  
  output$image_cnv_tree <- renderUI({
    selected_id <- getScTum()
    
    image_src <- paste0("cnvTree/",selected_id,"_tree.png?version=0")
    if (!is.null(image_src)) {
      img_tag <- tags$img(src = image_src, width = "180px")
    } else {
      img_tag <- tags$p("No image selected.")
    }
    
    # Return the generated img tag
    img_tag
  })
  
  
  output$image_cnv_snRNA <- renderUI({
    selected_id <- getScTum()
    dataType <- getDataType()
    if (dataType == "snRNA-seq") {
        image_src <- paste0("cnvSnRNAseq/",selected_id,"_RNA.png")
    } else {
        image_src <- paste0("cnvSnATACseq/",selected_id,"_ATAC.png")
    }
    if (!is.null(image_src)) {
      img_tag <- tags$img(src = image_src, width = "700px")
    } else {
      img_tag <- tags$p("No image selected.")
    }
    
    # Return the generated img tag
    img_tag
  })
  
 
  
    
  output$dimplot <- renderPlotly({
      tId <- getScTum()
      selGrp <- getAnnGroup()
      dataType <- getDataType()
      #print(paste(tId,selGrp,dataType))
      if (dataType == "snRNA-seq") {
        targRes <- dmRes
      } else {
        targRes <- atacRes
      }
      #print(dim(targRes))
      if(selGrp == "Subclone") {
           selColors = subcloneCols
      } else if (grepl("Group",selGrp) ){
          selColors = varColors
      } else {
           selColors = "YlGnBu"
      }
      #print(selColors)
      p <- plot_ly( x = targRes[,"UMAP_1"], y = targRes[,"UMAP_2"],
          color=targRes[,selGrp],
          #text = rownames(dmRes), what to show on move mouse
          colors=selColors,
          #marker=list(colorscale='Viridis'),
          size = I(3),
          mode="markers",
          #source = "subset",
          type="scattergl")
      
  })
  
  ### WGS FOCUS
  
  observeEvent(input$EvoAna == "initiation",{
    curId <<- input$inputWGSData
    if (curId %in% names(WGSMutTimeR)){
      updateRadioButtons(session, "MutDens", choices = c("this paper", "MutationTimeR"))
    } 
  })
  
  getWGSTum <- eventReactive(input$inputWGSData, {
    curId <<- input$inputWGSData
    ## load in all data for this tumor
    cnvsegments <<- readRDS(WGSCNVs[[curId]])
    vafFit <<- read.delim(WGSVafFits[[curId]]) 
    vaf <<- read.delim(WGSVafs[[curId]]) 
    sigs <<- read.delim(WGSSigs[[curId]]) 
    mutdens <<- read.delim(WGSMutdensities[[curId]]) 

    timeline <<- read.delim(EventTimelines[[curId]])
    
    analyses <- "initiation"
    if (curId %in% names(NeutralFits)){
      mobsterfit <<- readRDS(Mobster[[curId]]) 
      neutrfit <<- read.delim(NeutralFits[[curId]]) 
      analyses  <- c(analyses,"growth")
    } 
    updateRadioButtons(session, "EvoAna", choices = analyses)
    
    if (curId %in% names(WGSMutTimeR)){
      muttimer <<- read.delim(WGSMutTimeR[[curId]]) 
      muttimerevo <<- read.delim(WGSEvoMutTimeR[[curId]])
      updateRadioButtons(session, "MutDens", choices = c("this paper", "MutationTimeR"))
    } 
    cns <- sort(as.numeric(unique(vaf$CN)))
    sel.cn <- ifelse(2 %in% cns, 2, cns[1])
    updateSelectInput(session, "selCN", choices = cns,selected = sel.cn)
    
    message("Change dataset to ",curId)
    curId
  })
  
  
  getPurity <- reactive( {
    curId <<- input$inputWGSData
    purity <<- wgsAnn[curId,"Purity"]
    purity
  })
  
  output$mnp11 <- renderText({ 
    tId <- getWGSTum()
    type <- wgsAnn[tId,"mnp11"]
    type <- gsub(".*MB, ", "", type)
    HTML(paste0("Group: <b>",type,"</b>"))
  })
  
  output$mnp12 <- renderText({ 
    tId <- getWGSTum()
    type <- wgsAnn[tId,"mnp12"]
    type <- gsub(".*MB_", "", type)
    HTML(paste0("Subgroup: <b>",type, "</b>")) 
  })
  
  output$ageWGS <- renderText({ 
    tId <- getWGSTum()
    paste0("Age (years): ",wgsAnn[tId,"Age"]) 
  })
  
  output$genderWGS <- renderText({
    tId <- getWGSTum()
    paste0("Gender: ",wgsAnn[tId,"Gender"])
  })
  
  output$cnvplot <- renderPlot({
    tId <- getWGSTum()
    print(cnvsegments + theme(text = element_text(size = 12)))
  })
  
  output$image_cnvplot <- renderUI({
    selected_id <-  getWGSTum()
    
    image_src <- paste0("WGScnv/CNV_plot_",selected_id,".png")
    if (!is.null(image_src)) {
      img_tag <- tags$img(src = image_src, width = "550px")
    } else {
      img_tag <- tags$p("No image selected.")
    }
    
    # Return the generated img tag
    img_tag
  })
  
  output$mutsigplot <- renderPlot({
    tId <- getWGSTum()
    p <- ggplot(sigs, aes(x=Type, y=Fraction.of.SNVs, fill=Signature)) + geom_col() + scale_fill_manual(values=sig.colors) +
      scale_x_discrete(name="") + scale_y_continuous("Fraction of SNVs") + coord_polar("y", start = 0) +
      theme_void() + theme(legend.position = "bottom")
    print(p)
  })
  
  output$vafplot <- renderPlot({
    tId <- getWGSTum()
    purity <- as.numeric(getPurity())
    CN <- as.numeric(input$selCN)
    
    p <- ggplot(vaf[vaf$CN == input$selCN,,drop=F], aes(x=VAF)) + geom_histogram(binwidth = 1/vaf[vaf$CN == input$selCN,,drop=F]$avg.depth[1]) + 
      scale_x_continuous(limits=c(0,1), name="Variant allele frequency") + ggtitle(paste0("CN = ", CN)) +
      geom_vline(xintercept = seq(1,CN)*purity/(CN*purity + 2*(1-purity)), linetype = 2, size = 2.5) +
      scale_y_continuous(name="Number of SNVs") +
      geom_line(data=vafFit[vafFit$CN == input$selCN,,drop=F],
                aes(x=x, y=p, col = Peak, group = B), size=2.5) + scale_color_manual(values = time.colors) +
      guides(color = guide_legend(title = "Fitted distribution of clonal SNVs")) + theme(legend.position = "bottom")
    
    print(p)
  })
  
  output$nonampdensplot <- renderPlot({
    
    tId <- getWGSTum()
    binwidth = (max(mutdens$Mean) - min(mutdens$Mean))/20
    
    p <- ggplot(mutdens[mutdens$Timing=="Post-CNV",], aes(x=Mean)) + 
      geom_histogram( binwidth = binwidth, fill=time.colors["MRCA"])+ 
      scale_y_continuous(name="Number of genomic segments") + 
      scale_x_continuous(name="SNVs per Mb", limits = c(-0.05*(max(mutdens$Mean)), max(mutdens$Mean)*1.1))
    
    print(p)
   })
  
  output$ampdensplot <- renderPlot({
    
    tId <- getWGSTum()
    binwidth = (max(mutdens$Mean) - min(mutdens$Mean))/20
    
    p <- ggplot(mutdens[mutdens$Timing=="Pre-CNV",], aes(x=Mean)) + 
      geom_histogram( binwidth = binwidth, fill=time.colors["ECA"])+ 
      scale_y_continuous(name="Number of genomic segments") + 
      scale_x_continuous(name="SNVs per Mb", limits =  c(-0.05*(max(mutdens$Mean)), max(mutdens$Mean)*1.1)) 
    
    print(p)
  })
  
  output$timelineplot <- renderPlot({
    
    tId <- getWGSTum()
    mrca <- MRCAdensities[MRCAdensities$Sample==tId,]
    eca <- ECAdensities[ECAdensities$Sample==tId,]
    timeline$Segment <- factor(timeline$Segment, levels = timeline$Segment)
    timeline$Time <- factor(timeline$Time, levels = c("< ECA", "ECA", "> ECA, < MRCA", "ECA/MRCA", "MRCA", "> MRCA"))
    p <- ggplot(timeline, aes(x=Mean, xmin=Min, xmax=Max, y=as.numeric(as.factor(Segment)) - 
                                0.5, col = Segment)) + scale_shape_manual(values = time.shapes, drop = F) +
      geom_point(aes(shape=Time), size = 3, show.legend = T) + geom_errorbarh(height=0, size = 1)+ 
      geom_vline(data=data.frame(x=mrca["Mean"]), aes(xintercept=Mean/3.3/10^3),  col=time.colors["MRCA"],
                 linetype=2, size = 1) + 
      annotate(geom = "text", label = "MRCA \n(mean; 95% CI)", x = unlist(mrca["Max"]/3.3/10^3), hjust = 0, 
               y = nrow(timeline)+1, col=time.colors["MRCA"]) +
      geom_ribbon(data=data.frame(xmin=rep(unlist(mrca["Min"]),2),
                                  xmax=rep(unlist(mrca["Max"]),2),
                                  y=c(-2, 2+max(1,nrow(timeline)))), 
                  aes(xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3,y=y),
                  inherit.aes = F, fill=time.colors["MRCA"], col=NA, alpha=0.5) +
      geom_vline(data=data.frame(x=eca["Mean"]), aes(xintercept=Mean/3.3/10^3), col=time.colors["ECA"],
                 linetype=2, size = 1) + 
      geom_ribbon(data=data.frame(xmin=rep(unlist(eca["Min"]), 2),
                                  xmax=rep(unlist(eca["Max"]),2),
                                  y=c(-2, 2+nrow(timeline))), 
                  aes( xmin=xmin/3.3/10^3, xmax=xmax/3.3/10^3, y=y),
                  inherit.aes = F, fill=time.colors["ECA"], col=NA, alpha=0.5) +
      annotate(geom = "text", label = "ECA \n(mean; 95% CI)", x = unlist(eca["Max"]/3.3/10^3), hjust = 0, 
               y = -1, col=time.colors["ECA"]) +
      scale_x_continuous(name="SNVs per Mb", limits =  c(-0.05*(max(mutdens$Mean)), max(mutdens$Max)*1.1)) +
      scale_y_continuous(limits=c(-2, 2+max(1,nrow(timeline)))) + 
      theme(axis.line.y=element_blank(), axis.text.y=element_blank(),
            axis.ticks.y=element_blank(), axis.title.y=element_blank()) +
      guides(col=guide_legend(ncol=2, title = "Genomic segment \n(chromosome - copy number - \nmajor allele count)"),
             shape = guide_legend(title = "Time of copy \nnumber change"))
    
    print(p)
  })
  
  output$timelineplotMutTimeR <- renderPlot({
    tId <- getWGSTum()
    muttimerevo$SegId <- factor(muttimerevo$SegId, levels = unique(muttimerevo$SegId))
    ggplot(muttimerevo, aes(x = Time, xmin = Time.lo, xmax = Time.up, y = SegNo, col = SegId)) + 
      geom_point(size = 3) + geom_errorbarh(height = 0, size = 1.5) + theme(axis.line.y = element_blank(), axis.title.y = element_blank(), 
                                                        axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
      geom_vline(xintercept = 1, col = time.colors["MRCA"], linetype = 2, size = 1.5) + 
      scale_x_continuous(name = "Mutation time relative to MRCA", limits = c(0, 2)) +
      scale_y_continuous(limits = c(-2, nrow(muttimerevo) +2)) +
      guides(col=guide_legend(ncol=2))+
      annotate(geom = "text", label = "MRCA \n(mean; 95% CI)", x = 1.1, hjust = 0, 
               y = nrow(muttimerevo)+1, col=time.colors["MRCA"]) +
      guides(col=guide_legend(ncol=2, title = "Genomic segment \n(chromosome - copy number - \nmajor allele count)"))
  })
  
  output$neutralevolplot <- renderPlot({
    tId <- getWGSTum()
    if(nrow(neutrfit[neutrfit$CN==input$selCN,,drop=F])>0){
      ggplot(neutrfit[neutrfit$CN==input$selCN,,drop=F], aes(x=VAF, y=Data, 
                                                             ymin = Data - sqrt(Data),
                                                             ymax = Data + sqrt(Data)
      )) + 
        geom_ribbon(aes(ymin=Mmin, 
                        ymax=Mmax),
                    fill="lightslateblue")+ geom_point() + geom_errorbar(width=0.01) + 
        scale_y_continuous(name="Cumulative number of SNVs") +
        ggtitle(paste("CN = ", input$selCN))
    }

  })
  
  output$MutationTimeRplot <- renderPlot({
    tId <- getWGSTum()
    p <- ggplot(muttimer, aes(x = Lachesis_mean, xmin = Lachesis_min, xmax = Lachesis_max,
                              y = MutationTimer_mean, ymin = MutationTimer_min, ymax = MutationTimer_max)) +
      geom_errorbar(width = 0) + geom_errorbarh(height = 0) + geom_point() + 
      geom_abline(slope = 1, intercept = 0, linetype = 2) + theme(aspect.ratio = 1) +
      scale_x_continuous(limits = c(0, max(c(muttimer$Lachesis_max, muttimer$MutationTimer_max))), name = "This paper")+
      scale_y_continuous(limits = c(0, max(c(muttimer$Lachesis_max, muttimer$MutationTimer_max))), name = "MutationTimeR")
    print(p)
  })
  
  output$mobsterplot <- renderPlot({
    tId <- getWGSTum()
    plot(mobsterfit$best) + ggtitle("") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=11, color="black"),
                                               panel.background = element_blank(), axis.line = element_line(colour = "black", size=1),
                                               axis.ticks = element_line(color = "black", size=1),
                                               axis.text = element_text(size=11, color="black"),
            line = element_line(size = 1.2))
  })
  
  output$image_mobsterplot <- renderUI({
    selected_id <-  getWGSTum()
    
    image_src <- paste0("mobster/Mobster_",selected_id,".png")
    if (!is.null(image_src)) {
      img_tag <- tags$img(src = image_src, width = "250px")
    } else {
      img_tag <- tags$p("No image selected.")
    }
    
    # Return the generated img tag
    img_tag
  })
  
  
  output$MutDensPlots <- renderUI( if(input$MutDens == "this paper"){
    fluidRow(
      column(width = 3, 
             div(
               h5("Densities of non-amplified SNVs", align = "left"),
               tags$div (
                 plotOutput("nonampdensplot", height = "250px")
                 
               )
             )
            ),
      column(width = 3,
             div(
               h5("Densities of amplified SNVs", align = "left"),
               tags$div (
                 plotOutput("ampdensplot", height = "250px")
               )
             ),
      ),
      column(width = 6,
             div(
               h5("Evolutionary timeline", align = "left"),
               tags$div (
                 plotOutput("timelineplot", height = "250px")
               )
             )
      )
    )
    
  }else if(input$MutDens == "MutationTimeR"){
    fluidRow(
      column(width = 4, 
             div(
               h5("Mutation time relative to MRCA", align = "left"),
               tags$div (
                 plotOutput("MutationTimeRplot", height = "250px")
               )
             )
        ),
      column(width = 6,
             div(
               h5("Evolutionary timeline (MutationTimeR)", align = "left"),
               tags$div (
                 plotOutput("timelineplotMutTimeR", height = "250px")
               )
             )
        )
    )

  })
  
  output$EvoChoices <- renderUI({
    if(input$EvoAna == "initiation"){
      radioButtons('MutDens', 'Mutation density method:', choices = "this paper")
    }
  }                
  )
  output$EvoPlots <- renderUI( if(input$EvoAna == "initiation"){
    fluidRow(
      uiOutput("MutDensPlots")
    )
   }else if(input$EvoAna == "growth"){
    fluidRow(
      column(width = 4,
             div(
               h5("Mobster fit", align = "left"),
               tags$div (
                # plotOutput("mobsterplot", height = "250px")
                 uiOutput("image_mobsterplot", height = "250px")
               )
             )
      ),
      column(width = 4,
             div(
               h5("Subclonal evolution model this paper", align = "left"),
               tags$div (
                 plotOutput("neutralevolplot", height = "250px")
               )
             )
      )
    )
  }          
  )
  
}

shinyApp(ui = ui, server = server)