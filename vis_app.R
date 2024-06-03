library(shiny)
library(shinyjs)
library(shinythemes)
library(plotly)
library(stringr)
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



# default start ID
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
      combProps <- summary(as.factor(dmRes$Subclone)) *100 / nrow(dmRes)
    } else {
      combProps <- summary(as.factor(atacRes$Subclone)) *100 / nrow(atacRes)
    }
    subcloneInfo <- c()
    for (i in 1:(length(combProps)-1)) {
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
  
  
}

shinyApp(ui = ui, server = server)