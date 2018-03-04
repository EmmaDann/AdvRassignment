### SHINY SERVER ###
library(shiny)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(tcgatools)

#### Define UI ----
ui <- fluidPage(
  titlePanel("EXPLORE TCGA"),
  sidebarLayout(
    sidebarPanel(
      # Select TCGA data file
      fileInput("file", h6("Choose file")),
      #--- Horizontal line
      tags$hr(),
      # Select expression data of interest
      radioButtons("XlevelType",
                   label = "Expression level type",
                   choices = list("RNA expression"="RNA",
                                  "RNA expression (z-transformed)"="RNAZ",
                                  "Protein expression"="RPPA")
                   )
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Expression level",
                           plotOutput("XpressBoxplots"),
                           checkboxInput("normalize",
                                        label = "Show normalized data",
                                        value=FALSE),
                           dataTableOutput("datatable")
                           ),
                  tabPanel("Copy number alterations",
                           plotOutput("scatterplots"),
                           dataTableOutput("cnaTable")
                           ),
                  tabPanel("Metadata",
                           textInput("selectedMetadataColumn", "Select clinical data", placeholder = "AGE,OS_MONTHS,RACE..."),
                           dataTableOutput("metadatatable"))
        
      )
    )
  )
)

### Define server logic ----
server <- function(input, output, session) {
  # --- REACTIVE FUNCTIONS: LOAD AND RESHAPE DATA
  load_data <- reactive({ # Load data filling empty spaces with NAs
    req(input$file)
    df <- read.csv(input$file$datapath, na.strings=c('','NA')) 
    return(df)
  })
  
  normalize_data <- reactive({ 
    tcga.data <- load_data()
    if (input$normalize) {
      tcga.data = norm.columns(tcga.data, cols.to.normalize = extract.cols.ix(tcga.data, input$XlevelType))
    } 
    return(tcga.data)
  })
  
  reshape4xpressBoxplot <- reactive({
    tcga.table <- normalize_data()
    # Extract columns of expression data of interest
    selected.data <- extract.subtype(tcga.table, subtype=input$XlevelType, with_metadata = FALSE)
    # Make long df with information about sample, gene and expression level (to plot)
    long.tcga.table = melt(selected.data, variable.name = "gene", value.name = "xpress.level") %>%
      mutate(gene=sapply(gene, function(x) gsub(pattern = '.+\\.', replacement = '',x))) # Remove subtype name
    return(long.tcga.table)
  })
  
  selectCNAdata <- reactive({
    tcga.table <- normalize_data()
    # Extract columns of CNA and expression data of interest
    cna.data <- cbind(
      extract.subtype(tcga.table, paste0('CNA|', input$XlevelType), with_metadata = FALSE),
      OS_MONTHS=tcga.table$OS_MONTHS
    )
    return(cna.data)
  })
  
  reshape4CNAplot <- reactive({
    # Reshaping dataframe to yield information on expression and CNA for one gene in one sample (to plot)
    long_table <- melt(selectCNAdata(), id.vars = c('sample', 'OS_MONTHS')) %>% 
      mutate(subtype = sapply(variable, function(var) strsplit(x = as.character(var),split= '\\.')[[1]][1]), 
             gene = sapply(variable, function(var) strsplit(x = as.character(var),split= '\\.')[[1]][2] ))
    wide_table <- dcast(long_table, formula = sample + gene + OS_MONTHS ~ subtype, value.var = 'value')
    colnames(wide_table) <- c('sample', 'gene', 'OS_MONTHS', 'CNA', 'Expression')  
    return(wide_table)
  })    

  load_metadata <- reactive({
    metadata <- load_data()[- extract.cols.ix(load_data(),prefix = 'RNA|CNA|MUT|RNAZ|PPA')]
    return(metadata)
  })
  
  # --- TAB 1: EXPRESSION DATA --- 
  output$XpressBoxplots <- renderPlot({
    data <- reshape4xpressBoxplot()
    ggplot(data, aes(gene, xpress.level)) +
      geom_boxplot(varwidth = TRUE, fill='lightblue') +
      ylab('Expression level') 
  })

  output$datatable <- renderDataTable(
    # Display table of selected expression data
    extract.subtype(normalize_data(), subtype=input$XlevelType, with_metadata = FALSE)
  )
 
  # --- TAB 2: CNA DATA --- 

  output$scatterplots <- renderPlot({
      # Plot relationship between copy number and expression for available genes
      data <- reshape4CNAplot()
      ggplot(data, aes(CNA, Expression, color=OS_MONTHS)) +
        facet_wrap( ~ gene , scales='free') +
        geom_point() +
        # Using scale_color_gradientn to obtain high resolution on lower values of OS_MONTHS
        scale_color_gradientn(colors=brewer.pal(n = 9, name = 'YlOrRd'), 
                              values=c(1.0,0.75,0.5,0.45,0.4,0.35,0.3,0.25, 0.2,0.15,0.1,0),
                              name="Months of survival") +
        geom_vline(xintercept=0, color='red') 
    })
  
  output$cnaTable <- renderDataTable(
    # Display table of selected CNA data
    selectCNAdata()
  )

  # --- TAB 3: METADATA --- 
  output$metadatatable <- renderDataTable(
     # Display table of metadata
     if (input$selectedMetadataColumn!='') {
      load_metadata()[,c('sample',input$selectedMetadataColumn)]
     } else {
      load_metadata()
     }
  )
}

### Run the app ----
shinyApp(ui = ui, server = server)
