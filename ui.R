# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("scRNA-seq cluster visualization"),
  
  # Sidebar layout with input and output definitions ----
  fluidRow(
    
    # Sidebar panel for inputs ----
    column(9,   
         textInput(inputId = "genelist",
         label = "Genes to view", width="100%",
         #value = "SNAP25,ENO2,SLC17A7,SLC17A6,GAD1,GAD2,SLC32A1,LAMP5,SST,CHODL,PVALB,VIP,CUX2,RORB,RBP4,GJA1,FGFR3,GFAP,OLIG1,OPALIN,PDGFRA,AIF1,TYROBP,NOSTRIN")
         value = "Gene1,Gene2")
    ),
    column(2,radioButtons("modeval","Mode:",
                          c("Display" = "displayval",
                            "Find DE Genes"= "findgenesval"))
    ),
    column(1,
           actionButton("go", "Plot",style="color: #ffffff; background-color: #008800; border-color: #333333")
    )
  ),
  fluidRow(
    
    column(6,   
           textInput(inputId = "clustlist",
                     label = "Clusters to view (cluster numbers or 'all'/'neuron'/'glia')", width="100%",
                     value = "all")
    ),
    column(2,
           radioButtons("scaletyp", "Scale:",
                        c("Z-score" = "zscore",
                          "Log10(CPM+1)" = "logcpm",
                          "Z-score of Log10(CPM+1)" ="zscorelogcpm"))
    )
  ),
  fluidRow(
    column(5,   
           textInput(inputId = "clustonlist",
                     label = "'On' clusters (numbers or 'neuron'/'glia')", width="100%",
                     value = "1,2,3")
    ),
    column(5,   
           textInput(inputId = "clustofflist",
                     label = "'Off' clusters (numbers or 'all'/'neuron'/'glia')", width="100%",
                     value = "all")
    ),
    column(2,   
           textInput(inputId = "numgeneval",
                     label = "Number of genes to show", width="100%",
                     value = "15")
    )
  ),
  fluidRow(
      column(12,
        h3(verbatimTextOutput("caption"))
      )),
  fluidRow(
    column(12,
           h3(verbatimTextOutput("caption2"))
    )
  ),
  fluidRow(
      column(12,
        plotOutput(outputId = "dotPlot",height="800px")
      )
  )
)
