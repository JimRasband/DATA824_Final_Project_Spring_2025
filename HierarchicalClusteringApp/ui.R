library(shiny)

ui <- fluidPage(
  titlePanel("Differential Gene Expression by Ensembl ID :: Hierarchical Clustering Heatmap and KEGG Pathway Generator"),
  #note--i selected a max of 50 genes because I was concerned about crashing the app
  #i think this is a relatively safe amount were the app to use proper p-values and log2 fold change
  #as it stands, it might be more handy to increase that number.
  #call it future proofing as i continue to work on this. Same basic deal r/t number of clusters.  
  #Maybe more in future, but 10 was definitely plenty for a set of 50 genes.
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Choose Gene Expression Excel File", accept = ".xlsx"),
      numericInput("num_genes", "Number of Genes to Display (10–50)", value = 50, min = 10, max = 50),
      numericInput("k", "Number of Clusters", value = 3, min = 1, max = 10),
      checkboxInput("run_kegg", "Run KEGG Pathway Analysis", value = FALSE),
      
      conditionalPanel(
        condition = "input.run_kegg == true",
        radioButtons("kegg_output_type", "KEGG Output Type:",
                     choices = c("Table" = "table", "Barplot" = "barplot"),
                     selected = "table")
      ),
      
      actionButton("run_clustering", "Run Clustering")
    ),
    
    mainPanel(
      plotOutput("dendrogram_plot", height = "700px"),
      conditionalPanel(
        condition = "input.run_kegg == true",
        h3("KEGG Pathway Results"),
        uiOutput("kegg_output_ui")
      )
    )
  )
)
