library(shiny)
library(readxl)
library(dplyr)
library(ggdendro)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(KEGGREST)
library(grid)
library(biomaRt)

server <- function(input, output, session) {
  
  # Takes ensembl codes and goes to biomart hsapiens dataset to grab ensembl ids
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Read and preprocess gene expression data
  gene_data <- reactive({
    req(input$file)
    #allows for upload of excel, reads uploads as excel
    df <- readxl::read_excel(input$file$datapath)
    #takes the first column (genename) and converts to a numeric vector
    #omits any na values (typically cells of 0 in the data I used--just gets rid of things which don't fit)
    df <- df %>%
      mutate(across(where(is.character), as.character)) %>%
      mutate(across(-GeneName, ~ as.numeric(as.character(.)))) %>%
      na.omit()
    
    # For each gene, calculates max absolute differential expression as how far from zero
    df$difference <- apply(dplyr::select(df, -GeneName), 1, function(x) max(abs(x), na.rm = TRUE))
    
    # Select most differentially expressed genes (i.e., furthest from 0 for expression values
    # selects the number of genes set by the user (e.g., the top 20, the top 50, etc.)
    top_n <- input$num_genes
    top_genes <- df %>%
      arrange(desc(difference)) %>%
      slice(1:top_n)
    
    # creates a dataframe using those top genes and their DE values
    mat <- as.data.frame(dplyr::select(top_genes, -GeneName, -difference))
    rownames(mat) <- top_genes$GeneName
    
    return(list(mat = mat, top_genes = top_genes))
  })
  
  # Convert Ensembl → Entrez → KEGG Pathways, store them for use in KEGG pathway analysis
  # Had there been more time, probably something that would be better to move up in the app
  # Because it would likely make it possible for either of Ensembl or Entrez IDs to be used in user uploaded data
  # And also personally I just prefer Entrez IDs as far as making things easy to understand.
  get_kegg_pathways_from_ensembl <- function(ensembl_ids) {
    mapping <- getBM(
      filters = "ensembl_gene_id",
      attributes = c("ensembl_gene_id", "entrezgene_id"),
      values = ensembl_ids,
      mart = ensembl
    )
    # In case there's not a match--mostly paranoia about typos in dataset
    entrez_ids <- na.omit(unique(mapping$entrezgene_id))
    if (length(entrez_ids) == 0) return(NULL)
    #stores those as KEGG id's
    kegg_ids <- paste0("hsa:", entrez_ids)
    #tryCAtch() is a very handy function for making mistakes like i did throughout
    #allows for "robust errors" such that the app doesn't crash if there's an error. 
    pathways <- tryCatch({
      keggLink("pathway", kegg_ids)
    }, error = function(e) {
      return(NULL)
    })
    
    return(pathways)
  }
  
  observeEvent(input$run_clustering, {
    req(gene_data())
    #this was a handy bit--gives an alert it's working.
    #with larger sets of genes, particularly from large files where it had to choose some genes
    #it often took a while--i think a user might think it's failing when it's just going slowly.
    withProgress(message = 'Generating heatmap...', value = 0.1, {
      data <- gene_data()$mat
      k <- input$k
      #ditto here to keep the user informed, basically, and not too impatient (me)
      incProgress(0.5, detail = "Clustering and plotting...")
      #Gives the clustered heatmap of DE genes by value, keeps participant ID at the bottom
      #such that it's easier to understand where data is coming from and think about participant characteristics
      #which i think is especially important given the goal was adjustable cluster and gene set sizes.
      output$dendrogram_plot <- renderPlot({
        ph <- pheatmap(data,
                       cluster_rows = TRUE,
                       cluster_cols = TRUE,
                       show_rownames = TRUE,
                       show_colnames = TRUE,
                       cutree_rows = k,
                       scale = "row",
                       color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
                       fontsize_row = 6,
                       main = paste("Top", input$num_genes, "Most Differentially Expressed Genes -", k, "Clusters"),
                       silent = TRUE)
        grid::grid.newpage()
        grid::grid.draw(ph$gtable)
      })
      
      # KEGG Analysis
      if (input$run_kegg) {
        top_genes <- gene_data()$top_genes
        ensembl_ids <- top_genes$GeneName
        kegg_paths <- get_kegg_pathways_from_ensembl(ensembl_ids)
        
        kegg_summary <- if (!is.null(kegg_paths)) {
          table(sub("path:", "", kegg_paths))
        } else {
          NULL
        }
        
        #Chatgpt assisted here--the || is interesting and maybe helpful
        #chatgpt alleges that it's safer r/t errors
        #it's basically functiong as an "or" type function so that if either condition is true, 
        #you get the message that there's no_kegg message, below "No KEGG pathways found."
        #in any case, apparently you can actually use a single line | as well, but probably shouldn't.
        output$kegg_output_ui <- renderUI({
          if (is.null(kegg_summary) || length(kegg_summary) == 0) {
            verbatimTextOutput("no_kegg")
          } else if (input$kegg_output_type == "table") {
            tableOutput("kegg_table")
          } else {
            plotOutput("kegg_barplot")
          }
        })
        
        #In the event KEGG pathways aren't present, so that the app doesn't break and so that users
        #can understand what's going on.  This bit prints output to say "none found."  Couple of safeguards
        #baked into the last few lines r/t possible missing data, too, with the idea of avoiding app crashing.
        output$no_kegg <- renderPrint({
          "No KEGG pathways found."
        })
        #creates the table--here I was a bit sloppy, but I ran out of time to make the table pretty
        # I would like for it to have included the participants who were associated with a particular pathway
        # as well as the ensembl ids.
        # that, or adjust such that entrez id is used consistently in all instances.
        #This would have made it much easier to draw conclusions, I think. 
        output$kegg_table <- renderTable({
          as.data.frame(kegg_summary, responseName = "Gene Count") %>%
            dplyr::rename(Pathway_ID = Var1)
        })
        #barplot turned out a bit nicer--not that it's saying much.  I set the limit to 15 pathways.
        output$kegg_barplot <- renderPlot({
          top_kegg <- sort(kegg_summary, decreasing = TRUE)[1:min(15, length(kegg_summary))]
          df <- as.data.frame(top_kegg)
          colnames(df) <- c("Pathway_ID", "Gene_Count")
          #this gives a bit more pizazz and a label for the barplot:
          ggplot(df, aes(x = reorder(Pathway_ID, -Gene_Count), y = Gene_Count)) +
            geom_bar(stat = "identity", fill = "#F1C338") +
            labs(title = "Top 15 KEGG Pathways", x = "KEGG Pathway ID", y = "Gene Count") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        })
      }
      #this lone warrior is for the speed with which the progress bars (e.g., r/t generating heatmaps)
      #move along.  It's really only relevant when you use larger data sets, I think, 
      #but it looks quite a lot nicer regardless (in my view).
      incProgress(.4)
    })  
  })    
}      
