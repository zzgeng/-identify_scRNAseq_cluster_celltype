##Using the most variable genes from each clusters to identify the potential cell identity for each clusters

##Using chi-square test

##Cell type-specific markers list is based on the data from https://panglaodb.se/markers.html?cell_type=%27choose%27 

## contingency table
##################    DEG  | not DEG
#with GO         |    A    |   B
#without GO      |    C    |   D 


library(dplyr)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)

gene_list <- as.data.frame(readxl::read_xlsx("Desktop/HBO/PanglaoDB_markers_27_Mar_2020.xlsx"))
cell_type <- unique(gene_list$`cell type`)

Marker_cluster <- read.delim("Desktop/HBO/cluster_marker.txt")
cluster <- unique(Marker_cluster$cluster)


results <- matrix(nrow = length(cell_type))
rownames(results) <-  cell_type

for(x in cluster){
  a <- c()
  for(y in 1:length(cell_type)){
  # Get the matrix for chi-square test
  ##up & GO
  N11 <- sum(Marker_cluster[Marker_cluster$cluster == x, 7] %in% gene_list[gene_list$`cell type` == cell_type[y],2])
  
  ##GO !up
  N21 <-  nrow(gene_list[gene_list$`cell type` == cell_type[y],]) - N11
  
  ##!GO uo
  N12 <- nrow(Marker_cluster[Marker_cluster$cluster == x, ]) - N11
  
  N22 <- 24000 - N21
  
  chi_matrix <- matrix(c(N11, N21, N12, N22),2,2)
  
  #pvalue
  a[y] <- chisq.test(chi_matrix)$p.value
  }
results <- cbind(results,a)
}

results <- results[,-1]
colnames(results) <- paste0("Cluster_", c(0:9))

Marker_plot <-  function(input, top_n = 5){
  output <- matrix(ncol = ncol(input))
  name_ROW <- c()               
  for(x in 1:ncol(input)){
    A <- input[order(input[,x], decreasing = F), ]
    A <- -log10(A[1:top_n,])
    name_ROW <- c(name_ROW, rownames(A))
    output <- rbind(output, A)
  }
  output <- output[-1,]
  rownames(output) <- c(name_ROW)
  return(output)
  
}

marker_results <- Marker_plot(results, top_n = 5)


Heatmap(t(scale(t(marker_results))), cluster_rows = F, cluster_columns = F, column_split = paste0("C", 0:(ncol(marker_results)-1)), row_split = rep(paste0("C", 0:(ncol(marker_results)-1)), each = 5))