state_expression_matrix <- function(expression_matrix, monocle_object, heatmap_genes, state_branch_data){
  #compute matrix from monocle_object
  average_matrix <- as.numeric()
  pseudotimes <- as.numeric()
  clusters <- as.character()
  individuals <- as.character()
  
  for(i in 1:dim(state_branch_data)[[1]]){
    
    #get data frame for each State i in the phenoData of our monocle_object
    temp <- subset(pData(monocle_object), State == as.character(state_branch_data[[i,1]]))
    
    #Compute the average pseudotime at State i
    pseudotimes <- c(pseudotimes, mean(temp$Pseudotime))
    
    #Establish the major cell type label for each State (which cell type is most common across cells)
    temp_clust <- table(temp$general_cluster)
    #print(which.max(temp_clust))
    clusters <- c(clusters, which.max(temp_clust))
    
    temp_ind <- table(temp$individual_label)
    #print(which.max(temp_ind))
    individuals <- c(individuals, which.max(temp_ind))
    
    #Compute average expression of heatmap genes for all cells in State i
    temp_mat <- expression_matrix[heatmap_genes,rownames(temp)]
    
    #If ncells > 1 in matrix, compute average, otherwise, take raw data
    if(dim(temp)[[1]] > 1){
      means_for_state <- rowMeans(temp_mat)
    }else{
      means_for_state <- temp_mat
    }	
    
    #append average expressions for each gene to average matrix
    average_matrix <- cbind(average_matrix, means_for_state)
  }
  
  #reformat clusters into a data frame object
  clusters <- names(clusters)
  clusters <- data.frame(clusters)
  rownames(clusters) <- state_branch_data$State
  
  individuals <- names(individuals)
  individuals <- data.frame(individuals)
  rownames(individuals) <- state_branch_data$State
  
  #reformat average_matrix into a data frame object
  colnames(average_matrix) <- state_branch_data[,1]
  plot_matrix <- as.matrix(average_matrix)
  
  #reformat pseudtimes into a data frame object
  pseudotimes <- data.frame(pseudotimes)
  rownames(pseudotimes) <- state_branch_data$State
  
  #combine initial dataset with new dataframe objects
  state_branch_data <- cbind(state_branch_data, pseudotimes,clusters, individuals)
  
  #set all metadata as factors
  state_branch_data$Branch <-factor(state_branch_data$Branch)
  state_branch_data$State <-factor(state_branch_data$State) 
  state_branch_data$clusters <- factor(state_branch_data$clusters)
  state_branch_data$individuals <- factor(state_branch_data$individuals)
  
  #generate final annotation dataframe and put in logical order for annotation rows
  anno_df <- as.data.frame(state_branch_data)
  rownames(anno_df) <- state_branch_data$State
  anno_df <- anno_df[order(anno_df$clusters),]
  anno_df <- anno_df[order(anno_df$pseudotimes),]
  anno_df <- anno_df[order(anno_df$Branch),]
  colnames(anno_df) <- c("State","Branch","Pseudotime","Cluster","Individual")
  
  return(list(plot_matrix,anno_df))
  
}	

