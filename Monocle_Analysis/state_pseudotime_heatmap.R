
#state_pseudotime_heatmap takes in the output of state_expression_matrix (seperated into object 1 and object 2)
#output gives a heatmap of the state-based gene expressions with state, cell type (cluster) and individual 
#annotations including gaps for major visual branches

library('pheatmap')

state_pseudotime_heatmap <- function(plot_matrix,anno_df){
  
  #set colors for heatmap
  cols <- colorRampPalette(c("blue","white","red"))(256)
  
  #set annotation colors that match State colors from ggplot2
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  #set colors for each annotation row individually
  state_cols <- gg_color_hue(length(levels(anno_df$State)))
  names(state_cols) <- levels(anno_df$State)
  
  clust_cols <- gg_color_hue(length(levels(anno_df$Cluster)))
  names(clust_cols) <- levels(anno_df$Cluster)
  
  ind_cols <- gg_color_hue(length(levels(anno_df$Individual)))
  names(ind_cols) <- levels(anno_df$Individual)
  
  anno_cols  = list(State = state_cols,  Cluster=clust_cols, Individual=ind_cols)
  
  #set breakpoints based on manual branch labels
  gap_locations <- cumsum(as.numeric(table(anno_df$Branch)))[-length(cumsum(as.numeric(table(anno_df$Branch))))]
  
  #remove manual branch data for plotting
  anno_df <- anno_df[,c("State","Cluster","Individual")]
  
  #create heatmap
  p_s_heatmap <- pheatmap(
    log(plot_matrix[,rownames(anno_df)]+1),
    scale="row",
    color= cols,
    cluster_cols=F,
    gaps_col = gap_locations,
    annotation_col=anno_df,
    annotation_colors=anno_cols)
  
  #return heatmap
  return(p_s_heatmap)
}
