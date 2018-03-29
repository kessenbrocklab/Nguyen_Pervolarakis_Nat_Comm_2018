library('monocle')
library('mclust')

#all data curated contains a matrix consisting for all 4 droplet based dataset, 
#trimmed for genes kept in the Seurat analysis
#each of ind4_data, ind5_data, ind6_data, and ind_7 data give the names of 
#the cells kept from each individaul in the Seurat analysis

########################Ind4 Only Ordering##############################
	
	#individual 4 only matrix (gene and cell trimmed based on Seurat analysis)
	expression.data1 <- all_data_curated[,ind4_data]
	
	pheno.data1 <- colnames(expression.data1)
	pheno.data.df1 <- data.frame(pheno.data1)
	rownames(pheno.data.df1) <- colnames(expression.data1) 
	pd1 <- new('AnnotatedDataFrame', data = pheno.data.df1)

	feature.data1 <- rownames(expression.data1)
	feature.data.df1 <- data.frame(feature.data1)
	rownames(feature.data.df1) <- rownames(expression.data1)
	fd1 <- new('AnnotatedDataFrame', data = feature.data.df1) frame

	I4_only <- newCellDataSet(expression.data1, phenoData = pd1, featureData = fd1, lowerDetectionLimit=0.5, expressionFamily = negbinomial.size())
	colnames(fData(I4_only)) <- "gene_short_name"
	
	#extract Seurat cluster designations
	pData(I4_only)$ClusterID <- factor(pData(I4_only)$ClusterID)
	full_cluster_subset <- full_cluster[rownames(pData(I4_only)),]
	general_cluster <- full_cluster_subset$ClusterID
	pData(I4_only) <- cbind(pData(I4_only),general_cluster)
	
	I4_only <- estimateSizeFactors(I4_only)
	I4_only <- estimateDispersions(I4_only)
	I4_only <- detectGenes(I4_only, min_expr = 0.1) 

	I4_only <- setOrderingFilter(I4_only, marker_genes_i4)
	I4_only <- reduceDimension(I4_only, max_components = 2) 
	I4_only <- orderCells(I4_only, reverse = F)

########################Ind5 Only Ordering##############################
	
	#individual 5 only matrix (gene and cell trimmed based on Seurat analysis)
	expression.data1 <- all_data_curated[,ind5_data]
	
	pheno.data1 <- colnames(expression.data1)
	pheno.data.df1 <- data.frame(pheno.data1)
	rownames(pheno.data.df1) <- colnames(expression.data1)
	pd1 <- new('AnnotatedDataFrame', data = pheno.data.df1)

	feature.data1 <- rownames(expression.data1)
	feature.data.df1 <- data.frame(feature.data1)
	rownames(feature.data.df1) <- rownames(expression.data1)
	fd1 <- new('AnnotatedDataFrame', data = feature.data.df1)

	I5_only <- newCellDataSet(expression.data1, phenoData = pd1, featureData = fd1, lowerDetectionLimit=0.5, expressionFamily = negbinomial.size()) #generate CellDataSet used by Monocle
	colnames(fData(I5_only)) <- "gene_short_name"
	
	#extract Seurat cluster designations
	pData(I5_only)$ClusterID <- factor(pData(I5_only)$ClusterID)
	full_cluster_subset <- full_cluster[rownames(pData(I5_only)),]
	general_cluster <- full_cluster_subset$ClusterID
	pData(I5_only) <- cbind(pData(I5_only),general_cluster)
	
	I5_only <- estimateSizeFactors(I5_only)
	I5_only <- estimateDispersions(I5_only)
	I5_only <- detectGenes(I5_only, min_expr = 0.1)

	I5_only <- setOrderingFilter(I5_only, marker_genes_i5)
	I5_only <- reduceDimension(I5_only, max_components = 2) 
	I5_only <- orderCells(I5_only, reverse = F)
	
########################Ind6 Only Ordering##############################
	
	#individual 6 only matrix (gene and cell trimmed based on Seurat analysis)
	expression.data1 <- all_data_curated[,ind6_data]
	
	pheno.data1 <- colnames(expression.data1)
	pheno.data.df1 <- data.frame(pheno.data1)
	rownames(pheno.data.df1) <- colnames(expression.data1)
	pd1 <- new('AnnotatedDataFrame', data = pheno.data.df1)

	feature.data1 <- rownames(expression.data1)
	feature.data.df1 <- data.frame(feature.data1)
	rownames(feature.data.df1) <- rownames(expression.data1)
	fd1 <- new('AnnotatedDataFrame', data = feature.data.df1)

	I6_only <- newCellDataSet(expression.data1, phenoData = pd1, featureData = fd1, lowerDetectionLimit=0.5, expressionFamily = negbinomial.size())
	colnames(fData(I6_only)) <- "gene_short_name"
	
	#extract Seurat cluster designations
	pData(I6_only)$ClusterID <- factor(pData(I6_only)$ClusterID)
	full_cluster_subset <- full_cluster[rownames(pData(I6_only)),]
	general_cluster <- full_cluster_subset$ClusterID
	pData(I6_only) <- cbind(pData(I6_only),general_cluster)
	
	I6_only <- estimateSizeFactors(I6_only)
	I6_only <- estimateDispersions(I6_only)
	I6_only <- detectGenes(I6_only, min_expr = 0.1)

	I6_only <- setOrderingFilter(I6_only, marker_genes_i6)
	I6_only <- reduceDimension(I6_only, max_components = 2) 
	I6_only <- orderCells(I6_only, reverse = F)
	
########################Ind7 Only Ordering##############################
	
	#individual 7 only matrix (gene and cell trimmed based on Seurat analysis)
	expression.data1 <- all_data_curated[,ind7_data]
	
	pheno.data1 <- colnames(expression.data1)
	pheno.data.df1 <- data.frame(pheno.data1)
	rownames(pheno.data.df1) <- colnames(expression.data1)
	pd1 <- new('AnnotatedDataFrame', data = pheno.data.df1)

	feature.data1 <- rownames(expression.data1)
	feature.data.df1 <- data.frame(feature.data1)
	rownames(feature.data.df1) <- rownames(expression.data1)
	fd1 <- new('AnnotatedDataFrame', data = feature.data.df1)

	I7_only <- newCellDataSet(expression.data1, phenoData = pd1, featureData = fd1, lowerDetectionLimit=0.5, expressionFamily = negbinomial.size()) #generate CellDataSet used by Monocle
	colnames(fData(I7_only)) <- "gene_short_name"
	
	#extract Seurat cluster designations
	pData(I7_only)$ClusterID <- factor(pData(I7_only)$ClusterID)
	full_cluster_subset <- full_cluster[rownames(pData(I7_only)),]
	general_cluster <- full_cluster_subset$ClusterID
	pData(I7_only) <- cbind(pData(I7_only),general_cluster)
	
	I7_only <- estimateSizeFactors(I7_only)
	I7_only <- estimateDispersions(I7_only)
	I7_only <- detectGenes(I7_only, min_expr = 0.1)

	I7_only <- setOrderingFilter(I7_only, marker_genes_i7)
	I7_only <- reduceDimension(I7_only, max_components = 2) 
	I7_only <- orderCells(I7_only, reverse = F)

########################Correlation Analysis Begins##############################
		
  	#Differential expression test on each individual ordering
	diff_test_res.I4 <- differentialGeneTest(I4_only, fullModelFormulaStr = "~State")
	diff_test_res.I5 <- differentialGeneTest(I5_only, fullModelFormulaStr = "~State")
	diff_test_res.I6 <- differentialGeneTest(I6_only, fullModelFormulaStr = "~State")
	diff_test_res.I7 <- differentialGeneTest(I7_only, fullModelFormulaStr = "~State")

	#FDR < 0.01
	i5_sig <- subset(diff_test_res.I5, qval < 0.01)
	i6_sig <- subset(diff_test_res.I6, qval < 0.01)
	i7_sig <- subset(diff_test_res.I7, qval < 0.01)
	i4_sig <- subset(diff_test_res.I4, qval < 0.01)
	
	#Find genes that had FDR < 0.01 across all individuals
	shared_genes <- Reduce(intersect, list(rownames(i4_sig),rownames(i5_sig),rownames(i6_sig),rownames(i7_sig)))
	
	#Generate average state expression matrix for shared_genes 
	i4_mat_data <- state_expression_matrix(all_data_curated,I4_only,shared_genes,state_branch_i4)
	i5_mat_data <- state_expression_matrix(all_data_curated,I5_only,shared_genes,state_branch_i5)
	i6_mat_data <- state_expression_matrix(all_data_curated,I6_only,shared_genes,state_branch_i6)
	i7_mat_data <- state_expression_matrix(all_data_curated,I7_only,shared_genes,state_branch_i7)
	
	#Compute gene correlation matrix for each individual
	i4_corr <- cor(t(i4_mat_data[[1]]))
	i5_corr <- cor(t(i5_mat_data[[1]]))
	i6_corr <- cor(t(i6_mat_data[[1]]))
	i7_corr <- cor(t(i7_mat_data[[1]]))
	
	#Compute average correlation across individual correlation matricies
	my.list <- list(i4_corr,i5_corr, i6_corr, i7_corr)
	averages <- Reduce("+", my.list) / length(my.list)

	#Take all genes that had a correlation > 0.8 with at least one other gene (genes that have module potential)
	correlation_sets_avg <- as.character()
	for(i in 1:dim(averages)[1]){
		correlation_sets_avg <- c(correlation_sets_avg,names(averages[i,averages[i,] > 0.8]))
	}
	
	#Reduce average correlation matrix to only genes that have potential to be in modules
	unique_genes_of_set_avg <- unique(correlation_sets_avg)
	reduced_correlation_matrix_averages <- averages[unique_genes_of_set_avg,unique_genes_of_set_avg]
	
	#Cluster on averaged correlation matrix
	corr_clust_avg <- Mclust(reduced_correlation_matrix_averages)
	cluster_labels_avg <- data.frame(corr_clust_avg$classification[order(corr_clust_avg$classification)])
	
	#Write out mclust (EM algorithm) gene/cluster labels to a table
	write.table(cluster_labels_avg, "~/Desktop/Correlation_Sets_Avg.txt")
