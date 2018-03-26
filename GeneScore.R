GeneScore<-function(matrix,idents,genes,plot.title){
  set.seed(1)
  #Matrix
  ind1_trimmed_log<-matrix
  print(dim(ind1_trimmed_log))
  #Cell IDs
  cell.id.trimmed<-idents
  gene_agg<-rowSums(ind1_trimmed_log)
  head(gene_agg)
  gene_agg_sort<-sort(gene_agg)
  range(gene_agg_sort) #0.7280198 2152.8409068
  gene_agg_sort_binned<-cut(gene_agg_sort,25,labels=FALSE)
  #gene_agg_sort_binned
  #hist(gene_agg_sort_binned)
  Bin1_genes<-which(gene_agg_sort_binned == "1")
  Bin2_genes<-which(gene_agg_sort_binned == "2")
  Bin3_genes<-which(gene_agg_sort_binned == "3")
  Bin4_genes<-which(gene_agg_sort_binned == "4")
  Bin5_genes<-which(gene_agg_sort_binned == "5")
  Bin6_genes<-which(gene_agg_sort_binned == "6")
  Bin7_genes<-which(gene_agg_sort_binned == "7")
  Bin8_genes<-which(gene_agg_sort_binned == "8")
  Bin9_genes<-which(gene_agg_sort_binned == "9")
  Bin10_genes<-which(gene_agg_sort_binned == "10")
  Bin11_genes<-which(gene_agg_sort_binned == "11")
  Bin12_genes<-which(gene_agg_sort_binned == "12")
  Bin13_genes<-which(gene_agg_sort_binned == "13")
  Bin14_genes<-which(gene_agg_sort_binned == "14")
  Bin15_genes<-which(gene_agg_sort_binned == "15")
  Bin16_genes<-which(gene_agg_sort_binned == "16")
  Bin17_genes<-which(gene_agg_sort_binned == "17")
  Bin18_genes<-which(gene_agg_sort_binned == "18")
  Bin19_genes<-which(gene_agg_sort_binned == "19")
  Bin20_genes<-which(gene_agg_sort_binned == "20")
  Bin21_genes<-which(gene_agg_sort_binned == "21")
  Bin22_genes<-which(gene_agg_sort_binned == "22")
  Bin23_genes<-which(gene_agg_sort_binned == "23")
  Bin24_genes<-which(gene_agg_sort_binned == "24")
  Bin25_genes<-which(gene_agg_sort_binned == "25")
  #head(gene_agg_sort)
  
  tmp<-rbind(gene_agg_sort,rep(0,length(gene_agg_sort)))
  tmp<-as.data.frame(tmp)
  tmp[2,Bin1_genes]<-"Bin1"
  tmp[2,Bin2_genes]<-"Bin2"
  tmp[2,Bin3_genes]<-"Bin3"
  tmp[2,Bin4_genes]<-"Bin4"
  tmp[2,Bin5_genes]<-"Bin5"
  tmp[2,Bin6_genes]<-"Bin6"
  tmp[2,Bin7_genes]<-"Bin7"
  tmp[2,Bin8_genes]<-"Bin8"
  tmp[2,Bin9_genes]<-"Bin9"
  tmp[2,Bin10_genes]<-"Bin10"
  tmp[2,Bin11_genes]<-"Bin11"
  tmp[2,Bin12_genes]<-"Bin12"
  tmp[2,Bin13_genes]<-"Bin13"
  tmp[2,Bin14_genes]<-"Bin14"
  tmp[2,Bin15_genes]<-"Bin15"
  tmp[2,Bin16_genes]<-"Bin16"
  tmp[2,Bin17_genes]<-"Bin17"
  tmp[2,Bin18_genes]<-"Bin18"
  tmp[2,Bin19_genes]<-"Bin19"
  tmp[2,Bin20_genes]<-"Bin20"
  tmp[2,Bin21_genes]<-"Bin21"
  tmp[2,Bin22_genes]<-"Bin22"
  tmp[2,Bin23_genes]<-"Bin23"
  tmp[2,Bin24_genes]<-"Bin24"
  tmp[2,Bin25_genes]<-"Bin25"
  
  of_interest_genes<-genes
  
  #Run Scores####
  #of_interest_genes
  of_interest_bin_ID<-colnames(tmp) %in% of_interest_genes[,1] #Identify colnames of genes in tmp
  of_interest_bin_ID<-grep("TRUE",of_interest_bin_ID) #Identify column numbers of genes in tmp
  of_interest_binned<-tmp[,of_interest_bin_ID] #Trim binned matrix to only cell cycle genes
  of_interest_sampling<-(of_interest_binned[2,]) #Determine bins that need to be sampled from
  
  #Count number of genes in each bin
  counts<-list()
  for(i in 1:25)
  {
    counts[[i]]<-length(grep(paste("Bin",i,"$",sep=''),of_interest_sampling[1,]))
  }
  counts<-unlist(counts)
  #Loop through and make a list of each of the genes for each sample, with each index being the list of genes for sample
  gene_samples<-list()
  k=1
  bin_names<-list()
  for(i in 1:25)
  {
    bin_names[i]<-unlist(noquote(paste("Bin",i,"_genes",sep='')))
  }
  
  for(i in 1:25)
  {
    if(counts[i]>0)
    {
      for(j in 1:counts[i])
      {
        if (length(get(bin_names[[i]]))>100)
        {
          gene_samples[[k]]<-sample(colnames(tmp)[get(bin_names[[i]])], 100)
          #print(paste("Bin",i,"_genes",sep=''))
        }
        else
        {
          gene_samples[[k]]<-sample(colnames(tmp)[get(bin_names[[i]])], 100,replace=TRUE)
        }
        
        k = k + 1
      }
    }
    else
    {}
  }
  
  #Loop through length of cell cycle sampling and create sample of genes from its corresponding bin
  #of_interest_sampling
  
  
  of_interest_inds<-match(colnames(of_interest_sampling),rownames(ind1_trimmed_log))
  of_interest_matrix<-ind1_trimmed_log[of_interest_inds,]
  of_interest_gene_scores<-colMeans(of_interest_matrix)
  
  #Loop through all cell cycle sample numbers and rbind together forever
  of_interest_control_matrix<-
    ind1_trimmed_log[gene_samples[[1]],]
  #rownames(of_interest_control_matrix)<-1:100
  for(i in 2:length(gene_samples))
  {
    of_interest_control_matrix<-rbind(of_interest_control_matrix,ind1_trimmed_log[gene_samples[[i]],],make.row.names=TRUE)
  }
  
  of_interest_control_scores<-colMeans(of_interest_control_matrix)
  of_interest_scores<-of_interest_gene_scores-of_interest_control_scores
  of_interest_score_matrix<-as.matrix(of_interest_scores)
  trimmed_cells_clusterID<-cell.id.trimmed[,1] %in% rownames(of_interest_score_matrix)
  trimmed_cells_clusterID_IDs<-grep("TRUE",trimmed_cells_clusterID) #Identify column numbers of genes in tmp
  cell_id_trimmed_trimmed<-cell.id.trimmed[trimmed_cells_clusterID_IDs,]
  of_interest_score_plotting_matrix<-cbind.data.frame(of_interest_score_matrix, cell_id_trimmed_trimmed[,2])
  colnames(of_interest_score_plotting_matrix)<-c("Score","ClusterID")
  #Plotting  Scores CHANGE TITLE FOR DIFFERENT GRAPHS####
  q<-ggplot(of_interest_score_plotting_matrix,aes(x=ClusterID, y=Score,color=ClusterID))
  q+
    geom_violin(aes(fill=ClusterID))+
    geom_jitter(height=0,color="black",size=0.1)+
    ggtitle(plot.title)+
    theme_bw()
  
}
