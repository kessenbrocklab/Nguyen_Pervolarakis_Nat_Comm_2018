# Nguyen_Pervolarakis_Nat_Comm_2018
Supplementary software for "Profiling human breast epithelial cells using single cell RNA sequencing identifies cell diversity"

Table of Contents
=================
* [Monocle_Analysis ](https://github.com/kessenbrocklab/Nguyen_Pervolarakis_Nat_Comm_2018/edit/master/Monocle_Analysis)
    * [Correlation_Analysis_for_Combined_10x_Ordering.R](https://github.com/kessenbrocklab/Nguyen_Pervolarakis_Nat_Comm_2018/edit/master/Monocle_Analysis/Correlation_Analysis_for_Combined_10x_Ordering.R): Analysis used to choose Droplet based ordering genes from Supplementary Data 2
    * [state_expression_matrix.R](https://github.com/kessenbrocklab/Nguyen_Pervolarakis_Nat_Comm_2018/edit/master/Monocle_Analysis/state_expression_matrix.R): Custom function built to generate a gene expression matrix by State (used in the combined Droplet based Monocle ordering - Figure 6)
      * [Combined_4000_State_Branch.txt](https://github.com/kessenbrocklab/Nguyen_Pervolarakis_Nat_Comm_2018/edit/master/Monocle_Analysis/Combined_4000_State_Branch.txt): example of the "state_branch_data" data frame structure required for input into state_expression_matrix.R
      * [state_expression_heatamp.R](https://github.com/kessenbrocklab/Nguyen_Pervolarakis_Nat_Comm_2018/edit/master/Monocle_Analysis/state_expression_heatmap.R): Custom function used to generate a formatted heatmap of the expression matrix and annotation data frame output by state_expression_matrix.R (not directly used in manuscript)
 


