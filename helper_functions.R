calc_pls <- function(data_1,data_2,n=2 ){
  
  #calculate PLS
  pls <- mixOmics::spls(data_1, data_2, ncomp = n, keepX = 5, keepY=5)
  
  #calc correlation between PC1
  x_axis_comp <- as_tibble(pls$variates$X, rownames = "sample_id")%>%
    rename_with(~str_c("x_", .x))%>%
    rename(sample_id = x_sample_id)
  y_axis_comp <- as_tibble(pls$variates$Y, rownames = "sample_id")%>%
    rename_with(~str_c("y_", .x))%>%
    rename(sample_id = y_sample_id)
  
  comp_table <- left_join(x_axis_comp,y_axis_comp, by = "sample_id")
  cor_test_results <- cor.test(comp_table$x_comp1, comp_table$y_comp1, method = "spearman")
  
  # write nice title to plot
  rho <- round(cor_test_results$estimat, digits = 3)
  p_value <- formatC(as.numeric(cor_test_results$p.value), digits=3)
  subtitle <- str_c("Spearman ??=",rho , ", p=", p_value)
  
  
  return(list(comp_table = comp_table,
              subtitle= subtitle,
              # pls = pls, 
              rho = rho,
              p_value = p_value))
}  

randomize_df <- function(df) {
  # set seed for reproducibility
  # set.seed(123)
  set.seed(NULL)
  
  # randomize values within each column
  df[] <- lapply(df, base::sample)
  
  # return randomized data frame
  return(df)
}

remove_rare_features <- function(feat_table, 
                                 filter_method = 'both', 
                                 prevalence_cutoff = 0.1, 
                                 avg_abundance_cutoff = 0.05) {
  # Argument verifications 
  if (prevalence_cutoff < 0 | prevalence_cutoff > 1) 
    error('Provided an invalid prevalence_cutoff value')
  if (avg_abundance_cutoff < 0) 
    error('Provided an invalid avg_abundance_cutoff value')
  
  # Initialize table with all current features
  new_feat_table <- feat_table
  n_samples <- nrow(feat_table)
  n_taxa_before_filter <- ncol(feat_table)-1
  
  # Prevalence calculations (number of non-zero values per feature)
  if (filter_method %in% c('prevalence', 'both')) {
    frequencies <- colSums(new_feat_table[,-1]>0) / n_samples
    new_feat_table <- new_feat_table[,c(TRUE, frequencies > prevalence_cutoff)]
  }
  
  # Average abundance calculations
  if (filter_method %in% c('avg_abundance', 'both')) {
    avg_abundances <- colSums(new_feat_table[,-1]) / n_samples
    new_feat_table <- new_feat_table[,c(TRUE, avg_abundances > avg_abundance_cutoff)]
  }
  
  n_taxa_after_filter <- ncol(new_feat_table)-1
  message(n_taxa_after_filter, ' of ',n_taxa_before_filter, ' taxa were removed.')
  
  return(new_feat_table)
}

normalize_rows <- function(df) {
  # loop through each row of the dataframe
  for (i in 1:nrow(df)) {
    # calculate the sum of values in the current row
    row_sum <- sum(df[i,])
    
    # divide each value in the current row by the row_sum
    df[i,] <- df[i,]/row_sum
  }
  
  # return the normalized dataframe
  return(df)
}