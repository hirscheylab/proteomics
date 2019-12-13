# dat <- top_n(final_peptide_stats[[i]], 100, `Relative LOG2_FOLD_CHANGE`)
# dat <- subset(dat, `Relative LOG2_FOLD_CHANGE` > 0)
# dat <- subset(dat, `Relative P_VALUE` < 0.1)
# if(nrow(dat) > 0){
#   dat <- enrichr(genes = dat$GN, databases = "KEGG_2019_Human")$KEGG_2019_Human
#   dat <- subset(dat, Adjusted.P.value < 0.3)
#   
#   print(ggplot(data = dat, aes(x=reorder(Term, -Adjusted.P.value), y = Adjusted.P.value)) +
#           geom_col(aes(fill = Adjusted.P.value)) +
#           coord_flip() + 
#           scale_fill_gradient(low = "green4",
#                               high = "grey95") +
#           xlab("Pathway") +
#           ylab("Adjusted P-Value") +
#           theme_light())
# }
# 
# 
# dat <- top_n(final_peptide_stats[[i]], -100, `Relative LOG2_FOLD_CHANGE`)
# dat <- subset(dat, `Relative LOG2_FOLD_CHANGE` < 0)
# dat <- subset(dat, `Relative P_VALUE` < 0.1)
# if(nrow(dat) > 0){
#   dat <- enrichr(genes = dat$GN, databases = "KEGG_2019_Human")$KEGG_2019_Human
#   dat <- subset(dat, Adjusted.P.value < 0.3)
#   
#   dat <- enrichr(genes = c(), databases = "KEGG_2019_Human")$KEGG_2019_Human
#   print(ggplot(data = dat, aes(x=reorder(Term, -Adjusted.P.value), y = Adjusted.P.value)) +
#           geom_col(aes(fill = Adjusted.P.value)) +
#           coord_flip() + 
#           scale_fill_gradient(low = "green4",
#                               high = "grey95") +
#           xlab("Pathway") +
#           ylab("Adjusted P-Value") +
#           labs("Adjusted P-Value") +
#           theme_light())
# }

enrich <- function(data, topn, p_cutoff = 1, enrichr_p_cutoff = 1){
  data <- na.omit(data)
  names(data)[2] <- "LOG2_FOLD_CHANGE"
  names(data)[3] <- "P_VALUE"
  
  dat <- subset(data, `P_VALUE` < p_cutoff)
  dat <- top_n(dat, topn, `LOG2_FOLD_CHANGE`)
  dat <- subset(dat, `LOG2_FOLD_CHANGE` > 0)

  View(dat)
  
  if(nrow(dat) > 0){
    dat <- enrichr(genes = dat$GN, databases = "KEGG_2019_Human")$KEGG_2019_Human
    dat <- subset(dat, Adjusted.P.value < enrichr_p_cutoff)
    
    print(ggplot(data = dat, aes(x=reorder(Term, -Adjusted.P.value), y = Adjusted.P.value)) +
            geom_col(aes(fill = Adjusted.P.value)) +
            coord_flip() + 
            scale_fill_gradient(low = "green4",
                                high = "grey95") +
            xlab("Pathway") +
            ylab("Adjusted P-Value") +
            theme_light())
  }
}
