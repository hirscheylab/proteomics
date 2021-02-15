

protein <- function(protein_data, peptide_data){
  protein_data <- na.omit(protein_data)
  names(protein_data)[2] <- "LOG2_FOLD_CHANGE"
  names(protein_data)[3] <- "P_VALUE"
  
  peptide_data <- na.omit(peptide_data)
  names(peptide_data)[2] <- "LOG2_FOLD_CHANGE"
  names(peptide_data)[3] <- "P_VALUE"
  
  dat <- left_join(protein_data, select(peptide_data, Accession, LOG2_FOLD_CHANGE), by = "Accession")
  dat <- melt(dat, id = "Accession", measure = c("LOG2_FOLD_CHANGE.y", "LOG2_FOLD_CHANGE.x"))
  dat <- left_join(dat, select(protein_data, Accession, LOG2_FOLD_CHANGE), by = "Accession")
  
  print(ggplot(data = dat, aes(x = reorder(Accession, LOG2_FOLD_CHANGE), y = value, color = variable)) +
          geom_point(alpha = 0.3) + 
          scale_color_manual(values = c("red2", "grey40")) +
          theme(axis.ticks.x = element_blank(),
                axis.text.x = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(color = "grey"),
                panel.background = element_rect(fill = "white", color = "grey")) + 
          geom_text_repel(data = top_n(dat, 20, abs(value)), 
                          aes(label=`Accession`),
                          segment.size = 0.2,
                          segment.color = "grey50",
                          point.padding = 1,
                          color = "black") +
          labs(color = "Modifications and Relative Adbundance") +
          xlab("Proteins") +
          ylab("Relative" ~Log[2]~ "Fold Change"))
}
