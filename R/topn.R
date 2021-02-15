

topn <- function(data, topn){
  data <- na.omit(data)
  names(data)[2] <- "REL_LOG2_FOLD_CHANGE"
  names(data)[3] <- "REL_P_VALUE"
 
  dat <- bind_rows(top_n(data, topn, REL_LOG2_FOLD_CHANGE), top_n(data, -topn, REL_LOG2_FOLD_CHANGE))
  dat$id <- 1:nrow(dat)
  dat$sigstars <- NA
  for(j in 1:nrow(dat)){
    if(dat$REL_P_VALUE[j] < 0.001){
      dat$sigstars[j] = "***"
    }
    if(0.01 > dat$REL_P_VALUE[j] & dat$REL_P_VALUE[j] > 0.01){
      dat$sigstars[j] = "**"
    }
    if(0.1 > dat$REL_P_VALUE[j] & dat$REL_P_VALUE[j] > 0.05){
      dat$sigstars[j] = "*"
    }
    if(dat$REL_P_VALUE[j] > 0.1){
      dat$sigstars[j] = ""
    }
  }
  dat <- arrange(dat, REL_LOG2_FOLD_CHANGE)
  
  topn <- ggplot(data=dat,
               aes(x=reorder(factor(id), REL_LOG2_FOLD_CHANGE),
                   y=REL_LOG2_FOLD_CHANGE)) + 
          geom_bar(aes(fill = REL_LOG2_FOLD_CHANGE < 0),
                   stat = "identity",
                   width = 0.5) + 
          scale_fill_manual(guide = FALSE,
                            breaks = c(TRUE, FALSE),
                            values=c("green4","red3")) +
          scale_x_discrete(labels=dat$GN) +
          coord_flip() +
          geom_text(data = dat,
                    aes(label = sigstars),
                    nudge_y = (ifelse(dat$REL_LOG2_FOLD_CHANGE<0, (max(abs(dat$REL_LOG2_FOLD_CHANGE))*-0.04), (max(abs(dat$REL_LOG2_FOLD_CHANGE))*0.04)))) +
          ylab("Relative" ~Log[2]~ "Fold Change") + 
          xlab("Gene ID") +
          theme_light()
  
  return(topn)
}
