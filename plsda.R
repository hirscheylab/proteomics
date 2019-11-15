library(mixOmics)
else if (method == "plsda"){
  
  X <- as.matrix(sapply(peptide_data_formatted[3: ncol(peptide_data_formatted)], as.numeric))
  Y <- peptide_data_formatted$Group #data_multi$Group
  
  plsda_res <- mixOmics::plsda(X, Y, ncomp = 4) #mixOmics::plsda(X, Y, ncomp = components)
  
  PLSDAi <- data.frame(plsda_res$variates$X, Groups = Y)
  
  scoresplot <- ggplot(PLSDAi, aes(x = comp.1, y = comp.2, col = Groups))+
    geom_point(size=3,alpha=0.5) +
    xlab("Component 1") +
    ylab("Component 2") +
    stat_ellipse(aes(x = comp.1, y = comp.2, col = Groups), type = "norm") +
    theme_minimal()
  
  print(scoresplot)
  #####
  
  set.seed(69)
  perf_plsda <- mixOmics::perf(plsda_res, validation = "Mfold", folds = 3,#mixOmics::perf(plsda_res, validation = validation, folds = folds,
                               progressBar = TRUE, auc = TRUE, nrepeat = 10)#progressBar = TRUE, auc = TRUE, nrepeat = nrepeat)
  
  overall <- data.frame(perf_plsda$error.rate[1]) %>%
    round(4) %>%
    rownames_to_column("Component")
  
  ber <- data.frame(perf_plsda$error.rate[2]) %>%
    round(4) %>%
    rownames_to_column("Component")
  
  errors_plsda1 <- reshape2::melt(ber, id.vars=c("Component"))
  errors_plsda2 <- reshape2::melt(overall, id.vars=c("Component"))
  errors_plsda <- rbind(errors_plsda1, errors_plsda2)
  
  errors_plsda_plot <- ggplot(data = errors_plsda, aes(x = Component, y = value, group = variable)) +
    geom_line(aes(color=variable)) +
    geom_point(aes(color=variable)) +
    theme_minimal() +
    geom_point(size=3,alpha=0.5)
  
  ####
  
  plsda_vip <- data.frame(mixOmics::vip(plsda_res))
  plsda_vip <- plsda_vip[order(plsda_vip[,1], decreasing = T) ,]
  
  plsda_vip <- plsda_vip %>% rownames_to_column("Variable")
  plsda_vip_top <- plsda_vip[1:15 ,]
  
  plsda_vip_top <- plsda_vip_top %>%
    mutate(Variable = factor(Variable, levels = Variable[order(comp.1)]))
  
  vip_plsda_plot <- ggplot(plsda_vip_top, aes(x = Variable, y = comp.1, fill = NULL)) +
    geom_bar(stat="identity", fill = rep(c("lightblue"), nrow(plsda_vip_top))) +
    coord_flip() +
    ylab("VIP") +
    geom_label(data = plsda_vip_top, aes(label = round(comp.1, 2))) +
    theme_minimal()
  print(vip_plsda_plot)
  ####
  
  scores_plsda <- PLSDAi %>% dplyr::select(-Groups) %>% round(4)
  
  return(list(scoresplot = scoresplot, errors_plsda = errors_plsda, errors_plsda_plot = errors_plsda_plot, plsda_vip_table = plsda_vip,
              vip_plsda_plot = vip_plsda_plot, scores_plsda = scores_plsda))
}