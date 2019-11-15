 

pool_normalization <- function(data, num_sample) {
    head <- select(data, (1:(ncol(data) - num_sample)))
    data <- bind_cols(select(data, "unique_id"),
                      select(data, tail(1:ncol(data), num_sample)))
    pool_index <- grep("Group = NA", names(data))
    pool_avg <- rowMeans(data[pool_index])
    
    pool_index <- grep("Group = NA", names(data))
    pool_avg <- rowMeans(data[pool_index])
     <- select(data, grep("LOG2 - AVG", names(data))) - data$Pool_AVG
    names(LOG2subAVG_poolnorm) <- gsub("LOG2 - AVG", "Pool Normalized LOG2-AVG", names(LOG2subAVG_poolnorm))
    data <- bind_cols(data, LOG2subAVG_poolnorm)
    
    # Initialize data structures and calculate the Log2 abundances
    head <- select(data, (1:(ncol(data) - num_sample)))
    data <- bind_cols(select(data, "unique_id"),
                      select(data, tail(1:ncol(data), num_sample)))
    data_log2 <- log2(data[2:ncol(data)])
    names(data_log2) <- paste("LOG2 ", names(data_log2))
    
    # Calculate the average of the Log2 abundances
    data_log2$LOG2_AVG <- rowMeans(select(data_log2, grep("LOG2 ", names(data_log2))))
    
    # Calculate the Log2 normalization of the abundances (Log2 Abundance - Log2 Average)
    data_log2_normalized <- select(data_log2, grep("LOG2 ", names(data_log2))) - data_log2$LOG2_AVG
    names(data_log2_normalized) <- gsub("LOG2", "Log2 Normalized", names(data_log2_normalized))
    data_log2_normalized <- bind_cols(select(data, "unique_id"),
                                      data_log2_normalized)
    data_log2_normalized <- left_join(head, data_log2_normalized, by = "unique_id")
    return(data_log2_normalized)
}