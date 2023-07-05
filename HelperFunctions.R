formatData <- function(data_Path, labels, dataset) {
  # This function formats data that is used for feature selection methods
  
  # Create the fullpath 
  path <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\data\\", data_Path)
  
  # Read the data from the CSV file 
  data <- read.csv(file = path, header = FALSE)
  
  # Extract the row names.
  row_names <- data[, 1]
  
  # Extract the column names.
  col_names <- data[1, ]
  
  # Transpose the data to interchange rows and columns.
  tdata <- t(data)
  
  # Loop through each label and assign it to the corresponding row in the transposed data.
  for (i in 1:dim(labels[1])[1]) {
    tdata[i] <- labels[i, 1]
  }
  
  # Set the column names of the transposed data as the first row.
  colnames(tdata) <- tdata[1, ]
  
  # Create the output path for the labeled data.
  outputPath1 <- paste0("Labeled_", data_Path)
  output_Path2 <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\Labeled_data\\", outputPath1)
  outputPath <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\Labeled_Data\\", data_Path)
  
  # Write the transposed data with labels to a CSV file.
  write.csv(tdata, output_Path2, row.names = FALSE)
  
  # Return the  data.
  return(tdata)
  
}

formatData_Other_KBest <- function(data_Path, labels, dataset) {
  #This function is also used for formatting data for feature selection, but for different datasets
  
  # Create the full path.
  path <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\data\\", data_Path)
  
  # Read the data from the CSV file.
  data <- read.csv(file = path, header = FALSE)
  
  # Extract the row names.
  row_names <- data[, 1]
  
  # Extract the column names.
  col_names <- data[1, ]
  
  # Transpose the data to interchange rows and columns.
  tdata <- t(data)
  
  # Loop through each label and assign it to the corresponding row in the transposed data.
  for (i in 1:length(labels)) {
    tdata[i, 1] <- labels[[i]]
  }
  
  # Set the column names of the transposed data as the first row of the transposed data.
  colnames(tdata) <- tdata[1, ]
  
  # Create the output path 
  output_Path <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\Labeled_data\\Labeled_", data_Path)
  
  # Write the transposed data with labels to a CSV file.
  write.csv(tdata, output_Path, row.names = FALSE)
  
  # Return data.
  return(tdata)
}

get_Labels_From_PKL <- function(dataset){
  #Returns labels from a PKL file
  
  #Set the datapath to get PKL files from
  datapath <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\cancer_", dataset)
  datapath <- paste0(datapath, ".pkl")
  
  #Load the labels
  label_Dataset <- py_load_object(datapath) 
  label <- c(label_Dataset[[1]], label_Dataset[[2]])
  
  #Format the labels to get correct ground labels.
  if (dataset == "BRCA"){
    label <- as.list(ifelse(label == 1, 0, 1))
  }
  if(dataset == "KIRC"){
    label <- as.list(ifelse(label == 0, 1, 0))
  }
  if (dataset == "LIHC"){
    label <- as.list(ifelse(label == 11, 0, 1))
  }
  if (dataset == "LUAD"){
    label <- as.list(ifelse(label == 1, 1, 0))
  }
  label <- c("class", label)
  
  #returns the labels
  return (label)
}

basic_MDICC_LB <- function(dataset, evaltype){
  # This function runs MDICC with Boruta and Lasso. 
  
  # Set up path and initialization of parameters
  path <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\Lasso_Boruta_Selection\\", dataset, "\\data_", evaltype,"\\")
  setwd(path) #data path
  list <- list.files()
  data <- data.frame()
  data1 <- list()
  X <- list()
  
  # Read data
  for(i in list){
    path <- i
    data <- read.csv(file = path, header = TRUE)
    # rownames(data) <- data$X
    data <- data[-1]
    data11 <- as.matrix(data)
    data1[[i]] = scale(data11, center=TRUE, scale=TRUE) 
    data2 = t(data1[[i]])
    d1 = dist(data2)
    d1 = as.matrix(d1)
    X[[i]] <- d1
  }
  
  # parameter setting
  k1 = 18 # the neighbor of affinity matrix
  k2 = 42 # 
  k3 = 2  # number of cluster
  c  = 3  # c = k3(c>2) or c = 3(c=2)
  
  # calculate affinity matrix
  aff = list()
  for(i in 1:3){
    a = as.matrix(X[[i]])
    xxx = testaff(a,k1)
    aff[[i]] = xxx
  }
  
  # network fusion
  test = MDICC(aff,c = c,k = k2)
  test_S = as.matrix(test)
  
  #Get scores
  score <- MDICCscore(test_S, k3, paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\label.csv"), 'label2')
  names(score) = c('RI','ARI','NMI','Accu','F1')
  label = MDICClabel(test_S,k3)
  MDICCresult = list()
  MDICCresult[['score']] = score
  MDICCresult[['label']] = label
  print(MDICCresult)
  return(MDICCresult)
}

basic_MDICC <- function(dataset){
  #This function is copied and implemented straight from MDICC

  setwd(paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\data")) #data path
  list <- list.files()
  data <- data.frame()
  data1 <- list()
  X <- list()
  
  for(i in list){ba
    path <- i
    ndata <- read.csv(file = path, header = TRUE)
    # rownames(data) <- data$X
    data <- ndata[-1]
    data11 <- as.matrix(data)
    data1[[i]] = scale(data11, center=TRUE, scale=TRUE) 
    data2 = t(data1[[i]])
    d1 = dist(data2)
    d1 = as.matrix(d1)
    X[[i]] <- d1
  }
  
  
  # parameter setting
  k1 = 18 # the neighbor of affinity matrix
  k2 = 42 # 
  k3 = 2  # number of cluster
  c  = 3  # c = k3(c>2) or c = 3(c=2)
  
  # calculate affinity matrix
  aff = list()
  for(i in 1:3){
    a = as.matrix(X[[i]])
    xxx = testaff(a,k1)
    aff[[i]] = xxx
  }
  
  # network fusion
  test = MDICC(aff,c = c,k = k2)
  test_S = as.matrix(test)
  

  score <- MDICCscore(test_S, k3, paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\label.csv"), 'label2')
  names(score) = c('RI','ARI','NMI','Accu','F1')
  label = MDICClabel(test_S,k3)
  MDICCresult = list()
  MDICCresult[['score']] = score
  MDICCresult[['label']] = label
  return(MDICCresult)
}

basic_MDICC_RFECV <- function(dataset){
  #Function to run MDICC for RFECV datasets
  
  get_Results <- function(dataset, j){
    #Helper function to run MDICC on dataset number j for dataset.
    path <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\fs_", dataset, "\\RS", j, "\\")
    print(path)
    setwd(path) #data path
    list <- list.files()
    data <- data.frame()
    data1 <- list()
    X <- list()
    
    for(i in list){
      path <- i
      data <- read.csv(file = path, header = TRUE)
      data <- data[-1]
      data <- as.data.frame(t(data))
      # rownames(data) <- data$X

      data11 <- as.matrix(data)
      data1[[i]] = scale(data11, center=TRUE, scale=TRUE) 
      data2 = t(data1[[i]])
      d1 = dist(data2)
      d1 = as.matrix(d1)
      X[[i]] <- d1
    }
    
    
    # parameter setting
    k1 = 18 # the neighbor of affinity matrix
    k2 = 42 # 
    k3 = 2  # number of cluster
    c  = 3  # c = k3(c>2) or c = 3(c=2)
    
    # calculate affinity matrix
    aff = list()
    for(i in 1:3){
      a = as.matrix(X[[i]])
      xxx = testaff(a,k1)
      aff[[i]] = xxx
    }
    
    # network fusion
    test = MDICC(aff,c = c,k = k2)
    test_S = as.matrix(test)
    
    
    score <- MDICCscore(test_S, k3, paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\label.csv"), 'label2')
    names(score) = c('RI','ARI','NMI','Accu','F1')
    label = MDICClabel(test_S,k3)
    MDICCresult = list()
    MDICCresult[['score']] = score
    MDICCresult[['label']] = label
    return(MDICCresult)
  }
  output = list()
  
  #Run MDICC for RFECV
  for(iteration in seq(1, 2, by = 1)){
    output[[iteration]] <- get_Results(dataset, iteration)
  }
  return(output)
}

getKBestScores <- function(dataset, division){
  #Function to run MDICC for Kbest scores in a loop.
  
  #Get data and set working directory
  wd_Labeled_data = paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\Labeled_data")
  # wd_Labeled_data = paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\data")
  setwd(wd_Labeled_data)
  
  #Initialize parameters
  list <- list.files()
  data <- data.frame()
  data1 <- list()
  X <- list()
  
  #Import datasets to run Kbest on
  dataframes <- list()
  for(i in list){
    print(paste0("Loading dataset ", i, " from ", wd_Labeled_data))
    path <- i
    dataframe <- read.csv(file = path, header = TRUE)
    dataframes[[i]] <- dataframe
  }
  # dataframes <- list(labeled_Gene, labeled_Methyl, labeled_miRNA)
  
  # parameter setting
  k1 = 2 # the neighbor of affinity matrix
  k2 = 42 #
  k3 = 2 # number of cluster
  c = 3 # c = k3(c>2) or c = 3(c=2)
  
  
  start <- 1
  #Number of divisions
  if((as.integer(1800/division)) <= 2){
    start <- 3
  }
  
  # Create an empty list to store the results
  results_list <- list()
  
  # Loop the number of divisions
  for (iter in seq(start, division, by = 1)) {
    print(paste0("Iteration: ", iter - start + 1))
    data <- data.frame()
    data1 <- list()
    ndata <- data.frame()
    X <- list()
    
    # Loop over the 3 datasets
    for (i in list) {
      
      # Store dataframe in ndata
      ndata <- dataframes[[i]]
      
      #Calculate number of features based on size of dataframe
      nbr_col <- ncol(ndata)
      k_value <- (nbr_col/division) * iter
      
      # Apply Kbest features from sklearn
      data <- feature_SelectionKBest(ndata,k_value-1, i, iter, dataset)
      
      # MDICC magic
      data11 <- as.matrix(data)
      data1[[i]] = scale(data11, center=TRUE, scale=TRUE)
      data2 = t(data1[[i]])
      d1 = dist(data2)
      d1 = as.matrix(d1)
      X[[i]] <- d1
      # print(i)
      # print(as.integer(k_value -1))
    }
    # print(d1[1,])
    print(dim(d1))
    
    # Calculate affinity matrix
    aff <- list()
    for (i in 1:3) {
      a <- as.matrix(X[[i]])
      xxx <- testaff(a, k1)
      aff[[i]] <- xxx
    }
    
    # Network fusion
    test <- MDICC(aff, c = c, k = k2)
    test_S <- as.matrix(test)
    
    # Calculate results
    score <- MDICCscore(test_S, k3, paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", dataset, "\\label.csv"), 'label2')
    names(score) <- c('RI', 'ARI', 'NMI', 'Accu', 'F1')
    label <- MDICClabel(test_S, k3)
    MDICCresult <- list()
    MDICCresult[['score']] <- score
    MDICCresult[['label']] <- label
    
    # Store the results in the list
    results_list[[as.character(iter)]] <- MDICCresult
    print(MDICCresult)
    print(d1[1,])
  }
  return(results_list)
}

output_Single_Result <- function(results_list, division, score, dataset){
  #Function used to make graphical representations for single results
  
  #Get scores
  first_scores <- sapply(results_list, function(result) result[1])
  first_scores_RI <- sapply(first_scores, function(result) result[score])
  meassure <- ""
  
  #Get score
  if (score == 1) {
    meassure <- "AR Score"
  } else if (score == 2) {
    meassure <- "ARI Score"
  } else if (score == 3) {
    meassure <- "NMI Score"
  } else if (score == 4) {
    meassure <- "Accuracy Score"
  } else if (score == 5) {
    meassure <- "F1 Score"
  } else {
    meassure <- "Invalid score"
  }
  
  # plot graph
  plot(seq(1, division, by = 1), first_scores_RI, type = "l", col = "blue", xlab = "Iteration (K)", ylab = meassure)
  legend("bottomright", legend = c("y1"), col = c("blue"), lty = 1)
  title(paste0(meassure, " for ", toupper(dataset), " dataset"))
  
  #Order results
  top_10_positions <- order(first_scores_RI, decreasing = TRUE)[1:division]
  
  # Get the maximum values of the top 10
  top_10_max_values <- first_scores_RI[top_10_positions]
  
  # Print the positions and maximum values
  # print(top_10_positions)
  # print(top_10_max_values)
  
  #Show Results
  # View(results_list)
}

output_Result <- function(result_BRCA, result_LIHC, score, division){
  #Function used for displaying scores of Kbest
  
  # result_BRCA <- brca_K_Best
  # result_LIHC <- lihc_K_Best
  # score <- 2
  # division <- division
  
  #Get best scores for BRCA and LIHC
  brca_scores <- sapply(result_BRCA, function(result) result[1])
  brca_scores_first <- sapply(brca_scores, function(result) result[score])
  
  lihc_scores <- sapply(result_LIHC, function(result) result[1])
  lihc_scores_first <- sapply(lihc_scores, function(result) result[score])
  
  #get top features and positions
  brca_top_features <- get_Top_Features(brca_scores_first, division)
  lihc_top_features <- get_Top_Features(lihc_scores_first, division)
  

  #Set the right label for the graph
  if (score == 1) {
    measure <- "AR Score"
  } else if (score == 2) {
    measure <- "ARI Score"
  } else if (score == 3) {
    measure <- "NMI Score"
  } else if (score == 4) {
    measure <- "Accuracy Score"
  } else if (score == 5) {
    measure <- "F1 Score"
  } else {
    measure <- "Invalid score"
  }
  
  #Set up dataframe and graph
  for(dev in dev.list()){dev.off(dev)}
  data <- data.frame(x = seq(1, division-2, by = 1), y1 = brca_scores_first, y2 = lihc_scores_first)
  
  #Set limit in case the number of division was smaller than 50.
  limit <- division
  if(division >50){
    limit <- 50
  }
  
  #Get the highest score
  pos_brca <- which(brca_scores_first == max(brca_scores_first))[1]
  pos_lihc <- which(lihc_scores_first == max(lihc_scores_first))[1]
  library(ggplot2)
  
  # Plot graph
  plot <- ggplot(data, aes(x = x)) +
    geom_line(aes(y = y1, color = "BRCA"), linetype = "solid", show.legend = T) +
    geom_line(aes(y = y2, color = "LIHC"), linetype = "solid", show.legend = T) +
    labs(title = paste0(measure, " Performance of BRCA and LIHC datasets"), 
         subtitle = measure, x = "Iteration (K)", y = measure) +
    theme_bw() +
    theme(axis.title = element_text(), legend.position = "bottom") +
    scale_x_continuous(breaks = seq(0, division, 100), limits = c(1, division)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    geom_vline(aes(xintercept = pos_brca, color = "BRCA"), linetype = "dashed" , show.legend = F) +
    geom_vline(aes(xintercept = pos_lihc, color ="LIHC"), linetype = "dashed", show.legend = F) +
    scale_color_manual(values = c("BRCA" = "blue", "LIHC" = "red"), name = "Dataset: ") +  
    guides(color = guide_legend(title = "Datasets:")) +  
    geom_point(aes(x=x, y = NaN)) +
    annotate("text", x = lihc_top_features$index[1] + 50, y = 0 , label = paste0("K = ", lihc_top_features$index[1] ), color = "red") +
    annotate("text", x = brca_top_features$index[1] - 75, y = 0 , label = paste0("K = ", brca_top_features$index[1]), color = "blue")

  #Print graph (print needed for ggplot to work)
  print(plot)
  
}

output_Result_subset <- function(result_BRCA, result_LIHC, score, division, subset){
  #Helper function used to display results for a subset of original data.
  output_Result(result_BRCA[1:subset], result_LIHC[1:subset], score, subset + 2)
}

plot_Bars_RFECV <- function(brca_RFECV, lihc_RFECV, plot, lihcvalue){
  #Function used to display results of RFECV
  
  #Get best scores for BRCA and LIHC
  first_elements <- lapply(brca_RFECV, function(x) x[[1]])
  first_elementsX <- lapply(first_elements, function(x) x[[plot]])
  
  first_elements_lihc <- lapply(lihc_RFECV, function(x) x[[1]])
  first_elements_lihc_X <- lapply(first_elements_lihc, function(x) x[[plot]])
  
  #Get values to plot
  values = c()
  for(value in first_elementsX){
    values <- c(values, value)
  }
  for(value in first_elements_lihc_X){
    values <- c(values, value)
  }
  
  #Set up dataframe and vertical line.
  df <- data.frame(
    category = factor(c("BRCA", "LIHC", "BRCA 1", "BRCA 2", "LIHC 1", "LIHC 2")),
    value = c(1, lihcvalue,values)
  )
  sorted_column <- sort(df$value, decreasing = TRUE)
  second_highest <- sorted_column[2]
  
  plotvalue <- sprintf("%0.3f", second_highest)

  
  # Create the plot
  ggplot(df, aes(x=category, y=value, fill=category)) +
    geom_bar(stat="identity", color = "black") +
    labs(x = "RFECV with random state", y = paste0(plot, "score"), fill = "Category", title = paste0("Bar plot values of ", plot)) +
    theme_minimal() + scale_fill_viridis_d(option = 'D', direction = -1) +
    geom_hline(yintercept = second_highest, linetype = "dashed") + 
    geom_point(aes(x=category, y=value)) + 
    annotate("text", x = 3, y = second_highest + 0.05, label = plotvalue)
}

plot_Bars_LB <- function(brca_boruta, lihc_boruta,plot, type){
  #Function used to plot results of Lasso and Boruta

  #Get values from results
  values <- c(brca_boruta$score[plot], lihc_boruta$score[plot])
  
  #Set up dataframe and vertical line.
  df <- data.frame(
    category = factor(c("BRCA", "LIHC")),
    value = values
  )
  
  plotvalue <- sprintf("%0.3f", max(values))
  
  # Create the plot
  ggplot(df, aes(x=category, y=value, fill=category)) +
    geom_bar(stat="identity", color = "black") +
    labs(x = type, y = paste0(plot, "score"), fill = "Category", title = paste0("Bar plot values of ", plot)) +
    theme_minimal() + scale_fill_viridis_d(option = 'D', direction = -1) +
    geom_hline(yintercept = max(values), linetype = "dashed") + 
    geom_point(aes(x=category, y=value)) + 
    annotate("text", x = 1.5, y = max(values) + 0.05, label = plotvalue)
}

plot_Bars_ALL <- function(names, values, plot){
  #Function used for a generic plot. This function is mainly used during the research to get a better understanding of results.
  
  #Set up dataframe
  df <- data.frame(
    category = factor(names),
    value = values
  )
  
  # Create the plot
  ggplot(df, aes(x=category, y=value, fill=category)) +
    geom_bar(stat="identity", color = "black") +
    labs(x = paste0(plot, " scores"), y = paste0(plot, "score"), fill = "Category", title = paste0("Bar plot values of ", plot)) +
    theme_minimal() + scale_fill_viridis_d(option = 'D', direction = -1) 
}

plot_Bars_size <- function(values, plot){
  #Function used to plot number of features per dataset. These values are hard coded in comments due to a lack of time.
  
  names <- c("Gene", "Methyl", "miRNA")
  ##  BRCA
  # RFECV1
  # Gene:15661
  # Methyl: 1133
  # miRNA:1496
  # 
  # RFECV2
  # Gene: 15676
  # Methyl: 13849
  # miRNA:168
  # 
  # Lasso
  # Gene: 23
  # Methyl: 21
  # miRNA: 23
  # 
  # Boruta
  # Gene:55
  # Methyl: 65
  # miRNA:37
  
  ##  LIHC
  # RFECV1
  # Gene:11662
  # Methyl: 904
  # miRNA:1007
  # 
  # RFECV2
  # Gene: 20110
  # Methyl: 2212
  # miRNA:1043
  # 
  # Lasso
  # Gene: 24
  # Methyl: 33
  # miRNA: 7
  # 
  # Boruta
  # Gene:74
  # Methyl: 80
  # miRNA:65
  
  # plot_Bars_size(c(15661, 1133, 1496), " RFECV with random state 1 for BRCA")
  # plot_Bars_size(c(15676, 13849, 168), " RFECV with random state 2 for BRCA")
  # plot_Bars_size(c(11662, 904, 1007), " RFECV with random state 1 for LIHC")
  # plot_Bars_size(c(20110, 2212, 1043), " RFECV with random state 2 for LIHC")
  # plot_Bars_size(c(55, 65, 37), " Boruta for BRCA")
  # plot_Bars_size(c(74, 80, 65), " Boruta for LIHC")
  # plot_Bars_size(c(23, 21, 23), " Lasso for BRCA")
  # plot_Bars_size(c(24, 33, 7), " Lasso for LIHC")
  
  #Create dataframe
  df <- data.frame(
    category = factor(names),
    value = values
  )
  
  # Create the plot
  ggplot(df, aes(x=category, y=value, fill=category)) +
    geom_bar(stat="identity", color = "black") +
    labs(x = paste0(plot, " dataset"), y = "Features", fill = "Category", title = paste0("Number of features per dataset of", plot)) +
    theme_minimal() + scale_fill_viridis_d(option = 'D', direction = -1) 
}

get_Top_Features <- function(scores, division){
  #Helper function to get top features.
  
  #order features
  top_positions <- order(scores, decreasing = TRUE)[1:division]
  
  # Get the maximum values of the top 10
  max_values <- scores[top_positions]
  
  # Print the positions and maximum values
  print("BRCA Scores and positions")
  print(top_positions)
  
  print("LIHC Scores and positions")
  print(max_values)
  return(list(index = top_positions, value = max_values))
}

lasso_Boruta_Fs <- function(datapath){
  #Function used to create datasets after Boruta and Lasso have been applied.
  
  #Get files
  setwd(paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\", datapath, "\\data"))
  dirs <- list.files()
  
  #Read in original data
  gene <- read.csv(file = dirs[1], header = TRUE)
  methyl <- read.csv(file = dirs[2], header = TRUE)
  miRNA  <- read.csv(file = dirs[3], header = TRUE)
  
  #Create subsets with Boruta/Lasso applied
  create_Subset(gene,dirs[1], "Boruta", datapath)
  create_Subset(methyl,dirs[2], "Boruta", datapath)
  create_Subset(miRNA,dirs[3], "Boruta", datapath)
  
  create_Subset(gene,dirs[1], "Lasso", datapath)
  create_Subset(methyl,dirs[2], "Lasso", datapath)
  create_Subset(miRNA,dirs[3], "Lasso", datapath)
}

create_Subset <- function(dataset, file, method, datapath){
  #Helper function to create a subset after Boruta or Lasso determined features.
  
  #Set working directory
  # setwd(paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\Lasso_Boruta_Selection\\", method, "\\"))
  
  #Set rownames to be the first row
  rownames(dataset) <- dataset[,1]
  
  #Delete first row
  dataset <- dataset[,-1]
  
  #Set Path 
  path <- paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\Lasso_Boruta_Selection\\", datapath, "\\", method, "\\", method, "_", file)
  # print(path)
  
  #Read in headers after feature selection has been applied
  miRNA_headers <- read.csv(file = path, header = TRUE)
  # print(miRNA_headers)
  rownames(miRNA_headers) <- miRNA_headers[,1]
  miRNA_headers <- miRNA_headers[,-1]
  
  #Methylation changes header names. These lines set them back to match the headers of the original data.
  if(file == "methyl.csv")
    rownames(miRNA_headers) <- gsub("\\.", "-", rownames(miRNA_headers))
  
  #Get common rows of subset and original data
  common_rows <- rownames(dataset) %in% rownames(miRNA_headers)
  
  #Create subset with features selected and create csv file
  subset_miRNA <- dataset[common_rows, ]
  write.csv(subset_miRNA, paste0("C:\\Users\\matth\\anaconda3\\envs\\data1\\data1\\Lasso_Boruta_Selection\\", datapath, "\\data_", tolower(method), "\\", "\\Fs_",method , "_", file))
}