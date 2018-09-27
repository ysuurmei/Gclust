# Loading required packages ----
library(reshape)
library(dplyr)
library(Matrix)
library(fpc)
library(igraph)
library(RColorBrewer)
library(ggplot2)
library(ggforce)
library(GGally)
library(intergraph)
library(data.table)
library(arules)

# Function definitions ----

percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

ClustPlot <- function(n, x = 2, results, df, colcutoff = 0.3){
 
  tbl <- sort(table(results$cluster)[table(results$cluster)>x], decreasing = T)
  temp1 <- results[results$cluster == names(tbl)[[n]],]
  temp2 <- df[df$Order_number %in% temp1$Order_number,]
  temp3 <- subset(temp2, select = c("Order_number", "Item_description"))
  temp4a <- data.table::melt(temp3, id = "Order_number", variable.name = "Var", value.factor = F)
  
  temp4 <- dcast(temp4a, Order_number ~ value, fun.aggregate = length)
  
  rownames(temp4) <- temp4$Order_number
  temp4 <- t(matrix(as.numeric(unlist(temp4)),nrow=nrow(temp4), dimnames = list(rownames(temp4), colnames(temp4))))
  print(ncol(temp4[-1,]))
  temp4 <- temp4[-1,]
  
  
  temp4 <- temp4[, order(colSums(temp4), decreasing = T)]
  temp4[!rowSums(temp4)>=colcutoff*ncol(temp4) & temp4 > 0] <- -1
  
  
  
  
  
  temp5 <- data.table::melt(as.matrix(temp4), id = temp4$Order_number, variable.name = "Var")
  
  srt <- temp5 %>% group_by(Var1) %>% summarize(sum = sum(value))
  
  temp5$Var1 <- factor(temp5$Var1, levels=(temp5$Var1)[order(srt$sum, decreasing = F)])
  
  temp5 <- temp5
  
  
  if(nrow(temp4)> 30){
    warning("Large number of rows, only plotting 20 most significant items", call. = nrow(temp4) > 30)
    print("4")
    if(!is.null(nrow(temp4))){
      temp <- rowSums(temp4)[order(rowSums(temp4), decreasing = T)]
    }
    keeps <- c(names(temp[temp>0]), names(temp[(length(temp)-20):length(temp)]))
    temp5 <- temp5[temp5$Var1 %in% keeps,]
  }
  
  
  plt <- ggplot(temp5, aes(as.numeric(as.factor(Var2)), Var1, fill = value)) + 
    theme_minimal()+
    geom_raster() + 
    geom_tile(colour="white" ,alpha = 0.1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title.x =  element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "grey50", fill=NA, size=1)
    )+ 
    scale_fill_gradient2(low="tomato", mid = "white",  high="#00A600FF")+
    scale_x_discrete(breaks = 1:length(unique(temp5$Var2)),
                     labels = unique(temp5$Var2))
  
  
  tblout <- temp2 %>% group_by(Order_number) %>%
    summarize(OrderSize = n()
    )
  
  print(tblout, n = 50)
  
  output <- list(clustertable = tbl, infotable = tblout, plot = plt)
  print(plt)
  return(output)
}

jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( J )
}

clustdist <- function(results, df){
  tbl <- table(results$cluster)
  
  tblout <- data.frame(matrix(nrow = 1, ncol = 2))
  names(tblout) <- c("Cluster", "Distance")
  
  for (i in names(tbl)) {
    temp <- results[results$cluster == i,]
    
    pivot <- subset(df, select = c("Order_number", "Item"))[df$Order_number %in% temp$Order_number,]
    pivot <- data.table(pivot)
    pivot <- melt(pivot, id = "Order_number")
    pivot <- dcast(pivot, Order_number~value, fun = length)
    dist <- as.matrix(pivot[,2:ncol(pivot)])
    dist[dist>0] <- 1
    dist <- jaccard(dist)
    temp2 <- c(i, signif(mean(dist@x[upper.tri(dist)], na.rm =T), digits = 2))
    tblout <- rbind(tblout, temp2)
  }

  if (nrow(tblout) <= 1) stop("No significant clusters, consider reducing 'minMember'")
  
  
  tblout <- tblout[-1,]
  tblout <- cbind(tblout, as.vector(tbl))
  names(tblout) <- c("Cluster", "Distance", "n")
  tblout <- tblout[order(tblout$n, decreasing = T),]
  tblout$cumsum <- cumsum(tblout$n)
  rownames(tblout) <- 1:nrow(tblout)
  
  u_order <- length(unique(df$Order_number))
  n_order <- sum(tblout$n) 
  prct <- sum(tblout$n)/length(unique(df$Order_number))
  n_clust <- length(tblout$Cluster)
  avg_dist <- as.numeric(tblout$Distance) %*% as.numeric(tblout$n/sum(tblout$n))
  
  param <- data.frame(u_order = u_order, n_order = n_order, prct = prct, n_clust = n_clust, avg_dist = avg_dist)
  
  cat("\n", n_order, "Out of", u_order, "Orders were clustered (", percent(prct, digits = 1),")" )
  cat("\nThe Orders are distributed over ", n_clust, "clusters, with an average within cluster similarity of", percent(avg_dist, digits = 2), "\n")
  
  
  return(list(results = tblout, summarystats = param))
}

to_adj <- function(df, start){
  pivot <- subset(df, select = c("Order_number", "Item"))
  pivot <- data.table(pivot)
  pivot <- melt(pivot, id = "Order_number")
  pivot <- dcast(pivot, Order_number~value, fun = length)
  dist <- pivot[,2:ncol(pivot)]
  
  ind <- dist > 1
  dist[ind] <- 1
  
  dist <- to_sparse(dist)
  dist <- jaccard(dist)
  
  ind <- dist@x <= start
  dist@x[ind] <- 0
  
  return(list(dist = dist, names = pivot[,1]))
}

louvain <- function(adj_matrix, names, minMember = 0){
  
  network <- graph_from_adjacency_matrix(adj_matrix, diag = F, mode = "max", weighted = T) #Note undirected was default
  clst <- cluster_louvain(network, weights = NA)
  
  output <- data.frame(Order_number = names, cluster = clst$membership, stringsAsFactors = F)
  
  keeps <- table(output$cluster)[table(output$cluster)>= minMember]
  output <- output[output$cluster %in% names(keeps),]
  return(output)
  
}

optimizer <- function(df, min, max, step, minMember = 5){
    
  adj <- to_adj(df, start = min)
  test <- data.frame()
    
  for (i in seq(min, max, step)){
    cat("\n i = ", i)
    ind <- adj$dist@x <= i
    adj$dist@x[ind] <- 0
    results <- louvain(adj_matrix = adj$dist, names = adj$names, minMember = minMember)
    summary <- clustdist(results, df = df)
    test <- rbind(test, cbind(summary$summarystats, i))

  }
  
  optim = test[match(max(test$n_clust), test$n_clust),6]
  mx <- 1.75*max(test$n_clust)
  
  plot <- ggplot(test, aes(x = i)) + 
    geom_line(aes(y = prct, col = "blue")) +
    geom_line(aes(y = avg_dist, col = "darkgreen")) +
    scale_y_continuous(sec.axis = ~.*mx) +
    geom_line(aes(y = n_clust/mx, col = "red")) +
    geom_vline(xintercept = optim, color = "grey") +
    geom_text(aes(x=optim+0.07, label= paste("Maximum @ ", optim), y=0.9), colour="black")+
    scale_colour_manual(name = "Legend",
                        values =c('blue'='blue','red'='red', 'darkgreen' ='darkgreen'), labels = c('%Orders clustered','Avg Clust simil', 'N clusters')) + 
    theme(legend.position = c(0.8, 0.98), 
          legend.justification = c(0, 1))
  print(plot)
  return(list(table = test, plot = plot))
}

to_sparse <- function(d_table){
  
  i_list <- lapply(d_table, function(x) which(x != 0))
  counts <- unlist(lapply(i_list, length), use.names = F)
  
  sparseMatrix(
    i = unlist(i_list, use.names = F),
    j = rep(1:ncol(d_table), counts),
    x = unlist(lapply(d_table, function(x) x[x != 0]), use.names = F),
    dims = dim(d_table),
    dimnames = list(NULL, names(d_table)))
}

# Loading the dataset ----

  df <- read.csv("OnlineRetail.csv", encoding = "UTF-8_BOM")
  
  names(df)[1:3] <- c("Order_number", "ItemA", "Item_description")
  
  df$Item <- gsub("[^0-9]", "", df$ItemA)
  df[nchar(df$Item)!=5,]$Item <- df[nchar(df$Item)!=5,]$ItemA
  
  df <- df %>% group_by(Order_number) %>% mutate(n = n())
  
  keeps <- df[df$n >= 15 & df$n <= 50,]
  df <- df[df$Order_number %in% keeps$Order_number,]
  
  
# Running the code for the Louvain algorithm ----
  test <- optimizer(df = df, min = 0.1, max = 0.35, step = 0.05, minMember = 50)
  
  adj <- to_adj(df, start = 0.20) 
  
  results <- louvain(adj_matrix = adj$dist, names = adj$names, minMember = 50)
  summary <- clustdist(results, df = df)
  summary
  
# Plotting the clustering outcomes ----
  outcome <- ClustPlot(5, 4, df = df, results = results, colcutoff = 0.2)

