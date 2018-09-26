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
  temp3 <- subset(temp2, select = c("Order_number", "Description"))
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
    temp2 <- c(i, signif(mean(dist[upper.tri(dist)]), digits = 2))
    tblout <- rbind(tblout, temp2)
  }
  
  tblout <- tblout[-1,]
  #rownames(tblout) <- 1:nrow(tblout)
  tblout <- cbind(tblout, as.vector(tbl))
  names(tblout) <- c("Cluster", "Distance", "n")
  tblout$cumsum <- cumsum(tblout$n)
  return(tblout)
}

gclust <- function(df, start = 1, finish = 0.5, step = 0.1, minMember = 10 ){
  clusters <- data.frame()
  i <- 1
  n <- (start-finish)/step
  output <- list()
  
  while (start >= finish){
    pivot <- subset(df, select = c("Order_number", "Item"))
    
    pivot <- pivot[!pivot$Order_number %in% clusters$Order_number,]
    
    pivot <- data.table(pivot)
    pivot <- melt(pivot, id = "Order_number")
    pivot <- dcast(pivot, Order_number~value, fun = length)
    
    dist <- as.matrix(pivot[,2:ncol(pivot)])
    dist[dist>0] <- 1
    
    dist <- jaccard(dist)
    
    ind <- dist@x <= start
    
    dist@x[ind] <- 0
    dist <- as.matrix(dist)
    rownames(dist) <- as.matrix(pivot[,1])
    
    network <- graph_from_adjacency_matrix(dist, diag = F, weighted = T, add.rownames = T)
    results <- data.frame(Order_number = pivot[,1], cluster = components(network)$membership)
    results$cluster <- paste0(LETTERS[[i]], ":", results$cluster)
    
    tbl <- table(results$cluster)[table(results$cluster)>minMember]
    clusters <- rbind(clusters, results[results$cluster %in% names(tbl),])
    
    cat("\nIteration", i, "is completed,", n-i+1, "iterations remaining")
    
    i <- i+1
    start <- start - step
  }
  
  summary <- clustdist(clusters, df)
  keeps <- summary[summary$Distance >= finish,1]
  
  summary <- summary[summary$Cluster %in% keeps,]
  clusters <- clusters[clusters$cluster %in% keeps,]
  
  
  output$clusters <- clusters
  output$summary <- summary
  
  return(output)
}

checkclust <- function (check, check2, var1, var2){
  par(mar = c(5, 5, 3, 5))
  plot(check[[var1]], type ="p",
       col = "grey", ylim = c(0,1))
  abline(lm(check[[var1]]~as.vector(1:nrow(check))), col = "grey")
  par(new = TRUE)
  plot(check[[var2]], type = "l", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "grey", lty = 2)
  axis(side = 4)
  mtext("Cumulative observations", side = 4, line = 3)
  par(new = TRUE)
  plot(check2[[var1]], type ="p", ylab = "Within Cluster Distance (Jaccard)",
       main = "Cluster Distance Plot", xlab = "Cluster",
       col = "blue", ylim = c(0,1))
  abline(lm(check[[var1]]~as.vector(1:nrow(check))), col = "blue")
  par(new = TRUE)
  plot(check2[[var2]], type = "l", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", col = "red", lty = 2)
  axis(side = 4)
  mtext("Cumulative observations", side = 4, line = 3)
  
  
}

clustoptim <- function(results, check, df, cutoff1 = 0.7, cutoff2 = 0.7){
  templist <- data.frame()
  results2 <- results
  
  for (i in unique(check$Cluster)){
    temp1 <- results[results$cluster  == i,]
    temp2 <- df[df$Order_number %in% temp1$Order_number,]
    temp2 <- subset(temp2, select = c("Item", "Order_number"))
    temp2 <- temp2[!duplicated(temp2),]
    temp3 <- (table(temp2$Item)[table(temp2$Item) > cutoff1*length(unique(temp2$Order_number))])
    templist <- bind_rows(templist, temp3)
    
  }
  
  templist[is.na(templist)] <- 0
  templist[templist > 0] <- 1
  rownames(templist) <- unique(check$Cluster)
  templist <- as.matrix(templist)
  dist <- jaccard(templist)
  
  dist[dist < cutoff2] = 0
  
  
  rownames(dist) <- unique(check$Cluster)
  network <- graph_from_adjacency_matrix(as.matrix(dist), diag = F, weighted = T, add.rownames = T)
  member <- data.frame(Name = unique(check$Cluster), cluster = components(network)$membership)
  
  names <- names(table(member$cluster)[table(member$cluster)>1])
  names <- member[member$cluster %in% names,]
  
  for (i in unique(names$cluster)){
    rename <- names[names$cluster == i,]
    
    results2[results2$cluster %in% as.character(rename$Name),2] <- as.character(i)
    
  }
  results2$cluster <- as.factor(results2$cluster)
  levels(results2$cluster) <- 1:length(unique(results2$cluster))
  results2 %>% group_by(cluster) %>% summarise(n=n()) %>% ungroup() %>% arrange(-n)
  return(results2)
}

removeWords <- function(str, stopwords) {
  x <- unlist(strsplit(str, " "))
  paste(x[!x %in% stopwords], collapse = " ")
}


# Loading the dataset ----

#Load the online retail order dataset
  
  df <- read.csv("OnlineRetail.csv", encoding = "UTF-8_BOM")
  head(df)
  
  df <- df[,c(1,2,3)]
  names(df)[1:2] <- c("Order_number", "ItemA", "Item_description")
  
  df$Item <- gsub("[^0-9]", "", df$ItemA)
  
  df <- df %>% group_by(Order_number) %>% mutate(n = n())
  
  keeps <- df[df$n >= 3 & df$n <= 25,]
  df <- df[df$Order_number %in% keeps$Order_number,]
  
  
  # Running the code for the Louvain algorithm ----
  pivot <- subset(df, select = c("Order_number", "Item"))
  pivot <- data.table(pivot)
  pivot <- melt(pivot, id = "Order_number")
  pivot <- dcast(pivot, Order_number~value, fun = length)
  dist <- pivot[,2:ncol(pivot)]
  
  ind <- dist > 1
  dist[ind] <- 0
  gc()
  dist <- jaccard(as.matrix(dist))
  
  start <- 0.35
  ind <- dist@x <= start
  dist@x[ind] <- 0
  dist <- as.matrix(dist)
  rownames(dist) <- as.matrix(pivot[,1])
  
  network <- graph_from_adjacency_matrix(dist, diag = F, weighted = T, add.rownames = T, mode = "undirected")
  test <- cluster_louvain(network, weights = NA)
  
  results <- data.frame(Order_number = rownames(dist), cluster = test$membership, stringsAsFactors = F)
  head(results)
  
  keeps <- table(results$cluster)[table(results$cluster)>=3]
  results <- results[results$cluster %in% names(keeps),]
  
  summary <- clustdist(results, df = df)
  
  
  sum(summary$n)
  summary
  as.numeric(summary$Distance) %*% as.numeric(summary$n/sum(summary$n))
  cat(sum(summary$n), "Out of", length(unique(df$Order_number)), "Orders were clustered (", percent(sum(summary$n)/length(unique(df$Order_number)), digits = 1),")" )
  length(summary$Cluster)

outcome <- ClustPlot(2, 4, df = df, results = results, colcutoff = 0.3)
outcome$plot
outcome$clustertable


# Running the code for the Gclust algorithm ----

results <- gclust(df, start = 1, finish = 0.2, step = 0.05, minMember = 3)

results2 <- clustoptim(results$clusters, results$summary, df, cutoff1 = 0.7, cutoff2 = 0.7)

check2 <- clustdist(results2, df)

checkclust(check, check2, "Distance", "cumsum")

keeps <- check2[check2$Distance > 0.3, 1]
results2 <- results2[results2$cluster %in% keeps,]

outcome <- ClustPlot(15, 1, results = results2, df = df, colcutoff = 0.5)

outcome$plot
sum(outcome$clustertable)
outcome$clustertable

