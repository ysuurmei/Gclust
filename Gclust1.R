# Loading required packages ----
library(reshape)
library(dplyr)
library(Matrix)
library(igraph)
library(ggplot2)
library(data.table)
library(countrycode)


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

  temp4 <- temp4[-1,]
  
  
  temp4 <- temp4[, order(colSums(temp4), decreasing = T)]
  temp4[temp4>1] <- 1
  temp4[!rowSums(temp4)>=colcutoff*ncol(temp4) & temp4 > 0] <- -1
  
  
  
  
  
  temp5 <- data.table::melt(as.matrix(temp4), id = temp4$Order_number, variable.name = "Var")
  
  srt <- temp5 %>% group_by(Var1) %>% summarize(sum = sum(value))
  
  temp5$Var1 <- factor(temp5$Var1, levels=(temp5$Var1)[order(srt$sum, decreasing = F)])

  #temp6 <- merge(temp5, subset(temp2, select = c("Item_description", "Order_number", "Usage Category")), 
  #               by.x = c("Var1", "Var2"), by.y = c("Material_description", "Order_number"), all.x = T)
  
  temp6a <- subset(temp2, select = c("Order_number", "CountryCode" ,"InvoiceDate"))
  temp6a <- temp6a[temp6a$Order_number %in% temp5$Var2,]
  temp6a <- temp6a[!duplicated(temp6a),]
  temp6 <- merge(temp5, temp6a, by.x = "Var2",
                 by.y = "Order_number", all.x = T, all.y = F)
  temp6 <- temp6[!duplicated(temp6[,c(1,2,3)]),]
  
  temp6$lbls <- as.factor(paste(temp6$Country, temp6$InvoiceDate))
  temp6$lbls1 <- as.factor(temp6$Country)
  temp6$lbls2 <- as.factor(paste(year(temp6$InvoiceDate), "-", months.Date(temp6$InvoiceDate, abbreviate = T)))
  
  
  if(nrow(temp4)> 30){
    warning("Large number of rows, only plotting 20 most significant items", call. = nrow(temp4) > 30)

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
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
          axis.title.x =  element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(colour = "grey50", fill=NA, size=1)
    )+ 
    scale_fill_gradient2(low="tomato", mid = "white",  high="#00A600FF")+
    scale_x_continuous(breaks = 1:length(unique(temp6$Var2)),
                       labels = temp6$lbls1[seq(1, length(temp6$lbls1), nrow(temp6)/length(unique(temp6$Var2)))],
                       sec.axis = sec_axis(~.,
                                           breaks = 1:length(unique(temp6$Var2)),
                                           labels = temp6$lbls2[seq(1, length(temp6$lbls2), nrow(temp6)/length(unique(temp6$Var2)))]))
  
  
  words <- table(strsplit(paste(temp2$Item_description, collapse = " "), " +"))/nrow(temp2)
  words <- words[!is.na(words)]
  words <- words[order(words, decreasing = T)]
  words <- data.frame(Word = names(words)[1:20], Prct = percent(words[1:20], digits = 1))
  
  tblout <- temp2 %>% group_by(Order_number) %>%
    summarize(
      OrderSize = n(),
      OrderValue = sum(UnitPrice),
      Country = first(CountryCode),
      Date = first(InvoiceDate) 
    )
  sumtable <- data.frame(N_Orders = length(unique(temp2$Order_number)),
                         Avg_value = mean(tblout$OrderValue),
                         Avg_size = mean(tblout$OrderSize))
  
  print(tblout, n = 30)
  
  output <- list(clustertable = tbl, infotable = tblout, plot = plt, words = words, summary = sumtable)
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
  return(list(results = output, network = network, cluster = clst))
  
}

optimizer <- function(df, min, max, step, minMember = 5){
    
  adj <- to_adj(df, start = min)
  test <- data.frame()
    
  for (i in seq(min, max, step)){
    cat("\n i = ", i)
    ind <- adj$dist@x <= i
    adj$dist@x[ind] <- 0
    results <- louvain(adj_matrix = adj$dist, names = adj$names, minMember = minMember)
    summary <- clustdist(results$results, df = df)
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
  
  df$InvoiceDate <- as.Date(as.character(df$InvoiceDate), format = "%d-%m-%y")
  df$CountryCode <- countrycode(df$Country, origin = "country.name", destination = "genc3c")
  
  df <- df[!is.na(df$CountryCode),]
  
  df$Item <- gsub("[^0-9]", "", df$ItemA)
  df[nchar(df$Item)!=5,]$Item <- df[nchar(df$Item)!=5,]$ItemA
  
  
  df <- df %>% group_by(Order_number) %>% mutate(n = n())
  keeps <- df[df$n >= 15 & df$n <= 50,]
  df <- df[df$Order_number %in% keeps$Order_number,]
  
  
# Running the code for the Louvain algorithm ----
  test <- optimizer(df = df, min = 0.1, max = 0.35, step = 0.05, minMember = 50)
  
  adj <- to_adj(df, start = 0.20) 
  
  results <- louvain(adj_matrix = adj$dist, names = adj$names, minMember = 50)
  summary <- clustdist(results$results, df = df)
  summary
  
# Plotting the clustering outcomes ----
  outcome <- ClustPlot(5, 20, df = df, results = results, colcutoff = 0.2)
  outcome$summary
  outcome$words
  plot(as.factor(outcome$infotable$Country))
