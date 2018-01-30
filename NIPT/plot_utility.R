### Scatter plot (1 probe, 1 graph)
# input: AB data.table
# input: probe_id
# output: sactter plot

AB_scatter_plot <- function(plot_dt = AB.dt, pid = "PH_NIPT_210017") {
  plot(x = unlist(plot_dt[probe_id == pid, 2:(samples_Num + 1)]), 
       y = unlist(plot_dt[probe_id == pid, (samples_Num+2):(2*samples_Num+1)]),
       main = pid,
       xlab = "A",
       ylab = "B",
       asp = 1)#,
  # xlim = c(0, 4500),
  # ylim = c(0, 1000))
}


AB_scatter_plot_color_by_clustering <- function(plot_dt, cluster_dt, pid = "PH_NIPT_210018") {
  #clustering 先reshape data.table
  tmp <- subset(cluster_dt, probe_id == pid, select = c(2:(1+2*samples_Num)))
  tmp_A <- tmp[, 1:(length(tmp)/2)]
  tmp_B <- tmp[, (length(tmp)/2+1):length(tmp)]
  names(tmp_B) <- names(tmp_A) # 為了配合rbind時需要對應的column names
  tmp1 <- rbind(tmp_A, tmp_B)
  tmp1 <- t(tmp1)
  # tmp_modeling <- cbind(tmp1, (tmp1[,1] - tmp1[,2]))
  tmp_modeling <- as.data.table(tmp1[,1] - tmp1[,2])
  # colnames(tmp_modeling) <- c("A", "B", "ratio")
  colnames(tmp_modeling) <- c("ratio")
  
  #str(tmp1)
  # md <- Mclust(tmp_modeling, G = 3)
  # md <- kmeans(tmp_modeling, 3)
  # cat("group #:", md$G, "\n")
  #----
  E.dist <- dist(tmp_modeling, method="euclidean")
  # h.E.cluster <- hclust(E.dist, method="centroid")
  # h.E.cluster <- hclust(E.dist, method="ward.D2")
  h.E.cluster <- hclust(E.dist, method="complete")
  
  # plot(h.E.cluster, xlab="歐式距離")
  cut.h.cluster <- cutree(h.E.cluster, k=3)
  #----
  
  
  
  # tmp1 <- cbind(tmp1, md$classification)
  tmp1 <- cbind(tmp1, cut.h.cluster)
  # tmp1 <- cbind(tmp1, md$cluster)
  colnames(tmp1) <- c("A", "B", "cluster")
  cluster_index <- ncol(tmp1)
  
  
  #plot(tmp1[, 1], tmp1[,2], col=tmp1[,3], asp = TRUE)
  
  tmp_ori <- subset(plot_dt, probe_id == pid, select = c(2:(1+2*samples_Num)))
  tmp_ori_A <- tmp_ori[, 1:(length(tmp_ori)/2)]
  tmp_ori_B <- tmp_ori[, (length(tmp_ori)/2+1):length(tmp_ori)]
  names(tmp_ori_B) <- names(tmp_ori_A)
  tmp_ori <- rbind(tmp_ori_A, tmp_ori_B)
  tmp_ori <- t(tmp_ori)
  tmp_ori <- as.data.table(tmp_ori)
  tmp_ori <- cbind(tmp_ori, tmp1[,cluster_index])
  colnames(tmp_ori) <- c("A", "B", "cluster")
  
  # specify the order of the clusters
  c1 <- tmp_ori[cluster==1, ]
  c2 <- tmp_ori[cluster==2, ]
  c3 <- tmp_ori[cluster==3, ]
  
  c1 <- as.numeric(apply(c1[,-3], 2, mean))
  c2 <- as.numeric(apply(c2[,-3], 2, mean))
  c3 <- as.numeric(apply(c3[,-3], 2, mean))
  
  c1 <- c1[1] - c1[2]
  c2 <- c2[1] - c2[2]
  c3 <- c3[1] - c3[2]
  
  clust_order <- order(c(c1,c2,c3), decreasing = TRUE)
  color_flag <- palette()[1:3]
  frac_str <- c("AA","AB","BB")
  # layout(rbind(1,2), heights=c(7,1))
  plot(tmp_ori[,A], tmp_ori[,B], pch = 20, col=tmp_ori[,cluster], asp = TRUE, xlab = "A", ylab = "B", main = pid)
  
  
  legend('topright','groups',legend = frac_str, cex = 0.8, pch = 20 ,col=c(color_flag[clust_order[1]], color_flag[clust_order[2]], color_flag[clust_order[3]]))
  
}




### Plotting for new sampel judgement
#  New sample, judge the A, B percentages
newSampleJudgment_plot <- function(newPoint, pid = "PH_NIPT_210018") {
  tmp <- subset(AB.log2.scale.dt, probe_id == pid, select = c(2:(1+2*samples_Num)))
  tmp_A <- tmp[, 1:(length(tmp)/2)]
  tmp_B <- tmp[, (length(tmp)/2+1):length(tmp)]
  names(tmp_B) <- names(tmp_A) # 為了配合rbind時需要對應的column names
  tmp1 <- rbind(tmp_A, tmp_B)
  tmp1 <- t(tmp1)
  tmp1 <- cbind(tmp1, (tmp1[,1] - tmp1[,2]))
  colnames(tmp1) <- c("A", "B", "ratio")
  
  #str(tmp1)
  
  md <- Mclust(tmp1, G = 3)
  
  # 用組間距離判定是否需改為2群
  # dist_vec <- c(dist(t(mdg3$parameters$mean)))
  # if (sum(dist_vec < 0.5)) md <- Mclust(tmp1, G = 2)
  
  tmp1 <- cbind(tmp1[, -3], md$classification)
  colnames(tmp1) <- c("A", "B", "cluster")
  
  #plot(tmp1[, 1], tmp1[,2], col=tmp1[,3], asp = TRUE)
  
  tmp_ori <- subset(AB.dt, probe_id == pid, select = c(2:(1 + samples_Num*2)))
  tmp_ori_A <- tmp_ori[, 1:(length(tmp_ori)/2)]
  tmp_ori_B <- tmp_ori[, (length(tmp_ori)/2+1):length(tmp_ori)]
  names(tmp_ori_B) <- names(tmp_ori_A)
  tmp_ori <- rbind(tmp_ori_A, tmp_ori_B)
  tmp_ori <- t(tmp_ori)
  tmp_ori <- as.data.table(tmp_ori)
  tmp_ori <- cbind(tmp_ori, tmp1[,3])
  colnames(tmp_ori) <- c("A", "B", "cluster")
  
  plot(tmp_ori[,A], tmp_ori[,B], col=tmp_ori[,cluster], asp = TRUE, xlab = "A", ylab = "B", main = pid, pch = 19)
  points(newPoint[1], newPoint[2], cex = 5, pch = 3)
  
  
  # specify the order of the clusters
  c1 <- tmp_ori[cluster==1, ]
  c2 <- tmp_ori[cluster==2, ]
  c3 <- tmp_ori[cluster==3, ]
  
  c1 <- as.numeric(apply(c1[,-3], 2, mean))
  c2 <- as.numeric(apply(c2[,-3], 2, mean))
  c3 <- as.numeric(apply(c3[,-3], 2, mean))
  
  c1 <- c1[1] - c1[2]
  c2 <- c2[1] - c2[2]
  c3 <- c3[1] - c3[2]
  
  clust_order <- order(c(c1,c2,c3), decreasing = TRUE)
  color_flag <- palette()[1:3]
  frac_str <- c("AA","AB","BB")
  legend('topright','groups', legend=frac_str, cex = 0.8, pch = 20 ,col=c(color_flag[clust_order[1]], color_flag[clust_order[2]], color_flag[clust_order[3]]))
  
  reg1 <-lm(tmp_ori[cluster==clust_order[1], B] ~ tmp_ori[cluster==clust_order[1], A])
  abline(reg1$coefficients[[1]], reg1$coefficients[[2]], col = color_flag[clust_order[1]])
  
  reg2 <-lm(tmp_ori[cluster==clust_order[2], B] ~ tmp_ori[cluster==clust_order[2], A])
  abline(reg2$coefficients[[1]], reg2$coefficients[[2]], col = color_flag[clust_order[2]])
  
  reg3 <-lm(tmp_ori[cluster==clust_order[3], B] ~ tmp_ori[cluster==clust_order[3], A])
  abline(reg3$coefficients[[1]], reg3$coefficients[[2]], col = color_flag[clust_order[3]])
  
  # slope to angle
  
  angle1 <- atan(reg1$coefficients[[2]]) * 180 /pi
  angle2 <- atan(reg2$coefficients[[2]]) * 180 /pi
  angle3 <- atan(reg3$coefficients[[2]]) * 180 /pi
  
  slope1 <- reg1$coefficients[[2]]
  slope2 <- reg2$coefficients[[2]]
  slope3 <- reg3$coefficients[[2]]
  
  intercept1 <- reg1$coefficients[[1]]
  intercept2 <- reg2$coefficients[[1]]
  intercept3 <- reg3$coefficients[[1]]
  
  intersection12 <- twoLines2intersection(intercept1 = intercept1,intercept2 = intercept2,slope1 = slope1,slope2 = slope2)
  intersection23 <- twoLines2intersection(intercept1 = intercept3,intercept2 = intercept2,slope1 = slope3,slope2 = slope2)
  
  result <- point2AB_percent(newPoint, intercept1 = intercept1, intercept2 = intercept2, intercept3 = intercept3, slope1 = slope1, slope2 = slope2, slope3 = slope3)
  resultStr <- paste("A:",result[1],"% ; B:", result[2],"%")
  
  return(resultStr)
}

# AB_adjust_scatter_plot <- function(dt = AB.dt, pid = pid) {
#   
#   data.m <- matrix(NA, nrow = samples_Num, ncol = 2)
#   # first rotation, let the AB type to horizontal (-45 degree)
#   for (i in 1:samples_Num) {
#     tmp <- AB.dt[probe_id == pid, c(1+get("i"), 1+samples_Num+get("i")), with = FALSE]
#     data.m[i,] <- Rotation(tmp, - probe_frac_slope(pid = pid, frac_lev = 1))
#   }
#   
#   # second, amplify the distance of the points(by magtitude)
#   for (i in 1:samples_Num) {
#     polar_tmp <- cart2polar(data.m[i,1], data.m[i,2])
#     mag <- slope2mag(probe_frac_slope(pid = pid, frac_lev = 2), probe_frac_slope(pid = pid, frac_lev = 0))
#     
#     cart_tmp <- polar2cart(polar_tmp[1], mag*polar_tmp[2])
#     data.m[i,1] <- cart_tmp[1]
#     data.m[i,2] <- cart_tmp[2]
#   }
#   if (mag == 1) cat("Can't calculate mag \n")
#   
#   # Third, re-rotate the AB type to 45 degree
#   data.m <- Rotation(data.m, pi/4)
#   plot(x = data.m[,1], 
#        y = data.m[,2],
#        main = pid,
#        xlab = "A",
#        ylab = "B",
#        asp = 1)
#   # xlim = c(0, 20000),
#   # ylim = c(0, 20000))
# }

### General scatter plot
# input: 
general_scatter_plot_save <- function(newPoint, plot_dt, cluster_dt, plot_f = plot_func(plot_dt, ...), outputDir = "/originDir/", outputFile = ".origin.png") {
  for (i in 1:plot_dt[,.N]) {
    pid <- as.character(plot_dt[i, 1][[1]])
    cat("ID:", pid, "\n")
    png(paste0(plot_outputPath, outputDir, pid, outputFile))
    plot_f(plot_dt, cluster_dt, pid = pid)
    #plot_f(newPoint, pid = pid)
    dev.off()
  }
}



## palette() 可以選顏色