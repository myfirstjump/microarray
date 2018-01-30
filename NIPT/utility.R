
### Calculate the median from the same probe_id.
# input: RGlist.
# input: Suffix vector, to determine the target(A or B).
# input: Sample number.
# output: Numeric vector, with suffix length.
probeMedian <- function(RGlist, suffix, sample=1) {
  tmp = numeric()
  for (i in 1:length(suffix)) {
    v1 <- median(RGlist$Diff[grep(suffix[i], RGlist$genes$ID), sample])
    tmp <- c(tmp, v1)
  }
  return(tmp)
}


### remove the records(probe_id), mean of AB signal intensity less than 250.
# input: data.table, default is AB.dt.
# input: Remove threshold, default is 250.
# output: A data.table without the specific records.

dt_remove_by_threshold <- function(dt=AB.dt, threshold=250){
  dt$probe_id <- as.character(dt$probe_id)
  reserveList <- c() #list want to be reserved, boolean vector.
  for (row in 1: length(dt$probe_id)) {
    for (j in 1:samples_Num) {
      target <- dt[row, c(paste0("sample.", get("j"),".A"), paste0("sample.", get("j"),".B"))]
      target <- mean(as.numeric(target))
      if (target < threshold) {
        reserveList <- c(reserveList, FALSE)
        break
      } else {
        reserveList <- c(reserveList, TRUE)
        break
      }
    }
  }
  return(dt[reserveList, ])
}




### Calculate AA, AB, BB mean points, and output the slope(group by A,B fraction, in each probe_id).
# input: pid
# input: fraction level, 0(fraction under 0.33) or 1(fraction between 0.33 and 0.66) or 2(fraction over 0.66)
# output: c(mean(A), mean(B))
# probe_frac_slope <- function(pid = "PH_NIPT_210676", frac_lev) {
#   pointList <- list()
#   for (i in 1:length(sampleList)) {
#     tmp <- sampleList[[i]][probe_id==pid & fraction==frac_lev, c(2:3)]
#     pointList[[paste0("sample.", i)]] <- tmp
#   }
#   sum.A <- 0
#   sum.B <- 0
#   cnt <- 0
#   for (i in 1:samples_Num) {
#     if (length(pointList[[i]][[1]]) != 0) {
#       cnt <- cnt + 1
#       sum.A <- sum.A + pointList[[i]][[1]]
#       sum.B <- sum.B + pointList[[i]][[2]]
#     }
#   }
#   mean.A <- sum.A / cnt
#   mean.B <- sum.B / cnt
#   #points(x = mean0.A, y = mean0.B, type = "p", col = "red")
#   if (cnt != 0) {
#     return(atan(mean.B/mean.A))  
#   } else {
#     #at("one-type missing: ", pid, " level:", frac_lev)
#     return(-1)
#   }
# }

### Calculation: Magnitude of the angle should be amplified.
### slope -> degree -> mag
# imput: 2 slope
# output: magnitude
slope2mag <- function(slope1, slope2) {
  if (slope1 == -1 | slope2 == -1) {
    return(1)
  } else {
    deg1 <- slope1
    deg2 <- slope2
    d <- deg1 - deg2
    mag <- (pi/2) / d
    return(mag)
  }
}

### Polar Transformation functions
cart2polar <- function(x,y) {
  r <- (x**2+y**2)**0.5
  theta = atan2(y, x)   # same as atan(y/x)
  return(c(r, theta))
}

polar2cart <- function(r, theta) {
  return(c(r*cos(theta), r*sin(theta)))
}

  
### pheno_fraction_vs_probe_plot (1 sample, 1 graph)
# input: A/B proportion table(data.frame), belongs to [0, 1]
# input: sample number
# output: distribution plot
# pheno_fraction_vs_probe_plot <- function(proportion_df, sample_num) {
#   
# }

### two lines -> intersection point
twoLines2intersection <- function(intercept1, intercept2, slope1, slope2) {
  x = (intercept2 - intercept1) / (slope1 - slope2)
  y = slope1 * x + intercept1
  return(c(x,y))
}

### 判斷點在線的左右
determin_LR_by_pt_and_line <- function(point, lineSlope, lineIntercept) {
  x <- point[1]
  y <- point[2]
  tmp <- y - lineIntercept - lineSlope*x
  if (tmp > 0 ) return("Left")
  else return("Right")
}

twoPoints2slope <- function(point1, point2) {
  x1 <- point1[1]
  y1 <- point1[2]
  x2 <- point2[1]
  y2 <- point2[2]
  m <- (y2-y1) / (x2-x1)
  return(m)
}

point2AB_percent <- function(point, slope1 = slope1, slope2 = slope2, slope3 = slope3, intercept1 = intercept1, intercept2 = intercept2, intercept3 = intercept3) {
  side <- determin_LR_by_pt_and_line(point = point, lineSlope = slope2, lineIntercept = intercept2)
  if (side == "Right") offset <- twoLines2intersection(intercept1, intercept2, slope1, slope2)
  else offset <- twoLines2intersection(intercept3, intercept2, slope3, slope2)
  
  m <- twoPoints2slope(point1 = point, point2 = offset)
  if (side == "Right") {
    A_percent <- 50 + 50*((slope2-m)/(slope2-slope1))
    B_percent <- 50 - 50*((slope2-m)/(slope2-slope1))
  } else {
    A_percent <- 50 - 50*((m-slope2)/(slope3-slope2))
    B_percent <- 50 + 50*((m-slope2)/(slope3-slope2))
  }
  
  if (A_percent > 100) { A_percent = 100; B_percent = 0 }
  if (B_percent > 100) { B_percent = 100; A_percent = 0 }
  
  
  return(c(A_percent, B_percent))
}


#  New sample, judge the A, B percentages
newSampleJudgment <- function(newPoint, pid = "PH_NIPT_210018") {
  tmp <- subset(AB.log2.scale.dt, probe_id == pid, select = c(2:(1+2*samples_Num)))
  tmp_A <- tmp[, 1:(length(tmp)/2)]
  tmp_B <- tmp[, (length(tmp)/2+1):length(tmp)]
  names(tmp_B) <- names(tmp_A) # 為了配合rbind時需要對應的column names
  tmp1 <- rbind(tmp_A, tmp_B)
  tmp1 <- t(tmp1)
  #tmp1 <- tmp1[,1] - tmp1[,2]
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
  
  #plot(tmp_ori[,A], tmp_ori[,B], col=tmp_ori[,cluster], asp = TRUE, xlab = "A", ylab = "B", main = pid, pch = 20)
  #points(newPoint[1], newPoint[2], cex = 5, pch = 3)
  
  
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
  
  # regression line
  reg1 <-lm(tmp_ori[cluster==clust_order[1], B] ~ tmp_ori[cluster==clust_order[1], A])
  #abline(reg1$coefficients[[1]], reg1$coefficients[[2]], col = color_flag[clust_order[1]])
  
  reg2 <-lm(tmp_ori[cluster==clust_order[2], B] ~ tmp_ori[cluster==clust_order[2], A])
  #abline(reg2$coefficients[[1]], reg2$coefficients[[2]], col = color_flag[clust_order[2]])
  
  reg3 <-lm(tmp_ori[cluster==clust_order[3], B] ~ tmp_ori[cluster==clust_order[3], A])
  #abline(reg3$coefficients[[1]], reg3$coefficients[[2]], col = color_flag[clust_order[3]])
  
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
  A_percent <- result[1]
  B_percent <- result[2]
  resultStr <- paste("A:",round(A_percent,2),"% ; B:", round(B_percent,2),"%")
  
  return(c(A_percent, B_percent))
}
  
newSample_probeProp_dist <- function(newPoint = c(6000, 550)) {
  vec <- numeric()
  for (i in 1:AB.dt[,.N]) {
    pid <- as.character(AB.dt[i, 1][[1]])
    vec <- c( vec, newSampleJudgment(newPoint = newPoint, pid = pid)[1])
  }
  hist(vec, breaks = seq(0,100), freq = FALSE, density = 50, col = "red", border = "blue",
       main = paste("intensity A:", newPoint[1], "; intensity B:", newPoint[2]),
       xlab = "A signal proportion")
}


### 根據子佔母之比例，計算AB fraction
fetal_fraction_2_intensity_fraction <- function(fetal_frac, fetal = "AA", maternal = "AA") {
  str1 <- factor(strsplit(fetal, "")[[1]], levels = c("A", "B"))
  str2 <- factor(strsplit(maternal, "")[[1]], levels = c("A", "B"))
  
  tab1 <- table(str1)
  tab2 <- table(str2)
  
  A <- tab1["A"][[1]] * fetal_frac + tab2["A"][[1]] * (1-fetal_frac)
  B <- tab1["B"][[1]] * fetal_frac + tab2["B"][[1]] * (1-fetal_frac)
  
  return( A / (A+B) )
}

### input = A_suf, and spot index; output = A/B ratio list("PHXXXX"=1234, "PHXXXX"=5678)
# probeMedian_ratio <- function(RGlist, A_suf=A_suf, spot=1) {
#   tmp <- list()
#   B_suf <- ATransB(RGlist,A_suf)
#   for (i in 1:length(A_suf)) {
#     v1 = median(RGlist$Diff[grep(A_suf[i], RGlist$genes$ID), spot]) 
#     v2 = median(RGlist$Diff[grep(B_suf[i], RGlist$genes$ID), spot])
#     str = substr(A_suf[i], start = 1, stop = nchar(A_suf[i])-1) # name of list entries
#     tmp[[str]] = v1/(v1+v2)
#   }
#   return(tmp)
# }

###
# getSpotsRatio <- function(RGlist, spots) {
#   result <- list()
#   tmp_str <- c()
#   for (i in spots) {
#     tmp <- probeMedian_ratio(RGlist, A_suf, spot=i)
#     result[[paste("spot", i)]] <- tmp
#   }
#   return(result)
# }


### Produce phenotype numbers of each probe_id (RS_ID).
# input: sample_genotype data.table.
# output: sample_genotype data.table with an extra column typeNum
add_pheno_num_col <- function(sample_genotype.dt = sample_genotype) {
  typeNum <- numeric()
  for (i in 1:sample_genotype.dt[,.N]) {
    vec <- as.matrix(sample_genotype.dt)[i,-1]
    AB_exist_vec <- sapply(vec, function(x) substr(x, 1,1)==substr(x,3,3))
    AB_exist <- FALSE %in% AB_exist_vec
    # AA_BB_exist <- TRUE %in% AB_exist_vec
    # if (AA_BB_exist)  AA_BB_num <- length(table(vec[AB_exist_vec]))
    AA_BB_num <- length(table(vec[AB_exist_vec]))
    num <- ifelse (AB_exist, AA_BB_num+1, AA_BB_num)
    typeNum <- c(typeNum, num)
  }
  dt <- cbind(sample_genotype.dt[,1], typeNum)
  return(dt)
}

