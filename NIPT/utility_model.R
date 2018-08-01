### allele pair to AB type
# A_type: G, B_type:T -> G|G: AA type
# input: data.table with columns A_type, B_type and target sample allele for each rsID
# ex. data_chr[,  ]
# output: sample columns indicates type level (AA, AB, BB)
allele_2_ABtype <- function(x) {
  tmpA <- data_chr[, A_type]
  tmpB <- data_chr[, B_type]
  y <- rep(0, length(x))
  boo1 <- substr(x, 1,1) == tmpA
  boo2 <- substr(x, 3,3) == tmpA
  for (i in 1: length(y)) {
    if (boo1[i] & boo2[i] ) y[i] = "AA"
    else if ( !boo1[i] & !boo2[i] ) y[i] = "BB"
    else y[i] = "AB"
  }
  return(y)
}


### Calculate the median from the same probe_id.
# input: RGlist.
# input: Suffix vector, to determine the target(A or B).
# input: Sample number.
# output: Numeric vector, with suffix length.
probe_Median <- function(RGlist, suffix, sample) {
  tmp = numeric()
  for (i in 1:length(suffix)) {
    # cat(is.null(RGlist$Diff[grep(suffix[i], RGlist$genes$ID), sample], "\n"))
    v1 <- median(RGlist$E[grep(suffix[i], RGlist$genes$ID), sample])
    tmp <- c(tmp, v1)
  }
  return(tmp)
}


### Build intensity table
build_intensity_table <- function(array = RGlist, suffix = A_suffix18, sample2gpr) {
  
  A.dt <- data.table(probe_id = sapply(suffix, function(chr) substr(chr, start = 1, stop = nchar(chr)-1)))
  type_chr <- substr(suffix[1], start = nchar(suffix[1]), stop = nchar(suffix[1]))
  gpr_name <- colnames(array$E)
  sample_number <- dim(sample2gpr)[1]
  sample2gpr_col_name <- colnames(sample2gpr)[2] # should be like 'GPR'
  
  for (i in 1:sample_number) {
    tmp_gpr <- c(as.matrix(sample2gpr[i, sample2gpr_col_name]))
    cat("sample", i,"\n")
    cat(c(as.matrix(sample2gpr[i, 1])), "\n")
    if (tmp_gpr %in% gpr_name) {
      tmp_gpr <- as.character(tmp_gpr)
      cat(tmp_gpr, "\n")
      A.dt$placeholder_name <- probe_Median(array, suffix = suffix, sample=tmp_gpr)
      tmp_sample <- as.character(sample2gpr[i, 1])
      names(A.dt)[names(A.dt) == "placeholder_name"] <- paste0(tmp_sample, paste0(".", type_chr))
    } else cat("No gpr file \n")
  }
  return(A.dt)
}


### Intensity data QC threshold AB值一起篩選 20180611測試規則: A*B > 2e07 or A*B < 5e03 refuse
apply_threshold_inter <- function(dt, upper = 2e07, lower = 5e03) {
  
  row_length <- nrow(dt)
  col_name <- gsub(".A$", "", grep(".A$", names(dt), value = TRUE))
  col_length <- length(col_name)
  
  # build new table (combine A & B)
  new.dt <- as.data.table(matrix(rep(TRUE, row_length*col_length), nrow = row_length))
  colnames(new.dt) <- col_name
  
  for (i in 1:row_length) {
    for (j in 1:col_length) {
      sample <- col_name[j]
      target <- as.matrix(dt[i, c(paste0(sample,".A"), paste0(sample ,".B")), with = FALSE])
      A <- as.numeric(target)[1]
      B <- as.numeric(target)[2]
      
      if ( A*B > upper | A*B < lower) {
        new.dt[i,j] <- FALSE
        cat(sample, '\n')
        cat('A:', target[1], 'B:', target[2], '\n')
      }
      
      # if ( (A > upper & B > upper) | ( A < lower & B < lower ) ) {
      #   new.dt[i,j] <- FALSE
      #   cat(sample, "\n")
      #   cat("A:", target[1],"B:", target[2],"\n")
      # }
    }
  }
  return(new.dt)
}

### get A B indensity substraction table
AB_substract_apply <- function(dt_list) {
  sample_name <- dt_list$sample_name
  log_dt <- dt_list$log_table
  
  dt <- log_dt[, probe_id]
  for (name in sample_name) {  # 可以簡化，用data.table的 := 創造新columns
    tmp_dt <- log_dt[, paste0(name,".A"),with=FALSE] - log_dt[, paste0(name,".B"),with=FALSE]
    dt <- cbind(dt, tmp_dt)
  }
  colnames(dt) <- c("probe_id", sample_name)
  return(dt)  
}


### New sample allele fraction reference from old samples  180426
target_allele_fraction <- function(sample){
  
  cat("Chr18 start \n" )
  factor_1 <- target_allele_fraction_factor_by_ref_data(data_list_target=sample_connection_target, data_list_ref=sample_connection, sample = sample, chr = 18)
  fraction_vec_18 <- factor_1$fraction
  # type_vec_18 <- factor_1$type
  
  cat("Chr21 start \n" )
  factor_2 <- target_allele_fraction_factor_by_ref_data(data_list_target=sample_connection_target, data_list_ref=sample_connection, sample = sample, chr = 21)
  fraction_vec_21_e1 <- factor_2$fraction
  # type_vec_21_e1 <- factor_2$type
  chr_gap <- 10
  a_allele_vec <- c(fraction_vec_18, rep(0, chr_gap), fraction_vec_21_e1)
  # type_vec <- c(type_vec_18, rep("BLANK", 10), type_vec_21_e1)
  
  # color_flag
  # color_flag <- c("black", "BLANK")#palette()[c(1,2,3,5,6,7)]
  # 
  # f <- function(x) {
  #   switch(x, AA = color_flag[1], AB = color_flag[2], BB = color_flag[3], BLANK = "white") 
  # }
  
  
  
  num_18 <- length(fraction_vec_18)
  num_21 <- length(fraction_vec_21_e1)
  color_vec <- c(rep('black', num_18), rep('white', chr_gap), rep('black', num_21))
  plot(x = seq(1, num_18+num_21+chr_gap),
       y = a_allele_vec,
       col = color_vec,
       pch = 19,
       main = sample,
       xlab = "Probe",
       ylab = "Fraction of A allele"
  )
  
  abline(v = num_18, col = 'black')
  abline(v = num_18+chr_gap, col = 'black')
  mtext(side = 3, text = 'Chr21', adj = 1)
  mtext(side = 3, text = 'Chr18', adj = 0)
}



target_allele_fraction_factor_by_ref_data <- function(data_list_target, data_list_ref, sample, chr) {
  
  if (!(sample %in% data_list_target$sample_name)) {
    cat("No this sample in chromosome \n")
    fraction_vec <- rep(0, 250)
    #type_vec <- rep("BLANK", 250)
    comp <- list(fraction = fraction_vec)#, type = type_vec)
    return(comp)
  }
  
  pid_ordered_vec <- c(as.matrix(data_chr[order(Probe_Start)][Chromosome==chr][, .(probe_id)]))
  # pid_ordered_vec <- intersect(pid_ordered_vec, data_list_target$filter[Check==1 & whether_A_minor==0, probe_id])
  # pid_ordered_vec <- intersect(pid_ordered_vec, data_list_target$filter[Check==1 & maternal_fetal_same==FALSE, probe_id])
  
  pid_ordered_vec <- intersect(pid_ordered_vec, data_list_target$filter[Check==1 | Check==2, probe_id])
  # pid_ordered_vec <- intersect(pid_ordered_vec, data_list_target$filter[Check==1, probe_id])
  fraction_vec <- c() # 準備收集 allele fraction 的向量
  
  for (pid in pid_ordered_vec) {
    cat('pid: ', pid, '\n')
    point <- c(as.matrix(data_list_target$log_substract[probe_id == pid, sample, with = FALSE]))  # point 來自 target sample
    result <- fraction_of_A_allele_by_pid(pid = pid, point_value = point, data_list = data_list_ref)  # 根據 ref 計算 target 的 allele fraction
    fraction <- result[["fraction"]]
    fraction_vec <- c(fraction_vec, fraction)
    threshold <- result[["threshold"]]
  }
  comp <- list(fraction = fraction_vec)#, type = type_vec)
  return(comp)
}

### determine allele fraction
fraction_of_A_allele_by_pid <- function(pid, point_value, data_list) {
  # point_value <- 2.55
  # pid <- "PH_NIPT_180422"
  # data_list <- sample_connection
  fraction <- -1
  threshold <- -1
  
  log_sub <- data_list$log_substract
  sample_name <- data_list$sample_name
  
  threshold <- c(as.matrix(data_list$threshold[probe_id==pid, sample_name, with = FALSE]))
  if(length(threshold[threshold==TRUE]) < 10) {
    # cat("Too few samples meets the threshold.",pid, "\n")
    fraction <- NULL
    too_few_sample_bool <- TRUE
    result <- list(fraction=fraction, too_few_sample_flag = too_few_sample_bool)
    return(result)
  } #cat("Too few samples meets the threshold. \n"))
  
  value <- log_sub[probe_id == pid, sample_name, with = FALSE]
  value <- value[,j = names(value)[threshold], with = FALSE] # apply threshold
  type <- sample_type[probe_id == pid, sample_name, with = FALSE]
  type <- type[,j = names(type)[threshold], with = FALSE] # apply threshold
  type_len <- length(levels(factor(as.matrix(type))))
  type_levels <- levels(factor(as.matrix(type)))
  
  common_cols <- intersect(colnames(value), colnames(type))
  pid_dt <- rbind(value[, common_cols, with = FALSE], type[, common_cols, with = FALSE])
  pid_dt <- transpose(pid_dt)
  
  names(pid_dt) <- c("value", "type")
  set(x = pid_dt, j = 1, value = as.numeric(pid_dt[,value]))
  
  bound <- pid_dt[, .(quan.2 = quantile(value, probs = 0.2)[[1]], median = median(value), quan.8 = quantile(value, probs = 0.8)[[1]]), by = .(type)]
  
  if (type_len == 3) { 
    AA <- bound[type == "AA", quan.2]
    AB <- bound[type == "AB", median]
    BB <- bound[type == "BB", quan.8]
    # cat("point:", point_value, "AA:", AA, " AB:", AB, " BB:", BB, "\n")
    if (point_value >= AA) {
      tmp <- 100
    } else if (point_value < AA & point_value >= AB){
      tmp <- 50 + (point_value-AB)/(AA-AB)*50
    } else if (point_value < AB & point_value >= BB) {
      tmp <- (point_value-BB)/(AB-BB)*50
      # tmp <- 0
    } 
    else {
      tmp <- 0
    } 
  } else if (type_len == 2) { ## only 2 cluster
    k1 <- type_levels[1]
    k2 <- type_levels[2]
    # k1,k2 組合為 AA,AB or AB,BB -> check by second char
    base <- ifelse (substr(k1, 2,2) == "A", 50, 0)  
    
    m1 <- ifelse (base == 50, bound[type == k1, quan.2], bound[type == k1, median] )
    m2 <- ifelse (base == 50, bound[type == k2, median], bound[type == k2, quan.8] )
    if (point_value > m1) tmp <- base + 50 
    else if (point_value < m1 & point_value > m2) tmp <- base + (point_value-m2)/(m1-m2)*50
    else tmp <- base
    # tmp <- -200
  } else { ## maybe only 1 cluster **************** This should be modify 判斷這個class是屬於哪一個type
    tmp <- 50
  }
  
  too_few_sample_bool <- FALSE
  result <- list(fraction = tmp, too_few_sample_flag = too_few_sample_bool)
  return(result)
  
}


apply_threshold <- function(dt, upper = 10000, lower = 200) {
  # input: 
  # output:
  row_length <- nrow(dt)
  col_name <- gsub(".A$", "", grep(".A$", names(dt), value = TRUE))
  col_length <- length(col_name)
  
  # build new table (combine A & B)
  new.dt <- as.data.table(matrix(rep(TRUE, row_length*col_length), nrow = row_length))
  colnames(new.dt) <- col_name
  
  for (i in 1:row_length) {
    for (j in 1:col_length) {
      sample <- col_name[j]
      target <- as.matrix(dt[i, c(paste0(sample,".A"), paste0(sample ,".B")), with = FALSE])
      A <- as.numeric(target)[1]
      B <- as.numeric(target)[2]
      
      if ( (A > upper & B > upper) | ( A < lower & B < lower ) ) {
        new.dt[i,j] <- FALSE
        cat(sample, "\n")
        cat("A:", target[1],"B:", target[2],"\n")
      }
    }
  }
  return(new.dt)
}





### scatter plot

# ex. 20180419, 20180426, 20180620 output
scatter_plot_for_pid <- function(pid, file) {
  # pid <- 'PH_NIPT_180002'
  # file <- 
  cat("probe:",pid,"begin\n")
  rsid <- data_chr[probe_id==pid, RS_ID]
  # e1
  # parameters
  plot.dt <- sample_connection$log_table
  sample_name <- sample_connection$sample_name
  filter_vec <- c(as.matrix(sample_connection$threshold[probe_id==pid, -1]))
  if(length(filter_vec[filter_vec==TRUE]) < 10) return(cat("Too few samples meets threshold. \n"))
  A_vec <- c()
  B_vec <- c()
  sample_vec <- c()
  for (i in 1:length(sample_name)) {
    judge <- filter_vec[i]
    if (judge) {
      sample <- sample_name[i]
      a <- c(as.matrix(plot.dt[probe_id == pid, paste0(sample, ".A"), with = FALSE]))
      b <- c(as.matrix(plot.dt[probe_id == pid, paste0(sample, ".B"), with = FALSE]))
      
      sample_vec <- c(sample_vec, sample)
      A_vec <- c(A_vec, a)
      B_vec <- c(B_vec, b)
    }
  }
  allele_type <- c(as.matrix(sample_type[probe_id==pid, sample_vec, with = FALSE]))
  
  # color_flag
  color_flag <- c("red", "green", "cyan")#palette()[c(1,2,3,5,6,7)]
  f <- function(x) {
    switch(x, AA = color_flag[1], AB = color_flag[2], BB = color_flag[3])
  }
  
  col <- sapply(allele_type, f)
  
  # plot
  png(file, res = 100, width = 1080, height = 1080)  ##################################################
  plot(x = A_vec, y = B_vec, pch = 19, col = col, cex = 2, main = paste(rsid, "/", pid), xlab = "A", ylab = "B", xlim=c(5, 13), ylim=c(5,13))
  for (i in 1:length(A_vec)) {
    text(x = A_vec[i], y = B_vec[i], labels = sample_vec[i], cex = 0.8)
  }
  
  # this is for the target sample 180426, 180620
  for (target in sample_name_target) {
    target <- as.character(target)
    judge <- c(as.matrix(sample_connection_target$threshold[probe_id==pid, target, with=FALSE]))
    cat(judge, '\n')
    if (judge) {
      col1 <- paste0(target, '.A')
      col2 <- paste0(target, '.B')
      target_value_1 <- c(as.matrix(AB_log.dt_target[probe_id==pid, col1, with=FALSE]))
      target_value_2 <- c(as.matrix(AB_log.dt_target[probe_id==pid, col2, with=FALSE]))
      # cat(target_value_1,' ',target_value_2, '\n')
      points(x=target_value_1, y=target_value_2, pch=19, cex=3, col='grey')
      text(x=target_value_1, y=target_value_2, labels=target, cex = 0.8)
      # cat('target points on!', target, '\n')
    }
  }
  dev.off()
}
