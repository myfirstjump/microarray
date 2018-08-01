

get_allele_fraction_cv <- function(sample){
  
  cat("Chr18 start \n" )
  factor_1 <- allele_fraction_factor_by_data(data_list = sample_connection_target, sample = sample, chr = 18)
  fraction_vec_18 <- factor_1$fraction
  type_vec_18 <- factor_1$type
  
  
  cat("Chr21 start \n" )
  factor_2 <- allele_fraction_factor_by_data(data_list = sample_connection_target, sample = sample, chr = 21)
  fraction_vec_21 <- factor_2$fraction
  type_vec_21_e1 <- factor_2$type
  
  # 計算CV
  allele_fraction.dt <- data.table(frac=c(factor_1$fraction, factor_2$fraction), type=c(factor_1$type, factor_2$type))
  
  cv_table <- allele_fraction.dt[, .(sd = sd(frac)), by = type]
  
  a_allele_vec_with_gap <- c(fraction_vec_18, rep(0, 10), fraction_vec_21)
  type_vec <- c(type_vec_18, rep('BLANK', 10), type_vec_21_e1)
  
  # color_flag
  color_flag <- c("red", "green", "cyan")#palette()[c(1,2,3,5,6,7)]
  # color_flag <- c('BLANK', 'green', 'BLANK')
  
  f <- function(x) {
    switch(x, AA = color_flag[1], AB = color_flag[2], BB = color_flag[3], BLANK = "white") 
  }
  
  color_vec <- sapply(type_vec, f)
  
  num_18 <- length(fraction_vec_18)
  num_21 <- length(fraction_vec_21)
  plot(x = seq(1, num_18+num_21+10),
       y = a_allele_vec_with_gap,
       col = color_vec,
       pch = 19,
       main = sample,
       xlab = "Probe",
       ylab = "Fraction of A allele"
  )
  
  abline(v = num_18, col = "black")
  abline(v = num_18+10, col = "black")
  mtext(side = 3, text = "Chr21", adj = 1)
  mtext(side = 3, text = "Chr18", adj = 0)
  
  return(cv_table)
}
#===========================================================================
allele_fraction_factor_by_data <- function(data_list, sample, chr) {
  # 若40個樣本，於單一probe有超過10個樣本通過threshold，則此probe計算
  #
  #
  
  # Simulation arguments
  # data_list <- sample_connection
  # sample <- 'NA06994_1'
  # chr <- 18
  
  # Attributes
  s_name <- data_list$sample_name
  
  
  if (!(sample %in% data_list$sample_name)) {
    cat("No this sample in the sample_list. \n")
    fraction_vec <- rep(0, 250)
    type_vec <- rep('BLANK', 250)
    comp <- list(fraction = fraction_vec, type = type_vec)
    return(comp)
  }
  
  pid_ordered_vec <- data_chr[order(Probe_Start)][Chromosome==chr][, probe_id]
  # if apply probe pick
  pid_ordered_vec <- intersect(pid_ordered_vec, data_list$filter[Check==1 | Check==2, probe_id])
  
  fraction_vec <- c()
  type_vec <- c()
  
  for (pid in pid_ordered_vec) {
    # pid <- 'PH_NIPT_180002'
    point <- c(as.matrix(data_list$log_substract[probe_id == pid, sample, with = FALSE]))
    type <- c(as.matrix(sample_type[probe_id == pid, sample, with = FALSE]))
    
    probe_pass_threshold_bool <- c(as.matrix(data_list$threshold[probe_id == pid, sample, with = FALSE]))
    
    
    ## 計算時，範圍設定 0~100
    # result <- fraction_of_A_allele_by_pid(pid = pid, point_value = point, data_list = data_list)
    ## 計算時，分開給值0~1
    result <- get_allele_value_type_mapping(pid = pid, pt = point, data_list = data_list, type)
    
    # bool <- result[['check_flag']]
    # if (bool) { #(result[['check_flag']]) {
    #   cat("Probe id:", pid, "\n")
    #   value_check <- c(as.matrix(data_list$log_substract[probe_id == pid, sample_name, with = FALSE]))
    #   cat(value_check, '\n')
    # }
    
    fraction <- result[["fraction"]]
    
    
    too_few_sample_flag <- result[['too_few_sample_flag']]
    if (!too_few_sample_flag & probe_pass_threshold_bool) { # sample & probe條件下，本身threshold為True，且此probe上通過threshold的探針數量 > 10
      if (fraction > 3.001) {
        cat('PID', pid, '\n', 'Origin:', point, '  Mapping:', fraction, 'Type:', type, '\n')
        cat('*******************************\n')
        next # *** 跳過排除所有例外後，mapping值仍大於1的情形
      }
      fraction_vec <- c(fraction_vec, fraction)
      type_vec <- c(type_vec, type)
    } #else { cat("*****", pid, "******", "\n")}
  }
  comp <- list(fraction = fraction_vec, type = type_vec)
  return(comp)
}

#===============================================================================================
### 已知Allele type的情況下，AA的點，在Allele fraction plot中，的散布程度；AB, BB也是一樣
### 若Allele type 不知，此function就沒用意義了
get_allele_value_type_mapping <- function(pid, pt, data_list, type) {
  # pid <- 'PH_NIPT_210263'
  # pt <- -0.1943217
  # data_list <- sample_connection
  
  fraction <- -1
  threshold <- -1
  
  log_sub <- data_list$log_substract
  s_name <- data_list$sample_name
  
  threshold <- c(as.matrix(data_list$threshold[probe_id==pid, s_name, with = FALSE]))
  if(length(threshold[threshold==TRUE]) < 10) {
    # cat("Too few samples meets the threshold.",pid, "\n")
    fraction <- NULL
    too_few_sample_bool <- TRUE # 傳給後面function使用者，告知threshold之後，此probe剩餘sample少於10
    result <- list(fraction=fraction, too_few_sample_flag = too_few_sample_bool)
    return(result)
  } #cat("Too few samples meets the threshold. \n"))
  
  value <- log_sub[probe_id == pid, s_name, with = FALSE]
  value <- value[,j = names(value)[threshold], with = FALSE] # apply threshold
  type <- sample_type[probe_id == pid, s_name, with = FALSE]
  type <- type[,j = names(type)[threshold], with = FALSE] # apply threshold
  type_len <- length(levels(factor(as.matrix(type)))) # how many types
  type_levels <- levels(factor(as.matrix(type)))
  
  common_cols <- intersect(colnames(value), colnames(type))
  pid_dt <- rbind(value[, common_cols, with = FALSE], type[, common_cols, with = FALSE])
  pid_dt <- transpose(pid_dt)
  
  names(pid_dt) <- c("value", "type")
  set(x = pid_dt, j = 1, value = as.numeric(pid_dt[,value]))
  
  ### 與之前判斷是差異之處。
  bound <- pid_dt[, .(median = median(value), sd=sd(value)), by = .(type)]
  
  if (type_len == 3) {
    if (type == 'AA') {
      AA_median <- bound[type == 'AA', median]
      AA_sd <- bound[type == 'AA', sd]
      
    } else if (type == 'AB') {
      AB_median <- bound[type == 'AB', median]
      AB_sd <- bound[type == 'AB', sd]
      
    } else if (type == 'BB') {
      BB_median <- bound[type == 'BB', median]
      BB_sd <- bound[type == 'BB', sd]
      
    }
    
    
    
  }
  
  too_few_sample_bool <- FALSE # 此probe上，通過threshold的probe多於10個
  result <- list(fraction = mapping_value, too_few_sample_flag = too_few_sample_bool)
  return(result)
}
