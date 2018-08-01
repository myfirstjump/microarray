### 需求: NIPT Statistics procedure
### 日期: 180731
### 細節: NIPT 的新Samples進來時，計算所有可能需要的統計量
### 

library(limma)
library("data.table")
library("readxl")
library(xlsx)
library(dplyr)

data_path_gpr_ref <- "//banana/檔案交換區/for EdwardChen/SNP Array base NIPT 20180419/GPR"
data_path_gpr_new <- "//banana/檔案交換區/For Albus/Edward/NIPTG2資料20180727"
data_path_probe_design <- "C:/Users/edwardchen/Desktop/DataSets/NIPT/probe_design"
data_path_sample <- "C:/Users/edwardchen/Desktop/DataSets/NIPT/sample"
output_path <- "C:/Users/edwardchen/Desktop/Mission/201803_NIPT/output"

# ====================================== 與建立Model部分一致 =================================================
# ====================================== 統計量計算在下方    =================================================
# gpr data import
setwd(data_path_gpr_new)
gpr_file <- dir(pattern = "*\\.gpr$", path = getwd())
sample_number_target <- length(gpr_file)
raw_data_target <- read.maimages(gpr_file , source="genepix", columns=list("ID", Rf="F635 Median", Fcv='F635 CV'))

# probe design tables
setwd(data_path_probe_design)
probe_design_file <- dir()
sheet <- excel_sheets(path = "NIPT測試晶片探針.xlsx")
probe_design <- read_xlsx(path = "NIPT測試晶片探針.xlsx", sheet = sheet)

# sample relative file
setwd(data_path_sample)
sample_relative_files <- dir()
sample_genotype <- read.table(file = "1000 genome sample genotype_result.txt", header = T)
sheet <- excel_sheets(path = "gpr_mapping_20180418.xlsx")
sample2gpr <- read_xlsx(path = "gpr_mapping_20180418.xlsx", sheet = sheet[1], col_names = TRUE)
sample2gpr <- sample2gpr[,-1]



### Build data table

# probe_design add column probe_id
col_index <- which(colnames(probe_design) == "rs_id")
colnames(probe_design)[col_index] <- "RS_ID"
probe_design <- as.data.table(probe_design)
probe_design <- probe_design[!grep("B$", probe_design[,`Phalanx ID`])]
probe_design[,"probe_id"] <- gsub("A", "", probe_design[,`Phalanx ID`])

# 截取出探針設計的資訊: RS_ID,       probe_id, allele type, Chromosome, Probe_Start, A_type, B_type
data_chr <- probe_design[,c("RS_ID", "probe_id", "allele type", "Chromosome", "Probe_Start")]
data_chr[, "A_type"] <- apply(data_chr[, "allele type"], 2, function(x) substr(x,1,1))
data_chr[, "B_type"] <- apply(data_chr[, "allele type"], 2, function(x) substr(x,3,3))

# A,B intensity table
A_suffix <- unique(grep("A$", raw_data$genes$ID, perl=TRUE, value = TRUE))
A_suffix <- grep('NIPT', A_suffix, value=TRUE)# length(A_suffix) # 1417
B_suffix <- unique(grep("B$", raw_data$genes$ID, perl=TRUE, value = TRUE))
B_suffix <- grep('NIPT', B_suffix, value=TRUE)# length(B_suffix) # 1417
# Add extra columns to sample_genotype  #加入重複實驗的allele到sample_genotype
sample_name <- c(as.matrix(sample2gpr[, 1]))
tmp <- intersect(colnames(sample_genotype), sample_name)
sample_copy_name <- setdiff(sample_name, tmp) # 2
sample_genotype[, sample_copy_name] <- sample_genotype[, substr(sample_copy_name, start = 1, stop = 7)]
data_chr <- merge(data_chr, sample_genotype[,c("RS_ID", sample_name)])

# allele type in each sample
sample_type <- cbind(data_chr[, c("RS_ID", "probe_id")], apply(data_chr[,sample_name, with = FALSE], 2 , allele_2_ABtype))

### 

A.dt <- build_intensity_table(array = raw_data, suffix = A_suffix, sample2gpr = sample2gpr)
B.dt <- build_intensity_table(array = raw_data, suffix = B_suffix, sample2gpr = sample2gpr)
AB.dt <- merge(A.dt, B.dt, by = "probe_id")
AB_log.dt <- cbind(AB.dt[,1], log2(AB.dt[, -1] + 1))


# filter: 之後畫圖時，probe_id的索引
AB_filtered <- merge(AB.dt[,1], data_chr[,1:2], by = "probe_id")

setwd(data_path_sample)
data_picked <- read_xlsx(path = "點位挑選_20180423.xlsx", col_names = TRUE)
colnames(data_picked) <- c("probe_id", "Check")
AB_filtered <- merge(AB_filtered, data_picked, by = "probe_id")

# threshold table: upper = 10000, lower = 150 boolean table 由AB.dt判斷是否內入分析，捨去訊號過低或過高的情形
# Coriell -> A,B各自的threshold; Real-sample -> A*B的threshold
threshold.dt <- apply_threshold(dt = AB.dt, upper = 10000, lower = 150)
threshold.dt <- apply_threshold_inter(dt = AB.dt, upper = 1e07, lower = 1e04)

threshold.dt <- cbind(AB.dt[,probe_id], threshold.dt)
colnames(threshold.dt)[1] <- "probe_id"
sample_remove <- c()

### table connection: (build without sample_remove)
sample_connection <- list(table = AB.dt, log_table = AB_log.dt, filter = AB_filtered, threshold = threshold.dt, sample_name = setdiff(sample_name, sample_remove))

AB_log_substract.dt <- AB_substract_apply(sample_connection)
sample_connection$log_substract <- AB_log_substract.dt







#==================================== 以下為加入各項欲觀測之數據 ======================================

# Use:
#   raw_data_target
#   sample2gpr_target
# 數據 #1: 5 repetition 的 cv
raw_df <- cbind(raw_data_target$genes, raw_data_target$E)

# 更換column names 例如: gpr 3170709828 -> NA09446之類的
for (i in 1:dim(sample2gpr_target)[1]) { 
  idx <- which(colnames(raw_df) == sample2gpr_target[i,][[2]])
  cat(idx, '\n')
  colnames(raw_df)[idx] <- sample2gpr_target[i, 1]
}

raw_df <- as.data.table(raw_df)
raw_df <- raw_df[, c('ID', sample_name_target), with=FALSE]

# 改cv_df <- raw_df[ID %in% c(A_suffix, B_suffix),]
raw_df <- raw_df[ID %in% c(A_suffix, B_suffix),]
raw_df <- raw_df[order(ID)]

# 加入 threshold
# 一個探針複製十份 因為 A,B各5重複 (目的:用point wise multiply 探針intensity和boolean threshold)
threshold_copy.dt <- sample_connection_target$threshold
for (i in 1:9) { 
  threshold_copy.dt <- rbindlist(list(threshold_copy.dt, sample_connection_target$threshold), use.names=TRUE, fill=TRUE)
}
threshold_copy.dt <- threshold_copy.dt[order(probe_id)]
# 把threshold判斷FALSE的值改為NA
threshold_copy.dt[,2:ncol(threshold_copy.dt)] <- threshold_copy.dt[,2:ncol(threshold_copy.dt)][, lapply(.SD, function(x) replace(x, which(x==FALSE), NA))]

sd_df <- raw_df
sd_df[,2:ncol(sd_df)] <- raw_df[,2:ncol(raw_df)] * threshold_copy.dt[, 2:ncol(threshold_copy.dt)]

sd_df <- sd_df[, lapply(.SD, (function (x) sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)) ), by=ID] # .SD為by以外其他的columns


# 加入統計量
sd_mean <- colMeans(sd_df[, 2:ncol(sd_df)], na.rm = TRUE)
mean_holder_dt <- data.table(rbind(NULL,c(0, sd_mean)))
colnames(mean_holder_dt)[1] <- 'ID'
# *** 5 rep cv mean -> #0

median_holder_dt <- sd_df[,2:ncol(sd_df)][, lapply(.SD, median, na.rm=TRUE)]
median_holder_dt <- cbind(median_holder_dt, list(1))
colnames(median_holder_dt)[which(colnames(median_holder_dt)=='V2')] <- 'ID'
# *** 5 rep cv median -> #1

# rbind 統計量
sd_df <- rbindlist(list(sd_df, mean_holder_dt), use.names=TRUE, fill=TRUE)
sd_df <- rbindlist(list(sd_df, median_holder_dt), use.names=TRUE, fill=TRUE)

sd_df[which(sd_df[,ID]==0), ID := '5rep_mean']
sd_df[which(sd_df[,ID]==1), ID := '5rep_median']




### allele_frac_sd
# 沒有Allele type 則無法計算
# cv_df <- rbindlist(list(cv_df, as.list(c(2, allele_frac_cv))), use.names = FALSE, fill = FALSE)#, idcol = TRUE)
# *** allele fraction cv -> #2


### 平均訊號
# median_df[,2:ncol(median_df)] <- median_df[,2:ncol(median_df)] * threshold_copy.dt[, 2:ncol(threshold_copy.dt)]
# *** 訊號先不使用threshold過濾 ex.全部訊號都進入平均計算，看出平均高低。

avg_sig_df <- raw_df
avg_sig_df <- avg_sig_df[, lapply(.SD, (function (x) median(x, na.rm = TRUE))), by=ID]
# median 取平均
avg_sig_vec <- colMeans(avg_sig_df[, 2:ncol(avg_sig_df)], na.rm = TRUE)
mean_holder_dt <- data.table(rbind(NULL,c(3, avg_sig_vec)))
# *** signal mean -> #3

colnames(mean_holder_dt)[1] <- 'ID'
avg_sig_df <- rbindlist(list(raw_df, mean_holder_dt), use.names=TRUE, fill=TRUE)
avg_sig_df[which(avg_sig_df[,ID]==3), ID := 'Average_signal']



###########
write.xlsx(avg_sig_df, 'C:/Users/edwardchen/Desktop/NIPT_Sample_statistics_0731.xlsx', sheetName = 'raw_data', row.names = FALSE)
write.xlsx(sd_df, 'C:/Users/edwardchen/Desktop/NIPT_Sample_statistics_0731.xlsx', sheetName = '5rep', append = TRUE, row.names = FALSE)
###########






















### pixel 裡面的Variation
tmp_df <- cbind(raw_data$genes, raw_data$Fcv)
for (i in 1:dim(sample2gpr)[1]) { # 更換column names
  idx <- which(colnames(tmp_df) == sample2gpr[i,][[2]])
  cat(idx, '\n')
  colnames(tmp_df)[idx] <- sample2gpr[i, 1]
}
raw_df <- as.data.table(tmp_df)
raw_df <- raw_df[, c('ID', sample_name), with=FALSE]
pixel_cv_df <- raw_df[ID %in% c(A_suffix, B_suffix),]
pixel_cv_df <- pixel_cv_df[order(ID)]
pixel_cv_df[, 2:ncol(pixel_cv_df)] <- pixel_cv_df[, 2:ncol(pixel_cv_df)] * threshold_copy.dt[, 2:ncol(threshold_copy.dt)]
pixel_cv_df <- pixel_cv_df[, lapply(.SD, (function (x) median(x, na.rm = TRUE)) ), by=ID]
### 
pixel_mean <- colMeans(pixel_cv_df[, 2:ncol(pixel_cv_df)], na.rm = TRUE)
holder_dt <- data.table(rbind(NULL,c(4, pixel_mean)))
# *** mean of pixel cv -> #4
colnames(holder_dt)[1] <- 'ID'
cv_df <- rbindlist(list(cv_df, holder_dt), use.names=TRUE, fill=TRUE)


### Target: Allele fraction 的散布程度
allele_frac_cv <- c()
for (s in sample_name) {
  cat(s, '\n')
  tmp <- get_allele_fraction_cv(sample = s)
  allele_frac.dt <- data.table(frac=c(tmp[[1]]$fraction, tmp[[2]]$fraction), type=c(tmp[[1]]$type, tmp[[2]]$type))
  AB_frac_vec <- allele_frac.dt[type=='AB',][, frac]
  allele_frac_cv <- c(allele_frac_cv, sd(AB_frac_vec)/mean(AB_frac_vec))
}

get_allele_fraction_cv <- function(sample){
  
  cat("Chr18 start \n" )
  factor_1 <- allele_fraction_factor_by_data(data_list = sample_connection, sample = sample, chr = 18)
  fraction_vec_18 <- factor_1$fraction
  type_vec_18 <- factor_1$type
  
  
  cat("Chr21 start \n" )
  factor_2 <- allele_fraction_factor_by_data(data_list = sample_connection, sample = sample, chr = 21)
  fraction_vec_21_e1 <- factor_2$fraction
  type_vec_21_e1 <- factor_2$type
  
  a_allele_vec <- c(fraction_vec_18, rep(0, 10), fraction_vec_21_e1)
  type_vec <- c(type_vec_18, rep("BLANK", 10), type_vec_21_e1)
  
  # color_flag
  color_flag <- c("red", "green", "cyan")#palette()[c(1,2,3,5,6,7)]
  
  f <- function(x) {
    switch(x, AA = color_flag[1], AB = color_flag[2], BB = color_flag[3], BLANK = "white") 
  }
  
  color_vec <- sapply(type_vec, f)
  
  num_18 <- length(fraction_vec_18)
  num_21 <- length(fraction_vec_21_e1)
  plot(x = seq(1, num_18+num_21+10),
       y = a_allele_vec,
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
  # return(list(factor_1,factor_2)) # 計算Allele fraction分散程度所需
}






write.xlsx(cv_df, 'C:/Users/edwardchen/Desktop/5_repetitions_cv_df_0702.xlsx', row.names = FALSE)








































sample_name_target <- as.character(sample_name_target)
N_statistics <- 3
sd_table <- matrix(data = rep(0, length(sample_name_target)*N_statistics), nrow = N_statistics, byrow = TRUE)
for (i in 1:length(sample_name_target)) {
  each_sample <- sample_name_target[i]
  cat('Sample', each_sample, 'begin...\n')
  each_cv_table <- get_allele_fraction_cv(each_sample)
  cv_table[1, i] <- each_cv_table[type=='AA', sd]
  cv_table[2, i] <- each_cv_table[type=='AB', sd]
  cv_table[3, i] <- each_cv_table[type=='BB', sd]
}

cv.dt <- as.data.table(cv_table)
row.names(cv.dt) <- c('AA', 'AB', 'BB')
colnames(cv.dt) <- sample_name_target
cv.dt

# get_allele_fraction_cv('NA06994_2')
write.xlsx(cv.dt, 'C:/Users/edwardchen/Desktop/NIPT_allele_frac_SD_table.xlsx')


