### 需求: NIPT Allele Fraction Plot
### 日期: 180727
### 細節: 使用四月份建立的Coriell Sample Model，在每一個probe計算新樣本allele fraction的數值
### 

library(limma)
library("data.table")
library("readxl")
#source("C:/Users/edwardchen/Documents/git-workspaces/NIPT/utility.R")

data_path_gpr_ref <- "//banana/檔案交換區/for EdwardChen/SNP Array base NIPT 20180419/GPR"
data_path_gpr_new <- "//banana/檔案交換區/For Albus/Edward/NIPTG2資料20180727"
data_path_probe_design <- "C:/Users/edwardchen/Desktop/DataSets/NIPT/probe_design"
data_path_sample <- "C:/Users/edwardchen/Desktop/DataSets/NIPT/sample"
output_path <- "C:/Users/edwardchen/Desktop/Mission/201803_NIPT/output"

# reference gpr data import
setwd(data_path_gpr_ref)
gpr_file <- dir(pattern = "*\\.gpr$", path = getwd())
sample_number <- length(gpr_file)
raw_data <- read.maimages(gpr_file , source="genepix", columns=list("ID", Rf="F635 Median"))

# probe design tables
setwd(data_path_probe_design)
probe_design_file <- dir()
sheet <- excel_sheets(path = 'NIPT測試晶片探針.xlsx')
probe_design <- read_xlsx(path = 'NIPT測試晶片探針.xlsx', sheet = sheet)

# sample relative file
setwd(data_path_sample)
sample_relative_files <- dir()
sample_genotype <- read.table(file = "1000 genome sample genotype_result.txt", header = T)
sheet <- excel_sheets(path = "gpr_mapping_20180418.xlsx")
sample2gpr <- read_xlsx(path = "gpr_mapping_20180418.xlsx", sheet = sheet[1], col_names = TRUE)
sample2gpr <- sample2gpr[,-1]
#sample2gpr <- sample2gpr[sample2gpr['GPR']!=3170709403,]  # 這個檔案壞了

### Build data table

# probe_design add column probe_id
col_index <- which(colnames(probe_design) == 'rs_id')
colnames(probe_design)[col_index] <- 'RS_ID'
probe_design <- as.data.table(probe_design)
probe_design <- probe_design[!grep("B$", probe_design[,`Phalanx ID`])]
probe_design[,"probe_id"] <- gsub("A", "", probe_design[,`Phalanx ID`])

# 截取出探針設計的資訊: RS_ID,       probe_id, allele type, Chromosome, Probe_Start, A_type, B_type
data_chr <- probe_design[,c("RS_ID", "probe_id", "allele type", "Chromosome", "Probe_Start")]
data_chr[, "A_type"] <- apply(data_chr[, "allele type"], 2, function(x) substr(x,1,1))
data_chr[, "B_type"] <- apply(data_chr[, "allele type"], 2, function(x) substr(x,3,3))

# A,B intensity table
A_suffix <- unique(grep("A$", raw_data$genes$ID, perl=TRUE, value = TRUE)) # length(A_suffix) # 1455
A_suffix <- grep('NIPT', A_suffix, value=TRUE)
B_suffix <- unique(grep("B$", raw_data$genes$ID, perl=TRUE, value = TRUE)) # length(B_suffix) # 1455
B_suffix <- grep('NIPT', B_suffix, value=TRUE)

# Add extra columns to sample_genotype  #加入重複實驗的allele到sample_genotype
sample_name <- c(as.matrix(sample2gpr[, 1]))
tmp <- intersect(colnames(sample_genotype), sample_name)
sample_copy_name <- setdiff(sample_name, tmp) # 2
sample_genotype[, sample_copy_name] <- sample_genotype[, substr(sample_copy_name, start = 1, stop = 7)]
data_chr <- merge(data_chr, sample_genotype[,c("RS_ID", sample_name)])

# allele type in each sample
sample_type <- cbind(data_chr[, c("RS_ID", "probe_id")], apply(data_chr[,sample_name, with = FALSE], 2 , allele_2_ABtype))

A.dt <- build_intensity_table(array = raw_data, suffix = A_suffix, sample2gpr = sample2gpr)
B.dt <- build_intensity_table(array = raw_data, suffix = B_suffix, sample2gpr = sample2gpr)
AB.dt <- merge(A.dt, B.dt, by = "probe_id")
AB_log.dt <- cbind(AB.dt[,1], log2(AB.dt[, -1] + 1))


# filter: 之後畫圖時，probe_id的索引
AB_filtered <- merge(AB.dt[,1], data_chr[,1:2], by = "probe_id")

setwd(data_path_sample)
data_picked <- read_xlsx(path = "點位挑選_20180423.xlsx", col_names = TRUE)
colnames(data_picked) <- c('probe_id', 'Check')
AB_filtered <- merge(AB_filtered, data_picked, by = "probe_id") # 666 3

# threshold table: upper = 10000, lower = 150 boolean table 由AB.dt判斷是否內入分析，捨去訊號過低或過高的情形
threshold.dt <- apply_threshold(dt = AB.dt, upper = 10000, lower = 150)
threshold.dt <- cbind(AB.dt[,probe_id], threshold.dt)
colnames(threshold.dt)[1] <- "probe_id"
sample_remove <- c() # If needs

### table connection: (build without sample_remove)
sample_connection <- list(table = AB.dt, log_table = AB_log.dt, filter = AB_filtered, threshold = threshold.dt, sample_name = setdiff(sample_name, sample_remove))

AB_log_substract.dt <- AB_substract_apply(sample_connection)
sample_connection$log_substract <- AB_log_substract.dt


### Target sample import

# new gpr import
setwd(data_path_gpr_new)
gpr_file <- dir(pattern = "*\\.gpr$", path = getwd())
sample_number_target <- length(gpr_file)
raw_data_target <- read.maimages(gpr_file , source="genepix", columns=list("ID", Rf="F635 Median"))

# sample relative file
setwd(data_path_sample)
sample_relative_files <- dir()
sheet <- excel_sheets(path = "gpr_mapping_20180729.xlsx")
sample2gpr_target <- read_xlsx(path = "gpr_mapping_20180729.xlsx", sheet = sheet[1], col_names = TRUE)

### Build data table

# Add extra columns to sample_genotype  #加入重複實驗的allele到sample_genotype
sample_name_target <- c(as.matrix(sample2gpr_target[, 1]))
# allele type in each sample
sample_type <- cbind(data_chr[, c("RS_ID", "probe_id")], apply(data_chr[,sample_name, with = FALSE], 2 , allele_2_ABtype))

A.dt_target <- build_intensity_table(array = raw_data_target, suffix = A_suffix, sample2gpr = sample2gpr_target)
B.dt_target <- build_intensity_table(array = raw_data_target, suffix = B_suffix, sample2gpr = sample2gpr_target)
AB.dt_target <- merge(A.dt_target, B.dt_target, by = "probe_id")
AB_log.dt_target <- cbind(AB.dt_target[,1], log2(AB.dt_target[, -1] + 1))

# threshold table: upper = 10000, lower = 150 boolean table 由AB.dt判斷是否內入分析，捨去訊號過低或過高的情形
# threshold.dt_target <- apply_threshold_inter(dt = AB.dt_target, upper = 10000, lower = 150)
threshold.dt_target <- apply_threshold_inter(dt = AB.dt_target, upper = 1e07, lower = 1e04)
threshold.dt_target <- cbind(AB.dt_target[,probe_id], threshold.dt_target)
colnames(threshold.dt_target)[1] <- "probe_id"
sample_remove <- c()

### table connection: (build without sample_remove)
sample_connection_target <- list(table = AB.dt_target, log_table = AB_log.dt_target, filter = AB_filtered, threshold = threshold.dt_target, sample_name = setdiff(sample_name_target, sample_remove))
sample_connection_target$sample_name <- as.character(sample_connection_target$sample_name)
AB_log_substract.dt_target <- AB_substract_apply(sample_connection_target)
sample_connection_target$log_substract <- AB_log_substract.dt_target

### scatter plot

general_scatter_plot_by_pid_save <- function(output_dir = "/自取/", outputFile = ".scatter.png") {
  pid_vec <- sample_connection$filter[Check==1 | Check==2 ,probe_id]
  for (i in 1:length(pid_vec)) {
    pid <- pid_vec[i]
    file <- paste0(output_path, output_dir, pid, outputFile)
    scatter_plot_for_pid(file = file, pid = pid)
  }
}

general_scatter_plot_by_pid_save()


### allele plot
target_allele_fraction('3170709828')
target_allele_fraction('3170709829')
target_allele_fraction('3170709830')

general_allele_fraction_plot_by_sample_save <- function(output_dir = "/自取/", outputFile = ".allele_frac.png") {
  sample_name_all <- setdiff(sample_connection_target$sample_name, sample_remove)
  sample_name_all <- '3170709828'
  for (sample in sample_name_all) {
    cat("sample:", sample, "initialize... (function:general_allele_fraction_plot_by_sample_save)\n")
    cat(class(sample), '\n')
    png(paste0(output_path, output_dir, sample, outputFile))
    target_allele_fraction(sample)
    dev.off()
  }
}

general_allele_fraction_plot_by_sample_save()

setwd(output_path)
write.csv(sample_connection_target$table, file='sample_connection_target.csv')


sample_connection_target_bk <- sample_connection_target

sample_connection_target <- sample_connection

target_allele_fraction("NA06994_1")
