### NIPT project

### packages import
library(limma)
library(xlsx)
library(mclust)
library("data.table")
## library("spdep")
library("readxl")
source("utility.R")
source("plot_utility.R")

gprFileDir1 <- "C:/Users/edwardchen/Desktop/DataSets/Cy3 GPR/SNP Array base NIPT 20180108/gpr/"
gprFileDir2 <- "C:/Users/edwardchen/Desktop/DataSets/Cy3 GPR/SNP Array base NIPT 20180109/Gpr/"
Probe2rsID_Dir <- "C:/Users/edwardchen/Desktop/Resource"
probe_pick_Dir <- "C:/Users/edwardchen/Desktop/DataSets/Cy3 GPR"
plot_outputPath <- "C:/Users/edwardchen/Desktop/Mission/4. NIPD/output"
txt_outputPath2 <- "C:\\Users\\edwardchen\\Documents\\git-workspaces\\NIPT\\"
true_label_Dir <- "C:\\Users\\edwardchen\\Desktop\\Mission\\4. NIPD"

# gprFileDir1 <- "C:\\Users\\CJSCOPE\\Desktop\\pbg\\NIPT\\dataset\\SNP Array base NIPT 20180108\\gpr"
# gprFileDir2 <- "C:\\Users\\CJSCOPE\\Desktop\\pbg\\NIPT\\dataset\\SNP Array base NIPT 20180109\\Gpr"
# Probe2rsID_Dir <- "C:\\Users\\CJSCOPE\\Desktop\\pbg\\NIPT\\dataset"
# probe_pick_Dir <- "C:\\Users\\CJSCOPE\\Desktop\\pbg\\NIPT\\dataset"
# plot_outputPath <- "C:/Users/edwardchen/Desktop/Mission/4. NIPD/output"
# txt_outputPath2 <- "C:\\Users\\edwardchen\\Documents\\git-workspaces\\NIPT\\"

### Data import
# gpr files
setwd(gprFileDir1)
files1 <- dir(pattern="*\\.gpr$", full.names = T, path = getwd())
setwd(gprFileDir2)
files2 <- dir(pattern="*\\.gpr$", full.names = T, path = getwd())
files <- c(files1,files2)

# smaple list files
setwd(gprFileDir1)
sample_list_file1 <- dir(path = "..//", pattern = "\\.xlsx$")
sample_sheets_name1 <- excel_sheets(path = paste0("..//", sample_list_file1))
sample_df1 <- read_xlsx(path = paste0("..//", sample_list_file1), sheet = sample_sheets_name1)

setwd(gprFileDir2)
sample_list_file2 <- dir(path = "..//", pattern = "\\.xlsx$")
sample_sheets_name2 <- excel_sheets(path = paste0("..//", sample_list_file2))
sample_df2 <- read_xlsx(path = paste0("..//", sample_list_file2), sheet = sample_sheets_name2[1])

sample_df <- rbind(sample_df1, sample_df2)
colnames(sample_df)[1] <- "SNP_sample"

# True label file
setwd(true_label_Dir)
genotype_result <- read.table("1000 genome sample genotype_result.txt", header = TRUE)
genotype_result <- as.data.table(genotype_result)

# take out the genotype of the samples
geno_index <- match(sample_df[,1], colnames(genotype_result))
geno_id <- c(na.omit(colnames(genotype_result)[geno_index]))
sample_genotype <- subset(genotype_result, select = c("RS_ID" ,geno_id))
rsID2typeNum.dt <- add_pheno_num_col(sample_genotype.dt = sample_genotype)
colnames(rsID2typeNum.dt) <- c("rsID", "typeNum")

# sample amounts
samples_Num <- length(files)

# read gpr file
RGlist <- read.maimages(files , source="genepix", columns=list("ID", Rf="F635 Median", Rb="B635 Median", Diff="F635 Median - B635"))

# probeID <-> rsID 轉換表
setwd(Probe2rsID_Dir)
rsID2Probe_df <- read.table("Probe2rsID.txt", header = TRUE)

# 點位篩選
setwd(probe_pick_Dir)
sheets_name <- excel_sheets(path = "NIPT SNP挑選20180117.xlsx")
data_pick <- read_xlsx(path = "NIPT SNP挑選20180117.xlsx", sheet = sheets_name[2])
colnames(data_pick) <- c("rsID", "Index")
data_pick <- merge(data_pick, rsID2Probe_df, by = "rsID")
colnames(data_pick)[3] <- "probe_id"
data_pick.dt <- as.data.table(data_pick)
rm(data_pick)



### Build A,B data.table
A_suffix <- unique(grep("A$", RGlist$genes$ID, perl=TRUE, value = TRUE)) # length(A_suffix) # 751
B_suffix <- unique(grep("B$", RGlist$genes$ID, perl=TRUE, value = TRUE)) # length(B_suffix) # 750

A.dt <- data.table(probe_id = sapply(A_suffix, function(chr) substr(chr, start = 1, stop = nchar(chr)-1)))
for (i in 1:length(files)) {
  A.dt$placeholder_name <- probeMedian(RGlist, suffix = A_suffix, sample=i)
  names(A.dt)[names(A.dt) == "placeholder_name"] <- paste0("sample.", i, ".A")
}
B.dt <- data.table(probe_id = sapply(B_suffix, function(chr) substr(chr, start = 1, stop = nchar(chr)-1)))
for (i in 1:length(files)) {
  B.dt$placeholder_name <- probeMedian(RGlist, B_suffix, sample=i)
  names(B.dt)[names(B.dt) == "placeholder_name"] <- paste0("sample.", i, ".B")
}

# AB table
AB.dt <- merge(A.dt, B.dt, by = "probe_id")  # dim(AB.dt) [1] 750 117
AB.dt_bk <- AB.dt
AB.dt <- dt_remove_by_threshold(AB.dt, 250)  # dim(AB.dt) [1] 428 117

# AB ratio table
# ratio.dt <- data.table(
#   probe_id = AB.dt$probe_id
# )
# for (i in 1:samples_Num) {
#   c1 <- AB.dt[[paste0("sample.",get("i"),".A")]]
#   c2 <- AB.dt[[paste0("sample.",get("i"),".B")]]
#   ratio.dt$placeholder_name <-  c1 / (c1 + c2)
#   names(ratio.dt)[names(ratio.dt) == "placeholder_name"] <- paste0("sample.", get("i"), ".ratio")
# }
# 
# AB.dt <- merge(AB.dt, ratio.dt, by = "probe_id")
dim(AB.dt) # [1] 428 105
AB.dt <- merge(AB.dt, data_pick.dt, by = "probe_id")
dim(AB.dt) # [1] 315 107 加入 rsID and Index(=1 or 2)
AB.l1.dt <- AB.dt[Index==1, ]
dim(AB.l1.dt) # [1] 265 107
AB.dt <- AB.l1.dt
AB.dt <- merge(AB.dt, rsID2typeNum.dt, by = "rsID")[, union(names(AB.dt), names(rsID2typeNum.dt)), with = FALSE]

# adjusted table (raw data tranformation)
AB.log2.dt <- log2(1 + subset(AB.dt, select = c(2:(samples_Num*2+1))))
AB.log2.dt <- cbind( AB.dt[, .(probe_id)], AB.log2.dt)
AB.log2.scale.dt <- cbind( AB.dt[, .(probe_id)], apply(AB.log2.dt, 2, scale))
AB.log2.scale.dt <- cbind( AB.log2.scale.dt, AB.dt[, .(typeNum)])

# save available probe ID
#write.table(AB.dt[,1], paste0(txt_outputPath2, "availableProbeID.txt"), sep="\t")

# For each sample, assign to a new data.table, wrapped by a list.
# sampleList <- list()
# for (i in 1:samples_Num) {
#   name <- paste0("sample.", i, ".dt")
#   cols <- c(1, 1+i, 1+i+samples_Num, 1+i+samples_Num*2)
#   tmp <- subset(AB.dt, select = cols)
#   tmp[, group_1 := tmp[, 4]<1/3][, group_2 := tmp[, 4]<2/3]
#   tmp[, fraction := tmp[,5]+tmp[,6]]
#   sampleList[[paste0("sample.", i, ".dt")]] <- tmp[, c(1:4, 7)]
# }


### AB_scatter_plot (1 probe, 1 graph)
AB_scatter_plot(plot_dt = AB.dt, pid = "PH_NIPT_210676")

### Calculate A, B allele proportion
AB_scatter_plot_color_by_clustering(cluster_dt = AB.log2.scale.dt, plot_dt = AB.dt, pid = "PH_NIPT_210014")
# AB_scatter_plot_color_by_clustering(cluster_dt = AB.log2.scale.dt, plot_dt = AB.log2.scale.dt, pid = "PH_NIPT_210224")

### AB_scatter_plot_with_color_by_Mcluster
# general_scatter_plot_save(newPoint = c(500, 5000), plot_dt = AB.dt, cluster_dt = AB.log2.scale.dt, plot_f = newSampleJudgment_plot, outputDir = "/RegressionLine/", outputFile = ".png")
# general_scatter_plot_save(plot_dt = AB.log2.scale.dt, cluster_dt = AB.log2.scale.dt, plot_f = AB_scatter_plot_color_by_clustering, outputDir = "/kmeans_origin~A-B/", outputFile = ".png")
general_scatter_plot_save(plot_dt = AB.log2.scale.dt, cluster_dt = AB.log2.scale.dt, plot_f = AB_scatter_plot_color_by_clustering, outputDir = "/hcluster_transfrom(log.scale.complete)~A,B,A-B/", outputFile = ".png")
general_scatter_plot_save(plot_dt = AB.log2.scale.dt, cluster_dt = AB.log2.scale.dt, plot_f = AB_scatter_plot_color_by_clustering, outputDir = "/Mcluster_transform(log.scale)~A-B/", outputFile = ".png")
### 
newSampleJudgment(newPoint = c(6000, 6000), pid = "PH_NIPT_210014")
newSampleJudgment_plot(newPoint = c(6000, 6000), pid = "PH_NIPT_210014")

### 
newSample_probeProp_dist(newPoint = c(5400, 5500))




