


pid = "PH_NIPT_210014"
newPoint = c(6000, 3000)
cluster_dt = AB.log2.scale.dt

tmp <- subset(cluster_dt, probe_id == pid, select = c(2:(1+2*samples_Num)))
tmp_A <- tmp[, 1:(length(tmp)/2)]
tmp_B <- tmp[, (length(tmp)/2+1):length(tmp)]
names(tmp_B) <- names(tmp_A) # 為了配合rbind時需要對應的column names
tmp1 <- rbind(tmp_A, tmp_B)
tmp1 <- t(tmp1)
tmp_modeling <- cbind(tmp1, (tmp1[,1] - tmp1[,2]))
# tmp_modeling <- as.data.frame(tmp1[,1] - tmp1[,2])
colnames(tmp_modeling) <- c("A", "B", "ratio")
# colnames(tmp_modeling) <- c("ratio")

#----
E.dist <- dist(tmp_modeling, method="euclidean")
h.E.cluster <- hclust(E.dist)
plot(h.E.cluster, xlab="歐式距離")
cut.h.cluster <- cutree(h.E.cluster, k=3)
#----

#str(tmp1)
#md <- Mclust(tmp_modeling, G = 3)
#dist_vec <- c(dist(t(md$parameters$mean)))
#if (max(dist_vec)/min(dist_vec) > 3.6) md <- Mclust(tmp1, G = 2)
#cat("group #:", md$G, "\n")

tmp1 <- cbind(tmp1, md$classification)
tmp1 <- cbind(tmp1, cut.h.cluster)

colnames(tmp1) <- c("A", "B", "cluster")
cluster_index <- ncol(tmp1)





#------

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

setwd(true_label_Dir)
genotype_result <- read.table("1000 genome sample genotype_result.txt", header = TRUE)
genotype_result


sample_df[,1]
colnames(genotype_result)



