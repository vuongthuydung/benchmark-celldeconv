setwd("path/to/your/data") # change to your folder containing data files

data1 <- read.csv("GSM2230757_human1_umifm_counts.csv", stringsAsFactors = FALSE)
data2 <- read.csv("GSM2230758_human2_umifm_counts.csv", stringsAsFactors = FALSE)
data3 <- read.csv("GSM2230759_human3_umifm_counts.csv", stringsAsFactors = FALSE)
data4 <- read.csv("GSM2230760_human4_umifm_counts.csv", stringsAsFactors = FALSE)

# data <- data[match(data$barcode,unique(data$barcode)),]
data <- rbind(data1,data2,data3,data4)
exprs_data <- data[, -c(1:3)]
exprs_data <- t(exprs_data)
rownames(exprs_data) <- colnames(data)[4:ncol(data)]
colnames(exprs_data) <- paste0("cell_", seq(nrow(data)))
metadata <- paste0("cell_", seq(nrow(data)))
metadata <- cbind(metadata,data[,3])
metadata <- cbind(metadata,c(rep("human1",nrow(data1)),rep("human2",nrow(data2)),rep("human3",nrow(data3)),rep("human4",nrow(data4))))
colnames(metadata) <- c("cellID", "cellType","sampleID")
write.table(metadata, "baron_phenoData.txt", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(as.data.frame(exprs_data), file = "sc_baron.rds")