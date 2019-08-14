library(corrplot)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)

base_path <- "/Users/kristen/Documents/sshfs_dir/"

# Load in count table and name based on sample
count_file <-
  "/Users/kristen/Documents/sshfs_dir/analysis_outs/accreta_countTable.txt"
count_table <- read.table(count_file, sep = "\t", row.names = 1, header = TRUE)

sample_file <- paste0(base_path, "sample_info.csv")

sample_table <- read.table(sample_file, sep = ",", header = TRUE,
                           row.names = 1)

count_table_names <- sample_table$Group.Name.B
names(count_table_names) <- rownames(sample_table)
new_colnames <- count_table_names[colnames(count_table)]

count_table_groupb <- count_table
colnames(count_table_groupb) <- new_colnames

# Correlations between all samples
pdf(paste0(base_path, "analysis_outs/images/correlation.pdf"))

correlation_matrix_groupb <- cor(count_table_groupb)

corrplot(as.matrix(correlation_matrix_groupb),
         method = "circle", order = "hclust")

correlation_matrix <- cor(count_table)

corrplot(as.matrix(correlation_matrix),
         method = "circle", order = "hclust")

dev.off()

count_table <- count_table[, rownames(sample_table)]

dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = sample_table,
                              design = ~Group.Name.B)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Group.Name.B, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup="Group.Name.B")
