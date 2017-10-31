library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)


rna <- read_tsv("RawData/GSE80744_ath1001_tx_norm_2016-04-21-UQ_gNorm_normCounts_k4.tsv.gz")
accessions <- read_tsv("RawData/accession_sra.tsv", col_names = c("SRA", "ID", "accession"))

# need to add group info
accessions$group <- c('Low', 'Low', 'Low', 'High', 'High', 'High', 'High', 'Low', 'Low',
                      'High', 'Low', 'High', 'Low', 'Low', 'Low', 'High')

colnames(rna) <- gsub("X", "", colnames(rna))
rna <- as.data.frame(rna)
rownames(rna) <- rna$gene_id

have_accessions <- accessions$ID %in% colnames(rna)
accessions <- accessions[have_accessions,]

gene <- rna["AT1G44575" , as.character(accessions$ID)]
gene <- gather(gene, ID, expression, 1:ncol(gene))
gene$ID <- as.numeric(gene$ID)

dat <- left_join(accessions, gene)

p <- ggplot(dat, aes(group, expression)) +
    geom_jitter(width = 0.1) +
    theme_bw() +
    ggtitle("NPQ4")

ggsave("Plots/rna.png", p)

pdf("Plots/npq4_hist.pdf")
hist(as.numeric(rna['AT1G44575',]), breaks = 50,
     col = 'lightblue', xlab = "Normalized counts", main = "AT1G44575 Expression 1001G")
dev.off()
