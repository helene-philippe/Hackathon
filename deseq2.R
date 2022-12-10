#if (!requireNamespace("BiocManager", quietly = TRUE))
#  + install.packages("BiocManager") BiocManager::install("DESeq2")


library("DESeq2")
library("eggCounts")
library("stats")
library("ggplot2")
library("dplyr")


count1<-read.table("SRR628582Counts",stringsAsFactors = TRUE, header=T, sep="\t")
count2<-read.table("SRR628583Counts",stringsAsFactors = TRUE, header=T, sep='\t')
count3<-read.table("SRR628584Counts",stringsAsFactors = TRUE, header=T, sep="\t")
count4<-read.table("SRR628585Counts",stringsAsFactors = TRUE, header=T, sep="\t")
count5<-read.table("SRR628586Counts",stringsAsFactors = TRUE, header=T, sep="\t")
count6<-read.table("SRR628587Counts",stringsAsFactors = TRUE, header=T, sep="\t")
count7<-read.table("SRR628588Counts",stringsAsFactors = TRUE, header=T, sep="\t")
count8<-read.table("SRR628589Counts",stringsAsFactors = TRUE, header=T, sep="\t")


newcount1<-count1[c('Geneid','SRR628582Aligned.sortedByCoord.out.bam')]
newcount2<-count2[c('Geneid','SRR628583Aligned.sortedByCoord.out.bam')]
newcount3<-count3[c('Geneid','SRR628584Aligned.sortedByCoord.out.bam')]
newcount4<-count4[c('Geneid','SRR628585Aligned.sortedByCoord.out.bam')]
newcount5<-count5[c('Geneid','SRR628586Aligned.sortedByCoord.out.bam')]
newcount6<-count6[c('Geneid','SRR628587Aligned.sortedByCoord.out.bam')]
newcount7<-count7[c('Geneid','SRR628588Aligned.sortedByCoord.out.bam')]
newcount8<-count8[c('Geneid','SRR628589Aligned.sortedByCoord.out.bam')]

counts<-merge(newcount1, newcount2, by=c('Geneid'))
counts<-merge(counts, newcount3, by=c('Geneid'))
counts<-merge(counts, newcount4, by=c('Geneid'))
counts<-merge(counts, newcount5, by=c('Geneid'))
counts<-merge(counts, newcount6, by=c('Geneid'))
counts<-merge(counts, newcount7, by=c('Geneid'))
countData<-merge(counts, newcount8, by=c('Geneid'))

geneid<-subset(counts, select=c(1))

condition <- DataFrame(condition=c(rep('Mutant',3), rep('Wildtype',5)))
rownames(condition) = colnames(countData)[2:9]
condition$condition <- as.factor(condition$condition)
condition$condition <- relevel(condition$condition,ref="Wildtype")

dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition, tidy=T)
dds <- DESeq(dds)

results <- results(dds)
summary(results, alpha=0.05)
results <- results[order(results$padj),] # organisé par padj pour avoir les gènes les plus importants

genes_interets = c("ENSG00000115524", "ENSG00000101019", "ENSG00000088256" ,"ENSG00000148848",
                   "ENSG00000245694", "ENSG00000114770", "ENSG00000156052", "ENSG00000163930")

code = c("SF3B1", "UQCC1", "GNA11", "ANKHD1", "CRNDE", "ABCC5", "GNAQ", "BAP1")

results_df = as.data.frame(results)
resultats_interets = results_df[genes_interets,]

png("volcano_all.png")
par(mfrow=c(1,1))
# Make a basic volcano plot
with(results, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-10,10)))
with(subset(results, padj<.05 & log2FoldChange>2), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(results, padj<.05 & log2FoldChange<(-2)), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()

vsdata <- vst(dds, blind=FALSE)
png("PCA_plot.png")
PCA_plot <- plotPCA(vsdata, intgroup="condition")
dev.off()

png("Disp.png")
Disp <- plotDispEsts(dds)
dev.off()

png("volcano_interet.png")
volcano_interet <- ggplot(data=na.omit(resultats_interets), aes(x=log2FoldChange, y=-log10(pvalue), label=code)) + 
  geom_point() + 
  theme_minimal() +
  geom_text() + 
  geom_vline(xintercept=c(-2, 2), col="blue", ) +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  xlim(-10,10)
dev.off()



write.csv(as.data.frame(results[order(results$padj),] ), file="condition_cancer_vs_control_all.csv")
write.csv(resultats_interets[order(results$padj),], file="condition_cancer_vs_control_genes.csv")

