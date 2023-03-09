library(ggplot2)
library(pheatmap)

normCounts_ <- read.csv("normalized_counts.csv", row.names = 1)
deseqRes_ <- read.csv("deseq-results-ordered.csv", row.names = 1)

#add another column to the results to tell if yes or no the padjis<0.05
deseqRes_$sig <- ifelse(deseqRes_$padj <= 0.05, "yes", "no")
#this column is important for ggplot coloring

#visualizations
#remove the outliers(NAs)
deseqRes_ <- na.omit(deseqRes_)
#1. ggplot plotMA
ggplot(deseqRes_, aes(x = log10(baseMean), y = log2FoldChange, color = sig)) + geom_point()

#2. volcano plot
ggplot(deseqRes_, aes(x = log2FoldChange, y = -log10(padj), color = sig)) + geom_point()

#3. heatmaps: only consider differentailly expressed genes: the yes
#but we know the yes are not in the normalized_counts file
#so we need to merge the two files into another
signi <- subset(deseqRes_, padj <= 0.05) #gives total row minus omitted 6824 rows
#merge signi and normcounts_ by rownames
allsig <- merge(normCounts_, signi, by = 0)

sigCounts <- allsig[, 2:127]
row.names(sigCounts) <- allsig$Row.names

pheatmap(log2(sigCounts + 1), scale = 'row')