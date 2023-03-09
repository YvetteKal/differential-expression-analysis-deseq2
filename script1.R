library(DESeq2)
library(apeglm)

#load our normlized count matrix
# Load TSV file
mydata <- read.csv("data/output_counts.csv", header = T, row.names = 1)
info <- read.table("data/condition.txt", header = T, sep = '\t') #rownmes in the file must have the same order than in the count matrix.

#link the two files to DESeq
dds <- DESeqDataSetFromMatrix(mydata, info, ~condition)

#Do not consider lowly expressed genes (lowly count reads) since they are hard to measure
#for differential expression (we cant tell the difference beween one read and another)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#main DESeq
ddsDE <- DESeq(dds)

#extract the normalized read counts produced by DESeq
normCounts <- counts(ddsDE, normalized = T)
write.csv(normCounts, "normalized_counts.csv")

#DESeq results
res <- results(ddsDE, alpha = 0.05) #alpha is the corrected p-value

#look at the summary of our res or summarise the output of the result file: %up and downs regulations
summary(res)

#to see more details
res

#order res by the most significant gene with adjusted p-value (increasing)
#output DESeq results
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, "deseq-results-ordered.csv")

#look at the conditions: the first one is the first showed
resultsNames(ddsDE)
#negative LFC:more expression in second condition and positive: more expression in the first condition

#visualization
#1. plotMA puts the mean expression on the x and the lfc in y
#and the color code gives the significance based on the p-value in blue
pdf("myplot.pdf", width=8, height=6)  # Set the figure width and height to 8 and 6 inches, respectively
plotMA(ddsDE, ylim = c(-5, 5))
dev.off()  # Close the PDF device. a pdf called myplot.pdf has been created



