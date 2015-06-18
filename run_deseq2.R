#

library(DESeq2)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(biomaRt )
library(org.Dm.eg.db)
library(ReportingTools)

# Provide working directory and read files
args <- commandArgs(TRUE)
#directory <- "/csc/analysis/Miguel-Aliaga/JJ/DESeq_Project_JJ/deseq2"
directory <- args[1]
#info<-read.table("/csc/analysis/Miguel-Aliaga/JJ/DESeq_Project_JJ/deseq2/metaData.csv",sep=",",header=T)
info<-read.table(args[2],sep=",",header=T)
sampleFiles=info$filename
sampleCondition=sub("(.gene_counts).*","",sampleFiles)

# Construct deseqdataset
sampleTable<-data.frame(sampleName=sampleCondition,fileName= info$filename,type=info$type,sample= info$sample)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design= ~type)

ddsColl <- ddsHTSeq

pdf("DESeq.pdf")
colData(ddsColl)$type <- factor( colData(ddsColl)$type, levels=c("virgin","mated"))

# Peform Differential expression analysis
ddsColl <- DESeq(ddsColl)

normcds <- round(as.data.frame(counts( ddsColl, normalized=TRUE )))
write.table(normcds, file="normalized.countstable.txt", quote=FALSE, sep="\t", row.names=TRUE)


# Transform data using rlog and vst
rld<-rlog(ddsColl)
vsd<-varianceStabilizingTransformation(ddsColl)

rlogMat<-assay(rld)
vstMat<-assay(vsd)

# Draw heatmaps on raw and transformed count data
select<-order(rowMeans(counts(ddsColl,normalized=TRUE)),decreasing=TRUE)[1:40]
hmcol<-colorRampPalette(brewer.pal(9,"GnBu"))(100)
heatmap.2(counts(ddsColl,normalized=TRUE)[select,],col= hmcol,Rowv=FALSE,Colv=FALSE,scale="none",dendrogram="none",trace="none",margin=c(10,9.5),main=" Heatmap of normalized counts")
heatmap.2(assay(rld)[select,],col= hmcol,Rowv=FALSE,Colv=FALSE,scale="none",dendrogram="none",trace="none",margin=c(10,9.1),main=" Heatmap rlog-transformed data")
heatmap.2(assay(vsd)[select,],col= hmcol,Rowv=FALSE,Colv=FALSE,scale="none",dendrogram="none",trace="none",margin=c(10,9.1),main=" Heatmap vsd-transformed data")

# Draw heatmap to look at sample to sample distances
distsRL<-dist(t(assay(rld)))

mat<-as.matrix(distsRL)
#rownames(mat)<-colnames(mat)<-with(colData(ddsColl),"type")
heatmap.2(mat,trace="none",col=rev(hmcol),margin=c(13,10),main=" Heatmap of sample distances")

# Plot PCA using rlog data
print(plotPCA(rld,intgroup=c("type","sample")))


print(plotPCA(vsd,intgroup=c("type","sample")))

# Differential expression analysis results
res <- results(ddsColl)

plotMA(res,main="MA Plot",ylim=c(-4,4))

## Write result to html using reportingTools--------------------


ensembl = useMart( "ensembl", dataset = "dmelanogaster_gene_ensembl" ) #speicific to drosophilla

report<-function(ensembl,name,out,ddsColl,contrast,factor)
{
res <- results(ddsColl,contrast)
res <- res[order(res$padj),]
genemap <- getBM( attributes = c("ensembl_gene_id", "flybasename_gene","flybasecgid_gene","entrezgene"), filters = "ensembl_gene_id", values = rownames(res), mart = ensembl )
idx <- match( rownames(res), genemap$ensembl_gene_id )
res$name <- genemap$flybasename_gene[ idx ]
res$cgid<-genemap$flybasecgid_gene[ idx ]
res$entrez<-genemap$entrezgene[ idx ]

write.csv(as.data.frame(res),file=paste0(name,".csv"))
results2<-res[!duplicated(res[,"entrez"]),]
results2<-results2[!is.na(results2[,"entrez"]),]
rownames(results2) <- results2$entrez

genemap2 <- getBM( attributes = c("ensembl_gene_id", "entrezgene"), filters = "ensembl_gene_id", values = rownames(ddsColl), mart = ensembl,uniqueRows=TRUE )
genemap2<-genemap2[!is.na(genemap2[,"entrezgene"]),]

idx <- match( rownames(ddsColl), genemap2$ensembl_gene_id )

rownames(ddsColl) = make.names(genemap2$entrezgene[ idx ], unique=TRUE)
trim.leading <- function (x)  sub("^\\X", "", x)
rownames(ddsColl)<-trim.leading(rownames(ddsColl))

des2Report <- HTMLReport(shortName = name, title = name, reportDirectory = out)
 publish(results2,des2Report, pvalueCutoff=0.05, annotation.db="org.Dm.eg.db", DataSet =ddsColl,n = 5000,factor = factor,  reportDir=out,make.plots = TRUE)
 finish(des2Report)

}


report(ensembl,"MatedVsVirgin","MatedVsVirgin",ddsColl,contrast=c("type","mated","virgin"),colData(ddsColl)$type)


#######################################################################################

dev.off()                             
