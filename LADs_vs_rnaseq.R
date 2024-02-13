###################################
##
## Correlation of LADs with DE genes with the absence of TOP2 DKO (transcription factories)
##
###################################
library(parallel)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg38)
library(gridExtra)

PROJECT <- "/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2021_01_meta_HCT116_breakome_prediction"
CORES <- 8
options(mc.cores=CORES)

setwd(PROJECT)
pdf("results/LADs_vs_rnaseq.pdf")

##
## data location
##
tracks <- c(lad="extdata/LADs/LaminB1_antibody.bw")
counts <- c(top2ab_dko="extdata/transcription_factory_rnaseq/3.TOP2A_ctrl_vs_TOP2AB_aux_irna_intron_GRCh38_100.csv")

# annotate with coordinates
#ann <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_name", "chromosome_name",
#                            "start_position", "end_position", "strand"),
###             filters    = "ensembl_gene_name",
###             values     = rownames(gro_counts),
#             mart       = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
#write.table(ann, file="ann/gene_ann.txt")
ann <- read.delim("ann/gene_ann.txt")

##
##  Correlation of LADs with DE genes with the absence of TOP2 DKO (transcription factories)
##
lad <- import.bw(tracks["lad"])

tf_rnaseq <- read.csv(counts["top2ab_dko"])
tf_rnaseq$padj[is.na(tf_rnaseq$padj)] <- 1          # fix NAs
tf_rnaseq$log2FoldChange[is.na(tf_rnaseq$log2FoldChange)] <- 0
tf_rnaseq$dge <- tf_rnaseq$padj < .05  # seems that Giuseppe used this threshold

# annotate
tf_rnaseq <- data.frame(gene  =tf_rnaseq$gene_name,
                        dge   =tf_rnaseq$dge,
                        log2fc=tf_rnaseq$log2FoldChange,
                        chr   =ann$Chromosome.scaffold.name[match(tf_rnaseq$gene_name,ann$Gene.name)],
                        start =ann$Gene.start..bp.[match(tf_rnaseq$gene_name,ann$Gene.name)],
                        end   =ann$Gene.end..bp.[match(tf_rnaseq$gene_name,ann$Gene.name)],
                        strand=ann$Strand[match(tf_rnaseq$gene_name,ann$Gene.name)])
tf_rnaseq <- tf_rnaseq[tf_rnaseq$chr %in% c(as.character(1:22), "X", "Y"), ]
tf_rnaseq$chr <- paste0("chr", tf_rnaseq$chr)
tf_rnaseq$strand <- ifelse(tf_rnaseq$strand < 0, "-", "+")
tf_rnaseq <- makeGRangesFromDataFrame(tf_rnaseq, keep.extra.columns=TRUE)

# calculate average LAD signal along the gene body
x <- decode(coverage(lad, weight="score")[tf_rnaseq])
tf_rnaseq$avg_lad_signal <- unlist(mclapply(x, mean, na.rm=TRUE))
tf_rnaseq$lam_B <- unlist(mclapply(x, function(x) { sum(x > 0) / length(x) }))
tf_rnaseq$lam_A <- 1-tf_rnaseq$lam_B

# test if average of FC in DEG is different in Lamin A (positive signal) or Lamin B (negative signal)
x <- ifelse(tf_rnaseq$avg_lad_signal[tf_rnaseq$dge] > 0, "B", "A")
y <- tf_rnaseq$log2fc[tf_rnaseq$dge]
p <- t.test(y ~ x)$p.value
boxplot(y ~ x, ylab="logFC", xlab="compartment", main="logFC of DEG, depending on the Lamina compartment")
points(y ~ jitter(as.numeric(as.factor(x))), pch=16, col="#000000AA")
legend("topright", legend=paste0("pval=", format.pval(p)))

# test if ratio of DGE is different in Lamin A (positive signal) or Lamin B (negative signal)
x <- table(dge=tf_rnaseq$dge, lad=ifelse(tf_rnaseq$avg_lad_signal > 0, "B", "A"))
f <- fisher.test(x)
x <- reshape2::melt(x)

ggplot(x, aes(x=lad, y=value, fill=dge)) +
  geom_bar(stat="identity", position="stack") +
  labs(title="ratio of DEG depending on the Lamina compartment",
       subtitle=paste0("pval=", format.pval(f$p.value)), x="compartment", y="") +
  scale_fill_manual("Diff. expressed", values=palette()) +
  theme_bw()

# test if ratio of up/down DGE is different in Lamin A (positive signal) or Lamin B (negative signal)
x <- table(dge=ifelse(tf_rnaseq$log2fc[tf_rnaseq$dge] > 0, "Up", "down"),
           lad=ifelse(tf_rnaseq$avg_lad_signal[tf_rnaseq$dge] > 0, "B", "A"))
f <- fisher.test(x)
x <- reshape2::melt(x)

p1 <- ggplot(x, aes(x=lad, y=value, fill=dge)) +
        geom_bar(stat="identity", position="stack") +
        labs(title="ratio of up/down DEG depending on the Lamina compartment",
             subtitle=paste0("pval=", format.pval(f$p.value)), x="compartment", y="") +
        scale_fill_manual("Diff. expressed", values=palette()) +
        theme_bw() +
        theme(legend.position="bottom")

p2 <- ggplot(x, aes(x=lad, y=value, fill=dge)) +
        geom_bar(stat="identity", position="fill") +
        labs(title="ratio of up/down DEG depending on the Lamina compartment",
             subtitle=paste0("pval=", format.pval(f$p.value)), x="compartment", y="") +
        scale_fill_manual("Diff. expressed", values=palette()) +
        theme_bw() +
        theme(legend.position="bottom")

grid.arrange(p1, p2, ncol=2)

dev.off()
