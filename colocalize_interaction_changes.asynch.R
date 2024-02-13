###############################
##
## NOTE: Same script as colocalize_interaction_changes.R, but with a new batch
##       of Hi-C data done in asynchronous population
##
##-----------------------------
## Signal co-localization analysis of increased interaction (Hi-C)seen in
## top2KOs vs. control, with all sorts of marks and other features:
##   - Lamina associated domains (LADs)
##   - Replication phase (repliseq)
##   - Chromatin accessibility
##   - R-loops
##   - Cohesin dependant insulation
##   - DSB
##
## Depending on the mark, could be simply colocalization (the typical 2-lines
## plots centered around a genomic feature), or an enrichment analysis compared
## to the regions around or outside.
##
## Split the analysis depending on the gene expression (GRO-seq).
##
## Things to do:
##   1-Compare new asynch pop experiment vs. old with cell cycle not arrested Hi-C analysis: #sites, width, etc.
##   2-Correlation with transcription (GRO-seq)?
##   3-Correlation with DE genes with the absence of TOP2 DKO (transcription 
##     factories: the excel VR sent that I previously did the bidirectional
##     analsysis) <-- Venn + HyperG test
##   4-Correlation with the boundaries of Repli-seq S-phases
##   5-Correlation with LADs (Lamin A/B compartments)
##   6-Correlation with histone marks
##   6-Correlation with enhancers
##
## NOTE: some experiments have 2 replicates, but for the sake of simplicity
## sometimes I use only the first replicate. Ideally, we'll merge (+average) the
## 2 replicates in 1 single track and repeat the analysis.
##
## RUN WITH:
## 
## sbatch --time=5:00:00 --nodes=1 --ntasks=8 --mem=250G -o colocalize_interaction_changes.async.out -J colocalize <<EOF
## #!/bin/bash
## ml R/Bioconductor_3.13_singularity
## Rscript colocalize_interaction_changes.asynch.R
## EOF
##
#################################
library(parallel)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(venneuler)
library(biomaRt)
library(BSgenome.Hsapiens.UCSC.hg38)

PROJECT <- "/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2021_01_meta_HCT116_breakome_prediction"
CORES <- 8
options(mc.cores=CORES)

setwd(PROJECT)
palette("Classic Tableau")   # change R's default palette to Tableau, which is a bit nicer (close to D3 also)
pdf("results/colocalize_interaction_changes.asynch.pdf")

##
## data location
##
tracks <- c(hic_asynch="extdata/increased_interaction_top2KOs.vs.control/Un_DKOWT_250_w10_g2_0.2.bw",    # asynch pop
            hic_na="extdata/increased_interaction_top2KOs.vs.control/intersect_R2_R1_av1.5.hg38.bw", # not arrested
            gro="extdata/groseq/GSE86165_HD1HD2.bw",
            lad="extdata/LADs/LaminB1_antibody.bw",
            repl_s1="extdata/repliseq/S1_sorted.bw",
            repl_s2="extdata/repliseq/S2_sorted.bw",
            repl_s3="extdata/repliseq/S3_sorted.bw",
            repl_s4="extdata/repliseq/S4_sorted.bw",
            repl_s5="extdata/repliseq/S5_sorted.bw",
            repl_s6="extdata/repliseq/S6_sorted.bw",
            repl_s7="extdata/repliseq/S7_sorted.bw",
            repl_s8="extdata/repliseq/S8_sorted.bw",
            repl_s9="extdata/repliseq/S9_sorted.bw",
            repl_s10="extdata/repliseq/S10_sorted.bw",
            repl_s11="extdata/repliseq/S11_sorted.bw",
            repl_s12="extdata/repliseq/S12_sorted.bw",
            repl_s13="extdata/repliseq/S13_sorted.bw",
            repl_s14="extdata/repliseq/S14_sorted.bw",
            repl_s15="extdata/repliseq/S15_sorted.bw",
            repl_s16="extdata/repliseq/S16_sorted.bw",
            atac="extdata/chromatin_accessibility/GSM2719724_Sample1.ATAC-seq-WT-Rep1.peaks.bw",
            rloop="extdata/rloops/imb_roukos_2020_01_5_Top2B_Mock_DRIPc_S5.unique.dupmarked.bw",
            insulation="extdata/cohesin_dependant_insulation/HCT116_cohesin_dependent_insulation_hg38.bw",
            dsb="extdata/dsb/imb_roukos_2020_07_9_longo_hct116_std_dox_1_hct116stddox1_CGTGTGAG_chr-loc-countDifferentUMI_ext.bw",
            H3K4me3="extdata/histone_marks/GSE101966_GSM2719740_H3K4me3_1.bw",
            H3K27ac="extdata/histone_marks/GSE101966_GSM2719748_H3K27ac_1.bw",
            H3K27me3="extdata/histone_marks/GSE101966_GSM2719756_H3K27me3_1.bw",
            H2AZ="extdata/histone_marks/GSE58638_GSM1415873_H2AZ_1.bw",
            H3K4me1="extdata/histone_marks/GSE58638_GSM1415875_H3K4me1_1.bw",
            H3K9me3="extdata/histone_marks/GSE58638_GSM1415878_H3K9me3_1.bw",
            H3K36me3="extdata/histone_marks/GSE58638_GSM1415879_H3K36me3_1.bw")

enhancers <- c(hct116="extdata/enhancers/hct116.bed",
               caco2="extdata/enhancers/caco2.bed",
               colo320="extdata/enhancers/colo320.bed",
               colo829="extdata/enhancers/colo829.bed",
               small_intestine="extdata/enhancers/small_intestine.bed",
#               hct116_starr_r1="extdata/enhancers/hct116_starr_seq.rep1.bed",
               hct116_starr_r2="extdata/enhancers/hct116_starr_seq.rep2.bed")

counts <- c(groseq="extdata/groseq/counts.txt",
            top2ab_dko="extdata/transcription_factory_rnaseq/3.TOP2A_ctrl_vs_TOP2AB_aux_irna_intron_GRCh38_100.csv")

##
##   1-Compare asynch vs. not arrested Hi-C analysis: #sites, width, etc.
##
# import and merge contiguous bins (Hi-C analysis was done in bins of 250k)
# constrains: min gap between regions 1kb, min region size 1kb
hic_asynch <- reduce(import.bw(tracks["hic_asynch"]), drop.empty.ranges=FALSE, min.gapwidth=1000)
hic_asynch <- hic_asynch[width(hic_asynch) > 1000]
hic_na <- reduce(import.bw(tracks["hic_na"]), drop.empty.ranges=FALSE, min.gapwidth=1000)
hic_na <- hic_na[width(hic_na) > 1000]

boxplot(list(Asynch=width(hic_asynch), Not_arrested=width(hic_na)), bty="n", main="Hi-C region width")

genome_len <- sum(seqlengths(hic_asynch))
x <- data.frame(y=c(length(hic_asynch), length(hic_na), sum(width(hic_asynch)) / genome_len, sum(width(hic_na)) / genome_len),
                x=c("Asynch", "Not_arrested", "Asynch", "Not_arrested"),
                what=c(rep("number of regions of enriched interaction", 2),
                       rep("fraction of the genome covered", 2)))
ggplot(x, aes(x=x, y=y, fill=x)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  ggsci::scale_fill_d3() +
  labs(x="", y="") +
  theme_bw() +
  facet_wrap(.~what, scales="free")

##
##   2-Correlation with transcription (GRO-seq)?
##
# read and transform the counts table (get mean log2 RPM for the WT samples)
gro_counts <- read.delim(counts["groseq"])
gro_counts <- do.call(cbind, lapply(gro_counts[, -1], function(x) tapply(x, gro_counts[, 1], sum))) # add counts from duplicated gene names
gro_counts <- gro_counts[, grepl("^HD", colnames(gro_counts))]
gro_counts <- apply(gro_counts, 2, function(x) log2(1 + x * 1e6 / sum(x)))  # log2 RPM
gro_counts <- apply(gro_counts, 1, mean, na.rm=TRUE)

# annotate with coordinates
#ann <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_name", "chromosome_name",
#                            "start_position", "end_position", "strand"),
###             filters    = "ensembl_gene_name",
###             values     = rownames(gro_counts),
#             mart       = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl"))
#write.table(ann, file="ann/gene_ann.txt")
ann <- read.delim("ann/gene_ann.txt")
gro_counts <- data.frame(gene  =names(gro_counts),
                         counts=gro_counts,
                         chr   =ann$Chromosome.scaffold.name[match(names(gro_counts),ann$Gene.name)],
                         start =ann$Gene.start..bp.[match(names(gro_counts),ann$Gene.name)],
                         end   =ann$Gene.end..bp.[match(names(gro_counts),ann$Gene.name)],
                         strand=ann$Strand[match(names(gro_counts),ann$Gene.name)])
gro_counts <- gro_counts[gro_counts$chr %in% c(as.character(1:22), "X", "Y"), ]
gro_counts$chr <- paste0("chr", gro_counts$chr)
gro_counts$strand <- ifelse(gro_counts$strand < 0, "-", "+")
gro_counts <- makeGRangesFromDataFrame(gro_counts, keep.extra.columns=TRUE)

# control plot to show the low/mid/high cutoffs
d <- density(gro_counts$counts)
gro_low_threshold  <- 2  # log2 RPM to consider as gene as not expressed
gro_high_threshold <- d$x[which.max(d$y)] # threshold mid/high expression (it ~5.5 log2 RPM)
gro_low_1000_threshold  <- sort(gro_counts$counts)[1000]
gro_high_1000_threshold <- rev(sort(gro_counts$counts))[1000]

plot(d, main="expression probability density function\nGRO-seq", xlab="")
abline(v=c(gro_low_threshold, gro_high_threshold, gro_low_1000_threshold, gro_high_1000_threshold), col=1:4, lty=2)
legend("topright", fill=1:4, legend=c("low expressed", "mid/high expressed", "bottom 1000", "top 1000"), title="thresholds")

# test if interacting regions are enriched with highly expressed genes
gro_counts$expression <- ifelse(gro_counts$counts < gro_low_threshold , "low",
                         ifelse(gro_counts$counts < gro_high_threshold, "mid","high"))
gro_counts$hic_asynch <- gro_counts %over% hic_asynch
gro_counts$hic_na <- gro_counts %over% hic_na

x <- list(Asynch=table(gro_counts$expression, gro_counts$hic_asynch),
          Not_arrested  =table(gro_counts$expression, gro_counts$hic_na))
p <- chisq.test(x[["Asynch"]])$p.value

df <- reshape2::melt(x)
df$Var1 <- factor(df$Var1, levels=c("low", "mid", "high"))

ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
  geom_bar(stat="identity", position="fill") +
  geom_text(aes(label=as.character(value)), position="fill", size=3) +
  ggsci::scale_fill_d3(name="in region of\nincreased\ninteraction") +
  labs(x="expression", y="", title="Genes on regions of increased interaction",
       subtitle=paste("by expression quantile, p-value", format.pval(p))) +
  theme_bw() +
  facet_wrap(.~L1)

#
# test now from the other side: on regions of increased interaction, how many genes in top/bottom 1000 we have?
#
gro_counts$expression2 <- ifelse(gro_counts$counts <= gro_low_1000_threshold , "bottom 1000",
                          ifelse(gro_counts$counts >= gro_high_1000_threshold, "top 1000","mid"))

x <- with(subset(gro_counts, expression2 %in% c("bottom 1000", "top 1000")), {
  list(Asynch=table(expression2[hic_asynch]),
       Not_arrested  =table(expression2[hic_na]))
})
bg <- length(gro_counts)     # all genes
A  <- sum(gro_counts$expression2 %in% "top 1000")    # 1000
B  <- sum(gro_counts$hic_asynch)                            # genes on regions of increased interaction
AB <- sum(gro_counts$expression2 %in% "top 1000" & gro_counts$hic_asynch)
p  <- sum(dhyper(AB:A, B, bg - B, A))

df <- reshape2::melt(x)

ggplot(df, aes(x=L1, y=value, fill=Var1)) +
  geom_bar(stat="identity", position="fill") +
  geom_text(aes(label=as.character(value)), position="fill", size=3) +
  ggsci::scale_fill_d3(name="expression quantile") +
  labs(x="Hi-C filter", y="", title="Top/bottom 1000 genes on regions of increased interaction",
       subtitle=paste("How many of the top/bottom genes are on regions of increased interaction? p-val", format.pval(p))) +
  theme_bw()

#
# test 3: instead of genes and counts, segment the genome and compare signal (bw)
#
BIN_SIZE <- round(mean(width(hic_asynch)) / 100)   # to have at least, 100 bins of ~10K per region, on average
bins <- tileGenome(seqlengths(BSgenome.Hsapiens.UCSC.hg38), tilewidth=BIN_SIZE, cut.last.tile.in.chrom=TRUE)
bins <- keepSeqlevels(bins, seqlevels(hic_asynch), pruning.mode="coarse")
bins$hic_asynch <- bins %within% hic_asynch   # completely contained within

gro <- import.bw(tracks["gro"])
gro <- keepSeqlevels(gro, seqlevels(hic_asynch), pruning.mode="coarse")
bins <- subsetByOverlaps(bins, gro)

#x <- coverage(gro, weight="score")[bins]
#y <- sapply(x, mean)  # takes 1h!! (and `unlist(mclapply(x, mean))` seems not to work...)
#bins$gro <- y
# use binnedAverage() instead!
x <- coverage(gro, weight="score")
bins <- binnedAverage(bins, x, "gro")

set.seed(666)
df <- as.data.frame(c(bins[bins$hic_asynch], sample(bins[!bins$hic_asynch], sum(bins$hic_asynch))))
df$l2gro <- log2(1 + df$gro)

ggpubr::ggboxplot(df, x="hic_asynch", y="l2gro", color="hic_asynch", legend="none", palette="d3") +
  ggpubr::stat_compare_means(method="t.test") +
  labs(x="region of increased interaction",
       y="average log2 GRO signal per 10KB bin",
       title="GRO signal on regions of increased interaction")

##
##   3-Correlation with DE genes with the absence of TOP2 DKO (transcription 
##     factories: the excel VR sent that I previously did the bidirectional
##     analsysis) <-- Venn + HyperG test
##
tf_rnaseq <- read.csv(counts["top2ab_dko"])
tf_rnaseq$padj[is.na(tf_rnaseq$padj)] <- 1          # fix NAs
tf_rnaseq$log2FoldChange[is.na(tf_rnaseq$log2FoldChange)] <- 0
tf_rnaseq$dge <- tf_rnaseq$padj < .05  # seems that Giuseppe used this threshold

# annotate
tf_rnaseq <- data.frame(gene  =tf_rnaseq$gene_name,
                        dge   =tf_rnaseq$dge,
                        chr   =ann$Chromosome.scaffold.name[match(tf_rnaseq$gene_name,ann$Gene.name)],
                        start =ann$Gene.start..bp.[match(tf_rnaseq$gene_name,ann$Gene.name)],
                        end   =ann$Gene.end..bp.[match(tf_rnaseq$gene_name,ann$Gene.name)],
                        strand=ann$Strand[match(tf_rnaseq$gene_name,ann$Gene.name)])
tf_rnaseq <- tf_rnaseq[tf_rnaseq$chr %in% c(as.character(1:22), "X", "Y"), ]
tf_rnaseq$chr <- paste0("chr", tf_rnaseq$chr)
tf_rnaseq$strand <- ifelse(tf_rnaseq$strand < 0, "-", "+")
tf_rnaseq <- makeGRangesFromDataFrame(tf_rnaseq, keep.extra.columns=TRUE)

# test if there's and overrepresentation of DE genes in increased interaction regions
invisible(
  Map(list(hic_asynch, hic_na), list("Asynch", "Not_arrested"), f=function(hic_asynch, subtitle) {
    tf_rnaseq$hic_asynch <- tf_rnaseq %over% hic_asynch

    bg <- length(tf_rnaseq)     # all genes
    A  <- sum(tf_rnaseq$dge)    # DE genes in the rna-seq experiment
    B  <- sum(tf_rnaseq$hic_asynch)    # genes on regions of increased interaction
    AB <- sum(tf_rnaseq$dge & tf_rnaseq$hic_asynch)  # DE genes on regions of increased interaction
    p  <- sum(dhyper(AB:A, B, bg - B, A))

    vd <- venneuler(reshape2::melt(list(`genes DE in TOP2AB DKO`=tf_rnaseq$gene[tf_rnaseq$dge],
                                        `genes on regions of increased interaction`=tf_rnaseq$gene[tf_rnaseq$hic_asynch])))
    vd$labels <- c("genes DE in TOP2AB DKO", "genes on regions of\nincreased interaction")
    plot(vd, col.txt=NA,
         main="probability of DE gene in TOP2AB DKO\non a region of increased interaction",
         sub=paste(subtitle, "--> p-value", format.pval(p), "(hypergeometric test of overrepresentation)"))
    x <- rbind(vd$centers,
               c(mean(vd$centers[, 1]), 0.4),
               c(0.1, 0.9))
    text(x, col=1, adj=0.5, cex=2/3,
         c(paste0(vd$labels[1], " (", A - AB, " genes)"),
           paste0(vd$labels[2], " (", B - AB, " genes)"),
           paste0("intersection (", AB, " genes)"),
           paste0("all genes tested:", bg)))
  })
)

##
##   4-Correlation with the boundaries of Repli-seq S-phases
##
# tile hic_asynch regions in 100 bins each, and extend them 250K on both ends
EXTEND <- 1e6  # since we normalize the signal to the max(signal), 1e6 works better for
BINS   <- 100  # maximizing differences between the center and the boudaries of the region
bins <- mclapply(tile(hic_asynch, BINS), function(x) {
  region5 <- unlist(tile(GRanges(seqnames(x)[1], IRanges(min(start(x)) - EXTEND, min(start(x)))), BINS))
  region3 <- unlist(tile(GRanges(seqnames(x)[1], IRanges(max(end(x)), max(end(x)) + EXTEND)), BINS))
  #seqinfo(x) <- seqinfo(region5) <- seqinfo(region3) <- Seqinfo(genome="hg38")
  seqinfo(x) <- seqinfo(region5) <- seqinfo(region3) <- seqinfo(hic_asynch)
  c(region5, x, region3)
})

# calculate the normalized average signal per repliseq track around those regions
co <- mclapply(tracks[grepl("^repl_s", names(tracks))], function(f) {
  x <- import.bw(f)
  seqlevels(x, pruning="coarse") <- seqlevels(hic_asynch)
  seqinfo(x) <- seqinfo(hic_asynch)
  x <- coverage(x, weight="score")
  # for each region, calculate the average signal in every bin
  co <- lapply(bins, function(bin) {
    #unlist(lapply(decode(x[trim(bin)]), mean, na.rm=TRUE))  # trim removes out-of-bounds regions
    binnedAverage(trim(bin), x, "a")$a
  })
  # aggregate the signal of all bins in a matrix with dimensions: regions x 300
  co <- do.call(rbind, co)
  co <- apply(co, 2, mean, na.rm=TRUE)
  data.frame(signal=co / max(co), bin=1:length(co))
})

# and plot them together
df <- reshape2::melt(co, id.vars="bin")
df$L1 <- factor(df$L1, levels=unique(df$L1)[order(as.integer(sub("^repl_s", "", unique(df$L1))))])
df$stage <- factor(ifelse(df$L1 %in% c("repl_s1", "repl_s2", "repl_s3"), "early",
                   ifelse(df$L1 %in% c("repl_s4", "repl_s5", "repl_s6", "repl_s7", "repl_s8"), "mid", "late")),
                   levels=c("early", "mid", "late"))

ggplot(df, aes(x=bin, y=value, color=L1)) +
  geom_line() +
  geom_vline(xintercept=3/2*BINS, col="red", lty=2) +
  geom_vline(xintercept=c(BINS, 2*BINS), col="blue", lty=2) +
  scale_color_manual(name="", values=colorRampPalette(c("blue", "red"))(length(levels(df$L1)))) +
  scale_x_continuous(breaks=c(0, BINS, 3/2*BINS, 2*BINS, 3*BINS),
                     labels=c(-EXTEND, "start", "center", "end", EXTEND)) +
  labs(title="co-localization with replication in S phase",
       x="center of the region of increased interaction",
       y="normalized signal") +
  theme_bw() +
  facet_wrap(.~stage, ncol=1)

# what a mess, let's plot 1 candidate for early, mid and late
ggplot(subset(df, L1 %in% c("repl_s1", "repl_s8", "repl_s16")), aes(x=bin, y=value, color=L1)) +
  geom_line() +
  geom_vline(xintercept=3/2*BINS, col="red", lty=2) +
  geom_vline(xintercept=c(BINS, 2*BINS), col="blue", lty=2) +
  scale_color_manual(name="", values=colorRampPalette(c("blue", "red"))(3)) +
  scale_x_continuous(breaks=c(0, BINS, 3/2*BINS, 2*BINS, 3*BINS),
                     labels=c(-EXTEND, "start", "center", "end", EXTEND)) +
  labs(title="co-localization with replication in S phase",
       x="center of the region of increased interaction",
       y="normalized signal") +
  theme_bw()

#
# plot as the overlap of the sites with early middle and late sites
# let's simplifying by looking at the top 10% windows of repliseq signal
#
Map(list(hic_asynch, hic_na), c("Asynch", "Not_arrested"), f=function(hic_asynch, tit) {
  overlap <- reshape2::melt(lapply(tracks[grepl("^repl_s", names(tracks))], function(f) {
    x <- import.bw(f)
    seqlevels(x, pruning="coarse") <- seqlevels(hic_asynch)
    seqinfo(x) <- seqinfo(hic_asynch)
    x <- x[x$score > quantile(ecdf(x$score), .9)]   # take the top 10% of replication peaks
    table(hic_asynch %over% x)
  }))
  overlap$L1 <- factor(overlap$L1, levels=unique(overlap$L1)[order(as.integer(sub("^repl_s", "", unique(overlap$L1))))])

  ggplot(overlap, aes(x=L1, y=value, fill=Var1)) +
    geom_bar(stat="identity", position="fill") +
    geom_text(aes(label=as.character(value)), position="fill", size=3) +
    ggsci::scale_fill_d3(name="on a window with high\nreplication activity") +
    labs(x="", y="", title="Regions of increased interaction on windows with high replication activity",
         subtitle=paste("Overlapping the top 10% windows (50Kb) of replication.", tit, "hic_asynch filters")) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
})

##
##   5-Correlation with LADs (Lamin A/B compartments)
##
# NOTE: binning hic_asynch regions we have done already in the previous chunk (4-Repliseq)
#       This is here just in case we want to use different EXTEND and BINS
# tile hic_asynch regions in 100 bins each, and extend them 1M on both ends
EXTEND <- 1e6  # since we normalize the signal to the max(signal), 1e6 works better for
BINS   <- 100  # maximizing differences between the center and the boudaries of the region
bins <- mclapply(tile(hic_asynch, BINS), function(x) {
  region5 <- unlist(tile(GRanges(seqnames(x)[1], IRanges(min(start(x)) - EXTEND, min(start(x)))), BINS))
  region3 <- unlist(tile(GRanges(seqnames(x)[1], IRanges(max(end(x)), max(end(x)) + EXTEND)), BINS))
  seqinfo(x) <- seqinfo(region5) <- seqinfo(region3) <- seqinfo(hic_asynch)
  c(region5, x, region3)
})

# calculate the normalized average signal around those regions
lad <- import.bw(tracks["lad"])
seqlevels(lad, pruning="coarse") <- seqlevels(hic_asynch)
seqinfo(lad) <- seqinfo(hic_asynch)
x <- coverage(lad, weight="score")
co <- mclapply(bins, function(bin) {
  #unlist(lapply(decode(coverage(lad, weight="score")[trim(bin)]), mean, na.rm=TRUE))  # trim removes out-of-bounds regions
  binnedAverage(trim(bin), x, "a")$a
})
# aggregate the signal of all bins in a matrix with dimensions: regions x 300
co <- do.call(rbind, co)
co <- apply(co, 2, mean, na.rm=TRUE)
co <- data.frame(signal=co, bin=1:length(co))

# and plot the LAD signal centered around the hic_asynch regions
ggplot(co, aes(x=bin, y=signal)) +
  geom_line() +
  geom_smooth() +
  geom_vline(xintercept=3/2*BINS, col="red", lty=2) +
  geom_vline(xintercept=c(BINS, 2*BINS), col="blue", lty=2) +
  scale_x_continuous(breaks=c(0, BINS, 3/2*BINS, 2*BINS, 3*BINS),
                     labels=c(-EXTEND, "start", "center", "end", EXTEND)) +
  labs(title="co-localization with Lamin B1 signal",
       subtitle="negative signal means existence of A-compartment instead of B-compartment",
       x="center of the region of increased interaction",
       y="average signal") +
  theme_bw()

#
# plot as the overlap of the sites with early middle and late sites
# let's simplifying by looking at the top 10% windows of repliseq signal
#
Map(list(hic_asynch, hic_na), c("Asynch", "Not_arrested"), f=function(hic_asynch, tit) {
  overlap <- reshape2::melt(table(hic_asynch %over% lad[lad$score < 0]))
  overlap$compartment <- ifelse(overlap$Var1, "A", "B")

  p <- ggplot(overlap, aes(x=compartment, y=value)) +
    geom_bar(stat="identity", position="dodge") +
    geom_text(aes(label=as.character(value)), size=3) +
    annotate("text", x=2, y=max(overlap$value),
             label=paste0("Chisqr p-val", format.pval(chisq.test(table(hic_asynch %over% lad[lad$score < 0]))$p.value))) +
    labs(x="compartment", y="regions of increased interation",
         title="Regions of increased interaction on LADs",
         subtitle=paste("Overlapping the A- and B-compartments.", tit, "hic_asynch filters")) +
    theme_bw()

  print(p)
})

#
# and yet another approach, this time seen from the other side: hic_asynch signal around A/B compartment transitions
#
REGION <- 1e5   # region to look at to identify Lamin B1 change of sign
EXTEND <- 1e6   # extend the region around Lamin B1 change of sign in both sides

# calculate coverage for every base of the chromosome (including gaps between regions)
x.gaps <- gaps(lad)   # add the regions with no info, and we'll add them a 0 score
x.gaps <- x.gaps[strand(x.gaps) == "*"]  # remove the +/- strand gaps (as thy cover the whole genome)
x.gaps$score <- 0
lad_all <- coverage(c(lad, x.gaps), weight="score")   # coverage of the lad track, per chr, with gaps scored to 0

# identify regions with signal changing sign for REGION bp, and from here create the A-B transition regions
ab_regions <- unlist(GRangesList(mcMap(lad_all, names(lad_all), f=function(co, chr) {  # parallelize per chr
  x  <- decode(co)   # signal per base
  x1 <- c(x[-1], 0)  # signal of the next base
  i  <- which(sign(x) != sign(x1))  # candidate positions
  i  <- i[i > REGION & i < length(x) - REGION]   # remove regions too close to the begining/end to avoid problems down
  i5 <- sapply(i, function(i) {    # identify true positions, where the trend is for at least REGION_kb up/down stream
    A <- mean(x[(i-REGION):i], na.rm=TRUE)
    B <- mean(x[(i+1):(i+REGION)], na.rm=TRUE)
    (A <= 0 & B >= 0)
  })
  i3 <- sapply(i, function(i) {    # identify true positions, where the trend is for at least REGION_kb up/down stream
    A <- mean(x[(i-REGION):i], na.rm=TRUE)
    B <- mean(x[(i+1):(i+REGION)], na.rm=TRUE)
    (A >= 0 & B <= 0)
  })
  tryCatch(GRanges(chr, IRanges(i[i3 | i5] - EXTEND, i[i3 | i5] + EXTEND), strand=ifelse(i5[i3 | i5], "+", "-")),
           error=function(e) GRanges(NULL))
})))

# get rid of regions in other chr that we have no hic_asynch data, and trim regions outside the boundaries
seqlevels(ab_regions, pruning.mode="coarse") <- seqlevels(hic_asynch)
seqinfo(ab_regions) <- seqinfo(hic_asynch)
ab_regions <- trim(ab_regions)
ab_regions <- ab_regions[width(ab_regions) == max(width(ab_regions))]# discard the few shorter regions after trimming, in edges of chrs
ab_regions <- ab_regions[ab_regions %over% hic_asynch]   # purge un-interesting regions that already don't have any hic_asynch signal

# merge regions and recenter
ab_regions <- reduce(ab_regions)
center <- mid(ab_regions)
start(ab_regions) <- center - EXTEND
end(ab_regions)   <- center + EXTEND

# count for the A-B regions, base by base, if it overlaps with a region of increased interaction
hic_asynch_ab <- coverage(hic_asynch)[ab_regions]
hic_asynch_ab <- t(sapply(hic_asynch_ab, decode))   # built a matrix with ab_regions rows by region_with columns, with 0/1 if it overlaps with hic_asynch
hic_asynch_ab[as.character(strand(ab_regions)) == "-", ] <- t(apply(hic_asynch_ab[as.character(strand(ab_regions)) == "-", ], 1, rev))

# and plot, normalized to the highest count
x <- colSums(hic_asynch_ab)
x <- tapply(x, cut(seq_len(length(x)), 1000), mean)
plot(x/max(x), type="l", axes=FALSE, #ylim=c(0,1),
     main="A/B compartment boundaries containing a region of increased interaction",
     xlab="compartment", ylab="overlaps with regions of increased interaction")
abline(v=length(x)/2, col="red", lty=2)
axis(1, c(1, length(x) / 2, length(x)), c(paste0("-", EXTEND), "A <----- + -----> B", EXTEND))
axis(2)

# summarizing all regions in 1 line
fields::image.plot(matrix(x / max(x)), axes=FALSE,
                   col=colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
                   main="A/B compartment boundaries containing a region of increased interaction")
abline(v=0.5, col="green", lty=2, lwd=2)
axis(1, c(0, 0.5, 1), c(paste0("-", EXTEND), "A <----- + -----> B", EXTEND))

##
## yeat another approach (for the plot): one line per hic_asynch_region, centered on ab_region
##
hic_asynch_ab <- hic_asynch[hic_asynch %over% ab_regions]
hic_asynch_ab <- do.call(rbind, mclapply(seq_along(hic_asynch_ab), function(i) {
  x_ab <- coverage(hic_asynch_ab[i])[ab_regions]
  x_ab <- t(sapply(x_ab, decode))   # built a matrix with ab_regions rows by region_with columns, with 0/1 if it overlaps with hic_asynch
  x_ab[as.character(strand(ab_regions)) == "-", ] <- t(apply(x_ab[as.character(strand(ab_regions)) == "-", ], 1, rev))
#  colSums(x_ab)
  x_ab[which.max(apply(x_ab, 1, sum)), ]  # return only the most representative ab_region for the hic_asynch region
}))

# summarize into bins to make the matrix manageable
ncols  <- seq_len(ncol(hic_asynch_ab))
bins   <- cut(ncols, breaks=1000)
hic_asynch_ab <- do.call(cbind, mclapply(split(as.data.frame(t(hic_asynch_ab)), bins), colSums))
hic_asynch_ab_col <- ifelse(hic_asynch_ab > 0, "#FFFF00", "#000000")  # yellow if signal, black if not

# plot the matrix
hic_asynch_ab <- ifelse(hic_asynch_ab > 0, 1, 0)
plot(c(0, 1), c(0, 1), type="n", bty="n", axes=FALSE,
     main="A/B compartment boundaries containing a region of increased interaction",
     xlab="area of A-->B transition", ylab="region of increased interaction")
Axis(side=1, at=c(0, 0.5, 1), labels=c(-EXTEND, "A|B", EXTEND))
h <- hclust(dist(hic_asynch_ab))#, method="manhattan"))
r <- as.raster(hic_asynch_ab_col[h$order,])
rasterImage(r, 0, 0, 1, 1, interpolate=FALSE)
abline(v=0.5, col="red", lty=2)

##
##   6-Correlation with histone marks
##
h3 <- c("H3K4me3", "H3K27ac", "H3K27me3", "H2AZ", "H3K4me1", "H3K9me3", "H3K36me3")

# tile hic_asynch regions in pseudo regions of 3x100 bins
EXTEND <- 1e5  # since we normalize the signal to the max(signal), 1e6 works better for
BINS   <- 100  # maximizing differences between the center and the boudaries of the region
hic_asynch_pseudo <- mclapply(tile(hic_asynch, BINS), function(x) {
  region5 <- unlist(tile(GRanges(seqnames(x)[1], IRanges(min(start(x)) - EXTEND, min(start(x)))), BINS))
  region3 <- unlist(tile(GRanges(seqnames(x)[1], IRanges(max(end(x)), max(end(x)) + EXTEND)), BINS))
  c(region5, x, region3)
})

# calculate the coverage of each mark
hic_asynch_pseudo_h3 <- mclapply(setNames(h3, h3), function(h3) {
  cat(h3, fill=T)   # some feedback, as it runs rather slow
  h3 <- import.bw(tracks[h3])
  seqlevels(h3, pruning="coarse") <- seqlevels(hic_asynch)
  seqinfo(h3) <- seqinfo(hic_asynch)
  x <- coverage(h3, weight="score")
  co <- do.call(rbind, lapply(hic_asynch_pseudo, function(bin) {
    #sapply(decode(coverage(h3, weight="score")[bin]), mean)   # vector of 3x100 bins for region x with the mean coverage
    binnedAverage(bin, x, "a")$a
  }))
  co <- apply(co, 2, sum, na.rm=TRUE)  # aggregate signal per bin
  data.frame(signal=co, bin=seq_len(length(co)))
})
  
# melt data and plot
x <- reshape2::melt(hic_asynch_pseudo_h3, id.vars="bin")
colnames(x) <- c("bin", "var", "value", "h3")

p <- ggplot(x, aes(x=bin, y=value)) +
      geom_line(color="black") +
      geom_smooth(color="red") +
      geom_vline(xintercept=3/2*BINS, col="red", lty=2) +
      geom_vline(xintercept=c(BINS, 2*BINS), col="blue", lty=2) +
      scale_x_continuous(breaks=c(0, BINS, 3/2*BINS, 2*BINS, 3*BINS),
                         labels=c(-EXTEND, "start", "center", "end", EXTEND)) +
      scale_color_manual(name="", values=palette()) +
      labs(title="histone marks within the regions of increased interaction", x="", y="aggregated signal") +
      theme_bw() +
      facet_wrap(~h3, scales="free")

print(p)

#
# zoom in at the start, end and center of the region
# 
hic_asynch$center <- start(hic_asynch) + width(hic_asynch) / 2
WIDTH <- 1e5

p2 <- Map(list(start(hic_asynch), end(hic_asynch), hic_asynch$center),   # position of interest
          c("start", "end", "center"),
          f=function(position, position.name) {

  # build the regions around the position of interest, and calculate the signal
  hic_asynch_pseudo <- GRanges(seqnames(hic_asynch), IRanges(position - WIDTH, position + WIDTH))

  hic_asynch_pseudo_h3 <- mclapply(setNames(h3, h3), function(h3) {
    cat(h3, fill=T)   # some feedback, as it runs rather slow
    h3 <- import.bw(tracks[h3])
    seqlevels(h3, pruning="coarse") <- seqlevels(hic_asynch)
    seqinfo(h3) <- seqinfo(hic_asynch)
    co <- do.call(rbind, lapply(coverage(h3, weight="score")[hic_asynch_pseudo], decode))
    co <- apply(co, 2, sum, na.rm=TRUE)  # aggregate signal per bin
    data.frame(signal=co, bin=seq_len(length(co)))
  })
    
  # melt data and plot
  x <- reshape2::melt(hic_asynch_pseudo_h3, id.vars="bin")
  colnames(x) <- c("bin", "var", "value", "h3")
  BINS <- max(x$bin)

  p <- ggplot(x, aes(x=bin, y=value, color=h3)) +
        geom_line(color="black") +
        geom_smooth(color="red") +
        geom_vline(xintercept=BINS/2, col="red", lty=2) +
        scale_x_continuous(breaks=c(0, BINS/2, BINS),
                           labels=c(as.character(-WIDTH), position.name, as.character(WIDTH))) +
        scale_color_manual(name="", values=palette()) +
        labs(title="histone marks within the regions of increased interaction",
             subtitle=paste(position.name, "of the region of increased interaction"),
             x="", y="aggregated signal") +
        theme_bw() +
        facet_wrap(~h3, scales="free")
  p
})

print(p2[[1]])
print(p2[[2]])
print(p2[[3]])

##
##   7-Correlation with enhancers
##
# tile hic_asynch regions in pseudo regions of 3x100 bins
EXTEND <- 1e5  # since we normalize the signal to the max(signal), 1e6 works better for
BINS   <- 100  # maximizing differences between the center and the boudaries of the region
hic_asynch_pseudo <- mclapply(tile(hic_asynch, BINS), function(x) {
  region5 <- unlist(tile(GRanges(seqnames(x)[1], IRanges(min(start(x)) - EXTEND, min(start(x)))), BINS))
  region3 <- unlist(tile(GRanges(seqnames(x)[1], IRanges(max(end(x)), max(end(x)) + EXTEND)), BINS))
  c(region5, x, region3)
})

# calculate the coverage of each mark
hic_asynch_pseudo_enhancers <- lapply(enhancers, function(enhancer) {
  cat(enhancer, fill=T)   # some feedback, as it runs rather slow
  enhancer <- import.bed(enhancer)
  enhancer$score <- as.numeric(enhancer$name)
  seqlevels(enhancer, pruning.mode="coarse") <- seqlevels(hic_asynch)
  seqinfo(enhancer) <- seqinfo(hic_asynch)
  enhancer <- trim(enhancer)
  co <- do.call(rbind, mclapply(hic_asynch_pseudo, function(x) {
    sapply(decode(coverage(enhancer)[x]), function(x) 1e3 * sum(x) / length(x))   # vector of 3x100 bins for region x with the mean coverage
  }))
  co <- apply(co, 2, sum, na.rm=TRUE)  # aggregate signal per bin
  data.frame(signal=co, bin=seq_len(length(co)))
})
  
# melt data and plot
x <- reshape2::melt(hic_asynch_pseudo_enhancers, id.vars="bin")
colnames(x) <- c("bin", "var", "value", "enhancers")

p <- ggplot(x, aes(x=bin, y=value)) +
      geom_line(color="black") +
      geom_smooth(color="red") +
      geom_vline(xintercept=3/2*BINS, col="red", lty=2) +
      geom_vline(xintercept=c(BINS, 2*BINS), col="blue", lty=2) +
      scale_x_continuous(breaks=c(0, BINS, 3/2*BINS, 2*BINS, 3*BINS),
                         labels=c(-EXTEND, "start", "center", "end", EXTEND)) +
      scale_color_manual(name="", values=palette()) +
      labs(title="Enhancers within the regions of increased interaction", x="", y="aggregated signal (RPK)") +
      theme_bw() +
      facet_wrap(~enhancers, scales="free")

print(p)

dev.off()
