####################################
##
## Investigate the role of supercoiling at the LAD borders and MicroC loops
##
##   1- Does supercoiling accumulate at LAD borders?
##   2- Separate MicroC loops into small/medium/large and answer if:
##      2.1- Supercoiling accumulates in the contacts?
##      2.2- the number of supercoiling peaks within the loop differ?
##      2.3- the length of the supercoiling peaks within the loop differ, depending on the size of the loop?
##   3- Classifiy supercoiling at the loop contacts taking into account the
##      directionality of the transcription (TTchem-seq) and whether the contacts have CTCF+Rad21
##   4- How many loop contacts overlap a Rad21+CTCF peak? <-- determine anchored loop contacts
##   5- Now, look at the loops not anchored (!Rad21 !CTCF) if they have supercoiling or not
##   6- Same as 5- but with a different categorization suggested by Akis:
##      a-CTCF binding at both anchors without TSS and H3K27ac
##        (these are the long strong loops which are not transcirption dependent 2)  
##      b-no CTCF at any anchor
##      c-all the rest
##   7- VR asked for:
##      a-Rad21 colocalization in loop contacts for the same categorization suggested by Akis
##      b-number of Rad21 %over% CTCF peaks, and viceversa
##      c-supercoiling signal at TSS - body - TES (as separate panels)
##
####################################
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(ggplot2)
library(cowplot)
library(parallel)
library(Rcpp)
library(grid)
library(pheatmap)
library(RColorBrewer)

CORES <- 16
PROJECT <- "/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2021_01_meta_HCT116_breakome_prediction"
setwd(PROJECT)
options(mc.cores=CORES)
theme_set(theme_cowplot(font_size=8))

pdf("results/supercoiling.pdf")

##
## read data
##
# Blacklisted regions
bl <- import.bed("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz")  # ENCODE's blacklisted regions

# GapR
x <- read.delim("data/imb_roukos_2022_03_longo_GapR_CnTSMC1/GapR/results/idr/GFP_WT.vs.GFP_MT.bed", head=FALSE, comment.char="#")
gapr.peaks <- with(x, GRanges(seqnames=x$V1, IRanges(start=x$V2, end=x$V3), signal=x$V7))
gapr.peaks <- gapr.peaks[!(gapr.peaks %over% bl)]

gapr.r1 <- import("data/imb_roukos_2022_03_longo_GapR_CnTSMC1/GapR/tracks/filtered/imb_roukos_2022_03_01_gapr_GFP_WT_1.unique.duprm.bw")
gapr.r2 <- import("data/imb_roukos_2022_03_longo_GapR_CnTSMC1/GapR/tracks/filtered/imb_roukos_2022_03_02_gapr_GFP_WT_2.unique.duprm.bw")
gapr.signal <- (coverage(gapr.r1, weight="score") + coverage(gapr.r2, weight="score")) / 2
gapr.signal <- gapr.signal[standardChromosomes(Hsapiens)]

# Cohesin
rad21.peaks <- import.bb("extdata/cohesin_dependant_insulation/ENCFF568PEO.bigBed")
rad21.peaks <- rad21.peaks[!(rad21.peaks %over% bl)]

rad21.signal <- coverage(import("extdata/cohesin_dependant_insulation/ENCFF027QAE.bigWig"), weight="score")
rad21.signal <- rad21.signal[standardChromosomes(Hsapiens)]

# CTCF
x <- read.delim("extdata/ctcf_loop_anchoring/ENCFF273SBR.bed", head=FALSE, comment.char="#")
ctcf.peaks <- with(x, GRanges(seqnames=x$V1, IRanges(start=x$V2, end=x$V3), signal=x$V7))
ctcf.peaks <- ctcf.peaks[!(ctcf.peaks %over% bl)]

ctcf.signal <- coverage(import("extdata/ctcf_loop_anchoring/ENCFF447IHY.bigWig"), weight="score")
ctcf.signal <- ctcf.signal[standardChromosomes(Hsapiens)]

# TOP2A/B
top2a.signal  <- coverage(import("extdata/top2AB/top2a.bw"), weight="score")
top2b.signal  <- coverage(import("extdata/top2AB/top2b.bw"), weight="score")
endseq.signal <- import("extdata/top2AB/end_seq.bw")
endseq.signal$score <- log2(1 + endseq.signal$score)
endseq.signal <- coverage(endseq.signal, weight="score")

# H3K27AC
x <- read.delim("extdata/histone_marks/GSE101966_GSM2719748_H3K27ac_1.bed", head=FALSE, comment.char="#")
h3k27ac.peaks <- with(x, GRanges(seqnames=x$V1, IRanges(start=x$V2, end=x$V3), signal=x$V5))

# LADs
# merge regions separated by no more than 250k)
# C version is fast like hell compared to anything that could be done in R
cppFunction('List mergeRegionsC(NumericVector start, NumericVector end, NumericVector score, int maxgap) {
  int n = start.size();
  LogicalVector merged(n);
  for(int i = 0; i < n-1; i++) {
    if((start[i+1] - end[i] <= maxgap) &&
       ((score[i] >= 0 && score[i+1] >= 0) ||
        (score[i] <  0 && score[i+1] <  0))) {
      //merge i with the next tile
      start[i+1] = start[i];
      score[i+1] = score[i] > 0 ? max(score[Rcpp::Range(i, i+1)]) : min(score[Rcpp::Range(i, i+1)]);
      merged[i] = true;
    }
  }
  return List::create(_["start"] = start, _["score"] = score, _["merged"]=merged);
}')

importLADs <- function(f) {
  lad <- import.bed(f)
  lad <- lad[!(lad %over% bl)]
  lad <- lapply(split(lad, seqnames(lad)), function(lad) { # to avoid comparing last of chr_i with first of chr_i+1
    x <- mergeRegionsC(start(lad), end(lad), lad$score, 250000)
    start(lad) <- x$start   # fetch new region start
    lad$score  <- x$score   # fetch new region score
    lad[!x$merged]
  })
  unlist(GRangesList(lad))
}

LADs <- list(dko=importLADs("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.IGG.LADs.bed"),
             wt =importLADs("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/TOP2AAID_NoAUX_LaminB1.vs.IGG.LADs.bed"))

diffLADs.gaps  <- unstrand(import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.TOP2AAID_NoAUX_LaminB1.LADs_diff.only_1e5_diffs.bed"))

# ttseq
ann <- getBM(attributes=c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand"),
             filters   ="biotype", values    ="protein_coding",
             mart      =useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=110))

read.ttseq <- function(f) {
  f <- list.files("data/sfb_roukos_2022_09_longo_ttseq_top2/results/subread-count", pattern=f, full.names=TRUE)
  x <- do.call(cbind, lapply(f, read.delim, head=FALSE, row.names=1))
  x <- apply(x, 2, function(x) 1e6 * x / sum(x))   # normalize to RPM
  x <- apply(x, 1, mean)                           # average replicates
  names(x) <- sub("\\.\\d+$", "", names(x))        # remove version number from gene name
  x <- x[names(x) %in% ann$ensembl_gene_id]        # subset protein coding genes only
  i <- match(names(x), ann$ensembl_gene_id)
  GRanges(ann$chromosome_name[i], IRanges(ann$start_position[i], ann$end_position[i]),
          strand=ann$strand[i], id=ann$ensembl_gene_id[i], gene=ann$external_gene_name[i], rpm=x)
}

ttseq <- list(dko=read.ttseq("sfb_roukos_2022_09_[12]_WT_[12].readcounts.tsv"),
              wt =read.ttseq("sfb_roukos_2022_09_[78]_DKO_[12].readcounts.tsv"))
ttseq <- lapply(ttseq, function(x) { x$rpkm <- 1e3 * x$rpm / width(x); x })   # normalize expression by gene length
seqlevelsStyle(ttseq$dko) <- "UCSC"
seqlevelsStyle(ttseq$wt)  <- "UCSC"

read.dge <- function(f) {
  x <- read.csv(f)
  x$BH.adjusted.p.values[is.na(x$BH.adjusted.p.values)] <- 1  # correct p-value for non-tested genes
  x$gene_id <- sub("\\.\\d+$", "", x$gene_id)
  makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)
}
dge <- read.dge("data/sfb_roukos_2022_09_longo_ttseq_top2/results/DE_DESeq2/DKO.vs.WT.csv")
dge <- dge[dge$gene_id %in% ttseq$wt$id]
seqlevelsStyle(dge)  <- "UCSC"

# MicroC loops
microc.loops <- list(dko=import("data/micro_c_top2_mutants/analysis_papantonis/DKO_g1.loops",  format="bedpe"),
                     wt =import("data/micro_c_top2_mutants/analysis_papantonis/ctrl_g1.loops", format="bedpe"))

seqlevelsStyle(first(microc.loops$wt))   <- "UCSC"
seqlevelsStyle(second(microc.loops$wt))  <- "UCSC"
seqlevelsStyle(first(microc.loops$dko))  <- "UCSC"
seqlevelsStyle(second(microc.loops$dko)) <- "UCSC"

##
##   1- Does supercoiling accumulate at LAD borders?
##
# function to meta-plot signal around a Pairs GRanges object
plotPairsSignal <- function(regions, signal, tiles=100, what="", main="", sub="", extend=c("", "")) {

  seqlevels(first(regions)) <-
  seqlevels(second(regions)) <-
  seqlevels(mcols(regions)$within) <-
  seqlevels(signal, pruning.mode="coarse") <-
    standardChromosomes(Hsapiens)

  # bin regions + calculate aggregated average signal per bin
  get_vector <- function(regions, signal) {
    # define regions
    regions  <- keepSeqlevels(regions , names(signal), pruning.mode="coarse")
    regions.binned <- unlist(tile(regions, n=tiles))
    regions.binned$region <- rep(seq_along(regions), each=tiles)

    # get binned signal average
    x <- binnedAverage(regions.binned, signal, "score")
    x <- matrix(x$score, ncol=tiles, byrow=TRUE)
    colSums(x, na.rm=TRUE) / length(regions.binned)   # aggregate signal from multiple regions, normalized to the number of regions
  }

  # get_vector for the 3 subregions: start of the loop, body of the loop, end of the loop
  x <- mclapply(list(left=first(regions), within=mcols(regions)$within, right=second(regions)), get_vector, signal)
  x <- c(x$left, x$within, x$right)

  # plot something quick and dirty
  plot(x, type="n", xlab=NA, ylab="arbitrary units", axes=FALSE, main=main, sub=sub, cex=2/3, cex.lab=2/3, cex.main=2/3, cex.sub=2/3)
  Axis(side=1, at=c(1, tiles/2, 3*tiles/2, 3*tiles-tiles/2, 3*tiles), cex.axis=2/3,
       labels=c(extend[1], "left\nboundary", paste(what, "center"), "right\nboundary", extend[2]))
  Axis(side=2, cex.axis=2/3)
  lines(x)
  abline(v=c(tiles/2, 3*tiles-tiles/2), lty=2, col='red')
}

# divide LADs into small (<1e5), medium (<1e6) and large (>1e6)
# define lad borders and call the plot
opar <- par(mfrow=c(2, 2))
EXTEND <- 5000
x <- LADs$wt[width(LADs$wt) < 1e5 & width(LADs$wt) > 2 * EXTEND]
x <- Pairs(first =GRanges(seqnames(x), IRanges(start(x) - EXTEND, start(x) + EXTEND)),
           within=GRanges(seqnames(x), IRanges(start(x) + EXTEND, end(x)   - EXTEND)),
           second=GRanges(seqnames(x), IRanges(end(x)   - EXTEND, end(x)   + EXTEND)))
plotPairsSignal(x, gapr.signal, what="LAD", extend=c("-5KB", "+5KB"),
                main="Supercoiling accumulation at LAD boundaries", sub="Small LADs <100KB")

EXTEND <- 50000
x <- LADs$wt[width(LADs$wt) > 1e5 & width(LADs$wt) < 1e6]
x <- Pairs(first =GRanges(seqnames(x), IRanges(start(x) - EXTEND, start(x) + EXTEND)),
           within=GRanges(seqnames(x), IRanges(start(x) + EXTEND, end(x)   - EXTEND)),
           second=GRanges(seqnames(x), IRanges(end(x)   - EXTEND, end(x)   + EXTEND)))
plotPairsSignal(x, gapr.signal, what="LAD", extend=c("-50KB", "+50KB"),
                main="Supercoiling accumulation at LAD boundaries", sub="Medium LADs <1MB")

EXTEND <- 500000
x <- LADs$wt[width(LADs$wt) > 1e6]
x <- Pairs(first =GRanges(seqnames(x), IRanges(start(x) - EXTEND, start(x) + EXTEND)),
           within=GRanges(seqnames(x), IRanges(start(x) + EXTEND, end(x)   - EXTEND)),
           second=GRanges(seqnames(x), IRanges(end(x)   - EXTEND, end(x)   + EXTEND)))
plotPairsSignal(x, gapr.signal, what="LAD", extend=c("-500KB", "+500KB"),
                main="Supercoiling accumulation at LAD boundaries", sub="Large LADs >1MB")
par(opar)

##
##   2- Separate MicroC loops into small/medium/large and answer if:
##

##
##      2.1- Supercoiling accumulates in the contacts?
##
# divide loops into small (<1e5), medium (<1e6) and large (>1e6)
# define lad borders and call the plot
w <- start(second(microc.loops$wt)) - end(first(microc.loops$wt))   # loop width

opar <- par(mfrow=c(2, 2))
EXTEND <- 10000
x <- microc.loops$wt[w < 1e5 & w > 2 * EXTEND + 100] # +100, to guarantee we have at least 100 tiles in the body (default in plotPairsSignal)
x <- Pairs(first =GRanges(seqnames(first(x)), IRanges(start(first(x))  - EXTEND, end(first(x))    + EXTEND)),
           within=GRanges(seqnames(first(x)), IRanges(end(first(x))    + EXTEND, start(second(x)) - EXTEND)),
           second=GRanges(seqnames(first(x)), IRanges(start(second(x)) - EXTEND, end(second(x))   + EXTEND)))
plotPairsSignal(x, gapr.signal, what="loop", extend=c("-10KB", "+10KB"), 
                main="Supercoiling accumulation at loop boundaries", sub="Small loops <100KB")

EXTEND <- 10000
x <- microc.loops$wt[w > 1e5 & w < 1e6 & w > 2 * EXTEND + 100]
x <- Pairs(first =GRanges(seqnames(first(x)), IRanges(start(first(x))  - EXTEND, end(first(x))    + EXTEND)),
           within=GRanges(seqnames(first(x)), IRanges(end(first(x))    + EXTEND, start(second(x)) - EXTEND)),
           second=GRanges(seqnames(first(x)), IRanges(start(second(x)) - EXTEND, end(second(x))   + EXTEND)))
plotPairsSignal(x, gapr.signal, what="loops", extend=c("-10KB", "+10KB"),
                main="Supercoiling accumulation at loop boundaries", sub="Medium loops <1MB")

EXTEND <- 10000
x <- microc.loops$wt[w >= 1e6 & w > 2 * EXTEND + 100]
x <- Pairs(first =GRanges(seqnames(first(x)), IRanges(start(first(x))  - EXTEND, end(first(x))    + EXTEND)),
           within=GRanges(seqnames(first(x)), IRanges(end(first(x))    + EXTEND, start(second(x)) - EXTEND)),
           second=GRanges(seqnames(first(x)), IRanges(start(second(x)) - EXTEND, end(second(x))   + EXTEND)))
plotPairsSignal(x, gapr.signal, what="loops", extend=c("-10KB", "+10KB"),
                main="Supercoiling accumulation at loop boundaries", sub="Large loops >=1MB")
par(opar)

##
##      2.2- the number of supercoiling peaks within the loop differ?
##
EXTEND <- 50000
x <- unlist(GRangesList(lapply(list(`small <100KB`=microc.loops$wt[w < 1e5],
                                    `medium <1MB` =microc.loops$wt[w > 1e5 & w < 1e6],
                                    `large >=1MB` =microc.loops$wt[w >= 1e6]),
                               function(x) GRanges(seqnames(first(x)), IRanges(start(first(x)) - EXTEND, end(second(x)) + EXTEND)))))
x$class <- names(x)   # names(x) is small, medium, large
x$gapr.peaks <- countOverlaps(x, gapr.peaks)

# number of loops with at least one GapR peak
df <- reshape2::melt(lapply(split(x$gapr.peaks, x$class), function(x) table(x > 0)))
p1 <- ggplot(df, aes(x=L1, y=value, fill=Var1)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual("At least 1 GapR peak", values=palette()) +
  labs(title="Loops with at least one GapR peak", subtitle="Including contacts and loop body",
       x="loop size", y="fraction of loops") +
  theme(legend.position="bottom")

# unnormalized
df <- as.data.frame(mcols(x))
p2 <- ggplot(df, aes(x=class, y=gapr.peaks)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title="GapR peaks within loops", subtitle="Including contacts and loop body",
       x="loop size", y="number of GapR peaks (log10)")

# normalized by loop size (GapR peaks per KB of loop)
df$gapr.peaks.norm <- 1e3 * df$gapr.peaks / width(x)
p3 <- ggplot(df, aes(x=class, y=gapr.peaks.norm)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title="GapR peaks within loops", subtitle="Including contacts and loop body\nnormalized by KB of loop length",
       x="loop size", y="number of GapR peaks\nper KB of loop (log10)")

cowplot::plot_grid(p1, p2, p3, nrow=2, ncol=2)

##
## Repeat the same, but excluding the loop body
##
x.loop.body <- unlist(GRangesList(lapply(list(`small <100KB`=microc.loops$wt[w < 1e5],
                                              `medium <1MB` =microc.loops$wt[w > 1e5 & w < 1e6],
                                              `large >=1MB` =microc.loops$wt[w >= 1e6]),
                                         function(x) GRanges(seqnames(first(x)), IRanges(end(first(x)), start(second(x)))))))
x$gapr.peaks.body    <- countOverlaps(x.loop.body, gapr.peaks)
x$gapr.peaks.contacts <- x$gapr.peaks - x$gapr.peaks.body

# number of loops with at least one GapR peak
df <- reshape2::melt(lapply(split(x$gapr.peaks.contacts, x$class), function(x) table(x > 0)))
p1 <- ggplot(df, aes(x=L1, y=value, fill=Var1)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual("At least 1 GapR peak", values=palette()) +
  labs(title="Loops with at least one GapR peak", subtitle="Within loop contacts",
       x="loop size", y="fraction of loops") +
  theme(legend.position="bottom")

# unnormalized
df <- as.data.frame(mcols(x))
p2 <- ggplot(df, aes(x=class, y=gapr.peaks.contacts)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title="GapR peaks within loops", subtitle="Within loop contacts",
       x="loop size", y="number of GapR peaks (log10)")

cowplot::plot_grid(p1, p2, nrow=2, ncol=2)

##
## Repeat the same, but for loop bodies excluding peaks on contacts of nested loops
##
x.loop.contacts <- unlist(GRangesList(lapply(list(`small <100KB`=microc.loops$wt[w < 1e5],
                                                 `medium <1MB` =microc.loops$wt[w > 1e5 & w < 1e6],
                                                 `large >=1MB`  =microc.loops$wt[w >= 1e6]),
                                         function(x) c(GRanges(seqnames(first(x)), IRanges(start(first(x))  - EXTEND, end(first(x)))),
                                                       GRanges(seqnames(first(x)), IRanges(start(second(x)), end(second(x)) + EXTEND))))))
gapr.peaks.body        <- gapr.peaks[countOverlaps(gapr.peaks, x.loop.contacts) == 0]  # keep GapR peaks not on contacts
x$gapr.peaks.body.only <- countOverlaps(x.loop.body, gapr.peaks.body)

# number of loops with at least one GapR peak
df <- reshape2::melt(lapply(split(x$gapr.peaks.body.only, x$class), function(x) table(x > 0)))
p1 <- ggplot(df, aes(x=L1, y=value, fill=Var1)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual("At least 1 GapR peak", values=palette()) +
  labs(title="Loops with at least one GapR peak", subtitle="Within bodies exluding contacts",
       x="loop size", y="fraction of loops") +
  theme(legend.position="bottom")

# unnormalized
df <- as.data.frame(mcols(x))
p2 <- ggplot(df, aes(x=class, y=gapr.peaks.body.only)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title="GapR peaks within loops", subtitle="Within bodies exluding contacts",
       x="loop size", y="number of GapR peaks (log10)")

# normalized by loop size (GapR peaks per KB of loop)
df$gapr.peaks.norm <- 1e3 * df$gapr.peaks.body.only / width(x)
p3 <- ggplot(df, aes(x=class, y=gapr.peaks.norm)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(title="GapR peaks within loops", subtitle="Within bodies exluding contacts\nnormalized by KB of loop length",
       x="loop size", y="number of GapR peaks\nper KB of loop (log10)")

cowplot::plot_grid(p1, p2, p3, nrow=2, ncol=2)

##
##      2.3- the length of the supercoiling peaks within the loop differ, depending on the size of the loop?
##
get_signal_vector <- function(regions, signal, tiles=200) {
  seqlevels(regions$first, pruning.mode="coarse") <-
  seqlevels(regions$second, pruning.mode="coarse") <-
  seqlevels(regions$within, pruning.mode="coarse") <-
  seqlevels(signal, pruning.mode="coarse") <-
    standardChromosomes(Hsapiens)

  get_vector <- function(regions, signal) {
    # define regions
    regions  <- keepSeqlevels(regions , names(signal), pruning.mode="coarse")
    regions.binned <- unlist(tile(regions, n=tiles))
    regions.binned$region <- rep(seq_along(regions), each=tiles)

    # get binned signal average
    x <- binnedAverage(regions.binned, signal, "score")
    x <- matrix(x$score, ncol=tiles, byrow=TRUE)
    colSums(x, na.rm=TRUE) / length(regions.binned)   # aggregate signal from multiple regions, normalized to the number of regions
  }

  # get_vector for the 3 subregions: start of the loop, body of the loop, end of the loop
  x <- mclapply(regions, get_vector, signal)
  c(x$first, x$within, x$second)
}

w <- end(second(microc.loops$wt)) - start(first(microc.loops$wt))
loop.widths <- factor(ifelse(w < 1e5, "<100KB",
                      ifelse(w < 1e6, "<1MB", ">=1MB")))

EXTEND <- 1e4
TILES  <- 200
df <- do.call(rbind, mclapply(split(microc.loops$wt, loop.widths), function(x) {
  regions <- list(first =GRanges(seqnames(first(x)), IRanges(start(first(x))  - EXTEND, end(first(x))    + EXTEND)),
                  second=GRanges(seqnames(first(x)), IRanges(start(second(x)) - EXTEND, end(second(x))   + EXTEND)),
                  within=GRanges(seqnames(first(x)), IRanges(end(first(x))    + EXTEND, start(second(x)) - EXTEND)))
  get_signal_vector(regions, gapr.signal, tiles=TILES)
}))

pheatmap::pheatmap(df, cluster_rows=FALSE, cluster_cols=FALSE, legend=FALSE, main="Aggregated GapR signal depending on loop size")
downViewport("matrix.4-3-4-3")
grid.lines(x=1 * (TILES / 2) / 600, y=c(0,1), gp=gpar(col="black", lwd=2, lty=2))  # annotate middle of the first contact region
grid.lines(x=5 * (TILES / 2) / 600, y=c(0,1), gp=gpar(col="black", lwd=2, lty=2))  # annotate middle of the second contact region
popViewport()

##
##   3- Classifiy supercoiling at the loop contacts taking into account the
##      directionality of the transcription (TTchem-seq)
##
# How I'll do it? 3 scenarios (ok, maybe 4 if we consider no transcription):
#
#     anchor           anchor              anchor
# 5'--> | <--3'    5'--> |                   | <--3'
# ttseq = ttseq    ttseq > 4*ttseq   4*ttseq < ttseq
#
# Do it separately for left and right contacts
#
EXTEND <- 1e6   # region around the anchor boundaries to look at genes
contacts <- lapply(list(left=first(microc.loops$wt), right=second(microc.loops$wt)), function(x) {
  Pairs(first =GRanges(seqnames(x), IRanges(mid(x) - EXTEND, mid(x))),
        second=GRanges(seqnames(x), IRanges(mid(x)         , mid(x) + EXTEND)))
})

# look first at the right anchor (no idea how to do it for the right (second)  anchor yet)
# pick only genes transcribing 5'-->3' toward the loop
i <- as.data.frame(findOverlaps(first(contacts$left), ttseq$wt[decode(strand(ttseq$wt)) == "+"]))
i$rpkm <- ttseq$wt$rpkm[i$subjectHits]
rpkm <- reshape2::melt(lapply(split(i$rpkm, i$queryHits), sum))
mcols(contacts$left)$rpkm_first_5 <- 0
mcols(contacts$left)$rpkm_first_5[as.numeric(rpkm$L1)] <- rpkm$value

# pick only genes transcribing 3'-->5' toward the loop
i <- as.data.frame(findOverlaps(second(contacts$left), ttseq$wt[decode(strand(ttseq$wt)) == "-"]))
i$rpkm <- ttseq$wt$rpkm[i$subjectHits]
rpkm <- reshape2::melt(lapply(split(i$rpkm, i$queryHits), sum))
mcols(contacts$left)$rpkm_first_3 <- 0
mcols(contacts$left)$rpkm_first_3[as.numeric(rpkm$L1)] <- rpkm$value

# classify directionality of the transcription (see classes in above comment)
mcols(contacts$left)$class <- ifelse(mcols(contacts$left)$rpkm_first_5 > 4 * mcols(contacts$left)$rpkm_first_3, "left only",
                              ifelse(mcols(contacts$left)$rpkm_first_3 > 4 * mcols(contacts$left)$rpkm_first_5, "right only", "both sides"))
mcols(contacts$left)$class <- ifelse(with(mcols(contacts$left), (class == "both sides" & rpkm_first_5 == 0 & rpkm_first_3 == 0)),
                                    "none", mcols(contacts$left)$class)

# now, plot the supercoiling accumulation around the (left) anchor
# function to meta-plot signal around a Pairs GRanges object
plotSignal <- function(regions, signal, tiles=100, arrow_direction="none", main="", sub="", what="", extend=c("", "")) {

  seqlevels(regions) <-
  seqlevels(signal, pruning.mode="coarse") <-
    standardChromosomes(Hsapiens)

  # bin regions + calculate aggregated average signal per bin
  get_vector <- function(regions, signal) {
    # define regions
    regions  <- keepSeqlevels(regions , names(signal), pruning.mode="coarse")
    regions.binned <- unlist(tile(regions, n=tiles))
    regions.binned$region <- rep(seq_along(regions), each=tiles)

    # get binned signal average
    x <- binnedAverage(regions.binned, signal, "score")
    x <- matrix(x$score, ncol=tiles, byrow=TRUE)
    colSums(x, na.rm=TRUE) / length(regions.binned)   # aggregate signal from multiple regions, normalized to the number of regions
  }

  # get_vector for the 3 subregions: start of the loop, body of the loop, end of the loop
  x <- get_vector(regions, signal)

  # plot something quick and dirty
  plot(x, type="n", xlab=NA, ylab="arbitrary units", axes=FALSE, main=main, sub=sub, cex=2/3, cex.lab=2/3, cex.main=2/3, cex.sub=2/3)
  Axis(side=1, at=c(1, tiles/2, tiles), labels=c(extend[1], what, extend[2]), cex.axis=2/3)
  Axis(side=2, cex.axis=2/3)
  lines(x)
  abline(v=tiles/2, lty=2)

  # add left and right average
  mean.left  <- mean(x[1:(tiles/2)])
  mean.right <- mean(x[(tiles/2+1):tiles])
  segments(x0=1      , x1=tiles/2, y0=mean.left , y1=mean.left , lty=2, col="blue")
  segments(x1=tiles/2, x0=tiles  , y0=mean.right, y1=mean.right, lty=2, col="red")
  text(x=1    , y=mean.left , col="blue", labels=paste0("mean=", signif(mean.left , digits=2)), adj=c(0, -0.2), cex=2/3)
  text(x=tiles, y=mean.right, col="red" , labels=paste0("mean=", signif(mean.right, digits=2)), adj=c(1, -0.2), cex=2/3)

  # add transcription arrows
  switch(arrow_direction,
         "both sides"={ arrows(x0=2*tiles/8, x1=3*tiles/8, y0=max(x), y1=max(x), length=.08);
                        arrows(x1=5*tiles/8, x0=6*tiles/8, y0=max(x), y1=max(x), length=.08);
                        text(x=2*tiles/8, y=max(x), adj=c(0, -0.2), labels="transcription", cex=2/3);
                        text(x=6*tiles/8, y=max(x), adj=c(1, -0.2), labels="transcription", cex=2/3) },
         "left only" ={ arrows(x0=2*tiles/8, x1=3*tiles/8, y0=max(x), y1=max(x), length=.08);
                        text(x=2*tiles/8, y=max(x), adj=c(0, -0.2), labels="transcription", cex=2/3) },
         "right only"={ arrows(x1=5*tiles/8, x0=6*tiles/8, y0=max(x), y1=max(x), length=.08);
                        text(x=6*tiles/8, y=max(x), adj=c(1, -0.2), labels="transcription", cex=2/3) })
}

transcribed <- with(mcols(contacts$left), (class == "both sides" & rpkm_first_5 > 1 & rpkm_first_3 > 1) |
                                         (class == "left only"  & rpkm_first_5 > 1) |
                                         (class == "right only" & rpkm_first_3 > 1))
opar <- par(mfrow=c(2, 2))
lapply(split(contacts$left[transcribed], mcols(contacts$left[transcribed])$class), function(x) {
  regions <- GRanges(seqnames(first(x)), IRanges(start(first(x)), end(second(x))))
  plotSignal(regions, gapr.signal, what="left anchor", extend=c("-1MB", "+1MB"), tiles=1000,
             arrow_direction=mcols(x)$class[1],
             main="GapR signal around the left anchor",
             sub=paste("Direction of transcription", mcols(x)$class[1], "(", length(x), "sites)"))
#  readline("...")
})
par(opar)

##
##   4- How many loop contacts overlap a CTCF+Rad21 peak? <-- determine anchored loop contacts
##
df <- list(rad21=list(first =sum( first(microc.loops$wt) %over% rad21.peaks & !second(microc.loops$wt) %over% rad21.peaks),
                      second=sum(!first(microc.loops$wt) %over% rad21.peaks &  second(microc.loops$wt) %over% rad21.peaks),
                      both  =sum( first(microc.loops$wt) %over% rad21.peaks &  second(microc.loops$wt) %over% rad21.peaks),
                      none  =sum(!first(microc.loops$wt) %over% rad21.peaks & !second(microc.loops$wt) %over% rad21.peaks)),
           ctcf =list(first =sum( first(microc.loops$wt) %over% ctcf.peaks  & !second(microc.loops$wt) %over% ctcf.peaks),
                      second=sum(!first(microc.loops$wt) %over% ctcf.peaks  &  second(microc.loops$wt) %over% ctcf.peaks),
                      both  =sum( first(microc.loops$wt) %over% ctcf.peaks  &  second(microc.loops$wt) %over% ctcf.peaks),
                      none  =sum(!first(microc.loops$wt) %over% ctcf.peaks  & !second(microc.loops$wt) %over% ctcf.peaks)),
           `ctcf+rad21`=list(first =sum(  first(microc.loops$wt) %over% rad21.peaks &   first(microc.loops$wt) %over% ctcf.peaks &
                                        !second(microc.loops$wt) %over% rad21.peaks & !second(microc.loops$wt) %over% ctcf.peaks),
                             second=sum(! first(microc.loops$wt) %over% rad21.peaks & ! first(microc.loops$wt) %over% ctcf.peaks &
                                         second(microc.loops$wt) %over% rad21.peaks &  second(microc.loops$wt) %over% ctcf.peaks),
                             both  =sum(  first(microc.loops$wt) %over% rad21.peaks &   first(microc.loops$wt) %over% ctcf.peaks &
                                         second(microc.loops$wt) %over% rad21.peaks &  second(microc.loops$wt) %over% ctcf.peaks),
                             none  =sum(! first(microc.loops$wt) %over% rad21.peaks & ! first(microc.loops$wt) %over% ctcf.peaks &
                                        !second(microc.loops$wt) %over% rad21.peaks & !second(microc.loops$wt) %over% ctcf.peaks))) |>
      reshape2::melt()
df$L1 <- factor(df$L1, levels=c("ctcf", "rad21", "ctcf+rad21"))
df$L2 <- factor(df$L2, levels=c("none", "first", "second", "both"),
                labels=c("none", "first contact end", "second contact end", "both contact ends"))

p1 <-  ggplot(df, aes(x=L1, y=value, fill=L2, label=value)) +
    geom_bar(stat="identity", position="stack") +
    geom_text(position="stack", vjust=1) +
    scale_fill_viridis_d("Loop contact end (10KB)\nhas CTCF+Rad21 peaks", direction=-1) +
    labs(subtitle="How many loop contacts overlap a CTCF + Rad21 peaks", title="Determine anchored loop contacts", x="", y="") +
    theme(legend.position="bottom")

##
##   5- Now, look at the loops not anchored (!Rad21 !CTCF) if they have supercoiling or not
##
anchored <- ifelse(  first(microc.loops$wt) %over% rad21.peaks &   first(microc.loops$wt) %over% ctcf.peaks &
                   !second(microc.loops$wt) %over% rad21.peaks & !second(microc.loops$wt) %over% ctcf.peaks , "first",
            ifelse(! first(microc.loops$wt) %over% rad21.peaks & ! first(microc.loops$wt) %over% ctcf.peaks &
                    second(microc.loops$wt) %over% rad21.peaks &  second(microc.loops$wt) %over% ctcf.peaks , "second",
            ifelse(  first(microc.loops$wt) %over% rad21.peaks &   first(microc.loops$wt) %over% ctcf.peaks &
                    second(microc.loops$wt) %over% rad21.peaks &  second(microc.loops$wt) %over% ctcf.peaks , "both",
            ifelse(! first(microc.loops$wt) %over% rad21.peaks & ! first(microc.loops$wt) %over% ctcf.peaks &
                   !second(microc.loops$wt) %over% rad21.peaks & !second(microc.loops$wt) %over% ctcf.peaks , "none", "single binding"))))
anchored <- factor(anchored, levels=c("both", "first", "second", "none", "single binding"))

# size of loops
widths <- start(second(microc.loops$wt)) - end(first(microc.loops$wt))

df <- data.frame(w=widths, a=anchored)
p2 <- ggplot(df, aes(x=a, y=w)) +
  geom_boxplot() +
  labs(title="size of loops depending on Rad21+CTCF occupancy on the contacts", x="CTCF+Rad21 contact occupancy", y="loop width") +
  scale_y_log10()

# size of GapR peaks
df <- reshape2::melt(lapply(split(microc.loops$wt, anchored), function(x) {
  sum(width(gapr.peaks[gapr.peaks %over% x])) / length(x)
}))
df$L1 <- factor(df$L1, levels=c("both", "first", "second", "none", "single binding"))

p3 <- ggplot(df, aes(y=value, x=L1)) +
  geom_bar(stat="identity") +
  labs(title="size of GapR peaks depending on Rad21+CTCF occupancy on the contacts", x="CTCF+Rad21 contact occupancy",
       y="Average GapR peak size per loop")

cowplot::plot_grid(p1, p2, p3, ncol=2, nrow=2)

# shape of the supercoiling (metaloop plot)
plot_signal_vector <- function(signal, classes, EXTEND=1e4, TILES=200, ...) {
  # function to aggregate in a vector the signal per tile in a vector of tiles
  get_signal_vector <- function(regions, signal, tiles=200) {
    seqlevels(regions$first) <-
    seqlevels(regions$second) <-
    seqlevels(regions$within) <-
    seqlevels(signal, pruning.mode="coarse") <-
      standardChromosomes(Hsapiens)

    get_vector <- function(regions, signal) {
      # define regions
      regions  <- keepSeqlevels(regions , names(signal), pruning.mode="coarse")
      regions.binned <- unlist(tile(regions, n=tiles))
      regions.binned$region <- rep(seq_along(regions), each=tiles)

      # get binned signal average
      x <- binnedAverage(regions.binned, signal, "score")
      x <- matrix(x$score, ncol=tiles, byrow=TRUE)
      colSums(x, na.rm=TRUE) / length(regions.binned)   # aggregate signal from multiple regions, normalized to the number of regions
    }

    # get_vector for the 3 subregions: start of the loop, body of the loop, end of the loop
    x <- mclapply(regions, get_vector, signal)
    c(x$first, x$within, x$second)
  }

  # split the regions and get the signal vector per group
  df <- do.call(rbind, mclapply(split(microc.loops$wt, classes), function(x) {
    regions <- list(first =GRanges(seqnames(first(x)), IRanges(start(first(x))  - EXTEND, end(first(x))    + EXTEND)),
                    second=GRanges(seqnames(first(x)), IRanges(start(second(x)) - EXTEND, end(second(x))   + EXTEND)),
                    within=GRanges(seqnames(first(x)), IRanges(end(first(x))    + EXTEND, start(second(x)) - EXTEND)))
    get_signal_vector(regions, signal, tiles=TILES)
  }))

  # plot the summary as a heatmap
  pheatmap(df, cluster_rows=FALSE, cluster_cols=FALSE, legend=FALSE, ...)
  downViewport("matrix.4-3-4-3")
  grid.lines(x=1 * (TILES / 2) / 600, y=c(0,1), gp=gpar(col="black", lwd=2, lty=2))  # annotate middle of the first contact region
  grid.lines(x=5 * (TILES / 2) / 600, y=c(0,1), gp=gpar(col="black", lwd=2, lty=2))  # annotate middle of the second contact region
  popViewport()
}

plot_signal_vector(gapr.signal, classes=anchored, main="Aggregated GapR signal")

# shape of top2A/B for the very same regions
plot_signal_vector(top2a.signal , classes=anchored, main="Aggregated Top2A signal")
plot_signal_vector(top2b.signal , classes=anchored, main="Aggregated Top2B signal")
plot_signal_vector(endseq.signal, classes=anchored, main="Aggregated END-seq signal")

##
##   6- Same as 5- but with a different categorization suggested by Akis:
##      a-CTCF binding at both anchors without TSS and H3K27ac
##        (these are the long strong loops which are not transcirption dependent
##      b-no CTCF at any anchor
##      c-all the rest
##
# get all annotated TSS (limit to those with evidence of the isoform (==TSL1))
ann <- getBM(attributes=c("ensembl_gene_id", "chromosome_name", "transcript_start", "transcript_end", "strand", "transcript_tsl"),
             mart      =useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=110))
ann$tss <- ifelse(ann$strand == 1, ann$transcript_start, ann$transcript_end)
ann$tes <- ifelse(ann$strand == 1, ann$transcript_end  , ann$transcript_start)
ann <- ann[grepl("^tsl1", ann$transcript_tsl), ]   # limit to transcripts with high evidence
tss <- GRanges(ann$chromosome_name, IRanges(start=ann$tss, end=ann$tss))
seqlevelsStyle(tss) <- "UCSC"

# classify loops
akis <- ifelse(  first(microc.loops$wt) %over% ctcf.peaks    &  second(microc.loops$wt) %over% ctcf.peaks   &
                !first(microc.loops$wt) %over% tss           & !second(microc.loops$wt) %over% tss          &
                !first(microc.loops$wt) %over% h3k27ac.peaks & !second(microc.loops$wt) %over% h3k27ac.peaks, "a",
         ifelse(!first(microc.loops$wt) %over% ctcf.peaks    & !second(microc.loops$wt) %over% ctcf.peaks, "b", "c"))
akis <- factor(akis, levels=c("a", "b", "c"), labels=c("long strong loop\nno transcription dep.", "no CTCF at any anchor", "other"))


p1 <- ggplot(data.frame(classes=akis), aes(x=akis)) +
    geom_bar() +
    geom_text(stat='count', aes(label=..count..), vjust=-0.5, size=3) +
    labs(title="Number of loops ~ Akis' classes", x="", y="") 

# size of loops
widths <- start(second(microc.loops$wt)) - end(first(microc.loops$wt))

df <- data.frame(w=widths, a=akis)
p2 <- ggplot(df, aes(x=a, y=w)) +
  geom_boxplot() +
  labs(title="Size of loops ~ Akis' classes", x="", y="loop width") +
  scale_y_log10()

# size of GapR peaks
df <- reshape2::melt(lapply(split(microc.loops$wt, akis), function(x) {
  sum(width(gapr.peaks[gapr.peaks %over% x])) / length(x)
}))

p3 <- ggplot(df, aes(y=value, x=L1)) +
  geom_bar(stat="identity") +
  labs(title="size of GapR peaks ~ Akis' classes", x="", y="Average GapR peak size per loop")

cowplot::plot_grid(p1, p2, p3, ncol=2, nrow=2)

# shape of the supercoiling and top2AB (metaloop plot)
plot_signal_vector(gapr.signal  , classes=akis, main="Aggregated GapR signal")
plot_signal_vector(top2a.signal , classes=akis, main="Aggregated Top2A signal")
plot_signal_vector(top2b.signal , classes=akis, main="Aggregated Top2B signal")
plot_signal_vector(endseq.signal, classes=akis, main="Aggregated END-seq signal")

##
##   7- VR asked for:
##

##
##      a-Rad21 colocalization in loop contacts for the same categorization suggested by Akis
##
plot_signal_vector(rad21.signal , classes=akis, main="Aggregated Rad21 signal")

##
##      b-number of Rad21 %over% CTCF peaks, and viceversa
##
all.peaks <- reduce(c(rad21.peaks, ctcf.peaks))
x <- reshape2::melt(table(ifelse(all.peaks %over% rad21.peaks & all.peaks %over% ctcf.peaks, "Rad21 + CTCF",
                          ifelse(all.peaks %over% rad21.peaks, "Rad21 only", "CTCF only"))))

ggplot(x, aes(x=0, y=value, fill=Var1)) +
  geom_bar(stat="identity", position="stack") +
  geom_text(aes(label=value), position=position_stack(vjust=0.5)) +
  scale_fill_manual("", values=colorRampPalette(palette()[1:2])(3)) +
  scale_x_continuous(breaks=NULL) +
  labs(title="Colocalization of Rad21 and CTCF peaks", y="number of peaks", x="")

##
##      c-supercoiling signal at TSS - body - TES (as separate panels)
##
EXT  <- 2000
# NOTE: instead of using naive annotation, take the gene classification that GP provided
#----------
#ann2 <- ann[ann$transcript_end - ann$transcript_start > 4 * EXT, ]  # look at transcripts >8kb
ann <- do.call(rbind, list(read.delim("data/sfb_roukos_2022_09_longo_ttseq_top2.additional_data/high_expressed_genes_WT.bed", head=FALSE) |>
                             dplyr::mutate(class="high"),
                           read.delim("data/sfb_roukos_2022_09_longo_ttseq_top2.additional_data/medium_expressed_genes_WT.bed", head=FALSE) |>
                             dplyr::mutate(class="medium"),
                           read.delim("data/sfb_roukos_2022_09_longo_ttseq_top2.additional_data/low_expressed_genes_WT.bed", head=FALSE) |>
                             dplyr::mutate(class="low"),
                           read.delim("data/sfb_roukos_2022_09_longo_ttseq_top2.additional_data/noexp_expressed_genes_WT.bed", head=FALSE) |>
                             dplyr::mutate(class="noexp")))
colnames(ann) <- c("chromosome_name", "transcript_start", "transcript_end", "gene_name", "counts", "strand", "class")
ann2 <- ann[ann$transcript_end - ann$transcript_start > 4 * EXT, ]  # look at transcripts >8kb
ann2$tss <- ifelse(ann2$strand == "+", ann2$transcript_start, ann2$transcript_end)
ann2$tes <- ifelse(ann2$strand == "+", ann2$transcript_end  , ann2$transcript_start)
#----------
tss  <- GRanges(ann2$chromosome_name, IRanges(start=ann2$tss - EXT, end=ann2$tss + EXT), strand=ann2$strand, class=ann2$class)
tes  <- GRanges(ann2$chromosome_name, IRanges(start=ann2$tes - EXT, end=ann2$tes + EXT), strand=ann2$strand, class=ann2$class)
body <- GRanges(ann2$chromosome_name, IRanges(start=ann2$transcript_start + EXT, end=ann2$transcript_end - EXT), strand=ann2$strand, class=ann2$class)
i    <- tss %over% tes  # remove transcripts which its start may overlap other's end.
tss  <- tss[!i]
tes  <- tes[!i]
body <- body[!i]

get_vector <- function(regions, signal, tiles=201) {
  # define regions
  seqlevels(regions, pruning.mode="coarse") <-
    seqlevels(signal, pruning.mode="coarse") <-
    standardChromosomes(Hsapiens)
  regions  <- keepSeqlevels(regions , names(signal), pruning.mode="coarse")
  regions.binned <- unlist(tile(regions, n=tiles))
  regions.binned$region <- rep(seq_along(regions), each=tiles)

  # get binned signal average
  x <- binnedAverage(regions.binned, signal, "score")
  x <- matrix(x$score, ncol=tiles, byrow=TRUE)
  i <- which(decode(strand(regions)) == "+")   # reverse transcripts in the minus strand
  x[i, ] <- t(apply(x[i, ], 1, rev))
  #colSums(x, na.rm=TRUE)                       # aggregate signal from multiple regions
  colMeans(x, na.rm=TRUE)                       # aggregate signal from multiple regions
}

signal <- mclapply(setNames(unique(ann$class), unique(ann$class)), function(cl) {
  mclapply(list(tss =tss[tss$class == cl],
                body=body[body$class == cl],
                tes =tes[tes$class == cl]),
           get_vector, gapr.signal)
})

par(mfrow=c(3, 3), mar=c(2, 2, 2, 2))
Map(list("tss", "body", "tes"),
    list(c("-2kb", "TSS", "+2kb"), c("TSS+2kb", "gene\ncenter", "TES-2kb"), c("-2kb", "TES", "+2kb")),
    list(TRUE, FALSE, FALSE),
  f=function(feature, labels, show.legend) {
    x <- lapply(signal, function(signal) signal[[feature]])
    plot(x[[1]], main=feature, xlab="", ylab="", type="n", xaxt="n", ylim=range(unlist(x)))
    Map(x, palette()[1:length(x)], f=function(x, col) lines(x, col=col))
    if(show.legend) {
      legend("topright", names(x), fill=palette()[1:length(x)])
    }
    axis(1, at=c(1, 101, 201), labels=labels)
    axis(2)
})

##
## Let's do the same with gene length. Also provided by GP
## Since windows in body will be of variable length, let's look at TSS and TES only
##
ann <- do.call(rbind, list(read.delim("data/sfb_roukos_2022_09_longo_ttseq_top2.additional_data/hg38.knownGene.very_large.bed", head=FALSE) |>
                             dplyr::mutate(class="very large"),
                           read.delim("data/sfb_roukos_2022_09_longo_ttseq_top2.additional_data/hg38.knownGene.large.bed", head=FALSE) |>
                             dplyr::mutate(class="large"),
                           read.delim("data/sfb_roukos_2022_09_longo_ttseq_top2.additional_data/hg38.knownGene.medium.bed", head=FALSE) |>
                             dplyr::mutate(class="medium"),
                           read.delim("data/sfb_roukos_2022_09_longo_ttseq_top2.additional_data/hg38.knownGene.small.bed", head=FALSE) |>
                             dplyr::mutate(class="small")))
colnames(ann) <- c("chromosome_name", "transcript_start", "transcript_end", "gene_name", "counts", "strand", "class")
ann2 <- ann
ann2$tss <- ifelse(ann2$strand == "+", ann2$transcript_start, ann2$transcript_end)
ann2$tes <- ifelse(ann2$strand == "+", ann2$transcript_end  , ann2$transcript_start)
#----------
tss  <- GRanges(ann2$chromosome_name, IRanges(start=ann2$tss - EXT, end=ann2$tss + EXT), strand=ann2$strand, class=ann2$class)
tes  <- GRanges(ann2$chromosome_name, IRanges(start=ann2$tes - EXT, end=ann2$tes + EXT), strand=ann2$strand, class=ann2$class)
#i    <- tss %over% tes  # remove transcripts which its start may overlap other's end.
#tss  <- tss[!i]
#tes  <- tes[!i]

signal <- mclapply(setNames(unique(ann$class), unique(ann$class)), function(cl) {
  mclapply(list(tss =tss[tss$class == cl],
                tes =tes[tes$class == cl]),
           get_vector, gapr.signal)
})

Map(list("tss", "tes"),
    list(c("-2kb", "TSS", "+2kb"), c("-2kb", "TES", "+2kb")),
    list(TRUE, FALSE),
  f=function(feature, labels, show.legend) {
    x <- lapply(signal, function(signal) signal[[feature]])
    plot(x[[1]], main=feature, xlab="", ylab="", type="n", xaxt="n", ylim=range(unlist(x)))
    Map(x, palette()[1:length(x)], f=function(x, col) lines(x, col=col))
    if(show.legend) {
      legend("topright", names(x), fill=palette()[1:length(x)])
    }
    axis(1, at=c(1, 101, 201), labels=labels)
    axis(2)
})

dev.off()
