#################################
##
## A less naive analysis of Cohesin (Rad21) - Supercoiling (GapR)
## NOTE: for Rad21 we take the subset of Rad21 peaks overlapping CTCF (Rad21 seems to have other functions otherwise)
##
#################################
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(ggplot2)
library(eulerr)
library(parallel)

PROJECT="/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2021_01_meta_HCT116_breakome_prediction"
GTF="/fsimb/common/genomes/homo_sapiens/gencode/release-43_GRCh38.p13/canonical/annotation/gencode.v43.annotation.gtf"
setwd(PROJECT)
pdf("results/colocalize_supercoiling_cohesin.pdf")

##
## download peaks/signal tracks
##
# Blacklisted regions
bl <- import.bed("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz")  # ENCODE's blacklisted regions

# CTCF
x <- read.delim("extdata/ctcf_loop_anchoring/ENCFF273SBR.bed", head=FALSE, comment.char="#")
ctcf.peaks <- with(x, GRanges(seqnames=x$V1, IRanges(start=x$V2, end=x$V3), signal=x$V7))
ctcf.peaks <- ctcf.peaks[!(ctcf.peaks %over% bl)]

ctcf.signal <- coverage(import("extdata/ctcf_loop_anchoring/ENCFF447IHY.bigWig"), weight="score")
ctcf.signal <- ctcf.signal[standardChromosomes(Hsapiens)]

# RAD21
rad21.peaks <- import.bb("extdata/cohesin_dependant_insulation/ENCFF568PEO.bigBed")
rad21.peaks <- rad21.peaks[!(rad21.peaks %over% bl)]
rad21.peaks <- rad21.peaks[rad21.peaks %over% ctcf.peaks]

rad21.signal <- coverage(import("extdata/cohesin_dependant_insulation/ENCFF027QAE.bigWig"), weight="score")
rad21.signal <- rad21.signal[standardChromosomes(Hsapiens)]

# GapR
gapr.peaks <- makeGRangesFromDataFrame(read.delim("data/imb_roukos_2022_03_longo_GapR_CnTSMC1/GapR/results/idr/GFP_WT.vs.GFP_MT.bed", head=FALSE),
                                       seqnames.field="V1", start.field="V2", end.field="V3", keep=TRUE)
gapr.peaks <- gapr.peaks[!(gapr.peaks %over% bl)]

gapr.r1 <- import("data/imb_roukos_2022_03_longo_GapR_CnTSMC1/GapR/tracks/filtered/imb_roukos_2022_03_01_gapr_GFP_WT_1.unique.duprm.bw")
gapr.r2 <- import("data/imb_roukos_2022_03_longo_GapR_CnTSMC1/GapR/tracks/filtered/imb_roukos_2022_03_02_gapr_GFP_WT_2.unique.duprm.bw")
gapr.signal <- (coverage(gapr.r1, weight="score") + coverage(gapr.r2, weight="score")) / 2
gapr.signal <- gapr.signal[standardChromosomes(Hsapiens)]

# annotation
genes <- import(GTF)
genes <- genes[genes$type == "gene" & genes$gene_type == "protein_coding"]

seqlevels(rad21.peaks, pruning.mode="coarse") <-
  seqlevels(rad21.signal, pruning.mode="coarse") <-
  seqlevels(gapr.peaks, pruning.mode="coarse") <-
  seqlevels(gapr.signal, pruning.mode="coarse") <-
  standardChromosomes(Hsapiens)

##
## Colocalization of the 2 peaks
##
metaplot <- function(regions, signal, signal.name="", ext=3000, bins=201, add.random=FALSE, title="", subtitle="", doplot=FALSE) {
  # define regions
  regions <- keepSeqlevels(regions, names(signal), pruning.mode="coarse")
  middle <- mid(regions)
  start(regions) <- middle - ext
  end(regions)   <- middle + ext
  regions.binned <- unlist(tile(regions, n=bins))
  regions.binned$region <- rep(seq_along(regions), each=bins)

  # get binned signal average
  x <- binnedAverage(regions.binned, signal, "score")
  x <- matrix(x$score, ncol=bins, byrow=TRUE)
  x <- list(colSums(x, na.rm=TRUE) / length(regions.binned))   # aggregate signal from multiple regions, normalized to the number of regions
  names(x) <- signal.name

  # add some random regions to show the bg
  if(add.random) {
    regions.random <- shift(regions, shift=runif(length(regions), -10000, 10000)) # randomly shift region +/- 1M
    regions.random.binned <- unlist(tile(regions.random, n=bins))                   # tile into bins
    regions.random.binned$region <- rep(seq_along(regions.random), each=bins)
    x.random <- binnedAverage(regions.random.binned, signal, "score") # get binned signal average
    x.random <- matrix(x.random$score, ncol=bins, byrow=TRUE)
    x$random <- colSums(x.random, na.rm=T) / length(regions.random)   # aggregate signal + normalize to the number of regions
  }

  # plot something
  df    <- reshape2::melt(x)
  df$x  <- rep(1:bins, length(x))

  if(doplot) {
    pal   <- c(palette()[1], "grey30")
    ggplot(df, aes(x=x, y=value, color=L1)) +
      geom_line() +
      geom_vline(xintercept=(1+bins)/2, lty=2) +
      labs(title=title, subtitle=subtitle, y="normalized signal (arbitrary units)", x="") +
      scale_color_manual("", values=pal) +
      scale_x_continuous(breaks=c(0, (1+bins)/2, bins), labels=c(-ext, "center", ext)) +
      theme(legend.position="bottom")
  } else {
    df
  }
}

# GapR singal on Cohesin peaks
x <- list()
rad21.x <- rad21.peaks
ranges(rad21.x) <- start(rad21.x) + rad21.x$peak  # center the region around the summit
x$all <- metaplot(rad21.x, gapr.signal, signal.name="GapR", doplot=FALSE,
                  title="Supercoiling signal on Cohesin peaks",
                  subtitle=paste0("all ", length(rad21.x), " rad21 peaks)"))

tss.plus <- genes[decode(strand(genes)) == "+"]
end(tss.plus)   <- start(tss.plus) + 500
start(tss.plus) <- start(tss.plus) - 5000
x$fwd <- metaplot(rad21.x[rad21.x %over% tss.plus], gapr.signal, signal.name="GapR", doplot=FALSE,
                  title="Supercoiling signal on Cohesin peaks",
                  subtitle=paste0("on promoters of forward genes (", sum(rad21.x %over% tss.plus), " rad21 peaks)"))

tss.minus <- genes[decode(strand(genes)) == "-"]
start(tss.minus) <- end(tss.minus) - 500
end(tss.minus)   <- end(tss.minus) + 5000
x$rev <- metaplot(rad21.x[rad21.x %over% tss.minus], gapr.signal, signal.name="GapR", doplot=FALSE,
                  title="Supercoiling signal on Cohesin peaks",
                  subtitle=paste0("on promoters of reverse genes", sum(rad21.x %over% tss.minus), " rad21 peaks)"))

genes.tss <- genes
start(genes.tss) <- ifelse(decode(strand(genes)) == "+", start(genes) - 5000, start(genes))
end(genes.tss)   <- ifelse(decode(strand(genes)) == "+", end(genes), end(genes) + 5000)
x$opensea <- metaplot(rad21.x[!(rad21.x %over% genes.tss)], gapr.signal, signal.name="GapR", doplot=FALSE,
                     title="Supercoiling signal on Cohesin peaks",
                     subtitle=paste0("_not_ on genes", sum(!(rad21.x %over% genes.tss)), " rad21 peaks)"))

# now plot all together
x$all$L1     <- paste0("all ", length(rad21.x) ," peaks")
x$fwd$L1     <- paste0(sum(rad21.x %over% tss.plus), " peaks on gene promoters (forward)")
x$rev$L1     <- paste0(sum(rad21.x %over% tss.minus), " peaks on gene promoters (reverse)")
x$opensea$L1 <- paste0(sum(!(rad21.x %over% genes.tss)), " peaks not on genes (TSS+body)")
df <- do.call(rbind, x)

ext  <- 3000
bins <- 201
bin_width <- 2 * ext / bins
max_pos <- x$all$x[order(x$all$value, decreasing=TRUE)[1:2]]
ggplot(df, aes(x=x, y=value, color=L1)) +
  geom_line() +
  geom_vline(xintercept=(1+bins)/2, lty=2) +
  geom_vline(xintercept=max_pos, lty=2, color="red") +
  annotate("text", x=max(max_pos[2]) + 1, y=max(df$value), vjust=0, hjust=0,
           label=paste("peaks at", paste(round(max_pos * bin_width - bins * bin_width / 2 - bin_width / 2), collapse=" and "), "bp from center")) +
  labs(title="Supercoiling signal on Cohesin peaks", y="normalized signal (arbitrary units)", x="") +
  scale_color_discrete("Rad21+CTCF peaks") +
  scale_x_continuous(breaks=c(0, (1+bins)/2, bins), labels=c(-ext, "center", ext)) +
  theme(legend.position="bottom") +
  guides(color=guide_legend(nrow=2, byrow=TRUE))

##
## Venn diagram with the overlap
##
peaks <- reduce(c(rad21.peaks, gapr.peaks), drop.empty.ranges=TRUE)

peaks$`rad21+ctcf` <- peaks %over% rad21.peaks
peaks$gapr         <- peaks %over% gapr.peaks
peaks$`gene TSS`   <- peaks %over% tss.minus | peaks %over% tss.plus

plot(euler(as.data.frame(mcols(peaks)), shape="ellipse"), quantities=TRUE,
     main="Co-localization of supercoiling and\ncohesin insulation (Rad21 + CTCF)")

dev.off()
