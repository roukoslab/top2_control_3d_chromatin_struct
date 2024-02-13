##################################
##
## Some metrics to describe the colocalization between supercoiling and cohesin
## with microC loops
##
##################################
set.seed(666)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(eulerr)
library(ggplot2)
library(parallel)

PROJECT="/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2021_01_meta_HCT116_breakome_prediction"
GTF="/fsimb/common/genomes/homo_sapiens/gencode/release-43_GRCh38.p13/canonical/annotation/gencode.v43.annotation.gtf"
CORES=16
setwd(PROJECT)
options(mc.cores=CORES)
pdf("results/supercoiling_cohesin_vs_microC_loops.pdf")

##
## download peaks/signal tracks
##
# Blacklisted regions
bl <- import.bed("https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz")  # ENCODE's blacklisted regions

# MicroC loops
microc.loops <- list(ctrl=import("data/micro_c_top2_mutants/analysis_papantonis/ctrl_g1.loops", format="bedpe"),
                     dko =import("data/micro_c_top2_mutants/analysis_papantonis/DKO_g1.loops" , format="bedpe"))
seqlevelsStyle(first(microc.loops$ctrl))  <- "UCSC"
seqlevelsStyle(second(microc.loops$ctrl)) <- "UCSC"
seqlevelsStyle(first(microc.loops$dko))   <- "UCSC"
seqlevelsStyle(second(microc.loops$dko))  <- "UCSC"

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

seqlevels(first(microc.loops$ctrl), pruning.mode="coarse") <-
  seqlevels(second(microc.loops$ctrl), pruning.mode="coarse") <-
  seqlevels(first(microc.loops$dko), pruning.mode="coarse") <-
  seqlevels(second(microc.loops$dko), pruning.mode="coarse") <-
  seqlevels(rad21.peaks, pruning.mode="coarse") <-
  seqlevels(rad21.signal, pruning.mode="coarse") <-
  seqlevels(gapr.peaks, pruning.mode="coarse") <-
  seqlevels(gapr.signal, pruning.mode="coarse") <-
  standardChromosomes(Hsapiens)

##
## Metaplot with the RAD21/GapR signal in the region of increased interaction
##
metaplot <- function(regions, signal,
                     signal.random=Reduce(`+`, signal) / length(signal),
                     ext=1e4, bins=200, add.random=TRUE, title="", subtitle="") {
  # extend regions
  start(first(regions)) <- start(first(regions)) - ext
  end(second(regions))  <- end(second(regions)) + ext
  within <- GRanges(seqnames(first(regions)), IRanges(end(first(regions) + 1), start(second(regions) - 1)))

  # bin regions + calculate aggregated average signal per bin
  get_vector <- function(regions, signal) {
    # define regions
    regions  <- keepSeqlevels(regions , names(signal), pruning.mode="coarse")
    regions.binned <- unlist(tile(regions, n=bins))
    regions.binned$region <- rep(seq_along(regions), each=bins)

    # get binned signal average
    x <- binnedAverage(regions.binned, signal, "score")
    x <- matrix(x$score, ncol=bins, byrow=TRUE)
    colSums(x, na.rm=TRUE) / length(regions.binned)   # aggregate signal from multiple regions, normalized to the number of regions
  }

  # get_vector for the 3 subregions: start of the loop, body of the loop, end of the loop
  x <- mclapply(signal, function(signal) {
    x <- mclapply(list(left=first(regions), within=within, right=second(regions)), get_vector, signal)
    c(x$left, x$within, x$right)
  })

  # add some random regions to show the bg
  if(add.random) {
    regions.random <- list(left  =shift(first(regions) , runif(length(first(regions)) , -1e5, -1e4 )),  # shift regions
                           right =shift(second(regions), runif(length(second(regions)), +1e4, +1e5 )))
    regions.random$within <- shift(regions.random$left, -1e4)   # take the contiguous region upstream of the left anchor
    x.random <- mclapply(regions.random, get_vector, signal.random)
    x$random <- c(x.random$left, x.random$within, x.random$right)
  }

  # scale down to have all tracks starting at y=0
  x <- lapply(x, function(x) x - (min(x) - 0))

  # plot something
  df   <- reshape2::melt(x)
  df$x <- rep(1:(3*bins), length(x))
  df$L1<- factor(df$L1, levels=c(names(signal), "random"))
  pal  <- c(palette()[1:length(signal)], "grey30")
  ymin <- min(df$value)
  ymax <- ymin - .05*diff(range(df$value))

  p <- ggplot(df, aes(x=x, y=value, color=L1)) +
         geom_line() +
         geom_vline(xintercept=c(bins-bins/2, bins, 2*bins, 2*bins+bins/2), lty=2) +  # assume `ext` == width(first(regions)) == width(second(regions))
         labs(title=title, subtitle=subtitle, y="normalized signal (arbitrary units)", x="") +
         annotate("rect", xmin=bins-bins/2, xmax=bins         , ymin=ymin, ymax=ymax, alpha=.2, fill="green") +
         annotate("rect", xmin=bins       , xmax=2*bins       , ymin=ymin, ymax=ymax, alpha=.2, fill="grey") +
         annotate("rect", xmin=2*bins     , xmax=2*bins+bins/2, ymin=ymin, ymax=ymax, alpha=.2, fill="green") +
         scale_color_manual("", values=pal) +
         scale_x_continuous(breaks=c(1, bins-bins/2, 3/2*bins, 2*bins+bins/2, 3*bins),
                            labels=c(-ext, "contact start window", "within the loop", "contact end window", ext)) +
         theme(legend.position="bottom")
  p
}

##
## How many loops have cohesin/gapr?
##
pairs_over <- function(x, y, f=`|`) {
  do.call(f, list(first(x) %over% y, second(x) %over% y))
}

df <- data.frame(value=c(sum(!(pairs_over(microc.loops$ctrl, rad21.peaks)) & !(pairs_over(microc.loops$ctrl, gapr.peaks))),
                         sum(  pairs_over(microc.loops$ctrl, rad21.peaks)  & !(pairs_over(microc.loops$ctrl, gapr.peaks))),
                         sum(  pairs_over(microc.loops$ctrl, rad21.peaks)  &   pairs_over(microc.loops$ctrl, gapr.peaks)),
                         sum(!(pairs_over(microc.loops$ctrl, rad21.peaks)) &   pairs_over(microc.loops$ctrl, gapr.peaks))),
                         
                 class=factor(1:4, labels=c("nothing",
                                            "cohesin (Rad21+CTCF)",
                                            "cohesin+supercoiling",
                                            "supercoiling")))

pal <- c("grey", colorRampPalette(palette()[1:2])(3))
ggplot(df, aes(x=0, y=value, fill=class)) +
  geom_bar(stat="identity", position="stack") +
  geom_text(aes(label=value), position="stack", vjust=1) +
  scale_fill_manual("loop contact and", values=pal) +
  scale_x_discrete(breaks=NULL) +
  labs(title="co-ocupancy of cohesin and supercoiling with loops", x="", y="contact loops") +
  theme_bw()

##
## now do metaplots for all loops / loop+cohesin / loop+supercoiling / loop+both / loop+none
##
# scale the 2 signal tracks, otherwise we cannot plot them together
# scaling goes like this: gapr signal is scaled up to the average signal per peak divided by the number of bases on peaks
s <- mean(mean(rad21.signal[rad21.peaks])) / mean(mean(gapr.signal[gapr.peaks])) * sum(width(rad21.peaks)) / sum(width(gapr.peaks))
gapr.signal.scaled <- s * gapr.signal
names(gapr.signal.scaled) <- names(gapr.signal)
s <- mean(mean(rad21.signal[rad21.peaks])) / mean(mean(ctcf.signal[ctcf.peaks])) * sum(width(rad21.peaks)) / sum(width(ctcf.peaks))
ctcf.signal.scaled <- s * ctcf.signal
names(ctcf.signal.scaled) <- names(ctcf.signal)
signal <- list(cohesin=rad21.signal, ctcf=ctcf.signal.scaled, supercoiling=gapr.signal.scaled)
signal.random <- Reduce(`+`, signal) / length(signal)

p1 <- metaplot(microc.loops$ctrl,
         signal, signal.random, title="Cohesin insulation and Supercoiling on MicroC loop contacts",
         subtitle="All loop contacts")
p2 <- metaplot(microc.loops$ctrl[  pairs_over(microc.loops$ctrl, rad21.peaks)  &   pairs_over(microc.loops$ctrl, gapr.peaks) ],
         signal, signal.random, title="Cohesin insulation and Supercoiling on MicroC loop contacts",
         subtitle="Loop contacts +cohesin +supercoiling")
p3 <- metaplot(microc.loops$ctrl[  pairs_over(microc.loops$ctrl, rad21.peaks)  & !(pairs_over(microc.loops$ctrl, gapr.peaks))],
         signal, signal.random, title="Cohesin insulation and Supercoiling on MicroC loop contacts",
         subtitle="Loop contacts +cohesin -supercoiling")
p4 <- metaplot(microc.loops$ctrl[!(pairs_over(microc.loops$ctrl, rad21.peaks)) &   pairs_over(microc.loops$ctrl, gapr.peaks) ],
         signal, signal.random, title="Cohesin insulation and Supercoiling on MicroC loop contacts",
         subtitle="Loop contacts -cohesin +supercoiling")
p5 <- metaplot(microc.loops$ctrl[!(pairs_over(microc.loops$ctrl, rad21.peaks)) & !(pairs_over(microc.loops$ctrl, gapr.peaks))],
         signal, signal.random, title="Cohesin insulation and Supercoiling on MicroC loop contacts",
         subtitle="Loop contacts -cohesin -supercoiling")

plot(p1)
nolabs  <- labs(title="", y="")
noguide <- theme(legend.position="none")
noscale <-scale_x_continuous(breaks=NULL, labels=NULL)
cowplot::plot_grid(p2 + nolabs + noguide + noscale,
                   p3 + nolabs + noguide + noscale,
                   p4 + nolabs + noguide + noscale,
                   p5 + nolabs + noguide + noscale)

##
## plot assessing the significance of GapR enrichemnt at (all) loop anchors
##
# a binnedAverage-like function to show the max, not the average
# why? because the windows are too broad and the signal looks diminished in the boxplot below
binnedView <- function(bins, numvar, varname, na.rm=FALSE, fun=IRanges::viewMeans) {
  if (!is(bins, "GRanges")) 
    stop("'x' must be a GRanges object")
  if (!is(numvar, "RleList")) 
    stop("'numvar' must be an RleList object")
  if (!identical(seqlevels(bins), names(numvar))) 
    stop("'seqlevels(bin)' and 'names(numvar)' must be identical")
  fun2 <- function(v, na.rm = FALSE) {
    if (!isTRUEorFALSE(na.rm)) 
    stop("'na.rm' must be TRUE or FALSE")
    result <- fun(v, na.rm = na.rm)
    w0 <- width(v)
    v1 <- trim(v)
    w1 <- width(v1)
    if (na.rm) {
      na_count <- sum(is.na(v1))
      w0 <- w0 - na_count
      w1 <- w1 - na_count
    }
    result <- result * w1/w0
    result[w0 != 0L & w1 == 0L] <- 0
    result
  }
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  result_list <- lapply(names(numvar), function(seqname) {
    v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
    fun2(v, na.rm = na.rm)
  })
  new_mcol <- unsplit(result_list, as.factor(seqnames(bins)))
  mcols(bins)[[varname]] <- new_mcol
  bins
}

# select the regions 10kb up/down -stream of the loop anchors
regions <- c(GRanges(first(microc.loops$ctrl) , class="left anchor"),
             GRanges(second(microc.loops$ctrl), class="right anchor"))
random.outside <- GRanges(shift(regions, runif(length(microc.loops$ctrl), -1e5, -1e4)), class="random outside loop")
random.inside  <- resize(GRanges(seqnames(first(microc.loops$ctrl)),
                                 IRanges(end(first(microc.loops$ctrl)), start(second(microc.loops$ctrl))),
                                 class="random inside loop"), width=width(first(microc.loops$ctrl)), fix="center")

regions <- c(regions, sample(random.outside, length(microc.loops$ctrl)), random.inside)

# compute average signal inside the range
regions$signal <- binnedView(regions, gapr.signal.scaled, "score", fun=IRanges::viewMaxs)$score
df <- as.data.frame(mcols(regions))
comp <- list(c("left anchor", "random inside loop"),
             c("left anchor", "random outside loop"),
             c("left anchor", "right anchor"),
             c("right anchor", "random inside loop"),
             c("right anchor", "random outside loop"),
             c("random inside loop", "random outside loop"))

ggplot(df, aes(x=class, y=1+signal, color=class)) +
  geom_violin(fill=NA) +
  stat_summary(fun.data=median_hilow, mult=1, geom="pointrange", color="red") +
  ggpubr::stat_compare_means(method="anova")+#, label.y=40) +
  ggpubr::stat_compare_means(label="p.signif", method="t.test", ref.group="left anchor")   +
  scale_color_manual("", values=c(palette()[2:1], palette()[1:2])) +
  scale_y_log10() +
  labs(title="GapR signal in different regions", y="scaled signal (arbitrary unit, log)", x="") +
  theme(legend.position="none")

dev.off()

