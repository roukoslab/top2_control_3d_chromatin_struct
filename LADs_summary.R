###########################
##
## So far, very simple script to describe some global stats for LADs:
##   - number and length of the LADs
##   - bases covered
##   - bases gained/lost upon TOP2 DKO
##
#############################
library(rtracklayer)
library(ggplot2)

FDR <- .01
PROJECT <- "/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2021_01_meta_HCT116_breakome_prediction"
setwd(PROJECT)

pdf("results/LADs_summary.pdf")

##
## read data
##

# LADs
LADs <- list(dko=import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.IGG.LADs.bed"),
             wt =import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/TOP2AAID_NoAUX_LaminB1.vs.IGG.LADs.bed"))
diffLADs.gaps <- unstrand(import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.TOP2AAID_NoAUX_LaminB1.LADs_diff.only_1e5_diffs.bed"))

##
## Plots
##

##   - bases covered
df <- data.frame(value=c(sum(width(LADs$dko)), sum(width(LADs$wt))) / 1e6,
                 class=c("DKO", "WT"))

ggplot(df, aes(x=class, y=value, fill=class)) +
  geom_bar(stat="identity") +
  scale_fill_manual("", values=palette()) +
  labs(title="Total bases covered", x="", y="bases covered (millions)") +
  theme_classic() +
  theme(legend.position="none")

##   - number and length of the LADs
df <- data.frame(value=c(length(LADs$dko), length(LADs$wt)),
                 class=c("DKO", "WT"))

ggplot(df, aes(x=class, y=value, fill=class)) +
  geom_bar(stat="identity") +
  scale_fill_manual("", values=palette()) +
  labs(title="Number of LADs detected", x="", y="") +
  theme_classic() +
  theme(legend.position="none")

df <- data.frame(value=c(width(LADs$dko), width(LADs$wt)),
                 class=c(rep("DKO", length(LADs$dko)), rep("WT", length(LADs$wt))))

ggplot(df, aes(x=class, y=value, color=class)) +
  ggbeeswarm::geom_beeswarm() +
  geom_violin(fill=NA) +
  stat_summary(fun.data=median_hilow, geom="pointrange", color="red") +
  ggpubr::stat_compare_means(label.x=1.5) +
  scale_color_manual("", values=palette()) +
  scale_y_log10() +
  labs(title="Length of LADs", y="bases covered (log10)", x="") +
  theme_classic() +
  theme(legend.position="none")

##   - bases gained/lost upon TOP2 DKO
df <- data.frame(value=sapply(split(width(diffLADs.gaps), diffLADs.gaps$score > 0), sum) / 1e6,
                 class=c("lost", "gained"))   # lost is score == -1, gained score == +1. Gaps were calculated as `dko - wt`

ggplot(df, aes(x=class, y=value, fill=class)) +
  geom_bar(stat="identity") +
  scale_fill_manual("", values=palette()) +
  labs(title="Bases gained/lost upon TOP2 DKO", x="", y="bases covered (millions)") +
  theme_classic() +
  theme(legend.position="none")

dev.off()
