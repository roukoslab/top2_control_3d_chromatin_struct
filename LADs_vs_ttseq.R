###########################
##
## Second attempt to relate changes in LADs to changes in expression.
## This time, using our own datasets:
##   -LADs: imb_roukos_2021_30_longo_cutnrun_laminb_top2
##   -expression: sfb_roukos_2022_09_longo_ttseq_top2
##
## ToDo:
##   -hyperG to test if number of DGE genes is overrepresented in diff-LADs.
##    As diffLADs, let's try firts epic-df standard results
##
#############################
library(rtracklayer)
library(venneuler)
library(ggplot2)

FDR <- .01
PROJECT <- "/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2021_01_meta_HCT116_breakome_prediction"
setwd(PROJECT)

pdf("results/LADs_vs_ttseq.pdf")

##
## read data
##
LADs           <- import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.IGG.LADs.bed")
diffLADs.epic  <- import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.TOP2AAID_NoAUX_LaminB1.epic-df.LADs_diff.bed")
diffLADs.epic_250k_merged <- import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.TOP2AAID_NoAUX_LaminB1.epic-df.LADs_diff.250k_merged.bed")
diffLADs.edger <- import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k_edger.10k.bed")
diffLADs.gaps  <- unstrand(import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_NoAUX_LaminB1.vs.TOP2AAID_NoAUX_LaminB1.LADs_diff.bed"))

read.dge <- function(f) {
  x <- read.csv(f)
  x$BH.adjusted.p.values[is.na(x$BH.adjusted.p.values)] <- 1  # correct p-value for non-tested genes
  makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)
}
dge.exons <- read.dge("data/sfb_roukos_2022_09_longo_ttseq_top2/STD_analysis/results/DE_DESeq2/DKO.vs.WT.csv")
dge.exons_introns <- read.dge("data/sfb_roukos_2022_09_longo_ttseq_top2/results/DE_DESeq2/DKO.vs.WT.csv")

##
##   -hyperG to test if number of DGE genes is overrepresented in diff-LADs
##
test_overlap <- function(x, diffLADs, main="") {
  x <- x[x$mean.of.normalized.counts.for.all.samples > 0, ]     # remove not expressed genes

  x$DGE     <- x$BH.adjusted.p.values < FDR
  x$LAD     <- x %over% LADs
  x$diffLAD <- x %over% diffLADs

  # test (only genes on LADs)
  bg <- sum(x$LAD)
  A  <- sum(x$LAD & x$diffLAD)
  B  <- sum(x$LAD & x$DGE)
  AB <- sum(x$LAD & x$diffLAD & x$DGE)
  pval <- sum(dhyper(AB:A, B, bg - B, A))

  # and plot something
  df <- data.frame(diffLAD=c(TRUE, TRUE, FALSE, FALSE),
                   dge    =c(TRUE, FALSE, TRUE, FALSE),
                   value  =c(sum(x$LAD &  x$diffLAD &  x$DGE), # TRUE , TRUE
                             sum(x$LAD &  x$diffLAD & !x$DGE), # TRUE , FALSE
                             sum(x$LAD & !x$diffLAD &  x$DGE), # FALSE, TRUE
                             sum(x$LAD & !x$diffLAD & !x$DGE)))# FALSE, FALSE

  p <- ggplot(df, aes(x=diffLAD, y=value, fill=dge)) +
         geom_bar(stat="identity", position="fill") +
         scale_fill_manual("differentially expressed", values=palette()) +
         labs(title=main,
              subtitle=paste0("p-value=", format.pval(pval)),
              x="differential LAD", y="fraction of genes on LADs") +
         theme_bw() +
         theme(legend.position="bottom")
  print(p)

  # more plots (as a Venn)
  df <- as.data.frame(x[x$LAD])   # select tested genes in the HyperG test
  df <- df[, c("LAD", "DGE", "diffLAD")]
  v  <- venneuler::venneuler(df)
  names(v$labels) <- v$labels
  v$labels["LAD"]     <- paste(bg, "expressed genes\non LADs")
  v$labels["diffLAD"] <- paste(A , "expressed genes\non diff LADs")
  v$labels["DGE"]     <- paste(B , "DGE genes\non LADs")
  plot(v, main=main, sub=paste("genes overlap=", AB, "\np-value=", format.pval(pval)))

  # plot FC of DGE genes on LADs
  df <- as.data.frame(x[x$DGE])
  df2 <- data.frame(LAD    =c(TRUE, TRUE, FALSE, FALSE),    # boxplot annotation with number of genes
                    diffLAD=c(TRUE, FALSE, TRUE, FALSE),
                    label  =c(sum( df$LAD &  df$diffLAD),
                              sum( df$LAD & !df$diffLAD),
                              sum(!df$LAD &  df$diffLAD),
                              sum(!df$LAD & !df$diffLAD)),
                    y      =min(df$log2.fold.change..MLE...group.DKO.vs.WT) - .5)
  df2$label <- paste(df2$label, "genes")

  p <- ggplot(df, aes(x=LAD, y=log2.fold.change..MLE...group.DKO.vs.WT, fill=diffLAD)) +
         geom_boxplot(position=position_dodge(0.8)) +
         ggpubr::stat_compare_means() +
         scale_fill_manual("on differential LAD", values=palette()) +
         geom_text(data=df2, aes(y=y, label=label, color=diffLAD), position=position_dodge(0.8), color="black") +
         labs(title=main, subtitle="Fold Change of differentially expressed genes", x="on LAD", y="log2 FC") +
         theme_bw() +
         theme(legend.position="bottom")
  print(p)

  df
}

# call the function
l <- list("DE genes (exons) on LADs vs. Epic diff LADs"              =list(dge=dge.exons        , diffLADs=diffLADs.epic),
          "DE genes (exons) on LADs vs. Epic 250k diff LADs"         =list(dge=dge.exons        , diffLADs=diffLADs.epic_250k_merged),
         #"DE genes (exons) on LADs vs. EdgeR diff LADs"             =list(dge=dge.exons        , diffLADs=diffLADs.edger),
         #"DE genes (exons) on LADs vs. Epic+EdgeR diff LADs"        =list(dge=dge.exons        , diffLADs=reduce(c(diffLADs.epic, diffLADs.edger))),
          "DE genes (exons) on LADs vs. LAD gaps"                    =list(dge=dge.exons        , diffLADs=diffLADs.gaps),
          "DE genes (exons+introns) on LADs vs. Epic diff LADs"      =list(dge=dge.exons_introns, diffLADs=diffLADs.epic),
          "DE genes (exons+introns) on LADs vs. Epic 250k diff LADs" =list(dge=dge.exons_introns, diffLADs=diffLADs.epic_250k_merged),
         #"DE genes (exons+introns) on LADs vs. EdgeR diff LADs"     =list(dge=dge.exons_introns, diffLADs=diffLADs.edger),
         #"DE genes (exons+introns) on LADs vs. Epic+EdgeR diff LADs"=list(dge=dge.exons_introns, diffLADs=reduce(c(diffLADs.epic, diffLADs.edger))))
          "DE genes (exons+introns) on LADs vs. LAD gaps"            =list(dge=dge.exons_introns, diffLADs=diffLADs.gaps))

res <- Map(l, names(l), f=function(x, x.name) test_overlap(x$dge, x$diffLADs, main=x.name))
names(res) <- gsub("DE genes", "DEG",    # shorten sheet name
              gsub("exons", "e",
              gsub("introns", "i",
              gsub("vs.", "vs",
              gsub("diff LADs", "dLAD", names(res))))))
writexl::write_xlsx(res, "results/LADs_vs_ttseq.xlsx")

dev.off()
