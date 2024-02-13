###########################
##
## Second attempt to relate changes in LADs to changes in chromatin conformation and expression.
## This time, using our own datasets:
##   -LADs: imb_roukos_2021_30_longo_cutnrun_laminb_top2
##   -expression: sfb_roukos_2022_09_longo_ttseq_top2
##   -microC: ../data/micro_c_top2_mutants/
##
## ToDo:
##   -hyperG to test if number of DGE genes is overrepresented in diff-LADs.
##    As diffLADs, let's simple differential DKO-WT LADs instead of epic-df results
##
#############################
library(rtracklayer)
library(venneuler)
library(grid)
library(VennDiagram)
library(ggplot2)
library(biomaRt)

FDR <- .01
PROJECT <- "/fsimb/groups/imb-bioinfocf/projects/sfb/roukos/sfb_roukos_2021_01_meta_HCT116_breakome_prediction"
setwd(PROJECT)

pdf("results/LADs_vs_microC.pdf")

##
## read data
##

# LADs
LADs <- reduce(c(import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.IGG.LADs.bed"),
                 import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/TOP2AAID_NoAUX_LaminB1.vs.IGG.LADs.bed")))
#diffLADs.epic  <- import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.TOP2AAID_NoAUX_LaminB1.epic-df.LADs_diff.bed")
#diffLADs.epic_250k_merged <- import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.TOP2AAID_NoAUX_LaminB1.epic-df.LADs_diff.250k_merged.bed")
#diffLADs.edger <- import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k_edger.10k.bed")
diffLADs.gaps  <- unstrand(import.bed("data/imb_roukos_2021_30_longo_cutnrun_laminb_top2/results/epic_1k/DOUBLE_PlusAUX_LaminB1.vs.TOP2AAID_NoAUX_LaminB1.LADs_diff.only_1e5_diffs.bed"))

# ttseq
ann <- getBM(attributes=c("ensembl_gene_id", "gene_biotype"),
             filters   ="biotype",
             values    ="protein_coding",
             mart      =useEnsembl(biomart="ensembl", 
                                   dataset="hsapiens_gene_ensembl"))

read.dge <- function(f) {
  x <- read.csv(f)
  x$BH.adjusted.p.values[is.na(x$BH.adjusted.p.values)] <- 1  # correct p-value for non-tested genes
  x$gene_id <- sub("\\.\\d+$", "", x$gene_id)
  makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)
}
#dge.exons <- read.dge("data/sfb_roukos_2022_09_longo_ttseq_top2/STD_analysis/results/DE_DESeq2/DKO.vs.WT.csv")
dge.exons_introns <- read.dge("data/sfb_roukos_2022_09_longo_ttseq_top2/results/DE_DESeq2/DKO.vs.WT.csv")

# microC
f <- list(up  ="data/micro_c_top2_mutants/analysis_papantonis/INS2951_INS2948_250_w10_g1.0.9up.bed",
          down="data/micro_c_top2_mutants/analysis_papantonis/INS2951_INS2948_250_w10_g1.0.02down.bed")
microC <- unlist(as(Class="GRangesList", object=Map(f, names(f), f=function(x, x.name) {
  x <- import.bed(x)
  x$direction <- x.name
  x
})))
names(microC) <- NULL

##
##   -hyperG to test if number of DGE genes is overrepresented in diff-LADs
##
test_overlap <- function(dge, diffLADs, main="") {
  x <- dge[dge$mean.of.normalized.counts.for.all.samples > 0]   # remove not expressed genes
  #criteria <- "expressed"
  x <- dge[dge$gene_id %in% ann$ensembl_gene_id]  # take all protein coding genes
  criteria <- "protein coding"
  #x <- dge[dge$mean.of.normalized.counts.for.all.samples > 0 & dge$gene_id %in% ann$ensembl_gene_id]   # both (expressed + protein coding)
  #criteria <- "expressed & protein coding"

  x$DGE     <- x$BH.adjusted.p.values < FDR
  x$LAD     <- x %over% LADs
  x$diffLAD <- x %over% diffLADs
  x$microC  <- x %over% microC
  x$DGE_sign     <- ifelse(x$DGE & x$log2.fold.change..MLE...group.DKO.vs.WT < 0, "down",
                    ifelse(x$DGE & x$log2.fold.change..MLE...group.DKO.vs.WT > 0, "up", "non-DGE"))
  x$diffLAD_sign <- ifelse(x %over% diffLADs[diffLADs$score > 0], "gain",
                    ifelse(x %over% diffLADs[diffLADs$score < 0], "loss",
                    ifelse(x$LAD, "no-change", "non-LAD")))
  x$microC_sign  <- ifelse(x %over% microC[microC$direction == "up"], "up",
                    ifelse(x %over% microC[microC$direction == "down"], "down", "no-change"))

  # basic plot comparing average expression of DKO vs. WT in gain/loss of LADs
  x$DKO <- log2(1 + rowMeans(as.matrix(mcols(x)[, grepl("^DKO", colnames(mcols(x)))])))
  x$WT  <- log2(1 + rowMeans(as.matrix(mcols(x)[, grepl("^WT", colnames(mcols(x)))])))
  plot(x$DKO ~ x$WT, main="average expression of DKO vs. WT in gain/loss of LADs", xlab="WT", ylab="DKO", type="n")
  points(x$DKO[x$diffLAD_sign == "gain"] ~ x$WT[x$diffLAD_sign == "gain"], col="red")
  points(x$DKO[x$diffLAD_sign == "loss"] ~ x$WT[x$diffLAD_sign == "loss"], col="blue")
  abline(0, 1)
  legend("bottomright", fill=c("red", "blue"), legend=c("gene on LAD gain in DKO", "gene on LAD loss in DKO"))

  # or put it differently: is DKO-WT expression different in gained/lost LADs?
  # I tried all imaginable things, but differences between DKO and WT are mostly inexistent
  # finally, came up with this silly barplot showing the fraction of genes with expression increased by 1/3 (33%) in DKO
  xx <- table(x$DKO - x$WT > log2(1 + 1/3), x$diffLAD_sign)
  df <- reshape2::melt(apply(xx, 2, function(x) x[2] / sum(x)))  # calculate fraction of upregulated genes
  p  <- fisher.test(xx[, c("gain", "loss")])$p.value   # test if fraction of upregulated genes is significant in LAD loss compared gain
  pval <- data.frame(group1="gain", group2="loss", p=p, p.format=format.pval(p))
  p <- ggplot(df, aes(x=rownames(df), y=value)) +
         geom_bar(aes(fill=rownames(df)), stat="identity") +
         ggpubr::stat_pvalue_manual(pval, y.position=1.05*max(df$value), label="fisher-test p={p.format}") +
         labs(title="Fraction of genes with expression increased by 1/3 (33%) in DKO", x="", y="") +
         scale_fill_manual("", values=c("red", "blue", "grey", "grey"), guide=FALSE)
  print(p)

  # test (only genes on LADs)
  bg <- sum(x$LAD)
  A  <- sum(x$LAD & x$diffLAD)
  B  <- sum(x$LAD & x$DGE)
  AB <- sum(x$LAD & x$diffLAD & x$DGE)
  pval <- sum(dhyper(AB:A, B, bg - B, A))

  # plot as an Euler
  df <- as.data.frame(x[x$LAD])   # select tested genes in the HyperG test
  df <- df[, c("LAD", "DGE", "diffLAD")]
  v  <- venneuler::venneuler(df)
  names(v$labels) <- v$labels
  v$labels["LAD"]     <- paste(bg, paste(criteria, "genes\non LADs"))
  v$labels["diffLAD"] <- paste(A , paste(criteria, "genes\non diff LADs"))
  v$labels["DGE"]     <- paste(B , "DGE genes\non LADs")
  plot(v, main=main, sub=paste("genes overlap=", AB, "\np-value=", format.pval(pval)))

  # plot as a Venn, separating up/down genes, gain/loss LADs
  df <- x[x$LAD]   # select tested genes in the HyperG test
  df <- list(gained_LAD=df$gene_id[df$diffLAD_sign == "gain"],
             lost_LAD  =df$gene_id[df$diffLAD_sign == "loss"],
             DE_up     =df$gene_id[df$DGE_sign == "up"],
             DE_down   =df$gene_id[df$DGE_sign == "down"])
  n <- c(paste(criteria, "genes\non gained LADs"),
         paste(criteria, "genes\non lost LADs"),
         "significantly upregulated",
         "significantly downregulated")
  grid.newpage()
  grid.draw(VennDiagram::venn.diagram(df, filename=NULL, fill=palette()[1:4], alpha=.3, category.names=n))

  # again, as an Euler (seems to not work, doesn't overlap some of the circles which should indeed overlap)
  df <- x[x$LAD]   # select tested genes in the HyperG test
  df <- data.frame(#LAD=TRUE,
                   gained_LAD=df$diffLAD_sign == "gain",
                   lost_LAD  =df$diffLAD_sign == "loss",
                   DE_up     =df$DGE_sign == "up",
                   DE_down   =df$DGE_sign == "down")
  v  <- venneuler::venneuler(df)
  names(v$labels) <- v$labels
  plot(v)

  # plot FC of DGE genes on LADs, separating by gain/loss lads
  df <- reshape2::melt(table(x$diffLAD_sign[x$DGE], x$DGE_sign[x$DGE]))
  p <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
         geom_bar(stat="identity", position="fill") +
         geom_text(aes(label=as.character(value)), position="fill") +
         scale_fill_manual("FC in DKO vs. WT", values=palette()) +
         labs(title="fraction of DE genes on LADs", x="LAD changes", y="fraction of DE genes") +
         theme_bw() +
         theme(legend.position="bottom")
  print(p)

  # plot differential LADs vs. differential microC
  df <- list(gain=list(up  =sum(width(intersect(diffLADs[diffLADs$score > 0], microC[microC$direction == "up"]))),
                       down=sum(width(intersect(diffLADs[diffLADs$score > 0], microC[microC$direction == "down"])))),
             loss=list(up  =sum(width(intersect(diffLADs[diffLADs$score < 0], microC[microC$direction == "up"]))),
                       down=sum(width(intersect(diffLADs[diffLADs$score < 0], microC[microC$direction == "down"])))))
  df <- reshape2::melt(df) 
  df$value <- df$value / 1e6
  
  p <- ggplot(df, aes(x=L1, y=value, fill=L2)) +
         geom_bar(stat="identity") +
         scale_fill_manual("microC", values=palette()) +
         labs(title="Overlap between changes in LADs and regions of increased interaction",
              subtitle=paste("Significance (chisqr test):", format.pval(chisq.test(matrix(df$value, ncol=2))$p.value)),
              x="LAD changes", y="Millions of bases") +
         theme_bw() +
         theme(legend.position="bottom")
  print(p)

  # plot FC of DGE genes on regions of microC increased interaction, separating by up/down
  df <- reshape2::melt(table(x$microC_sign[x$DGE], x$DGE_sign[x$DGE]))
  p <- ggplot(df, aes(x=Var1, y=value, fill=Var2)) +
         geom_bar(stat="identity", position="fill") +
         geom_text(aes(label=as.character(value)), position="fill") +
         scale_fill_manual("FC in DKO vs. WT", values=palette()) +
         labs(title="fraction of DE genes on microC differential regions", x="microC changes", y="fraction of DE genes") +
         theme_bw() +
         theme(legend.position="bottom")
  print(p)

  # plot as an Euler
  x$criteria <- TRUE
  df <- as.data.frame(x[x$criteria])   # select tested genes in the HyperG test
  df <- df[, c("criteria", "DGE", "microC")]
  bg <- sum(x$criteria)
  A  <- sum(x$criteria & x$microC)
  B  <- sum(x$criteria & x$DGE)
  AB <- sum(x$criteria & x$microC & x$DGE)
  pval <- sum(dhyper(AB:A, B, bg - B, A))

  v  <- venneuler::venneuler(df)
  names(v$labels) <- v$labels
  v$labels["criteria"] <- paste(bg, paste(criteria, "genes"))
  v$labels["microC"]   <- paste(A , paste(criteria, "genes\non diff microC"))
  v$labels["DGE"]      <- paste(B , "DE genes")
  plot(v, main=main, sub=paste("genes overlap=", AB, "\np-value=", format.pval(pval)))
}

# call the function
l <- list(#"DE genes (exons) on LADs vs. Epic diff LADs"              =list(dge=dge.exons        , diffLADs=diffLADs.epic),
          #"DE genes (exons) on LADs vs. Epic 250k diff LADs"         =list(dge=dge.exons        , diffLADs=diffLADs.epic_250k_merged),
          #"DE genes (exons) on LADs vs. EdgeR diff LADs"             =list(dge=dge.exons        , diffLADs=diffLADs.edger),
          #"DE genes (exons) on LADs vs. Epic+EdgeR diff LADs"        =list(dge=dge.exons        , diffLADs=reduce(c(diffLADs.epic, diffLADs.edger))),
          #"DE genes (exons) on LADs vs. LAD gaps"                    =list(dge=dge.exons        , diffLADs=diffLADs.gaps),
          #"DE genes (exons+introns) on LADs vs. Epic diff LADs"      =list(dge=dge.exons_introns, diffLADs=diffLADs.epic),
          #"DE genes (exons+introns) on LADs vs. Epic 250k diff LADs" =list(dge=dge.exons_introns, diffLADs=diffLADs.epic_250k_merged),
          #"DE genes (exons+introns) on LADs vs. EdgeR diff LADs"     =list(dge=dge.exons_introns, diffLADs=diffLADs.edger),
          #"DE genes (exons+introns) on LADs vs. Epic+EdgeR diff LADs"=list(dge=dge.exons_introns, diffLADs=reduce(c(diffLADs.epic, diffLADs.edger))))
          "DE genes (exons+introns) on LADs vs. LAD gaps"            =list(dge=dge.exons_introns, diffLADs=diffLADs.gaps))

res <- Map(l, names(l), f=function(x, x.name) test_overlap(x$dge, x$diffLADs, main=x.name))

dev.off()
