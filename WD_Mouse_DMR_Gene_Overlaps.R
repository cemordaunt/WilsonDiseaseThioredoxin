# WD Mouse DMR Gene Overlaps ####
# Charles Mordaunt
# 4/10/18

setwd()

# Packages ####
library(VennDiagram)
library(ggplot2)
library(reshape2)
library(grid)
library(scales)
library(R.utils)
library(rtracklayer)
library(reshape)
library(devtools)
library(plyr)
library(sm)
library(GeneOverlap)
library(biomaRt)

# Gene Overlap ####
# Genes from GREAT
WTctrl_vs_WTchol_hyper <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_WTchol_Hyper_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTctrl_vs_WTchol_hypo <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_WTchol_Hypo_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTctrl_vs_txJctrl_hyper <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_txJctrl_Hyper_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTctrl_vs_txJctrl_hypo <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_txJctrl_Hypo_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
txJctrl_vs_txJchol_hyper <- sort(unique(as.character(unlist(read.delim("Tables/txJctrl_vs_txJchol_Hyper_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
txJctrl_vs_txJchol_hypo <- sort(unique(as.character(unlist(read.delim("Tables/txJctrl_vs_txJchol_Hypo_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTchol_vs_txJchol_hyper <- sort(unique(as.character(unlist(read.delim("Tables/WTchol_vs_txJchol_Hyper_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTchol_vs_txJchol_hypo <- sort(unique(as.character(unlist(read.delim("Tables/WTchol_vs_txJchol_Hypo_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTctrl_vs_txJchol_hyper <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_txJchol_Hyper_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTctrl_vs_txJchol_hypo <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_txJchol_Hypo_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTchol_vs_txJctrl_hyper <- sort(unique(as.character(unlist(read.delim("Tables/WTchol_vs_txJctrl_Hyper_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTchol_vs_txJctrl_hypo <- sort(unique(as.character(unlist(read.delim("Tables/WTchol_vs_txJctrl_Hypo_DMR_GeneList.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WTctrl_vs_WTchol_all <- sort(unique(c(WTctrl_vs_WTchol_hyper, WTctrl_vs_WTchol_hypo)))
WTctrl_vs_txJctrl_all <- sort(unique(c(WTctrl_vs_txJctrl_hyper, WTctrl_vs_txJctrl_hypo)))
txJctrl_vs_txJchol_all <- sort(unique(c(txJctrl_vs_txJchol_hyper, txJctrl_vs_txJchol_hypo)))
WTchol_vs_txJchol_all <- sort(unique(c(WTchol_vs_txJchol_hyper, WTchol_vs_txJchol_hypo)))
WTctrl_vs_txJchol_all <- sort(unique(c(WTctrl_vs_txJchol_hyper, WTctrl_vs_txJchol_hypo)))
WTchol_vs_txJctrl_all <- sort(unique(c(WTchol_vs_txJctrl_hyper, WTchol_vs_txJctrl_hypo)))

all_genes <- sort(unique(c(WTctrl_vs_WTchol_all, WTctrl_vs_txJctrl_all, txJctrl_vs_txJchol_all, WTchol_vs_txJchol_all, WTctrl_vs_txJchol_all, WTchol_vs_txJctrl_all)))

gene_list <- list(WTctrl_vs_WTchol_all, WTctrl_vs_WTchol_hyper, WTctrl_vs_WTchol_hypo,
                  WTctrl_vs_txJctrl_all, WTctrl_vs_txJctrl_hyper, WTctrl_vs_txJctrl_hypo,
                  txJctrl_vs_txJchol_all, txJctrl_vs_txJchol_hyper, txJctrl_vs_txJchol_hypo, 
                  WTchol_vs_txJchol_all, WTchol_vs_txJchol_hyper, WTchol_vs_txJchol_hypo, 
                  WTctrl_vs_txJchol_all, WTctrl_vs_txJchol_hyper, WTctrl_vs_txJchol_hypo, 
                  WTchol_vs_txJctrl_all, WTchol_vs_txJctrl_hyper, WTchol_vs_txJctrl_hypo)
names(gene_list) <- c("WTctrl_vs_WTchol_all", "WTctrl_vs_WTchol_hyper", "WTctrl_vs_WTchol_hypo",
                      "WTctrl_vs_txJctrl_all", "WTctrl_vs_txJctrl_hyper", "WTctrl_vs_txJctrl_hypo",
                      "txJctrl_vs_txJchol_all", "txJctrl_vs_txJchol_hyper", "txJctrl_vs_txJchol_hypo", 
                      "WTchol_vs_txJchol_all", "WTchol_vs_txJchol_hyper", "WTchol_vs_txJchol_hypo", 
                      "WTctrl_vs_txJchol_all", "WTctrl_vs_txJchol_hyper", "WTctrl_vs_txJchol_hypo", 
                      "WTchol_vs_txJctrl_all", "WTchol_vs_txJctrl_hyper", "WTchol_vs_txJctrl_hypo")

# Choline Gene Overlap Venns
venn.diagram(list("WTctrl_vs_txJctrl"=WTctrl_vs_txJctrl_hyper, "WTctrl_vs_txJchol"=WTctrl_vs_txJchol_hyper),
             "Figures/WD Mouse Liver DMR Choline Hyper Gene Overlap.png", height=8, width=10, imagetype="png", units="in", fontfamily="sans", cat.fontfamily="sans", 
             fill=c("lightblue", "lightpink"), cex=3, lwd=4, cat.cex=2, cat.pos=c(0,0), cat.dist=c(0.05, 0.05), rotation.degree=180, ext.text=FALSE)
# "Cebpb"  "Egflam" "Gdnf"   "Ptpn1" 

venn.diagram(list("WTctrl_vs_txJctrl"=WTctrl_vs_txJctrl_hypo, "WTctrl_vs_txJchol"=WTctrl_vs_txJchol_hypo),
             "Figures/WD Mouse Liver DMR Choline Hypo Gene Overlap.png", height=8, width=10, imagetype="png", units="in", fontfamily="sans", cat.fontfamily="sans", 
             fill=c("lightblue", "lightpink"), cex=3, lwd=4, cat.cex=2, cat.pos=c(0,0), cat.dist=c(0.05, 0.05), rotation.degree=180, ext.text=FALSE)

intersect(WTctrl_vs_txJctrl_hyper, WTchol_vs_txJchol_hyper)
#[1] "Rpl38" "Sdk2" 

intersect(WTctrl_vs_txJctrl_hypo, WTchol_vs_txJchol_hypo)
# None

# Thioredoxin gene overlap ####
txn_genes <- read.delim("Tables/thioredoxin_gene_list_mm10.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
txn_genes <- sort(unique(as.character(unlist(txn_genes))))
intersect_txn_gene_list <- matrix(nrow = length(gene_list), ncol = 1, dimnames=list(names(gene_list), "txn_genes"))
for(i in 1:length(gene_list)){
        intersect_txn_gene_list[i,1] <- length(intersect(gene_list[[i]], txn_genes))
}
rownames(intersect_txn_gene_list)[intersect_txn_gene_list[,1] > 0]
# "txJctrl_vs_txJchol_all"  "txJctrl_vs_txJchol_hypo" 
# "WTchol_vs_txJchol_all"   "WTchol_vs_txJchol_hyper" 
# "WTctrl_vs_txJchol_all"   "WTctrl_vs_txJchol_hyper"

intersect(gene_list$txJctrl_vs_txJchol_hypo, txn_genes) # Prdx2
intersect(gene_list$WTchol_vs_txJchol_hyper, txn_genes) # Gpx4
intersect(gene_list$WTctrl_vs_txJchol_hyper, txn_genes) # Hif1a

# Mouse and Human DMR Gene Overlap ####
# Mouse genes
names(gene_list)
# "WTctrl_vs_WTchol_all"     "WTctrl_vs_WTchol_hyper"   "WTctrl_vs_WTchol_hypo"    "WTctrl_vs_txJctrl_all"    "WTctrl_vs_txJctrl_hyper" 
# "WTctrl_vs_txJctrl_hypo"   "txJctrl_vs_txJchol_all"   "txJctrl_vs_txJchol_hyper" "txJctrl_vs_txJchol_hypo"  "WTchol_vs_txJchol_all"   
# "WTchol_vs_txJchol_hyper"  "WTchol_vs_txJchol_hypo"   "WTctrl_vs_txJchol_all"    "WTctrl_vs_txJchol_hyper"  "WTctrl_vs_txJchol_hypo"  
# "WTchol_vs_txJctrl_all"    "WTchol_vs_txJctrl_hyper"  "WTchol_vs_txJctrl_hypo"  

# Human genes
# WD liver specific DMRs all, hyper, hypo
WD_liver_hyper <- sort(unique(as.character(unlist(read.delim("C:/Users/Booboo/Charles/Documents/Programming/Wilson's Disease Liver/GREAT/WD Liver Hyper DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WD_liver_hypo <- sort(unique(as.character(unlist(read.delim("C:/Users/Booboo/Charles/Documents/Programming/Wilson's Disease Liver/GREAT/WD Liver Hypo DMR GREAT Distal Genes.txt", sep="\t", header=FALSE, stringsAsFactors=FALSE)))))
WD_liver_all <- sort(unique(c(WD_liver_hyper, WD_liver_hypo)))
WD_liver_gene_list <- list(WD_liver_all, WD_liver_hyper, WD_liver_hypo)
names(WD_liver_gene_list) <- c("WD_liver_all", "WD_liver_hyper", "WD_liver_hypo")

# Convert from hg19 to mouse mm10 genes with ensembl biomart
WD_liver_mouse_hyper <- sort(unique(as.character(unlist(read.delim("Tables/WD Liver Hyper DMR Genes to Mouse.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)))))
WD_liver_mouse_hypo <- sort(unique(as.character(unlist(read.delim("Tables/WD Liver Hypo DMR Genes to Mouse.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)))))
WD_liver_mouse_list <- list(WD_liver_mouse_hyper, WD_liver_mouse_hypo)
names(WD_liver_mouse_list) <- c("WD_liver_mouse_hyper", "WD_liver_mouse_hypo")

# Mouse Background genes from GREAT
WTctrl_vs_txJctrl_back <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_txJctrl_background_genes.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,1]))))
WTctrl_vs_WTchol_back <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_WTchol_background_genes.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,1]))))
txJctrl_vs_txJchol_back <- sort(unique(as.character(unlist(read.delim("Tables/txJctrl_vs_txJchol_background_genes.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,1]))))
WTchol_vs_txJchol_back <- sort(unique(as.character(unlist(read.delim("Tables/WTchol_vs_txJchol_background_genes.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,1]))))
WTctrl_vs_txJchol_back <- sort(unique(as.character(unlist(read.delim("Tables/WTctrl_vs_txJchol_background_genes.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,1]))))
WTchol_vs_txJctrl_back <- sort(unique(as.character(unlist(read.delim("Tables/WTchol_vs_txJctrl_background_genes.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,1]))))
mouse_background_genes <- sort(unique(c(WTctrl_vs_txJctrl_back, WTctrl_vs_WTchol_back, txJctrl_vs_txJchol_back, WTchol_vs_txJchol_back, WTctrl_vs_txJchol_back, WTchol_vs_txJctrl_back)))

# Get Genes
gom <- newGOM(gene_list, WD_liver_mouse_list, genome.size = length(mouse_background_genes)) # genome.size = union of all mouse background genes
intersects_genes <- getNestedList(gom, "intersection")

WD_liver_mouse_hyper_intersect <- intersects_genes$WD_liver_mouse_hyper
WD_liver_mouse_hyper_intersect <- lapply(WD_liver_mouse_hyper_intersect, sort)
maxLength <- max(sapply(WD_liver_mouse_hyper_intersect, length))
WD_liver_mouse_hyper_intersect <- lapply(WD_liver_mouse_hyper_intersect, function(x){c(x, rep(NA, maxLength - length(x)))})
WD_liver_mouse_hyper_intersect <- as.data.frame(WD_liver_mouse_hyper_intersect)
write.table(WD_liver_mouse_hyper_intersect, "Tables/WD Liver to Mouse hyper DMR gene Mouse DMR gene overlap.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

WD_liver_mouse_hypo_intersect <- intersects_genes$WD_liver_mouse_hypo
WD_liver_mouse_hypo_intersect <- lapply(WD_liver_mouse_hypo_intersect, sort)
maxLength <- max(sapply(WD_liver_mouse_hypo_intersect, length))
WD_liver_mouse_hypo_intersect <- lapply(WD_liver_mouse_hypo_intersect, function(x){c(x, rep(NA, maxLength - length(x)))})
WD_liver_mouse_hypo_intersect <- as.data.frame(WD_liver_mouse_hypo_intersect)
write.table(WD_liver_mouse_hypo_intersect, "Tables/WD Liver to Mouse hypo DMR gene Mouse DMR gene overlap.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Human WD Liver and Mouse liver DMR Gene Overlap ####
# Overlap
gom <- newGOM(WD_liver_gene_list_to_mouse, gene_list, genome.size = length(mouse_background_genes)) # genome.size = union of all mouse background genes
oddsRatio <- getMatrix(gom, "odds.ratio") 
pValue <- getMatrix(gom, "pval")
fdr <- matrix(p.adjust(pValue, method="fdr"), nrow=nrow(pValue), ncol=ncol(pValue), byrow=FALSE)
dimnames(fdr) <- dimnames(pValue)
intersects <- getMatrix(gom, "intersection")
overlapResults <- as.data.frame(cbind(intersects, oddsRatio, pValue, fdr))
colnames(overlapResults) <- paste(rep(c("intersect", "odds_ratio", "p_value", "q_value"), each=ncol(intersects)), colnames(overlapResults), sep="_")
write.table(overlapResults, "Tables/Mouse Human WD Liver DMR Gene Overlap Results.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

intersects_genes <- getNestedList(gom, "intersection")
intersects_genes$WTctrl_vs_txJctrl_hyper$WD_liver_hyper
# "Alkbh5"   "Myo15"    "Mlxipl"   "Cebpb"    "Ptpn1"    "Grid2ip"  "Zfp750"   "Kdelr2"   "Mink1"    "Tmem184a" "Sh3pxd2a"
intersects_genes$WTctrl_vs_txJctrl_hypo$WD_liver_hypo
# "B4galnt4" "Ykt6"     "Camk2b"   "Lonrf2"   "Afap1"    "Aff3"     "Wnt11"    "Zfp644"   "Cacna1h"  "Psapl1"   "Sept9"    "Fam129b"  "Trim11"   "Cacnb1"  
# "Mad1l1"   "Ppcdc"    "Myt1l"    "Prkcz"    "Hlcs"     "Tab1"     "Elfn1"    "Tpsg1"    "Sim2"    

# logqvalue Heatmap hyper hypo
logqValue <- -log10(fdr)
logqValue <- melt(logqValue)
colnames(logqValue) <- c("Human", "Mouse", "logqValue")
logqValue <- subset(logqValue, grepl("hyp", Human) & grepl("hyp", Mouse))

logqValue$Mouse <- factor(logqValue$Mouse, levels=rev(c("WTctrl_vs_txJctrl_hyper","WTctrl_vs_txJctrl_hypo",
                                                        "WTctrl_vs_WTchol_hyper","WTctrl_vs_WTchol_hypo",
                                                        "txJctrl_vs_txJchol_hyper","txJctrl_vs_txJchol_hypo",
                                                        "WTchol_vs_txJchol_hyper","WTchol_vs_txJchol_hypo",
                                                        "WTctrl_vs_txJchol_hyper","WTctrl_vs_txJchol_hypo",
                                                        "WTchol_vs_txJctrl_hyper","WTchol_vs_txJctrl_hypo")), ordered=TRUE)
logqValue$Human <- as.character(logqValue$Human)
logqValue$Human[logqValue$Human == "WD_liver_hyper"] <- "Hyper"
logqValue$Human[logqValue$Human == "WD_liver_hypo"] <- "Hypo"
logqValue$Human <- factor(logqValue$Human, levels=c("Hyper", "Hypo"), ordered=TRUE)

gg <- ggplot(data = logqValue)
gg +
        geom_tile(aes(x = Human, y = Mouse, fill = logqValue)) +
        scale_fill_gradientn("-log(q-value)\n", colors = c("Black", "#FF0000"), values = c(0,1), na.value = "#FF0000", limits=c(0,max(logqValue$logqValue))) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25), 
              axis.ticks = element_line(size = 1), legend.key = element_blank(), legend.text=element_text(size=24),
              panel.grid.minor = element_blank(), legend.position = c(1.35, 0.805), legend.key.size = unit(2, "lines"),
              legend.background = element_blank(),
              plot.margin = unit(c(1,11,1,7), "lines"), 
              axis.text.x = element_text(size = 20, color = "Black", angle = 0, hjust = 0.5, vjust = 0.5),
              axis.text.y = element_blank(),
              axis.title.x = element_text(size=24), 
              axis.title.y=element_blank(),
              legend.title = element_text(size = 24)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))
ggsave("Figures/Mouse Human WD Liver Hyper Hypo DMR Gene Overlap Heatmap.png", dpi = 600, width = 7, height = 8, units = "in")