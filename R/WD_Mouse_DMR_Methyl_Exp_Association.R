# WD Mouse Liver DMR Methylation and Gene Expression Association ####
# Charles Mordaunt
# 5/3/18

setwd()

# Packages ####
library(reshape2)
library(ggplot2)
library(scales)

# Functions ####
combineFiles <- function(chroms, prefix, suffix){
# Combine chromosome-level files together
        DMRs <- NULL
        for(i in 1:length(chroms)){
                temp <- NULL
                if(file.exists(paste(prefix,chroms[i],suffix, sep=""))){
                        temp <- read.delim(paste(prefix,chroms[i],suffix, sep=""), header = TRUE, sep = "\t", stringsAsFactors=FALSE)
                        DMRs <- rbind(DMRs, temp)
                }
        }
        DMRs
}

exp_anova <- function(formula){
# Run ANOVA stats on RNA-seq gene expression data
        means <- aggregate(formula, FUN=mean)
        sds <- aggregate(formula, FUN=sd)
        test_aov <- aov(formula)
        p_aov <- summary(test_aov)[[1]]$'Pr(>F)'[1]
        test_tukey <- TukeyHSD(test_aov)
        stats <- c(means[,2], sds[,2], p_aov, test_tukey[[1]][1,], test_tukey[[1]][2,], test_tukey[[1]][3,], test_tukey[[1]][4,], test_tukey[[1]][5,], test_tukey[[1]][6,])
        names(stats) <- c(paste("mean", as.character(means[,1]), sep="_"), paste("sd", as.character(means[,1]), sep="_"), "ANOVA_p",
                          paste(rep(dimnames(test_tukey[[1]])[[1]], each=4), c("diff", "lwr", "upr", "padj"), sep="_"))
        stats
}

# Data ####
# Methylation Data
chroms <- c(paste("chr",1:19,sep=""), "chrM")
conc_dmr_meth <- combineFiles(chroms=chroms, prefix="Concatenated_DMR_Methyl/Concatenated_DMR_Methyl_", suffix = "_DMR_methylation.txt")
table(is.na(conc_dmr_meth$VMNS001A))
# FALSE  TRUE 
# 642    30 
# 30 DMRs failed coverage requirement of 1 read in all samples, removed
conc_dmr_meth <- subset(conc_dmr_meth, !is.na(VMNS001A))
conc_dmr_meth$DMRid <- as.character(paste("DMR", 1:nrow(conc_dmr_meth), sep="_"))
conc_dmr_meth_bed <- conc_dmr_meth[,c("chr", "start", "end", "DMRid")]
#write.table(conc_dmr_meth_bed, "UCSC Tracks/Concatenated DMRs with ID for GREAT.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# DMR Genes
conc_dmr_genes <- read.delim("GREAT/Concatenated DMR Genes by Gene.txt", sep="\t", header=FALSE, stringsAsFactors = FALSE)
conc_dmr_genes <- conc_dmr_genes[2:nrow(conc_dmr_genes), 1:2]
colnames(conc_dmr_genes) <- c("Gene", "DMR")
DMRid <- strsplit(conc_dmr_genes$DMR, " ")
table(sapply(DMRid, length))
#   2   4   6   8 
# 835 140  17   5 
DMRid <- lapply(DMRid, function(x){
        if(length(x) == 2){return(x[1])}
        if(length(x) == 4){return(c(x[1], x[3]))}
        if(length(x) == 6){return(c(x[1], x[3], x[5]))}
        if(length(x) == 8){return(c(x[1], x[3], x[5], x[7]))}})
names(DMRid) <- conc_dmr_genes$Gene
DMRid_m <- melt(DMRid)
colnames(DMRid_m) <- c("DMR", "Gene")
DMRid_m$DMR <- as.character(DMRid_m$DMR)
conc_dmr_meth_gene <- merge(DMRid_m, conc_dmr_meth, by.x="DMR", by.y="DMRid", all.x=TRUE, all.y=FALSE)
colnames(conc_dmr_meth_gene)[6:21]
# "VMNS001A" "VMNS002A" "VMNS003A" "VMNS004A" "VMNS001B" "VMNS002B" "VMNS003B" "VMNS004B" "VMNS001C" "VMNS002C" "VMNS003C" "VMNS004C" "VMNS001D" "VMNS002D" "VMNS003D" "VMNS004D"
colnames(conc_dmr_meth_gene)[6:21] <- c(paste("WTctrl", 1:4, sep=""), paste("txJctrl", 1:4, sep=""), paste("WTchol", 1:4, sep=""), paste("txJchol", 1:4, sep=""))
dmr_genes <- sort(unique(as.character(conc_dmr_meth_gene$Gene)))
rm(conc_dmr_genes, conc_dmr_meth, conc_dmr_meth_bed, DMRid_m, chroms, DMRid)

# Expression Data
WTctrl_vs_txJctrl_exp <- read.delim("RNAseq/WTctrl_vs_txJctrl_all_RNAseq.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
WTctrl_vs_WTchol_exp <- read.delim("RNAseq/WTctrl_vs_WTchol_all_RNAseq.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
txJctrl_vs_txJchol_exp <- read.delim("RNAseq/txJctrl_vs_txJchol_all_RNAseq.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
WTchol_vs_txJchol_exp <- read.delim("RNAseq/WTchol_vs_txJchol_all_RNAseq.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)

table(duplicated(WTctrl_vs_txJctrl_exp$Gene.Name)) # TRUE 78, Duplicate genes, don't know which transcript is which

table(WTctrl_vs_txJctrl_exp$Gene.Name == WTctrl_vs_WTchol_exp$Gene.Name &
              WTctrl_vs_txJctrl_exp$Gene.Name == txJctrl_vs_txJchol_exp$Gene.Name &
              WTctrl_vs_txJctrl_exp$Gene.Name == WTchol_vs_txJchol_exp$Gene.Name)
#  TRUE 41128, Gene names are all in the same order

table(WTctrl_vs_txJctrl_exp$WTctrl1 == WTctrl_vs_WTchol_exp$WTctrl1) #TRUE
table(WTctrl_vs_txJctrl_exp$txJctrl1 == txJctrl_vs_txJchol_exp$txJctrl1) #TRUE
table(WTctrl_vs_WTchol_exp$WTchol1 == WTchol_vs_txJchol_exp$WTchol1) #TRUE
table(txJctrl_vs_txJchol_exp$txJchol1 == WTchol_vs_txJchol_exp$txJchol1) #TRUE
# RPKMs are the same between different comparisons

table(WTctrl_vs_txJctrl_exp$Gene.Name == WTchol_vs_txJchol_exp$Gene.Name) # TRUE 
all_exp <- cbind(WTctrl_vs_txJctrl_exp, WTchol_vs_txJchol_exp[,2:9])
dmr_exp <- subset(all_exp, Gene.Name %in% dmr_genes)
table(duplicated(dmr_exp$Gene.Name)) # TRUE 1
dmr_exp$Gene.Name[duplicated(dmr_exp$Gene.Name)] # "Map3k4"
dmr_exp[dmr_exp$Gene.Name == "Map3k4",]
#       Gene.Name    WTctrl1    WTctrl2   WTctrl3    WTctrl4  txJctrl1   txJctrl2   txJctrl3   txJctrl4   WTchol1    WTchol2   WTchol3   WTchol4   txJchol1   txJchol2    txJchol3   txJchol4
# 25652    Map3k4 4.34525285 4.18282996 4.2076879 4.61754160 3.4275772 4.12563213 4.39071730 4.13132798 5.1177298 4.58474934 4.7052378 4.9896934 3.40648621 4.37041953  3.72316266 4.08745489
# 25653    Map3k4 0.05013771 0.07066345 0.1489366 0.02964279 0.1322626 0.03803318 0.08943846 0.05429337 0.2949503 0.03845071 0.1313294 0.1159047 0.04929328 0.06938399  0.04585701 0.01888994
# Use Map3k4 with higher exp
dmr_exp <- dmr_exp[!row.names(dmr_exp) == "25653",]
table(duplicated(dmr_exp$Gene.Name)) # FALSE 

zero_test <- apply(dmr_exp[,2:ncol(dmr_exp)], 1, function(x){sum(x > 0)}) # How many non-zero values
table(zero_test)
#  0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16 
# 81  16  22  12   9   8  11   8   7   7  10   9   7   8  13  28 721 
dmr_exp <- subset(dmr_exp, zero_test >= 4) # Remove genes with less than 4 non-zero values
rm(all_exp, txJctrl_vs_txJchol_exp, WTchol_vs_txJchol_exp, WTctrl_vs_txJctrl_exp, WTctrl_vs_WTchol_exp, zero_test)

# Associate Methylation and Expression ####
colnames(conc_dmr_meth_gene) <- paste("meth", colnames(conc_dmr_meth_gene), sep="_")
colnames(dmr_exp) <- paste("exp", colnames(dmr_exp), sep="_")
dmr_exp[,2:ncol(dmr_exp)] <- lapply(dmr_exp[,2:ncol(dmr_exp)], function(x){log2(x+1)}) # log2 transform, add 1 to RPKM to prevent -Inf values

meth_exp <- merge(conc_dmr_meth_gene, dmr_exp, by.x="meth_Gene", by.y="exp_Gene.Name", all=FALSE)
colnames(meth_exp)[1:5] <- c("Gene", "DMR", "chr", "start", "end")
meth <- as.matrix(meth_exp[,6:21])
exp <- as.matrix(meth_exp[,22:37])

meth_exp_stats <- NULL
for(i in 1:nrow(meth_exp)){
        if(i %% 50 == 0){cat(paste(rownames(meth_exp)[i], "\t", sep=""))}
        assoc <- NULL
        stats <- NULL
        assoc <- summary(lm(exp[i,] ~ meth[i,]))
        stats <- as.character(assoc$coefficients[2,1:4])
        meth_exp_stats <- rbind(meth_exp_stats, stats)
}
meth_exp_stats <- as.data.frame(meth_exp_stats, stringsAsFactors=FALSE)
meth_exp_stats <- as.data.frame(apply(meth_exp_stats, 2, as.numeric))
colnames(meth_exp_stats) <- c("estimate", "std_error", "tstat", "pvalue")
meth_exp <- cbind(meth_exp, meth_exp_stats)
meth_exp$qvalue <- p.adjust(meth_exp$pvalue, "fdr")
meth_exp$top <- meth_exp$pvalue < 0.05
meth_exp$top <- factor(meth_exp$top, levels=c("TRUE", "FALSE"), ordered=TRUE)
table(meth_exp$top) # 67 associations
meth_exp$logp <- -log10(meth_exp$pvalue)
write.table(meth_exp, "Tables/DMR Methylation Expression Association Results.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Add DMR comparison to DMR-Gene Associations ####
dmrs <- read.delim("Tables/Concatenated_DMRs_WD_Mouse_Liver.bed", sep="\t", header=FALSE, stringsAsFactors = FALSE)
colnames(dmrs) <- c("chr", "start", "end", "Comparison")
assoc <- read.delim("Tables/DMR Methylation Expression Association Results.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
assoc <- merge(assoc, dmrs, by=c("chr", "start", "end"), all.x=TRUE, all.y=FALSE)
write.table(assoc, "Tables/DMR Methylation Expression Association Results with Comparison.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names = FALSE)

# Meth and Expression Volcano Plot ####
meth_exp_label <- subset(meth_exp, pvalue < 0.005)
gg <- ggplot()
gg +
        geom_point(data=meth_exp, aes(x = estimate, y = logp, color = top), size = 2.5) +
        geom_text(data=meth_exp_label, aes(x = estimate, y=logp, label=Gene), color="black", hjust=0, nudge_x=0.2, size=7, check_overlap = TRUE, fontface="italic") +
        annotate(geom="text", x=c(0.45, 1.1), y=c(3.98, 3.75), label="*", size=12, hjust=0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(2, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.85, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black"),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 24),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 26), 
              legend.title = element_text(color = "black", size = 24)) +
        scale_color_manual(name = "DMR-Gene\nAssociation", values = c("TRUE" = "#FF3366", "FALSE" = "#3366CC")) +
        xlab(expression(log[2]*"(RPKM+1) / Methylation (%)")) +
        ylab(expression(-log[10]*"(p-value)")) +
        coord_cartesian(xlim=c(-6,6), ylim=c(0,4.5)) +
        scale_y_continuous(expand=c(0.01,0.01))
ggsave("Figures/Meth Exp Association Volcano Plot.png", dpi = 600, width = 8, height = 8, units = "in")

# Sept9 Methylation and Expression
sept9 <- subset(meth_exp, Gene == "Sept9" & qvalue < 0.05)
sept9 <- data.frame(meth = unlist(sept9[1,grep("meth", colnames(sept9))]),
                    exp = unlist(sept9[1,grep("exp", colnames(sept9))]))
sept9$group <- factor(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), levels=c("WTctrl", "txJctrl", "WTchol", "txJchol"), ordered=TRUE)

gg <- ggplot(sept9)
gg +
        geom_smooth(aes(x=meth, y=exp), color="black", method="lm", se=FALSE) +
        geom_point(aes(x = meth, y = exp, color=group), size = 4) +
        annotate(geom="text", x=0.435, y=3.8, label="Sept9 ~ DMR_405", size=10) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(2, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.85, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black", vjust=0),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 24),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 26), 
              legend.title = element_blank()) +
        xlab("Methylation (%)") +
        ylab(expression(log[2]*"(RPKM+1)")) +
        coord_cartesian(ylim=c(3,3.8)) +
        scale_x_continuous(labels = function(x)x*100) +
        scale_color_manual(breaks = c("WTctrl", "txJctrl", "WTchol", "txJchol"), 
                           values = c("WTctrl"="#3366CC", "txJctrl"="#FF3366", "WTchol"="#009933","txJchol"="#FF66CC"))
ggsave("Figures/Sept9 Methylation and Expression Plot.png", dpi = 600, width = 8, height = 8, units = "in")

# Tbc1d1 Methylation and Expression
Tbc1d1 <- subset(meth_exp, Gene == "Tbc1d1" & qvalue < 0.05)
Tbc1d1 <- Tbc1d1[1,] # Remove Duplicate
Tbc1d1 <- data.frame(meth = unlist(Tbc1d1[1,grep("meth", colnames(Tbc1d1))]),
                    exp = unlist(Tbc1d1[1,grep("exp", colnames(Tbc1d1))]))
Tbc1d1$group <- factor(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), levels=c("WTctrl", "txJctrl", "WTchol", "txJchol"), ordered=TRUE)

gg <- ggplot(Tbc1d1)
gg +
        geom_smooth(aes(x=meth, y=exp), color="black", method="lm", se=FALSE) +
        geom_point(aes(x = meth, y = exp, color=group), size = 4) +
        annotate(geom="text", x=0.45, y=2.05, label="Tbc1d1 ~ DMR_164", size=10, hjust=0) +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              legend.key = element_blank(), legend.key.size = unit(2, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(0.85, 0.87), 
              legend.background = element_blank(), legend.text = element_text(size = 24, color = "Black", vjust=0),
              plot.margin = unit(c(1,1,1,1), "lines"), 
              axis.text = element_text(color = "black", size = 24),
              axis.ticks = element_line(color = "black", size = 1.25),
              axis.title = element_text(color = "black", size = 26), 
              legend.title = element_blank()) +
        xlab("Methylation (%)") +
        ylab(expression(log[2]*"(RPKM+1)")) +
        coord_cartesian(ylim=c(1.57,2.05)) +
        scale_x_continuous(labels = function(x)x*100) +
        scale_color_manual(breaks = c("WTctrl", "txJctrl", "WTchol", "txJchol"), 
                           values = c("WTctrl"="#3366CC", "txJctrl"="#FF3366", "WTchol"="#009933","txJchol"="#FF66CC"))
ggsave("Figures/Tbc1d1 Methylation and Expression Plot.png", dpi = 600, width = 8, height = 8, units = "in")

# Meth and Expression Heatmap ####
meth_exp_heat <- subset(meth_exp, top == "TRUE", select=-c(chr, start, end))
meth_exp_heat <- meth_exp_heat[order(meth_exp_heat$pvalue),]
meth_exp_heat$Assoc <- paste("Assoc", 1:nrow(meth_exp_heat), sep="_")
meth_exp_heat <- meth_exp_heat[,c(42, 1:34)]

meth_heat <- meth_exp_heat[,grep("meth", colnames(meth_exp_heat))]
meth_means <- rowMeans(meth_heat)
meth_heat <- (meth_heat - meth_means)*100
meth_exp_heat[,grep("meth", colnames(meth_exp_heat))] <- meth_heat

exp_heat <- meth_exp_heat[,grep("exp", colnames(meth_exp_heat))]
exp_means <- rowMeans(exp_heat)
exp_heat <- exp_heat - exp_means
meth_exp_heat[,grep("exp", colnames(meth_exp_heat))] <- exp_heat

meth_exp_heat <- melt(meth_exp_heat, id.vars=c("Assoc", "Gene", "DMR"))
meth_exp_heat$data <- rep(NA, nrow(meth_exp_heat))
meth_exp_heat$data[grepl("meth", meth_exp_heat$variable)] <- "Methylation"
meth_exp_heat$data[grepl("exp", meth_exp_heat$variable)] <- "Expression"
meth_exp_heat$data <- factor(meth_exp_heat$data, levels=c("Methylation", "Expression"), ordered=TRUE)
meth_exp_heat$variable <- as.character(meth_exp_heat$variable)
meth_exp_heat$variable <- gsub(".*_", "", meth_exp_heat$variable)
meth_exp_heat[,c("Assoc", "Gene", "DMR")] <- lapply(meth_exp_heat[,c("Assoc", "Gene", "DMR")], as.factor)
meth_exp_heat$variable <- factor(meth_exp_heat$variable, levels=c(paste(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), 1:4, sep="")), ordered=TRUE)
row_order <- hclust(dist(as.matrix(meth_heat)), "ward.D")$order # Cluster DMRs by methylation

# Run exp_anova by group on expression of top genes to subset genes actually changing by group, and correlating with DMR methylation
pheno <- data.frame("Sample"=paste(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), 1:4, sep=""),
                    "Group"=rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4))
exp_stats <- cbind(exp_heat, t(apply(exp_heat, 1, function(x){exp_anova(x ~ pheno$Group)})))
exp_stats$Assoc <- paste("Assoc", 1:nrow(exp_stats), sep="_")
exp_stats$Gene <- meth_exp_heat$Gene[match(exp_stats$Assoc, meth_exp_heat$Assoc)]
write.table(exp_stats, "Tables/DMR-Gene Expression ANOVA Stats.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

table(exp_stats$ANOVA_p < 0.05) # 25 with main effect of group on expression
exp_heat_ANOVA <- subset(exp_stats, ANOVA_p < 0.05)
exp_heat_ANOVA <- exp_heat_ANOVA[,c(grep("exp", colnames(exp_heat_ANOVA)), ncol(exp_heat_ANOVA)-1, ncol(exp_heat_ANOVA))]
exp_heat_ANOVA <- subset(exp_heat_ANOVA, !duplicated(Gene)) # Remove genes with more than one association, same expression values

# Expression Heatmap, clustered by expression, ANOVA p < 0.05
exp_heat_ANOVA_m <- melt(exp_heat_ANOVA, id.vars=c("Assoc", "Gene"))
row_order <- hclust(dist(as.matrix(exp_heat_ANOVA[,1:16])), "ward.D")$order # Cluster genes by expression
exp_heat_ANOVA_m$variable <- gsub(".*_", "", as.character(exp_heat_ANOVA_m$variable))
exp_heat_ANOVA_m$variable <- factor(exp_heat_ANOVA_m$variable, levels=c(paste(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), 1:4, sep="")), ordered=TRUE)
exp_heat_ANOVA_m$Gene <- factor(exp_heat_ANOVA_m$Gene, levels=unique(exp_heat_ANOVA_m$Gene)[row_order], ordered=TRUE)

gg <- ggplot(data = exp_heat_ANOVA_m)
gg +
        geom_tile(aes(x = variable, y = Gene, fill = value)) +
        scale_fill_gradientn(expression(log[2]*"(RPKM+1)"), colors = c("#0000FF", "black", "#FF0000"), values = c(0,0.7,1), na.value = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(color="black", size=1), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.23, 0.86), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,9,1,1), "lines"), 
              axis.text.y = element_text(color="black", size=18, hjust=1, face="italic"),
              axis.text.x=element_text(color="black", size=18, angle=90, hjust=1, vjust=0.5),
              axis.title = element_blank(), 
              legend.title = element_text(size = 18)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))
ggsave("Figures/Methylation Expression Associations Expression Heatmap ANOVA genes.png", dpi = 600, width = 8, height = 7, units = "in")
