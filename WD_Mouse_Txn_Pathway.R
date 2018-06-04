# Thioredoxin Pathway Analysis ####
# Charles Mordaunt
# 5/14/18

setwd()

# Packages ####
library(ggplot2)
library(scales)
library(reshape2)
library(dplyr)

# Functions ####
meth_anova <- function(formula){
# Run ANOVA stats on mouse liver promoter methylation data
        means <- aggregate(formula, FUN=mean)
        sds <- aggregate(formula, FUN=sd)
        test_aov <- aov(formula)
        p_aov <- summary(test_aov)[[1]]$'Pr(>F)'[3]
        test_tukey <- TukeyHSD(test_aov)
        stats <- c(means[,3], sds[,3], test_tukey[[1]], test_tukey[[2]], p_aov, test_tukey[[3]][1,], test_tukey[[3]][2,], test_tukey[[3]][3,], test_tukey[[3]][4,],
                   test_tukey[[3]][5,], test_tukey[[3]][6,])
        names(stats) <- c(paste("mean", paste(as.character(means[,1]), as.character(means[,2]), sep=""), sep="_"), 
                          paste("sd", paste(as.character(means[,1]), as.character(means[,2]), sep=""), sep="_"),
                          paste(rep(c("Strain", "Treatment"), each=4), c("diff", "lwr", "upr", "padj"), sep="_"), "p_Strain_Treatment",
                          paste(rep(c(dimnames(test_tukey[[3]])[[1]]), each=4), c("diff", "lwr", "upr", "padj"), sep="_"))
        stats
}

human_meth_ttest <- function(formula){
# Run ttest on human liver promoter methylation data
        test <- t.test(formula)
        means <- as.numeric(test$estimate)
        diff <- means[1] - means[2]
        conf <- as.numeric(test$conf.int)
        tstat <- as.numeric(test$statistic)
        pval <- as.numeric(test$p.value)
        stats <- as.numeric(c(means, diff, conf, tstat, pval))
        names(stats) <- c("mean_WD", "mean_HC", "mean_diff", "conf_lwr", "conf_upr", "tstat", "pvalue")
        stats
}

exp_anova <- function(formula){
# Run ANOVA stats on RNAseq gene expression data
        means <- aggregate(formula, FUN=mean)
        sds <- aggregate(formula, FUN=sd)
        test_aov <- aov(formula)
        p_aov <- summary(test_aov)[[1]]$'Pr(>F)'[3]
        test_tukey <- TukeyHSD(test_aov)
        stats <- c(means[,3], sds[,3], test_tukey[[1]], test_tukey[[2]], p_aov, test_tukey[[3]][1,], test_tukey[[3]][2,], test_tukey[[3]][3,], test_tukey[[3]][4,],
                   test_tukey[[3]][5,], test_tukey[[3]][6,])
        names(stats) <- c(paste("mean", paste(as.character(means[,1]), as.character(means[,2]), sep=""), sep="_"), 
                          paste("sd", paste(as.character(means[,1]), as.character(means[,2]), sep=""), sep="_"),
                          paste(rep(c("Strain", "Treatment"), each=4), c("diff", "lwr", "upr", "padj"), sep="_"), "p_Strain_Treatment",
                          paste(rep(c(dimnames(test_tukey[[3]])[[1]]), each=4), c("diff", "lwr", "upr", "padj"), sep="_"))
        stats
}

# Data ####
mouse_permeth <- read.delim("Tables/WD_Mouse_Thioredoxin_Promoter_AvgMeth2col.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
mouse_2col <- read.delim("Tables/WD_Mouse_Thioredoxin_Promoter_AvgMeth2col.txt.2col", sep="\t", header=TRUE, stringsAsFactors = FALSE)
# human thioredoxin genes mapped to mouse with OrthoRetriever http://lighthouse.ucsf.edu/orthoretriever/
# promoter regions -5kb to +1kb
# methylation data if region has at least 1 read in 1 sample

# Check Coverage ####
coverage <- mouse_2col[,c(4, grep("total", colnames(mouse_2col)))]
minreads <- apply(coverage[,2:ncol(coverage)], 1, min)
summary(minreads)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.0    56.0    92.0   100.9   133.0   286.0 
table(minreads < 10)
# FALSE  TRUE 
# 40     1 
# Gm10639
rm(coverage, minreads, mouse_2col)

# Stats ####
meth <- mouse_permeth[,grep("VMNS", colnames(mouse_permeth))]
pheno <- data.frame("Sample"=colnames(meth), "Strain"=c(rep("WT", 8), rep("txJ", 8)), 
                    "Treatment"=c(rep("ctrl", 4), rep("chol", 4), rep("ctrl", 4), rep("chol", 4)))
meth_stats <- cbind(mouse_permeth[,c("Chromosome", "Start", "End", "Name")], t(apply(meth, 1, function(x){meth_anova(x ~ pheno$Strain * pheno$Treatment)})))
colnames(meth_stats)[grep("p", colnames(meth_stats))]
#"Strain_padj" "Treatment_padj" "p_Strain_Treatment" "WT:chol-txJ:chol_padj" "txJ:ctrl-txJ:chol_padj" "WT:ctrl-txJ:chol_padj" "txJ:ctrl-WT:chol_padj" "WT:ctrl-WT:chol_padj" "WT:ctrl-txJ:ctrl_padj" 
meth_stats$Name[meth_stats$Strain_padj < 0.05] #"Gstp1"  "Hif1an"
meth_stats$Name[meth_stats$Treatment_padj < 0.05] #"Txn2"
meth_stats$Name[meth_stats$p_Strain_Treatment < 0.05] #"Gm10639"
meth_stats$Name[meth_stats$'WT:chol-txJ:chol_padj' < 0.05] # None
meth_stats$Name[meth_stats$'txJ:ctrl-txJ:chol_padj' < 0.05] # None
meth_stats$Name[meth_stats$'WT:ctrl-txJ:chol_padj' < 0.05] # None
meth_stats$Name[meth_stats$'txJ:ctrl-WT:chol_padj' < 0.05] # None
meth_stats$Name[meth_stats$'WT:ctrl-WT:chol_padj' < 0.05] # None
meth_stats$Name[meth_stats$'WT:ctrl-txJ:ctrl_padj' < 0.05] # None

meth_stats$Strain_fdrq <- p.adjust(meth_stats$Strain_padj, method="fdr")
meth_stats$Treatment_fdrq <- p.adjust(meth_stats$Treatment_padj, method="fdr")
meth_stats$Strain_Treatment_fdrq <- p.adjust(meth_stats$p_Strain_Treatment, method="fdr")

meth_stats$Name[meth_stats$Strain_fdrq < 0.05] # None
meth_stats$Name[meth_stats$Treatment_fdrq < 0.05] # None
meth_stats$Name[meth_stats$Strain_Treatment_fdrq < 0.05] # None
write.table(meth_stats, file="Tables/Mouse Thioredoxin Pathway Genes Promoter Methylation Stats.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Methylation Heatmap ####
meth_heat <- mouse_permeth[,4:ncol(mouse_permeth)]
colnames(meth_heat)[2:ncol(meth_heat)] <- c(paste("WTctrl", 1:4, sep=""), paste("WTchol", 1:4, sep=""), paste("txJctrl", 1:4, sep=""), paste("txJchol", 1:4, sep=""))
meth_means <- rowMeans(meth_heat[,2:ncol(meth_heat)])
meth_heat[,2:ncol(meth_heat)] <- (meth_heat[,2:ncol(meth_heat)] - meth_means)*100
row_order <- hclust(dist(as.matrix(meth_heat[,2:ncol(meth_heat)])), "ward.D")$order # Cluster DMRs by methylation
meth_heat <- melt(meth_heat, id.vars="Name")
meth_heat$variable <- factor(meth_heat$variable, levels=c(paste(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), 1:4, sep="")), ordered=TRUE)
meth_heat$Name <- factor(meth_heat$Name, levels=unique(meth_heat$Name)[row_order], ordered=TRUE)

gg <- ggplot(data = meth_heat)
gg +
        geom_tile(aes(x = variable, y = Name, fill = value)) +
        scale_fill_gradientn("Promoter\nMethylation (%)", colors = c("#0000FF", "black", "#FF0000"), values = c(0,0.55,1), na.value = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(color="black", size=1), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.25, 0.88), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,10,1,1), "lines"), 
              axis.text.y = element_text(color="black", size=12, hjust=1, face="italic"),
              axis.text.x=element_text(color="black", size=12, angle=90, hjust=1, vjust=0.5),
              axis.title=element_blank(), legend.title = element_text(size = 17)) +
        scale_x_discrete(expand=c(0,0))
ggsave("Figures/WD Mouse Thioredoxin Pathway Promoter Methylation Heatmap.png", dpi = 600, width = 7, height = 8, units = "in")

# Methylation Heatmap Top ####
top_genes <- c("Gpx1", "Gpx4", "Gstcd", "Gstk1", "Gsto1", "Gstp1", "Hif1an", "Mgst3", "Prdx1", "Prdx2", "Prdx3", "Prdx4", "Txn1", "Txn2", "Txnrd1", "Txnrd2")
mouse_permeth_top <- subset(mouse_permeth, Name %in% top_genes)
meth_heat <- mouse_permeth_top[,4:ncol(mouse_permeth_top)]
colnames(meth_heat)[2:ncol(meth_heat)] <- c(paste("WTctrl", 1:4, sep=""), paste("WTchol", 1:4, sep=""), paste("txJctrl", 1:4, sep=""), paste("txJchol", 1:4, sep=""))
meth_means <- rowMeans(meth_heat[,2:ncol(meth_heat)])
meth_heat[,2:ncol(meth_heat)] <- (meth_heat[,2:ncol(meth_heat)] - meth_means)*100
meth_heat <- melt(meth_heat, id.vars="Name")
meth_heat$variable <- factor(meth_heat$variable, levels=c(paste(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), 1:4, sep="")), ordered=TRUE)
meth_heat$Name <- factor(meth_heat$Name, levels=levels(rel_exp$Name), ordered=TRUE)

gg <- ggplot(data = meth_heat)
gg +
        geom_tile(aes(x = variable, y = Name, fill = value)) +
        scale_fill_gradientn("Promoter\nMethylation (%)", colors = c("#0000FF", "black", "#FF0000"), values = c(0,0.25,1), na.value = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(color="black", size=1), legend.key = element_blank(), legend.key.size=unit(1.75, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(1.26, 0.8), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,10,1,1), "lines"), 
              axis.text.y = element_text(color="black", size=16, hjust=1, face="italic"),
              axis.text.x=element_text(color="black", size=16, angle=90, hjust=1, vjust=0.5),
              axis.title=element_blank(), legend.title = element_text(size = 17)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))
ggsave("Figures/WD Mouse Top Thioredoxin Pathway Promoter Methylation Heatmap.png", dpi = 600, width = 7, height = 7, units = "in")

# Mouse Thioredoxin Pathway Gene Expression ####
# Data ####
WTctrl_vs_txJctrl_exp <- read.delim("RNAseq/WTctrl_vs_txJctrl_all_RNAseq.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
WTctrl_vs_WTchol_exp <- read.delim("RNAseq/WTctrl_vs_WTchol_all_RNAseq.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
txJctrl_vs_txJchol_exp <- read.delim("RNAseq/txJctrl_vs_txJchol_all_RNAseq.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
WTchol_vs_txJchol_exp <- read.delim("RNAseq/WTchol_vs_txJchol_all_RNAseq.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
all_exp <- cbind(WTctrl_vs_txJctrl_exp, WTchol_vs_txJchol_exp[,2:9])
mouse_txn_genes <- sort(unique(as.character(mouse_permeth$Name)))
txn_exp <- subset(all_exp, Gene.Name %in% mouse_txn_genes)
table(duplicated(txn_exp$Gene.Name)) # All FALSE
zero_test <- apply(txn_exp[,2:ncol(txn_exp)], 1, function(x){sum(x > 0)}) # How many non-zero values
table(zero_test)
# 0  2 12 15 16 
# 4  2  1  1 33 
txn_exp <- subset(txn_exp, zero_test > 0) # Remove genes with no expression
txn_exp[,2:ncol(txn_exp)] <- lapply(txn_exp[,2:ncol(txn_exp)], function(x){log2(x+1)}) # log2 transform, add 1 to RPKM to prevent -Inf values
rm(all_exp, txJctrl_vs_txJchol_exp, WTchol_vs_txJchol_exp, WTctrl_vs_txJctrl_exp, WTctrl_vs_WTchol_exp, zero_test)

# Stats ####
exp <- txn_exp[,2:ncol(txn_exp)]
pheno <- data.frame("Sample"=paste(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), 1:4, sep=""),
                    "Strain"=rep(c("WT", "txJ", "WT", "txJ"), each=4),
                    "Treatment"=rep(c("ctrl", "chol"), each=8))
exp_stats <- cbind(txn_exp, t(apply(exp, 1, function(x){exp_anova(x ~ pheno$Strain * pheno$Treatment)})))
exp_stats$Strain_q <- p.adjust(exp_stats$Strain_padj, method="fdr")
exp_stats$Treatment_q <- p.adjust(exp_stats$Treatment_padj, method="fdr")
exp_stats$Strain_Treatment_q <- p.adjust(exp_stats$p_Strain_Treatment, method="fdr")
exp_stats$`txJ:ctrl-txJ:chol_q` <- p.adjust(exp_stats$`txJ:ctrl-txJ:chol_padj`, method="fdr")
exp_stats$`WT:chol-txJ:chol_q` <- p.adjust(exp_stats$`WT:chol-txJ:chol_padj`, method="fdr")
exp_stats$`WT:ctrl-txJ:chol_q` <- p.adjust(exp_stats$`WT:ctrl-txJ:chol_padj`, method="fdr")
exp_stats$`txJ:ctrl-WT:chol_q` <- p.adjust(exp_stats$`txJ:ctrl-WT:chol_padj`, method="fdr")
exp_stats$`WT:ctrl-txJ:ctrl_q` <- p.adjust(exp_stats$`WT:ctrl-txJ:ctrl_padj`, method="fdr")
exp_stats$`WT:ctrl-WT:chol_q` <- p.adjust(exp_stats$`WT:ctrl-WT:chol_padj`, method="fdr")
write.table(exp_stats, "Tables/Mouse Txn Gene Expression ANOVA Stats.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Expression Heatmap ####
exp_means <- rowMeans(exp)
rel_exp <- exp - exp_means
row_order <- hclust(dist(as.matrix(rel_exp)), "ward.D")$order # Cluster genes by expression
rel_exp$Name <- exp_stats$Gene.Name
rel_exp <- melt(rel_exp, id.vars="Name")
rel_exp$variable <- factor(rel_exp$variable, levels=c(paste(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), 1:4, sep="")), ordered=TRUE)
rel_exp$Name <- factor(rel_exp$Name, levels=unique(rel_exp$Name)[row_order], ordered=TRUE)

gg <- ggplot(data = rel_exp)
gg +
        geom_tile(aes(x = variable, y = Name, fill = value)) +
        scale_fill_gradientn(expression(log[2]*"(Expression)"), colors = c("#0000FF", "black", "#FF0000"), values = c(0,0.6,1), na.value = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(color="black", size=1), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.25, 0.88), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,10,1,1), "lines"), 
              axis.text.y = element_text(color="black", size=12, hjust=1, face="italic"),
              axis.text.x=element_text(color="black", size=12, angle=90, hjust=1, vjust=0.5),
              axis.title=element_blank(), legend.title = element_text(size = 17)) +
        scale_x_discrete(expand=c(0,0))
ggsave("Figures/WD Mouse Thioredoxin Pathway Gene Expression Heatmap.png", dpi = 600, width = 7, height = 8, units = "in")

# Expression Heatmap Top ####
top_genes <- c("Gpx1", "Gpx4", "Gstcd", "Gstk1", "Gsto1", "Gstp1", "Hif1an", "Mgst3", "Prdx1", "Prdx2", "Prdx3", "Prdx4", "Txn1", "Txn2", "Txnrd1", "Txnrd2")
exp_top <- subset(exp, exp_stats$Gene.Name %in% top_genes)
exp_means <- rowMeans(exp_top)
rel_exp <- exp_top - exp_means
exp_row_order <- hclust(dist(as.matrix(rel_exp)), "ward.D")$order # Cluster genes by expression
rel_exp$Name <- exp_stats$Gene.Name[exp_stats$Gene.Name %in% top_genes]
rel_exp <- melt(rel_exp, id.vars="Name")
rel_exp$variable <- factor(rel_exp$variable, levels=c(paste(rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4), 1:4, sep="")), ordered=TRUE)
rel_exp$Name <- factor(rel_exp$Name, levels=unique(rel_exp$Name)[exp_row_order], ordered=TRUE)

gg <- ggplot(data = rel_exp)
gg +
        geom_tile(aes(x = variable, y = Name, fill = value)) +
        scale_fill_gradientn(expression(log[2]*"(RPKM+1)"), colors = c("#0000FF", "black", "#FF0000"), values = c(0,0.6,1), na.value = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(color="black", size=1), legend.key = element_blank(), legend.key.size=unit(1.75, "lines"),
              panel.grid.minor = element_blank(), legend.position = c(1.25, 0.82), 
              legend.background = element_blank(),
              plot.margin = unit(c(1,10,1,1), "lines"), 
              axis.text.y = element_text(color="black", size=16, hjust=1, face="italic"),
              axis.text.x=element_text(color="black", size=16, angle=90, hjust=1, vjust=0.5),
              axis.title=element_blank(), legend.title = element_text(size = 17)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0))
ggsave("Figures/WD Mouse Top Thioredoxin Pathway Gene Expression Heatmap.png", dpi = 600, width = 7, height = 7, units = "in")

# Human WD Liver Thioredoxin Pathway Gene Methylation ####

# Data ####
human_permeth <- read.delim("Tables/WD_Human_Thioredoxin_Promoter_AvgMeth2col.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
human_2col <- read.delim("Tables/WD_Human_Thioredoxin_Promoter_AvgMeth2col.txt.2col", sep="\t", header=TRUE, stringsAsFactors = FALSE)

# Check Coverage ####
coverage <- human_2col[,c(4, grep("total", colnames(human_2col)))]
minreads <- apply(coverage[,2:ncol(coverage)], 1, min)
summary(minreads)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 15.0    70.0   115.0   138.3   170.0   492.0 
table(minreads < 10) # All False
rm(coverage, minreads, human_2col)

# Stats ####
meth <- human_permeth[,grep("VMDK", colnames(human_permeth))]
pheno <- data.frame("Sample"=colnames(meth), "Diagnosis"=c(rep("Healthy_control", 6), rep("Wilson_disease", 10), rep("Disease_control", 5)))
meth <- meth[,!pheno$Diagnosis == "Disease_control"]
pheno <- subset(pheno, !Diagnosis == "Disease_control")
pheno$Diagnosis <- factor(as.character(pheno$Diagnosis), levels=c("Wilson_disease", "Healthy_control"), ordered=TRUE)
human_meth_stats <- cbind(human_permeth[,c("Chromosome", "Start", "End", "Name")], t(apply(meth, 1, function(x){human_meth_ttest(x ~ pheno$Diagnosis)})))
human_meth_stats$fdr_qvalue <- p.adjust(human_meth_stats$pvalue, method="fdr")

human_meth_stats$Name[human_meth_stats$pvalue < 0.05] # "PRDX1" "GSTM2" "GSTM5" "GSTM5" "GSTM5" "MGST3" "GSTP1" "MGST1" "GSTZ1" "GSTZ1" "TXN2"  "TXN2"  "TXN2" 
human_meth_stats$Name[human_meth_stats$fdr_qvalue < 0.05] # None
write.table(human_meth_stats, file="Tables/Human Thioredoxin Pathway Genes Promoter Methylation Stats No DC.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
human_meth_sig <- subset(human_meth_stats, pvalue < 0.05)

# Dot plots ####
# GSTM2
GSTM2 <- cbind(pheno, as.numeric(meth[human_meth_stats$Name == "GSTM2" & human_meth_stats$pvalue < 0.05,]))
GSTM2$Diagnosis <- as.character(GSTM2$Diagnosis)
GSTM2$Diagnosis[GSTM2$Diagnosis == "Healthy_control"] <- "Healthy\nControl"
GSTM2$Diagnosis[GSTM2$Diagnosis == "Wilson_disease"] <- "Wilson\nDisease"
GSTM2$Diagnosis <- factor(GSTM2$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)
colnames(GSTM2)[3] <- "GSTM2"
means <- data.frame("Diagnosis"= c("Healthy\nControl", "Wilson\nDisease"),
                    "mean"=as.numeric(human_meth_stats[human_meth_stats$Name == "GSTM2" & human_meth_stats$pvalue < 0.05, c("mean_HC", "mean_WD")])) 
means$Diagnosis <- factor(means$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)

g <- ggplot()
g + 
        geom_dotplot(data=GSTM2, aes(x=Diagnosis, y=GSTM2, fill=Diagnosis, color=Diagnosis),
                     binaxis = "y", stackdir = "center", binwidth = 0.02, stackratio = 1.5, dotsize = 1) +
        geom_errorbar(data=means, aes(ymin=mean, ymax=mean, x=Diagnosis), width=0.8, color="black", size=2) +
        theme_bw(base_size = 28) +
        theme(legend.direction = 'vertical', legend.position = "none", legend.key.size = unit(1.5, "line"), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=26), legend.text=element_text(size=24),
              axis.text = element_text(color = "black", size=26), legend.background = element_blank(), axis.title.y = element_text(size=28),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4), labels = function(x) x*100) +
        coord_cartesian(ylim=c(0.3,0.8)) +
        scale_fill_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                          values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        scale_color_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                           values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        ylab("Promoter Methylation (%)") +
        annotate("text", x=0.45, y=0.8, label="GSTM2", fontface="italic", size=10, hjust=0) +
        annotate("text", x=1.5, y=max(means$mean)+0.07, label="*", size=16, fontface="bold", hjust=0.5) +
        annotate("segment", x=1.2, xend=1.8, y=max(means$mean)+0.05, yend=max(means$mean)+0.05, size=1.5)
ggsave("Figures/GSTM2 Human Promoter Methylation Dotplot no DC.png", dpi = 600, width = 6, height = 6, units = "in")

# GSTM5
GSTM5 <- cbind(pheno, as.numeric(meth[human_meth_stats$Name == "GSTM5" & human_meth_stats$pvalue < 0.05,][1,])) # Take first TSS out of 3
GSTM5$Diagnosis <- as.character(GSTM5$Diagnosis)
GSTM5$Diagnosis[GSTM5$Diagnosis == "Healthy_control"] <- "Healthy\nControl"
GSTM5$Diagnosis[GSTM5$Diagnosis == "Wilson_disease"] <- "Wilson\nDisease"
GSTM5$Diagnosis <- factor(GSTM5$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)
colnames(GSTM5)[3] <- "GSTM5"
means <- data.frame("Diagnosis"= c("Healthy\nControl", "Wilson\nDisease"),
                    "mean"=as.numeric(human_meth_stats[human_meth_stats$Name == "GSTM5" & human_meth_stats$pvalue < 0.05, c("mean_HC", "mean_WD")][1,])) 
means$Diagnosis <- factor(means$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)

g <- ggplot()
g + 
        geom_dotplot(data=GSTM5, aes(x=Diagnosis, y=GSTM5, fill=Diagnosis, color=Diagnosis),
                     binaxis = "y", stackdir = "center", binwidth = 0.02, stackratio = 1.5, dotsize = 1) +
        geom_errorbar(data=means, aes(ymin=mean, ymax=mean, x=Diagnosis), width=0.8, color="black", size=2) +
        theme_bw(base_size = 28) +
        theme(legend.direction = 'vertical', legend.position = "none", legend.key.size = unit(1.5, "line"), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=26), legend.text=element_text(size=24),
              axis.text = element_text(color = "black", size=26), legend.background = element_blank(), axis.title.y = element_text(size=28),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4), labels = function(x) x*100) +
        coord_cartesian(ylim=c(0.45,0.95)) +
        scale_fill_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                          values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        scale_color_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                           values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        ylab("Promoter Methylation (%)") +
        annotate("text", x=0.45, y=0.95, label="GSTM5", fontface="italic", size=10, hjust=0) +
        annotate("text", x=1.5, y=max(means$mean)+0.07, label="*", size=16, fontface="bold", hjust=0.5) +
        annotate("segment", x=1.2, xend=1.8, y=max(means$mean)+0.05, yend=max(means$mean)+0.05, size=1.5)
ggsave("Figures/GSTM5 Human Promoter Methylation Dotplot no DC.png", dpi = 600, width = 6, height = 6, units = "in")

# GSTP1
GSTP1 <- cbind(pheno, as.numeric(meth[human_meth_stats$Name == "GSTP1" & human_meth_stats$pvalue < 0.05,]))
GSTP1$Diagnosis <- as.character(GSTP1$Diagnosis)
GSTP1$Diagnosis[GSTP1$Diagnosis == "Healthy_control"] <- "Healthy\nControl"
GSTP1$Diagnosis[GSTP1$Diagnosis == "Wilson_disease"] <- "Wilson\nDisease"
GSTP1$Diagnosis <- factor(GSTP1$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)
colnames(GSTP1)[3] <- "GSTP1"
means <- data.frame("Diagnosis"= c("Healthy\nControl", "Wilson\nDisease"),
                    "mean"=as.numeric(human_meth_stats[human_meth_stats$Name == "GSTP1" & human_meth_stats$pvalue < 0.05, c("mean_HC", "mean_WD")])) 
means$Diagnosis <- factor(means$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)

g <- ggplot()
g + 
        geom_dotplot(data=GSTP1, aes(x=Diagnosis, y=GSTP1, fill=Diagnosis, color=Diagnosis),
                     binaxis = "y", stackdir = "center", binwidth = 0.02, stackratio = 1.5, dotsize = 1) +
        geom_errorbar(data=means, aes(ymin=mean, ymax=mean, x=Diagnosis), width=0.8, color="black", size=2) +
        theme_bw(base_size = 28) +
        theme(legend.direction = 'vertical', legend.position = "none", legend.key.size = unit(1.5, "line"), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=26), legend.text=element_text(size=24),
              axis.text = element_text(color = "black", size=26), legend.background = element_blank(), axis.title.y = element_text(size=28),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4), labels = function(x) x*100) +
        coord_cartesian(ylim=c(0.1,0.6)) +
        scale_fill_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                          values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        scale_color_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                           values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        ylab("Promoter Methylation (%)") +
        annotate("text", x=0.45, y=0.6, label="GSTP1", fontface="italic", size=10, hjust=0) +
        annotate("text", x=1.5, y=max(means$mean)+0.07, label="*", size=16, fontface="bold", hjust=0.5) +
        annotate("segment", x=1.2, xend=1.8, y=max(means$mean)+0.05, yend=max(means$mean)+0.05, size=1.5)
ggsave("Figures/GSTP1 Human Promoter Methylation Dotplot no DC.png", dpi = 600, width = 6, height = 6, units = "in")

# GSTZ1
GSTZ1 <- cbind(pheno, as.numeric(meth[human_meth_stats$Name == "GSTZ1" & human_meth_stats$pvalue < 0.05,][1,])) # Take 1st TSS out of 2
GSTZ1$Diagnosis <- as.character(GSTZ1$Diagnosis)
GSTZ1$Diagnosis[GSTZ1$Diagnosis == "Healthy_control"] <- "Healthy Control"
GSTZ1$Diagnosis[GSTZ1$Diagnosis == "Wilson_disease"] <- "Wilson Disease"
GSTZ1$Diagnosis <- factor(GSTZ1$Diagnosis, levels=c("Healthy Control", "Wilson Disease"), ordered=TRUE)
colnames(GSTZ1)[3] <- "GSTZ1"
means <- data.frame("Diagnosis"= c("Healthy Control", "Wilson Disease"),
                    "mean"=as.numeric(human_meth_stats[human_meth_stats$Name == "GSTZ1" & human_meth_stats$pvalue < 0.05, c("mean_HC", "mean_WD")][1,])) 
means$Diagnosis <- factor(means$Diagnosis, levels=c("Healthy Control", "Wilson Disease"), ordered=TRUE)

g <- ggplot()
g + 
        geom_dotplot(data=GSTZ1, aes(x=Diagnosis, y=GSTZ1, fill=Diagnosis, color=Diagnosis),
                     binaxis = "y", stackdir = "center", binwidth = 0.02, stackratio = 1.5, dotsize = 1) +
        geom_errorbar(data=means, aes(ymin=mean, ymax=mean, x=Diagnosis), width=0.8, color="black", size=2) +
        theme_bw(base_size = 28) +
        theme(legend.direction = 'vertical', legend.position = c(0.65, 0.92), legend.key.size = unit(2, "line"), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_blank(), legend.text=element_text(size=26),
              axis.text = element_text(color = "black", size=26), legend.background = element_blank(), axis.title.y = element_text(size=28),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4), labels = function(x) x*100) +
        coord_cartesian(ylim=c(0,0.5)) +
        scale_fill_manual(breaks = c("Healthy Control", "Wilson Disease"), 
                          values = c("Healthy Control"="#3366CC", "Wilson Disease"="#FF3366")) +
        scale_color_manual(breaks = c("Healthy Control", "Wilson Disease"), 
                           values = c("Healthy Control"="#3366CC", "Wilson Disease"="#FF3366")) +
        ylab("Promoter Methylation (%)") +
        annotate("text", x=0.45, y=0.5, label="GSTZ1", fontface="italic", size=10, hjust=0) +
        annotate("text", x=1.5, y=max(means$mean)+0.07, label="*", size=16, fontface="bold", hjust=0.5) +
        annotate("segment", x=1.2, xend=1.8, y=max(means$mean)+0.05, yend=max(means$mean)+0.05, size=1.5)
ggsave("Figures/GSTZ1 Human Promoter Methylation Dotplot no DC.png", dpi = 600, width = 6, height = 6, units = "in")

# MGST1
MGST1 <- cbind(pheno, as.numeric(meth[human_meth_stats$Name == "MGST1" & human_meth_stats$pvalue < 0.05,])) 
MGST1$Diagnosis <- as.character(MGST1$Diagnosis)
MGST1$Diagnosis[MGST1$Diagnosis == "Healthy_control"] <- "Healthy\nControl"
MGST1$Diagnosis[MGST1$Diagnosis == "Wilson_disease"] <- "Wilson\nDisease"
MGST1$Diagnosis <- factor(MGST1$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)
colnames(MGST1)[3] <- "MGST1"
means <- data.frame("Diagnosis"= c("Healthy\nControl", "Wilson\nDisease"),
                    "mean"=as.numeric(human_meth_stats[human_meth_stats$Name == "MGST1" & human_meth_stats$pvalue < 0.05, c("mean_HC", "mean_WD")])) 
means$Diagnosis <- factor(means$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)

g <- ggplot()
g + 
        geom_dotplot(data=MGST1, aes(x=Diagnosis, y=MGST1, fill=Diagnosis, color=Diagnosis),
                     binaxis = "y", stackdir = "center", binwidth = 0.02, stackratio = 1.5, dotsize = 1) +
        geom_errorbar(data=means, aes(ymin=mean, ymax=mean, x=Diagnosis), width=0.8, color="black", size=2) +
        theme_bw(base_size = 28) +
        theme(legend.direction = 'vertical', legend.position = "none", legend.key.size = unit(1.5, "line"), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=26), legend.text=element_text(size=24),
              axis.text = element_text(color = "black", size=26), legend.background = element_blank(), axis.title.y = element_text(size=28),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4), labels = function(x) x*100) +
        coord_cartesian(ylim=c(0.4,0.9)) +
        scale_fill_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                          values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        scale_color_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                           values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        ylab("Promoter Methylation (%)") +
        annotate("text", x=0.45, y=0.9, label="MGST1", fontface="italic", size=10, hjust=0) +
        annotate("text", x=1.5, y=max(means$mean)+0.07, label="*", size=16, fontface="bold", hjust=0.5) +
        annotate("segment", x=1.2, xend=1.8, y=max(means$mean)+0.05, yend=max(means$mean)+0.05, size=1.5)
ggsave("Figures/MGST1 Human Promoter Methylation Dotplot no DC.png", dpi = 600, width = 6, height = 6, units = "in")

# MGST3
MGST3 <- cbind(pheno, as.numeric(meth[human_meth_stats$Name == "MGST3" & human_meth_stats$pvalue < 0.05,])) 
MGST3$Diagnosis <- as.character(MGST3$Diagnosis)
MGST3$Diagnosis[MGST3$Diagnosis == "Healthy_control"] <- "Healthy\nControl"
MGST3$Diagnosis[MGST3$Diagnosis == "Wilson_disease"] <- "Wilson\nDisease"
MGST3$Diagnosis <- factor(MGST3$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)
colnames(MGST3)[3] <- "MGST3"
means <- data.frame("Diagnosis"= c("Healthy\nControl", "Wilson\nDisease"),
                    "mean"=as.numeric(human_meth_stats[human_meth_stats$Name == "MGST3" & human_meth_stats$pvalue < 0.05, c("mean_HC", "mean_WD")])) 
means$Diagnosis <- factor(means$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)

g <- ggplot()
g + 
        geom_dotplot(data=MGST3, aes(x=Diagnosis, y=MGST3, fill=Diagnosis, color=Diagnosis),
                     binaxis = "y", stackdir = "center", binwidth = 0.02, stackratio = 1.5, dotsize = 1) +
        geom_errorbar(data=means, aes(ymin=mean, ymax=mean, x=Diagnosis), width=0.8, color="black", size=2) +
        theme_bw(base_size = 28) +
        theme(legend.direction = 'vertical', legend.position = "none", legend.key.size = unit(1.5, "line"), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=26), legend.text=element_text(size=24),
              axis.text = element_text(color = "black", size=26), legend.background = element_blank(), axis.title.y = element_text(size=28),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4), labels = function(x) x*100) +
        coord_cartesian(ylim=c(0,0.5)) +
        scale_fill_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                          values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        scale_color_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                           values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        ylab("Promoter Methylation (%)") +
        annotate("text", x=0.45, y=0.5, label="MGST3", fontface="italic", size=10, hjust=0) +
        annotate("text", x=1.5, y=max(means$mean)+0.07, label="*", size=16, fontface="bold", hjust=0.5) +
        annotate("segment", x=1.2, xend=1.8, y=max(means$mean)+0.05, yend=max(means$mean)+0.05, size=1.5)
ggsave("Figures/MGST3 Human Promoter Methylation Dotplot no DC.png", dpi = 600, width = 6, height = 6, units = "in")

# PRDX1
PRDX1 <- cbind(pheno, as.numeric(meth[human_meth_stats$Name == "PRDX1" & human_meth_stats$pvalue < 0.05,])) 
PRDX1$Diagnosis <- as.character(PRDX1$Diagnosis)
PRDX1$Diagnosis[PRDX1$Diagnosis == "Healthy_control"] <- "Healthy\nControl"
PRDX1$Diagnosis[PRDX1$Diagnosis == "Wilson_disease"] <- "Wilson\nDisease"
PRDX1$Diagnosis <- factor(PRDX1$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)
colnames(PRDX1)[3] <- "PRDX1"
means <- data.frame("Diagnosis"= c("Healthy\nControl", "Wilson\nDisease"),
                    "mean"=as.numeric(human_meth_stats[human_meth_stats$Name == "PRDX1" & human_meth_stats$pvalue < 0.05, c("mean_HC", "mean_WD")])) 
means$Diagnosis <- factor(means$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)

g <- ggplot()
g + 
        geom_dotplot(data=PRDX1, aes(x=Diagnosis, y=PRDX1, fill=Diagnosis, color=Diagnosis),
                     binaxis = "y", stackdir = "center", binwidth = 0.02, stackratio = 1.5, dotsize = 1) +
        geom_errorbar(data=means, aes(ymin=mean, ymax=mean, x=Diagnosis), width=0.8, color="black", size=2) +
        theme_bw(base_size = 28) +
        theme(legend.direction = 'vertical', legend.position = "none", legend.key.size = unit(1.5, "line"), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=26), legend.text=element_text(size=24),
              axis.text = element_text(color = "black", size=26), legend.background = element_blank(), axis.title.y = element_text(size=28),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4), labels = function(x) x*100) +
        coord_cartesian(ylim=c(0.3,0.8)) +
        scale_fill_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                          values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        scale_color_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                           values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        ylab("Promoter Methylation (%)") +
        annotate("text", x=0.45, y=0.8, label="PRDX1", fontface="italic", size=10, hjust=0) +
        annotate("text", x=1.5, y=max(means$mean)+0.07, label="*", size=16, fontface="bold", hjust=0.5) +
        annotate("segment", x=1.2, xend=1.8, y=max(means$mean)+0.05, yend=max(means$mean)+0.05, size=1.5)
ggsave("Figures/PRDX1 Human Promoter Methylation Dotplot no DC.png", dpi = 600, width = 6, height = 6, units = "in")

# TXN2
TXN2 <- cbind(pheno, as.numeric(meth[human_meth_stats$Name == "TXN2" & human_meth_stats$pvalue < 0.05,][1,])) # Take 1st TSS out of 3 
TXN2$Diagnosis <- as.character(TXN2$Diagnosis)
TXN2$Diagnosis[TXN2$Diagnosis == "Healthy_control"] <- "Healthy\nControl"
TXN2$Diagnosis[TXN2$Diagnosis == "Wilson_disease"] <- "Wilson\nDisease"
TXN2$Diagnosis <- factor(TXN2$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)
colnames(TXN2)[3] <- "TXN2"
means <- data.frame("Diagnosis"= c("Healthy\nControl", "Wilson\nDisease"),
                    "mean"=as.numeric(human_meth_stats[human_meth_stats$Name == "TXN2" & human_meth_stats$pvalue < 0.05, c("mean_HC", "mean_WD")][1,])) 
means$Diagnosis <- factor(means$Diagnosis, levels=c("Healthy\nControl", "Wilson\nDisease"), ordered=TRUE)

g <- ggplot()
g + 
        geom_dotplot(data=TXN2, aes(x=Diagnosis, y=TXN2, fill=Diagnosis, color=Diagnosis),
                     binaxis = "y", stackdir = "center", binwidth = 0.02, stackratio = 1.5, dotsize = 1) +
        geom_errorbar(data=means, aes(ymin=mean, ymax=mean, x=Diagnosis), width=0.8, color="black", size=2) +
        theme_bw(base_size = 28) +
        theme(legend.direction = 'vertical', legend.position = "none", legend.key.size = unit(1.5, "line"), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(size=26), legend.text=element_text(size=24),
              axis.text = element_text(color = "black", size=26), legend.background = element_blank(), axis.title.y = element_text(size=28),
              axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        scale_y_continuous(breaks=pretty_breaks(n=4), labels = function(x) x*100) +
        coord_cartesian(ylim=c(0.3,0.8)) +
        scale_fill_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                          values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        scale_color_manual(breaks = c("Healthy\nControl", "Wilson\nDisease"), 
                           values = c("Healthy\nControl"="#3366CC", "Wilson\nDisease"="#FF3366")) +
        ylab("Promoter Methylation (%)") +
        annotate("text", x=0.45, y=0.8, label="TXN2", fontface="italic", size=10, hjust=0) +
        annotate("text", x=1.5, y=max(means$mean)+0.07, label="*", size=16, fontface="bold", hjust=0.5) +
        annotate("segment", x=1.2, xend=1.8, y=max(means$mean)+0.05, yend=max(means$mean)+0.05, size=1.5)
ggsave("Figures/TXN2 Human Promoter Methylation Dotplot no DC.png", dpi = 600, width = 6, height = 6, units = "in")

# Methylation Heatmap ####
meth_heat <- human_permeth[,4:20] # Removed DC samples
colnames(meth_heat)[2:ncol(meth_heat)] <- c(paste("HC", 1:6, sep=""), paste("WD", 1:10, sep=""))
meth_heat$Name <- as.factor(meth_heat$Name)
meth_heat$pvalue <- human_meth_stats$pvalue
meth_heat <- meth_heat %>% group_by(Name) %>% filter(pvalue == min(pvalue))
meth_heat <- unique(meth_heat)
meth_heat <- meth_heat[,1:ncol(meth_heat)-1]

meth_means <- rowMeans(meth_heat[,2:ncol(meth_heat)])
meth_heat[,2:ncol(meth_heat)] <- (meth_heat[,2:ncol(meth_heat)] - meth_means)*100
row_order <- hclust(dist(as.matrix(meth_heat[,2:ncol(meth_heat)])), "ward.D")$order # Cluster DMRs by methylation
meth_heat <- melt(meth_heat, id.vars="Name")
meth_heat$variable <- factor(meth_heat$variable, levels=c(paste("HC", 1:6, sep=""), paste("WD", 1:10, sep="")), ordered=TRUE)
meth_heat$Name <- factor(meth_heat$Name, levels=unique(meth_heat$Name)[row_order], ordered=TRUE)

gg <- ggplot(data = meth_heat)
gg +
        geom_tile(aes(x = variable, y = Name, fill = value)) +
        scale_fill_gradientn("Promoter\nMethylation (%)\n", colors = c("#0000FF", "black", "#FF0000"), values = c(0,0.55,1), na.value = "black") +
        theme_bw(base_size = 24) +
        theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", size = 1.25),
              axis.ticks = element_line(color="black", size=1), legend.key = element_blank(), 
              panel.grid.minor = element_blank(), legend.position = c(1.25, 0.84), legend.key.size = unit(1.5, "lines"),
              legend.background = element_blank(),
              plot.margin = unit(c(1,10,1,1), "lines"), 
              axis.text.y = element_text(color="black", size=12, hjust=1, face="italic"),
              axis.text.x=element_text(color="black", size=12, angle=90, hjust=1, vjust=0.5),
              axis.title=element_blank(), legend.title = element_text(size = 17)) +
        scale_x_discrete(expand=c(0,0))
ggsave("Figures/WD Human Thioredoxin Pathway Promoter Methylation Heatmap No DC.png", dpi = 600, width = 7, height = 8, units = "in")
