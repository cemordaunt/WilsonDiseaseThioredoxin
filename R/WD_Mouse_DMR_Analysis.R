# WD Mouse DMR Analysis
# Charles Mordaunt
# 4/10/18

setwd()

# Packages ####
library(GenomicRanges)
library(VennDiagram)
library(ggdendro)
library(ggplot2)
library(reshape2)
library(grid)
library(scales)
library(ggbiplot)
library(R.utils)
library(rtracklayer)
library(LOLA)
library(reshape)
library(devtools)
library(plyr)
library(sm)
library(ChIPpeakAnno)
library(biomaRt)

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

covCutoff <- function(numCtrl, numExp, minDiff){
# Determine minimum coverage threshold
        precisionCtrl <- (1/numCtrl)*2
        precisionExp <- (1/numExp)*2
        precisionSum <- precisionCtrl + precisionExp
        coverage <- ceiling(precisionSum/minDiff)
        coverage
}

DMRttest <- function(meth_info, numCtrl, numExp, reads){
# Run ttest on DMR whole region methylation
        meth <- meth_info[,16:ncol(meth_info)]
        DMRs <- 1:nrow(meth)
        ttest <- matrix(nrow = length(DMRs), ncol = 5)
        for(i in DMRs){
                association <- NULL
                association <- t.test(meth[i,(numCtrl+1):(numCtrl+numExp)], meth[i,1:numCtrl])
                ttest[i,1] <- association$estimate[1] - association$estimate[2]
                ttest[i,2] <- association$conf.int[1]
                ttest[i,3] <- association$conf.int[2]
                ttest[i,4] <- association$statistic
                ttest[i,5] <- association$p.value
        }
        colnames(ttest) <- c("meanDiff", "confIntL", "confIntR", "tstat", "pValue")
        ttest <- as.data.frame(ttest)
        ttest$pValueFDR <- p.adjust(ttest$pValue, "fdr")
        ttest$pValueBonf <- p.adjust(ttest$pValue, "bonf")
        
        minreads <- apply(reads, 1, min)
        ttest$minReads <- minreads
        meanreads <- apply(reads, 1, mean)
        ttest$meanReads <- meanreads
        
        ttest <- cbind(meth_info[,1:15], ttest[,2:ncol(ttest)])
        ttest
}

makeGRange <- function(ttest, direction=c("all", "hyper", "hypo")){
# Convert ttest data.frame to GRange
        if("all" == direction){
                GR <- GRanges(seqnames = ttest$chr, ranges=IRanges(start=ttest$start, end=ttest$end))
        }
        if("hyper" == direction){
                ttest <- subset(ttest, meanDiff > 0)
                GR <- GRanges(seqnames = ttest$chr, ranges=IRanges(start=ttest$start, end=ttest$end))
        }
        if("hypo" == direction){
                ttest <- subset(ttest, meanDiff < 0)
                GR <- GRanges(seqnames = ttest$chr, ranges=IRanges(start=ttest$start, end=ttest$end))
        }
        GR
}

GRangetoBED <- function(GR, writeFile=TRUE, fileName){
# Convert GRange to BED and write
        GR <- as.data.frame(GR)
        BED <- GR[,c("seqnames", "start", "end")]
        BED <- BED[order(BED$seqnames, BED$start),]
        if(writeFile){
                write.table(BED, fileName , sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        }
        GR
}

ttestToBED <- function(ttest, direction, prefix){
# Convert ttest data.frame to BED
        if("all" %in% direction){
                BED <- ttest[,c("chr", "start", "end", "DMRid")]
                write.table(BED, paste(prefix, "DMRs.bed", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        }
        if("hyper" %in% direction){
                BED <- subset(ttest, meanDiff > 0)
                BED <- BED[,c("chr", "start", "end", "DMRid")]
                write.table(BED, paste(prefix, "hypermethylated_DMRs.bed", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        }
        if("hypo" %in% direction){
                BED <- subset(ttest, meanDiff < 0)
                BED <- BED[,c("chr", "start", "end", "DMRid")]
                write.table(BED, paste(prefix, "hypomethylated_DMRs.bed", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)        
        }
}

prepGREAT_mm10 <- function(ttest, direction, background, prefix, chroms=c(paste("chr", 1:19, sep=""), "chrM")){
# Make DMR and background BED for GREAT input
        # Background
        GR_background <- GRanges(seqnames = background$chr, ranges=IRanges(start=background$start, end=background$end))
        seqlevelsStyle(GR_background) <- "UCSC"
        if(!isDisjoint(GR_background)){GR_background <- disjoin(GR_background)} 
        if(length(GR_background) > 1000000){cat("\nWarning: Need to reduce background to < 1M regions\n")}
        background <- as.data.frame(GR_background)[,c("seqnames", "start", "end")]
        colnames(background) <- c("chr", "start", "end")
        background$chr <- as.character(background$chr)
        background <- background[order(background$chr, background$start),]
        background <- unique(subset(background, chr %in% chroms))
        write.table(background, paste(prefix, "Background_for_GREAT.bed", sep="_"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
        
        # All DMRs
        if("all" %in% direction){
                GR_DMRs <- GRanges(seqnames = ttest$chr, ranges=IRanges(start=ttest$start, end=ttest$end), DMRid=ttest$DMRid)
                seqlevelsStyle(GR_DMRs) <- "UCSC"
                if(!isDisjoint(GR_DMRs)){GR_DMRs <- disjoin(GR_DMRs)} 
                GR_DMRs2 <- redefineUserSets(GR_DMRs, GR_background) # metadata lost
                GR_DMRs2 <- unlist(GRangesList(GR_DMRs2))
                overlaps <- findOverlaps(GR_DMRs, GR_DMRs2, type="any", select="first") # Match regions
                elementMetadata(GR_DMRs2) <- elementMetadata(GR_DMRs)[[1]][overlaps] # Add metadata back in
                DMRs <- as.data.frame(GR_DMRs2)[,c("seqnames", "start", "end")]#, "X")]
                colnames(DMRs) <- c("chr", "start", "end")#, "DMRid")
                DMRs$chr <- as.character(DMRs$chr)
                DMRs <- DMRs[order(DMRs$chr, DMRs$start),]
                DMRs <- unique(subset(DMRs, chr %in% chroms))
                write.table(DMRs, paste(prefix, "DMRs_for_GREAT.bed", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        }
        
        # Hyper DMRs
        if("hyper" %in% direction){
                hyper <- subset(ttest, meanDiff > 0)
                GR_hyper <- GRanges(seqnames = hyper$chr, ranges=IRanges(start=hyper$start, end=hyper$end), metadata=hyper$DMRid)
                seqlevelsStyle(GR_hyper) <- "UCSC"
                if(!isDisjoint(GR_hyper)){GR_hyper <- disjoin(GR_hyper)} 
                GR_hyper2 <- redefineUserSets(GR_hyper, GR_background) # metadata lost
                GR_hyper2 <- unlist(GRangesList(GR_hyper2))
                overlaps <- findOverlaps(GR_hyper, GR_hyper2, type="any", select="first") # Match regions
                elementMetadata(GR_hyper2) <- elementMetadata(GR_hyper)[[1]][overlaps] # Add metadata back in
                hyper <- as.data.frame(GR_hyper2)[,c("seqnames", "start", "end", "X")]
                colnames(hyper) <- c("chr", "start", "end", "DMRid")
                hyper$chr <- as.character(hyper$chr)
                hyper <- hyper[order(hyper$chr, hyper$start),]
                hyper <- unique(subset(hyper, chr %in% chroms))
                write.table(hyper, paste(prefix, "hypermethylated_DMRs_for_GREAT.bed", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        }
        
        # Hypo DMRs
        if("hypo" %in% direction){
                hypo <- subset(ttest, meanDiff < 0)
                GR_hypo <- GRanges(seqnames = hypo$chr, ranges=IRanges(start=hypo$start, end=hypo$end), metadata=hypo$DMRid)
                seqlevelsStyle(GR_hypo) <- "UCSC"
                if(!isDisjoint(GR_hypo)){GR_hypo <- disjoin(GR_hypo)} 
                GR_hypo2 <- redefineUserSets(GR_hypo, GR_background) # metadata lost
                GR_hypo2 <- unlist(GRangesList(GR_hypo2))
                overlaps <- findOverlaps(GR_hypo, GR_hypo2, type="any", select="first") # Match regions
                elementMetadata(GR_hypo2) <- elementMetadata(GR_hypo)[[1]][overlaps] # Add metadata back in
                hypo <- as.data.frame(GR_hypo2)[,c("seqnames", "start", "end", "X")]
                colnames(hypo) <- c("chr", "start", "end", "DMRid")
                hypo$chr <- as.character(hypo$chr)
                hypo <- hypo[order(hypo$chr, hypo$start),]
                hypo <- unique(subset(hypo, chr %in% chroms))
                write.table(hypo, paste(prefix, "hypomethylated_DMRs_for_GREAT.bed", sep="_"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        }
}

# Heatmap with Pheno Data Functions
mydplot_pheno <- function(ddata, row=!col, col=!row, labels=col) {
# Make dendrogram
        ## plot a dendrogram
        yrange <- range(ddata$segments$y)
        yd <- yrange[2] - yrange[1]
        nc <- max(nchar(as.character(ddata$labels$label)))
        tangle <- if(row) { 0 } else { 90 }
        tshow <- col
        p <- ggplot() +
                geom_segment(data=segment(ddata), aes(x=x, y=y, xend=xend, yend=yend), lwd = 0.45) +
                labs(x = NULL, y = NULL) + theme_dendro()
        if(row) {
                p <- p +
                        scale_x_continuous(expand=c(0.5/length(ddata$labels$x),0)) +
                        coord_flip()
        } else {
                p <- p +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, color = "black"))
        }
        return(p)
}

g_legend_pheno <-function(a.gplot){
# Get phenotype legend
        ## from
        ## http://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
}

ggheatmap.show_pheno <- function(L, widths=c(0.02,0.85,0.11,0.02), heights=c(0.02,0.2,0.06,0.7,0.02)){
# plot heatmap with phenotype data
        grid.newpage()
        top.layout <- grid.layout(nrow = 5, ncol = 4,
                                  widths = unit(widths, "null"),
                                  heights = unit(heights, "null"))
        pushViewport(viewport(layout=top.layout))
        print(L$col, vp=viewport(layout.pos.col=2, layout.pos.row=2))
        print(L$row, vp=viewport(layout.pos.col=3, layout.pos.row=4))
        ## print centre without legend
        print(L$centre +
                      theme(axis.line=element_blank(),
                            axis.text.x=element_blank(),axis.text.y=element_blank(),
                            axis.ticks=element_blank(),
                            axis.title.x=element_blank(),axis.title.y=element_blank(),
                            legend.position="none",
                            panel.background=element_blank(),
                            panel.border=element_blank(),panel.grid.major=element_blank(),
                            panel.grid.minor=element_blank(),plot.background=element_blank()),
              vp=viewport(layout.pos.col=2, layout.pos.row=4))
        
        print(L$phenoData +
                      theme_bw(base_size = 24) +
                      theme(panel.grid.major = element_blank(), panel.border = element_blank(),
                            legend.key = element_blank(), legend.key.size = unit(1, "lines"),
                            panel.grid.minor = element_blank(), legend.position = "none", 
                            legend.background = element_blank(), legend.text = element_text(size = 12, color = "Black"),
                            plot.margin = unit(c(0,0,0,-0.45), "lines"), 
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.title = element_blank(), 
                            legend.title = element_blank(),
                            plot.title = element_blank()), vp=viewport(layout.pos.col=2, layout.pos.row=3))
        
        ## add heatmap legend
        legend <- g_legend_pheno(L$centre +
                                         theme(legend.title = element_blank(), 
                                               legend.text = element_text(size = 15),
                                               legend.background = element_blank(),
                                               legend.position = c(0.7, -0.7)))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=2))
        grid.draw(legend)
        upViewport(0)
        
        ## add pheno legend
        phenoLegend <- g_legend_pheno(L$phenoData +
                                              theme(legend.title = element_blank(), 
                                                    legend.text = element_text(size = 15),
                                                    legend.direction = "vertical",
                                                    legend.position = c(0.94, 0.93),
                                                    legend.background = element_blank()))
        pushViewport(viewport(layout.pos.col=3, layout.pos.row=3))
        grid.draw(phenoLegend)
        upViewport(0)
}

ggheatmap2_pheno <- function(x, phenoData, hm.colours=my.colours, my.values, low, high) {
# Make heatmap object with phenotype data
        if(is.null(colnames(x)))
                colnames(x) <- sprintf("col%s",1:ncol(x))
        if(is.null(rownames(x)))
                rownames(x) <- sprintf("row%s",1:nrow(x))
        ## plot a heatmap
        ## x is an expression matrix
        row.hc <- hclust(dist(x), "ward.D")
        col.hc <- hclust(dist(t(x)), "ward.D")
        row.dendro <- dendro_data(as.dendrogram(row.hc),type="rectangle")
        col.dendro <- dendro_data(as.dendrogram(col.hc),type="rectangle")
        
        ## dendro plots
        col.plot <- mydplot_pheno(col.dendro, col=TRUE, labels=FALSE) +
                theme(plot.margin = unit(c(0,-0.5,0,-0.6), "lines"),
                      axis.text.x = element_blank())
        row.plot <- mydplot_pheno(row.dendro, row=TRUE, labels=FALSE) +
                theme(plot.margin = unit(c(0,2,0,0), "lines"))
        
        ## order of the dendros
        col.ord <- match(col.dendro$labels$label, colnames(x))
        row.ord <- match(row.dendro$labels$label, rownames(x))
        xx <- x[row.ord,col.ord]
        dimnames(xx) <- NULL
        xx <- melt(xx)
        
        # Heatmap
        centre.plot <- ggplot(xx, aes(X2,X1)) + 
                geom_tile(aes(fill=value, colour=value)) +
                scale_fill_gradientn(colours = hm.colours, values = my.values, limits = c(low, high), na.value = "black") +
                scale_colour_gradientn(colours = hm.colours, values = my.values, limits = c(low, high), na.value = "black") +
                labs(x = NULL, y = NULL) +
                scale_x_continuous(expand=c(0,0)) +
                scale_y_continuous(expand=c(0,0),breaks = NULL) +
                theme(plot.margin = unit(rep(0, 4), "lines"))
        
        # phenoData
        sample.ord <- match(col.dendro$labels$label, as.character(phenoData$Sample))
        phenoData$Sample <- factor(as.character(phenoData$Sample), levels = as.character(phenoData$Sample)[sample.ord], ordered = TRUE)
        phenoData_m <- melt(phenoData, id.vars = "Sample")
        phenoData_m$variable <- factor(phenoData_m$variable, levels = rev(unique(phenoData_m$variable)), ordered = TRUE)
        phenoData.plot <- ggplot(phenoData_m, aes(Sample, variable)) +
                geom_tile(aes(fill=value, color=value)) +
                scale_x_discrete(expand=c(0,0)) +
                scale_y_discrete(expand=c(0,0)) +
                scale_color_manual(breaks = c("WTctrl", "txJctrl", "WTchol", "txJchol"), 
                                   values = c("WTctrl"="#3366CC", "txJctrl"="#FF3366", "WTchol"="#009933","txJchol"="#FF66CC")) +
                scale_fill_manual(breaks = c("WTctrl", "txJctrl", "WTchol", "txJchol"), 
                                  values = c("WTctrl"="#3366CC", "txJctrl"="#FF3366", "WTchol"="#009933","txJchol"="#FF66CC"))
        ret <- list(col=col.plot,row=row.plot,centre=centre.plot, phenoData=phenoData.plot)
        invisible(ret)
}

# Combine Files ####
chroms <- c(paste("chr",1:19,sep=""), "chrM")

# DMRs
suffix <- "_silver_DMR_info.txt"
DMRs_WTctrl_vs_WTchol <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_WTchol/WTctrl_vs_WTchol_", suffix = suffix)
DMRs_WTctrl_vs_txJctrl <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_txJctrl/WTctrl_vs_txJctrl_", suffix = suffix)
DMRs_txJctrl_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_txJctrl_vs_txJchol/txJctrl_vs_txJchol_", suffix = suffix)
DMRs_WTchol_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_WTchol_vs_txJchol/WTchol_vs_txJchol_", suffix = suffix)
DMRs_WTctrl_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_txJchol/WTctrl_vs_txJchol_", suffix = suffix)
DMRs_WTchol_vs_txJctrl <- combineFiles(chroms=chroms, prefix="DMRs_WTchol_vs_txJctrl/WTchol_vs_txJctrl_", suffix = suffix)

# Methylation
suffix <- "_silver_DMR_methylation.txt"
meth_WTctrl_vs_WTchol <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_WTchol/WTctrl_vs_WTchol_", suffix = suffix)
meth_WTctrl_vs_txJctrl <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_txJctrl/WTctrl_vs_txJctrl_", suffix = suffix)
meth_txJctrl_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_txJctrl_vs_txJchol/txJctrl_vs_txJchol_", suffix = suffix)
meth_WTchol_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_WTchol_vs_txJchol/WTchol_vs_txJchol_", suffix = suffix)
meth_WTctrl_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_txJchol/WTctrl_vs_txJchol_", suffix = suffix)
meth_WTchol_vs_txJctrl <- combineFiles(chroms=chroms, prefix="DMRs_WTchol_vs_txJctrl/WTchol_vs_txJctrl_", suffix = suffix)

# Coverage
suffix <- "_silver_DMR_coverage.txt"
cov_WTctrl_vs_WTchol <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_WTchol/WTctrl_vs_WTchol_", suffix = suffix)
cov_WTctrl_vs_txJctrl <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_txJctrl/WTctrl_vs_txJctrl_", suffix = suffix)
cov_txJctrl_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_txJctrl_vs_txJchol/txJctrl_vs_txJchol_", suffix = suffix)
cov_WTchol_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_WTchol_vs_txJchol/WTchol_vs_txJchol_", suffix = suffix)
cov_WTctrl_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_txJchol/WTctrl_vs_txJchol_", suffix = suffix)
cov_WTchol_vs_txJctrl <- combineFiles(chroms=chroms, prefix="DMRs_WTchol_vs_txJctrl/WTchol_vs_txJctrl_", suffix = suffix)

# Background
suffix <- "_background_DMR_coverage.txt"
backCov_WTctrl_vs_WTchol <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_WTchol/WTctrl_vs_WTchol_", suffix = suffix)
backCov_WTctrl_vs_txJctrl <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_txJctrl/WTctrl_vs_txJctrl_", suffix = suffix)
backCov_txJctrl_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_txJctrl_vs_txJchol/txJctrl_vs_txJchol_", suffix = suffix)
backCov_WTchol_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_WTchol_vs_txJchol/WTchol_vs_txJchol_", suffix = suffix)
backCov_WTctrl_vs_txJchol <- combineFiles(chroms=chroms, prefix="DMRs_WTctrl_vs_txJchol/WTctrl_vs_txJchol_", suffix = suffix)
backCov_WTchol_vs_txJctrl <- combineFiles(chroms=chroms, prefix="DMRs_WTchol_vs_txJctrl/WTchol_vs_txJctrl_", suffix = suffix)

# Subset WTctrl_vs_WTchol DMRs ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_WTctrl_vs_WTchol_info <- merge(meth_WTctrl_vs_WTchol, DMRs_WTctrl_vs_WTchol, by = c("chr", "start", "end"))
meth_WTctrl_vs_WTchol_info <- meth_WTctrl_vs_WTchol_info[,c(1:3,(4+numSamples):ncol(meth_WTctrl_vs_WTchol_info),4:(3+numSamples))]
cov_WTctrl_vs_WTchol_info <- merge(cov_WTctrl_vs_WTchol, DMRs_WTctrl_vs_WTchol, by = c("chr", "start", "end"))
cov_WTctrl_vs_WTchol_info <- cov_WTctrl_vs_WTchol_info[,c(1:3,(4+numSamples):ncol(cov_WTctrl_vs_WTchol_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_WTctrl_vs_WTchol_info[,16:ncol(cov_WTctrl_vs_WTchol_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample

# T-test and subset
ttest_WTctrl_vs_WTchol <- DMRttest(meth_info=meth_WTctrl_vs_WTchol_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_WTctrl_vs_WTchol <- subset(ttest_WTctrl_vs_WTchol, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #1186 -> 120
ttest_WTctrl_vs_WTchol <- ttest_WTctrl_vs_WTchol[order(ttest_WTctrl_vs_WTchol$chr, ttest_WTctrl_vs_WTchol$start),]
ttest_WTctrl_vs_WTchol$DMRid <- paste("DMR", 1:nrow(ttest_WTctrl_vs_WTchol), sep="_")

# Subset WTctrl_vs_txJctrl DMRs ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_WTctrl_vs_txJctrl_info <- merge(meth_WTctrl_vs_txJctrl, DMRs_WTctrl_vs_txJctrl, by = c("chr", "start", "end"))
meth_WTctrl_vs_txJctrl_info <- meth_WTctrl_vs_txJctrl_info[,c(1:3,(4+numSamples):ncol(meth_WTctrl_vs_txJctrl_info),4:(3+numSamples))]
cov_WTctrl_vs_txJctrl_info <- merge(cov_WTctrl_vs_txJctrl, DMRs_WTctrl_vs_txJctrl, by = c("chr", "start", "end"))
cov_WTctrl_vs_txJctrl_info <- cov_WTctrl_vs_txJctrl_info[,c(1:3,(4+numSamples):ncol(cov_WTctrl_vs_txJctrl_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_WTctrl_vs_txJctrl_info[,16:ncol(cov_WTctrl_vs_txJctrl_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample

# T-test and subset
ttest_WTctrl_vs_txJctrl <- DMRttest(meth_info=meth_WTctrl_vs_txJctrl_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_WTctrl_vs_txJctrl <- subset(ttest_WTctrl_vs_txJctrl, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #1076 -> 97
ttest_WTctrl_vs_txJctrl <- ttest_WTctrl_vs_txJctrl[order(ttest_WTctrl_vs_txJctrl$chr, ttest_WTctrl_vs_txJctrl$start),]
ttest_WTctrl_vs_txJctrl$DMRid <- paste("DMR", 1:nrow(ttest_WTctrl_vs_txJctrl), sep="_")

# Subset txJctrl_vs_txJchol DMRs ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_txJctrl_vs_txJchol_info <- merge(meth_txJctrl_vs_txJchol, DMRs_txJctrl_vs_txJchol, by = c("chr", "start", "end"))
meth_txJctrl_vs_txJchol_info <- meth_txJctrl_vs_txJchol_info[,c(1:3,(4+numSamples):ncol(meth_txJctrl_vs_txJchol_info),4:(3+numSamples))]
cov_txJctrl_vs_txJchol_info <- merge(cov_txJctrl_vs_txJchol, DMRs_txJctrl_vs_txJchol, by = c("chr", "start", "end"))
cov_txJctrl_vs_txJchol_info <- cov_txJctrl_vs_txJchol_info[,c(1:3,(4+numSamples):ncol(cov_txJctrl_vs_txJchol_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_txJctrl_vs_txJchol_info[,16:ncol(cov_txJctrl_vs_txJchol_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample

# T-test and subset
ttest_txJctrl_vs_txJchol <- DMRttest(meth_info=meth_txJctrl_vs_txJchol_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_txJctrl_vs_txJchol <- subset(ttest_txJctrl_vs_txJchol, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #1299 -> 110
ttest_txJctrl_vs_txJchol <- ttest_txJctrl_vs_txJchol[order(ttest_txJctrl_vs_txJchol$chr, ttest_txJctrl_vs_txJchol$start),]
ttest_txJctrl_vs_txJchol$DMRid <- paste("DMR", 1:nrow(ttest_txJctrl_vs_txJchol), sep="_")

# Subset WTchol_vs_txJchol DMRs ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_WTchol_vs_txJchol_info <- merge(meth_WTchol_vs_txJchol, DMRs_WTchol_vs_txJchol, by = c("chr", "start", "end"))
meth_WTchol_vs_txJchol_info <- meth_WTchol_vs_txJchol_info[,c(1:3,(4+numSamples):ncol(meth_WTchol_vs_txJchol_info),4:(3+numSamples))]
cov_WTchol_vs_txJchol_info <- merge(cov_WTchol_vs_txJchol, DMRs_WTchol_vs_txJchol, by = c("chr", "start", "end"))
cov_WTchol_vs_txJchol_info <- cov_WTchol_vs_txJchol_info[,c(1:3,(4+numSamples):ncol(cov_WTchol_vs_txJchol_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_WTchol_vs_txJchol_info[,16:ncol(cov_WTchol_vs_txJchol_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample

# T-test and subset
ttest_WTchol_vs_txJchol <- DMRttest(meth_info=meth_WTchol_vs_txJchol_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_WTchol_vs_txJchol <- subset(ttest_WTchol_vs_txJchol, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #1458 -> 117
ttest_WTchol_vs_txJchol <- ttest_WTchol_vs_txJchol[order(ttest_WTchol_vs_txJchol$chr, ttest_WTchol_vs_txJchol$start),]
ttest_WTchol_vs_txJchol$DMRid <- paste("DMR", 1:nrow(ttest_WTchol_vs_txJchol), sep="_")

# Subset WTctrl_vs_txJchol DMRs ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_WTctrl_vs_txJchol_info <- merge(meth_WTctrl_vs_txJchol, DMRs_WTctrl_vs_txJchol, by = c("chr", "start", "end"))
meth_WTctrl_vs_txJchol_info <- meth_WTctrl_vs_txJchol_info[,c(1:3,(4+numSamples):ncol(meth_WTctrl_vs_txJchol_info),4:(3+numSamples))]
cov_WTctrl_vs_txJchol_info <- merge(cov_WTctrl_vs_txJchol, DMRs_WTctrl_vs_txJchol, by = c("chr", "start", "end"))
cov_WTctrl_vs_txJchol_info <- cov_WTctrl_vs_txJchol_info[,c(1:3,(4+numSamples):ncol(cov_WTctrl_vs_txJchol_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_WTctrl_vs_txJchol_info[,16:ncol(cov_WTctrl_vs_txJchol_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample

# T-test and subset
ttest_WTctrl_vs_txJchol <- DMRttest(meth_info=meth_WTctrl_vs_txJchol_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_WTctrl_vs_txJchol <- subset(ttest_WTctrl_vs_txJchol, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #1326 -> 121
ttest_WTctrl_vs_txJchol <- ttest_WTctrl_vs_txJchol[order(ttest_WTctrl_vs_txJchol$chr, ttest_WTctrl_vs_txJchol$start),]
ttest_WTctrl_vs_txJchol$DMRid <- paste("DMR", 1:nrow(ttest_WTctrl_vs_txJchol), sep="_")

# Subset WTchol_vs_txJctrl DMRs ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Merge Tables
meth_WTchol_vs_txJctrl_info <- merge(meth_WTchol_vs_txJctrl, DMRs_WTchol_vs_txJctrl, by = c("chr", "start", "end"))
meth_WTchol_vs_txJctrl_info <- meth_WTchol_vs_txJctrl_info[,c(1:3,(4+numSamples):ncol(meth_WTchol_vs_txJctrl_info),4:(3+numSamples))]
cov_WTchol_vs_txJctrl_info <- merge(cov_WTchol_vs_txJctrl, DMRs_WTchol_vs_txJctrl, by = c("chr", "start", "end"))
cov_WTchol_vs_txJctrl_info <- cov_WTchol_vs_txJctrl_info[,c(1:3,(4+numSamples):ncol(cov_WTchol_vs_txJctrl_info),4:(3+numSamples))]

# Coverage
reads <- as.matrix(cov_WTchol_vs_txJctrl_info[,16:ncol(cov_WTchol_vs_txJctrl_info)])
table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample

# T-test and subset
ttest_WTchol_vs_txJctrl <- DMRttest(meth_info=meth_WTchol_vs_txJctrl_info, numCtrl=numCtrl, numExp=numExp, reads=reads)
ttest_WTchol_vs_txJctrl <- subset(ttest_WTchol_vs_txJctrl, abs(meanDiff) >= minDiff & minReads >= coverage & pValue <= 0.05) #1212 -> 107
ttest_WTchol_vs_txJctrl <- ttest_WTchol_vs_txJctrl[order(ttest_WTchol_vs_txJctrl$chr, ttest_WTchol_vs_txJctrl$start),]
ttest_WTchol_vs_txJctrl$DMRid <- paste("DMR", 1:nrow(ttest_WTchol_vs_txJctrl), sep="_")

# Concatenate DMRs into one table ####
ttest_list <- list(ttest_WTctrl_vs_WTchol, ttest_WTctrl_vs_txJctrl, ttest_txJctrl_vs_txJchol, ttest_WTchol_vs_txJchol, ttest_WTctrl_vs_txJchol, ttest_WTchol_vs_txJctrl)
names(ttest_list) <- c("WTctrl_vs_WTchol", "WTctrl_vs_txJctrl", "txJctrl_vs_txJchol", "WTchol_vs_txJchol", "WTctrl_vs_txJchol", "WTchol_vs_txJctrl")

all_DMRs <- data.frame(NULL)
for(i in 1:length(ttest_list)){
        temp <- NULL
        tempHyper <- NULL
        tempHypo <- NULL
        temp <- ttest_list[[i]]
        tempHyper <- subset(temp, meanDiff > 0, select=c(chr, start, end))
        tempHyper$DMRlist <- rep(paste(names(ttest_list)[i], "hyper", sep="_"))
        tempHypo <- subset(temp, meanDiff < 0, select=c(chr, start, end))
        tempHypo$DMRlist <- rep(paste(names(ttest_list)[i], "hypo", sep="_"))
        all_DMRs <- rbind(all_DMRs, tempHyper, tempHypo)
}
write.table(all_DMRs, "Tables/Concatenated_DMRs_WD_Mouse_Liver.bed", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Make Bed files ####
ttestToBED(ttest=ttest_WTctrl_vs_WTchol, direction=c("all", "hyper", "hypo"), "UCSC Tracks/WTctrl_vs_WTchol_WD_Mouse_Liver")
ttestToBED(ttest=ttest_WTctrl_vs_txJctrl, direction=c("all", "hyper", "hypo"), "UCSC Tracks/WTctrl_vs_txJctrl_WD_Mouse_Liver")
ttestToBED(ttest=ttest_txJctrl_vs_txJchol, direction=c("all", "hyper", "hypo"), "UCSC Tracks/txJctrl_vs_txJchol_WD_Mouse_Liver")
ttestToBED(ttest=ttest_WTchol_vs_txJchol, direction=c("all", "hyper", "hypo"), "UCSC Tracks/WTchol_vs_txJchol_WD_Mouse_Liver")
ttestToBED(ttest=ttest_WTctrl_vs_txJchol, direction=c("all", "hyper", "hypo"), "UCSC Tracks/WTctrl_vs_txJchol_WD_Mouse_Liver")
ttestToBED(ttest=ttest_WTchol_vs_txJctrl, direction=c("all", "hyper", "hypo"), "UCSC Tracks/WTchol_vs_txJctrl_WD_Mouse_Liver")

# Subset WTctrl_vs_WTchol Background ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Coverage
reads <- as.matrix(backCov_WTctrl_vs_WTchol[,4:ncol(backCov_WTctrl_vs_WTchol)])
#table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample
backCov_WTctrl_vs_WTchol$minReads <- apply(reads, 1, min)

# Subset
backCov_WTctrl_vs_WTchol <- subset(backCov_WTctrl_vs_WTchol, minReads >= coverage)
background_WTctrl_vs_WTchol <- backCov_WTctrl_vs_WTchol[,c("chr", "start", "end")]
background_WTctrl_vs_WTchol <- background_WTctrl_vs_WTchol[order(background_WTctrl_vs_WTchol$chr, background_WTctrl_vs_WTchol$start),]
write.table(background_WTctrl_vs_WTchol, "UCSC Tracks/WTctrl_vs_WTchol_WD_Mouse_Liver_Background.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Subset WTctrl_vs_txJctrl Background ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Coverage
reads <- as.matrix(backCov_WTctrl_vs_txJctrl[,4:ncol(backCov_WTctrl_vs_txJctrl)])
#table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample
backCov_WTctrl_vs_txJctrl$minReads <- apply(reads, 1, min)

# Subset
backCov_WTctrl_vs_txJctrl <- subset(backCov_WTctrl_vs_txJctrl, minReads >= coverage)
background_WTctrl_vs_txJctrl <- backCov_WTctrl_vs_txJctrl[,c("chr", "start", "end")]
background_WTctrl_vs_txJctrl <- background_WTctrl_vs_txJctrl[order(background_WTctrl_vs_txJctrl$chr, background_WTctrl_vs_txJctrl$start),]
write.table(background_WTctrl_vs_txJctrl, "UCSC Tracks/WTctrl_vs_txJctrl_WD_Mouse_Liver_Background.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Subset txJctrl_vs_txJchol Background ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Coverage
reads <- as.matrix(backCov_txJctrl_vs_txJchol[,4:ncol(backCov_txJctrl_vs_txJchol)])
#table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample
backCov_txJctrl_vs_txJchol$minReads <- apply(reads, 1, min)

# Subset
backCov_txJctrl_vs_txJchol <- subset(backCov_txJctrl_vs_txJchol, minReads >= coverage)
background_txJctrl_vs_txJchol <- backCov_txJctrl_vs_txJchol[,c("chr", "start", "end")]
background_txJctrl_vs_txJchol <- background_txJctrl_vs_txJchol[order(background_txJctrl_vs_txJchol$chr, background_txJctrl_vs_txJchol$start),]
write.table(background_txJctrl_vs_txJchol, "UCSC Tracks/txJctrl_vs_txJchol_WD_Mouse_Liver_Background.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Subset WTchol_vs_txJchol Background ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Coverage
reads <- as.matrix(backCov_WTchol_vs_txJchol[,4:ncol(backCov_WTchol_vs_txJchol)])
#table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample
backCov_WTchol_vs_txJchol$minReads <- apply(reads, 1, min)

# Subset
backCov_WTchol_vs_txJchol <- subset(backCov_WTchol_vs_txJchol, minReads >= coverage)
background_WTchol_vs_txJchol <- backCov_WTchol_vs_txJchol[,c("chr", "start", "end")]
background_WTchol_vs_txJchol <- background_WTchol_vs_txJchol[order(background_WTchol_vs_txJchol$chr, background_WTchol_vs_txJchol$start),]
write.table(background_WTchol_vs_txJchol, "UCSC Tracks/WTchol_vs_txJchol_WD_Mouse_Liver_Background.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Subset WTctrl_vs_txJchol Background ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Coverage
reads <- as.matrix(backCov_WTctrl_vs_txJchol[,4:ncol(backCov_WTctrl_vs_txJchol)])
#table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample
backCov_WTctrl_vs_txJchol$minReads <- apply(reads, 1, min)

# Subset
backCov_WTctrl_vs_txJchol <- subset(backCov_WTctrl_vs_txJchol, minReads >= coverage)
background_WTctrl_vs_txJchol <- backCov_WTctrl_vs_txJchol[,c("chr", "start", "end")]
background_WTctrl_vs_txJchol <- background_WTctrl_vs_txJchol[order(background_WTctrl_vs_txJchol$chr, background_WTctrl_vs_txJchol$start),]
write.table(background_WTctrl_vs_txJchol, "UCSC Tracks/WTctrl_vs_txJchol_WD_Mouse_Liver_Background.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Subset WTchol_vs_txJctrl Background ####
# Parameters
numCtrl <- 4
numExp <- 4
numSamples <- numCtrl + numExp
minDiff <- 0.1

# Coverage
reads <- as.matrix(backCov_WTchol_vs_txJctrl[,4:ncol(backCov_WTchol_vs_txJctrl)])
#table(is.na(reads)) # No NA values
#if(sum(is.na(reads))>0){reads[is.na(reads)] <- 0} # replaced NA values with 0 (no reads)
coverage <- covCutoff(numCtrl=numCtrl, numExp=numExp, minDiff=minDiff) #10 reads per sample
backCov_WTchol_vs_txJctrl$minReads <- apply(reads, 1, min)

# Subset
backCov_WTchol_vs_txJctrl <- subset(backCov_WTchol_vs_txJctrl, minReads >= coverage)
background_WTchol_vs_txJctrl <- backCov_WTchol_vs_txJctrl[,c("chr", "start", "end")]
background_WTchol_vs_txJctrl <- background_WTchol_vs_txJctrl[order(background_WTchol_vs_txJctrl$chr, background_WTchol_vs_txJctrl$start),]
write.table(background_WTchol_vs_txJctrl, "UCSC Tracks/WTchol_vs_txJctrl_WD_Mouse_Liver_Background.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names=FALSE)

# Overlap DMRs ####
# Make DMR GRanges
GR_WTctrl_vs_WTchol <- makeGRange(ttest=ttest_WTctrl_vs_WTchol, direction="all")
GR_WTctrl_vs_txJctrl <- makeGRange(ttest=ttest_WTctrl_vs_txJctrl, direction="all")
GR_txJctrl_vs_txJchol <- makeGRange(ttest=ttest_txJctrl_vs_txJchol, direction="all")
GR_WTchol_vs_txJchol <- makeGRange(ttest=ttest_WTchol_vs_txJchol, direction="all")
GR_WTctrl_vs_txJchol <- makeGRange(ttest=ttest_WTctrl_vs_txJchol, direction="all")
GR_WTchol_vs_txJctrl <- makeGRange(ttest=ttest_WTchol_vs_txJctrl, direction="all")
GR_WTctrl_vs_txJctrl_hyper <- makeGRange(ttest=ttest_WTctrl_vs_txJctrl, direction="hyper")
GR_WTctrl_vs_txJctrl_hypo <- makeGRange(ttest=ttest_WTctrl_vs_txJctrl, direction="hypo")
GR_WTctrl_vs_txJchol_hyper <- makeGRange(ttest=ttest_WTctrl_vs_txJchol, direction="hyper")
GR_WTctrl_vs_txJchol_hypo <- makeGRange(ttest=ttest_WTctrl_vs_txJchol, direction="hypo")

# Make Background GRanges ####
GR_background_WTctrl_vs_WTchol <- GRanges(seqnames = background_WTctrl_vs_WTchol$chr, ranges=IRanges(start=background_WTctrl_vs_WTchol$start, end=background_WTctrl_vs_WTchol$end))
GR_background_WTctrl_vs_txJctrl <- GRanges(seqnames = background_WTctrl_vs_txJctrl$chr, ranges=IRanges(start=background_WTctrl_vs_txJctrl$start, end=background_WTctrl_vs_txJctrl$end))
GR_background_txJctrl_vs_txJchol <- GRanges(seqnames = background_txJctrl_vs_txJchol$chr, ranges=IRanges(start=background_txJctrl_vs_txJchol$start, end=background_txJctrl_vs_txJchol$end))
GR_background_WTchol_vs_txJchol <- GRanges(seqnames = background_WTchol_vs_txJchol$chr, ranges=IRanges(start=background_WTchol_vs_txJchol$start, end=background_WTchol_vs_txJchol$end))
GR_background_WTctrl_vs_txJchol <- GRanges(seqnames = background_WTctrl_vs_txJchol$chr, ranges=IRanges(start=background_WTctrl_vs_txJchol$start, end=background_WTctrl_vs_txJchol$end))
GR_background_WTchol_vs_txJctrl <- GRanges(seqnames = background_WTchol_vs_txJctrl$chr, ranges=IRanges(start=background_WTchol_vs_txJctrl$start, end=background_WTchol_vs_txJctrl$end))
GR_background_list <- list(GR_background_WTctrl_vs_WTchol, GR_background_WTctrl_vs_txJctrl, GR_background_txJctrl_vs_txJchol, GR_background_WTchol_vs_txJchol, GR_background_WTctrl_vs_txJchol, GR_background_WTchol_vs_txJctrl)
names(GR_background_list) <- c("WTctrl_vs_WTchol", "WTctrl_vs_txJctrl", "txJctrl_vs_txJchol", "WTchol_vs_txJchol", "WTctrl_vs_txJchol", "WTchol_vs_txJctrl")

# GREAT Prep ####
prepGREAT_mm10(ttest_WTctrl_vs_WTchol, direction=c("all", "hyper", "hypo"), background=background_WTctrl_vs_WTchol, prefix="UCSC Tracks/GREAT/WTctrl_vs_WTchol_WD_Mouse_Liver")
prepGREAT_mm10(ttest_WTctrl_vs_txJctrl, direction=c("all", "hyper", "hypo"), background=background_WTctrl_vs_txJctrl, prefix="UCSC Tracks/GREAT/WTctrl_vs_txJctrl_WD_Mouse_Liver")
prepGREAT_mm10(ttest_txJctrl_vs_txJchol, direction=c("all", "hyper", "hypo"), background=background_txJctrl_vs_txJchol, prefix="UCSC Tracks/GREAT/txJctrl_vs_txJchol_WD_Mouse_Liver")
prepGREAT_mm10(ttest_WTchol_vs_txJchol, direction=c("all", "hyper", "hypo"), background=background_WTchol_vs_txJchol, prefix="UCSC Tracks/GREAT/WTchol_vs_txJchol_WD_Mouse_Liver")
prepGREAT_mm10(ttest_WTctrl_vs_txJchol, direction=c("all", "hyper", "hypo"), background=background_WTctrl_vs_txJchol, prefix="UCSC Tracks/GREAT/WTctrl_vs_txJchol_WD_Mouse_Liver")
prepGREAT_mm10(ttest_WTchol_vs_txJctrl, direction=c("all", "hyper", "hypo"), background=background_WTchol_vs_txJctrl, prefix="UCSC Tracks/GREAT/WTchol_vs_txJctrl_WD_Mouse_Liver")

# Choline DMR Overlap Venns
pdf(file="Figures/Choline Effect Hyper Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_WTctrl_vs_txJctrl_hyper, GR_WTctrl_vs_txJchol_hyper), 
                        NameOfPeaks = c("WTctrl_vs_txJctrl", "WTctrl_vs_txJchol"), fontfamily="sans", cat.fontfamily="sans", 
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 180, margin = 0.02,
                        cat.cex = 2, cex = 3, fill = c("lightblue", "lightpink"), cat.pos = c(180,180), 
                        cat.dist = c(0.05, 0.05), ext.dist=-0.2, ext.length=0.9, ext.text=FALSE)
dev.off()

pdf(file="Figures/Choline Effect Hypo Overlap Venn.pdf", width=10, height=8, onefile = FALSE)
venn <- makeVennDiagram(Peaks = list(GR_WTctrl_vs_txJctrl_hypo, GR_WTctrl_vs_txJchol_hypo), 
                        NameOfPeaks = c("WTctrl_vs_txJctrl", "WTctrl_vs_txJchol"), fontfamily="sans", cat.fontfamily="sans", 
                        maxgap = 0, minoverlap = 1, connectedPeaks = "min", rotation.degree = 180, margin = 0.02,
                        cat.cex = 2, cex = 3, fill = c("lightblue", "lightpink"), cat.pos = c(180,180), 
                        cat.dist = c(0.05, 0.05), ext.dist=-0.2, ext.length=0.9, ext.text=FALSE)
dev.off()

# Concatenated DMR Methylation Analysis ####
# Methylation Data
chroms <- c(paste("chr",1:19,sep=""), "chrM")
conc_dmr_meth <- combineFiles(chroms=chroms, prefix="Concatenated_DMR_Methyl/Concatenated_DMR_Methyl_", suffix = "_DMR_methylation.txt")
table(is.na(conc_dmr_meth$VMNS001A))
# FALSE  TRUE 
# 642    30 
conc_dmr_meth <- subset(conc_dmr_meth, !is.na(VMNS001A))

# Heatmap with Pheno Data
phenoData <- data.frame("Sample"=colnames(conc_dmr_meth)[4:ncol(conc_dmr_meth)])
phenoData$Diagnosis <- rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4)
phenoData$Sample <- factor(phenoData$Sample, levels=unique(phenoData$Sample), ordered=TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels=c("WTctrl", "txJctrl", "WTchol", "txJchol"), ordered=TRUE)

meth <- conc_dmr_meth[,4:ncol(conc_dmr_meth)]
methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg
methdiff <- as.matrix(methdiff)
methplot <- ggheatmap2_pheno(x = methdiff, phenoData=phenoData, hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), 
                             my.values = c(0,0,0.5,1,1), low = min(methdiff), high = max(methdiff))
pdf(file="Figures/Concatenated DMRs All Samples Heatmap phenoData.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show_pheno(methplot, widths=c(0.02,0.81,0.12,0.06), heights=c(0.02,0.12,0.04,0.8,0.02))
dev.off()

# PCA Plot PC1 PC2
data <- t(as.matrix(meth))
diagnosis <- phenoData$Diagnosis
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)
# PC1 19.4%
# PC2 17.4%
# PC3 16.0%

g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, 
              var.axes = FALSE, varname.abbrev = FALSE, choices = 1:2,ellipse.prob = 0.95)
g + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.33, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), legend.text=element_text(size=18)) +
        coord_cartesian(xlim = c(-23, 23), ylim = c(-23,23)) +
        xlab("PC1 (19% of Variance)") +
        ylab("PC2 (17% of Variance)") +
        scale_color_manual(breaks = c("WTctrl", "txJctrl", "WTchol", "txJchol"), 
                           values = c("WTctrl"="#3366CC", "txJctrl"="#FF3366", "WTchol"="#009933","txJchol"="#FF66CC")) +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        geom_point(aes(color = diagnosis), size=2.5)
ggsave("Figures/Concatenated DMRs All Samples PC1 PC2 plot.png", dpi = 600, width = 8, height = 8, units = "in")

# Choline Effect DMR Heatmap ####
# WTctrl vs txJctrl only
chol_eff <- subset(all_DMRs, DMRlist %in% c("WTctrl_vs_txJctrl_hyper", "WTctrl_vs_txJctrl_hypo"))
chol_eff_meth <- merge(chol_eff, conc_dmr_meth, by=c("chr", "start", "end"))

phenoData <- data.frame("Sample"=colnames(chol_eff_meth)[5:ncol(chol_eff_meth)])
phenoData$Diagnosis <- rep(c("WTctrl", "txJctrl", "WTchol", "txJchol"), each=4)
phenoData$Sample <- factor(phenoData$Sample, levels=unique(phenoData$Sample), ordered=TRUE)
phenoData$Diagnosis <- factor(phenoData$Diagnosis, levels=c("WTctrl", "txJctrl", "WTchol", "txJchol"), ordered=TRUE)

meth <- chol_eff_meth[,5:ncol(chol_eff_meth)]
methavg <- rowMeans(meth, na.rm = TRUE)
methdiff <- meth - methavg
methdiff <- as.matrix(methdiff)
methplot <- ggheatmap2_pheno(x = methdiff, phenoData=phenoData, hm.colours = c("#0000FF", "#0000FF", "Black", "#FF0000", "#FF0000"), 
                             my.values = c(0,0,0.5,1,1), low = -0.34, high = 0.34)
pdf(file="Figures/Choline Effect Genotype only DMRs All Samples Heatmap phenoData.pdf", width=10, height=8, onefile = FALSE)
ggheatmap.show_pheno(methplot, widths=c(0.02,0.79,0.18,0.02), heights=c(0.02,0.12,0.04,0.8,0.02))
dev.off()

# WTctrl vs txJ ctrl only PCA Plot
data <- t(as.matrix(meth))
diagnosis <- phenoData$Diagnosis
data.pca <- prcomp(data, center = TRUE, scale. = TRUE) 
plot(data.pca, type = "l")
summary(data.pca)
# PC1 43.9%
# PC2 7.4%

g <- ggbiplot(data.pca, obs.scale = 1, var.scale = 1, groups = diagnosis, ellipse = TRUE, circle = FALSE, 
              var.axes = FALSE, varname.abbrev = FALSE, choices = 1:2,ellipse.prob = 0.95)
g + 
        theme_bw(base_size = 25) +
        theme(legend.direction = 'horizontal', legend.position = c(0.33, 0.96), panel.grid.major = element_blank(), 
              panel.border = element_rect(color = "black", size = 1.25), axis.ticks = element_line(size = 1.25), 
              legend.key = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(),
              axis.text = element_text(color = "black"), legend.background = element_blank(), legend.text=element_text(size=18)) +
        coord_cartesian(xlim = c(-12, 12), ylim = c(-12, 12)) +
        xlab("PC1 (44% of Variance)") +
        ylab("PC2 (7% of Variance)") +
        scale_color_manual(breaks = c("WTctrl", "txJctrl", "WTchol", "txJchol"), 
                           values = c("WTctrl"="#3366CC", "txJctrl"="#FF3366", "WTchol"="#009933","txJchol"="#FF66CC")) +
        scale_x_continuous(breaks=pretty_breaks(n=4)) +
        scale_y_continuous(breaks=pretty_breaks(n=4)) +
        geom_point(aes(color = diagnosis), size=2.5)
ggsave("Figures/WT vs txJ DMRs All Samples PC1 PC2 plot.png", dpi = 600, width = 8, height = 8, units = "in")
