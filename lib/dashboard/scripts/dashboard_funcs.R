

#' Extract Transcript Features
#'
#' Prepare Gene Model from GTF with Productivity information.
#' Update FLAIR gtf annotation with UTR and CDS 
#' features from productivity annotation for a single Gene.
#'
#' @param productivity data.table create with \link[prepProductivity]
#' @param gtf 
#' @param geneID 
#'
#' @return
#' @export
#'
#' @examples
extractTXFeatures <- function(productivity, gtf, geneID) {
    
    require(GenomicRanges)
    
    # skip NGO and NST as they are most likely not translated at all
    geneFeaturesDetailed <- productivity[gene_id == geneID #& productivity %in% c("PTC","PRO")
                                         , ]
    cds <- geneFeaturesDetailed[,.(chrom, "cdsStart" = thickStart, "cdsEnd" = thickEnd,strand,type = "cds"),
                                by = isoform_id]
    cdsGRL <- split(GenomicRanges::makeGRangesFromDataFrame(cds,keep.extra.columns = TRUE),
                    f = cds$isoform_id)
    
    transcript <- geneFeaturesDetailed[,.(chrom, "utrStart" = chromStart, "utrEnd" = chromEnd,strand),
                                       by = isoform_id]
    transcriptGRL <- split(GenomicRanges::makeGRangesFromDataFrame(transcript,keep.extra.columns = TRUE),
                           f = transcript$isoform_id)
    
    utrGRL <- mendoapply(FUN = function(x,y) {z = setdiff(x,y); z$type = "utr"; return(z) }, 
                         transcriptGRL,cdsGRL)
    
    geneFeatures <- gtf[gene_id == geneID & type != "transcript"]
    geneFeatures <- merge(geneFeatures,productivity[,c("isoform_id","productivity")],
                          by.x="transcript_id",
                          by.y="isoform_id")
    # geneFeatures <- geneFeatures[productivity %in% c("PRO","PTC")]
    
    geneFeaturesList <- split(as(geneFeatures,"GRanges"),
                              f = geneFeatures$transcript_id)
    
    txList <- GenomicRanges::GRangesList(mapply(function(exons,cds,utr) { 
        
        
        ## coding features
        cds <- intersect(exons,cds,ignore.strand=FALSE)
        if(length(cds)>0) cds$type <- "cds"
        ## utr features
        utr <- intersect(exons,utr,ignore.strand=FALSE)
        utr$type <- "utr"
        
        newTx <- sort(c(cds,utr))
        mcols(newTx) <- cbind(mcols(newTx),mcols(exons)[1,c("transcript_id","gene_id","productivity")])
        mcols(newTx)$group <- mcols(newTx)[1,"transcript_id"]
        
        return(newTx)
        
    }, exons = geneFeaturesList, cds = cdsGRL, utr = utrGRL))
    
    names(txList) <- sprintf("%s",substr(names(txList) ,0,13))
    names(txList) <- make.unique(names(txList),sep = "-")
    
    return(txList)
    
}

plotGeneModel <- function(txList,geneID, asEvents=NULL, tx2gene) {
    
    require(ggbio)
    
    geneName <- tx2gene[gene_id == geneID, gene_name][1]
    
    p_model <-    ggplot() + 
        ## first plot CDS
        geom_alignment(txList, aes(group = group,
                                   type = type, 
                                   fill = productivity) ,
                       cds.rect.h = 0.25,label = TRUE, stat = "identity",
                       names.expr = function(x) sprintf("%s",substr(x,0,13))) +
        ## then UTR
        ggtitle(label = geneName,
                subtitle = geneID) +
        scale_fill_manual(values = alpha(colour = c(c("PRO" = "darkgreen","PTC" = "orange",
                                     "NGO" = "black","NST" = "black"),
                                     c("ir" = "blue","es" = "red",
                                       "alt3" = "orange","alt3" = "purple")),
                                     alpha=c(rep(1,4),rep(0.1,4))),
                          labels = c(c(
                              "PRO" = "productive transcript",
                              "PTC" = "premature termination codon",
                              "NGO" = "no start codon",
                              "NST" = "start but no stop codon"
                          ),
                          c("ir" = "retained intron",
                            "es" = "skipped exon",
                            "alt3" = "alt 3' splice site",
                            "alt5" = "alt 5' splice site"
                            )),
                          name = "transcript productivity"
                          ) +
        scale_y_discrete(labels = names(txList),breaks = names(txList)) +
        theme_bw()
    
    if(! is.null(asEvents) ) {

        asFeatures <- na.omit(asEvents[gene_id == geneID & grepl("exclusion",feature_id)])

        p_model <- p_model +
            ggplot2::geom_rect(data = asFeatures, aes(xmin = start,
                                                      xmax = end,
                                                      ymin = -Inf,
                                                      ymax = Inf, fill = event),
                               alpha=0.1,show.legend = FALSE)

        if(nrow(asFeatures[event == "ir"])>=1) {
            p_model <- p_model +
                ggplot2::geom_rect(data = asFeatures[event == "ir"],
                                   aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = "ir"),
                                   alpha=0.1,fill="blue",show.legend = TRUE)
        }
        if(nrow(asFeatures[event == "es"])>=1) {
            p_model <- p_model +
                ggplot2::geom_rect(data = asFeatures[event == "es"],
                                   aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = "es"),
                                   alpha=0.1,fill="red",show.legend = TRUE) # +

        }
        if(nrow(asFeatures[event == "alt3"])>=1) {
            p_model <- p_model +
                ggplot2::geom_rect(data = asFeatures[event == "alt3"],
                                   aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = "alt3"),
                                   alpha=0.1,fill="orange",show.legend = TRUE) # +

        }
        if(nrow(asFeatures[event == "alt5"])>=1) {
            p_model <- p_model +
                ggplot2::geom_rect(data = asFeatures[event == "alt5"],
                                   aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = "alt5"),
                                   alpha=0.1,fill="purple",show.legend = TRUE) # +

        }

    }
    return( p_model) 
}

plotTxExpression <- function(norm.counts,geneID,gene2symbol,productivity,txOrder=NULL,
                             labels = c("N2" = "N2","mrg1" = bquote("mrg-1"^{"-/-"})), meta) {
    
    require(data.table)
    require(reshape2)
    require(ggplot2)
    
    
    txIds <- productivity[gene_id == geneID #& productivity %in% c("PTC","PRO")
                          ,transcript_id]
    txIds <- txIds[txIds %in% row.names(norm.counts)]
    
    expressionTable <- as.data.table(norm.counts[txIds,],keep.rownames = TRUE)
    expressionTable$transcript_id <- productivity[match(txIds,transcript_id),isoform_id]
    expressionTable$transcript_id_labels <- factor(
        make.unique(sprintf("%s",substr(expressionTable$transcript_id ,0,13)),sep = "-"),
        levels = txOrder)
    
    expressionLevels <- data.table::melt(expressionTable)
    setnames(expressionLevels, c("rn","variable"),c("tx_gene_id","sample"))
    expressionLevels <- cbind(expressionLevels,meta[expressionLevels$sample,c("group","batch")])
    
    # prod <- productivity[match(txIds,transcript_id), productivity]
    # textColor <- sapply(prod,switch, "PRO" = "grey","PTC" = "orange","NGO" = "black","NST" = "black")
    
    p_box <- ggplot(expressionLevels,
                    aes(x = transcript_id_labels, 
                        y = value, fill = group)) +
        geom_boxplot(position = "dodge") +
        # geom_jitter(position = position_dodge(width = 0.75)) +
        coord_flip() + 
        scale_fill_manual(labels = labels, 
                          breaks = names(labels),
                          values = setNames(c("#D55E00","#0072B2"),names(labels))) +
        labs(x = element_blank(), y = "Expression", fill = "Group") + 
        theme_bw() #+
    # theme(axis.text.y = element_text(angle = 0, hjust = 1, colour = textColor))
    # ggeasy::easy_text_color(teach = TRUE) +
    
    return(p_box)
    
}

plotDTU <- function(tx_exp_prod, highlight_genes, isGeneID=FALSE,
                    gene2symbol, title = bquote("N2 vs mrg-1"^{"-/-"})) {
    
    require(ggplot2)
    require(ggrepel)
    
    if(isGeneID) {
        highlight_genes <- gene2symbol[match(highlight_genes,gene_id), 
                                       gene_name]
    }
    
    tx_exp_prod <- tx_exp_prod[productivity %in% c("PTC","PRO"),]
    
    p_tx <- ggplot() + 
        geom_point(data = na.omit(tx_exp_prod),mapping =  aes(x=log2FoldChange, y = -log10(padj)),
                   color = "grey",show.legend = TRUE, size = 0.5) + 
        geom_point(data = na.omit(tx_exp_prod)[gene_name %in% highlight_genes],mapping = aes(x=log2FoldChange, y = -log10(padj),color = productivity),show.legend = TRUE, shape=21, size = 3) + 
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1,1), linetype = "dashed") +
        ggtitle(bquote("N2 vs mrg-1"^{"-/-"}),"Differential Transcript Usage") +
        scale_color_manual(values = c("PRO" = "darkgreen","PTC" = "orange",
                                      "NGO" = "black","NST" = "black"),
                           labels = c(
                               "PRO" = "productive transcript",
                               "PTC" = "premature termination codon",
                               "NGO" = "no start codon",
                               "NST" = "start but no stop codon"
                           )) +
        geom_label_repel(na.omit(tx_exp_prod)[gene_name %in% highlight_genes],#& altSpliceStat=="AS_DE"],
                         mapping = aes(x=log2FoldChange,
                                       y = -log10(padj),
                                       label = sprintf("%s",substr(isoform_id,0,13)),
                                       color = productivity),show.legend = FALSE,size = 3,force = 10) +
        # ylim(c(0,15)) +
        # xlim(c(-15,15)) +
        # labs(x = bquote("N2 Log2(Foldchange) mrg-1"^{"-/-"}) " Log2(Foldchange) ")
        theme_bw() + 
        theme(legend.title = element_blank())
    
    return(p_tx)
}


plotDGE <- function(gene_exp, highlight_genes, isGeneID=FALSE,gene2symbol, 
                    title = bquote("N2 vs mrg-1"^{"-/-"})) {
    
    require(ggplot2)
    require(ggrepel)
    
    if(isGeneID) {
        highlight_genes <- gene2symbol[match(highlight_genes,gene_id), 
                                       gene_name]
    }
    
    p_gene <- ggplot() + 
        geom_point(data = na.omit(gene_exp),mapping =  aes(x=log2FoldChange, y = -log10(padj)), color = "grey",show.legend = TRUE, size = 0.5) +
        geom_point(data = na.omit(gene_exp)[gene_name %in% highlight_genes],mapping = aes(x=log2FoldChange, y = -log10(padj)),color = "red",show.legend = TRUE,shape=21, size = 3) + 
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1,1), linetype = "dashed") +
        ggtitle(label = title,subtitle = "Differential Gene Expression") + 
        # scale_color_manual(values = c("PRO" = "grey","PTC" = "orange",
        #                               "NGO" = "black","NST" = "black"),
        #                    labels = c(
        #                        "PRO" = "productive transcript",
        #                        "PTC" = "premature termination codon",
        #                        "NGO" = "no start codon",
        #                        "NST" = "start but no stop codon"
        #                              )) +
        geom_label_repel(na.omit(gene_exp)[gene_name %in% highlight_genes],
                         mapping = aes(x=log2FoldChange,
                                       y = -log10(padj),
                                       label = gene_name), 
                         color = "black",show.legend = FALSE, size = 3) +
        theme_bw() + 
        theme(legend.title = element_blank())
    
    return(p_gene)
}


summaryPlot <- function(geneID,productivity,gtf,
                        asEvents, tx2gene, norm.counts,
                        gene2symbol ,meta ,
                        boxPlotLabels = c(N2 = "N2", mrg1 = bquote("mrg-1"^{"-/-"})), 
                        tx_exp_prod, 
                        isGeneID = FALSE, 
                        filterAsEvent = "ir",
                        filterProductivity = NULL,
                        gene_exp ,filename = NULL) {
    
    require(ggplot2)
    require(patchwork)
    
    if(!isGeneID) {
        geneID <- gene2symbol[gene_name == geneID, gene_id]
    } 
    
    geneName <- gene2symbol[geneID, gene_name]
    
    if(!is.null(filterProductivity) & 
       all(filterProductivity %in% c("PRO","PTC","NST","NGO"))) {
        productivity <- productivity[productivity %in% filterProductivity,]
    }
    
    txList <- extractTXFeatures(productivity = productivity,
                                gtf = gtf,
                                geneID = geneID)
    
    if(!is.null(filterAsEvent)) {
        if(filterAsEvent %in% c("ir","es","alt5","alt3")) {
            asEvents <- asEvents[event %in% filterAsEvent]
        }
    }
    
    p_model <- plotGeneModel(txList = txList,
                             asEvents = asEvents,
                             geneID = geneID,
                             tx2gene = tx2gene)
    
    ## extract the exact order of tx labels
    
    ## extract data from ggplot
    # pg <- ggplot_build(p_model)
    ## inspect data per layer
    # pg$data
    ## extract last label layer from p_model
    # pl <- layer_data(p_model, i = 6)
    txOrder <- with(layer_data(p_model, i = 6),label[order(y)])
    
    p_box <- plotTxExpression(norm.counts = norm.counts, 
                              geneID = geneID, 
                              gene2symbol = gene2symbol, 
                              productivity = productivity,
                              txOrder = txOrder,
                              meta = meta,
                              labels = boxPlotLabels
                             )
    
    p_tx <- plotDTU(tx_exp_prod = tx_exp_prod, 
                    highlight_genes = geneID, 
                    isGeneID = TRUE,
                    gene2symbol = gene2symbol, 
                    title = NULL)
    
    p_gene <- plotDGE(gene_exp = gene_exp, 
                      highlight_genes = geneID, 
                      isGeneID=TRUE, 
                      gene2symbol = gene2symbol, 
                      title = NULL)
    
    
    
    p_summary <- (
        p_gene + ggtitle(label = NULL) + coord_flip() + 
            theme(axis.title = element_text(size = 9)) + 
            p_tx + ggtitle(label = NULL)+ coord_flip() + 
            theme(axis.title = element_text(size = 9)) + 
            plot_layout(guides = "collect" ) & theme(legend.position = "right") 
    ) / (
        (p_model + theme(legend.position = "none",
                         axis.text.x = element_blank(),
                         axis.ticks = element_blank(),
                         text = element_text(size = 7)) + 
             ggtitle(label = NULL, subtitle = NULL)) +
            (p_box + 
                 scale_y_log10() +
                 theme(
                     axis.text.y = element_blank(),
                     axis.ticks.y = element_blank(),
                     axis.title = element_text(size = 9)))
    ) 
    
    if(!is.null(filename)) {
        for (f in filename) {
            ggsave(filename = filename, plot = p_summary,width = 8,height = 5)
        }
    }
    
    return(p_summary)
    
}

summarizeDiffSplice <- function(asCountsFile, asStatsFile) {
    
    
    ## resolve gene assignment of isoform_ids column
    resolveIsoforms <- function(x) {
        sp = strsplit(x,",")
        gene_id = sapply(sp,function(sp) tx2gene[gsub(pattern="_.*",replacement="",sp),gene_id])
        tab = lapply(gene_id, table)
        tab[lengths(tab) > 1] = lapply(tab[lengths(tab) > 1],
                                       function(x) {
                                           newTab = sort(x,decreasing=TRUE)
                                           if((unique(newTab) == 1) && (grepl(pattern = ":",x = names(newTab)[1]))) 
                                               newTab = c(newTab[-1],newTab[1])
                                           return(newTab)})
        uniq_genes = transpose(lapply(tab,names))
        iso_counts = transpose(tab)
        list(gene_id = uniq_genes[[1]],
             # alt_gene_id = uniq_genes[[-1]],
             n_isoforms=iso_counts[[1]]#,
             # alt_n_isoforms = iso_counts[[-1]]
        )}
    
    
    asCounts <- fread(asCountsFile)
    setkey(asCounts, feature_id)
    asCounts[, c("gene_id", "n_isoforms") := resolveIsoforms(isoform_ids)]
    asCounts[, isoform_ids := NULL]
    
    asStats <- fread(asStatsFile)
    setkey(asStats, feature_id)
    setnames(x = asStats, "gene_id", "coordinate")
    
    mDat <- merge(
        asStats,
        asCounts,
        by = c("feature_id", "coordinate"),
        suffixes = c(".norm", ".raw"),
        all = TRUE
    )
    setorder(mDat, pvalue, feature_id)
    
    fwrite(
        x = mDat,
        file = gsub(".tsv", ".genesAssigned.tsv", asStatsFile),
        sep = "\t",
        na = "NA"
    )
    
    
}


prepDiffSplice <- function(asListFiles) {
    
    asList <- lapply(setNames(asListFiles,
                              tstrsplit(basename(asListFiles), split = "\\.")[[2]]),
                     fread)
    
    as_dt <- rbindlist(asList,idcol = TRUE)
    setnames(as_dt,".id","event")
    as_dt[,log2FoldChange := log2(mrg11a_mrg1_batch1.norm/N21_N2_batch1.norm)]
    setorder(as_dt, pvalue, coordinate)
    
    extractCoords <- function(coord) {
        res <- tstrsplit(coord,":|-|_")
        chrom <- res[1]
        pos <- transpose(lapply(transpose(res[2:3]),function(x) as.integer(c(min(x),max(x)))))
        names(pos) <- c("start","end")
        strand <- transpose(lapply(transpose(res[2:3]),function(x){
            ifelse(x[1] > x[2],"+","-")
        } ))
        return(c(chrom,pos,strand))
    }
    
    as_dt[, c("chrom","start","end","strand") := extractCoords(coordinate) ]
    
    fwrite(as_dt,file = file.path(dataDir,"diffSplice_noShortReads","diffSplice_summary.tsv"),
           sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
    
    return(as_dt)
    
}

prepProductivity <- function(productivityBedFile) {
    
    require(data.table)
    
    tx_prod_full <- fread(productivityBedFile)
    names(tx_prod_full)[1:12] <- c("chrom","chromStart","chromEnd","name","score","strand","thickStart",
                                   "thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
    
    splitName <- function(name) {
        s <- strsplit(name,"_")
        l <- lengths(s)
        prod <- mapply(function(x,y) x[y],x = s, y = l)
        tx_id <- mapply(function(x,y) paste(x[-y],collapse = "_"),x = s, y = l)
        gene_id <- mapply(function(x,y) paste(x[y-1],collapse = "_"),x = s, y = l)
        iso_id <- mapply(function(x,y) paste(x[-c(y,y-1)],collapse = "_"),x = s, y = l)
        # s[lengths(s)>3] = lapply(s[lengths(s)>3], function (x) {
        #     c(paste(x[1:2],collapse="_"),
        #       x[c(1,2)])
        # }) 
        # s = transpose(s)
        return(list(
            "transcript_id" = tx_id,
            "isoform_id" = iso_id,
            "gene_id" = gene_id,
            "productivity" = prod ))
    }
    tx_prod_full[, c("transcript_id","isoform_id", "gene_id" ,"productivity" ) := splitName(name) ]
    
    return(tx_prod_full)
}

prepNormCounts <- function(diffExpFolder) {
    
    
    txCountsFile <- list.files(file.path(diffExpFolder),
                               pattern = "filtered_iso_counts_drim.tsv", 
                               recursive = TRUE,
                               full.names = TRUE)
    
    require(data.table)
    require(DESeq2)
    
    tx_counts <- fread(txCountsFile)
    tx_counts[,V1 := NULL]
    
    meta <- data.frame(tstrsplit(names(tx_counts)[-c(1,2)],"_")[1:3])
    names(meta) <-  c("sample","group","batch")
    rownames(meta) <- meta$sample
    
    rownames(tx_counts) <- tx_counts$feature_id
    
    designFormula <- "~ group"
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = tx_counts[,-c(1,2)],colData = meta,design = as.formula(designFormula))
    dds <- DESeq2::DESeq(dds)
    vst <- DESeq2::varianceStabilizingTransformation(dds)
    norm.counts <- DESeq2::counts(dds, normalized = TRUE)
    row.names(norm.counts) <- tx_counts$feature_id
    
    return(list(norm.counts = norm.counts, 
                meta = meta))
}


prepFlairGTF <- function(flairGtfFile) {
    
    data.table::as.data.table(rtracklayer::import.gff(con = flairGtfFile, format = 'gtf'))
    
}

prepGene2Symbol <- function(ensemblGtfFile) {

    gtfBaseRef <- data.table::as.data.table(rtracklayer::import.gff(con = ensemblGtfFile, format = 'gtf'))
    gene2symbol <- gtfBaseRef[type == "gene"][,c("gene_id", "gene_name")]
    setkey(gene2symbol, gene_id)
    
    return(gene2symbol)
}

prepTx2Gene <- function(gtfData) {
    
    tx2gene <- gtfData[type == "transcript",c("gene_id","transcript_id")]
    tx2gene[, tx_gene_id := sprintf("%s_%s",transcript_id,gene_id)]
    tx2gene[, gene_name := gene2symbol[tx2gene$gene_id, gene_name, on = "gene_id"]]
    setkey(tx2gene, transcript_id)
    
    return(tx2gene)
    
}

prepAsData <- function(diffSpliceFolder,gene2symbol) {
    
    asListFiles <- list.files(path = file.path(diffSpliceFolder),
                              pattern = "\\..*.genesAssigned.tsv", 
                              recursive = TRUE,
                              full.names = TRUE) 
    
    if(length(asListFiles) == 0 ) {
        asCountsFiles <- list.files(path = file.path(diffSpliceFolder),
                                    pattern = "\\..*quant.tsv",
                                    recursive = TRUE,
                                    full.names = TRUE) 
        
        asStatsFiles <- list.files(path = file.path(diffSpliceFolder),
                                   pattern = "\\..*drimseq2_results.tsv", 
                                   recursive = TRUE,
                                   full.names = TRUE) 
        
        mapply(summarizeDiffSplice,asCountsFiles, asStatsFiles)
        
        asListFiles <- list.files(path = file.path(diffSpliceFolder),
                                  pattern = "\\..*.genesAssigned.tsv", 
                                  recursive = TRUE,
                                  full.names = TRUE) 
    } 
    as_dt <- prepDiffSplice(asListFiles)
    as_dt[,gene_name := gene2symbol[as_dt$gene_id,gene_name, on ="gene_id"]]
    
    return(as_dt)
}

prepTxProductivity <- function(diffExpFolder,tx2gene,
                               tx_productivity,
                               Case,Control) {
    
    diuResultsFile <- list.files(file.path(diffExpFolder),
                                 pattern = "diu.*.drimseq2_results.tsv", 
                                 recursive = TRUE,
                                 full.names = TRUE)
    tx_exp <- fread(diuResultsFile)
    ### FIXME currently use my own colun labels, needs to be generalized
    Case_1 <- names(tx_exp_prod)[which(tstrsplit(names(tx_exp_prod),"_")[[2]] %in% Case)[1]]
    Control_1 <- names(tx_exp_prod)[which(tstrsplit(names(tx_exp_prod),"_")[[2]] %in% Control)[1]]

    tx_exp_prod[,log2FoldChange := log2(.SD[,1]/.SD[,2]) , .SDcols = c(Case_1,Control_1)]
    tx_exp[, gene_name := tx2gene[tx_exp$feature_id, gene_name, on = "tx_gene_id"]]
    setnames(tx_exp,"adj_pvalue","padj")
    tx_prod <- tx_productivity[,c("transcript_id","isoform_id", "gene_id" ,"productivity" )]
    setkey(tx_prod, "transcript_id")
    tx_exp_prod <- merge(tx_exp,tx_prod, by.x= c("feature_id","gene_id"),by.y = c("transcript_id","gene_id"))
    
    return(tx_exp_prod)
}


prepGeneExp <- function(diffExpFolder, gene2symbol) {
    
    dgeResultsFile <- list.files(file.path(diffExpFolder),
                                 pattern = "dge.*.deseq2_results.tsv", 
                                 recursive = TRUE,
                                 full.names = TRUE)
    gene_exp <- fread(dgeResultsFile)
    setnames(gene_exp,"V1","geneID")
    gene_exp[,gene_name := gene2symbol[gene_exp$geneID, gene_name ,on = "gene_id"]]
    
    return(gene_exp)
}


prepData <- function(productivityBedFile, 
                     flairGtfFile, 
                     ensemblGtfFile, 
                     diffSpliceFolder,
                     diffExpFolder,
                     dataDir,
                     readsManifest,
                     Case = "mrg1",
                     Control = "N2"
                     ) {

    require(data.table)
    
    #### prepare productivity     
    cat("#### prepare productivity\n")
    tx_productivity <- prepProductivity(productivityBedFile)
    
    #### prepare gtf
    cat("#### prepare gtf\n")
    gtfData <- prepFlairGTF(flairGtfFile)
    
    #### prepare gene2symbol
    cat("#### prepare gene2symbol\n")
    gene2symbol <- prepGene2Symbol(ensemblGtfFile)
    
    #### prepare tx2gene
    cat("#### prepare tx2gene\n")
    tx2gene <- prepTx2Gene(gtfData)
    
    #### prepare AS data
    cat("#### prepare AS data\n")
    as_dt <- prepAsData(diffSpliceFolder, gene2symbol)
    
    #### prepare norm.counts and meta
    cat("#### prepare norm.counts and meta\n")
    txCounts <-  prepNormCounts(diffExpFolder)
    
    #### prepare tx_exp and tx_prod data
    cat("#### prepare tx_exp and tx_prod data\n")
    tx_exp_prod <- prepTxProductivity(diffExpFolder,tx2gene,
                                      tx_productivity,Case,Control)
    
    #### prepare gene_exp
    cat("#### prepare gene_exp\n")
    gene_exp <- prepGeneExp(diffExpFolder, gene2symbol)
    
    
    #### save everything
    cat("#### save everything\n")
    saveRDS(list(
        productivity = tx_productivity,
        gtf = gtfData,
        asEvents = as_dt,
        tx2gene = tx2gene,
        norm.counts = txCounts$norm.counts,
        gene2symbol = gene2symbol,
        meta = txCounts$meta,
        tx_exp_prod = tx_exp_prod,
        gene_exp = gene_exp
    ),file = file.path(dataDir,"preppedData.RDS"))
    
    return(list(
        productivity = tx_productivity,
        gtf = gtfData,
        asEvents = as_dt,
        tx2gene = tx2gene,
        norm.counts = txCounts$norm.counts,
        gene2symbol = gene2symbol,
        meta = txCounts$meta,
        tx_exp_prod = tx_exp_prod,
        gene_exp = gene_exp
    ))
}
