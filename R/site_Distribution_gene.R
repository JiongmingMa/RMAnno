#' Visualization of RNA modification sites
#'
#' site_Distribution_gene visualize the RNA modification sites based on the
#' genomic region information and gene region information.
#'
#' @param gtf_file A data.frame or GRanges object containing gene annotations in GTF format.
#' @param geneName Character string specifying the gene symbol to plot. Only one gene can be visualized at a time.
#' @param seqInfo A Seqinfo object containing chromosome length information. If NULL, chromosome lengths will be automatically extracted from the GTF file.
#' @param gtf_path Character string specifying the path to a GTF file. Ignored if gtf_file is provided.
#' @param gr1 A GRanges object containing genomic positions to be plotted as data points or regions.
#' @param mergeSample Logical indicating whether to merge multiple samples when plotting. If TRUE, data from different samples will be aggregated.
#' @param rightmargin Numeric value specifying the right margin of the plot as a fraction of the total plot width. Default is 0.1 10.
#' @param ideogramheight Numeric value between 0 and 1 specifying the height of the chromosome ideogram track relative to the total plot height.
#' @param data1height Numeric value between 0 and 1 specifying the height of the main data track relative to the total plot height.
#' @param gene_r0 Numeric values between 0 and 1 defining the radial position inner and outer radius of the gene track in circular plots. gene_r0 is the inner radius, gene_r1 is the outer radius.
#' @param gene_r1 Numeric values between 0 and 1 defining the radial position inner and outer radius of the gene track in circular plots. gene_r0 is the inner radius, gene_r1 is the outer radius.
#' @param mergeTranscripts Logical indicating whether to merge transcript variants of the same gene into a single representation.
#' @param non.coding.exons.height Numeric value specifying the height of non-coding exons in the gene model representation.
#' @param non.coding.exons.col Numeric value specifying the height of non-coding exons in the gene model representation.
#' @param non.coding.exons.border.col Character string specifying the fill color and border color for non-coding exons.
#' @param coding.exons.col Character string specifying the fill color and border color for non-coding exons.
#' @param coding.exons.border.col Character string specifying the  fill color and border color for coding exons.
#' @param introns.col Character string specifying the color for intron lines.
#' @param region_col Character string specifying the fill color for genomic regions.
#' @param region_r0 Character string specifying the fill color for genomic regions.
#' @param region_r1 Numeric values between 0 and 1 defining the radial position inner and outer radius of region tracks in circular plots.
#' @param point1_col Character string specifying the color for primary data points.
#' @param point1_r0 Numeric values between 0 and 1 defining the radial position inner and outer radius for primary data points in circular plots.
#' @param point1_r1 Numeric values between 0 and 1 defining the radial position inner and outer radius for primary data points in circular plots.
#' @param point1_pch Integer or single character specifying the plotting character for primary data points see par"pch".
#' @param point1_cex Numeric value specifying the character expansion factor for primary data points.
#' @param point1_y Numeric value or vector specifying the y-coordinates for primary data points in linear plots.
#' @param point2_col Character string specifying the color for secondary data points.
#' @param point2_r0 Numeric values between 0 and 1 defining the radial position inner and outer radius for secondary data points in circular plots.
#' @param point2_r1 Numeric values between 0 and 1 defining the radial position inner and outer radius for secondary data points in circular plots.
#' @param point2_pch Integer or single character specifying the plotting character for secondary data points see par"pch".
#' @param point2_cex Numeric value specifying the character expansion factor for secondary data points.
#' @param point2_y Numeric value or vector specifying the y-coordinates for secondary data points in linear plots.
#' @param kpPoints A function or expression specifying custom point plotting using karyoploteR::kpPoints. Useful for adding additional data layers.
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomeInfoDb seqlevelsInUse
#' @importFrom karyoploteR getDefaultPlotParams
#' @importFrom karyoploteR plotKaryotype
#' @importFrom karyoploteR makeGenesDataFromTxDb
#' @importFrom karyoploteR mergeTranscripts
#' @importFrom karyoploteR kpPlotGenes
#' @importFrom karyoploteR kpAddLabels
#' @importFrom GenomicRanges GRangesList
#' @importFrom karyoploteR kpPlotRegions
#' @importFrom GenomicRanges reduce
#' @importFrom IRanges IRanges
#' @importFrom karyoploteR kpLines
#' @importFrom GenomeInfoDb seqnames
#' @export

site_Distribution_gene <- function(gtf_file, geneName, seqInfo, gtf_path, gr1, mergeSample = F,
                                   rightmargin = 0.2,ideogramheight = 5, data1height = 200, gene_r0 = 0, gene_r1 = 0.6, mergeTranscripts = T,
                                   non.coding.exons.height=0.5, non.coding.exons.col="#82BB22", non.coding.exons.border.col="#000000",
                                   coding.exons.col="#0050A4", coding.exons.border.col="#000000", introns.col="#000000",
                                   region_col = "#E50053", region_r1 = 0.8, region_r0 = 0,
                                   point1_col = "#E50053", point1_r1 = 0.8, point1_r0 = 0.6, point1_pch = 16, point1_cex = 2, point1_y = 1,
                                   point2_col = "#000000", point2_r1 = 0.8, point2_r0 = 0.6, point2_pch = 1, point2_cex = 2, point2_y = 1, kpPoints = T){

  if(seqInfo == "hg38"){

    seqInfo <- GenomicRanges::GRanges(
      seqnames = c("1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
                   "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                   "chrM", "chrX", "chrY"),
      ranges = IRanges::IRanges(
        start = 1,
        end = c(248956422, 133797422, 135086622, 133275309, 114364328, 107043718,
                101991189, 90338345, 83257441, 80373285, 58617616, 242193529,
                64444167, 46709983, 50818468, 198295559, 190214555, 181538259,
                170805979, 159345973, 145138636, 138394717, 16569, 156040895, 57227415)
      ),
      strand = "*"
    )

  }else if(seqInfo == "hg19"){

    seqInfo <- GenomicRanges::GRanges(
      seqnames = c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21",
                   "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
                   "chrM", "chrX", "chrY"),
      ranges = IRanges::IRanges(
        start = 1,
        end = c(249250621, 135534747, 135006516, 133851895, 115169878, 107349540,
                102531392, 90354753, 81195210, 78077248, 59128983, 243199373,
                63025520, 48129895, 51304566, 198022430, 191154276, 180915260,
                171115067, 159138663, 146364022, 141213431, 16569, 155270560, 59373566)
      ),
      strand = "*"
    )

  }else if(seqInfo == "mm10"){

    seqInfo <- GenomicRanges::GRanges(
      seqnames = c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                   "chr16", "chr17", "chr18", "chr19", "chr2", "chr3", "chr4",
                   "chr5", "chr6", "chr7", "chr8", "chr9", "chrM", "chrX", "chrY"),
      ranges = IRanges::IRanges(
        start = 1,
        end = c(195471971, 130694993, 122082543, 120129022, 120421639, 124902244,
                104043685, 98207768, 94987271, 90702639, 61431566, 182113224,
                160039680, 156508116, 151834684, 149736546, 145441459, 129401213,
                124595110, 16299, 171031299, 91744698)
      ),
      strand = "*"
    )

  }else if(seqInfo == "mm39"){

    seqInfo <- GenomicRanges::GRanges(
      seqnames = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
                   "chrX", "chrY", "chrM"),
      ranges = IRanges::IRanges(
        start = 1,
        end = c(195471971, 182113224, 160039680, 156508116, 151834684, 149736546,
                145441459, 129401213, 124595110, 130694993,
                122082543, 120129022, 120421639, 124902244, 104043685, 98207768,
                94987271, 90702639, 61431566, 171031299, 91744698, 16299)
      ),
      strand = "*"
    )

  }


  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path)
  geneRegion <- gtf_file$Gene_gr[which(gtf_file$Gene_gr$gene_name == geneName),]
  GenomeInfoDb::seqlevels(geneRegion) <- GenomeInfoDb::seqlevelsInUse(geneRegion)

  # 画图区域参数
  pp <- karyoploteR::getDefaultPlotParams(plot.type=1)
  pp$rightmargin <- rightmargin
  pp$ideogramheight <- ideogramheight
  pp$data1height <- data1height

  # 染色体参数
  kp <- karyoploteR::plotKaryotype(genome = seqInfo, zoom = geneRegion, plot.type = 1, plot.params = pp)
  genes.data <- karyoploteR::makeGenesDataFromTxDb(txdb,
                                      karyoplot = kp,
                                      plot.transcripts = T,
                                      plot.transcripts.structure = T)

  # 基因参数
  #genes.data$genes$name <- geneRegion$gene_name
  if(mergeTranscripts == T){genes.data <- karyoploteR::mergeTranscripts(genes.data)} # 转录本是否压缩

  karyoploteR::kpPlotGenes(kp, data = genes.data, r0 = gene_r0, r1 = gene_r1, data.panel = 1,
              non.coding.exons.height = non.coding.exons.height, non.coding.exons.col = non.coding.exons.col, non.coding.exons.border.col = non.coding.exons.border.col,
              add.strand.marks = TRUE, mark.height = 0.2, mark.width = 2, mark.distance = 5,
              coding.exons.col = coding.exons.col, coding.exons.border.col = coding.exons.border.col,
              introns.col = introns.col, marks.col = "#000000",
              border = NULL, avoid.overlapping = TRUE, clipping = TRUE, gene.names = NA)
  # 基因的lable
  karyoploteR::kpAddLabels(kp, labels = geneName, side = "left", data.panel = 1, cex = 1, r0 = gene_r0, r1 = gene_r1)

  gr2 <- lapply(gr1, function(x){
    gr3 <- x[which(x$Gene_Name == geneName)]
    GenomeInfoDb::seqlevels(gr3) <- GenomeInfoDb::seqlevelsInUse(gr3)
    return(gr3)
  })

  gr2 <- GenomicRanges::GRangesList(gr2)
  width_peak <- (region_r1 - region_r0)/length(gr1)
  sample_width <- 0.2 / length(gr1)

  print("gr2 Load")

  if(mergeSample == F){
    for(i in 1:length(gr1)){
      karyoploteR::kpPlotRegions(kp, data = GenomicRanges::reduce(gr2[[i]]), col = region_col, r1 = (region_r0 + i * width_peak - sample_width), r0 = (region_r0 + (i - 1) * width_peak + sample_width), data.panel = 1, border = region_col)
      karyoploteR::kpAddLabels(kp, labels = names(gr2[i]),
                  side = "right", data.panel = 1, cex = 1, r1 = (region_r0 + i * width_peak - sample_width), r0 = (region_r0 + (i - 1) * width_peak + sample_width))

      karyoploteR::kpLines(kp, chr = GenomeInfoDb::seqnames(geneRegion) ,x = (GenomicRanges::start(geneRegion): GenomicRanges::end(geneRegion)), y = 0.5, data.panel = 1,
              r1 = (region_r0 + i * width_peak - sample_width), r0 = (region_r0 + (i - 1) * width_peak + sample_width), col = "black", lty = 2)
      print(i)
    }
  }else if(mergeSample == T){
    gr3 <- unlist(gr2)
    gr3 <- GenomicRanges::reduce(gr3)
    karyoploteR::kpPlotRegions(kp, data = gr3, col = region_col, r1 = region_r1, r0 = region_r0, data.panel = 1, border = region_col)
    karyoploteR::kpAddLabels(kp, labels = "merged samples", side = "right", data.panel = 1, cex = 1, r0 = region_r0, r1 = region_r1)
  }
}


