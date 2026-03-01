#' Calculate Splice Sites
#'
#' This function can calculate splicing sites
#'
#' @param gtf_file A file path of GTF/GFF file
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges resize
#' @return A GRange object of splicing sites
#' @export
#'
calculate_splice_sites <- function(gtf_file) {
  splice_sites <- GenomicRanges::GRanges()
  gtf <- rtracklayer::import(gtf_file)
  exons <- gtf[gtf$type == "exon"]

  # 移除只有一个exon的转录本
  freq <- base::table(exons$transcript_id)
  unique_strings <- base::names(freq)[freq == 1]
  exons <- exons[which(!exons$transcript_id %in% unique_strings),]
  exons <- exons[order(exons$transcript_id, GenomicRanges::start(exons))]
  rm(freq)

  # 计算donor位置
  exons$Exon_rank <- paste0(exons$transcript_id,"_", exons$exon_number)
  freq <- table(exons$transcript_id)
  freq <- paste0(names(freq),"_", as.numeric(freq))

  # 移除末尾的Exon
  donor_pos <- exons[!exons$Exon_rank %in% freq]
  donor_pos <- donor_pos[,c("gene_id","gene_name","gene_biotype","transcript_id","transcript_name","transcript_biotype","exon_number","exon_id","ccds_id")]
  donor_pos <- GenomicRanges::resize(donor_pos, width = 2, fix = "end")
  donor_pos$splice_site_type  <- "5'donor"

  # 移除开头的Exon
  acceptor_pos <- exons[!exons$exon_number == "1"]
  acceptor_pos <- acceptor_pos[,c("gene_id","gene_name","gene_biotype","transcript_id","transcript_name","transcript_biotype","exon_number","exon_id","ccds_id")]
  acceptor_pos <- GenomicRanges::resize(acceptor_pos, width = 2, fix = "start")
  acceptor_pos$splice_site_type  <- "3'acceptor"

  splice_sites <- c(donor_pos, acceptor_pos)
  return(splice_sites)
}
