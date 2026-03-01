#' Extract genomic information from GTF/GFF file.
#'
#' ExtractGTFFile could prepare the genomic information for the annotation step.
#' It extracts the region information including 5'UTR, 3'UTR, exon and intron
#' information. Moreover, it could capture the information of genes including
#' gene name, gene id and biotype if GTF/GFF files have that column.
#'
#' @param gtf A file path of GTF/GFF file
#'
#' @return A list which contain genomic information
#' @export
#' @importFrom IRanges psetdiff
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom S4Vectors mcols
#' @importFrom S4Vectors split
#' @importFrom stringr str_remove
#' @importFrom stringr str_detect
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#' @importFrom GenomicFeatures fiveUTRsByTranscript
#' @importFrom GenomicFeatures threeUTRsByTranscript
#' @importFrom GenomicFeatures exons
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom GenomicFeatures genes
#' @importFrom GenomicFeatures transcripts
ExtractGTFFile <- function(gtf) {

  message("Importing GTF/GFF file.")
  gtf_gr_file <- rtracklayer::import(gtf)

  # GTF information check
  gtf_ensembl_GCA <- c("type", "gene_id", "gene_name", "gene_biotype")
  gff3_ensembl_GCA <- c("type", "gene_id", "Name", "biotype")

  gff3_ncbi_GCF <- c("type", "Name", "gbkey", "gene_biotype", "gene", "Dbxref")
  gtf_ncbi_GCF <- c("type", "gbkey", "gene_biotype", "gene")

  gtf_ucsc_GCA <- c("type", "gene_id", "gene_name")
  #gtf_ucsc_GCA_knownGene <- c("type", "gene_id", "gene_name")
  #gtf_ucsc_GCA_ncbiRefSeq <- c("type", "gene_id", "gene_name")
  #gtf_ucsc_GCA_refGene <- c("type", "gene_id", "gene_name")

  GTF_col <- base::colnames(S4Vectors::mcols(gtf_gr_file))
  message("Matching file type of GTF/GFF file.")

  if (all(gff3_ncbi_GCF %in% GTF_col)) {
    # 6 colunm names gff3_ncbi_GCF
    GTF_type <- "gff3_ncbi_GCF"

    # NCBI GTF does not have region information. It will be replaced by other function
    # within GenomicFeature package.

    # Prepare for region annotation information
    txdb_region <- GenomicFeatures::makeTxDbFromGRanges(gtf_gr_file)
    FiveUTR_gr <- unlist(GenomicFeatures::fiveUTRsByTranscript(txdb_region))
    ThreeUTR_gr <- unlist(GenomicFeatures::threeUTRsByTranscript(txdb_region))
    Exon_gr <- GenomicFeatures::exons(txdb_region)
    Intron_gr <- unlist(GenomicFeatures::intronsByTranscript(txdb_region))

    # Change column names
    gene_name <- which(names(S4Vectors::mcols(gtf_gr_file)) == "Name")
    names(S4Vectors::mcols(gtf_gr_file))[gene_name] <- "gene_name"
    gene_biotype <- which(names(S4Vectors::mcols(gtf_gr_file)) == "Dbxref")
    names(S4Vectors::mcols(gtf_gr_file))[gene_biotype] <- "gene_id"

    # Gene information: There are some duplicated row within. This part need to improve.
    Gene_gr <- gtf_gr_file[gtf_gr_file$type == "gene" | gtf_gr_file$type == "pseudogene"]
    Gene_gr$gene_id <- stringr::str_replace(sapply(Gene_gr$gene_id, "[[", 1), "GeneID:", "NCBI:")
    Transcript_gr <- GenomicFeatures::transcripts(txdb_region)

  } else if (all(gtf_ncbi_GCF %in% GTF_col)) {
    # 4 colunm names gtf_ncbi_GCF
    GTF_type <- "gtf_ncbi_GCF"

    # NCBI GTF does not have region information. It will be replaced by other function
    # within GenomicFeature package.

    # Prepare for region annotation information
    # txdb_region <- GenomicFeatures::makeTxDbFromGRanges(gtf_gr_file)
    # GFF3 file in NCBI might has some issue which need to check and improve.

    txdb_region <- GenomicFeatures::makeTxDbFromGFF(gtf)
    FiveUTR_gr <- unlist(GenomicFeatures::fiveUTRsByTranscript(txdb_region))
    ThreeUTR_gr <- unlist(GenomicFeatures::threeUTRsByTranscript(txdb_region))
    Exon_gr <- GenomicFeatures::exons(txdb_region)
    Intron_gr <- unlist(GenomicFeatures::intronsByTranscript(txdb_region))

    # Change column names
    gene_name <- which(names(S4Vectors::mcols(gtf_gr_file)) == "gene")
    names(S4Vectors::mcols(gtf_gr_file))[gene_name] <- "gene_name"

    # Gene information: There are some duplicated row within. This part need to improve.
    Gene_gr <- gtf_gr_file[gtf_gr_file$type == "gene"]
    Transcript_gr <- gtf_gr_file[gtf_gr_file$type == "transcript"]

  } else if (all(gtf_ensembl_GCA %in% GTF_col)) {
    # 4 colunm names gtf_ensembl_GCA
    GTF_type <- "gtf_ensembl_GCA"

    # Prepare for region annotation information
    FiveUTR_gr <- gtf_gr_file[gtf_gr_file$type == "five_prime_utr"][, c("type", "transcript_id")]
    ThreeUTR_gr <- gtf_gr_file[gtf_gr_file$type == "three_prime_utr"][, c("type", "transcript_id")]
    Exon_gr <- gtf_gr_file[gtf_gr_file$type == "exon"][, c("type", "transcript_id")]

    # Extract Intron information
    Transcript_gr <- gtf_gr_file[gtf_gr_file$type == "transcript"][, c("type", "transcript_id")]
    Transcript_gr <- Transcript_gr[base::order(Transcript_gr$transcript_id)]
    Transcript_id_select <- base::unique(Transcript_gr$transcript_id[Transcript_gr$transcript_id %in% Exon_gr$transcript_id])
    Transcript_gr <- Transcript_gr[Transcript_gr$transcript_id %in% Transcript_id_select]
    Exon_split <- Exon_gr[Exon_gr$transcript_id %in% Transcript_id_select]
    Exon_split <- S4Vectors::split(Exon_split, Exon_split$transcript_id)
    Intron_gr <- IRanges::psetdiff(Transcript_gr, Exon_split)
    Intron_gr <- base::unlist(Intron_gr)

    # Gene information
    Gene_gr <- gtf_gr_file[gtf_gr_file$type == "gene"]

  } else if (all(gff3_ensembl_GCA %in% GTF_col)) {
    # 4 colunm names gff3_ensembl_GCA
    GTF_type <- "gff3_ensembl_GCA"

    # Change column names
    gene_name <- which(names(S4Vectors::mcols(gtf_gr_file)) == "Name")
    names(S4Vectors::mcols(gtf_gr_file))[gene_name] <- "gene_name"
    gene_biotype <- which(names(S4Vectors::mcols(gtf_gr_file)) == "biotype")
    names(S4Vectors::mcols(gtf_gr_file))[gene_biotype] <- "gene_biotype"

    # Prepare for region annotation information
    FiveUTR_gr <- gtf_gr_file[gtf_gr_file$type == "five_prime_UTR"][, c("type", "Parent")]
    ThreeUTR_gr <- gtf_gr_file[gtf_gr_file$type == "three_prime_UTR"][, c("type", "Parent")]
    Exon_gr <- gtf_gr_file[gtf_gr_file$type == "exon"][, c("type", "Parent")]
    Exon_gr$Parent <- stringr::str_remove(unlist(Exon_gr$Parent), "transcript:")

    # Extract Intron information
    Transcript_gr <- gtf_gr_file[which(stringr::str_detect(gtf_gr_file$ID, "transcript:"))][, c("type", "ID", "transcript_id")]
    Transcript_gr <- Transcript_gr[order(Transcript_gr$transcript_id)]
    Exon_split <- S4Vectors::split(Exon_gr, Exon_gr$Parent)

    # Check whether there are duplicated transcript_id
    if (all(Transcript_gr$transcript_id == base::names(Exon_split))) {
      Intron_gr <- IRanges::psetdiff(Transcript_gr, Exon_split)
      Intron_gr <- base::unlist(Intron_gr)
    } else{
      Transcript_gr <- Transcript_gr[!duplicated(Transcript_gr$transcript_id) &
                        Transcript_gr$transcript_id %in% names(Exon_split)]
      Exon_split <- Exon_split[names(Exon_split) %in% Transcript_gr$transcript_id]
      Intron_gr <- IRanges::psetdiff(Transcript_gr, Exon_split)
      Intron_gr <- base::unlist(Intron_gr)
    }

    # Gene information
    Gene_gr <- gtf_gr_file[which(stringr::str_detect(gtf_gr_file$ID, "gene:"))]

  } else if (all(gtf_ucsc_GCA %in% GTF_col)) {
    # 3 colunm names gtf_ucsc_GCA
    GTF_type <- "gtf_ucsc_GCA"

    # UCSC GTF does not have region and gene information. It will be replaced by other function
    # within GenomicFeature package.

    # Prepare for region annotation information
    txdb_region <- GenomicFeatures::makeTxDbFromGRanges(gtf_gr_file)

    FiveUTR_gr <- unlist(GenomicFeatures::fiveUTRsByTranscript(txdb_region))
    ThreeUTR_gr <- unlist(GenomicFeatures::threeUTRsByTranscript(txdb_region))
    Exon_gr <- GenomicFeatures::exons(txdb_region)
    Intron_gr <- unlist(GenomicFeatures::intronsByTranscript(txdb_region))

    # Gene information: There are some duplicated row within. This part need to improve.
    Gene_gr <- GenomicFeatures::genes(txdb_region)
    Transcript_gr <- GenomicFeatures::transcripts(txdb_region)

  }

  gr_list <- base::list(GTF_type = GTF_type,
                        FiveUTR_gr = FiveUTR_gr,
                        ThreeUTR_gr = ThreeUTR_gr,
                        Exon_gr = Exon_gr,
                        Intron_gr = Intron_gr,
                        Gene_gr = Gene_gr,
                        Transcript_gr = Transcript_gr,
                        seqlevels = GenomeInfoDb::seqlevels(gtf_gr_file))

  return(gr_list)
}
