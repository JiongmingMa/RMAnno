#' Annotating region information
#'
#' This function can annotate region information, including 5'UTR, 3'UTR, Intron
#' and Exon.
#'
#' @param anno_info Region type (5'UTR, 3'UTR, Intron or Exon)
#' @param column_name Column name for region information
#' @param grange_file Sites/Peak information file in GRange format
#' @param selected_gr Region annotation information in GRange format
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors mcols
#' @return A GRange file annotated with the Region column
.AnnoType <-
  function(anno_info,
           column_name,
           grange_file,
           selected_gr) {

    # find overlapping region
    overlap_list <- GenomicRanges::findOverlaps(grange_file,selected_gr)
    qhit <- unique(S4Vectors::queryHits(overlap_list))

    # Annotate Type column and check whether the column is NA, if column is not
    # NA, paste anno_info and separate by "/".
    na_row <- is.na(S4Vectors::mcols(grange_file)[, column_name][qhit])
    S4Vectors::mcols(grange_file)[, column_name][qhit][na_row] <- anno_info

    if (any(!na_row)) {
      S4Vectors::mcols(grange_file)[, column_name][qhit][!na_row] <-
        paste(S4Vectors::mcols(grange_file)[, column_name][qhit][!na_row], anno_info, sep = "/")
    }

    return(grange_file)
  }


#' Annotating intergenic region
#'
#' This function can annotate intergenic region information.
#'
#' @param column_name Column name for region information
#' @param grange_file Sites/Peak information file in GRange format
#' @param selected_gr Region annotation information in GRange format
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @return A GRange file annotated with the Region column
.AnnoIntergenic <- function(column_name, grange_file, selected_gr) {

  # find overlapping region and not overlapping region
  overlap_list <- GenomicRanges::findOverlaps(grange_file, selected_gr)
  qhit <- unique(S4Vectors::queryHits(overlap_list))
  intergenic_overlap <- !(seq_along(grange_file) %in% qhit)

  # check whether intergenic region overlap with other region (5'UTR, 3'UTR, Exon, Intron)
  if (all(is.na(S4Vectors::mcols(grange_file)[, column_name][intergenic_overlap]))) {
    S4Vectors::mcols(grange_file)[, column_name][intergenic_overlap] <- "Intergenic"
  } else{
    message("Warning 1: Intergenic region overlap with other region. Please check your GTF/GFF or txdb file.")
    print(grange_file[intergenic_overlap][!is.na(S4Vectors::mcols(grange_file)[, column_name][intergenic_overlap])])
    grange_file$Region[which(is.na(S4Vectors::mcols(grange_file)[, column_name]))] <- "Intergenic"

  }
  return(grange_file)
}


#' Annotating gene information
#'
#' This function can annotate basic gene information, including Gene_id,
#' Gene_name and Bio_type
#'
#' @param grange_file Sites/Peak information file in GRange format
#' @param GTF_type Which database of GTF/GFF file from
#' @param selected_gr Gene annotation information in GRange format
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @return A GRange file annotated with the gene id, gene name and biotype.
.AnnoColumn <- function(grange_file, selected_gr, GTF_type){

  # find overlapping in gene level
  overlap_list <- GenomicRanges::findOverlaps(grange_file, selected_gr)

  # gene_id column
  split_list <- S4Vectors::split(selected_gr$gene_id[S4Vectors::subjectHits(overlap_list)], S4Vectors::queryHits(overlap_list))
  merge_list <- sapply(split_list, function(x) paste(x, collapse = "/"))
  grange_file$Gene_id[as.numeric(names(merge_list))] <- merge_list
  rm(split_list, merge_list)
  message("Gene_id column: Finish!")

  # UCSC GTF/GFF does not have other information
  if (GTF_type != "gtf_ucsc_GCA"){

  # gene_name column
  split_list <- split(selected_gr$gene_name[S4Vectors::subjectHits(overlap_list)], S4Vectors::queryHits(overlap_list))
  merge_list <- sapply(split_list, function(x) paste(x, collapse = "/"))
  grange_file$Gene_name[as.numeric(names(merge_list))] <- merge_list
  rm(split_list, merge_list)
  message("Gene_name column: Finish!")

  # bio_type column
  split_list <- split(selected_gr$gene_biotype[S4Vectors::subjectHits(overlap_list)], S4Vectors::queryHits(overlap_list))
  merge_list <- sapply(split_list, function(x) paste(x, collapse = "/"))
  grange_file$Bio_type[as.numeric(names(merge_list))] <- merge_list
  rm(split_list, merge_list)
  message("Bio_type column: Finish!")
  }

  return(grange_file)
}





#' Annotating other information
#'
#' This function can annotate other information, including RBP,
#' miRNA and SNP, which you need to provide the grange object for function.
#'
#' @param grange_file Sites/Peak information file in GRange format
#' @param other_gr The grange object contain other information
#' @param columns The colnames that you want to annotate from other_gr
#'
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom stringr str_detect
#' @return A GRange file annotated with the other information.
#' @export

infoAnno <- function(grange_file, other_gr, columns){

  if(all(GenomeInfoDb::seqlevels(grange_file) %in% GenomeInfoDb::seqlevels(other_gr))){
    message("chromosome check finished!")
  }else{
    message("Warning: Some chromosomes are not in the target dataset!
         You need to chech your chromosome names.")
  }

  if(all(!GenomeInfoDb::seqlevels(grange_file) %in% GenomeInfoDb::seqlevels(other_gr))){
    message("Warning: All the chromosomes don't correspond, try to change the format automatically.")

    if(all(stringr::str_detect(GenomeInfoDb::seqlevels(grange_file),"chr"))){
      if(all(!stringr::str_detect(GenomeInfoDb::seqlevels(other_gr),"chr"))){

        GenomeInfoDb::seqlevels(other_gr) <- paste0("chr", GenomeInfoDb::seqlevels(other_gr))
        GenomeInfoDb::seqlevels(other_gr)[which(GenomeInfoDb::seqlevels(other_gr)=="chrMT")] <- "chrM"

      }else{
        message("Chromosome corresponding Error! Please check data.")
      }

    }else if(all(!stringr::str_detect(GenomeInfoDb::seqlevels(grange_file),"chr"))){
      if(all(stringr::str_detect(GenomeInfoDb::seqlevels(other_gr),"chr"))){

        GenomeInfoDb::seqlevels(other_gr) <- sub("^chr", "", GenomeInfoDb::seqlevels(other_gr))
        GenomeInfoDb::seqlevels(other_gr)[which(GenomeInfoDb::seqlevels(other_gr)=="M")] <- "MT"

      }else{
        message("Chromosome corresponding Error! Please check data.")
      }
    }
  }

  if(length(columns) == 1){

    index1 <- GenomicRanges::findOverlaps(grange_file, other_gr)
    S4Vectors::mcols(grange_file)[,columns] <- NA
    S4Vectors::mcols(grange_file)[,columns][S4Vectors::queryHits(index1)] <-
      S4Vectors::mcols(other_gr)[,columns][S4Vectors::subjectHits(index1)]

  }else if(length(columns) > 1){

    index1 <- GenomicRanges::findOverlaps(grange_file, other_gr)
    S4Vectors::mcols(grange_file)[,columns] <- NA
    S4Vectors::mcols(grange_file)[,columns][S4Vectors::queryHits(index1),] <-
      S4Vectors::mcols(other_gr)[,columns][S4Vectors::subjectHits(index1),]

  }


  return(grange_file)
}




#' Gene annotation
#'
#' This function can perform corresponding gene annotations according to the
#' site/peak information and GTF files provided by users. Annotation information
#' includes region, gene id, gene name and biotype information.
#'
#' @param grange_file Sites/Peak information file in GRange format
#' @param gr_list A list contains genomic information extracted by ExtractGTFFile
#' @param region The vector of gene regions
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#' @importFrom GenomicFeatures intronsByTranscript
#' @return A annotated GRange file
#' @export
rmAnno <- function(grange_file,
                   gr_list,
                   region = c("5'UTR", "3'UTR", "Exon", "Intron", "Intergenic")){

  # Generate annotation file
  if (is.character(gr_list)){
    message("Import GTF file for annotation ...")
    gr_list <- ExtractGTFFile(gr_list)
  }else if (inherits(gr_list, "list") & length(gr_list) == 8){
    message("Import GRangeList for annotation ...")
    gr_list <- gr_list
  }else if (inherits(gr_list, "TxDb")){
    GTF_type <- "TxDb"
    FiveUTR_gr <- unlist(GenomicFeatures::fiveUTRsByTranscript(gr_list))
    ThreeUTR_gr <- unlist(GenomicFeatures::threeUTRsByTranscript(gr_list))
    Exon_gr <- GenomicFeatures::exons(gr_list)
    Intron_gr <- unlist(GenomicFeatures::intronsByTranscript(gr_list))
    Gene_gr <- unlist(GenomicFeatures::genes(gr_list, single.strand.genes.only = FALSE))
    Gene_gr$gene_id <- names(Gene_gr)
    seqlevels = GenomeInfoDb::seqlevels(gr_list)

    gr_list <- base::list(GTF_type = GTF_type,
                          FiveUTR_gr = FiveUTR_gr,
                          ThreeUTR_gr = ThreeUTR_gr,
                          Exon_gr = Exon_gr,
                          Intron_gr = Intron_gr,
                          Gene_gr = Gene_gr,
                          seqlevels = GenomeInfoDb::seqlevels(gr_list))
  }

  # Check the chromosome name format between GFF file and site information
  if (all(GenomeInfoDb::seqlevels(grange_file) %in% gr_list[[8]])){
    message("Checking the chromosome names: Identical.")
  }else{
    message("Error 1: The chromosome format between site/peak file and GTF/GFF file is not identical!
         Please revise chromosome name in site/peak file or GTF file.")
  }

  # Annotation for column
  grange_file$Region <- NA
  grange_file$Gene_id <- NA
  grange_file$Gene_name <- NA
  grange_file$Bio_type <- NA
  #grange_file$transcript_id <- NA

  # Annotate region column
  column_name = "Region"
  message("Annotating region column.")

  for (i in region){
    if (i == "5'UTR"){
      grange_file <- .AnnoType(i, column_name, grange_file, gr_list[[2]])

    }else if (i == "3'UTR"){
      grange_file <- .AnnoType(i, column_name, grange_file, gr_list[[3]])

    }else if (i == "Exon"){
      grange_file <- .AnnoType(i, column_name, grange_file, gr_list[[4]])

    }else if (i == "Intron"){
      grange_file <- .AnnoType(i, column_name, grange_file, gr_list[[5]])

    }else if (i == "Intergenic"){
      grange_file <- .AnnoIntergenic(column_name, grange_file, gr_list[[6]])

    }
  }
  rm(column_name, i)

  message("Annotating transcript.")
  # find overlapping in transcript. level
  overlap_list <- GenomicRanges::findOverlaps(grange_file, gr_list[[7]])

  # Transcript_id column
  grange_file$transcript_id <- NA
  split_list <- S4Vectors::split(gr_list[[7]]$transcript_id[S4Vectors::subjectHits(overlap_list)], S4Vectors::queryHits(overlap_list))
  merge_list <- sapply(split_list, function(x) paste(x, collapse = "/"))
  grange_file$transcript_id[as.numeric(names(merge_list))] <- merge_list
  rm(split_list, merge_list)
  message("Transcript_id column: Finish!")

  # Annotate Gene_id column
  message("Annotating gene_id, gene_name and gene_biotype columns ...")
  grange_file <- .AnnoColumn(grange_file, gr_list[[6]], gr_list[[1]])

  #all(is.na(grange_file$Gene_id[which(grange_file$Region == "Intergenic")]))   # check whether intergenic sites have gene id
  #all(is.na(grange_file$Gene_name[which(grange_file$Region == "Intergenic")]))
  #all(is.na(grange_file$Bio_type[which(grange_file$Region == "Intergenic")]))

  message("All steps have been completed!")
  return(grange_file)
}




