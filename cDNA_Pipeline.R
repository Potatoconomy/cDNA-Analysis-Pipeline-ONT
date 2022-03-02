rm(list=ls())

#library(tximportData) ## Later to convert transcripts into genes
library(seqinr)
library(dplyr)
library(tidyverse)
library(data.table)
library(Rsubread)


guppy_path <- "/Users/nes/Documents/SequencingPrograms/ont-guppy-cpu/bin/"
base_path <- "/Users/nes/Documents/SequencingPrograms/Practice_Data1/" # add check for / at end
fq_path <- "FastQ/"
ref_DNA <- "/Users/nes/Documents/SequencingPrograms/ReferenceGenomes/Rat/ENSEMBL_DNA_mRatBN7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa"
ref_cDNA <- "/Users/nes/Documents/SequencingPrograms/ReferenceGenomes/Rat/ENSEMBL_cDNA_mRatBN7/Rattus_norvegicus.mRatBN7.2.cdna.all.fa"
samtools_path <- "/opt/homebrew/bin/samtools"
ref_cDNA_df <- read.fasta(ref_cDNA)
GTF_mRat <- "/Users/nes/Documents/SequencingPrograms/ont-guppy-cpu/bin/Rat/ENSEMBL_GTF/Rattus_norvegicus.mRatBN7.2.105.gtf"

# Can find your samtools path in computer by typing in "which samtools" in terminal following homebrew installation


demultiplex <- function(){  # Later get fq path from sys
  kit <- "SQK-PCB109"
  barcoder <- paste(guppy_path, "/guppy_barcoder", sep="")
  data_path <- paste(base_path, fq_path, sep="")
  print(data_path)
  system(paste(barcoder, 
               "-i", data_path, 
               "-s", paste(base_path, "Demultiplexed/", sep=""),
               "--barcode_kits", kit,
               "--trim_barcodes",
               "--trim_adapters",
               "--trim_primers",
               "--num_extra_bases_trim", 10,
               sep=" "))
}

align_cDNA <- function(){
  aligner <- paste(guppy_path, "minimap2", sep="")
  align_path <- paste(base_path, "Minimap2_cDNA", sep="")
  system(paste("mkdir", align_path))  # add check for if exists
  system(paste("mkdir", paste(align_path, "/Counts_idx", sep="")))
  barcodes <- list.files(paste(base_path, "Demultiplexed", sep=""), pattern="^barcode")
  
  #print("Note: This only works if there is only ONE .fastq file in each barcode folder")
  
  for (i in barcodes){
    output <- paste(paste(align_path, "/", i, sep=""))
    counts_save_path_idxstats <- paste(align_path, "/Counts_idx/", i, "/", sep="")
    system(paste("mkdir", counts_save_path_idxstats))
    system(paste("mkdir", output))
    
    data <- list.files(paste(base_path, "Demultiplexed/", i, sep=""))
    data_ <- paste(base_path, "Demultiplexed/", i, "/", data, sep="")
    true_output_sam <- paste(output, "/", i, ".sam", sep="")
    true_output_bam <- paste(output, "/", i, ".bam", sep="")
    
    system(paste(aligner,
                 "-ax", "map-ont",
                 ref_cDNA,  
                 data_,
                 ">", true_output_sam))
    
    system(paste(samtools_path, "view", "-S", "-b", true_output_sam, ">", true_output_bam))
    system(paste(samtools_path, "sort", true_output_bam, "-o", true_output_bam))
    system(paste(samtools_path, "index", "-b", true_output_bam))
    system(paste("rm", true_output_sam))
    
    #make_saf(output, true_output_bam, i)
    idx_output <- paste(counts_save_path_idxstats, i, ".txt", sep="")
    system(paste(samtools_path, "idxstats", true_output_bam, ">", idx_output))
  }
}

align_DNA <- function(){
  aligner <- paste(guppy_path, "minimap2", sep="")
  align_path <- paste(base_path, "Minimap2_DNA", sep="")
  system(paste("mkdir", align_path))  
  system(paste("mkdir", paste(align_path, "/Counts", sep="")))
  system(paste("mkdir", paste(align_path, "/Counts_idx/", sep="")))

  barcodes <- list.files(paste(base_path, "Demultiplexed", sep=""), pattern="^barcode")
  
  #print("Note: This only works if there is only ONE .fastq file in each barcode folder")
  
  for (i in barcodes){
    counts_save_path <- paste(align_path, "/Counts/", i, "/", sep="")
    counts_save_path_idxstats <- paste(align_path, "/Counts_idx/", i, "/", sep="")
    
    system(paste("mkdir", counts_save_path))
    system(paste("mkdir", counts_save_path_idxstats))
    output <- paste(paste(align_path, "/", i, sep=""))
    system(paste("mkdir", output))
    data <- list.files(paste(base_path, "Demultiplexed/", i, sep=""))
    data_ <- paste(base_path, "Demultiplexed/", i, "/", data, sep="")
    true_output_sam <- paste(output, "/", i, ".sam", sep="")
    true_output_bam <- paste(output, "/", i, ".bam", sep="")
    
    system(paste(aligner,
                 "-ax", "splice",
                 ref_DNA,  
                 data_,
                 ">", true_output_sam))
    
    system(paste(samtools_path, "view", "-S", "-b", true_output_sam, ">", true_output_bam))
    system(paste(samtools_path, "sort", true_output_bam, "-o", true_output_bam))
    system(paste(samtools_path, "index", "-b", true_output_bam))
    system(paste("rm", true_output_sam))
    
    df <- get_counts(true_output_bam, GTF_mRat)
    write.table(df, file=paste(counts_save_path, i, ".csv", sep=""))
    
    idx_output <- paste(counts_save_path_idxstats, i, ".txt", sep="")
    system(paste(samtools_path, "idxstats", true_output_bam, ">", idx_output))
  }
}

make_saf <- function(output_folder, bam_file_name, barcode) {
  print(barcode)
  reads_file_name <- paste(output_folder, "/", barcode, "reads", ".txt", sep="")
  saf_save_path <- paste(output_folder, "/", barcode, "saf", ".saf", sep="")
  system(paste(samtools_path, "idxstats", bam_file_name, ">", reads_file_name ))
  reads <- head(read.table(reads_file_name), -1)
  #reads <- filter(a, a$V3 + a$V4 > 0)
  
  f <- function(i) {
    anno <- attr(ref_cDNA_df[[i]], which= "Annot")
    info <- str_split(anno, fixed(":"))
    #transcript <- str_extract(info[[1]][1], "E[^ ]+")
    chr <- paste("chr", info[[1]][3], sep="")
    start <- 1
    end <- length(ref_cDNA_df[[i]])
    #start <- info[[1]][4]
    #end <- info[[1]][5]
    data <- c(i, chr,start, end, "-") # c(paste("SN:",i,sep=""), chr, start, end, "-")
    return(data)
  }
  
  g <- lapply(reads$V1, f)
  gg <- as.data.frame(t(as.data.table(g)))
  colnames(gg) <- c("GeneID", "Chr", "Start", "End", "Strand")
  write.table(gg, file=saf_save_path,row.names = FALSE, sep="\t", quote=FALSE)
}

get_counts <- function(bam_file, GTF_file) {
  df <- featureCounts(bam_file, annot.ext = GTF_file, isGTFAnnotationFile = TRUE, useMetaFeatures = TRUE, countMultiMappingReads = TRUE)
  return(df$counts)
}

make_count_table_cDNA_idx <- function(count_folder) {
  # Get each column to be a sample and each row to be a gene in the gene name 
  # Need to get gene name from the cDNA fa file and convert the transcript id into that
  
  trans_to_ENSR <- function(i) {
    # Convert ENSRNOT to ENSRNOG using mRatBN7 GTF file
    anno <- attr(ref_cDNA_df[[i]], which= "Annot")
    info <- str_split(anno, fixed(":"))
    #transcript <- str_extract(info[[1]][1], "E[^ ]+")
    ENSR <- str_split(info[[1]][7], fixed(" "))[[1]][1]
    data <- c(i, ENSR)
    return(data)
  }
  
  files <- list.files(path = count_folder, pattern = "*.txt$", recursive = TRUE)
  print(files)
  for (file in files) {
    file_name <- paste(count_folder, "/", file, sep="")
    a <- read.table(file_name)
    a <- select(a_, c("V1","V3"))
    a <- a[a$V3 > 0 ,]
    
    g <- lapply(a$V1, trans_to_ENSR)
    gg <- as.data.frame(t(as.data.table(g)))
    colnames(gg) <- c("Transcript", "Gene")
    
    # Here I need to match up the transcripts to the Genes and get counts in that format
  }
}
make_count_table_cDNA_idx("/Users/nes/Documents/SequencingPrograms/Practice_Data1/Minimap2_cDNA/Counts_idx")


demultiplex()
align_cDNA()
align_DNA()
make_count_table_cDNA_idx("/Users/nes/Documents/SequencingPrograms/Practice_Data1/Minimap2_cDNA/Counts_idx")





