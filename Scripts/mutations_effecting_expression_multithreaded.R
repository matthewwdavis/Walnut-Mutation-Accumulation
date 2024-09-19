#Set up
## Libraries
library(tidyverse)
library(data.table)
library(future.apply)

## Set multithreading rules
plan(multicore, workers = 32)
options(future.globals.maxSize = 128 * 1024^3)

## Functions
parse_columns <- function(table, sample){
  
  # Parse the columns
  parse.table <- table %>%
    mutate(rna = sub("rna-([^=]+).*", "\\1", target_id),
           gene = sub(".*gene=([^=]+).*", "\\1", target_id),
           name = sub(".*name=([^=]+).*", "\\1", target_id),
           seq_id = sub(".*seq_id=([^=]+).*", "\\1", target_id),
           type = sub(".*type=([^=]+)", "\\1", target_id),
           rna = gsub("gene", "", rna),
           gene = gsub("gene-|name", "", gene),
           name = gsub("seq_id", "", name)) %>%
    select(-target_id) %>%
    select(seq_id, everything()) %>%
    mutate(Source = sample)
  
  # Make the chromosome column numeric
  parse.table$seq_id <- as.numeric(gsub("NC_0499.*?([0-9]+).*", "\\1", parse.table$seq_id))
  
  # Add groups for later coloring
  final.table <- parse.table %>%
    mutate(Group = case_when(
      seq_id == 4 ~ "chr_4",
      seq_id == 9 ~ "chr_9",
      TRUE ~ "Other"))
  
  final.table <- final.table %>%
    rename(CHROM = seq_id)
  
  return(final.table)
}

read_genes_gff <- function(file){
  
  gff <- fread(file, skip = 3)
  
  colnames(gff) <- c("CHROM", "SOURCE", "ID", "START", "STOP", "SCORE", "STRAND", "FRAME", "ATTRIBUTE")
  
  # Parse the columns
  parse.table <- gff %>%
    mutate(gene = sub(".*gene-([^=]+).*", "\\1", ATTRIBUTE),
           POS = (STOP + START) / 2) %>%
    select(-ATTRIBUTE) %>%
    filter(!grepl(";", gene))
  
  gff <- parse.table[grep("^NC_0499", parse.table$CHROM), ]
  gff$CHROM <- as.numeric(gsub("NC_0499(..).1_RagTag", "\\1.", gff$CHROM))
  
  gene_pos <- gff %>%
    select(CHROM, START, STOP, gene)
  
  genes <- gene_pos %>%
    mutate(gene = sub("_.*", "", gene))
  
  return(genes)  
}

## Read in data
# SNP data
dn_snps.df <- fread("./all_dn_snps.vcf")

# RNA data
cr11a <- fread("./11A_abundance.tsv")
cr12a <- fread("./12A_abundance.tsv")
cr13a <- fread("./13A_abundance.tsv")
cr15a <- fread("./15A_abundance.tsv")
cr15a1 <- fread("./15A1_abundance.tsv")
cr16a <- fread("./16A_abundance.tsv")
cr16a1 <- fread("./16A1_abundance.tsv")
cr16a2 <- fread("./16A2_abundance.tsv")
cr17a1 <- fread("./17A1_abundance.tsv")
cr17a2 <- fread("./17A2_abundance.tsv")
cr18a <- fread("./18A_abundance.tsv")
cr18a1 <- fread("./18A1_abundance.tsv")

# GFF
genes <- read_genes_gff("./ref_chandler_primary_default_scaffold_chr_only_liftoff_a99s99.gff")

## Filter data for analysis
# Get all unique snp values
uniq_snps <- unique(dn_snps.df$unique)

# filter for only biallelic sites in the stacks, remove largest dataframe
dn_stacks.df <- dn_snps.df %>%
  filter(!Source %in% c("ref_chandler", "tree2", "cr2_hifi", "cr85", "cr10", "cr13", "cr22", "cr21_1", "cr21_2")) %>%
  filter(GT %in% c("0/1", "1/1")) %>%
  filter(grepl("^([^,]*,[^,]*)$", AD)) %>%
  filter(str_length(REF) == 1) %>%
  filter(str_length(ALT) == 1)

rm(dn_snps.df)

# Parse and combine RNA data
cr11a <- parse_columns(cr11a, "cr2_11A")
cr12a <- parse_columns(cr12a, "cr2_12A")
cr13a <- parse_columns(cr13a, "cr2_13A")
cr15a <- parse_columns(cr15a, "cr2_15A")
cr15a1 <- parse_columns(cr15a1, "cr2_15A1")
cr16a <- parse_columns(cr16a, "cr2_16A")
cr16a1 <- parse_columns(cr16a1, "cr2_16A1")
cr16a2 <- parse_columns(cr16a2, "cr2_16A2")
cr17a1 <- parse_columns(cr17a1, "cr2_17A1")
cr17a2 <- parse_columns(cr17a2, "cr2_17A2")
cr18a <- parse_columns(cr18a, "cr2_18A")
cr18a1 <- parse_columns(cr18a1, "cr2_18A1")

all_rna <- rbind(cr11a, cr12a, cr13a, cr15a, cr15a1, cr16a, cr16a1, cr16a2, cr17a1, cr17a2, cr18a, cr18a1)

rm(cr11a, cr12a, cr13a, cr15a, cr15a1, cr16a, cr16a1, cr16a2, cr17a1, cr17a2, cr18a, cr18a1)

# Identify SNP presence in stacks
## Extract only unique SNP entries
# Create summary table of snp presence
snp_pres.df <- dn_stacks.df %>%
  filter(unique %in% uniq_snps) %>%
  group_by(unique) %>%
  summarize(
    CHROM = first(CHROM),
    POS = first(POS),
    unique = first(unique),
    REF = first(REF),
    ALT = first(ALT),
    QUAL = median(QUAL),
    sources_present = paste(unique(Source), collapse = ", ")
  )

# Remove snps present in every sample, in only 1 sample, in all samples but 1
snp_diff.df <- snp_pres.df %>%
  filter(str_count(sources_present, ",") != 11) %>%
  filter(str_count(sources_present, ",") != 10) %>%
  filter(str_count(sources_present, ",") > 0)

## Find differences in expression due to SNPs
# Merge gff and rna data
all_rna <- left_join(all_rna, genes, by = join_by(CHROM, gene))


# Look for difference in TPM where snps are present and where they are not present
all_rna <- data.table(all_rna)

snp_diff.df <- data.table(snp_diff.df)

r = unique(all_rna$rna)

association<-rbindlist(future_lapply(1:nrow(snp_diff.df), function(m){
  
  u = snp_diff.df$unique[m]
  
  samples <- snp_diff.df$sources_present[m]
  samples <- unlist(strsplit(samples, ", "))
  
  gene_assoc<-rbindlist(future_lapply(unique(all_rna$rna), function(r){
    
    rna_sub <- all_rna[rna==r]
    rna_sub$mut <- rna_sub$Source %in% samples
    if(uniqueN(rna_sub$mut)==1 | min(table(rna_sub$mut))==1) return(NULL)
    ttest <- t.test(rna_sub$tpm~rna_sub$mut)
    out <- data.table(unique=u, rna=r, pvalue=ttest$p.value)
    
    return(out)
  }))
  
  return(gene_assoc)
}))

  # Retrieve RNA locations
rna_locs <- all_rna[,.(CHROM=unique(CHROM), START=unique(START)), by=.(rna)]

  # Merge results from loop with RNA locations
gene_assoc_merge <- merge(rna_locs, association, by="rna")

## Write out the file
write_tsv(x = gene_assoc_merge, file = "./mutation_transcript_association.tsv")