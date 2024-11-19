# Setup
## Libraries

library(polymorphology2)
library(tidyverse)
library(pbapply)


## Functions

mutations_embryo_tree<-function(mutations_sub){
  mutations_sub<-mutations_sub[Type=="Embryo"]
  mutations_matrix<-data.table(table(UNIQUE=mutations_sub$UNIQUE,
                                     SOURCE=mutations_sub$Source))
  mutations_matrix<-dcast(mutations_matrix, SOURCE~UNIQUE, value.var = "N")
  distance_mat<-dist(mutations_matrix[,-1], method = "binary")
  clust<-hclust(distance_mat, method = "average")
  clust$labels<-mutations_matrix$SOURCE
  clustplot<-plot(clust)
  return(list(dist=distance_mat, plot=clustplot))
}

annot_mut_effect<-function(i){
  genic <- mutations_means$CDS[i]
  if(!genic) return("nonCDS")
  row <- mutations_means[i]
  POS <- row$POS
  REF <- row$REF
  ALT <- row$ALT
  CHROM = row$CHROM
  overlap <- sites_in_features(cds_good, row, mode = "any")[any==T]
  mod <- unique(cds_good[ID%in%overlap$ID]$Parent)[1]
  CDS_model <- cds_good[Parent%in%mod]
  DIR <- CDS_model$DIRECTION[1]
  CDS_POS <- unlist(lapply(1:nrow(CDS_model), function(j){
    CDS_model$START[j]:CDS_model$STOP[j]
  }))
  BPs <- CDS[[mod]]
  AA <- seqinr::translate(toupper(CDS[[mod]]))
  genome_BPs<-fasta[[CHROM]][CDS_POS]
  genome_AA<-seqinr::translate(toupper(genome_BPs))
  
  REF_genome<-toupper(fasta[[CHROM]][POS])
  if(REF_genome!=REF) return("REF_genome_mismatch")
  
  CDS_dt<-data.table(cPOS=CDS_POS, BP=toupper(cds[[mod]]))
  if(DIR=="-") {
    CDS_dt$cPOS<-rev(CDS_dt$cPOS)
    REF<-revcomp(REF)
    ALT<-revcomp(ALT)
    genome_AA<-seqinr::translate(unlist(strsplit(revcomp(toupper(fasta[[CHROM]][CDS_POS])), split="")))
  }
  if(!all(genome_AA==AA)) return("REF_genome_mismatch_AA")
  if(CDS_dt[cPOS==POS]$BP!=REF) return(CDS_dt[cPOS==POS]$BP)
  CDS_dt$AA<-rep(seqinr::translate(CDS_dt$BP), each=3)
  CDS_dt$BP2<-CDS_dt$BP
  CDS_dt$BP2[which(CDS_dt$cPOS==POS)]<-ALT
  CDS_dt[BP2!=BP]
  CDS_dt$newAA<-rep(seqinr::translate(CDS_dt$BP2), each=3)
  CDS_dt$AAPOS<-rep(1:length(seqinr::translate(CDS_dt$BP)), each=3)
  SYN=all(CDS_dt$AA==CDS_dt$newAA)
  message(i, mod, CDS_dt[which(CDS_dt$cPOS==POS)]$AA," to ", CDS_dt[which(CDS_dt$cPOS==POS)]$newAA)
  message(ifelse( SYN, "\tSyn","\tNon-Syn"))
  
  effect<-CDS_dt[which(CDS_dt$cPOS==POS)]
  effect$gene_model<-mod
  effect$CHROM<-CHROM
  ifelse(SYN, "Syn","Non-Syn")
}

vstrsplit<-function(x,...){
  unlist(strsplit(x,...))
}

CDS_codon_mutation_table <- function(CDS) {
  # Precompute mutations table
  muts <- data.table(table(REF = c("A", "T", "C", "G"), ALT = c("A", "T", "C", "G")))[REF != ALT]
  
  # Create codon table with all possible mutations
  codons <- data.table(SEQINR.UTIL$CODON.AA)
  codons[, CODON := toupper(CODON)]
  
  # Precompute the CODON_mutations_table
  CODON_mutations_table <- rbindlist(lapply(1:nrow(codons), function(i) {
    CODON <- codons$CODON[i]
    CODON_dt <- data.table(REF = strsplit(CODON, split = "")[[1]], POS = 1:3, REFCODON = CODON)
    CODON_dt <- unique(merge(CODON_dt, muts, by = "REF", allow.cartesian = TRUE))
    CODON_dt[, ALTCODON := sapply(1:.N, function(j) {
      REFCODON <- strsplit(CODON_dt$REFCODON[j], split = "")[[1]]
      REFCODON[CODON_dt$POS[j]] <- CODON_dt$ALT[j]
      paste0(REFCODON, collapse = "")
    })]
    CODON_dt[, `:=`(REFAA = codons$AA[match(REFCODON, codons$CODON)],
                    ALTAA = codons$AA[match(ALTCODON, codons$CODON)])]
    CODON_dt[, effect := ifelse(REFAA == ALTAA, "S", "N")]
    return(CODON_dt)
  }))
  
  # Convert each CDS into a mutation table
  pb <- txtProgressBar(min = 0, max = length(CDS), style = 3)
  CDS_dt_all <- rbindlist(lapply(1:length(CDS), function(i) {
    random_seq <- toupper(CDS[[i]])
    remainder <- length(random_seq) %% 3
    if (remainder != 0) {
      message("Gene not multiple of 3")
      return(NULL)
    }
    
    original_protein <- seqinr::translate(random_seq)
    CDS_dt <- data.table(POS = rep(1:3, times = length(random_seq) / 3),
                         DNAPOS = 1:length(random_seq),
                         REF = random_seq,
                         AA = rep(original_protein, each = 3))
    
    # Extract codon for each amino acid
    codons_list <- sapply(seq(1, length(random_seq), by = 3), function(start) {
      paste(random_seq[start:(start + 2)], collapse = "")
    })
    CDS_dt[, REFCODON := rep(codons_list, each = 3)]
    
    # Merge with mutations to add ALT column of all possible mutations at each codon
    CDS_dt_mutations_table <- unique(merge(CDS_dt, muts, by = "REF", allow.cartesian = TRUE))[order(DNAPOS)]
    
    # Merge with CODON_mutations_table to annotate effect
    CDS_dt_mutations_table <- merge(CDS_dt_mutations_table, CODON_mutations_table, by = c("REFCODON", "REF", "POS", "ALT"))
    setTxtProgressBar(pb, i)
    CDS_dt_mutations_table[, GENE := names(CDS)[i]]
    return(CDS_dt_mutations_table[, .(REF, ALT, effect, GENE)])
  }))
  close(pb)
  
  return(CDS_dt_all)
}

mutation_probs<-function(muts, genome, cds=NULL){
  message("Calculating base pair frequencies from genome...")
  bp_freq<-rbindlist(lapply(1:length(genome), function(i){
    message(paste("\tChr",i))
    bp_freq<-data.table(table(REF=genome[[i]]))[REF!="n"]
    bp_freq$REF<-toupper(bp_freq$REF)
    return(bp_freq)
  }))
  bp_freq<-bp_freq[REF %in% c("A","T","C","G"),.(N=sum(N)), by="REF"]
  
  mutations<-data.table((table(REF=muts$REF, ALT=muts$ALT)))[REF!=ALT]
  mutations<-merge(mutations, bp_freq, by="REF")
  mutations$rawPROB=(mutations$N.x/(mutations$N.y))
  
  if(!is.null(cds)){
    message("Calculating base pair frequencies from CDS...")
    bp_freqCDS<-rbindlist(lapply(1:length(cds), function(i){
      
      bp_freq<-data.table(table(REF=cds[[i]]))[REF!="n"]
      bp_freq$REF<-toupper(bp_freq$REF)
      return(bp_freq)
    }))
    bp_freqCDS<-bp_freqCDS[REF %in% c("A","T","C","G"),.(CDSN=sum(N)), by="REF"]
    
    mutations<-merge(mutations, bp_freqCDS, by="REF")
    mutations$rawPROB=(mutations$N.x/(mutations$N.y-mutations$CDSN))
  }
  return(mutations)
}

sample_neutrals<-function(neutrals, N, i){
  samples<-sapply(1:i, function(i){
    samp<-neutrals[sample(1:length(neutrals), N)]
    sum(samp=="N")/sum(samp=="S")
  })
}

GvsIG<-function(mutations_data, REGIONS){
  REGIONS$ID<-1:nrow(REGIONS)
  REGIONS_mut<-sites_in_features(REGIONS,mutations_data, mode = "counts")
  REGIONS$mut<-REGIONS_mut$counts[match(REGIONS$ID, REGIONS_mut$ID)]
  
  REGIONS_means<-REGIONS[,.(mut=sum(mut), LENGTH=sum(MAPP_LENGTH2), mutrate=sum(mut)/sum(MAPP_LENGTH2)), by=REGION]
  chitest<-chisq.test(REGIONS_means[,2:3])
  ratio<-REGIONS_means[REGION=="GENE"]$mutrate/REGIONS_means[REGION=="IG"]$mutrate
  return(data.table(p=chitest$p.value, ratio=ratio))
}

my_minimal<-function(){
  theme_minimal(base_size = 6)+
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.25, "cm"),
          legend.background = element_blank(),
          plot.title = element_text(hjust=0.5, size=6))
}

annot_mut_effect<-function(i){
  genic<-mutations_means$CDS[i]
  if(!genic) return("nonCDS")
  row<-mutations_means[i]
  POS<-row$POS
  REF<-row$REF
  ALT<-row$ALT
  CHROM=row$CHROM
  overlap<-sites_in_features(cds_good, row, mode = "any")[any==T]
  mod<-unique(cds_good[ID%in%overlap$ID]$Parent)[1]
  cds_model<-cds_good[Parent%in%mod]
  DIR<-CDS_model$DIRECTION[1]
  CDS_POS<-unlist(lapply(1:nrow(CDS_model), function(j){
    CDS_model$START[j]:CDS_model$STOP[j]
  }))
  BPs<-protein_cds[[mod]]
  AA<-seqinr::translate(toupper(protein_cds[[mod]]))
  genome_BPs<-fasta[[CHROM]][CDS_POS]
  genome_AA<-seqinr::translate(toupper(genome_BPs))
  
  REF_genome<-toupper(fasta[[CHROM]][POS])
  if(REF_genome!=REF) return("REF_genome_mismatch")
  
  CDS_dt<-data.table(cPOS=CDS_POS, BP=toupper(protein_cds[[mod]]))
  if(DIR=="-") {
    CDS_dt$cPOS<-rev(CDS_dt$cPOS)
    REF<-revcomp(REF)
    ALT<-revcomp(ALT)
    genome_AA<-seqinr::translate(unlist(strsplit(revcomp(toupper(fasta[[CHROM]][CDS_POS])), split="")))
  }
  if(!all(genome_AA==AA)) return("REF_genome_mismatch_AA")
  if(CDS_dt[cPOS==POS]$BP!=REF) return(CDS_dt[cPOS==POS]$BP)
  CDS_dt$AA<-rep(seqinr::translate(CDS_dt$BP), each=3)
  CDS_dt$BP2<-CDS_dt$BP
  CDS_dt$BP2[which(CDS_dt$cPOS==POS)]<-ALT
  CDS_dt[BP2!=BP]
  CDS_dt$newAA<-rep(seqinr::translate(CDS_dt$BP2), each=3)
  CDS_dt$AAPOS<-rep(1:length(seqinr::translate(CDS_dt$BP)), each=3)
  SYN=all(CDS_dt$AA==CDS_dt$newAA)
  message(i, mod, CDS_dt[which(CDS_dt$cPOS==POS)]$AA," to ", CDS_dt[which(CDS_dt$cPOS==POS)]$newAA)
  message(ifelse( SYN, "\tSyn","\tNon-Syn"))
  
  effect<-CDS_dt[which(CDS_dt$cPOS==POS)]
  effect$gene_model<-mod
  effect$CHROM<-CHROM
  ifelse(SYN, "Syn","Non-Syn")
}


# Genome and annotations
## Preparing data
### Read in data

# Mappability
mapp <- read.bedGraph("./ref_chandler_primary_default_scaffold.genmap.bedgraph")

# Genome
fasta <- read.fasta("./ref_chandler_primary_default_scaffold.fasta")

# Annotations
gff <- read.GFF("./ref_chandler_primary_default_scaffold_chr_only_liftoff_a99s99.gff")


### Rename chromosomes

gff$CHROM<-as.numeric(gsub("NC_0499(..).1_RagTag", "\\1", gff$CHROM))
unique(gff$CHROM)

mapp$CHROM<-as.numeric(gsub("NC_0499(..).1_RagTag", "\\1", mapp$CHROM))
mapp <- mapp[!is.na(CHROM)]
unique(gff$CHROM)

names(fasta)<-as.numeric(gsub("NC_0499(..).1_RagTag", "\\1", names(fasta)))
fasta <- fasta[!is.na(names(fasta))]
unique(names(fasta))

## Extract data
### Extract genes and cds sequences 

genes <- gff %>%
  filter(TYPE %in% "gene")

genes$GENE <- gsub("ID=gene-(.+);Dbxref.+", "\\1", genes$INFO)

sum(duplicated(genes$GENE))

cds <- gff %>%
  filter(TYPE %in% "CDS")

cds$Parent <- gsub(".+Parent=(.+);Dbxref.+", "\\1", cds$INFO)
sum(duplicated(genes$Parent))

cds$Product <- gsub(".+product=(.+);protein_id.+", "\\1", cds$INFO)
sum(duplicated(genes$Product))

cds$GENE <- gsub(".+CDS;gene=(.+);product.+", "\\1", cds$INFO)
sum(duplicated(cds$GENE))


cds$ID <- 1:nrow(cds)


### Protein amino acids

proteins <- pblapply(unique(cds$Parent), function(p){
  cds_p <- cds[Parent == p]
  dir <- unique(cds_p$DIRECTION)
  CHROM = cds_p$CHROM[1]
  genome_seq <- paste(sapply(1:nrow(cds_p), function(r){
    START = cds_p$START[r]
    STOP = cds_p$STOP[r]
    paste(fasta[[CHROM]][START:STOP], collapse="")
  }), collapse = "")
  if(dir == "-"){
    genome_seq <- revcomp(genome_seq)
  }
  genome_seq <- unlist(strsplit(genome_seq, split = ""))
  AA_seq<-seqinr::translate(genome_seq)
})

names(proteins) <- unique(cds$Parent)

#proteins
sum(duplicated(names(proteins)))

### Extract protein DNA sequences

protein_cds <- pblapply(unique(cds$Parent), function(p){
  cds_p <- cds[Parent == p]
  dir <- unique(cds_p$DIRECTION)
  CHROM = cds_p$CHROM[1]
  genome_seq <- paste(sapply(1:nrow(cds_p), function(r){
    START = cds_p$START[r]
    STOP = cds_p$STOP[r]
    paste(fasta[[CHROM]][START:STOP], collapse = "")
  }), collapse = "")
  if(dir == "-"){
    genome_seq<-revcomp(genome_seq)
  }
  genome_seq<-unlist(strsplit(genome_seq, split=""))
  return(genome_seq)
})

names(protein_cds)<-unique(cds$Parent)

##protein_cds
sum(duplicated(names(proteins)))

## Checking sequences
### Make sure protein sequences are correct (start with M, ends with *)

prot_check <- rbindlist(pblapply(names(proteins), function(p){
  codon1 <- proteins[[p]][1]
  codonend <- proteins[[p]][length(proteins[[p]])]
  return(data.table(p, codon1, codonend))
}), 
fill=T)

prot_check$CHECK <- prot_check$codon1 == "M" & prot_check$codonend == "*"
table(prot_check$CHECK)

### Only keep valid proteins in annotations

cds_good <- cds[Parent %in% prot_check[CHECK==T]$p]
genes_good <- genes[GENE %in% cds_good$GENE]

nrow(cds)
nrow(cds_good)

nrow(genes)
nrow(genes_good)

## Intergenic space and mappability
### Identify intergenic space

intergenic <- genes_good[, .(
  TYPE = "intergenic",
  START = STOP + 1L,
  STOP = data.table::shift(START, type = "lead") - 1L
), by = CHROM]

intergenic <- intergenic[!is.na(STOP)]

intergenic <- intergenic[START <= STOP]

intergenic$ID<-1:nrow(intergenic)

### Defining mappable regions

# Genic space
genes_good$ID<-1:nrow(genes_good)

genes_good$LENGTH <- genes_good$STOP - genes_good$START

mapp_genic<-features_in_features(genes_good, mapp[DEPTH==1], mode = "length")

# Intergenic space
intergenic$ID <- 1:nrow(intergenic)

intergenic$LENGTH <- intergenic$STOP - intergenic$START

mapp_genic <- features_in_features(intergenic, mapp[DEPTH==1], mode = "length")

### Ratio of gene to intergenic

sum(genes_good$LENGTH)/sum(intergenic$LENGTH)

## Identify Homopolymers
### Generate Homopolymer data

# Identify HPs in the genome
HPs <- make_homopolymer(fasta, 3)
fwrite(HPs,"./ref_chandler_primary_omnic.fasta.HPs.tsv", sep = "\t")

### Read homopolymer data in

#HPs <- fread("./Data/Homopolymers/ref_chandler_primary_omnic.fasta.HPs.tsv")

# Mutations
## Preparing data
### Read data

mutations <- fread("./Data/VCFs/Filtered/all_dn_snps.vcf")
# mutations_ns <-fread("./Data/VCFs/Filtered/dn_snps_no_share.vcf")

### Rename default column

mutations <- mutations %>%
  rename("SAMPLE"= "default", "SOURCE" = "Source") %>%
  dplyr::select(-c(ID, INFO, Hap, unique))

### Adjusting the unique column and identifying SNPs

# Mutations
mutations$UNIQUE <- paste(mutations$CHROM, mutations$POS, mutations$ALT, sep = "_")
mutations[, NUMIND := .N, by = UNIQUE] # NUMIND column is the number of Sources with that mutation
SNPS <- nchar(mutations$REF) == 1 & nchar(mutations$ALT) == 1
mutations$SNP <- SNPS

## Annotate the mutations with CDS and Genic regions
### Annotating with CDS

# Mutations
mutations$ID <- 1:nrow(mutations)
mutations_cds <- features_in_sites(cds_good, mutations)
mutations$CDS <- mutations_cds$overlaps[match(mutations$ID, mutations_cds$ID)]

table(mutations$CDS)

### Annotating with genes

# Mutations
mutations_gene <- features_in_sites(genes_good, mutations)
mutations$GENIC <- mutations_gene$overlaps[match(mutations$ID, mutations_gene$ID)]

table(mutations$GENIC)

## Split AD

mutations$ALTDP <- sapply(mutations$AD, function(x) as.numeric(vstrsplit(x, split=",")[2]))
mutations$REFDP <- sapply(mutations$AD, function(x) as.numeric(vstrsplit(x, split=",")[1]))


## Calculation of deviation from heterozygous expection

# Does the Chi square test warning matter?
mutations$CHIPHET <- pbsapply(1:nrow(mutations), function(i) {
  
  ALTDP <- mutations$ALTDP[1]
  REFDP <- mutations$REFDP[1]
  chitest <- chisq.test(c(ALTDP, REFDP))
  
  return(chitest$p.value)
})

## Distribution of shared SNPs

mutations %>%
  ggplot(aes(x = as.factor(NUMIND))) +
  geom_histogram(stat = "count") +
  theme_classic()

counts<-data.table(table(table(mutations$UNIQUE)))
ggplot(counts, aes(x=as.numeric(V1), y=N))+
  geom_bar(stat="identity")+
  scale_y_log10()

## Recalculating DP and VAF
### Calculate DP based on informative reads and normalize

mutations$DPREAL <- mutations$ALTDP + mutations$REFDP

mutations[, DPNORM:=log(DP/mean(DP)), by = SOURCE]
mutations[, DPNORM_REAL:=log(DPREAL/mean(DPREAL)), by = SOURCE]

### Calculate VAF based on informative reads

mutations$VAFREAL <- mutations$ALTDP/mutations$DPREAL

## Assign labels for samples

mutations <- mutations %>%
  mutate(CULTURETYPE = case_when(
    SOURCE %in% c("ref_chandler", "tree2") ~ "tree",
    SOURCE %in% c("cr2_hifi", "cr21_1", "cr21_2", "cr2_11A", "cr2_12A", "cr2_13A", "cr2_15A",
                  "cr2_15A1", "cr2_16A", "cr2_16A1", "cr2_16A2", "cr2_17A1", "cr2_17A2",
                  "cr2_18A", "cr2_18A1") ~ "embryo",
    SOURCE %in% c("cr85", "cr10", "cr13", "cr22") ~ "shoot")) %>%
  mutate(SEQTYPE = case_when(
    SOURCE %in% c("ref_chandler", "tree2", "cr2_hifi") ~ "long_read",
    SOURCE %in% c("cr21_1", "cr21_2", "cr2_11A", "cr2_12A", "cr2_13A", "cr2_15A", "cr2_15A1", 
                  "cr2_16A", "cr2_16A1","cr2_16A2", "cr2_17A1", "cr2_17A2", "cr2_18A",
                  "cr2_18A1", "cr85", "cr10", "cr13", "cr22") ~ "short_read"))

## Calculate summary stats for each SNP across all samples

# Reset object as data table to avoid warnings
mutations <- as.data.table(mutations)
# Mean VAF
mutations[, meanVAF:=mean(VAF), by = UNIQUE]
mutations[, meanVAFREAL:=mean(VAFREAL), by = UNIQUE]

# Median depth
mutations[, medDP:=median(DP), by = UNIQUE]
mutations[, medDPREAL:=median(DPREAL), by = UNIQUE]

# Median Quality
mutations[, medQUAL:=median(QUAL), by = UNIQUE]

## Annotate mappability to means

mutations_MAPP <- features_in_sites(mapp[DEPTH==1], mutations)
mutations$MAPP <- mutations_MAPP$overlaps[match(mutations$ID, mutations_MAPP$ID)]
table(mutations$MAPP, mutations$SEQTYPE)

## Annotate mutations near homopolymeric regions 
# Seems to be about 40 minutes 
mutations_neighbors <- rbindlist(pblapply((1:nrow(mutations)), function(i){ 
  CHROM = mutations$CHROM[i] 
  POS = mutations$POS[i] 
  UNIQUE = mutations$UNIQUE[i] 
  up = paste0(toupper(fasta[[CHROM]][(POS-10):(POS-1)]), collapse="") 
  down = paste0(toupper(fasta[[CHROM]][(POS+1):(POS+10)]), collapse="") 
  return(data.table(UNIQUE, up, down)) 
})) 

mutations <- merge(mutations, mutations_neighbors, by="UNIQUE") 

## Annotate trimer context for mutations

# mutations contexts are only created for SNPs with this function, not INDELs. INDELS are filled with NAs.
# Takes around 15 minutes
mutations$TRICONTEXT <- tricontexts(mutations, fasta)
mutations$MUT <- paste0(substr(mutations$TRICONTEXT, 2, 2), substr(mutations$TRICONTEXT, 4, 4), substr(mutations$TRICONTEXT, 6, 6))

mutations$MUT[mutations$MUT == "NANANA"] <- NA

## Create a means dataframe

mutations_means <- mutations[
                           .(CHROM=unique(CHROM),
                             POS=unique(POS),
                             REF=unique(REF),
                             ALT=unique(ALT),
                             CDS=unique(CDS),
                             GENIC=unique(GENIC)),
                           by=.(UNIQUE)]


# annotate effects of mutations
library(pbapply)
mutations_means$effect<-pbsapply(1:nrow(mutations_means), annot_mut_effect)
mutations_means$effect[mutations_means$effect=="nonCDS" & mutations_means$GENIC==T]<-"genic-nonCDS"

# add back to mutations data
mutations$effect<-mutations_means$effect[match(mutations$UNIQUE, mutations_means$UNIQUE)]

# Write out new files
fwrite(mutations, "./mutations.tsv", sep = "\t")
fwrite(cds_good, "./cds_good.tsv", sep = "\t")
fwrite(genes_good, "./genes_good.tsv", sep = "\t")


