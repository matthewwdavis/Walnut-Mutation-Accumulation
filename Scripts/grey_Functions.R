## Functions
mutations_embryo_tree<-function(mutations_sub){
  mutations_sub<-mutations_sub[Type=="Embryo"]
  mutations_matrix<-data.table(table(UNIQUE=mutations_sub$UNIQUE, SOURCE=mutations_sub$Source))
  mutations_matrix<-dcast(mutations_matrix, SOURCE~UNIQUE, value.var = "N")
  distance_mat<-dist(mutations_matrix[,-1], method = "binary")
  clust<-hclust(distance_mat, method = "average")
  clust$labels<-mutations_matrix$SOURCE
  clustplot<-plot(clust)
  return(list(dist=distance_mat, plot=clustplot))
}


annot_mut_effect<-function(i){
  genic<-mutations_means$CDS[i]
  if(!genic) return("nonCDS")
  row<-mutations_means[i]
  POS<-row$POS
  REF<-row$REF
  ALT<-row$ALT
  CHROM=row$CHROM
  overlap<-sites_in_features(CDS_good, row, mode = "any")[any==T]
  mod<-unique(CDS_good[ID%in%overlap$ID]$Parent)[1]
  CDS_model<-CDS_good[Parent%in%mod]
  DIR<-CDS_model$DIRECTION[1]
  CDS_POS<-unlist(lapply(1:nrow(CDS_model), function(j){
    CDS_model$START[j]:CDS_model$STOP[j]
  }))
  BPs<-cds[[mod]]
  AA<-seqinr::translate(toupper(cds[[mod]]))
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
    CDS_dt_mutations_table[, gene := names(CDS)[i]]
    return(CDS_dt_mutations_table[, .(REF, ALT, effect, gene)])
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
