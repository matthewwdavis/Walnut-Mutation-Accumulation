library(polymorphology2)
source("code/R/Functions.R")

# load reference genome ---------------------------------------------------
# Mappability (genmap map -I ref_chandler_primary_omnic.genmap -O genmap_results/ref_chandler_primary_omnic.genmap.K150E0 -bg -K 150 -E 0  )
mapp<-read.bedGraph("matt_data/all/genmap_results/ref_chandler_primary_omnic.genmap.K150E0.bedgraph")

# Primary assembly and annotation
fasta<-read.fasta("matt_data/all/ref_chandler_primary_omnic.fasta")

GFF<-read.GFF("matt_data/ref_chandler_primary_omnic.gff")

# GENEs and CDS
GENE<-GFF[TYPE=="gene"]
GENE$Gene<-gsub("ID=gene-(.+);Dbxref.+", "\\1", GENE$INFO)

CDS<-GFF[TYPE=="CDS"]
CDS$Parent<-gsub(".+Parent=(.+);Dbxref.+", "\\1", CDS$INFO)
CDS$Product<-gsub(".+product=(.+);protein_id.+", "\\1", CDS$INFO)

CDS$Gene<-gsub(".+CDS;gene=(.+);product.+", "\\1", CDS$INFO)
CDS$ID<-1:nrow(CDS)

# extract protein AA sequences
proteins<-lapply(unique(CDS$Parent), function(p){
  CDS_p<-CDS[Parent==p]
  DIR<-unique(CDS_p$DIRECTION)
  CHROM=CDS_p$CHROM[1]
  genome_seq<-paste(sapply(1:nrow(CDS_p), function(r){
    START=CDS_p$START[r]
    STOP=CDS_p$STOP[r]
    paste(fasta[[CHROM]][START:STOP], collapse="")
  }), collapse = "")
  if(DIR=="-"){
    genome_seq<-revcomp(genome_seq)
  }
  genome_seq<-unlist(strsplit(genome_seq, split=""))
  AA_seq<-seqinr::translate(genome_seq)
})
names(proteins)<-unique(CDS$Parent)

# extract protein DNA sequences
cds<-lapply(unique(CDS$Parent), function(p){
  CDS_p<-CDS[Parent==p]
  DIR<-unique(CDS_p$DIRECTION)
  CHROM=CDS_p$CHROM[1]
  genome_seq<-paste(sapply(1:nrow(CDS_p), function(r){
    START=CDS_p$START[r]
    STOP=CDS_p$STOP[r]
    paste(fasta[[CHROM]][START:STOP], collapse="")
  }), collapse = "")
  if(DIR=="-"){
    genome_seq<-revcomp(genome_seq)
  }
  genome_seq<-unlist(strsplit(genome_seq, split=""))
  return(genome_seq)
})
names(cds)<-unique(CDS$Parent)

# check that protein sequence is proper ORF (M....*)
prot_check<-rbindlist(lapply(names(proteins), function(p){
  codon1<-proteins[[p]][1]
  codonend<-proteins[[p]][length(proteins[[p]])]
  return(data.table(p, codon1, codonend))
}), fill=T)
prot_check$CHECK<-prot_check$codon1=="M" & prot_check$codonend =="*"

# subset CDS and GENE to include good ORF proteins only
CDS_good<-CDS[Parent %in% prot_check[CHECK==T]$p]
GENE_good<-GENE[Gene %in% CDS_good$Gene]


#Intergenic regions
intergenic <- GENE_good[, .(
  SOURCE = "TAIR10",
  TYPE = "intergenic",
  START = STOP + 1L,
  STOP = shift(START, type = "lead") - 1L
), by = CHROM]
intergenic <- intergenic[!is.na(STOP)]
intergenic <- intergenic[START <= STOP]
intergenic$ID<-1:nrow(intergenic)

# update chromosome names
CDS_good$CHROM<-as.numeric(gsub("NC_0499(..).1_RagTag", "\\1", CDS_good$CHROM))
GENE_good$CHROM<-as.numeric(gsub("NC_0499(..).1_RagTag", "\\1", GENE_good$CHROM))
mapp$CHROM<-as.numeric(gsub("NC_0499(..).1_RagTag", "\\1", mapp$CHROM))
mapp<-mapp[!is.na(CHROM)]
names(fasta)<-as.numeric(gsub("NC_0499(..).1_RagTag", "\\1", names(fasta)))

#Define mappable IG and Genic
GENE_good$ID<-1:nrow(GENE_good)
GENE_good$LENGTH<-GENE_good$STOP-GENE_good$START
mapp_GENIC<-features_in_features(GENE_good,mapp[DEPTH==1], mode = "length")
# GENE_good$MAPP_LENGTH<-mapp_GENIC$length[match(GENE_good$ID, mapp_GENIC$ID)]
# GENE_good$MAPP_LENGTH2<-ifelse(GENE_good$MAPP_LENGTH>GENE_good$LENGTH,GENE_good$LENGTH,  GENE_good$LENGTH-GENE_good$MAPP_LENGTH)

intergenic$ID<-1:nrow(intergenic)
intergenic$LENGTH<-intergenic$STOP-intergenic$START
mapp_GENIC<-features_in_features(intergenic,mapp[DEPTH==1], mode = "length")
# intergenic$MAPP_LENGTH<-mapp_GENIC$length[match(intergenic$ID, mapp_GENIC$ID)]
# intergenic$MAPP_LENGTH2<-ifelse(intergenic$MAPP_LENGTH>intergenic$LENGTH,intergenic$LENGTH,  intergenic$LENGTH-intergenic$MAPP_LENGTH)

sum(GENE_good$LENGTH)/sum(intergenic$LENGTH)


# Identify HPs in the genome
# HPs<-make_homopolymer(fasta, 3)
# fwrite(HPs,"tables/ref_chandler_primary_omnic.fasta.HPs.csv")
# Or
HPs<-fread("tables/ref_chandler_primary_omnic.fasta.HPs.csv")

# Mutations ---------------------------------------------------------------

mutations<-fread("matt_data/all_dn_snps.vcf")
mutations$SAMPLE<-mutations$default

mutations$UNIQUE<-paste(mutations$CHROM, mutations$POS, mutations$ALT)
mutations[,N:=.N, by=UNIQUE]
SNPS<-nchar(mutations$REF)==1 & nchar(mutations$ALT)==1
mutations$SNP<-SNPS

# annotate CDS and GENIC
mutations$ID<-1:nrow(mutations)
mutations_CDS<-features_in_sites(CDS_good, mutations)
mutations$CDS<-mutations_CDS$overlaps[match(mutations$ID, mutations_CDS$ID)]

mutations_gene<-features_in_sites(GENE_good, mutations)
mutations$GENIC<-mutations_gene$overlaps[match(mutations$ID, mutations_gene$ID)]

#parse SAMPLE
mutations$VAF<-sapply(mutations$SAMPLE, function(x) as.numeric(vstrsplit(x, split=":")[5]))
mutations$DP<-sapply(mutations$SAMPLE, function(x) as.numeric(vstrsplit(x, split=":")[3]))
mutations$GQ<-sapply(mutations$SAMPLE, function(x) as.numeric(vstrsplit(x, split=":")[2]))
mutations$GT<-sapply(mutations$SAMPLE, function(x) (vstrsplit(x, split=":")[1]))
mutations$AD<-sapply(mutations$SAMPLE, function(x) (vstrsplit(x, split=":")[4]))
mutations$ALTD<-sapply(mutations$AD, function(x) as.numeric(vstrsplit(x, split=",")[2]))
mutations$REFD<-sapply(mutations$AD, function(x) as.numeric(vstrsplit(x, split=",")[1]))
# calcualate deviation from hterozygous expectation (REFD:ALTD == 1:1)
mutations$ChiP<-pbsapply(1:nrow(mutations), function(i) {
  ALTD<-mutations$ALTD[i]
  REFD<-mutations$REFD[i]
  chitest<-chisq.test(c(ALTD, REFD))
  return(chitest$p.value)
})
# N=number of Source with mutation


# DP and VAF based on informative reads
mutations$DP2<-mutations$ALTD+mutations$REFD
mutations$VAF2<-mutations$ALTD/mutations$DP2
mutations[,DP_norm:=log(DP/mean(DP)), by=Source]
mutations[,DP_norm2:=log(DP2/mean(DP2)), by=Source]

mutations$Source_type<-ifelse(mutations$Source %in% c("ref_chandler","tree2"), "tree",
                              ifelse(mutations$Source %in% c("cr10","cr13","cr22","cr85"), "shoot", "embryo"))
mutations[,meanVAF:=mean(VAF), by=UNIQUE]
mutations[,medDP:=median(DP), by=UNIQUE]
mutations[,medQUAL:=median(QUAL), by=UNIQUE]

mutations_means<-mutations[QUAL>10 & SNP==T,
                           .(CHROM=unique(CHROM),
                             POS=unique(POS),
                             REF=unique(REF),
                             ALT=unique(ALT),
                             CDS=unique(CDS),
                             GENIC=unique(GENIC)),
                           by=.(UNIQUE)]

# annotate MAPPABILITY == 1
mutations_means$ID<-1:nrow(mutations_means)
mutations_means$MAPP<-features_in_sites(mapp[DEPTH==1], mutations_means)$overlaps

# annotate effects of mutations
library(pbapply)
mutations_means$effect<-pbsapply(1:nrow(mutations_means), annot_mut_effect)
mutations_means$effect[mutations_means$effect=="nonCDS" & mutations_means$GENIC==T]<-"genic-nonCDS"

# annotate Homopolymer (HP) adjacent mutations
mutations_means$HP<-sapply(1:nrow(mutations_means), function(i){
  CHROM=mutations_means$CHROM[i]
  POS=mutations_means$POS[i]
  REF=mutations_means$REF[i]
  ALT=mutations_means$ALT[i]
  fasta_ref<-toupper(fasta[[CHROM]][(POS)])
  if(fasta_ref!=REF) stop("REF")

  up=toupper(fasta[[CHROM]][(POS-3):(POS-1)])
  down=toupper(fasta[[CHROM]][(POS+1):(POS+3)])
  all(ALT==up) | all(ALT==down)
})

# annotate trimer context of mutations
mutations_means$tricontext<-tricontexts(mutations_means, fasta)
mutations_means$mut<-paste0(substr(mutations_means$tricontext, 2, 2), substr(mutations_means$tricontext, 4, 4), substr(mutations_means$tricontext, 6, 6))


# add annotations to mutations dataset
mutations_archive<-mutations
mutations<-mutations[QUAL>10 & SNP==T]
mutations$effect<-mutations_means$effect[match(mutations$UNIQUE, mutations_means$UNIQUE)]
mutations$HP<-mutations_means$HP[match(mutations$UNIQUE, mutations_means$UNIQUE)]
mutations$MAPP<-mutations_means$MAPP[match(mutations$UNIQUE, mutations_means$UNIQUE)]
mutations$tricontext<-mutations_means$tricontext[match(mutations$UNIQUE, mutations_means$UNIQUE)]
mutations$mut<-mutations_means$mut[match(mutations$UNIQUE, mutations_means$UNIQUE)]

mutations$dist<-pbsapply(1:nrow(mutations), function(i){
  C=mutations$CHROM[i]
  P=mutations$POS[i]
  S=mutations$Source[i]
  other<-mutations[CHROM==C & Source==S & POS!=P]
  min(abs(other$POS-P))
})

# write tables ------------------------------------------------------------
fwrite(mutations, "tables/mutations.csv")
fwrite(mutations_means, "tables/mutations_means.csv")
fwrite(CDS_good, "tables/CDS_good.csv")
fwrite(GENE_good, "tables/GENE_good.csv")


# Analyses ----------------------------------------------------------------

table(mutations$effect)
mutations$VAF2
summary(lm(VAF2~effect+mut, mutations[effect %in% c("Non-Syn","Syn") & !grepl("cr85|cr22|cr10|cr13|ref|tree", Source)]))

hist(mutations[effect %in% c("Non-Syn","Syn")]$QUAL)
table(CHROMmutations$Source)

CDS_muts_embryos<-mutations[effect %in% c("Non-Syn","Syn") & !grepl("cr85|cr22|cr10|cr13|ref|tree", Source)]
ggplot(CDS_muts_embryos, aes(x=VAF, fill=effect))+
  geom_histogram(alpha=0.5)+
  facet_grid(CHROM~mut)

CDS_muts_embryos$nonsyn<-CDS_muts_embryos$effect=="Non-Syn"
plot_bins(CDS_muts_embryos, yvar = "nonsyn", xvar = "VAF", bins = 10, xaxis = "numeric")

CDS_muts_embryos$nonsyn<-CDS_muts_embryos$effect=="Non-Syn"
plot_bins(CDS_muts_embryos, yvar = "nonsyn", xvar = "N", bins = 10, xaxis = "numeric")

CDS_muts_embryos[,.(mean(nonsyn)), by=N]


View(CDS_muts_embryos[VAF==1])

mutations_means2<-mutations[!grepl("cr85|cr22|cr10|cr13|ref|tree", Source),
                           .(CHROM=unique(CHROM),
                             POS=unique(POS),
                             REF=unique(REF),
                             ALT=unique(ALT),
                             CDS=unique(CDS),
                             GENIC=unique(GENIC),
                             DP=mean(DP),
                             QUAL=mean(QUAL),
                             QUAL_med=median(QUAL),
                             effect=unique(effect),
                             VAF=mean(VAF),
                             N=mean(N),
                             GQ=mean(GQ),
                             ALTD=mean(ALTD),
                             #ChiP=mean(-log10(ChiP)),
                             Count=.N,
                             mut=unique(mut),
                             HOM=sum(GT=="1/1")/.N,
                             # Type=unique(Type),
                             DP_norm=mean(DP_norm),
                             DP_norm2=mean(DP_norm2),
                             DP2=mean(DP2),
                             VAF2=mean(VAF2)),
                           by=.(UNIQUE)]

summary(lm(QUAL~effect+mut, mutations_means2[effect %in% c("Non-Syn","Syn")]))

table(mutations_means2$effect)

CDS_muts_embryos_mean<-mutations_means2[effect %in% c("Non-Syn","Syn") ]

View(CDS_muts_embryos_mean)

ggplot(CDS_muts_embryos_mean, aes(x=VAF, fill=effect))+
  geom_histogram(alpha=0.5)+
  facet_grid(CHROM~mut, scales="free_y")

table(CDS_muts_embryos_mean$VAF>.5, CDS_muts_embryos_mean$effect)



gene_windows_mutations<-sites_in_features(gene_windows, sites = muts_embryos,mode = "mean", value="QUAL")
gene_windows$VAF<-gene_windows_mutations$mean[match(gene_windows$ID, gene_windows_mutations$ID)]
plot_feature_windows(gene_windows, variable="VAF", mode="mean")

ggplot(gene_windows, aes(x=REGION, y=VAF))+geom_boxplot()+
  geom_smooth()
