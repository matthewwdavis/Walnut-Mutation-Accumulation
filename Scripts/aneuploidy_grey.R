library(data.table)
library(seqinr)
library(ggplot2)
library(tidyverse)
library(polymorphology)
library(polymorphology2)

# Functions ---------------------------------------------------------------
parsevcf_indel<-function(file){
  
  vcf<-read.vcfR(file)
  fix<-data.table(vcf@fix)
  gt<-data.table(vcf@gt)
  all<-cbind(fix, gt)
  all$SNP<-nchar(all$REF)==1 & nchar(all$ALT)==1
  gtsplit<-all[SNP==F] %>% separate(col=UnnamedSample, into=c("GT","GQ","DP","AD","VAF","PL"), sep=":")
  gtsplit<-data.table(gtsplit)
  gtsplit$VAF<-as.numeric(gtsplit$VAF)
  gtsplit$DP<-as.numeric(gtsplit$DP)
  gtsplit$QUAL<-as.numeric(gtsplit$QUAL)
  gtsplit$unique<-paste(gtsplit$CHROM, gtsplit$POS, gtsplit$ALT)
  return(gtsplit)
}

parsevcf<-function(file){
  
  vcf<-read.vcfR(file)
  fix<-data.table(vcf@fix)
  gt<-data.table(vcf@gt)
  all<-cbind(fix, gt)
  all$SNP<-nchar(all$REF)==1 & nchar(all$ALT)==1
  gtsplit<-all[SNP==T] %>% separate(col=UnnamedSample, into=c("GT","GQ","DP","AD","VAF","PL"), sep=":")
  gtsplit<-data.table(gtsplit)
  gtsplit$VAF<-as.numeric(gtsplit$VAF)
  gtsplit$DP<-as.numeric(gtsplit$DP)
  gtsplit$POS<-as.numeric(gtsplit$POS)
  gtsplit$QUAL<-as.numeric(gtsplit$QUAL)
  gtsplit$unique<-paste(gtsplit$CHROM, gtsplit$POS, gtsplit$ALT)
  gtsplit$POS<-as.numeric(gtsplit$POS)
  gtsplit$chr<-as.numeric(gsub("NC_0499.*?([0-9]+).*", "\\1", gtsplit$CHROM))
  return(gtsplit)
}

plotquals<-function(gtsplit){
  sum<-gtsplit[!is.na(VAF),.(QUAL=mean(QUAL), DP=mean(DP)), by=cut(VAF, 20)]
  
  ggplot(sum, aes(x=cut, y=QUAL))+
    geom_point()+
    theme(axis.text.x = element_text(angle=45, hjust=1))
}


filtergt<-function(gt, Q, Depth, chitest){
  filt<-gt[!grepl(",", ALT) & FILTER=="PASS" & denovo==T & GT=="0/1" & DP>Depth]
  denovo<-filt %>% separate(col=AD, into=c("RefD","AltD"), sep=",")
  denovo<-data.table(denovo)
  denovo$AltD<-as.numeric(denovo$AltD)
  denovo$RefD<-as.numeric(denovo$RefD)
  
  denovo$VAF2<-denovo$AltD/(denovo$AltD+denovo$RefD)
  denovo_filtered<-denovo[QUAL>Q]
  denovo_filtered$chr<-as.character(gsub("chr(.+)_RagTag","\\1", denovo_filtered$CHROM))
  denovo_filtered$CHROM<-denovo_filtered$chr
  denovo_filtered$ID<-1:nrow(denovo_filtered)
  denovo_filtered$POSITION<-denovo_filtered$POS
  denovo_filtered$chi<-unlist(apply(denovo_filtered, 1, function(r){
    chisq.test(c(as.numeric(r["AltD"]), as.numeric(r["RefD"])), p=c(.5,.5))$p.value
  }))
  denovo_filtered<-denovo_filtered[chi>chitest]
  return(denovo_filtered)
}

# Read in data ------------------------------------------------------------

cr2<-parsevcf("./Data/VCFs/cr2_primary_deepvariant.vcf.gz")

ref<-parsevcf("./Data/VCFs/ref_chandler_primary_deepvariant.vcf.gz")

term<-parsevcf("./Data/VCFs/term_chandler_primary_deepvariant.vcf.gz")

genes<-read.GFF("Data/GFFs/ref_chandler_primary_default_scaffold_liftoffv1.gff")
genes$CHROM<-as.numeric(gsub("NC_0499.*?([0-9]+).*", "\\1", genes$CHROM))
genes<-genes[TYPE=="CDS"]
genes$ID<-1:nrow(genes)


# annotate "de novo" ------------------------------------------------------

# in this case, because it compares the entire VCF (PASS or otherwise) denovo is only TRUE if the variant was not found in the other files at all, meaning even if that site was still called "Ref" rather than "PASS" in the other datasets. By doing so, we eliminate variants that were called at sites potentially prone to bad calling. The filter step below then only considers variants with FILTER==PASS and denovo==T
cr2$denovo<-!cr2$unique %in% c(ref$unique, term$unique)
ref$denovo<-!ref$unique %in% c(cr2$unique, term$unique)
term$denovo<-!term$unique %in% c(ref$unique, cr2$unique)

# analyze -----------------------------------------------------------------

ref$shared<-ref$unique %in% c(cr2$unique) & ref$unique %in% term$unique
ref_shared<-ref[shared==T & QUAL>40 & DP>30 & DP <100 & FILTER=="PASS"]

term$shared<-term$unique %in% c(cr2$unique) & term$unique %in% term$unique
term_shared<-term[shared==T & QUAL>40 & DP>30 & DP <100 & FILTER=="PASS"]

cr2$shared<-cr2$unique %in% c(ref$unique) & cr2$unique %in% term$unique
cr2_shared<-cr2[shared==T & QUAL>40 & DP>30 & DP <100 & FILTER=="PASS"]


ref_shared$genome<-"ref"
term_shared$genome<-"term"
cr2_shared$genome<-"cr2"

cr2_shared$MAF<-ifelse(cr2_shared$VAF>0.5, 1-cr2_shared$VAF, cr2_shared$VAF)

ref_shared$MAF<-ifelse(ref_shared$VAF>0.5, 1-ref_shared$VAF, ref_shared$VAF)

term_shared$MAF<-ifelse(term_shared$VAF>0.5, 1-term_shared$VAF, term_shared$VAF)

find_closest <- function(values) {
  sapply(values, function(x) {
    # Define the target values
    targets <- c(0.5, .33, 0)
    
    # Calculate the absolute differences between the value and the targets
    differences <- abs(targets - x)
    
    # Find the target that has the minimum difference
    closest_target <- targets[which.min(differences)]
    
    return(closest_target)
  })
}

cr2_shared$DP_scaled<-cr2_shared$DP/mean(cr2_shared$DP)
ref_shared$DP_scaled<-ref_shared$DP/mean(ref_shared$DP)
term_shared$DP_scaled<-term_shared$DP/mean(term_shared$DP)
shared_wide<-merge(merge(cr2_shared, ref_shared, by = "unique"), term_shared, by = "unique")
shared_wide$chi<-apply(shared_wide, 1, function(x){
  
  xsplit<-as.numeric(unlist(strsplit(x["AD.x"], split=",")))
  ysplit<-as.numeric(unlist(strsplit(x["AD.y"], split=",")))
  chi<-chisq.test(xsplit, p=prop.table(ysplit))
  chi$p.value
})

#pdf("~/Desktop/chitests.pdf", width=7, height=2)

shared_wide$closest<-find_closest(shared_wide$MAF.x)

ggplot(shared_wide, aes(x=POS.x, y=-log10(p.adjust(chi)), col=factor(closest)))+
  geom_point(size=0.1, alpha=.5)+
  facet_grid(~chr.x, space="free", scales = "free_x")+
  scale_color_manual(values=c("lightgreen","darkblue", "gray90"), name="MAF", labels=c("LOH","Triploid","Diploid"))+
  theme_classic(base_size = 6)+
  scale_x_continuous(breaks=(1:6)*10000000, labels=1:6, name="10e7 bp")


shared_long<-rbind(rbind(cr2_shared, ref_shared), term_shared)

ggplot(shared_long[!is.na(chr)], aes(x=factor(chr), y=DP_scaled, fill=(genome)))+
  geom_boxplot(outlier.size = 0, outlier.color = NA)+
  theme_classic()+
  theme(strip.background = element_blank(), legend.position = "bottom",
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_grid(~chr, scales = "free_x", space = "free")+
  scale_x_discrete("Chromosome")+
  scale_y_continuous(name="Mean scaled depth", limits=c(min(shared_long$DP_scaled), 2))+
  scale_fill_brewer(palette = "Dark2")
#ggsave("~/Desktop/depth_boxplot_cr2_long.pdf", height = 3, width = 10)
#ggsave("~/Desktop/depth_boxplot_cr2_long.png", height = 3, width = 8)

#dev.off()