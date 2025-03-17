library(polymorphology2)
library(vcfR)
library(data.table)
library(stringr)

read_vcf_data <- function(file) {
  message(file)
  vcf_data <- read.vcfR(file)
  vcf_data <- data.table(cbind(vcf_data@fix, vcf_data@gt))
  vcf_data<-vcf_data[FILTER=="PASS" & !grepl(",", ALT) & grepl("hap",CHROM)]
  # Extract the sample name using gsub
  sample_name <- gsub(".*\\/|_combohap_deepvariant.vcf.*", "", file)
  vcf_data[, SOURCE := sample_name]
  
  colnames(vcf_data)[10]<-"default"
  # Create a unique identifier
  vcf_data[, UNIQUE := paste(CHROM, POS, REF, ALT, sep = "_")]
  vcf_data[, SNP:=nchar(REF)==1 & nchar(ALT)==1]
  # Extract FORMAT field names
  format_fields <- strsplit(vcf_data$FORMAT[1], ":")[[1]]
  
  # Split the 'default' column into multiple columns
  format_values <- tstrsplit(vcf_data[[10]], ":", type.convert = TRUE)
  # Assign column names based on FORMAT field names
  names(format_values) <- format_fields
  format_values<- as.data.table(format_values)
  # Bind parsed columns back to vcf_data
  vcf_data <- cbind(vcf_data, format_values)
  
  return(vcf_data)
}

# Apply function to all VCF files and combine them
vcf_files <- list.files("matt_data/Combined_Haplotype/", pattern = "*.vcf", full.names = TRUE)
vcf_data_combined <- rbindlist(lapply(vcf_files, read_vcf_data))
vcf_data_combined[, SNP:=nchar(REF)==1 & nchar(ALT)==1]
vcf_data_combined[, UNIQUEN:=.N, by = UNIQUE]
vcf_data_combined$QUAL<-as.numeric(vcf_data_combined$QUAL)

ggplot(vcf_data_combined[SNP==T], aes(x=SOURCE, y=VAF)) + geom_violin()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(UNIQUEN~.)

ggplot(vcf_data_combined[SNP==T], aes(x=QUAL, y=VAF)) + geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(UNIQUEN~SOURCE)

ggplot(vcf_data_combined[SNP==T,.(N=.N), by=.(UNIQUEN, SOURCE, FIXED=VAF>0.95)], aes(x=SOURCE, y=UNIQUEN, fill=N))+geom_tile()+theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(FIXED~.)

hist(vcf_data_combined$QUAL, breaks=100)

summary(vcf_data_combined$DP)
ggplot(vcf_data_combined[SNP==T & SOURCE=="11A" & UNIQUEN==1 & QUAL>10], aes(x=DP, y=VAF)) + geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(UNIQUEN~SOURCE)

ggplot(vcf_data_combined[SNP==T & SOURCE=="ref_chandler" & UNIQUEN<2 & QUAL>10 & DP>12], aes(x=VAF)) + 
  geom_histogram()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(vcf_data_combined[SNP==T& UNIQUEN<2 & QUAL>0 & DP>10], aes(x=VAF, alpha=cut2(QUAL, g=10))) + 
  geom_histogram()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 1))+
  facet_grid(SOURCE~., scale="free_y")+
  theme_classic()

ggplot(vcf_data_combined[SOURCE=="CR85" & SNP==T& UNIQUEN<2 & QUAL>0 & DP>10], aes(x=VAF, alpha=cut2(QUAL, g=10))) + 
  geom_histogram()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.y = element_text(angle = 0, hjust = 1))+
  facet_grid(SOURCE~., scale="free_y")