

# functions and packages --------------------------------------------------

library(polymorphology2)

sim_model<-function(diploid=10000, triploid=3000, haploid=3000){
  # Define the ploidy categories
  set.seed(123)
  ploidy <- c(rep("haploid", haploid), rep("diploid", diploid), rep("triploid", triploid))
  n_sites <- length(ploidy)
  # Initialize the data.table
  dt <- data.table(
    Ploidy = ploidy,
    SCALED_DP = numeric(n_sites),
    MAF = numeric(n_sites)
  )

  # Function to simulate Gamma distribution close to 0 or 1, with bounding

  gamma_close_to_0_1_bounded <- function(n, mode, shape = 1, rate = 100) {
    if (mode == 0) {
      return(pmin(1, pmax(0, rgamma(n, shape = shape, rate = rate)))) # Values close to 0
    } else {
      return(pmin(1, pmax(0, 1 - rgamma(n, shape = shape, rate = rate)))) # Values close to 1
    }
  }

  # Simulate SCALED_DP and MAF based on Ploidy
  # For haploid

  dt[Ploidy == "haploid", `:=`(
    SCALED_DP = pmax(0, rnorm(.N, mean = 0.5, sd = 0.1)),  # Bounded by 0
    MAF = sample(c(gamma_close_to_0_1_bounded(.N/2, 0), gamma_close_to_0_1_bounded(.N/2, 1)))
  )]

  # For diploid
  dt[Ploidy == "diploid", `:=`(
    SCALED_DP = pmax(0, rnorm(.N, mean = 1, sd = 0.1)),  # Bounded by 0
    MAF = ifelse(runif(.N) < 0.5,
                 gamma_close_to_0_1_bounded(.N, 0),
                 pmax(0, pmin(0.5, rnorm(.N, mean = 0.5, sd = 0.1))))
  )]

  # For triploid
  dt[Ploidy == "triploid", `:=`(
    SCALED_DP = pmax(0, rnorm(.N, mean = 1.5, sd = 0.1)),  # Bounded by 0
    MAF = ifelse(runif(.N) < 0.5,
                 gamma_close_to_0_1_bounded(.N, 0),
                 pmax(0, pmin(0.5, rnorm(.N, mean = 0.33, sd = 0.1))))
  )]
View(dt)
  # Check the first few rows of the data.table

  library(randomForest)
  rf_model<-randomForest(factor(Ploidy)~SCALED_DP+MAF, data=dt[sample(1:nrow(dt), 10000)])
  print(rf_model)
  return(rf_model)
}

read_test<-function(filepath, subset=NULL, type="all"){
  test<-fread(filepath, col.names = c("CHROM","POS","ALT","default"))
  if(!is.null(subset)) test<-test[sample(1:nrow(test), subset)]
  if(type=="VAR") test<-test[grepl("1/1|0/1", default)]

  # Add DP and MAF columns by parsing the default column
  test[, `:=`(
    # Extract the second value of "default" for DP (convert to numeric)
    DP = as.numeric(sapply(strsplit(default, ":"), function(x) x[3])),

    # MAF calculation based on rules
    MAF = ifelse(
      grepl(",", ALT),  # Check if ALT has a ","
      {  # If there is a comma
        # Extract the third value of "default", split it by "," and calculate MAF
        sapply(strsplit(sapply(strsplit(default, ":"), function(x) x[4]), ","),
               function(values) {
                 numeric_values <- as.numeric(values)
                 non_zero_values <- numeric_values[numeric_values != 0]
                 if (length(non_zero_values) > 0) {
                   min(non_zero_values) / sum(numeric_values)
                 } else {
                   1  # If all values are zero, set MAF to 1
                 }
               })
      },
      0  # If no comma in ALT, set MAF to 1
    )
  )]

  # View the updated data table
  test$SCALED_DP<-test$DP/mean(test$DP)

  test$CHROM<-as.numeric(gsub("NC_0499(..).1_RagTag","\\1", test$CHROM))
  test$MAF<-ifelse(grepl("0/0",test$default),0, test$MAF)
  test$source<-gsub("matt_data/ploidy_pred/(.+)_deepvariant_selected.gvcf.gz","\\1",filepath)
  message(gsub("matt_data/ploidy_pred/(.+)_deepvariant_selected.gvcf.gz","\\1",filepath))
  return(test)
}

ploidy_windows<-function(test, windows=50){

  test_sum<-test[,.(Ha=sum(Ploidy=="haploid")/.N,
                    Di=sum(Ploidy=="diploid")/.N,
                    Tr=sum(Ploidy=="triploid")/.N,
                    SCALED_DP=mean(SCALED_DP),
                    DP=mean(DP),
                    MAF=mean(MAF),
                    sumHOM=sum(MAF==1)/sum(MAF>0),
                    sumCALL=sum(MAF>0)), by=.(CHROM, POS=as.numeric(cut(POS, breaks=windows)))]

  #add column Ploidy which returns "haploid", "diploid", or "triploid", depending on which of columns Ha, Di, or Tr are greatest
  test_sum[, Ploidy := fifelse(Ha == pmax(Ha, Di, Tr), "haploid",
                               fifelse(Di == pmax(Ha, Di, Tr), "diploid", "triploid"))]
  #add column Prob which gives the Ha, Di, or Tr value, depending on Ploidy
  test_sum[, Prob := fifelse(Ploidy == "haploid", Ha,
                             fifelse(Ploidy == "diploid", Di, Tr))]


  return(test_sum)
}

theme_chrom<-function(){
  theme_void(base_size = 6)+
    theme(panel.border = element_rect(fill=NA), legend.key.size = unit(0.25, "cm"),
          panel.spacing = unit(0.1, "lines"),
          strip.text.y.left = element_text(hjust = 1))

}

#NEW:

predict_ploidy <- function(data) {
  # Define a helper function to calculate ploidy based on a value and given thresholds
  get_ploidy <- function(value, thresholds, ploidies) {
    # Find the closest threshold
    closest <- which.min(abs(value - thresholds))
    return(ploidies[closest])
  }

  # Ploidy thresholds
  dp_thresholds <- c(0.5, 1, 1.5)
  maf_thresholds <- c(1, 0.5, 0.3)
  ploidies <- c("haploid", "diploid", "triploid")

  # Calculate SCALED_DP_ploidy for each row
  data$SCALED_DP_ploidy <- pbsapply(data$SCALED_DP, get_ploidy, dp_thresholds, ploidies)

  # Initialize final ploidy prediction
  data$final_ploidy <- data$SCALED_DP_ploidy

  # If MAF > 0, calculate MAF_ploidy and adjust final ploidy prediction
  maf_indices <- data$MAF > 0
  data$MAF_ploidy[maf_indices] <- sapply(data$MAF[maf_indices], get_ploidy, maf_thresholds, ploidies)

  data$final_ploidy[maf_indices] <- ifelse(
    data$SCALED_DP_ploidy[maf_indices] == "haploid" & data$MAF_ploidy[maf_indices] == "haploid", "haploid",
    ifelse(data$SCALED_DP_ploidy[maf_indices] == "triploid" & data$MAF_ploidy[maf_indices] == "triploid", "triploid", "diploid")
  )

  return(data$final_ploidy)
}

# run ---------------------------------------------------------------------


files<-list.files("matt_data/ploidy_pred", full.names = T)

datasets_sub<-rbindlist(lapply(files, function(filepath){
  test_sub<-read_test(filepath, subset=1000000, type="all")
  return(test_sub)
}))

# rf_model<-sim_model(diploid=10000, triploid=6000, haploid=6000)
# # 1M sites takes about 1 min to run
# datasets_sub$Ploidy<-predict(rf_model, datasets_sub)
datasets_sub$Ploidy<-predict_ploidy(datasets_sub)


datasets_sub_ploidy<-rbindlist(lapply(unique(datasets_sub$source), function(s){
  test_sub<-datasets_sub[source==s]
  ploidy<-ploidy_windows(test_sub, windows=100)
  ploidy$source<-s
  return(ploidy)
}))

mapping <- c(
  "11A" = "cr2_11A",
  "12A" = "cr2_12A",
  "13A" = "cr2_13A",
  "15A" = "cr2_15A",
  "15A1" = "cr2_15A1",
  "16A" = "cr2_16A",
  "16A1" = "cr2_16A1",
  "16A2" = "cr2_16A2",
  "17A1" = "cr2_17A1",
  "17A2" = "cr2_17A2",
  "18A" = "cr2_18A",
  "18A1" = "cr2_18A1",
  "CR10" = "cr10",
  "CR13" = "cr13",
  "cr2" = "cr2_hifi",
  "CR21-1" = "cr21_1",
  "CR21-2" = "cr21_2",
  "CR22" = "cr22",
  "CR85" = "cr85",
  "ref_chandler" = "ref_chandler",
  "term_chandler" = "tree2"
)

# Apply the mapping to rename the `source` column
datasets_sub_ploidy[, source := mapping[source]]

# plot --------------------------------------------------------------------

pdf("figures/sandbox/Ideograms.pdf", width=5, height=2)
datasets_sub_ploidy$type<-ifelse(datasets_sub_ploidy$source %in% c("ref_chandler","term_chandler"),"Tree", ifelse(datasets_sub_ploidy$source %in% c("CR10","CR13","CR22","CR85"),"Shoot", "Embryo"))

###NOTE tip_order from make_mutation_tree
datasets_sub_ploidy$source<-factor(datasets_sub_ploidy$source, levels=tip_order)
ggplot(datasets_sub_ploidy, aes(x=POS, y=source2, fill=Ploidy, alpha=(Prob)))+
  geom_tile()+
  facet_grid(source2~CHROM, scales="free",space="free", switch = "y")+
  theme_chrom()+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  scale_fill_manual(values=c("gray","purple2","green4"))

ggplot(datasets_sub_ploidy, aes(x=POS, y=source, fill=sumHOM, alpha=log(sumCALL)))+
  geom_tile()+
  facet_grid(source~CHROM, scales="free",space="free", switch = "y")+
  theme_chrom()+
  scale_fill_gradientn(colors=c("gray","red"))+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))

ggplot(datasets_sub_ploidy, aes(x=POS, y=source, fill=SCALED_DP))+
  geom_tile()+
  facet_grid(source~CHROM, scales="free",space="free", switch = "y")+  theme_chrom()+
  scale_fill_gradient2(low="purple2", mid="white", high="green4",midpoint = 1, limits=c(0, 3))+
  theme(panel.border = element_rect(fill=NA))+
  scale_y_discrete(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))

dev.off()

