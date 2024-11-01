library(polymorphology2)

# List of expression files
exp_files <- list.files("matt_data/expression_data/", full.names = TRUE)

# Read and combine all expression files
all_exp_data <- rbindlist(lapply(exp_files, function(file) {
  # Read each file
  exp_data <- fread(file)
  # Extract the source from the file name (assuming the format is consistent)
  exp_data$Source <- sub("_abundance.tsv$", "", basename(file))
  # Return the data frame
  return(exp_data)
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
all_exp_data[, Source := mapping[Source]]

# add columns
all_exp_data$Gene<-gsub(".+gene=gene-(.+)name=.+","\\1", all_exp_data$target_id)
all_exp_data$RNA<-gsub("(.+)gene=gene-.+","\\1", all_exp_data$target_id)
all_exp_data$CHROM<-as.numeric(gsub(".+=NC_0499(..).+","\\1", all_exp_data$target_id))
all_exp_data[,meantpm:=mean(tpm), by=RNA]
all_exp_data[,mintpm:=min(tpm), by=RNA]

all_exp_data[RNA==R]
tree_exp<-drop.tip(tree$tree, tip =tree$tree$tip.label[!tree$tree$tip.label %in% all_exp_data$Source] )

plot_exp_branch<-function(i, R){

  all_exp_data_tmp<-all_exp_data[mintpm!=0]

  node1 <- tree_exp$edge[i, 1]; node2 <- tree_exp$edge[i, 2]

  # Get the tips (samples) downstream of node2
  descendants <- Descendants(tree_exp, node2, type = "tips")[[1]]
  tips_in_branch <- tree_exp$tip.label[descendants]

  all_exp_data_tmp[, tip := Source %in% tips_in_branch]  # Add the tip column once
  all_exp_data_RNA <- all_exp_data_tmp[RNA == R]

  # Perform t-test based on random or not
  ttest <- if (random) {
    t.test((all_exp_data_RNA$tpm) ~ sample(all_exp_data_RNA$tip))
  } else {
    t.test(all_exp_data_RNA$tpm ~ all_exp_data_RNA$tip)
  }



  wilcoxtest <- if (random) {
    wilcox.test((all_exp_data_RNA$tpm) ~ sample(all_exp_data_RNA$tip))
  } else {
    wilcox.test(all_exp_data_RNA$tpm ~ all_exp_data_RNA$tip)
  }

  plot<-ggplot(all_exp_data_RNA, aes(x=tip, y=tpm))+
    geom_boxplot()

  tests<-data.table(
    RNA = R,
    Gene = unique(all_exp_data_RNA$Gene),
    p = ttest$p.value,
    stat = ttest$statistic,
    wilcoxp=wilcoxtest$p.value,
    w=wilcoxtest$statistic,
    tips = paste(tips_in_branch, collapse = " "),
    branch = i
  )
  return(list(plot, tests))
}


exp_tree_tests <- function(tree_exp, all_exp_data, Gene_sub = NULL, random = FALSE) {

  # Create a copy and subset once
  all_exp_data_tmp <- if (!is.null(Gene_sub)) {
    all_exp_data[Gene %in% Gene_sub]
  } else {
    all_exp_data
  }
  all_exp_data_tmp<-all_exp_data[mintpm!=0]
  unique_RNA <- unique(all_exp_data_tmp$RNA)  # Extract unique RNA once

  exp_tests <- rbindlist(lapply(seq_len(nrow(tree_exp$edge)), function(i) {
    node1 <- tree_exp$edge[i, 1]; node2 <- tree_exp$edge[i, 2]

    # Get the tips (samples) downstream of node2
    descendants <- Descendants(tree_exp, node2, type = "tips")[[1]]
    tips_in_branch <- tree_exp$tip.label[descendants]

    if (length(tips_in_branch) > 1) {
      all_exp_data_tmp[, tip := Source %in% tips_in_branch]  # Add the tip column once

      # Use pblapply to iterate over unique RNA values
      ttests <- rbindlist(pblapply(unique_RNA, function(R) {
        all_exp_data_RNA <- all_exp_data_tmp[RNA == R]

        # Perform t-test based on random or not
        ttest <- if (random) {
          t.test((all_exp_data_RNA$tpm) ~ sample(all_exp_data_RNA$tip))
        } else {
          t.test(all_exp_data_RNA$tpm ~ all_exp_data_RNA$tip)
        }

        wilcoxtest <- if (random) {
          wilcox.test((all_exp_data_RNA$tpm) ~ sample(all_exp_data_RNA$tip))
        } else {
          wilcox.test(all_exp_data_RNA$tpm ~ all_exp_data_RNA$tip)
        }

        # Return a data.table for each RNA test
        data.table(
          RNA = R,
          Gene = unique(all_exp_data_RNA$Gene),
          p = ttest$p.value,
          stat = ttest$statistic,
          wilcoxp=wilcoxtest$p.value,
          w=wilcoxtest$statistic,
          tips = paste(tips_in_branch, collapse = " "),
          branch = i
        )
      }))

      return(ttests)
    } else {
      return(NULL)
    }
  }))

  return(exp_tests)
}


make_exp_windows<-function(test, windows=50){

  expression_windows<-real_tests2[,.(stat=median(stat, na.rm=T),
                    sig=median(-log10(p), na.rm=T)
                   ), by=.(CHROM, POS=as.numeric(cut(START, breaks=windows)), branch)]

  ggplot(expression_windows, aes(x=POS, y=branch, fill=stat, alpha=sig))+
    geom_tile()+
    facet_grid(branch~CHROM, scales="free",space="free", switch = "y")+
    theme_chrom()+
    scale_fill_gradient2()+
    scale_y_discrete(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))
}


# Analyses: ---------------------------------------------------------------

## takes 1hour:
real_tests<-exp_tree_tests(tree_exp, all_exp_data, Gene_sub=NULL, random=F)
real_tests<-merge(real_tests,GENE_good, by=c("Gene"))
mintpm<-all_exp_data[,.(RNA, mintpm)]
real_tests$mintpm<-mintpm$mintpm[match(real_tests$RNA, mintpm$RNA)]
real_tests<-real_tests[mintpm>0]
real_tests[,padjust:=p.adjust(p, method="fdr"), by=branch]


real_tests$sig<- -log10(real_tests$p)
real_tests$sig<-ifelse(real_tests$stat<0, real_tests$sig, -real_tests$sig)
ggplot(real_tests[padjust<0.05], aes(x=START, col=stat>0, y=sig))+
  geom_point(size=0.5)+
  facet_grid(branch~CHROM, space="free_x", scales="free_x")+
  scale_color_manual(values=c("blue","red"))

sigcounts<-data.table(table(branch=real_tests[padjust<0.05]$branch, tips=real_tests[padjust<0.05]$tips))[N>0]
sigcounts

plot_exp_branch(i=15, R="rna-XM_018973053.2")

