library(ape)
library(phangorn)
library(polymorphology2)


source_type_colors<-c("embryo"="gold2", "shoot"="darkseagreen2","tree"="darkseagreen4")

# Define the function
gene_body_source_type <- function(source_type) {

  # Filter mutations for the given source type
  muts_filtered <- unique(mutations_filt[Source_type == source_type, .(CHROM, POS, effect)])
  muts_filtered$CHROM <- as.numeric(as.character(muts_filtered$CHROM))
  gene_windows$CHROM <- as.numeric(as.character(gene_windows$CHROM))

  # Calculate mutations in gene windows
  gene_windows_mutations <- sites_in_features(gene_windows, sites = muts_filtered, mode = "counts")
  gene_windows$muts <- gene_windows_mutations$counts[match(gene_windows$ID, gene_windows_mutations$ID)]

  # Plot and summarize the feature windows
  summary <- plot_feature_windows(gene_windows, variable = "muts", mode = "percent")$summary
  maxpos <- max(summary$RELATIVEPOS)

  # Chi-squared test
  chisq.test(summary[, 4:5])

  # Generate the plot
  plot <- ggplot(summary, aes(x = RELATIVEPOS, y = y, fill = REGION == "gene body")) +
    geom_bar(stat = "identity", col = "black", width = 0.75) +
    geom_vline(xintercept = c(maxpos / 3, maxpos / 3 * 2) + 0.5, linetype = "dashed") +
    theme_classic(base_size = 6) +
    scale_x_continuous(breaks = c(1, (maxpos / 3) + 0.5, (maxpos / 3 * 2) + 0.5, maxpos),
                       labels = c("-5kb", "START", "STOP", "+5kb")) +
    scale_fill_manual(values = c("gray", "green4"), guide="none") +
    ggtitle(paste("Ns =", sum(muts_filtered$effect == "Non-Syn"), "S =", sum(muts_filtered$effect == "Syn")))

  return(plot)
}

fasta_trimer_freq<-function(fasta, stop = NULL){
  fasta_windows<-rbindlist(lapply(1:length(fasta), function(i){
    data.table(CHROM=names(fasta)[i], START=1, STOP=length(fasta[[i]]), ID=i)
    data.table(CHROM=names(fasta)[i], START=1, STOP=ifelse(is.null(stop), length(fasta[[i]]), stop), ID=i)

  }))[!is.na(CHROM)]
  trimer_freq<-Nmerfrequency(features = fasta_windows, fasta=fasta, Nmer = 3, mode="counts")
  trimer_freq<-rbindlist(lapply(1:ncol(trimer_freq), function(i){
    data.table(TRI=names(trimer_freq)[i], N=sum(trimer_freq[[i]]))
  }))

  trimer_freq<-trimer_freq[!grepl("ID|N|V", TRI)]
}

make_tree<-function(mutations){
  mutations_matrix<-mutations
  mutations_matrix<-data.table(table(UNIQUE=mutations_matrix$UNIQUE, Source=mutations_matrix$Source))
  mutations_matrix<-dcast(mutations_matrix, Source~UNIQUE)
  mutations_matrix<-as.matrix(mutations_matrix)
  row.names(mutations_matrix)<-mutations_matrix[,1]
  mutations_matrix<-mutations_matrix[,-1]

  library(ape)
  library(phangorn)

  # Assume mutations_matrix is already defined (row names are Samples, colnames are UNIQUE)
  # Rows: samples, Columns: mutations

  # Step 1: Convert the binary mutation matrix into a phyDat object
  phy_data <- phyDat(mutations_matrix, type = "USER", levels = c(0, 1))

  # Step 2: Create a distance matrix using binary distances
  dist_matrix <- dist(mutations_matrix, method = "binary")

  # Step 3: Generate a UPGMA tree from the distance matrix
  tree <- upgma(dist_matrix)
  #plot(tree)

  # Step 4: Calculate branch lengths based on unique mutations
  # Initialize branch lengths
  branch_lengths <- numeric(nrow(tree$edge))

  branch_muts<-c()
  # Traverse each edge and calculate the unique mutations
  for (i in seq_len(nrow(tree$edge))) {
    message(i)
    node1 <- tree$edge[i, 1]
    node2 <- tree$edge[i, 2]

    # Get the tips (samples) downstream of node2
    descendants <- Descendants(tree, node2, type = "tips")[[1]]


    # Get the tip labels corresponding to these descendants
    tips_in_branch <- tree$tip.label[descendants]
    message(paste0(rownames(mutations_matrix)[descendants], collapse=" "))
    # Identify the columns in mutation_matrix that correspond to these tips
    tip_indices <- match(tips_in_branch, row.names(mutations_matrix))

    # Subset the mutation matrix for these tips
    relevant_mutations <- (mutations_matrix[tip_indices,])
    if(length(tip_indices)==1){
      relevant_mutations<-matrix(relevant_mutations, nrow=1)
    }
    relevant_mutations <- apply(relevant_mutations, 2, as.numeric)
    if(length(tip_indices)==1){
      relevant_mutations<-matrix(relevant_mutations, nrow=1)
    }
    out_relevant_mutations <- (mutations_matrix[-tip_indices,])
    if(length(tip_indices)==length(tree$tip.label)-1){
      out_relevant_mutations<-matrix(out_relevant_mutations, nrow=1)
    }
    out_relevant_mutations <- apply(out_relevant_mutations, 2, as.numeric)
    if(length(tip_indices)==length(tree$tip.label)-1){
      out_relevant_mutations<-matrix(out_relevant_mutations, nrow=1)
    }

    # Count mutations that are unique to this branch (present in at least one tip downstream, absent in others)
    group_mutations <- apply(relevant_mutations, 2, function(x) all(x==1))

    non_outgroup_mutations <- apply(out_relevant_mutations, 2, function(x) !any(x==1))
    group_mutations_names<-colnames(mutations_matrix)[which(group_mutations)]

    group_mutations_names<-colnames(mutations_matrix)[which(group_mutations & non_outgroup_mutations)]

    branch_muts<-c(branch_muts, group_mutations_names)
    unique_mutations<-sum(group_mutations & non_outgroup_mutations)
    # Set the branch length to the count of unique mutations
    branch_lengths[i] <- sum(unique_mutations)
  }
  branch_muts<-unique(branch_muts)
  # Step 5: Assign calculated branch lengths to the tree
  tree$edge.length <- branch_lengths
  tree$tip.color<-ifelse(tree$tip.label %in% c("ref_chandler","tree2"), "tree",
                         ifelse(tree$tip.label %in% c("cr10","cr13","cr22","cr85"), "shoot", "embryo"))

  # Ensure the tip labels and colors are stored correctly in the tree object
  tree$tip.label <- as.character(tree$tip.label)

  # Create a data frame from the tree object for easier mapping of colors
  tip_data <- data.frame(label = tree$tip.label, color = tree$tip.color)

  library(data.table)
  # Step 1: Create a ggtree object with tip labels as colored boxes
  p <- ggtree(tree, branch.length = "edge.length", size=0.25) %<+% tip_data +  # Use the `%<+%` operator to attach data
    geom_tiplab(aes(label = label, fill = color),
                geom = "label",
                size = 2,
                align = FALSE,
                offset = 20, linesize = 0.25) +  # Add tip labels as boxes
    #geom_text2(aes(label = branch.length), vjust = 1.5, hjust = 1.5, size = 2, col = "red") +  # Add branch length labels
    scale_fill_manual(values = source_type_colors, name=NULL, labels=c("Embryo","Shoot","Tree")) +  # Define fill colors
    theme_tree2() +
    theme(text = element_text(size = 6), legend.position = c(0.2, 0.8), legend.key.size = unit(0.25,"cm"))

  # Extract the maximum x position of the tips
  max_x <- max(p$data$x)

  # Expand the x-axis to fit all points
  p <- p + expand_limits(x = max_x * 1.2)

  # Step 2: Display the plot
  print(p)

  return(list(tree=tree, plot=p, branch_muts=branch_muts))
}

mutations_filt<-mutations[dist>300 & medDP>10 & MAPP==T]

pdf("figures/sandbox/mutations_tree_InDels.pdf", width=3.5, height=3)
tree<-make_tree(mutations_filt[])
tip_order <- get_taxa_name(tree$plot)
mutations_filt$branch<-mutations_filt$UNIQUE %in% tree$branch_muts
dev.off()

# Other analyses ----------------------------------------------------------

#trimer_freq<-fasta_trimer_freq(fasta, stop=1000000)
pdf("figures/sandbox/trimers_mutations.pdf", width=3.5, height=1)
plot_tricontexts(unique(mutations_filt[Source_type=="tree",.(UNIQUE, tricontext)])$tricontext, full=T, trimer_freq = trimer_freq, x_text = F)
plot_tricontexts(unique(mutations_filt[Source_type=="embryo",.(UNIQUE, tricontext)])$tricontext, full=T, trimer_freq = trimer_freq, x_text = F)
plot_tricontexts(unique(mutations_filt[Source_type=="shoot",.(UNIQUE, tricontext)])$tricontext, full=T, trimer_freq = trimer_freq, x_text = F)
dev.off()

gene_windows<-feature_windows(GENE_good, breaks=3,dist=5000, directed=T, IDcol="Gene")
gene_windows$CHROM<-as.numeric(as.character(gene_windows$CHROM))
mapp$CHROM<-as.numeric(as.character(mapp$CHROM))
gene_windows_mapp<-features_in_features(gene_windows, mapp, mode = "mean", value = "DEPTH")
gene_windows<-gene_windows[ID %in% gene_windows_mapp[mean==1]$ID]

pdf("figures/sandbox/gene_bodies.pdf", width=2, height=1)
gene_body_source_type("tree")
gene_body_source_type("shoot")
gene_body_source_type("embryo")
dev.off()

IG_G<-dcast(data.table(table(LOC=ifelse(mutations_filt$GENIC,"GENIC","INTERGENIC"), Source=mutations_filt$Source)), Source~LOC)
IG_G$Source_type<-ifelse(IG_G$Source %in% c("ref_chandler","tree2"), "tree",
                              ifelse(IG_G$Source %in% c("cr10","cr13","cr22","cr85"), "shoot", "embryo"))


genome_ratio<-sum(GENE_good$LENGTH)/sum(intergenic$LENGTH)
IG_G$chi<-sapply(1:nrow(IG_G), function(i){

 chitest<- chisq.test(c(IG_G$GENIC[i], IG_G$INTERGENIC[i]), p=prop.table(c(sum(GENE_good$LENGTH),sum(intergenic$LENGTH))))
 chitest$p.value
})

IG_G$Source<-factor(IG_G$Source, levels=rev(tip_order))
pdf("figures/sandbox/IG_G_SNPS.pdf", width=2, height=2)

ggplot(IG_G, aes(x=(GENIC/INTERGENIC), y=Source, fill=Source_type, size=-log10(chi)))+
  geom_segment(aes(x=(GENIC/INTERGENIC), xend=genome_ratio, yend=Source), size=0.25, linetype="dashed")+
  geom_point(shape=21)+
  geom_vline(xintercept = genome_ratio)+
  scale_size_continuous(range = c(1,2), guide="none")+
  my_minimal()+
  scale_x_continuous(limits=c(0,0.7))+
  theme(axis.line.x = element_line())+
  geom_text(aes(x=(GENIC/INTERGENIC)-0.05, label=ifelse(chi<2e-10,"**","*")), size=3)+
  scale_fill_manual(values=source_type_colors, guide="none")
dev.off()

# InDels ------------------------------------------------------------------

# ### ARCHIVE
# mutations_filt<-mutations[!grepl(",", ALT) & SNP==F & medQUAL>30 &QUAL>30 & medDP>20 & ALTD>5]
# mutations_filt[,N2:=.N, by=UNIQUE]
#
# mutations_filt<-mutations_filt[N2==N]
#
# pdf("figures/sandbox/mutations_tree_InDels.pdf", width=3.5, height=3)
# tree<-make_tree(mutations_filt[N2==N])
# tip_order <- get_taxa_name(tree$plot)
# mutations_filt$branch<-mutations_filt$UNIQUE %in% tree$branch_muts
# dev.off()
#
# mutations_filt$size<-nchar(mutations_filt$ALT)-nchar(mutations_filt$REF)
# mutations_filt_size_counts<-unique(mutations_filt[branch==T & N2==N & abs(size<50),.(CHROM, POS, REF, ALT, size, Source_type)])
# INS<-mutations_filt_size_counts[size==1]
# DEL<-mutations_filt_size_counts[size==-1]
# chisq.test(table(substr(DEL$REF, 2, 2), DEL$Source_type))
# chisq.test(table(substr(INS$ALT, 2, 2), INS$Source_type))
#
# table(substr(INS$ALT, 2, 2), INS$Source_type)
# table(substr(DEL$REF, 2, 2), DEL$Source_type)
#
# ggplot(mutations_filt_size_counts, aes(x=size, fill=Source_type))+
#   geom_histogram(bins=200)+
#   scale_y_log10()+
#   facet_grid(Source_type~.)+geom_vline(xintercept = 0)
#
# mutations_filt[,N2:=.N, by=UNIQUE]
#
