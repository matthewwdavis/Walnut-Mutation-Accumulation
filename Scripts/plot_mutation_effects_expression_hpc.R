# Setup
## Libraries
library(data.table)
library(tidyverse)

## Read in data
df <- fread("../mutation_transcript_association.tsv")

# Plotting
## Creating significance lines
alpha <- 0.01
tests <- length(1:nrow(df))

p_threshold <- -log10(alpha)

bon_threshold <- -log10(alpha / tests)

## Create color map
unique_chroms <- unique(sort(df$CHROM))
colors <- rep(c("skyblue4", "goldenrod"), length.out = length(unique_chroms))

color_map <- setNames(colors, unique_chroms)

## Plot
df %>%
  ggplot(aes(x=START/1e7, y=-log10(pvalue), color = factor(CHROM)))+
  geom_point(size = .5) +
  geom_hline(yintercept = p_threshold, linetype = "dotted", color = "grey30", linewidth = 0.5)+
  geom_hline(yintercept = bon_threshold, linetype = "dotted", color = "grey30", linewidth = 0.5)+
  geom_text(data = df %>% filter(-log10(pvalue) > bon_threshold), 
            aes(label = rna), 
            vjust = 0.75, hjust = -.1, size = 1) +
  facet_grid(~CHROM, space = "free", scales = "free_x") +
  scale_color_manual(values = color_map) +
  labs(x = "Position (10 Mb)", y = "- Log 10 P Value") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./mut_expression_effect_plot.pdf", height = 1, width = 4)
ggsave("./mut_expression_effect_plot.png", height = 1, width = 4)
