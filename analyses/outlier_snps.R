setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(pcadapt); library(qvalue)

# Function for converting VCF to bed (usable by pcadapt).
vcf2pc <- \(vcf_file) system(paste0("plink.exe --vcf ", vcf_file, 
                           " --make-bed --aec --out ", 
                           gsub("\\.vcf", "", vcf_file)))

# Make bed files.
vcf2pc("../data/vcfs/hdgbs_subset_134kSNPs_n362_imputed.vcf")


# Function for ocnducting pcadapt analyses.
pcadapt2 <- \(vcf, K, q.alpha) {
 
  # Read in bed file.
  input <- read.pcadapt(gsub("\\.vcf", "\\.bed", vcf), type = "bed")
  # Conduct PCadapt analysis.
  x <- pcadapt(input = input, K = K)
  # Output screeplot for determining best K value.
  plot(x, option = 'screeplot')
  
  # Build and format output dataframe.
  snps <- read.delim(gsub("\\.vcf", "\\.bim", vcf), 
                     header = FALSE)[,c(1,4)] %>% 
    `colnames<-`(., c("CHROM", "POS")) %>% 
    mutate(pval = x$pvalues,
           chi2 = x$chi2.stat,
           data = basename(gsub("\\.vcf", "", vcf)),
           # Adjust p-values for multiple comparisons.
           qval = qvalue(pval)$qvalues,
           out  = as.factor(case_when(
             qval >= q.alpha ~ 'N', # Non-outlier.
             qval <  q.alpha ~ 'Y'  # Outlier.
           )))
  print(summary(snps$out)); return(snps)
  
}

# Conduct pcadapt with specified parameters. Returns dataframe of all SNPs and their outlier scores/status.
hdgbs_imp <- pcadapt2(vcf = "../data/vcfs/hdgbs_subset_134kSNPs_n362_original.vcf", K = 5, q.alpha = 0.001)

# For better plotting, read in chromosome info.
chrs <- read.delim("../data/otsh_sequence_report.tsv") %>% 
  filter(!Chromosome.name %in% c("MT", "Un")) %>% 
  mutate(Ots = gsub("LG", "Ots", Chromosome.name),
         LGN = as.factor(str_sub(Ots, 4, 5)))

# "Custom" ggplot-based Manhattan plot function.
ggManhattan <- \(df) {
  
  # Add 'Ots##' style chromosome names.
  df <- merge(hdgbs_imp, 
              chrs[,c("RefSeq.seq.accession","Chromosome.name")], 
              by.x = "CHROM", by.y = "RefSeq.seq.accession") %>% 
    mutate(fillcol = case_when( # Awkward way to assign colours based on two simultaneous conditions.
      CHROM %in% unique(hdgbs_imp$CHROM)[seq(1, length(unique(hdgbs_imp$CHROM)), 2)] & out == "N" ~ "gray90",
      CHROM %in% unique(hdgbs_imp$CHROM)[seq(2, length(unique(hdgbs_imp$CHROM)), 2)] & out == "N" ~ "gray30",
      out == 'Y' ~ "red"
    ))
  
  ggplot() +
    geom_point(data = df,
               aes(x = POS,
                   y = -log10(pval),
                   colour = fillcol, 
                   fill = fillcol),
               shape = 21) +
    scale_fill_identity() +
    scale_colour_identity() +
    facet_grid(cols = vars(Chromosome.name),
               scales = "free_x",
               switch = 'x',
               space = "free_x") +
    theme_classic() +
    labs(x = NULL,
         y = expression(-log[10](p-value))) +
    # scale_colour_manual(values = rep(c("gray80", "gray30"), 17)) +
    # scale_fill_manual(values = rep(c("gray80", "gray30"),   17)) +
    theme(axis.text.x = element_blank(),
          axis.title.y = element_text(size = 16),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.spacing = unit(0.1, "cm"), 
          strip.background = element_blank(),
          strip.text.x = element_text(size = 8),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          strip.placement = "outside")

}

ggManhattan(hdgbs_imp)
