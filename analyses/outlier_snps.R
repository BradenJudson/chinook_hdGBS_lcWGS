setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(pcadapt); library(qvalue); library(cowplot)

# Function for converting VCF to bed (usable by pcadapt).
# Assumes VCF is gzipped.
vcf2pc <- \(vcf_file) system(paste0("plink.exe --vcf ", vcf_file, 
                           " --make-bed --aec --double-id --out ", 
                           gsub("\\.vcf.gz", "", vcf_file)))

# Make bed files.
vcf2pc("../data/vcfs/hdgbs_subset_134kSNPs_n362_imputed.vcf.gz")
vcf2pc("../data/vcfs/hdgbs_full_maf5_m15_original.vcf.gz")
vcf2pc("../data/vcfs/lcwgs_subset_134kSNPs_imputed.vcf.gz")

# Function for ocnducting pcadapt analyses.
pcadapt2 <- \(vcf, K, q.alpha) {
 
  # Read in bed file.
  input <- read.pcadapt(gsub("\\.vcf.gz", "\\.bed", vcf), type = "bed")
  # Conduct PCadapt analysis.
  x <- pcadapt(input = input, K = K)
  # Output screeplot for determining best K value.
  plot(x, option = 'screeplot')
  
  # Build and format output dataframe.
  snps <- read.delim(gsub("\\.vcf.gz", "\\.bim", vcf), 
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
hdgbs_imp <- pcadapt2("../data/vcfs/hdgbs_subset_134kSNPs_n362_original.vcf.gz", K = 5, q.alpha = 1/100)
hdgbs_ori <- pcadapt2("../data/vcfs/hdgbs_full_maf5_m15_original.vcf.gz", K = 5, q.alpha = 1/100)
lcwgs_imp <- pcadapt2("../data/vcfs/lcwgs_subset_134kSNPs_imputed.vcf.gz", K = 5, q.alpha = 1/100)

# For better plotting, read in chromosome info.
chrs <- read.delim("../data/otsh_sequence_report.tsv") %>% 
  filter(!Chromosome.name %in% c("MT", "Un")) %>% 
  mutate(Ots = gsub("LG", "Ots", Chromosome.name),
         LGN = as.factor(str_sub(Ots, 4, 5)))


Ots28 <- rbind(hdgbs_imp[hdgbs_imp$CHROM == "NC_056456.1",],
               hdgbs_ori[hdgbs_ori$CHROM == "NC_056456.1",],
               lcwgs_imp[lcwgs_imp$CHROM == "NC_056456.1",]) %>% 
  mutate(fillcol = case_when( # Awkward way to assign colours based on two simultaneous conditions.
    CHROM %in% unique(chrs$RefSeq.seq.accession)[seq(1, length(unique(chrs$RefSeq.seq.accession)), 2)] & out == "N" ~ "gray90",
    CHROM %in% unique(chrs$RefSeq.seq.accession)[seq(2, length(unique(chrs$RefSeq.seq.accession)), 2)] & out == "N" ~ "gray30",
    out == 'Y' ~ "red"
  ))

(chromplot <- ggplot(data = Ots28, 
                     aes(x = POS/1e6,
                         y = -log(pval),
                         colour = fillcol,
                         fill   = fillcol)) +
    geom_point(shape = 21) +
    theme_classic() +
    scale_fill_identity() +
    scale_colour_identity() +
    labs(x = "Position (Mbp)",
         y = expression(-log[10](p-value))) +
    theme(axis.title.y = element_text(size = 12),
          panel.grid   = element_blank(),
          legend.position  = "none",
          strip.background = element_blank(),
          strip.placement  = "inside",
          strip.text = element_text(size = 12)) +
    facet_wrap(~data, ncol = 1, scales = "free_x"))


ggsave("../plots/pcadapt_manhattan.tiff",
       dpi = 300, width = 10, height = 10)

