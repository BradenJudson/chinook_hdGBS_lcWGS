setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(pcadapt); library(qvalue); 
library(cowplot); library(bigutilsr); library(ggmagnify)

# Function for converting VCF to bed (usable by pcadapt).
# Assumes VCF is gzipped.
vcf2pc <- \(vcf_file) system(paste("plink.exe --vcf", vcf_file, 
                           " --make-bed --aec --double-id --out", 
                           gsub("\\.vcf.gz", "", vcf_file)))

# Make bed files.
# vcf2pc("../data/vcfs/hdgbs_subset_134kSNPs_n362_imputed.vcf.gz")
vcf2pc("../data/vcfs/hdgbs_n361_maf5_m70.vcf.gz")
vcf2pc("../data/vcfs/chinook_lcwgs_maf1_m15_n361_imputed.vcf.gz")
vcf2pc("../data/vcfs/lcWGS_full_8MSNPs_imputed.vcf.gz")

# Function for conducting pcadapt analyses.
pcadapt2 <- \(vcf, K, q.alpha) {
 
  # Read in bed file.
  input <- read.pcadapt(gsub("\\.vcf.gz", "\\.bed", vcf), type = "bed")
  # Conduct PCadapt analysis.
  x <- pcadapt(input = input, K = K, min.maf = 0.01)
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
hdgbs_ori <- pcadapt2("../data/vcfs/hdgbs_n361_maf5_m70.vcf.gz", K = 5, q.alpha = 1/100)
# lcwgs_imp <- pcadapt2("../data/vcfs/lcWGS_full_8MSNPs_imputed.vcf.gz", K = 5, q.alpha = 1/100)
lcwgs_imp2 <- pcadapt2("../data/vcfs/chinook_lcwgs_maf1_m15_n361_imputed.vcf.gz", K = 5, q.alpha = 1/100)

ggplot(data = lcwgs_imp2[lcwgs_imp2$CHROM == "NC_056456.1",], aes(x = POS, y = -log10(pval))) + geom_point() + facet_wrap(~CHROM)
ggplot(data = hdgbs_ori[hdgbs_ori$CHROM == "NC_056456.1",], aes(x = POS, y = -log10(pval))) + geom_point() + facet_wrap(~CHROM)

# For better plotting, read in chromosome info.
chrs <- read.delim("../data/otsh_sequence_report.tsv") %>% 
  filter(!Chromosome.name %in% c("MT", "Un")) %>% 
  mutate(Ots = gsub("LG", "Ots", Chromosome.name),
         LGN = as.factor(str_sub(Ots, 4, 5)))


# Isolate Ots28 for current purposes.
Ots28 <- rbind(hdgbs_imp[hdgbs_imp$CHROM == "NC_056456.1",],
               hdgbs_ori[hdgbs_ori$CHROM == "NC_056456.1",],
               lcwgs_imp[lcwgs_imp$CHROM == "NC_056456.1",]) 

(chromplot <- ggplot(data = Ots28, 
                     aes(x = POS/1e6,
                         y = -log(pval),
                         colour = out)) +
    scale_color_manual(values = c("gray70", "red1")) +
    geom_point() +
    theme_classic() +
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


# ------------------------------------------------------------------------------

# p-value calculations from https://github.com/Rosemeis/pcangsd/blob/master/scripts/pcadapt.R

lcwgs_outliers <- \(zscores, positions, K) {
  # Read in zscores.
  zscores <- read.table(zscores) %>% 
    `colnames<-`(., gsub("V", "PC", colnames(.)))
  K <- ncol(zscores)
  
  # Calcuate distances, p-values and q-values.
  dist <- dist_ogk(as.matrix(zscores))
  pval <- pchisq(dist, df = K, lower.tail = FALSE)
  qval <- qvalue(pval)$qvalues
  
  # Genomic positions for each locus.
  # positions <- read.table(positions, header = T) %>% 
  #   mutate(pos = as.numeric(gsub(".*\\.1\\_", "", marker)),
  #          chr = as.factor(gsub("\\.1\\_.*", "\\.1", marker))) %>% 
  #   dplyr::select(c(chr, pos))
  
  positions <- read_tsv(positions, col_names = c("chr", "pos"))
  
  # Bind it all together and return.
  data <- cbind(positions,
                data.frame(
                  stat = dist,
                  pval = pval,
                  qval = qval)
  )

}

lcwgs_full <- lcwgs_outliers("../data/fst/lcwgs_m15_maf005.pcadapt.zscores",
                             "../data/lcwgs_8MSNPs_positions.txt", K = 5)

lcwgs_full_n361 <- lcwgs_outliers("../data/fst/pcadapt_n361_lcwgs_pcadapt_wSites.pcadapt.zscores",
                                  "../data/chinook_lcwgs_maf1_m15_n361_snps.txt", K=5)

test <- read.delim("../data/pcadapt_n361_lcwgs_pcadapt_wSites.sites", col.names = "keep") %>% 
  cbind(., lcwgs_full_n361[1:5268574,]) %>% filter(keep == 1) %>% select(-c(keep))

ggplot(data = test, aes(x = pos, y = -log10(pval))) + geom_point() + facet_wrap(~chr, scales = "free_x")

# write.csv(lcwgs_full, "../data/fst/lcwgs_pcadapt_scores_pvals.csv", row.names = F)

lcwgs_full <- read.csv("../data/fst/lcwgs_pcadapt_scores_pvals.csv") 
 
# (ots28_full <- ggplot(data = lcwgs_full %>% 
#        filter(chr == "NC_056456.1"),
#        aes(x = pos/1e6, 
#            y = -log10(pval))) +
#   geom_point() +  theme_classic() +
#   theme(axis.title.y = element_text(size = 12),
#         panel.grid   = element_blank(),
#         legend.position  = "none",
#         strip.background = element_blank(),
#         strip.placement  = "inside",
#         strip.text = element_text(size = 12)) +
#   labs(x = "Position (Mbp)",
#        y = expression(-log[10](p-value)))) 


# greb1l <- ots28$qval<threshold & ots28$pos>=13278338 & ots28$pos<=13630075

manh_ins <- \(df, x) { 
  
  main <- ggplot(df, 
                 aes(x = {{x}}/1e6,
                     y = -log10(pval))) +
    geom_point() + theme_classic() +
    labs(x = "Position (Mbp)",
         y = expression(-log[10](p-value)))

  # ymax <- as.numeric(max(-log10(df$pval), na.rm = TRUE))

  mag <- main + coord_cartesian(clip = "off") +
    theme(plot.margin = ggplot2::margin(10,100,10,10)) +
    geom_magnify(from = c(13.8,14,0,40),
                 to   = c(48,56,0,80))

}

(lcwgs <- manh_ins(df = lcwgs_full_n361[lcwgs_full_n361$chr == "NC_056456.1",], x = pos))
(hdgbs <- manh_ins(df = hdgbs_imp[hdgbs_imp$CHROM == "NC_056456.1",], x = POS))

(stack <- cowplot::plot_grid(plotlist = list(lcwgs, hdgbs),
                             nrow = 2, align = "VH"))

cowplot::save_plot(plot = stack, base_width = 8, base_height = 6,
                   filename = "../plots/stacked_manhattans.png")

# 
# 
# 
# 
# j <- ots28_full +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = ggplot2::margin(10,80,10,10)) +
#   geom_magnify(from = c(37,39,0,75),
#                to   = c(48,56,0,75))
# 
# h <- ggplot(data = hdgbs_imp[hdgbs_imp$CHROM == "NC_056456.1",],
#        aes(x = POS/1e6, y = -log10(pval))) +
#   geom_point() +  theme_classic()
# 
# (h2 <- h +
#     coord_cartesian(clip = 'off', xlim = c(0, 40)) +
#     geom_magnify(from = c(37,39,0,18),
#                  to   = c(48,56,0,18)))
# 
# 
# 
# cowplot::plot_grid(plotlist = list(j,h2),
#                    ncol = 1, nrow = 2)
# 
# ggsave("../plots/stacked_manhattans.tiff", dpi = 300,
#        width = 12, height = 8)


