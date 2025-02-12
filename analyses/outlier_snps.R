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
# vcf2pc("../data/vcfs/hdgbs_n361_maf5_m70.vcf.gz")
# vcf2pc("../data/vcfs/chinook_lcwgs_maf1_m15_n361_imputed.vcf.gz")
# vcf2pc("../data/vcfs/lcWGS_full_8MSNPs_imputed.vcf.gz")
# vcf2pc("../data/vcfs_n361/lcwgs_maf5.gz")
vcf2pc("../data/vcfs_n361/hdgbs_maf5_m70.vcf.gz")


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
hdgbs  <- pcadapt2("../data/vcfs_n361/hdgbs_maf5_m70.vcf.gz", K = 5, q.alpha = 1/100)


# hdprune <- pcadapt2("../data/vcfs_n361/hdgbs_maf5_m70_pruned.vcf.gz", K = 5, q.alpha = 1/100)
#hdorig  <- pcadapt2("../data/vcfs_n361/hdgbs_maf5_m70.vcf.gz", K = 5, q.alpha = 1/100)
# hdgbs_ori <- pcadapt2("../data/vcfs/hdgbs_n361_maf5_m70.vcf.gz", K = 5, q.alpha = 1/100)
# lcwgs_imp <- pcadapt2("../data/vcfs/lcWGS_full_8MSNPs_imputed.vcf.gz", K = 5, q.alpha = 1/100)
lcwgs_ori <- pcadapt2("../../chinook_filtered_maf5.vcf.gz", K = 5, q.alpha = 1/100)
# lcwgs_imp2 <- pcadapt2("../data/vcfs/chinook_lcwgs_maf1_m15_n361_imputed.vcf.gz", K = 5, q.alpha = 1/100)

# BELOW NEED TO ADD SIGNIFICANCE THRESHOLD

# Define region of Ots28 containing GREB1L/ROCK1.
grebrock <- c(13278338:13630075)

Ots28man <- \(df) {
  
  # Maximum "height" on y-axis for points within GREB1L/ROCK1.
  gr <- max(-log10(df[df$CHROM == "NC_056456.1" & df$POS %in% grebrock, "pval"]))
  mp <- median(grebrock)/1e6 # Median of that region. For plotting a label.
  
  # First plot everything that isn't in GREB1L/ROCK1 for cleaner visualization.
  ggplot(data = df[df$CHROM == "NC_056456.1" & !df$POS %in% grebrock,],
         aes(x = POS/1e6, y = -log10(pval),
             colour = greb)) +
    labs(x = "Ots28 Position (Mb)",
         y = expression(-log[10](p-value))) +
    geom_point(color = "gray") +
    theme_classic() +
    # Two annotation calls for labelling GREB1L, etc.
    annotate("text", label="GREB1L/ROCK1",
             x = mp, y = gr*2.15,angle=45) +
    annotate("segment", x = mp, y = gr*2, xend = mp, yend = gr+gr*0.1,
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
    theme(legend.position = "none") +
    # Fill in points for GREB1L/ROCK1 in black so they are easier to see.
    geom_point(data = df[df$CHROM == "NC_056456.1" & df$POS %in% grebrock,],
               aes(x = POS/1e6, y = -log10(pval)), inherit.aes = F, color = "black")
}

(hdgbs28 <- Ots28man(hdgbs))
# (lcwgs   <- Ots28man(lcwgs_ori))

# ggplot(data = lcwgs_ori[lcwgs_ori$CHROM == "NC_056456.1",], aes(x = POS/1e6, y = -log10(pval))) + geom_point() + theme_bw() + geom_vline(xintercept = 13.6)
# ggplot(data = hdgbs[hdgbs$CHROM == "NC_056456.1",], aes(x = POS, y = -log10(pval))) + geom_point()
# ggplot(data = hdorig[hdorig$CHROM == "NC_056456.1",], aes(x = POS, y = -log10(pval))) + geom_point()

# # For better plotting, read in chromosome info.
# chrs <- read.delim("../data/otsh_sequence_report.tsv") %>% 
#   filter(!Chromosome.name %in% c("MT", "Un")) %>% 
#   mutate(Ots = gsub("LG", "Ots", Chromosome.name),
#          LGN = as.factor(str_sub(Ots, 4, 5)))
# 
# 
# # Isolate Ots28 for current purposes.
# Ots28 <- rbind(hdgbs_imp[hdgbs_imp$CHROM == "NC_056456.1",],
#                hdgbs_ori[hdgbs_ori$CHROM == "NC_056456.1",],
#                lcwgs_imp[lcwgs_imp$CHROM == "NC_056456.1",]) 
# 
# (chromplot <- ggplot(data = Ots28, 
#                      aes(x = POS/1e6,
#                          y = -log(pval),
#                          colour = out)) +
#     scale_color_manual(values = c("gray70", "red1")) +
#     geom_point() +
#     theme_classic() +
#     labs(x = "Position (Mbp)",
#          y = expression(-log[10](p-value))) +
#     theme(axis.title.y = element_text(size = 12),
#           panel.grid   = element_blank(),
#           legend.position  = "none",
#           strip.background = element_blank(),
#           strip.placement  = "inside",
#           strip.text = element_text(size = 12)) +
#     facet_wrap(~data, ncol = 1, scales = "free_x"))
# 
# 
# ggsave("../plots/pcadapt_manhattan.tiff",
#        dpi = 300, width = 10, height = 10)
# 

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

  # Read in genomic positions.
  positions <- read_tsv(positions, col_names = c("CHROM", "POS"))
  
  # Bind it all together and return.
  data <- cbind(positions,
                data.frame(
                  stat = dist,
                  pval = pval,
                  qval = qval)
                )

}

lcwgs_full <- lcwgs_outliers("../data/pcadapt/lcwgs_n361_pcadapt_7h2.pcadapt.zscores",
                             "../data/lcwgs_snp_positions_n361_7M.txt", K = 5)

(lcma28 <- Ots28man(lcwgs_full))

(mans <- cowplot::plot_grid(plotlist = list(hdgbs28, lcma), 
                            ncol = 1, labels = c("hdGBS", "lcWGS"), label_x = 0.025))

ggsave("../plots/ots28_manhattans.tiff", dpi = 300, width = 12, height = 8)

