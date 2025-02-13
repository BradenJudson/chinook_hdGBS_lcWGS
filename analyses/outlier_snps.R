setwd("~/ots_landscape_genetics/analyses")

library(tidyverse); library(pcadapt); library(qvalue); 
library(cowplot); library(bigutilsr); library(ggmagnify)

# Function for converting VCF to bed (usable by pcadapt).
# Assumes VCF is gzipped.
vcf2pc <- \(vcf_file) system(paste("plink.exe --vcf", vcf_file, 
                           " --make-bed --aec --double-id --out", 
                           gsub("\\.vcf.gz", "", vcf_file)))

# Make bed files necessary for pcadapt.
vcf2pc("../data/vcfs_n361/hdgbs_maf5_m70.vcf.gz")
vcf2pc("../data/vcfs_n361/chinook_filtered_maf5_imputed.vcf.gz")


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
lcimp  <- pcadapt2("../data/vcfs_n361/chinook_filtered_maf5_imputed.vcf.gz", K = 5, q.alpha = 1/100)
# write.csv(lcimp, "../data/lcwgs_n361_7Msnps_pcadapt.csv", row.names = F)

# BELOW NEED TO ADD SIGNIFICANCE THRESHOLD

# Define region of Ots28 containing GREB1L/ROCK1.
grebrock <- c(13278338:13630075)

Ots28man <- \(df) {
  
  # Maximum "height" on y-axis for points within GREB1L/ROCK1.
  gr <- max(-log10(df[df$CHROM == "NC_056456.1" & df$POS %in% grebrock, "qval"]))
  mp <- median(grebrock)/1e6 # Median of that region. For plotting a label.
  
  # First plot everything that isn't in GREB1L/ROCK1 for cleaner visualization.
  ggplot(data = df[df$CHROM == "NC_056456.1" & !df$POS %in% grebrock,],
         aes(x = POS/1e6, y = -log10(qval))) +
    labs(x = "Ots28 Position (Mb)",
         y = expression(-log[10](q-value))) +
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
               aes(x = POS/1e6, y = -log10(qval)), inherit.aes = F, color = "black")
}

(hdgbs28 <- Ots28man(hdgbs))
(lcimp28 <- Ots28man(lcimp))

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

(mans <- cowplot::plot_grid(plotlist = list(hdgbs28, lcma28, lcimp28), 
                            labels = c("hdGBS", "lcWGS", "Imputed lcWGS"), 
                            label_x = 0.025, ncol = 1))

ggsave("../plots/ots28_manhattans.tiff", dpi = 300, width = 12, height = 8)

