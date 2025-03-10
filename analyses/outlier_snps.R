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
  x <- pcadapt(input = input, K = K, min.maf = 0.05)
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
           )),
           loc = paste0(CHROM, "_", POS))
  
  print(summary(snps$out)); return(snps)
  
}

# Conduct pcadapt with specified parameters. Returns dataframe of all SNPs and their outlier scores/status.
hdgbs  <- pcadapt2("../data/vcfs_n361/hdgbs_maf5_m70.vcf.gz", K = 5, q.alpha = 1/100)
lcimp  <- pcadapt2("../data/vcfs_n361/chinook_filtered_maf5_imputed.vcf.gz", K = 5, q.alpha = 1/100)
# write.csv(lcimp, "../data/lcwgs_n361_7Msnps_pcadapt.csv", row.names = F)
# lcimp <- read.csv("../data/lcwgs_n361_7Msnps_pcadapt.csv")

genome_manhattan <- \(ds) {
  
  ggplot(data = ds,
         aes(x = POS/1e6, y = -log10(qval))) + 
    geom_point() + facet_wrap(~CHROM, scales= "free_x") +
    theme_bw() + labs(x = "Position (Mbp)")
  
}

(hdgbs_GW <- genome_manhattan(hdgbs))
ggsave("../plots/supp_figs/hdgbs_manhattan.tiff", dpi = 300, width = 15, height = 10)
(lcimp_GW <- genome_manhattan(lcimp))
ggsave("../plots/supp_figs/lcimp_manhattan.tiff", dpi = 300, width = 15, height = 10)

# ------------------------------------------------------------------------------

# p-value calculations from https://github.com/Rosemeis/pcangsd/blob/master/scripts/pcadapt.R

lcwgs_outliers <- \(zscores, positions, q.alpha) {
  
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
                ) %>% 
    mutate(out = case_when(qval >= q.alpha ~ "N",
                           qval <  q.alpha ~ "Y"),
           loc = paste0(CHROM, "_", POS))

}

lcwgs_full <- lcwgs_outliers("../data/pcadapt/lcwgs_n361_pcadapt_7h2.pcadapt.zscores",
                             "../data/lcwgs_snp_positions_n361_7M.txt", q.alpha = 1/1000)

(lcma28 <- Ots28man(lcwgs_full))

# (mans <- cowplot::plot_grid(plotlist = list(hdgbs28 + xlab(NULL),
#                                             lcma28 + xlab(NULL), lcimp28), 
#                             labels = c("hdGBS", "lcWGS", "Imputed lcWGS"), 
#                             label_x = 0.05, ncol = 1, hjust = 0))

ggsave("../plots/ots28_manhattans.tiff", dpi = 300, width = 12, height = 8)


# Regions under selection ------------------------------------------------------

library(GenomicRanges)

# Define region of Ots28 containing GREB1L/ROCK1.
grebrock <- c(13278338:13630075)

# Isolate Ots28 for this analysis and only consider "true" outliers.
hd_out <- hdgbs[hdgbs$out == "Y" & hdgbs$CHROM == "NC_056456.1",] %>% filter(!is.na(pval)) 
im_out <- lcimp[lcimp$out == "Y" & lcimp$CHROM == "NC_056456.1",] %>% filter(!is.na(pval))

nrow(hd_out); nrow(im_out)

# Compute genomic ranges of outliers with a given minimum gapwidth.
df2GR <- \(x, mgw) GRanges(seqnames = x$CHROM,
                      ranges = IRanges(start = x$POS,
                                       end   = x$POS)) %>% 
                  `names<-`(., x$loc) %>% 
                   reduce(., min.gapwidth = mgw,
                          with.revmap = F)


# Define outlier regions for each dataset.
hd_pval <- df2GR(hd_out, mgw = 100e3)
lcimp_p <- df2GR(im_out, mgw = 100e3)

# Isolate regions that overlap between the two datasets.
genOverlap <- subsetByOverlaps(hd_pval, lcimp_p, type = "any") 

nrow(as.data.frame(hd_pval@ranges))
nrow(as.data.frame(lcimp_p@ranges))
nrow(as.data.frame(genOverlap@ranges))

lcwgs_only <- subsetByOverlaps(lcimp_p, genOverlap, invert = T)
hdgbs_only <- subsetByOverlaps(hd_pval, genOverlap, invert = T)

reg_widths <- rbind(
  as.data.frame(lcwgs_only@ranges) %>% mutate(ds = "Imputed lcWGS only"),
  as.data.frame(hdgbs_only@ranges) %>% mutate(ds = "hdGBS only"),
  as.data.frame(genOverlap@ranges) %>% mutate(ds = "Shared regions")
)

waov <- aov(width ~ ds, data = reg_widths); summary(waov)
ggplot(data = reg_widths, aes(x = ds, y = log10(width))) + geom_boxplot()


# Create a vector of all possible SNPs that fall in the overlapping 
# regions defined above.
overlap_snps <- rowwise(as.data.frame(genOverlap@ranges)) %>%  
  mutate(seq = list(seq(start, end, by = 1))) %>% 
  select(seq) %>% unlist() %>% as.vector()

# Function for visualizing shared and non-shared outlier SNPs.
Ots28man <- \(df, shared_list, alpha) {
  
  df <- df[!is.na(df$pval),] %>% # Create a column for the SNP - is the outlier unique or not?
    mutate(outlier = factor(case_when(POS %in% shared_list & out == "Y" ~ "Shared outlier",
                                      POS %in% shared_list & out == "N" ~ "Non-outlier",
                                     !POS %in% shared_list & out == "Y" ~ "Unique outlier",
                                     !POS %in% shared_list & out == "N" ~ "Non-outlier"),
                            levels = c("Shared outlier", "Unique outlier", "Non-outlier")))

  ggplot(data = df[df$CHROM == "NC_056456.1",],
         aes(x = POS/1e6, y = -log10(qval), 
             fill = outlier, alpha = outlier)) +
    labs(x = "Ots28 Position (Mb)",
         y = expression(-log[10](q-value))) +
    geom_point(pch = 21, stroke = NA, size=2) +
    scale_fill_manual(values = c(`Non-outlier` = "gray", 
                                 `Shared outlier` = "black",
                                 `Unique outlier` = "skyblue2")) +
    scale_alpha_manual(values = c(1, alpha, alpha)) +
    geom_point(data = df[df$CHROM == "NC_056456.1" & df$outlier == "Shared outlier",],
               aes(x = POS/1e6, y = -log10(qval), fill = "black"), inherit.aes = F) +
    theme_classic() + 
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text  = element_text(size = 14),
          legend.box.background = element_rect(colour = "black")) +
    coord_cartesian(clip = "off", xlim = c(0, 43.2)) + 
    theme(plot.margin = ggplot2::margin(15, 175, 10, 10)) +
    geom_magnify(from = c(13.3, 13.6, 0, 8),
                 to   = c(47, 57, 0, 15),
                 proj.linetype = 0,
                 axes = "xy") +
    annotate("text", x = 52, y = 16.5, label = "GREB1L/ROCK1") +
    guides(fill = guide_legend(nrow = 1))
    
}

# Visualize each dataset. Adjust alpha as needed.
(hdgbs28 <- Ots28man(hdgbs, shared_list = overlap_snps, alpha = 1/2))
outlier_leg <- get_legend(hdgbs28)
(lcimp28 <- Ots28man(lcimp, shared_list = overlap_snps, alpha = 1/2))  

mans <- cowplot::plot_grid(plotlist = list(hdgbs28 + xlab(NULL) +
                                              theme(legend.position = "none"),
                                            lcimp28 +
                                              theme(legend.position = "none")), 
                            labels = c("hdGBS", "Imputed lcWGS"), 
                            label_x = 0.05, ncol = 1, hjust = 0) +
  theme(plot.margin = ggplot2::margin(20,0,0,0))
mans + draw_plot(outlier_leg, x = 1/8, y = 0.92, width = 0.6, height = 0.1)
# ConGenFunctions::insettr(mans, outlier_leg, "tl", height = 1, width = 1)

ggsave("../plots/shared_manhattans_100kb.tiff", dpi = 300, width = 12, height = 8)

 