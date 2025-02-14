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

# DON'T FORGET TO ADD SIGNIFICANCE THRESHOLD - MAYBE Q<0.05?

# Define region of Ots28 containing GREB1L/ROCK1.
grebrock <- c(13278338:13630075)

Ots28man <- \(df) {
  
  # Maximum "height" on y-axis for points within GREB1L/ROCK1.
  gr <- max(-log10(df[df$CHROM == "NC_056456.1" & df$POS %in% grebrock, "qval"]))
  mp <- median(grebrock)/1e6 # Median of that region. For plotting a label.
  
  # First plot everything that isn't in GREB1L/ROCK1 for cleaner visualization.
  ggplot(data = df[df$CHROM == "NC_056456.1" ,],
         aes(x = POS/1e6, y = -log10(qval), colour = out)) +
    labs(x = "Ots28 Position (Mb)",
         y = expression(-log[10](q-value))) +
    geom_point() +
    scale_colour_manual(values = c("gray80", "black")) +
    theme_classic() +
    # Two annotation calls for labelling GREB1L, etc.
    annotate("text", label="GREB1L/ROCK1",
             x = mp, y = gr*2.15,angle=45) +
    annotate("segment", x = mp, y = gr*2, xend = mp, yend = gr+gr*0.1,
             arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
    theme(legend.position = "none") 
    # # Fill in points for GREB1L/ROCK1 in black so they are easier to see.
    # geom_point(data = df[df$CHROM == "NC_056456.1" & df$POS %in% grebrock,],
    #            aes(x = POS/1e6, y = -log10(qval)), inherit.aes = F, color = "black")
}

(hdgbs28 <- Ots28man(hdgbs))
(lcimp28 <- Ots28man(lcimp))

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

(mans <- cowplot::plot_grid(plotlist = list(hdgbs28 + xlab(NULL),
                                            lcma28 + xlab(NULL), lcimp28), 
                            labels = c("hdGBS", "lcWGS", "Imputed lcWGS"), 
                            label_x = 0.05, ncol = 1, hjust = 0))

ggsave("../plots/ots28_manhattans.tiff", dpi = 300, width = 12, height = 8)


# Outlier comparisons ----------------------------------------------------------

# osnp <- \(x) paste0(x[x$out == "Y", "CHROM"], 
#                     x[x$out == "Y", "POS"])

common_snps <- data.frame(snp = c(paste0(lcimp$CHROM, lcimp$POS),
                                 paste0(hdgbs$CHROM, hdgbs$POS))) %>% 
  group_by(snp) %>% tally() %>% filter(n > 1) 


outliers <- list(
  lcimp = lcimp,
  lcwgs = lcwgs_full
)

ggVennDiagram(outliers)



# Regions under selection ------------------------------------------------------

library(GenomicRanges)

# gtf <- read.table("../data/GCF_018296145.1_Otsh_v2.0_genomic.gtf", 
#                   header=F, stringsAsFactors=F, sep="\t") %>% 
#   mutate(gene = gsub("gene_id ", "", gsub(";.*$", "", V9))) %>% 
#   filter(V3 == "gene") %>% 
#   column_to_rownames("gene") %>% 
#   select(c(1,4:5,7)) %>% 
#   `colnames<-`(., c("chromosome", "start", "end", "strand"))
#   
# # Assemble GRanges object.
# ggr <- GRanges(seqnames = gtf$chromosome, 
#                ranges = IRanges(start = gtf$start, end = gtf$end)) %>% 
#   `names<-`(., rownames(gtf)) 
# strand(ggr) <- gtf$strand
# 
# # Define promoters.
# gwpgr <- promoters(ggr, upstream = 2000, downstream = 200, use.names = T)

hd_out <- hdgbs[hdgbs$out == "Y" & hdgbs$CHROM == "NC_056456.1",] %>% filter(!is.na(pval))
im_out <- lcimp[lcimp$out == "Y" & lcimp$CHROM == "NC_056456.1",] %>% filter(!is.na(pval))

df2GR <- \(x, mgw) GRanges(seqnames = x$CHROM,
                      ranges = IRanges(start = x$POS,
                                       end   = x$POS)) %>% 
                  `names<-`(., x$loc) %>% 
                   reduce(., min.gapwidth = mgw,
                          drop.empty.ranges = F,
                          with.revmap = F, 
                          with.inframe.attrib = F)


hd_pval <- df2GR(hd_out, mgw = 2e4)

lcimp_p <- df2GR(im_out, mgw = 2e4)

test2 <- subsetByOverlaps(hd_pval, lcimp_p, type = "any") 

test <- Reduce(subsetByOverlaps, list(hd_pval, lcimp_p))
q <- as.data.frame(test2@ranges)

nrow(as.data.frame(hd_pval@ranges)); nrow(as.data.frame(lcimp_p@ranges)); nrow(q)

k <- data.frame(
  min_gapwidth = c(10000, 20000, 30000, 40000, 50000),
  hdgbs_regs   = c(105, 96, 87, 82, 77),
  lcwgs_regs   = c(618, 390, 285, 241, 205),
  overlaps     = c(75, 75, 72, 67, 66)
) 




