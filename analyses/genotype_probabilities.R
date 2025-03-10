setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(sqldf)

# files <- read.table("../data/bam_list_n453.txt", col.names = c("file")) %>%
#   mutate(fish_ID = gsub(".dedup.clip.bam", "",
#                         gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file))) 

# From: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018296145.1/ (June 2024).
chrs <- read_tsv("../data/otsh_Sequence_report.tsv", show_col_types = FALSE) %>% 
  filter(!`Chromosome name` %in% c("Un", "MT")) # Remove MT and unassembled contigs.

# Have to leverage sql for this otherwise we do not have enough memory.
gl_hist <- \(csv) {

  # To populate from for-loop.
  chr_list <- list()

  # For each chromosome, do the following:
  for (LG in c(unique(chrs$`RefSeq seq accession`))) {
  
    # For parsing with sql.
    chromosome <- toString(sprintf("'%s'", LG))
  
    # Because this takes a while, print progress.
    print(paste("Processing", noquote(gsub("'", "", chromosome))))
  
    df <- fn$read.csv.sql(file = csv, eol = "\n", 
               sql = "select * from file where chr == ($chromosome) ") %>% 
          pivot_longer(cols = -c("chr", "pos"), values_to = "GP") %>% select(-c(name)) 
  
    # Split the dataframe above and count genotype probabilities falling into 5% bins.
    # Convert the tabulated GP data and convert into a dataframe.
    bins <- as.data.frame(cut(df$GP, breaks = seq(0, 1, 1/20)) %>% table) %>% 
      `colnames<-`(., c("GP_bin", gsub("'", "", chromosome)))
  
    # Reassign to position in the parent list. 
    chr_list[[length(chr_list) + 1]] <- bins
  
    rm(df); gc() 
  
  }
  
  # Join list components into a dataframe by common bin values.
  # Sum values across each chromosome for genome-wide frequencies. 
  gp_bin_counts <- reduce(chr_list, full_join, by = "GP_bin") %>% 
    mutate(total  = rowSums(select(., contains("NC"))),
           rightB = as.numeric(gsub("]", "", gsub(".*,", "", GP_bin)))) %>% 
             relocate("rightB", 1)  
  
  return(gp_bin_counts)
  
}

imputed  <- gl_hist(csv = "../data/GLs/maximum_imputed_GLs.csv")
write.csv(imputed,  "../data/imputed_histogram_data2.csv",  row.names = F)

# UPDATED FEB 25 - WAITING FOR UPDATED IMPUTED NOW.
original <- gl_hist(csv = "../data/GLs/maximum_original_GLs.csv")
write.csv(original, "../data/original_histogram_data.csv", row.names = F)

hist2 <- \(df) { ggplot(data = df[df$NC_056429.1 != 0,],
                        aes(x = rightB, y = total/1e6)) +
    theme_classic() + scale_x_continuous(breaks = seq(0, 1, 1/20)) +
    geom_col(just = 1, colour = "black", fill = "grey90") +
    labs(y = "SNPs (millions)", x = "Genotype probability") +
    theme(panel.grid.minor = element_blank()) 
  }

(imp <- hist2(imputed))
ggsave("../plots/imputedGPs_hist.tiff", dpi = 300, width = 10, height = 8)

(ori <- hist2(original))
ggsave("../plots/originalGPs_hist.tiff", dpi = 300, width = 10, height = 8)


# imp <- read.csv("../data/imputed_histogram_data.csv")

dat <- rbind(imp[,c(1:2,37)] %>% mutate(ds = "Imputed"),
           original[,c(1:2,37)] %>% mutate(ds = "Original")) %>% 
  mutate(ds = factor(ds, levels = c("Original", "Imputed")))

ggplot(data = dat[dat$total > 0,], 
       aes(x = rightB, y = total/1e6, fill = ds)) +
  geom_col(position = "dodge",  just = 1, 
           colour = "black", size = 1/5) +
  scale_x_continuous(breaks = seq(0, 1, 1/20)) +
  labs(y = "SNPs (millions)", x = "Genotype probability") +
  scale_fill_manual(values = c("grey90", "grey30")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.9))
ggsave("../plots/gp_hist.tiff", dpi = 300, width = 8, height = 6)
