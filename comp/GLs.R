setwd("~/ots_landscape_genetics/comp")

library(tidyverse); library(sqldf)

files <- read.table("../data/bam_list_n453.txt", col.names = c("file")) %>%
  mutate(fish_ID = gsub(".dedup.clip.bam", "",
                        gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file))) 

# dat2 <- read.csv("../data/landgen_chinook_indvs2.csv")
# 
# 
# # Apply below framework to huge csv with genotype probabilities but do it chr-by-chr 
# # and calculate histogram percentages by group. Make an output file with new column for 
# # every chromosome and have a terminal column for totals. 
# 
# dflist <- list()
# 
# for (i in c(unique(dat2$site))) {
#   
#   population <- toString(sprintf("'%s'", i))
#   
#   df <- fn$read.csv.sql(file = "../data/landgen_chinook_indvs2.csv", 
#                             sql = paste0("select * from file where site == ($population)"))
#   
#   
#   # Calculate entries per 5% bin here, make summary and wrap in function to join at the end by chromosome
#   
#   dflist[[length(dflist)+1]] <- df
# }
# 
# # df <- read.csv.sql(file = "../../maximum_original_GLs_LF.csv", eol = "\n",
# #                       sql = paste0("select * from file where POS == 398"))
# # 
# # df2 <- read.csv.sql(file = "../data/test.csv",eol = "\n",
# #                     sql = "select * from file where POS == 523")
# # 
# # 

################################################################################

# From: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018296145.1/
chrs <- read_tsv("../data/otsh_Sequence_report.tsv", show_col_types = FALSE) %>% 
  filter(!`Chromosome name` %in% c("Un", "MT")) # Remove MT and unassembled contigs.

# To populate from for-loop.
chr_list <- list()

# For each chromosome, do the following:
for (LG in c(unique(chrs$`RefSeq seq accession`[1:2]))) {
  
  # For parsing with sql.
  chromosome <- toString(sprintf("'%s'", LG))
  
  # Because this takes a while, print progress.
  print(paste("Processing", noquote(gsub("'", "", chromosome))))
  
  # Use SQLDF to read in pre-filtered chunks of the original massive csv file.
  # Necessary to specify end of line character or this fails to parse header.
  df <- fn$read.csv.sql(file = "../../maximum_original_GLs_LF.csv", eol = "\n",
                        sql = paste0("select * from file where CHR == ($chromosome)"))
  
  # Split the dataframe above and count genotype probabilities falling into 5% bins.
  # Convert the tabulated GP data and convert into a dataframe.
  bins <- as.data.frame(cut(df$GP, breaks = seq(0, 1, 1/20)) %>% table) %>% 
    `colnames<-`(., c("GP_bin", gsub("'", "", chromosome)))
  
  # Reassign to position in the parent list. 
  chr_list[[length(chr_list) + 1]] <- bins
  
}

# Join list components into a dataframe by common bin values.
# Sum values across each chromosome for genome-wide frequencies. 
(gp_bin_counts <- reduce(chr_list, full_join, by = "GP_bin") %>% 
  mutate(total  = rowSums(select(., contains("NC"))),
         rightB = as.numeric(gsub("]", "", gsub(".*,", "", GP_bin))))) %>% 
  relocate("rightB", 1) # For plotting purposes, show right edge of each frequency bin.

# Manually create a histogram from left-aligned columns and count data.
(imp <- ggplot(data = gp_bin_counts[gp_bin_counts$NC_056429.1 != 0,], 
       aes(x = rightB, y = total/1e6)) + theme_classic() +
  geom_col(just = 1, colour = "black", fill = "grey90") +
  labs(y = "SNPs (millions)", x = "Genotype probability") +
  scale_x_continuous(breaks = seq(0, 1, 1/20)) +
  theme(panel.grid.minor = element_blank()))

ggsave("../plots/imputedGPs_hist.tiff", dpi = 300, width = 10, height = 8)




            