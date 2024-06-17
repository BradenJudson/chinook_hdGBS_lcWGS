setwd("~/ots_landscape_genetics/lcWGS/stats")

library(tidyverse); library(ggpmisc); library(ggExtra); library(cowplot)


# Reads and coverage -----------------------------------------------------------


cov <- read.table("bam_coverage.txt", sep = "",header = F, 
                  col.names = c("file", "coverage"))

reads <- read.table("reads.txt", sep = "", header = F,
                    col.names = c("file", "reads")) %>% 

d <- read.table("../../data/bam_list_n453.txt", header = F)
  
  
  
  mutate(sample = gsub(".dedup.clip.bam", "", 
                     gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file)))

dat <- merge(cov, reads, by = "file") %>% 
  file = gsub(".dedup.clip.bam", "", 
              gsub("^.*IDT_i\\d{1,3}_\\d{1,3}\\.", "", file))

(lre <- ggplot(data  = dat,
               aes(x = reads/1e6,
                   y = coverage)) + 
    scale_x_continuous(breaks = seq(0, 90, 10)) +
    scale_y_continuous(breaks = c(0, seq(1,3,1/2), 4)) +
    coord_cartesian(clip = 'on') +
    geom_segment(aes(y = mean(`coverage`), 
                     yend = mean(`coverage`), 
                     x = -Inf, xend = Inf),
                 linewidth = 1/2, linetype = 2) +
    geom_segment(aes(x = mean(`reads`)/1e6,
                     xend = mean(`reads`)/1e6,
                     y = -Inf, yend = Inf),
                 linewidth = 1/2, linetype = 2) +
    geom_smooth(method = "lm",
                colour = "red",
                alpha  = 1/6,
                linetype = 1) +
    geom_point(colour = "black",
               shape  = 21, 
               fill   = "gray80",
               size   = 2,
               alpha  = 4/5) +
    theme_bw() +
    theme(panel.grid = element_line(color = "gray95"),
          panel.grid.minor.y = element_blank()) +
    labs(x = "Reads (millions)", 
         y = "Average Individual Coverage") +
    stat_poly_eq(use_label(c("R2", "p")),
                 label.x = "left",
                 label.y = "top",
                 small.p = TRUE))

# Add density plots of marginal distributions.
(mds <- ggMarginal(lre, type = "density", colour = "black", 
                   linewidth = 1, fill = "gray90", alpha = 0.7))

save_plot("../../plots/chinook_lcWGS_reads_coverage.tiff", mds, ncol = 2,
          base_height = 6, base_asp = 1)


(sstat <- dat %>% 
    pivot_longer(cols = c("reads", "coverage")) %>% 
    group_by(name) %>% 
    summarise(mean_value = mean(value),
              med_value  = median(value),
              max_val    = max(value),
              min_val    = min(value),
              sd         = sd(value)))


# Individuals with lowest reads and/or coverage.
dat %>% arrange(reads) %>% head() %>% select(c(1:3L))



# Population-level stats -------------------------------------------------------

samples <- read.csv("../../data/landgen_chinook_indvs.csv")[,c(1,3)] %>% 
  group_by(fish_ID) %>% sample_n(1) %>% 
  mutate(fish_ID = gsub("\\.", "-", fish_ID),
         site_full = gsub("[[:space:]]Brood", "", site_full))

datpop <- merge(dat, samples, by.y = "fish_ID", by.x = "file") %>% 
  group_by(site_full) %>% 
  summarise(coverage_variance = var(coverage),
            mean_coverage = mean(coverage),
            min_coverage  = min(coverage))


write.csv(datpop, "ch_lcwgs_covvar.csv", row.names = F)










