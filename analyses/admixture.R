setwd("~/ots_landscape_genetics/analyses")

library(tidyverse)

qvals <- \(directory) {
 
  # List ngsAdmix output files.
  files <- paste0(directory, 
                  list.files(pattern = ".*.qopt",
                  path = directory))
  
  # Read in ngsAdmix output files and format them so they can be interpretted
  # as a "qlist" as per the pophelper workflow.
  data <- lapply(files, \(x) as.data.frame(readQ(x, filetype = "basic")) %>% 
                   `colnames<-`(., gsub(".*\\.", "", colnames(.))) )
  
  # Orders the dataframes such that they are listed with respect to the number of clusters.
  ss <- data[order(sapply(data, ncol), decreasing = F)]
  
}

# Function for reading in individual admixture values per dataset.
# qvals <- \(directory) {
#   
#   # List *.qopt files output from ngsAdmix.
#   files <- paste0(directory, 
#                   list.files(pattern = ".*.qopt",
#                   path = directory))
# 
#   # Read in each *.qopt file in a list.
#   # Files formatted to have Q1, Q2, Q3 column names and 
#   # extra whitespace is trimmed out of the file.
#   data <- lapply(files, 
#                  FUN = \(x) read_delim(file = x, 
#                             delim = " ", col_names = F) %>% 
#                    select_if(~ !any(is.na(.)) %>% 
#                    # `colnames<-`(., paste0("Q", seq(1:ncol(.)))) %>% 
#                    # rownames_to_column("id") %>% 
#                    # pivot_longer(cols = -c("id")) %>% 
#                    # group_by(id) %>% 
#                    # mutate(
#                    #   la = name[which.max(value)],
#                    #   ap = max(value)
#                    # ) %>% 
#                    # arrange(la, desc(ap)) %>% 
#                    # ungroup() %>% 
#                    # mutate(id = forcats::fct_inorder(factor(id)),
#                    #        name = factor(name, levels = paste0("Q", seq(2:14)))))
#   
#   # Orders list elements sequentially (q1, q2, q3, not q1, q10, etc.)
#   ss <- data[order(sapply(data, nrow), decreasing = F)]
#   
#   }



lcwgs_134k <- qvals(directory = "../data/admixture/lcWGS_134ksubset/")

k9 <- plotQ(as.qlist(lcwgs_134k[length(lcwgs_134k)]), returnplot = T, 
            exportplot = F, sortind = "all", ordergrp = T); k9$plot

j <- plotQ(as.qlist(lcwgs_134k[c(8,5)]), returnplot = T, 
          exportplot = F, clustercol = shiny, sortind = "all",
          imgoutput = "join", sharedindlab = F,
          splab = paste0("K = ", sapply(lcwgs_134k[c(8,5)], ncol)),
          basesize = 22, showindlab = T, ordergrp = T)
j$plot[[1]]


# j <- as.qlist(lcwgs_134k); plotQ(j, returnplot = T, exportplot = F)
# h <- alignK(j)
# plotQ(h, exportplot = FALSE, returnplot = T)
# 
# test <- readQ("../data/admixture/lcWGS_134ksubset/10.qopt", filetype = "basic")
# 
# adplot <- \(x) ggplot(data = x,
#                       aes(x = id,
#                           y = value,
#                           fill = name)) +
#   geom_col() + theme_classic() +
#   labs(x = NULL, y = "Ancestry proportion") +
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.title = element_blank(),
#         legend.position = "none") 
# 
# 
# p <- lapply(lcwgs_134k, adplot)
# 
# cowplot::plot_grid(plotlist = p, ncol = 1)
# 
# ggsave("../plots/lcwgs_134ksubset_admixture.png", dpi = 300, width = 7, height = 15)

optK <- \(directory) {
  files <- paste0(directory, 
                  list.files(pattern = ".*.log",
                             path = directory))
  data <- lapply(files, 
                 FUN = \(x) read_delim(file = x, 
                                       delim = " ", col_names = F)) %>% 
    lapply(., \(x) as.numeric(substr(x[nrow(x),2], start = 6, stop = 200)))
  
  df <- data.frame(q = as.numeric(gsub(".log", "",
                       gsub("../data/admixture/lcWGS_134ksubset/", "", files))),
                   loglike = unlist(data))
  
}

y <- optK(directory = "../data/admixture/lcWGS_134ksubset/")

ggplot(data = y, aes(x = q, y = -loglike)) +
  geom_point() + geom_line() + xlab("k") +
  scale_x_continuous(breaks = c(seq(1,15,2)))



