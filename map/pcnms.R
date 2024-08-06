args <- commandArgs(T)

# Confirm correct usage --------------------------------------------------------

if (length(args) != 2)
  stop("Usage: dbMEM.R <DISTMAT> <ALLELEFREQS>")

# Install and load necessary packages ------------------------------------------

for (p in c("vegan","tidyverse","data.table")) {
  if (!suppressMessages(require(p, character.only = T))) {
    message(paste("Installing:", p))
    install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)
    suppressMessages(require(p, character.only = T))}
  rm(p)
}

# Allele frequencies -----------------------------------------------------------

freqs <- as.data.frame(fread(file = args[2], header=T, sep="\t")) %>% 
  column_to_rownames(1)

message("read allele frequencies")

# Distance matrix --------------------------------------------------------------

rivdist <- as.data.frame(fread(file = args[1], header = T, sep = "\t")) %>% 
  column_to_rownames(1)

rivdist <- read.csv("../data/Otsh_distances_mat.csv", row.names = 1) %>% 
  `colnames<-`(., gsub("\\.", "", colnames(.))) %>% 
  `rownames<-`(., gsub(" ", "", rownames(.))) 

rivdist <- data.matrix(rivdist)
rivdist <- rivdist[rownames(rivdist) %in% rownames(freqs),
                   colnames(rivdist) %in% rownames(freqs)]

pcnms <- as.data.frame(pcnm(dis = rivdist)$vectors) 

message("PCNMs calculated")

# PCNM Selection ---------------------------------------------------------------

nmod  <- rda(freqs ~ 1, pcnms) # Null model.
fmod  <- rda(freqs ~ ., pcnms) # Full model.
optpcnm <- ordiR2step(nmod, scope = formula(fmod), direction = 'forward')
optcall <- summary(optpcnm)$call[2] # Prints formula with selected predictor variables.

message("Optimal PCNMS determined")

# Isolate optimal PCNMs and write to txt file.
(vars <- str_split(gsub(",(.*)","",gsub(".*\\~ ", "", optcall)), " \\+ ")[[1]])
write.table(vars, file = "optimal_pcnms.txt", quote = F, row.names = F, col.names = F)




