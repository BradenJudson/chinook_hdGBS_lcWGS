args <- commandArgs(T)

# Confirm correct usage --------------------------------------------------------

if (length(args) != 2)
  stop("Usage: pcnm.R <DISTMAT> <ALLELEFREQS>")

# Install and load necessary packages ------------------------------------------

for (p in c("vegan","tidyverse","data.table")) {
  if (!suppressMessages(require(p, character.only = T))) {
    message(paste("Installing:", p))
    install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)
    suppressMessages(require(p, character.only = T))}
  rm(p)
}

# Allele frequencies -----------------------------------------------------------

message("reading allele frequencies")

# fread corrupts some of the population text, so we use awk
# to isolate and re-assign the population names.
system(paste("awk '{print $1}'", args[2], "> populations.txt"))
freqs <- as.data.frame(fread(file = args[2]))[,-1]
print(str(freqs))
pops <- read.table("populations.txt", header = TRUE, quote = "")
rownames(freqs) <- c(pops[,1])
print(dim(freqs))
print(freqs[1:5, 1:5])
print(unique(rownames(freqs)))

message("read allele frequencies")

# Distance matrix --------------------------------------------------------------

rivdist <- read.csv(file = args[1], row.names = 1)

print(dim(rivdist))

message(paste0("Distance matrix for ", nrow(rivdist), " populations used"))

rivdist <- data.matrix(rivdist)
rivdist <- rivdist[rownames(rivdist) %in% rownames(freqs),
                   colnames(rivdist) %in% rownames(freqs)]

print(dim(rivdist))
pcnms <- as.data.frame(pcnm(dis = as.dist(rivdist))$vectors)

message("PCNMs calculated")

# PCNM Selection ---------------------------------------------------------------

nmod  <- rda(freqs ~ 1, pcnms) # Null model.
fmod  <- rda(freqs ~ ., pcnms) # Full model.
optpcnm <- ordiR2step(nmod, scope = formula(fmod), direction = 'forward')
optcall <- summary(optpcnm)$call[2] # Prints formula with selected predictor variables.

message("Optimal PCNMS determined")

# Isolate optimal PCNMs and write to txt file.
(vars <- str_split(gsub(",(.*)","",gsub(".*\\~ ", "", optcall)), " \\+ ")[[1]])
write.table(vars, file = "data/optimal_pcnms.txt", quote = F, row.names = F, col.names = F)
