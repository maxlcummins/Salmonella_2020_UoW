library(pheatmap)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(reshape2)
library(tidytree)
library(ggtree)

#### Path Definitions/Config ####
### IMPORTANT ###
### Check lines where there are hardcoded changes, such as line 21,
### if you plan to modify this script's inputs.

#Replace the variable below with the path to your SG17-135 repo
wrk_dir <-
  "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/AVC171/AVC171"

#Changes working directory
setwd(wrk_dir)

#source("scripts/gene_data.R")

#Paths to input files
path_to_tree <- "~/Desktop/core_gene_alignment_snp_sites.aln.tree"
path_to_abricate <- "analysis/abricate/ST95_all/pAVC171-IncF.txt"
plasrefname <- "pAnkS"
treerefname <- "AVC171"

#Minimum hit thresholds
#Minimum hit length (as a percentage [i.e. 0.5 = 0.5%])
min_hit_length <- 0.5
#Minimum nucleotide ID (also as a percentage [i.e. 90 = 90%])
min_hit_id <- 90

#### Read in, process and subset abricate data ####
#Read in the abricate genotype data sheet
#(small number of rows for colname reassignment)
#This is to reduce memory requirements
abricate_hits <-
  read_delim(
    path_to_abricate,
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    n_max = 10
  )

#Colname reassignment
colnames(abricate_hits)[c(1, 10:11)] <-
  c("name", "perc_coverage", "perc_identity")

#Extract column names for later reassignment
abricate_hits_colnames <- colnames(abricate_hits)

#Re-read in PAI abricate genotype data sheet
abricate_hits <-
  read_delim(
    path_to_abricate,
    "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    col_names = FALSE,
    skip = 1
  )

#Remove cases where there are multiple headers from concatenation of abricate reports
abricate_hits <- abricate_hits %>% filter(X2 != "SEQUENCE")


#Colname reassignment
colnames(abricate_hits) <- abricate_hits_colnames

#Convert percent coverage and identity to numeric type to allow filtering
abricate_hits$perc_coverage <- as.numeric(abricate_hits$perc_coverage)
abricate_hits$perc_identity <- as.numeric(abricate_hits$perc_identity)

#Filter to perc_identity > 95%
#abricate_hits <-
abricate_hits <- abricate_hits %>% filter(perc_identity > min_hit_id)
abricate_hits <- abricate_hits %>% filter(perc_coverage > min_hit_length)

#Trim excess characters the assembly names and reassign this to rownames
abricate_hits$name <- gsub("\\..*", "", abricate_hits$name)

#Read in the tree file
tree <-
  read.tree(file = path_to_tree)

#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

#Replace Reference with refname
tree$tip.label <- gsub("Reference", treerefname, tree$tip.label)

#Subset the hits to strains within the tree (this saves memory)
abricate_hits <- abricate_hits %>% filter(name %in% tree$tip.label)

#Extract from coverage column the coverage range and the length of the reference
abricate_hits$coverage_range <- gsub("\\/.*","", abricate_hits$COVERAGE)
abricate_hits$ref_length <- gsub(".*\\/","", abricate_hits$COVERAGE)

#Save length of plasmid reference as a variable for later use
ref_length <- as.numeric(unique(abricate_hits$ref_length))

#Replace the '-' with a ':' in the coverage
abricate_hits$coverage_range <- gsub("-",":", abricate_hits$coverage_range)

#Create a column for start and end coordinates of hits
abricate_hits$end <- gsub("[0-9]+:","", abricate_hits$coverage_range)
abricate_hits$start <- gsub(":[0-9]+","", abricate_hits$coverage_range)

#Select columns of interest
abricate_hits <- abricate_hits %>%
  select(name, gene = GENE, ref_length, start, end, percentage = perc_coverage)

#Convert start and end coordinates to numeric
abricate_hits$start <- as.numeric(abricate_hits$start)
abricate_hits$end <- as.numeric(abricate_hits$end)

#Create an empty matrix equal to length of ref plasmid
empty_plasrow <- rep(0, times = unique(abricate_hits$ref_length))

#Create an empty matrix with n rows (n = sample size) with ncol == length(ref plasmid)
empty_plasmatrix <- matrix(rep(empty_plasrow,
                               times = length(unique(abricate_hits$name))),
                           nrow = length(unique(abricate_hits$name)))

#Create a list of levels for sample names, a list of start coords and a list of end coords
#and bind these in a list of lists
start_ends <- list(as.list(as.integer(as.factor(abricate_hits$name))),
                   as.list(as.integer(abricate_hits$start)),
                   as.list(as.integer(abricate_hits$end)))

#Create a counter
counter <- 0

#Map the BLAST hits to our matrix of bp coordinates
for (i in 1:nrow(abricate_hits)){
  sample <- start_ends[[1]][[i]]
  start_coord <- start_ends[[2]][[i]]
  end_coord <- start_ends[[3]][[i]]
  empty_plasmatrix[sample, start_coord:end_coord] <- 1
  counter <- counter + 1
  if(counter %% 1000 == 0){
    message(paste(counter, "out of ", nrow(abricate_hits), "hits processed"))}
}

#Rename matrix
base_matrix <- empty_plasmatrix

#Remove old matrix
rm(empty_plasmatrix)

# Convert matrix to a dataframe
base_matrix <- as.data.frame(base_matrix, stringsAsFactors = FALSE)

# Assign sample names to rows
rownames(base_matrix) <- unique(abricate_hits$name)

# Get the length of the reference sequence
df_length <- length(abricate_hits)

# Generate indices that cover blocks of 100 columns in base_matrix
bin_ranges <- c(seq(from = 1, to = ref_length, by = 100))

bin_ranges2 <- c(seq(from = 100, to = ref_length, by = 100), ref_length)


# Split the ranges into two lists of  [1] start and [2] end indices 
bin_splits <- list(bin_ranges, bin_ranges2)

# Initialise empty vector for loop below
binned_hits <- vector()

#Create a counter
counter <- 0

# Binning loop
for (i in 1:length(bin_ranges2)){
    # If the last bin is only 1 base long then the as.matrix line won't work,
    #so we have to include the if statement below:
  if(i == length(bin_ranges) & bin_ranges[length(bin_ranges)] == bin_ranges2[length(bin_ranges2)]){
    row_sum <- 1
  }else{
    # Generate row sums (i.e. Number of matching bases) for 100 column chunks of the base_matrix
  row_sum <- as.matrix(rowSums(base_matrix[,bin_splits[[1]][i]:bin_splits[[2]][i]]))
  }
  # Bind them together in a new vector
  binned_hits <- cbind(binned_hits, row_sum)
  counter <- counter + 1
  if(counter %% 50 == 0){
    message(paste(counter, " samples out of ", length(unique(abricate_hits$name)), " processed" ))}
}

#Save rownames
nems <- rownames(binned_hits)

#Convert binned_hits to a data frame
binned_hits <- as.data.frame(binned_hits)

df6 <- as.data.frame(rowSums(binned_hits))

df6$working_names <- rownames(df6)

colnames(df6) <- c("Percent_hit","working_name")

df6$Percent_hit <- round((df6$Percent_hit/ref_length) * 100)


# Write the DF
write.csv(x = binned_hits, file = paste0("delims/",plasrefname, "_plasmid_coverage.csv"), row.names = TRUE)

write.csv(x = df6, file = paste0("delims/",plasrefname, "_plasmid_coverage_percentage.csv"), row.names = TRUE)
