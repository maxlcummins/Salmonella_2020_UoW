library(pheatmap)
library(ggplot)
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
path_to_tree <- "analysis/snippy/HC50_1106/fasttree/AVC171.clean.fullcore.tree"
path_to_abricate <- "analysis/abricate/AVC171_all/pEC14_114.txt"
plasrefname <- "pEC14_114"
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

ref_length <- as.numeric(unique(abricate_hits$ref_length))

#Replace the '-' with a ':' in the coverage
abricate_hits$coverage_range <- gsub("-",":", abricate_hits$coverage_range)

#Create a column for start and end coordinates of hits
abricate_hits$end <- gsub("[0-9]+:","", abricate_hits$coverage_range)
abricate_hits$start <- gsub(":[0-9]+","", abricate_hits$coverage_range)

#Select columns of interest
abricate_hits <- abricate_hits %>%
  select(name, gene = GENE, ref_length, start, end, percentage = perc_coverage)

abricate_hits$start <- as.numeric(abricate_hits$start)
abricate_hits$end <- as.numeric(abricate_hits$end)

#Assign empty vector
base_matrix <- c()

counter <- 1

# WE APOLOGISE FOR THE FOLLOWING NESTED FOR LOOP. IT WORKS BUT 'TIS SLOW M'LADY *TIPS FEDORA* 

# Separate the abricate hits DF into separate dfs for each sample
for(sample in unique(abricate_hits$name)){
  abricate_hits2 <- abricate_hits %>% filter(name == sample)
  #Create an empty matrix for each sample with a length = to that of the reference
  range_matrix <- rep(0, times = unique(abricate_hits$ref_length))
  
# Separate out individual hits for a given sample
  
  # reassign value to 1 for given BP coordinates based on start and end coords of each hit
  for(hit in 1:nrow(abricate_hits2)){
    start_ <- abricate_hits2[[hit,4]]
    end_ <- abricate_hits2[[hit,5]]
    range_matrix[start_:end_] <- 1
  }
  
  #bind together the hit matrix for all samples
  base_matrix <- rbind(base_matrix, range_matrix)
  print(paste(counter, " iterations"))
        counter <- counter + 1
}

base_matrix2 <- base_matrix2

base_matrix <- c()

empty_plasrow <- rep(0, times = unique(abricate_hits$ref_length))

empty_plasmatrix <- matrix(rep(empty_plasmatrix,
                               times = length(unique(abricate_hits$name))),
                           nrow = length(unique(abricate_hits$name)))

start_ends <- list(as.list(as.integer(as.factor(abricate_hits$name))),
                   as.list(as.integer(abricate_hits$start)),
                   as.list(as.integer(abricate_hits$end)))

#SOLN2
for (i in 1:nrow(abricate_hits)){
  # Generate row sums (i.e. Number of matching bases) for 100 column chunks of the base_matrix
  sample <- start_ends[[1]][[i]]
  start_coord <- start_ends[[2]][[i]]
  end_coord <- start_ends[[3]][[i]]
  empty_plasmatrix[sample, start_coord:end_coord] <- 1
  print("one iteration")
  # Bind them together in a new vector
}


# Convert matrix to a dataframe
base_matrix <- as.data.frame(base_matrix, stringsAsFactors = FALSE)

# Assign sample names to rows
rownames(base_matrix) <- unique(abricate_hits$name)

# Get the length of the reference sequence
df_length <- length(abricate_hits)

# Generate indices that cover blocks of 100 columns in base_matrix
bin_ranges <- c(seq(from = 1, to = ref_length, by = 99))

# Add the last two non-divisible co-ordinates for the remaining bases
bin_ranges1 <- c(bin_ranges, bin_ranges[length(bin_ranges)]+1, ref_length)

# Split the ranges into two lists of  [1] start and [2] end indices 
bin_splits <- split(bin_ranges1, 1+(seq_along(bin_ranges1)-1) %% 2)

# Initialise empty vector for loop below
binned_hits <- vector()

counter <- 0

# Binning loop
for (i in seq(1,length(bin_ranges1)/2)){
  # Generate row sums (i.e. Number of matching bases) for 100 column chunks of the base_matrix
  row_sum <- as.matrix(rowSums(base_matrix[,bin_splits[[1]][i]:bin_splits[[2]][i]]))
  # Bind them together in a new vector
  binned_hits <- cbind(binned_hits, row_sum)
  print(paste(counter, " iterations"))
  counter <- counter + 1
}

nems <- rownames(binned_hits)

binned_hits <- as.data.frame(binned_hits)



# Write the DF
write.csv(x = binned_hits, file = paste0("delims/",plasrefname, "_plasmid_coverage.csv"), row.names = TRUE)

write.csv(x = df6, file = paste0("delims/",plasrefname, "_plasmid_coverage_percentage.csv"), row.names = TRUE)

# Ignore this probs
tree1 <- ggtree(tree, branch.length = "none") %<+% metadata +
  geom_tippoint(aes(color = ColV, offset =100), size = .00000000001) +
  scale_color_manual(values = c("Yes" = "#df03fc", "No" = "white"))

              gheatmap(tree1, 
                      data = binned_hits,
                      font.size = 2,
                      hjust = 0,
                      colnames =FALSE,
                      width = 20,
                      offset = 0.1,
                      color = NULL) + 
  scale_fill_gradient(low = "white", high = "#8dd3c7", na.value = "white") +
  theme(legend.position = "right")
              

ggsave("FigureS10_pCERC4_map.tiff",
       figureS10, 
       path = "outputs/figures/", 
       device = "tiff", 
       width= 297, 
       height = 210, 
       unit ="mm", 
       dpi = "print")