library(pheatmap)
library(ggplot)
library(magrittr)
library(dplyr)
library(reshape2)
library(readr)
library(OneR)
library(microbenchmark)
library(tidytree)
library(ggtree)
library(magrittr)
library(dplyr)
library(readr)
library(reshape2)
library(ComplexHeatmap)
library(ggplot2)
library(ggimage)

#### Path Definitions/Config ####
### IMPORTANT ###
### Check lines where there are hardcoded changes, such as line 21,
### if you plan to modify this script's inputs.

#Replace the variable below with the path to your SG17-135 repo
path_to_repo <-
        "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/AVC171/AVC171"

#Changes working directory
setwd(path_to_repo)

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

#Replace the name of the reference (if a reference based tree) with the trees ref
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

#Separate the abricate hits DF into sample-wise dfs
for(sample in unique(abricate_hits$name)){
        abricate_hits2 <- abricate_hits %>% filter(name == sample)
        #Create an empty matrix for each sample with a length = to that of the reference
        range_matrix <- rep(0, times = unique(abricate_hits$ref_length))
        #Separate out individual hits for a given sample
        #reassign value to 1 for given BP coordinates based on start and end coords of each hit
        for(hit in 1:nrow(abricate_hits2)){
                start_ <- abricate_hits2[[hit,4]]
                end_ <- abricate_hits2[[hit,5]]
                range_matrix[start_:end_] <- 1
        }
        #bind together the hit matrix for all samples
        base_matrix <- rbind(base_matrix, range_matrix)
}

base_matrix <- as.data.frame(base_matrix, stringsAsFactors = FALSE)

rownames(base_matrix) <- unique(abricate_hits$name)


x <- 1

df <- as.data.frame(base_matrix[1:2,])

df_length <- length(df)

(as.numeric(unique(abricate_hits$ref_length)) - (as.numeric(unique(abricate_hits$ref_length)) %% 100) )/ 100 

bin_ranges <- c(seq(from = 1, to = ref_length, by = 100), ref_length)

sum(df[,bin_ranges[1]:bin_ranges[2]])


x <- 1
y <- 2

for(j in 1:nrow(df)){
        for(i in df){
                #print(j)
                a <- bin_ranges[x]
                b <- bin_ranges[y]
                test <- sum(df[j,a:b])
                x <- x + 1
                y <- y + 1
        }
        print(test)
        x <- 1
        y <- 2
}

#        bins <- lapply(bins, as.numeric)
#        binsums <- lapply(bins, sum)
#        listy_ <- rbind(listy_, unlist(binsums))
#        print(paste("another one:",x))
#        x <- x + 1
#        
#}

#listy_ <- rbind(listy_, rep(x = 100, times = ncol(listy_)))
#listy_ <- listy_[1:676,]
rownames(listy_) <- c(list_2$V1)#,"AVC171")

df6 <- as.data.frame(rowSums(listy_))

df6$working_names <- rownames(df6)

colnames(df6) <- c("ColV_percent_hit","working_name")

df6$ColV_percent_hit <- round((df6$ColV_percent_hit/length_of_gene) * 100)

write.csv(x = listy_, file = paste0("delims/",plasrefname, "_plasmid_coverage.csv"), row.names = TRUE)

write.csv(x = df6, file = paste0("delims/",plasrefname, "_plasmid_coverage_percentage.csv"), row.names = TRUE)

#df6$working_name <- rownames(df6)



#df6$working_name

#pheatmap(listy_, cluster_cols = FALSE, fontsize_col = 1)
#
#
#abc <- length(list_)/3
#
#df <- data.frame(matrix(unlist(list_), nrow = length(unique(abricate_hits$name)), byrow=T), stringsAsFactors = F)
#
#colnames(df) <- c("name","GENE","Coverage_percentage")
#
#df$Coverage_percentage[is.na(df$Coverage_percentage)] <- 0
#
#df$Coverage_percentage <- as.numeric(df$Coverage_percentage)
#
#final_table <- dcast(df, name ~ GENE)
#
#final_final_table <- final_table[1:nrow(final_table),2:ncol(final_table)]
#
#final_final_table_2 <- final_final_table
#
#final_final_table[final_final_table < 60] <- 0
#final_final_table[final_final_table >= 60] <- 1
#
#rownames(final_final_table) <- final_table$name
#
#final_table <- final_final_table
#
##write.csv(final_table, "analysis/PAIs_present_absent.csv")
#
#pheatmap(final_final_table, fontsize_row = 2)
#
#