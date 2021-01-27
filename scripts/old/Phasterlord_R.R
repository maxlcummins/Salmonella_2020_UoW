library(readr)
library(scales)
library(magrittr)
library(dplyr)

#### Reading in and preprocessing data ####
df <- read_delim("analysis/phagelord/ST95_all/phastaf_combined.txt", 
                               "\t", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)

colnames(df) <- c("scaffold", "start_coord","end_coord","gene","nucleotide_ID","strand","filename")

df <- df %>% filter(nucleotide_ID >= 90)

df$phage_ID <- gsub("-gi\\|.*","", df$gene)

df$filename <- gsub(".*ST95_all\\/","", df$filename)
df$filename <- gsub("\\/diamond.coords.bed","", df$filename)
df$accession <- gsub(".*_NC","NC", df$gene)
df$accession <- gsub("-gi.*","", df$accession)

headers <- read_delim("misc/prophage_headers_untrimmed.txt", 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE)

headers <- headers %>% arrange(X1) 

headers$X1 <- gsub(">","", headers$X1)

headers$phage_ID <- gsub("-gi\\|.*","", headers$X1)
headers$phage_name <- gsub(".*\\[","", headers$X1)
headers$phage_name <- gsub("\\].*","", headers$phage_name)
headers$accession <- gsub(".*_NC","NC", headers$X1)
headers$accession <- gsub("-gi.*","", headers$accession)


headers$protein_ID <- gsub(".*-gi\\|","", headers$X1)
headers$protein_ID <- gsub("\\|.*","", headers$protein_ID)

headers$protein_descriptor <- gsub(".*\\| ","", headers$X1)

seq_lengths <- read_delim("misc/taxid10239.tbl", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

seq_lengths_small <- seq_lengths %>% select(accession, length)



#### Scoring ####

proteins_per_phage <- headers %>% group_by(phage_ID) %>% summarise(counts = n())

protein_counts <- df %>% group_by(phage_ID, filename) %>% summarise(counts = n())

protein_counts <- left_join(protein_counts, proteins_per_phage, by = "phage_ID")

protein_counts$percentage <- percent(protein_counts$counts.x/protein_counts$counts.y, accuracy = 2)

protein_counts$percentage <- gsub("\\%","",protein_counts$percentage)

protein_counts$percentage <- as.numeric(protein_counts$percentage)




for(i in 1:nrow(protein_counts)){
        if(protein_counts[i,6] >= 100){
                protein_counts[i,7] <- 150
}else{
        protein_counts[i,7] <- "other"
}
        print(i)
}



