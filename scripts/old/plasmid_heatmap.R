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

### IMPORTANT ###
### Check lines where there are hardcoded changes, such as line 21,
### if you plan to modify this script's inputs.

#Replace the variable below with the path to your SG17-135 repo
path_to_repo <-
        "/Users/maxcummins/Dropbox/Doctorate/Manuscripts/AVC171/AVC171"

#Changes working directory
setwd(path_to_repo)

path_to_tree <- "analysis/snippy/AVC171_HC50/fasttree/AVC171_chromosome.clean.fullcore.tree"
path_to_pais <- "analysis/abricate/AVC171_all/pAVC171-F.txt"
refname <- "AVC171"

#Read in the abricate genotype data sheet (small number of rows for colname reassignment)
pais_df <-
        read_delim(
                path_to_pais,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE,
                n_max = 10
        )

#Colname reassignment
colnames(pais_df)[c(1, 10:11)] <-
        c("name", "perc_coverage", "perc_identity")
pais_df_colnames <- colnames(pais_df)

#Re-read in PAI abricate genotype data sheet
pais_df <-
        read_delim(
                path_to_pais,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE,
                col_names = FALSE,
                skip = 1
        )

#Remove cases where there are multiple headers from concatenation of abricate reports
pais_df <- pais_df %>% filter(X2 != "SEQUENCE")


#Colname reassignment
colnames(pais_df) <- pais_df_colnames

#Convert percent coverage and identity to numeric type to allow filtering
pais_df$perc_coverage <- as.numeric(pais_df$perc_coverage)
pais_df$perc_identity <- as.numeric(pais_df$perc_identity)

#Filter to perc_identity > 95%
#pais_df <-
pais_df <- pais_df %>% filter(perc_identity > 90)

#Trim excess characters the assembly names and reassign this to rownames
pais_df$name <- gsub("\\..*", "", pais_df$name)

#Read in the tree file
tree <-
        read.tree(file = path_to_tree)

#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

tree$tip.label <- gsub("Reference", refname, tree$tip.label)

pais_df <- pais_df %>% filter(name %in% tree$tip.label)

#Replace "SAL_HC4750AA_AS" with SG17-135
#pais_df$name <- gsub("SAL_HC4750AA_AS", "SG17-135", pais_df$name)

pais_df$newcov <- gsub("\\/.*","", pais_df$COVERAGE)
pais_df$length_gene <- gsub(".*\\/","", pais_df$COVERAGE)

pais_df$newcov <- gsub("-",":", pais_df$newcov)

new_df <- pais_df %>% group_by(name, GENE, length_gene) %>% filter(perc_coverage > 0.5) %>% summarise(start =paste(sort(unique(newcov)), collapse=","), end = paste(sort(unique(newcov)), collapse=",")) #%>% filter(grepl("SPI-1_", GENE))
#new_df <- pais_df %>% group_by(name, GENE, length_gene) %>% summarise(start =paste(sort(unique(newcov)), collapse=","), end = paste(sort(unique(newcov)), collapse=",")) #%>% filter(grepl("SPI-1_", GENE))

new_df$end <- gsub("[0-9]+:","", new_df$end)
new_df$start <- gsub(":[0-9]+","", new_df$start)


spl <-strsplit(as.character(new_df$start), ",")
start_coord <- data.frame(name = new_df$name, gene = new_df$GENE,
                          length_gene = new_df$length_gene,
                          chunk1 = sapply(spl, "[", 1),
                          chunk2 = sapply(spl, "[", 2),
                          chunk3 = sapply(spl, "[", 3),
                          chunk4 = sapply(spl, "[", 4),
                          chunk5 = sapply(spl, "[", 5),
                          chunk6 = sapply(spl, "[", 6),
                          chunk7 = sapply(spl, "[", 7),
                          chunk8 = sapply(spl, "[", 8),
                          chunk9 = sapply(spl, "[", 9),
                          chunk10 = sapply(spl, "[", 10),
                          chunk11 = sapply(spl, "[", 11),
                          chunk12= sapply(spl, "[", 12),
                          chunk13 = sapply(spl, "[", 13),
                          chunk14 = sapply(spl, "[", 14),
                          chunk15 = sapply(spl, "[", 15),
                          chunk16 = sapply(spl, "[", 16),
                          chunk17 = sapply(spl, "[", 17),
                          chunk18 = sapply(spl, "[", 18),
                          chunk19 = sapply(spl, "[", 19),
                          chunk20= sapply(spl, "[", 20),
                          chunk21 = sapply(spl, "[", 21),
                          chunk22 = sapply(spl, "[", 22),
                          chunk23 = sapply(spl, "[", 23),
                          chunk24 = sapply(spl, "[", 24),
                          chunk25 = sapply(spl, "[", 25),
                          chunk26 = sapply(spl, "[", 26),
                          chunk27 = sapply(spl, "[", 27),
                          chunk28= sapply(spl, "[", 28),
                          chunk29 = sapply(spl, "[", 29),
                          chunk30 = sapply(spl, "[", 30),
                          chunk31 = sapply(spl, "[", 31),
                          chunk32 = sapply(spl, "[", 32),
                          chunk33 = sapply(spl, "[", 33),
                          chunk34 = sapply(spl, "[", 34),
                          chunk35 = sapply(spl, "[", 35),
                          chunk36= sapply(spl, "[", 36),
                          chunk37 = sapply(spl, "[", 37),
                          chunk38 = sapply(spl, "[", 38),
                          chunk39 = sapply(spl, "[", 39),
                          chunk40 = sapply(spl, "[", 40),
                          chunk41 = sapply(spl, "[", 41),
                          chunk42 = sapply(spl, "[", 42),
                          chunk43 = sapply(spl, "[", 43),
                          chunk44= sapply(spl, "[", 44),
                          chunk45 = sapply(spl, "[", 45),
                          chunk46 = sapply(spl, "[", 46),
                          chunk47 = sapply(spl, "[", 47),
                          chunk48 = sapply(spl, "[", 48),
                          chunk49 = sapply(spl, "[", 49),
                          chunk50= sapply(spl, "[", 50),
                          chunk51 = sapply(spl, "[", 51),
                          chunk52 = sapply(spl, "[", 52),
                          chunk53 = sapply(spl, "[", 53),
                          chunk54 = sapply(spl, "[", 54),
                          chunk55 = sapply(spl, "[", 55),
                          chunk56 = sapply(spl, "[", 56),
                          chunk57 = sapply(spl, "[", 57),
                          chunk58 = sapply(spl, "[", 58),
                          chunk59 = sapply(spl, "[", 59),
                          chunk60 = sapply(spl, "[", 60),
                          chunk61 = sapply(spl, "[", 61),
                          chunk62 = sapply(spl, "[", 62),
                          chunk53 = sapply(spl, "[", 53),
                          chunk64 = sapply(spl, "[", 64),
                          chunk65 = sapply(spl, "[", 65),
                          chunk66 = sapply(spl, "[", 66),
                          chunk67 = sapply(spl, "[", 67),
                          chunk68 = sapply(spl, "[", 68),
                          chunk69 = sapply(spl, "[", 69)
                          )

start_coord <- melt(start_coord, id=1:3, value.name = "start")

start_coord <- start_coord %>% select(-starts_with("variable")) 

spl <-strsplit(as.character(new_df$end), ",")
end_coord <- data.frame(name = new_df$name, gene = new_df$GENE,
                        length_gene = new_df$length_gene,
                        chunk1 = sapply(spl, "[", 1),
                        chunk2 = sapply(spl, "[", 2),
                        chunk3 = sapply(spl, "[", 3),
                        chunk4 = sapply(spl, "[", 4),
                        chunk5 = sapply(spl, "[", 5),
                        chunk6 = sapply(spl, "[", 6),
                        chunk7 = sapply(spl, "[", 7),
                        chunk8 = sapply(spl, "[", 8),
                        chunk9 = sapply(spl, "[", 9),
                        chunk10 = sapply(spl, "[", 10),
                        chunk11 = sapply(spl, "[", 11),
                        chunk12= sapply(spl, "[", 12),
                        chunk13 = sapply(spl, "[", 13),
                        chunk14 = sapply(spl, "[", 14),
                        chunk15 = sapply(spl, "[", 15),
                        chunk16 = sapply(spl, "[", 16),
                        chunk17 = sapply(spl, "[", 17),
                        chunk18 = sapply(spl, "[", 18),
                        chunk19 = sapply(spl, "[", 19),
                        chunk20= sapply(spl, "[", 20),
                        chunk21 = sapply(spl, "[", 21),
                        chunk22 = sapply(spl, "[", 22),
                        chunk23 = sapply(spl, "[", 23),
                        chunk24 = sapply(spl, "[", 24),
                        chunk25 = sapply(spl, "[", 25),
                        chunk26 = sapply(spl, "[", 26),
                        chunk27 = sapply(spl, "[", 27),
                        chunk28= sapply(spl, "[", 28),
                        chunk29 = sapply(spl, "[", 29),
                        chunk30 = sapply(spl, "[", 30),
                        chunk31 = sapply(spl, "[", 31),
                        chunk32 = sapply(spl, "[", 32),
                        chunk33 = sapply(spl, "[", 33),
                        chunk34 = sapply(spl, "[", 34),
                        chunk35 = sapply(spl, "[", 35),
                        chunk36= sapply(spl, "[", 36),
                        chunk37 = sapply(spl, "[", 37),
                        chunk38 = sapply(spl, "[", 38),
                        chunk39 = sapply(spl, "[", 39),
                        chunk40 = sapply(spl, "[", 40),
                        chunk41 = sapply(spl, "[", 41),
                        chunk42 = sapply(spl, "[", 42),
                        chunk43 = sapply(spl, "[", 43),
                        chunk44= sapply(spl, "[", 44),
                        chunk45 = sapply(spl, "[", 45),
                        chunk46 = sapply(spl, "[", 46),
                        chunk47 = sapply(spl, "[", 47),
                        chunk48 = sapply(spl, "[", 48),
                        chunk49 = sapply(spl, "[", 49),
                        chunk50= sapply(spl, "[", 50),
                        chunk51 = sapply(spl, "[", 51),
                        chunk52 = sapply(spl, "[", 52),
                        chunk53 = sapply(spl, "[", 53),
                        chunk54 = sapply(spl, "[", 54),
                        chunk55 = sapply(spl, "[", 55),
                        chunk56 = sapply(spl, "[", 56),
                        chunk57 = sapply(spl, "[", 57),
                        chunk58 = sapply(spl, "[", 58),
                        chunk59 = sapply(spl, "[", 59),
                        chunk60 = sapply(spl, "[", 60),
                        chunk61 = sapply(spl, "[", 61),
                        chunk62 = sapply(spl, "[", 62),
                        chunk53 = sapply(spl, "[", 53),
                        chunk64 = sapply(spl, "[", 64),
                        chunk65 = sapply(spl, "[", 65),
                        chunk66 = sapply(spl, "[", 66),
                        chunk67 = sapply(spl, "[", 67),
                        chunk68 = sapply(spl, "[", 68),
                        chunk69 = sapply(spl, "[", 69)
)

end_coord <- melt(end_coord, id=1:3, value.name = "end")

end_coord <- end_coord %>% select(-starts_with("variable")) 

coords <- start_coord

coords$end <- end_coord$end

coords <- coords[complete.cases(coords),]

unique(coords$length_gene)

coords$start <- as.numeric(coords$start)
coords$end <- as.numeric(coords$end)
coords$length_gene <- as.numeric(coords$length_gene)


coords$percentage <-  (((coords$end-coords$start)+1)/coords$length_gene)*100

test <- coords# %>% filter(name == "SAL_AB7542AA_AS", gene == "SPI-12_NC_006905_P4") %>% arrange(desc(end))

list_ <- NULL

counter <- 1

for(sample in unique(test$name)){
        test2 <- test %>% filter(name == sample)
        for(gene_ in unique(test$gene)){
                test3 <- test2 %>% filter(gene == gene_)
                length_of_gene <- test3$length_gene[1]
                if(is.na(length_of_gene) == FALSE){
                        range_matrix <- rep(0, times = length_of_gene)
                        for(hit in 1:nrow(test3)){
                                start_ <- test3[hit,4]
                                end_ <- test3[hit,5]
                                range_matrix[start_:end_] <- 1
                                range_matrix[range_matrix > 1] <- 1
                        }
                }
                newline <- c(sample, range_matrix)
                list_ <- rbind(list_,newline)
                print(paste("another loop", counter))
                counter <- counter + 1
        }
        
}

list_2 <- as.data.frame(list_, stringsAsFactors = FALSE)

rm(list_)

base_ <- as.data.frame(list_2[,1:2])

cols <- c(2:ncol(list_2))

base_$V2 <- apply(list_2[ , cols ] , 1 , paste, collapse = "")

listy_ <- NULL

x <- 1

for(i in 1:nrow(new_df)){
        d <- unlist(list_2[i,2:ncol(list_2)])
        bins <- split(d, ceiling(seq_along(d)/100))
        bins <- lapply(bins, as.numeric)
        binsums <- lapply(bins, sum)
        listy_ <- rbind(listy_, unlist(binsums))
        print(paste("another one:",x))
        x <- x + 1
        
}

#listy_ <- rbind(listy_, rep(x = 100, times = ncol(listy_)))
listy_ <- listy_[1:676,]
rownames(listy_) <- c(list_2$V1)#,"AVC171")

df6 <- as.data.frame(rowSums(listy_))

df6$working_names <- rownames(df6)

colnames(df6) <- c("ColV_percent_hit","working_name")

df6$ColV_percent_hit <- round((df6$ColV_percent_hit/144998) * 100)

#df6$working_name <- rownames(df6)



#df6$working_name

pheatmap(listy_, cluster_cols = FALSE, fontsize_col = 1)


abc <- length(list_)/3

df <- data.frame(matrix(unlist(list_), nrow = length(unique(pais_df$name)), byrow=T), stringsAsFactors = F)

colnames(df) <- c("name","GENE","Coverage_percentage")

df$Coverage_percentage[is.na(df$Coverage_percentage)] <- 0

df$Coverage_percentage <- as.numeric(df$Coverage_percentage)

final_table <- dcast(df, name ~ GENE)

final_final_table <- final_table[1:nrow(final_table),2:ncol(final_table)]

final_final_table_2 <- final_final_table

final_final_table[final_final_table < 60] <- 0
final_final_table[final_final_table >= 60] <- 1

rownames(final_final_table) <- final_table$name

final_table <- final_final_table

#write.csv(final_table, "analysis/PAIs_present_absent.csv")

pheatmap(final_final_table, fontsize_row = 2)

