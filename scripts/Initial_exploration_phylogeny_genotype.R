library(ggtree)
library(phytools)
library(tidyr)
library(readr)
library(magrittr)
library(dplyr)

#Load paths to files needed for the script
#tree_path <- "output/Salmonella_only/pangenome/roary_99/accessory_binary_genes.fa.newick"
tree_path <- "output/Salmonella_only/pangenome/iqtree/snp_sites/core_gene_alignment_snp_sites.aln.treefile"
abricate_path <- "output/Salmonella_only/abricate/genotype.txt"
pointfinder_path <- "output/Salmonella_only/pointfinder/Pointfinder.txt"
ColV_data_path <- "output/Salmonella_only/abricate/ColV_Liu.txt"
pMLST_data <- "output/Salmonella_only/pMLST/pMLST.txt"
output_name <- "Salmonella_only"
#refname <- "AVC171"
#cgMLST_path <- "metadata/curated_metadata_all.txt"
#phylotype_path <- "metadata/ST95_Phylotypes.txt"
#serotype_path <- "metadata/ST95_Serotypes.txt"
serovar_data <- "output/Salmonella_only/sistr/serovars.csv"


#Load in abricateR script to amalgamate data
source("scripts/abricateR2.R")

#Run abricateR
abricateR(
        file = abricate_path,
        output = output_name,
        identity = 90,
        length = 90,
        writecsv = TRUE,
        pointfinder_data = pointfinder_path,
        ColV_Liu_data = ColV_data_path,
        pMLST_data = pMLST_data
)

tree <- read.tree(tree_path)

df <- as.data.frame(tree$tip.label)

colnames(df) <- "name"

df <- left_join(df, Salmonella_only_simple_summary_N90L90)

df[is.na(df)] <- 0

df0 <- df

df <- df %>% select(name, starts_with("card_"), starts_with("pointfinder"))
#df <- df %>% select(name, starts_with("vfdb_"))

#colnames(df) <- gsub("EC_custom_VGI-([0-9]_)","EC_custom_VGI-0\\1", colnames(df))

#df <- df[rowSums(df[2:ncol(df)]) > 15,]

df <- df[,order(colnames(df))]

save_names <- df$name

df <- df %>% select(-name)

df2 <- df

df <- df[, colSums(df != 0) > 0]

df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

rownames(df) <- save_names
rownames(df2) <- save_names


colnames(df) <- gsub("Escherichia_coli", "E.coli", colnames(df))
colnames(df) <- gsub("Shigella_flexneri", "S.flexneri", colnames(df))
colnames(df) <- gsub("chloramphenicol_acetyltransferase", "chlor_acetyltrans", colnames(df))

library(readr)
library(magrittr)
library(dplyr)

#### cgMLST data ####

cgMLST <- read_delim("delims/Salmonella_cgMLST.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)


human_cgMLST<- cgMLST %>% filter(`Source Niche` == "Human")

gull_cgMLST<- cgMLST %>% filter(`Source Niche` == "Wild Animal")

human_cgMLST$HC200[human_cgMLST$HC200 %in% gull_cgMLST$HC200]

human_cgMLST$HC100[human_cgMLST$HC100 %in% gull_cgMLST$HC100]

human_cgMLST$HC50[human_cgMLST$HC50 %in% gull_cgMLST$HC50]

human_cgMLST$HC20[human_cgMLST$HC20 %in% gull_cgMLST$HC20]

human_cgMLST$HC10[human_cgMLST$HC10 %in% gull_cgMLST$HC10]

HC5_overlap <- human_cgMLST$HC5[human_cgMLST$HC5 %in% gull_cgMLST$HC5]

cgMLST %>% filter(HC5 %in% HC5_overlap) %>% select(Name, HC5) %>% arrange(desc(HC5))

HC2_overlap <- human_cgMLST$HC2[human_cgMLST$HC2 %in% gull_cgMLST$HC2]

cgMLST %>% filter(HC2 %in% HC2_overlap) %>% select(Name, HC2) %>% arrange(desc(HC2))

human_cgMLST$HC2[human_cgMLST$HC2 %in% gull_cgMLST$HC2]

human_cgMLST$`HC0 (indistinguishable)`[human_cgMLST$HC5 %in% gull_cgMLST$`HC0 (indistinguishable)`]



#pheatmap(df, fontsize_row = 2, fontsize_col = 2, cluster_cols = TRUE)

#### Read in and clean tree file ####
tree <-
        read.tree(file = tree_path)

#pheatmap(df, fontsize_row = 6, fontsize_col = 6)

cols_needed <- sort(unique(c(serovars_simple$strain_source, serovars_simple$serovar_cgmlst, 0,1,2 ))) %>% as.data.frame()
#cols_needed <- sort(unique(c(serovars_simple$strain_source, serovars_simple$serovar_cgmlst, 0,1 ))) %>% as.data.frame()

colnames(cols_needed) <- "variable"

source_cols <- as.data.frame(unique(serovars_simple$strain_source))
serovar_cols <- as.data.frame(unique(serovars_simple$serovar_cgmlst))
genotype_cols <- as.data.frame(c(0,1,2))

colnames(source_cols) <- "variable"
colnames(serovar_cols) <- "variable"
colnames(genotype_cols) <- "variable"

source_cols$colour <- gsub("Gull", "red", source_cols$variable) 
source_cols$colour <- gsub("Human", "blue", source_cols$colour) 

serovar_cols$colour <- rainbow(n = length(unique(serovar_cols$variable)))
genotype_cols$colour <- c("white","purple", "purple2")
#genotype_cols$colour <- c("white","red")

col_list <- rbind(source_cols, serovar_cols, genotype_cols)

cols_list <- left_join(cols_needed, col_list)

#replace yellow colour with something else - it is illegible
cols_list$colour <- gsub("#F0FF00","black", cols_list$colour)

p <- ggtree(tree) %<+%
        serovars_simple + 
        geom_tiplab(
                aes(colour = strain_source),
                size = 1,
                align = TRUE
                )+ 
        geom_tiplab(
                aes(label = serovar_cgmlst,
                    colour = serovar_cgmlst
                    ),
                align = TRUE,
                linetype = NULL,
                offset = 0.3,
                size = 1
) + theme(legend.position = "none") 

gheatmap(
        p = p,
        data = df,
        color = 'grey',
        width = 5,
        offset = 0.6,
        font.size = 1.5,
        colnames_position = "top",
        hjust = 0,
        #colnames = FALSE,
        colnames_offset_y = -0.1,
        colnames_angle = 90,
) + ggplot2::ylim(NA, 180) + theme(legend.position = "none") +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = c(cols_list$colour),
                na.value = 'grey'
        ) 

serovars_simple <- serovars_simple %>% rename("name" = "genome")

df3 <- df2 %>% select(-card_mdsB, -card_mdsC)

amr_counts <- as.data.frame(rowSums(df2))

colnames(amr_counts) <- "amr_count"

amr_counts$name <- rownames(amr_counts)

all_data <- left_join(amr_counts, serovars_simple)

all_data <- left_join(all_data, df0)

all_data %>% group_by(strain_source) %>% summarise(average = mean(amr_count))

all_data %>% group_by(strain_source, amr_count) %>% summarise(counts = n()) %>% arrange(desc(amr_count))

all_data %>% group_by(strain_source, amr_count) %>% summarise(counts = n()) %>% arrange(desc(amr_count))

df4 <- df2 %>% select(starts_with("card_cat"),
               starts_with("card_CTX"),
               starts_with("card_floR"),
               starts_with("card_MCR"),
               starts_with("card_IMP"),
               starts_with("card_OXA"),
               starts_with("card_Qnr"),
               starts_with("card_SHV"),
               starts_with("pointfinder")
)

fqr_snp <- df4 %>% select(starts_with("pointfinder")) %>% rowSums()

df4 <- df4 %>% select(-starts_with("pointfinder"))

df4$fqr_snp <- fqr_snp

df4$fqr_snp[df4$fqr_snp > 1] <- 1

CIA <- as.data.frame(rowSums(df4))

colnames(CIA) <- "CIA"

CIA$name <- rownames(CIA)

all_data <- left_join(all_data, CIA)

all_data %>% group_by(strain_source, CIA) %>% summarise(counts = n()) %>% arrange(desc(CIA))

all_data %>% group_by(plasmidfinder_IncHI2_1, CIA) %>% summarise(counts = n()) %>% arrange(desc(CIA))

df2$name <- rownames(df2)

all_data %>% filter(name %in% c("Seagull_18_154_S34", "Seagull_18_89_S8", "W_1_F11")) %>% select(name, starts_with("plasmidfinder_Inc")) %>% View()

gull1 <- df2 %>% filter(name %in% c("Seagull_18_154_S34"))

gull <- gull1 %>% select(-name)

i <- (colSums(gull, na.rm=T) != 0) # 

gsub(".*_","",colnames(gull[, i]))

gull2 <- df2  %>% filter(name %in% c("W_1_F11"))

gull <- gull2 %>% select(-name)

i <- (colSums(gull, na.rm=T) != 0) # 

gsub(".*_","",colnames(gull[, i]))

gull3  <- df2 %>% filter(name %in% c("Seagull_18_89_S8"))

gull <- gull3 %>% select(-name)

i <- (colSums(gull, na.rm=T) != 0) # 

gsub(".*_","",colnames(gull[, i]))



#Salmonella_all_genotype_serovars<- left_join(serovars, Salmonella_simple_summary_N90L90)

#Salmonella_all_genotype_serovars %>% select(serovar_cgmlst, serovar_data, starts_with("EC"))

#write_csv(Salmonella_all_genotype_serovars, "processed_data/Salmonella_genotype_serovar.csv")



