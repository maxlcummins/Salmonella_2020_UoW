library(ggtree)
library(phytools)
library(tidyr)
library(readr)
library(magrittr)
library(dplyr)

#Load paths to files needed for the script
#tree_path <- "analysis/snippy/AVC171_all/fasttree/AVC171_chromosome.clean.fullcore.tree"
tree_path <- "output/pangenomelord/initial_pangenome_207/core_gene_alignment.aln_snp_sites.treefile"
abricate_path <- "output/abricate/Salmonella/Salmonella_abricate.txt"  
#pointfinder_path <- "analysis/pointfinder/AVC171_all_pointfinder.txt"
ColV_data_path <- "output/abricate/Salmonella/colV_zoetis.txt"
#pMLST_data <- "analysis/pMLST/pMLST_results.txt"
output_name <- "Salmonella"
#refname <- "AVC171"
#cgMLST_path <- "metadata/curated_metadata_all.txt"
#phylotype_path <- "metadata/ST95_Phylotypes.txt"
#serotype_path <- "metadata/ST95_Serotypes.txt"
serovar_data <- "output/sistr/Salmonella_serovars.csv"


#Load in abricateR script to amalgamate data
source("scripts/abricateR2.R")

#Run abricateR
abricateR(
        file = abricate_path,
        output = output_name,
        identity = 90,
        length = 90,
        writecsv = TRUE,
        #pointfinder_data = pointfinder_path,
        #ColV_Liu_data = ColV_data_path,
        #pMLST_data = pMLST_data
)

serovars <- read_csv("output/sistr/Salmonella_serovars.csv")

serovars$strain_source <- serovars$genome

serovars$strain_source <- gsub("SIML.*|ISLHD.*","Human", serovars$strain_source)
serovars$strain_source <- gsub("SG.*|W_.*|Seagull.*","Gull", serovars$strain_source)

serovars_simple <- serovars %>% select(genome, strain_source, qc_status, cgmlst_subspecies, serovar_cgmlst)

df <- as.data.frame(tree$tip.label)

colnames(df) <- "name"

df <- left_join(df, Salmonella_simple_summary_N90L90)

df[is.na(df)] <- 0

df <- df %>% select(name, starts_with("card_"))

#colnames(df) <- gsub("EC_custom_VGI-([0-9]_)","EC_custom_VGI-0\\1", colnames(df))

#df <- df[rowSums(df[2:ncol(df)]) > 15,]

df <- df[,order(colnames(df))]

save_names <- df$name

df <- df %>% select(-name)

df2 <- df

df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

rownames(df) <- save_names
rownames(df2) <- save_names


colnames(df) <- gsub("Escherichia_coli", "E.coli", colnames(df))
colnames(df) <- gsub("Shigella_flexneri", "S.flexneri", colnames(df))
colnames(df) <- gsub("chloramphenicol_acetyltransferase", "chlor_acetyltrans", colnames(df))

#pheatmap(df, fontsize_row = 2, fontsize_col = 2, cluster_cols = TRUE)

#### Read in and clean tree file ####
tree <-
        read.tree(file = tree_path)

#pheatmap(df, fontsize_row = 6, fontsize_col = 6)

cols_needed <- sort(unique(c(serovars_simple$strain_source, serovars_simple$serovar_cgmlst, 0,1,2 ))) %>% as.data.frame()

colnames(cols_needed) <- "variable"

source_cols <- as.data.frame(unique(serovars_simple$strain_source) )
serovar_cols <- as.data.frame(unique(serovars_simple$serovar_cgmlst))
genotype_cols <- as.data.frame(c(0,1,2))

colnames(source_cols) <- "variable"
colnames(serovar_cols) <- "variable"
colnames(genotype_cols) <- "variable"

source_cols$colour <- gsub("Gull", "red", source_cols$variable) 
source_cols$colour <- gsub("Human", "blue", source_cols$colour) 

serovar_cols$colour <- rainbow(n = length(unique(serovar_cols$variable)))
genotype_cols$colour <- c("white","purple", "purple2")

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
) + ggplot2::ylim(NA, 250) + theme(legend.position = "none") +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = c(cols_list$colour),
                na.value = 'grey'
        ) 

nodepie(df, cols = 1, )


serovars <- serovars %>% rename("name" = genome)

#Salmonella_all_genotype_serovars<- left_join(serovars, Salmonella_simple_summary_N90L90)

#Salmonella_all_genotype_serovars %>% select(serovar_cgmlst, serovar_data, starts_with("EC"))

#write_csv(Salmonella_all_genotype_serovars, "processed_data/Salmonella_genotype_serovar.csv")



