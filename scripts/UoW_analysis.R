library(ggtree)
library(phytools)
library(tidyr)
library(readr)
library(magrittr)
library(dplyr)


#Process enterobase data
if(!exists("metadata")){
        source("scripts/Enterobase_data_combine.R")
}

metadata <- enterobase_data %>% filter(Lab_Contact == "M. Cummins (UTS, Sydney)")

metadata$Name <- gsub("Ecoli", "E_coli", metadata$Name)

rm(assembly_stats)
rm(enterobase_data)

#Load paths to files needed for the script
#tree_path <- "output/Salmonella_only/pangenome/roary_99/accessory_binary_genes.fa.newick"
tree_path <- "output/new_output/all/accessory_tree_all.tree"
abricate_path <- "output/new_output/genotype.txt"
pointfinder_path <- "output/new_output/pointfinder.txt"
pMLST_data <- "output/new_output/pMLST.txt"

output_name <- "UoW_new"

output_dir <- "delims"


tree <- read.tree(tree_path)

tree <- midpoint.root(tree)

#Load in abricateR script to amalgamate data
source("scripts/abricateR2.R")

#Run abricateR
abricateR(
        abricate_in = abricate_path,
        output_directory = output_dir,
        output = output_name,
        identity = 90,
        length = 90,
        writecsv = FALSE,
        pointfinder_data = pointfinder_path,
        pMLST_data = pMLST_data
)

df <- UoW_new_simple_summary_N90L90

df <- df %>% rename(Name = "name")

df <- df %>% filter(Name %in% tree$tip.label)

metadata$working_name <- metadata$Assembly_Barcode

geno_meta <- left_join(metadata, df)

df_names <- df$Name

df <- df %>% select(-starts_with("Inc"), -Name)

df <- df[, colSums(df != 0) > 0]

df[sapply(df, is.integer)] <- lapply(df[sapply(df, is.integer)], 
                                     as.factor)
df$Name <- df_names

rownames(df) <- df$Name

geno_meta <- geno_meta %>% filter(Name %in% tree$tip.label)

geno_meta <- geno_meta %>% select(Name, everything())

rownames(geno_meta) <- geno_meta$Name

#### Colour definitions ####
genotype_vals <- as.factor(c(0, 1, 2, 3, 18, 29))

genotype_cols <- c("white","#fb8072", "blue","black", "red", "green")

names(genotype_cols) <- genotype_vals

Source_Niche_vals <- c("Human", "Wild Animal")

Source_Niche_cols<- c("blue","red")

names(Source_Niche_cols) <- Source_Niche_vals

var_cols <- c(Source_Niche_cols, genotype_cols)

p <- ggtree(tree) %<+%
        geno_meta +
        geom_tiplab(size = 1.5,
                    align = TRUE,
                    linesize = 0.15,
                    aes(color = Source_Niche)
        )

df_plot <- df %>% select(-Name,
                         starts_with("card"),
                         starts_with("colV_zoetis"),
                         starts_with("dfrA5_848"),
                         starts_with("EC_custom"),
                         starts_with("ecoh"),
                         starts_with("ISfinder_Feb_2020"),
                         starts_with("plasmidfinder"),
                         #starts_with("vfdb"),
                         starts_with("pointfinder")
)

heatmap <- gheatmap(
        p = p,
        data = df_plot,
        #colnames_offset_y = -0.1,
        font.size = 1.2,
        hjust = 0,
        colnames_position = "top",
        colnames = TRUE,
        colnames_angle = 90,
        #Accessory
        offset = 0.45,
        #CGA_snps
        #offset = 0.005,
        #Accessory
        width = 3.5,
        #CGA_snps
        #width = 3,
        color = "grey"
) +
        geom_tiplab(
                aes(label = IncF_RST),
                align = TRUE,
                linetype = NULL,
                #Accessory
                offset = 0.3,
                #CGA_snp_sites
                #offset = 0.006,
                size = 1.2
        ) +
        geom_tiplab(
                aes(label = paste0("ST",ST)),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -0.45,
                size = 1.2
        ) +
        geom_tiplab(
                aes(label = SISTR1_Serovar),
                align = TRUE,
                linetype = NULL,
                #Accessory
                #offset = 3.3,
                #CGA_snp_sites
                offset = -0.65,
                size = 1.2
        ) + scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = var_cols,
                na.value = 'grey'
        ) + ggplot2::ylim(NA, 170)#+

require(cowplot)
leg1 <- get_legend(heatmap)

plot(leg1)

heatmap + theme(legend.position="none")

heatmap

