library(readr)
library(magrittr)
library(dplyr)

meta_path <- "processed_data/metagenodata.csv"
refname1 <- "AVC171"
refname2 <- "ST95_all_core_gene_alignment"


snp_dists_path1 <- "analysis/snippy/AVC171_all/snp_dists/AVC171_chromosome.pairwise_snps.csv"
snp_dists_path2 <- "/Users/131785/core_gene_alignment_snp_dists.csv"


snp_name1 <- "AVC171_chr"
snp_name2 <- "ST95_all_core_gene_alignment"

#Read in cgMLST data
meta <- read_delim(cgMLST_path,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

colnames(meta) <- gsub("Assembly barcode", "working_name", colnames(meta))
colnames(meta) <- gsub(" ", "_", colnames(meta))


# This function cleans the file names and replaces "Reference" with the refname.
# It then melts the pairwise SNP table into a long list and generates a column
# we can later use to join multiple snp tables.
#
#It assigns these tables to the global environment based on the snp_name variable provided.
snp_process <- function(path_snp_data, snp_name, refname){
        snp_df <- read_csv(path_snp_data)
        colnames(snp_df)[1] <- "working_name"
        snp_df$working_name <- gsub(".out","",snp_df$working_name)
        snp_df$working_name <- gsub("Reference",refname, snp_df$working_name)
        colnames(snp_df) <- gsub(".out","",colnames(snp_df))
        colnames(snp_df) <- gsub("Reference", refname, colnames(snp_df))
        snp_df2 <- reshape2::melt(snp_df)
        colnames(snp_df2) <- c("working_name","working_name2", snp_name)
        snp_df2$pair <- paste(snp_df2$working_name, snp_df2$working_name2)
        assign(snp_name, snp_df2, envir = globalenv())
        
}

snp_process(snp_dists_path1, snp_name1, refname1)
snp_process(snp_dists_path2, snp_name2, refname1)

snp_df <- full_join(AVC171_chr, ST95_all_core_gene_alignment)

meta_small <- meta %>% select(working_name, Revised_Source_Niche, Country, Revised_Collection_Year, Classification)

snp_df <- left_join(snp_df, meta_small)

snp_df <- left_join(snp_df, meta_small, by = c("working_name2" = "working_name"))


colnames(snp_df) <- gsub("\\.x$", "1", colnames(snp_df))
colnames(snp_df) <- gsub("\\.y$", "2", colnames(snp_df))
colnames(snp_df) <- gsub("working_name$", "working_name1", colnames(snp_df))
colnames(snp_df) <- gsub("Revised_Source_Niche","Source", colnames(snp_df))
colnames(snp_df) <- gsub("Revised_Collection_Year","Year", colnames(snp_df))

cgMLST_data <- Metadata %>% select(working_name, starts_with("HC"), pSF_088_nores, pMLST)

colnames(cgMLST_data) <- gsub("ColV .*", "ColV", colnames(cgMLST_data))

snp_df <- left_join(snp_df, cgMLST_data, by = c("working_name1" = "working_name"))

snp_df <- left_join(snp_df, cgMLST_data, by = c("working_name2" = "working_name"))

colnames(snp_df) <- gsub("\\.x$", "_1", colnames(snp_df))
colnames(snp_df) <- gsub("\\.y$", "_2", colnames(snp_df))
colnames(snp_df) <- gsub("\\(indistinguishable\\)_", "", colnames(snp_df))

paired <- snp_df %>% select(HC400_1, HC200_1, HC100_1, HC50_1, HC20_1, HC10_1, HC5_1, HC2_1, HC0_1, pSF_088_nores_1, pMLST_1, working_name1, Source1, Country1, Year1, Classification1,
                            AVC171_chr, ST95_all_core_gene_alignment,
                            Classification2, Year2, Country2, Source2, working_name2, pMLST_2, pSF_088_nores_1, HC0_2, HC2_2, HC5_2, HC10_2, HC20_2, HC50_2, HC100_2, HC200_2, HC400_2)

paired <- paired %>% filter(working_name1 != working_name2)


human_non <- paired %>% filter(Source1 == "Human") %>%
        filter(Source2 != "Human") %>% filter(ST95_all_core_gene_alignment < 20)

human_non %>% group_by(Source2, HC20_1) %>% summarise(counts=n()) %>% View()
