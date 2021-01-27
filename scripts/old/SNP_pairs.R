library(readr)

meta_path <- "processed_data/metagenodata.csv"
refname1 <- "AVC171"
refname2 <- "ESC_NA4045AA_AS"
refname3 <- "core_gene_alignment"

snp_dists_path1 <- "analysis/snippy/AVC171_HC50/snp_dists/AVC171_chromosome.pairwise_snps.csv"
snp_dists_path2 <- "analysis/snippy/AVC171_HC50/snp_dists/pAVC171-IncF.pairwise_snps.csv"
snp_dists_path3 <- "analysis/snippy/PH/snp_dists/ESC_NA4045AA_AS.pairwise_snps.csv"
snp_dists_path4 <- "analysis/pangenome/HC50_1106/output/core_gene_alignment_snp_dists.csv"


snp_name1 <- "AVC171_chr"
snp_name2 <- "pAVC171-F"
snp_name3 <- "ESC_NA4045AA_AS"
snp_name4 <- "core_gene_alignment"

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
snp_process(snp_dists_path3, snp_name3, refname2)
snp_process(snp_dists_path4, snp_name4, refname3)

snp_df <- full_join(AVC171_chr, `pAVC171-F`)
snp_df <- full_join(snp_df, ESC_NA4045AA_AS)
snp_df <- full_join(snp_df, core_gene_alignment)


meta <- read_csv(meta_path)

meta_small <- meta %>% select(working_name, Revised_Source_Niche, Country, Revised_Collection_Year, Classification)

snp_df <- left_join(snp_df, meta_small)

snp_df <- left_join(snp_df, meta_small, by = c("working_name2" = "working_name"))


colnames(snp_df) <- gsub("\\.x$", "1", colnames(snp_df))
colnames(snp_df) <- gsub("\\.y$", "2", colnames(snp_df))
colnames(snp_df) <- gsub("working_name$", "working_name1", colnames(snp_df))
colnames(snp_df) <- gsub("Revised_Source_Niche","Source", colnames(snp_df))
colnames(snp_df) <- gsub("Revised_Collection_Year","Year", colnames(snp_df))

paired <- snp_df %>% select(working_name1, Source1, Country1, Year1, Classification1,
                            AVC171_chr, core_gene_alignment, `pAVC171-F`, ESC_NA4045AA_AS, 
                            Classification2, Year2, Country2, Source2, working_name2)

paired <- paired %>% filter(working_name1 != working_name2)


human_non <- paired %>% filter(Source1 == "Human") %>%
        filter(Source2 != "Human")

