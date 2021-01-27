library(ggtree)
library(phytools)
library(tidyr)

#Load paths to files needed for the script
#tree_path <- "analysis/snippy/AVC171_all/fasttree/AVC171_chromosome.clean.fullcore.tree"
tree_path <- "~/Desktop/core_gene_alignment_snp_sites.aln.tree"
abricate_path <- "analysis/abricate/ST95_all/Genotype.txt"  
pointfinder_path <- "analysis/pointfinder/AVC171_all_pointfinder.txt"
ColV_data_path <- "analysis/abricate/ST95_all/ColV_Liu.txt"
pMLST_data <- "analysis/pMLST/pMLST_results.txt"
output_name <- "ST95_all"
refname <- "AVC171"
cgMLST_path <- "metadata/curated_metadata_all.txt"
phylotype_path <- "metadata/ST95_Phylotypes.txt"
serotype_path <- "metadata/ST95_Serotypes.txt"


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

#### Read in and clean tree file ####
tree <-
        read.tree(file = tree_path)

#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (tree tip label)
tree$tip.label <- gsub("Reference", refname, tree$tip.label)

#### Read in and process cgMLST, O/H type and fimH type data ####
#Read in cgMLST data
Metadata <- read_delim(cgMLST_path,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

##Designate strains as HC5-4181 or Other
Metadata$HC50_or_other <-
        gsub("^1106$*", "HC50-1106", Metadata$HC50)
Metadata$HC50_or_other <-
        gsub("^[0-9].*", "Other", Metadata$HC50_or_other)

##Designate strains as HC200_Other
HC200_others <- Metadata$HC200 %>% table() %>% as.data.frame() %>% filter(Freq < 10)

#Designate strains as HC200_Other
Metadata$HC200_Other <- Metadata$HC200

Metadata$HC200_Other[Metadata$HC200_Other %in% HC200_others$.] <- "Other"

#Designate strains as HC100_Other
HC100_others <- Metadata$HC100 %>% table() %>% as.data.frame() %>% filter(Freq < 3)

Metadata$HC100_Other <- Metadata$HC100

Metadata$HC100_Other[Metadata$HC100_Other %in% HC100_others$.] <- "Other"

#Designate strains as HC50_other
HC50_others <- Metadata$HC50 %>% table() %>% as.data.frame() %>% filter(Freq <= 20)

#Designate strains as HC50_others
Metadata$HC50_Other <- Metadata$HC50

Metadata$HC50_Other[Metadata$HC50_Other %in% HC50_others$.] <- "Other"


#Read in serotype data
serotype <- read_delim(serotype_path,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

serotype <- serotype %>% select(Uberstrain, `O Antigen`, `H Antigen`)

#Read in fimH data
phylotype <- read_delim(phylotype_path,
                        "\t",
                        escape_double = FALSE,
                        trim_ws = TRUE)

phylotype <- phylotype %>% select(Uberstrain, `fimH (fimTyper)`)

Metadata <- left_join(Metadata, phylotype)

Metadata <- left_join(Metadata, serotype)

colnames(Metadata) <- gsub("fimH \\(fimTyper\\)", "fimH_type", colnames(Metadata))

colnames(Metadata) <- gsub(" Antigen", "_type", colnames(Metadata))

Metadata$O_type <- gsub("O50 or O2|O2 or O50","O2/O50", Metadata$O_type)

Metadata$O_type <- gsub("uncertain","O?", Metadata$O_type)

Metadata$O_type <- gsub("-","O-", Metadata$O_type)

Metadata$H_type <- gsub("uncertain","H?", Metadata$H_type)

Metadata$H_type <- gsub("-","H-", Metadata$H_type)

Metadata$fimH_type <- gsub("\\*","fimH?", Metadata$fimH_type)


#### General manipulation of Metadata ####

#Add new column working_name
Metadata$working_name <- Metadata$`Assembly barcode`

#replace the assembly barcode with strain name
Metadata$working_name <-
        gsub("ESC_SA8243AA_AS", refname, Metadata$working_name)

#Filter metadata table to only contain strains from the tree
Metadata <- Metadata %>% filter(working_name %in% tree$tip.label)

#Change the column order to put working_name first
Metadata <- Metadata %>% select(working_name, everything())

#set rowname to be equal to working name
rownames(Metadata) <- Metadata$working_name

#Remove spaces from column names for metadata - this usually causes issues
colnames(Metadata) <- gsub(" ", "_", colnames(Metadata))

#### Adding in of URLs for icons ####
#Generate a new column called Flag within which we will give URLs based on
#countries of origin, used by ggimage to render flags in our image.
Metadata$Flag <- Metadata$Country
Metadata$Flag <-
        gsub(
                "Australia",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Australia.jpg",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Denmark",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Denmark.jpg",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "United States",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/USA.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "^Ireland$",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Ireland.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Vietnam",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Vietnam.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "United Kingdom",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/United_Kindgom.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Mexico",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Mexico.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Ghana",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Ghana.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "United Arab Emirates",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/UAE.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Scotland",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Scotland.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Germany",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Germany.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Northern Ireland",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/Northern_Ireland.jpg",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Japan",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/japan.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Sweden",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/sweden.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Norway",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/norway.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "China",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/china.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Nepal",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/nepal.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Saudi Arabia",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/saudi-arabia.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Netherlands",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/netherlands.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "New Zealand",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/new-zealand.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Hungary",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/hungary.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "France",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/france.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Singapore",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/singapore.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Croatia",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/croatia.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Sri Lanka",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/sri-lanka.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Finland",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/finland.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Canada",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/canada.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "India",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/india.png",
                Metadata$Flag
        )
Metadata$Flag <-
        gsub(
                "Italy",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Flags/italy.png",
                Metadata$Flag
        )
#replace countries of origin that are unknown with a URL for a question mark
Metadata$Flag[is.na(Metadata$Flag)] <-
        "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/question.png"

#Generate a new column called Classification within which we will give URLs based on
#types of strain source, used by ggimage to render icons in our image.
Metadata$Classification_img <- Metadata$Classification
Metadata$Classification_img <-
        gsub(
                "SEPEC",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/blood.png",
                Metadata$Classification_img
        )
Metadata$Classification_img[is.na(Metadata$Classification_img)] <-
        "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/question.png"
Metadata$Classification_img <-
        gsub(
                "Faecal",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/colon.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "RMAE",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/meat.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "APEC",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/poultry.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "UPEC",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/urine.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "ExPEC",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/ExPEC.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "Environmental",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/earth.png",
                Metadata$Classification_img
        )
Metadata$Classification_img <-
        gsub(
                "Other",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/other.png",
                Metadata$Classification_img
        )

#Generate a new column called Pathogen within which we will use to give URLs based on
#pathotype of strains, used by ggimage to render icons in our image.
Metadata$Pathogen <- Metadata$Classification
Metadata$Pathogen <- gsub("ExPEC", "Systemic", Metadata$Pathogen)
Metadata$Pathogen <- gsub("APEC", "Systemic", Metadata$Pathogen)
Metadata$Pathogen <- gsub("UPEC", "Urine", Metadata$Pathogen)
Metadata$Pathogen <- gsub("SEPEC", "Systemic", Metadata$Pathogen)
Metadata$Pathogen <- gsub("RMAE", "Raw Chicken", Metadata$Pathogen)
Metadata$Pathogen <- gsub("Faecal", "Flora", Metadata$Pathogen)

#Here are the URLs for Pathotype images
Metadata$Pathogen_img <- Metadata$Pathogen
Metadata$Pathogen_img <-
        gsub(
                "Flora",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/colon.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Raw Chicken",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/raw_chicken2.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Systemic",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/blood.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Urine",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/urine.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Environmental",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/earth.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img <-
        gsub(
                "Other",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/other.png",
                Metadata$Pathogen_img
        )
Metadata$Pathogen_img[is.na(Metadata$Pathogen_img)] <-
        "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/question.png"

#Same as above for source data
Metadata$Revised_Source_Niche_img <- Metadata$Revised_Source_Niche
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Canine",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/dog.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Poultry",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/poultry.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Human",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/human.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Bovine",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/cow.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Environment",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/earth.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img <-
        gsub(
                "Other",
                "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/other.png",
                Metadata$Revised_Source_Niche_img
        )
Metadata$Revised_Source_Niche_img[is.na(Metadata$Revised_Source_Niche_img)] <-
        "https://raw.githubusercontent.com/maxlcummins/AVC171/master/images/Icons/question.png"


#### Processing of Genotypic data ####
#Filter the genotypic data to only contain strains in our tree
ST95_all_simple_summary_N90L90 <-
        ST95_all_simple_summary_N90L90 %>% filter(name %in% tree$tip.label)

#Rename our genotypic dataframe for brevity
df <- ST95_all_simple_summary_N90L90

#Reorder columns so that name is removed, ColV carriage (Liu) comes first, then every other column follows
df <- df %>% select(ColV, everything(),-name)

#Change ColV column to character
df$ColV <- as.character(df$ColV)

#Replace multiple hits for a given gene with a 1
df[df > 1] <- 1

#create a data frame for generating colsums for given genes (to determine prevelance of a given gene/trait)
colsum <- cbind(colnames(df), colSums(df)) %>% as.data.frame()
colsum$V2 <- as.numeric(colsum$V2)
df <- data.frame(lapply(df, as.numeric), stringsAsFactors = FALSE)
df <- df[, colSums(df != 0) > 0]

#Remove unwanted AMR genes from card
df <- df %>% select(
        -starts_with("card_acr"),
        -starts_with("card_bac"),
        -starts_with("card_bae"),
        -starts_with("card_cpxA"),
        -starts_with("card_CRP"),
        -starts_with("card_emr"),
        -starts_with("card_eptA"),
        -starts_with("card_Escherichia_coli_acrA"),
        -starts_with("card_Escherichia_coli_amp"),
        -starts_with("card_Escherichia_coli_emrE"),
        -starts_with("card_Escherichia_coli_mdfA"),
        -starts_with("card_evg"),
        -starts_with("card_gad"),
        -starts_with("card_H.NS"),
        -starts_with("card_kdpE"),
        -starts_with("card_marA"),
        -starts_with("card_mdt"),
        -starts_with("card_msbA"),
        -starts_with("card_pmrF"),
        -starts_with("card_tolC"),
        -starts_with("card_ugd"),
        -starts_with("card_yoj"),
        -starts_with("card_qacH"),
        -starts_with("card_mphA"),
        -starts_with("card_mphC"),
        -starts_with("card_linG"),
        -starts_with("card_determinant_of_bleomycin_resistance"),
        -starts_with("card_ErmB"),
        -starts_with("card_FosA3"),
        -starts_with("card_msrA"),
        -starts_with("card_SAT")
)

#Remove unwanted virulence genes - these are in our custom DB or members of operons we dont want to count additively
df <- df %>% select(
        -starts_with("EC_custom_malX"),
        -matches("EC_custom_eit[B-D]"),
        -starts_with("EC_custom_VGI"),
        -starts_with("EC_custom_malX"),
        -starts_with("EC_custom_yeeT"),
        -starts_with("EC_custom_iucD"),
        -starts_with("EC_custom_irp"),
        -starts_with("EC_custom_fim"),
        -starts_with("EC_custom_fyuA")
)

#Remove IS elements
df <- df %>% select(-starts_with("ISfinder_Feb_2020"))

#Remove unwanted virulence genes - these are in our custom DB or members of operons we dont want to count additively
df <- df %>% select(
        -starts_with("vfdb_chu"),
        -starts_with("vfdb_ent"),
        -starts_with("vfdb_chu"),
        -starts_with("vfdb_fepB|C|D"),
        -starts_with("vfdb_chu"),
        -starts_with("vfdb_fimA|B|C|D|E|F|G|I"),
        -starts_with("vfdb_gsp[D-M]"),
        -starts_with("vfdb_iro[B-D]"),
        -starts_with("vfdb_iuc[ABC]"),
        -starts_with("vfdb_iro[B-D]"),
        -starts_with("vfdb_pap[DEFHIJKX]"),
        -starts_with("vfdb_sfa[ABCDEFGHXY]"),
        -starts_with("vfdb_yag[WXYZ]"),
        -starts_with("vfdb_ybt[EPQSTUX]"),
        -starts_with("vfdb_ykgK")
)

#Remove unwanted virulence genes - these are in our custom DB or members of operons we dont want to count additively
df <- df %>% select(
        -matches("vfdb_chu"),
        -matches("vfdb_ent"),
        -matches("vfdb_chu"),
        -matches("vfdb_fep[BCD]"),
        -matches("vfdb_chu"),
        -matches("vfdb_fim[ABCDEFGI]"),
        -matches("vfdb_gsp[D-M]"),
        -matches("vfdb_iro[B-D]"),
        -matches("vfdb_iuc[ABC]"),
        -matches("vfdb_iro[B-D]"),
        -matches("vfdb_pap[BDEFHIJKX]"),
        -matches("vfdb_sfa[ABCDEFGHXY]"),
        -matches("vfdb_yag[WXYZ]"),
        -matches("vfdb_ybt[EPQSTUX]"),
        -matches("vfdb_ykgK")
)

#Save a copy of this df we use later
db_df <- df

#Change colnames to reflect gene function rather than database source
colnames(df) <- gsub("EC_custom_merA.*", "res_merA", colnames(df))
colnames(df) <- gsub("EC_custom_terA.*", "res_terA", colnames(df))
colnames(df) <- gsub("EC_custom_intI1.*", "res_intI1", colnames(df))
colnames(df) <- gsub("EC_custom_", "vir_", colnames(df))
colnames(df) <- gsub("card_", "res_", colnames(df))
colnames(df) <- gsub("vfdb_", "vir_", colnames(df))
colnames(df) <- gsub("plasmidfinder_", "plas_", colnames(df))

#creates a df to be used in generation of new colnames with gene sums
old_new_colnames <- rbind(colnames(df),colnames(df))

#Sum each hit for each gene
genesums <- as.data.frame(colSums(data.matrix(df)))

#rename genesum column to 'sum'
colnames(genesums) <- 'sum'

#paste together the new colnames and assign to our df with old and new names
genesums$rowsumcat <- paste(rownames(genesums), " (", genesums$sum , "/", nrow(df),")", sep="")

#duplicates our genotype table for a second figure containing gene sums
genotype_wsums <- df

#assigns new colnames to hit_table4
colnames(genotype_wsums) <- genesums$rowsumcat


df <- genotype_wsums

#Separate out our gene hits from different databases to allow us to change their cell contents for dataviz
r <- df %>% select(starts_with("res_"))
p <- df %>% select(starts_with("plas_"))
cv <- df %>% select(starts_with("ColV"))
v <- df %>% select(starts_with("vir_"))

#Set plasmid gene hits to be 2, virulence gene hits to be 3, colv genes to be 4
p[p == 1] <- 2
v[v == 1] <- 3
cv[cv == 1] <- 4

#Bind these columns back togeteer
df <- cbind(cv, r, p, v)

#replace multiple columns for sitA hits to be a column for a single sitA hits
df$vir_sitA <- df %>% select(starts_with("vir_sitA")) %>% rowSums()
df$vir_sitA[df$vir_sitA > 1] <- 3
df <- df %>% select(-starts_with("vir_sitA_"))

#replace multiple columns for incBOKZ hits to be a column for a single incBOKZ hit
df$plas_IncB.O.K.Z <-
        df %>% select(starts_with("plas_IncB")) %>% rowSums()
df$plas_IncB.O.K.Z[df$plas_IncB.O.K.Z > 1] <- 2
df <- df %>% select(-starts_with("plas_IncB.O.K.Z_"))

#Order the genes alphabetically (as they are prefixed by DB name they will be sorted alphabetically by DB name then by gene name)
df <- df[, order(names(df))]

#Trim database name from column names
colnames_save <- gsub("(res|vir|plas)_", "", colnames(df))

df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

colnames(df) <- colnames_save

rownames(df) <- ST95_all_simple_summary_N90L90$name

#### Reading and processing of plasmid mapping data ####

#Read in plasmid percent coveraage dataframe
plas_perc_cov <- read_csv("delims/pBCE049_1_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,2:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("pBCE049_1", "working_name")

Metadata <- left_join(Metadata, plas_perc_cov)

plas_perc_cov <- read_csv("delims/pUTI89_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,2:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("pUTI89", "working_name")

Metadata <- left_join(Metadata, plas_perc_cov)

plas_perc_cov <- read_csv("delims/pSF_088_nores_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,2:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("pSF_088_nores", "working_name")

Metadata <- left_join(Metadata, plas_perc_cov)

plas_perc_cov <- read_csv("delims/pAPEC_O2_ColV_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,2:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("pAPEC_O2_ColV", "working_name")

Metadata <- left_join(Metadata, plas_perc_cov)

plas_perc_cov <- read_csv("delims/pU1_51_B10_plasmid_coverage_percentage.csv")

plas_perc_cov <- plas_perc_cov[,2:ncol(plas_perc_cov)]

colnames(plas_perc_cov) <- c("pU1_51_B10", "working_name")

Metadata <- left_join(Metadata, plas_perc_cov)


#### Processing of pMLST data ####

pMLST_ <- ST95_all_simple_summary_N90L90 %>% select(name, pMLST)

pMLST_$pMLST <- gsub("^0$","", pMLST_$pMLST)
pMLST_$pMLST <- gsub("[A-Z]-:","", pMLST_$pMLST)
pMLST_$pMLST <- gsub(":K-","", pMLST_$pMLST)

colnames(pMLST_) <- c("working_name", "pMLST")

Metadata <- left_join(Metadata, pMLST_)

df_join <- df

df_join$working_name <- rownames(df_join)

Metadata <- left_join(Metadata, df_join)

metadata <- left_join(Metadata, ST95_all_simple_summary_N90L90, by = c("working_name" = "name"))

#### Designations of strains as APEC or ExPEC ####

APEC <- df %>% select(`iutA_KU578032.1 (330/668)`, `hylF_CP000836.1 (336/668)`, `iss_type1_pAPEC.O2.ColV (4/668)`, `iroN (388/668)`, `ompT_episomal_HM210637.1 (324/668)`)

APEC <- sapply(APEC, as.numeric) %>% as.data.frame()

APEC[APEC > 1] <- 1

APEC$APEC_val <- rowSums(APEC)

rownames(APEC) <- rownames(df)

APEC$working_name <- rownames(APEC)

APEC <- APEC %>% select(working_name, APEC_val)

APEC$APEC_val[APEC$APEC_val < 3] <- 0
APEC$APEC_val[APEC$APEC_val >= 3] <- 1

Metadata <- left_join(Metadata,APEC)


ExPEC <- df %>% select(`papC (585/668)`, `papA_CU928161 (504/668)`, `sfaS (101/668)`, `focG (1/668)`, starts_with("afaA"), starts_with("dra"), `kpsMT.II._K2.CP000468.1..APEC.O1.3376659.3378102 (655/668)`, `iutA_KU578032.1 (330/668)`)

ExPEC <- sapply(ExPEC, as.numeric) %>% as.data.frame()

ExPEC[ExPEC > 1] <- 1

ExPEC$V1 <- ExPEC$`papC (585/668)` + ExPEC$`papA_CU928161 (504/668)`

ExPEC$V2 <- ExPEC$`sfaS (101/668)` + ExPEC$`focG (1/668)`

ExPEC$V4 <- ExPEC$`kpsMT.II._K2.CP000468.1..APEC.O1.3376659.3378102 (655/668)`

ExPEC$V5 <- ExPEC$`iutA_KU578032.1 (330/668)`

ExPEC[ExPEC > 1] <- 1

ExPEC <- ExPEC %>% select(starts_with("V"))

ExPEC$ExPEC_val <- rowSums(ExPEC)

rownames(ExPEC) <- rownames(df)

ExPEC$working_name <- rownames(ExPEC)

ExPEC <- ExPEC %>% select(working_name, ExPEC_val)

ExPEC$ExPEC_val[ExPEC$ExPEC_val < 2] <- 0
ExPEC$ExPEC_val[ExPEC$ExPEC_val >= 2] <- 1

Metadata <- left_join(Metadata,ExPEC)

#### Counting of AMR and non-IncF repA genes (Inc only) ####

amr_df <- db_df %>% select(starts_with("card_"))

amr_count <- as.data.frame(rowSums(amr_df))

amr_count$working_name <- ST95_all_simple_summary_N90L90$name

colnames(amr_count) <- c("amr_counts","working_name")

Metadata <- left_join(Metadata, amr_count)


amr_gene_phenotype <- ST95_all.N90.L90.PASS %>% select(GENE, RESISTANCE)

amr_gene_phenotype <- amr_gene_phenotype %>% unique()

amr_gene_phenotype <- amr_gene_phenotype %>% filter(grepl("card", GENE))

amr_gene_phenotype$GENE <- gsub("\\(|\\)", ".", amr_gene_phenotype$GENE)
amr_gene_phenotype$GENE <- gsub("\\'", ".", amr_gene_phenotype$GENE)
amr_gene_phenotype$GENE <- gsub("\\-", ".", amr_gene_phenotype$GENE)

#Create a column for whether or not a strain is resistant to aminoglycosides
aminoglycosides <- amr_gene_phenotype %>% filter(grepl("aminoglycoside",RESISTANCE)) %>% select(GENE)
aminoglycosides <- aminoglycosides[[1]]
aminoglycosides <- aminoglycosides[aminoglycosides %in% colnames(amr_df)]

aminoglycoside <- amr_df %>% select(aminoglycosides)

aminoglycoside_res <- rowSums(aminoglycoside)

#Create a column for whether or not a strain is resistant to trimethoprim
trimethoprims <- amr_gene_phenotype %>% filter(grepl("diaminopyrimidine",RESISTANCE)) %>% select(GENE)
trimethoprims <- trimethoprims[[1]]
trimethoprims <- trimethoprims[trimethoprims %in% colnames(amr_df)]

trimethoprim <- amr_df %>% select(trimethoprims) %>% select(starts_with("card_dfrA"))

trimethoprim_res <- rowSums(trimethoprim)

#Create a column for whether or not a strain is resistant to fluoroquinolones
fluoroquinolones <- amr_gene_phenotype %>% filter(grepl("fluoroquinolone",RESISTANCE)) %>% select(GENE)
fluoroquinolones <- fluoroquinolones[[1]]
fluoroquinolones <- fluoroquinolones[fluoroquinolones %in% colnames(amr_df)]

fluoroquinolone <- amr_df %>% select(fluoroquinolones) %>% select(starts_with("card_Qnr"))

fluoroquinolone_snps <- ST95_all_simple_summary_N90L90 %>% select(starts_with("pointfinder_gyrA"), starts_with("pointfinder_parC"))

fluoroquinolone <- cbind(fluoroquinolone, fluoroquinolone_snps)

fluoroquinolone_res <- rowSums(fluoroquinolone)

#Create a column for whether or not a strain is resistant to sulfonamides
sulfonamides <- amr_gene_phenotype %>% filter(grepl("sulfonamide",RESISTANCE)) %>% select(GENE)
sulfonamides <- sulfonamides[[1]]
sulfonamides <- sulfonamides[sulfonamides %in% colnames(amr_df)]

sulfonamide <- amr_df %>% select(sulfonamides)

sulfonamide_snps <- ST95_all_simple_summary_N90L90 %>% select(starts_with("pointfinder_folP"))

sulfonamide <- cbind(sulfonamide, sulfonamide_snps)

sulfonamide_res <- rowSums(sulfonamide)

sul_tri <- cbind(sulfonamide_res, trimethoprim_res)

sul_tri_res <- rowSums(sulfonamide)



#Create a column for whether or not a strain is resistant to tetracyclines
tetracyclines <- amr_gene_phenotype %>% filter(grepl("tetracycline",RESISTANCE)) %>% select(GENE)
tetracyclines <- tetracyclines[[1]]
tetracyclines <- tetracyclines[tetracyclines %in% colnames(amr_df)]

tetracycline <- amr_df %>% select(tetracyclines) %>% select(starts_with("card_tet"))

tetracycline_res <- rowSums(tetracycline)

#Create a column for whether or not a strain is resistant to glycopeptides
glycopeptides <- amr_gene_phenotype %>% filter(grepl("peptide",RESISTANCE)) %>% select(GENE)
glycopeptides <- glycopeptides[[1]]
glycopeptides <- glycopeptides[glycopeptides %in% colnames(amr_df)]

glycopeptide <- amr_df %>% select(glycopeptides) %>% select(starts_with("card_MCR"))

glycopeptide_res <- rowSums(glycopeptide)

#Create a column for whether or not a strain is resistant to amphenicol
amphenicols <- amr_gene_phenotype %>% filter(grepl("phenicol",RESISTANCE)) %>% select(GENE)
amphenicols <- amphenicols[[1]]
amphenicols <- amphenicols[amphenicols %in% colnames(amr_df)]

amphenicol <- amr_df %>% select(amphenicols) %>% select(-starts_with("card_tet"), -starts_with("card_msrA"))

amphenicol_res <- rowSums(amphenicol)

#Create a column for whether or not a strain is resistant to third_gen_cephalosporin
beta_lactams <- amr_gene_phenotype %>% filter(grepl("ceph|penam|penem",RESISTANCE)) %>% select(GENE)
beta_lactams <- beta_lactams[[1]]
beta_lactams <- beta_lactams[beta_lactams %in% colnames(amr_df)]

beta_lactamase <- amr_df %>% select(beta_lactams) %>% select(-starts_with("card_tet"),
                                                            -starts_with("card_cmlA1"),
                                                            -starts_with("card_floR"))

ESBL <- amr_df %>% select(beta_lactams) %>% select(-starts_with("card_tet"),
                                                               -starts_with("card_cmlA1"),
                                                               -starts_with("card_floR"),
                                                               -starts_with("card_TEM"),
                                                               -starts_with("card_SHV"),
                                                               -starts_with("card_OXA"),
                                                               -starts_with("card_DHA"))


beta_lactam_res <- rowSums(beta_lactamase)
ESBL_res <- rowSums(ESBL)

ESBL_res[ESBL_res > 1] <- 1

ESBL_res <- as.data.frame(ESBL_res)

ESBL_res$working_name <- ST95_all_simple_summary_N90L90$name

class_res <- cbind(fluoroquinolone_res, aminoglycoside_res, sul_tri_res, tetracycline_res, glycopeptide_res, trimethoprim_res, amphenicol_res, beta_lactam_res)

class_res[class_res > 1] <- 1

class_res_count <- rowSums(class_res)

class_res <- as.data.frame(class_res)

class_res$working_name <- ST95_all_simple_summary_N90L90$name

class_res_count <- as.data.frame(class_res_count)

class_res_count$working_name <- ST95_all_simple_summary_N90L90$name

colnames(class_res_count)[1] <- "class_res_counts"
colnames(ESBL_res)[1] <- "ESBL_ress"

Metadata <- Metadata %>% left_join(class_res_count)
Metadata <- Metadata %>% left_join(class_res)
Metadata <- Metadata %>% left_join(ESBL_res)

Metadata$class_res_counts <- gsub("^","classes_res_0",Metadata$class_res_counts)
Metadata$ESBL_ress <- gsub("1","ESBL_pos",Metadata$ESBL_ress)
Metadata$ESBL_ress <- gsub("0","ESBL_neg",Metadata$ESBL_ress)

plas_df <- db_df %>% select(starts_with("plasmidfinder_Inc")) %>% select(-starts_with("plasmidfinder_IncF"))

plas_count <- as.data.frame(rowSums(plas_df))

plas_count$working_name <- ST95_all_simple_summary_N90L90$name

colnames(plas_count) <- c("plas_counts","working_name")

Metadata <- left_join(Metadata, plas_count)

Metadata$amr_counts <- gsub("^","count_res_0",Metadata$amr_counts)
Metadata$amr_counts <- gsub("count_res_010","count_res_10",Metadata$amr_counts)
Metadata$amr_counts <- gsub("count_res_011","count_res_11",Metadata$amr_counts)
Metadata$amr_counts <- gsub("count_res_012","count_res_12",Metadata$amr_counts)
Metadata$amr_counts <- gsub("count_res_013","count_res_13",Metadata$amr_counts)


Metadata$plas_counts <- gsub("^","count_plas_",Metadata$plas_counts)


#### Generation of small dataframe for Figure 1 ####

df$`ColV (328/668)`<- gsub("4","ColV_pos",df$`ColV (328/668)`)
df$`ColV (328/668)` <- gsub("0","ColV_neg",df$`ColV (328/668)`)

df_small <- as.data.frame(df$`ColV (328/668)`)
colnames(df_small) <- gsub("^df.*", "ColV", colnames(df_small))
df_small$working_name <- rownames(df)
df_small <- left_join(df_small, Metadata)
#df_small$Revised_Source_Niche <- replace_na(df_small$Revised_Source_Niche, "Unknown")
df_small$Revised_Source_Niche <- gsub("^", "Source_", df_small$Revised_Source_Niche)
df_small <- as.data.frame(df_small)
colnames(df_small) <- gsub("intI1.*", "intI1", colnames(df_small))
df_small <- df_small %>% select(Revised_Source_Niche, pU1_51_B10, pBCE049_1, pAPEC_O2_ColV, pSF_088_nores, pUTI89, intI1, class_res_counts, ESBL_ress, plas_counts)

for(i in 1:nrow(df_small)){
        rowmax <- max(df_small[i,2:6])
        print(rowmax)
        if(df_small[i,2] < rowmax){
                df_small[i,2] <- 0
        }
        if(df_small[i,3] < rowmax){
                df_small[i,3] <- 0
        }
        if(df_small[i,4] < rowmax){
                df_small[i,4] <- 0
        }
        if(df_small[i,5] < rowmax){
                df_small[i,5] <- 0
        }
        if(df_small[i,6] < rowmax){
                df_small[i,6] <- 0
        }
}

df_small2 <- df_small

#### Binning of percentage hits for plasmids into ranges #### 
df_small$pBCE049_1 <- paste0("pBCE049_1_",
                             as.character(
                                     cut(
                                             df_small$pBCE049_1,
                                             breaks = c(0, 50, 60, 70, 80, 90, 100),
                                             labels = c("0-50", "51-60", "61-70", "71-80", "81-90", "90-100"),
                                             include.lowest = TRUE
                                     )
                             ))

df_small$pSF_088_nores <- paste0("pSF_088_nores_",
                          as.character(
                                  cut(
                                          df_small$pSF_088_nores,
                                          breaks = c(0, 50, 60, 70, 80, 90, 100),
                                          labels = c("0-50", "51-60", "61-70", "71-80", "81-90", "90-100"),
                                          include.lowest = TRUE
                                  )
                          ))

df_small$pAPEC_O2_ColV <- paste0("pAPEC_O2_ColV_",
                          as.character(
                                  cut(
                                          df_small$pAPEC_O2_ColV,
                                          breaks = c(0, 50, 60, 70, 80, 90, 100),
                                          labels = c("0-50", "51-60", "61-70", "71-80", "81-90", "90-100"),
                                          include.lowest = TRUE
                                  )
                          ))

df_small$pUTI89 <- paste0("pUTI89_",
                          as.character(
                                  cut(
                                          df_small$pUTI89,
                                          breaks = c(0, 50, 60, 70, 80, 90, 100),
                                          labels = c("0-50", "51-60", "61-70", "71-80", "81-90", "90-100"),
                                          include.lowest = TRUE
                                  )
                          ))

df_small$pU1_51_B10 <- paste0("pU1_51_B10_",
                          as.character(
                                  cut(
                                          df_small$pU1_51_B10,
                                          breaks = c(0, 50, 60, 70, 80, 90, 100),
                                          labels = c("0-50", "51-60", "61-70", "71-80", "81-90", "90-100"),
                                          include.lowest = TRUE
                                  )
                          ))




rownames(df_small) <- rownames(df)

#### Determination of reference plasmids as present or absent ####
plasmid_data <- df_small2

plasmid_data$working_name <- rownames(df_small)

ColV_only <- Metadata %>% select(working_name, `ColV (328/668)`)

colV_comparisons_ <- plasmid_data %>% select(working_name, pSF_088_nores, pAPEC_O2_ColV, pUTI89, pU1_51_B10, pBCE049_1, Revised_Source_Niche)

colV_comparisons <- colV_comparisons_

colV_comparisons$sums<- rowSums(colV_comparisons[2:6])

colV_comparisons_ <- left_join(colV_comparisons_, ColV_only)

colV_comparisons_$pSF_088_nores[colV_comparisons_$pSF_088_nores < 80] <- 0
colV_comparisons_$pSF_088_nores[colV_comparisons_$pSF_088_nores >= 80] <- 1

colV_comparisons_$pAPEC_O2_ColV[colV_comparisons_$pAPEC_O2_ColV < 80] <- 0
colV_comparisons_$pAPEC_O2_ColV[colV_comparisons_$pAPEC_O2_ColV >= 80] <- 1

colV_comparisons_$pUTI89[colV_comparisons_$pUTI89 < 80] <- 0
colV_comparisons_$pUTI89[colV_comparisons_$pUTI89 >= 80] <- 1

colV_comparisons_$pU1_51_B10[colV_comparisons_$pU1_51_B10 < 80] <- 0
colV_comparisons_$pU1_51_B10[colV_comparisons_$pU1_51_B10 >= 80] <- 1

colV_comparisons_$pBCE049_1[colV_comparisons_$pBCE049_1 < 80] <- 0
colV_comparisons_$pBCE049_1[colV_comparisons_$pBCE049_1 >= 80] <- 1

colV_pos <- colV_comparisons_ %>% filter(`ColV (328/668)` == 4) %>% group_by(pSF_088_nores, pAPEC_O2_ColV) %>% summarise(counts = n()) 

colv_neg <- colV_comparisons_ %>% filter(`ColV (328/668)` == 0) %>% group_by(pUTI89, pU1_51_B10, pBCE049_1) %>% summarise(counts = n()) 

all_plas <- colV_comparisons_ %>% group_by(`ColV (328/668)`, pSF_088_nores, pAPEC_O2_ColV, pUTI89, pU1_51_B10, pBCE049_1) %>% summarise(counts = n()) 

plas_neg <- colV_comparisons %>% filter(sums == 0) %>% select(working_name) %>% as.data.frame()

plas_neg_pMLSTs <- Metadata %>% filter(working_name %in% plas_neg$working_name) %>% group_by(pMLST) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pAPEC_O2_ColV == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pSF_088_nores == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pU1_51_B10 == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pUTI89 == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

colV_comparisons_ %>% filter(pBCE049_1 == 1) %>% group_by(Revised_Source_Niche) %>% summarise(counts = n())

#### Generation of a separate df_small for Figure 1 which collapses plasmids into a single column ####
rownames(df_small2) <- rownames(df_small)

df_small3 <- df_small %>% select(starts_with("p"), -plas_counts)

df_small3$working_name <- rownames(df_small2)

df_small3 <- df_small3 %>% melt(id.vars = "working_name") %>% filter(!grepl("0-50$", value))

df_small3 <- df_small3 %>% filter(!grepl("variable", variable)) %>% select(working_name, value)

df_small3 <- df_small3 %>% filter(working_name != "ESC_KB9729AA_AS" & value  != "pSF_088_nores_81")

df_small3 <- df_small3 %>% rename(plasmid_map = value)

df_small4 <- df_small

df_small4$working_name <- rownames(df_small)

nems <- df_small4$working_name

df_small4 <- left_join(df_small4, df_small3)

extra_data <- Metadata %>% select(working_name, Pathogen)

df_small4 <- left_join(df_small4, extra_data)

df_small4 <- df_small4 %>% select(Pathogen, Revised_Source_Niche, plasmid_map, intI1, class_res_counts, ESBL_ress, plas_counts)

df_small4$Pathogen <- gsub("^","_pathogen_status_", df_small4$Pathogen)

rownames(df_small4) <- nems

pMLST_vs_plasmid_hit <- left_join(df_small3, pMLST_)

Metadata$`ColV (328/668)`<- gsub("4","ColV_pos",Metadata$`ColV (328/668)`)
Metadata$`ColV (328/668)` <- gsub("0","ColV_neg",Metadata$`ColV (328/668)`)

#### Generation of Scoary trait table ####

scoary_traits <- Metadata %>% select(working_name, Revised_Source_Niche, Pathogen, HC50, HC200)

scoary_traits$Revised_Source_Niche <- replace_na(scoary_traits$Revised_Source_Niche, replace = "Unknown")

scoary_traits$HC50 <- gsub("^","HC50_",scoary_traits$HC50)
scoary_traits$HC200 <- gsub("^","HC200_",scoary_traits$HC200)

scoary_traits$present <- rep(1, nrow(scoary_traits))

scoary_traits_2 <- dcast(scoary_traits, working_name ~ Revised_Source_Niche)

scoary_traits_2 <- dcast(scoary_traits, working_name ~ Pathogen) %>% left_join(scoary_traits_2, by = "working_name")

scoary_traits_2 <- dcast(scoary_traits, working_name ~ HC50) %>% left_join(scoary_traits_2, by = "working_name")

scoary_traits_2 <- dcast(scoary_traits, working_name ~ HC200) %>% left_join(scoary_traits_2, by = "working_name")

scoary_traits_2[is.na(scoary_traits_2)] <- 0

colnames(scoary_traits_2)[1] <- ""

#write.csv(scoary_traits_2, "delims/scoary_traits.csv", row.names = FALSE)
