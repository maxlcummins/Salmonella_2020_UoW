library(ggtree)
library(phytools)

#Load paths to files needed for the script
#tree_path <- "analysis/snippy/AVC171_all/fasttree/AVC171_chromosome.clean.fullcore.tree"
tree_path <- "~/Desktop/core_gene_alignment_snp_sites.aln.tree"
abricate_path <- "analysis/abricate/AVC171_all/abricate.txt"  
pointfinder_path <- "analysis/pointfinder/AVC171_all_pointfinder.txt"
ColV_data_path <- "analysis/abricate/AVC171_all/ColV_markers.txt"
pMLST_data <- "analysis/pMLST/pMLST_results.txt"
output_name <- "ST95_all"
refname <- "AVC171"
cgMLST_path <- "metadata/curated_metadata_all.txt"

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

#Read in the tree file
tree <-
        read.tree(file = tree_path)

#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (tree tip label)
tree$tip.label <- gsub("Reference", refname, tree$tip.label)

if(tree_path ==  "~/Desktop/core_gene_alignment_snp_sites.aln.tree"){
        tree <- phangorn::midpoint(tree = tree)
}

#Read in cgMLST data
Metadata <- read_delim(cgMLST_path,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

#Add new column working_name
Metadata$working_name <- Metadata$`Assembly barcode`

#replace the assembly barcode with strain name
Metadata$working_name <-
        gsub("ESC_SA8243AA_AS", "AVC171", Metadata$working_name)

#Filter metadata table to only contain strains from the tree
Metadata <- Metadata %>% filter(working_name %in% tree$tip.label)

#Change the column order to put working_name first
Metadata <- Metadata %>% select(working_name, everything())

#set rowname to be equal to working name
rownames(Metadata) <- Metadata$working_name

#Remove spaces from column names for metadata - this usually causes issues
colnames(Metadata) <- gsub(" ", "_", colnames(Metadata))

#Provides clade numbers to colour by.
#If you dont know what the node labels are use the line below under "get node labels"
#You may need to change the sizing..
#tree2 <- groupClade(tree, c(87, 86))

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

##Designate strains as HC5-4181 or Other
Metadata$HC50_or_other <-
        gsub("^1106$*", "HC50-1106", Metadata$HC50)
Metadata$HC50_or_other <-
        gsub("^[0-9].*", "Other", Metadata$HC50_or_other)

##Designate strains as HC200_Other
Metadata$HC200_Other <-
        gsub("^(8375|23651|28050|98164|1111|70040|80967|104870|121272|133658)$", "Other", Metadata$HC200)

HC100_others <- Metadata$HC100 %>% table() %>% as.data.frame() %>% filter(Freq < 3)

#Designate strains as HC100_Other
Metadata$HC100_Other <- Metadata$HC100

Metadata$HC100_Other[Metadata$HC100_Other %in% HC100_others$.] <- "Other"

#Designate strains as HC50_other
HC50_others <- Metadata$HC50 %>% table() %>% as.data.frame() %>% filter(Freq < 10)

#Designate strains as HC50_others
Metadata$HC50_Other <- Metadata$HC50

Metadata$HC50_Other[Metadata$HC50_Other %in% HC50_others$.] <- "Other"

ST95_all_simple_summary_N90L90 <-
        ST95_all_simple_summary_N90L90 %>% filter(name %in% tree$tip.label)

df <- ST95_all_simple_summary_N90L90

df <- df %>% select(ColV, everything(),-name)

df$ColV <- as.character(df$ColV)

df[df > 1] <- 1

colsum <- cbind(colnames(df), colSums(df)) %>% as.data.frame()

colsum$V2 <- as.numeric(colsum$V2)

df <- data.frame(lapply(df, as.numeric), stringsAsFactors = FALSE)

df <- df[, colSums(df != 0) > 0]

#Remove unwanted AMR genes from card
df <- df %>% select(
        -starts_with("card_acr"),-starts_with("card_bac"),-starts_with("card_bae"),-starts_with("card_cpxA"),-starts_with("card_CRP"),-starts_with("card_emr"),-starts_with("card_eptA"),-starts_with("card_Escherichia_coli_acrA"),-starts_with("card_Escherichia_coli_amp"),-starts_with("card_Escherichia_coli_emrE"),-starts_with("card_Escherichia_coli_mdfA"),-starts_with("card_evg"),-starts_with("card_gad"),-starts_with("card_H.NS"),-starts_with("card_kdpE"),-starts_with("card_marA"),-starts_with("card_mdt"),-starts_with("card_msbA"),-starts_with("card_pmrF"),-starts_with("card_tolC"),-starts_with("card_ugd"),-starts_with("card_yoj")
)

df <- df %>% select(
        -starts_with("EC_custom_malX"),-matches("EC_custom_eit[B-D]"),-starts_with("EC_custom_VGI"),-starts_with("EC_custom_malX"),-starts_with("EC_custom_yeeT"),-starts_with("EC_custom_pap"),-starts_with("EC_custom_iutA"),-starts_with("EC_custom_iucD"),-starts_with("EC_custom_irp"),-starts_with("EC_custom_pap"),-starts_with("EC_custom_fim"),-starts_with("EC_custom_fyuA")
)

df <- df %>% select(-starts_with("ISfinder_Feb_2020"))

df <- df %>% select(
        -starts_with("vfdb_chu"),-starts_with("vfdb_ent"),-starts_with("vfdb_chu"),-starts_with("vfdb_fepB|C|D"),-starts_with("vfdb_chu"),-starts_with("vfdb_fimA|B|C|D|E|F|G|I"),-starts_with("vfdb_gsp[D-M]"),-starts_with("vfdb_iro[B-D]"),-starts_with("vfdb_iuc[ABC]"),-starts_with("vfdb_iro[B-D]"),-starts_with("vfdb_pap[DEFHIJKX]"),-starts_with("vfdb_sfa[ABCDEFGHXY]"),-starts_with("vfdb_yag[WXYZ]"),-starts_with("vfdb_ybt[EPQSTUX]"),-starts_with("vfdb_ykgK")
)

df <- df %>% select(
        -matches("vfdb_chu"),-matches("vfdb_ent"),-matches("vfdb_chu"),-matches("vfdb_fep[BCD]"),-matches("vfdb_chu"),-matches("vfdb_fim[ABCDEFGI]"),-matches("vfdb_gsp[D-M]"),-matches("vfdb_iro[B-D]"),-matches("vfdb_iuc[ABC]"),-matches("vfdb_iro[B-D]"),-matches("vfdb_pap[BDEFHIJKX]"),-matches("vfdb_sfa[ABCDEFGHXY]"),-matches("vfdb_yag[WXYZ]"),-matches("vfdb_ybt[EPQSTUX]"),-matches("vfdb_ykgK")
)

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

r <- df %>% select(starts_with("res_"))
p <- df %>% select(starts_with("plas_"))
cv <- df %>% select(starts_with("ColV"))
v <- df %>% select(starts_with("vir_"))

p[p == 1] <- 2
v[v == 1] <- 3
cv[cv == 1] <- 4

df <- cbind(cv, r, p, v)

df$vir_sitA <- df %>% select(starts_with("vir_sitA")) %>% rowSums()

df$vir_sitA[df$vir_sitA > 1] <- 3

df <- df %>% select(-starts_with("vir_sitA_"), -starts_with("vir_kpsMT"))

df$plas_IncB.O.K.Z <-
        df %>% select(starts_with("plas_IncB")) %>% rowSums()

df$plas_IncB.O.K.Z[df$plas_IncB.O.K.Z > 1] <- 2

df <- df %>% select(-starts_with("plas_IncB.O.K.Z_"))

df <- df[, order(names(df))]

colnames_save <- gsub("(res|vir|plas)_", "", colnames(df))

df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

colnames(df) <- colnames_save

rownames(df) <- ST95_all_simple_summary_N90L90$name

df6 <- read_csv("delims/pCERC4_plasmid_coverage_percentage.csv")

df6 <- df6[,2:ncol(df6)]

Metadata <- left_join(Metadata, df6)

df6 <- read_csv("delims/pEC14_114_plasmid_coverage_percentage.csv")

df6 <- df6[,2:ncol(df6)]

plasref2 <- "pEC14_114"

colnames(df6) <- c(plasref2, "working_name")

Metadata <- left_join(Metadata, df6)

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

#write_csv(metadata, "processed_data/metagenodata.csv")

#Metadata$`ColV (329/669)`<- gsub("4","F35a",Metadata$`ColV (329/669)`)

p2 <- ggtree(tree, layout = 'circular',
             open.angle = 0) %<+%
        Metadata +
        geom_tiplab(size = 0.8,
                    align = TRUE,
                    linesize = 0.15,
                    aes(color = as.factor(`ColV (328/668)`))
                    ) +
        geom_tippoint(size = 0.2, aes(colour=as.factor(HC200_Other)))

#list_fields_circ_tree <- unique(c(Metadata$pMLST, Metadata$HC50_or_other, Metadata$`ColV (329/669)`, Metadata$Revised_Source_Niche)) %>% as.data.frame() %>% write.csv("~/names.csv")



colval <- c("white", #blank
            'orange', # no_colV (zero)
            "#80b1d3", #2
            "#fb8072", #3
            "black", #4
            #"brown", # Bovine
            "gold", # Canine
            # "green", # Environment
            "blue", # Human
            # "black", #Other
            "red" # Poultry
            # "grey" #?
)

colval <- c("white", #0
            "#bebada", #1
            "#80b1d3", #2
            "#fb8072", #3
            "black", #4
            #"brown", # Bovine
            "gold", # Canine
            # "green", # Environment
            "blue", # Human
            # "black", #Other
            "red" # Poultry
            # "grey" #?
)

df_small <- as.data.frame(df$`ColV (328/668)`)
colnames(df_small) <- gsub("^df.*", "ColV", colnames(df_small))
df_small$working_name <- rownames(df)
df_small <- left_join(df_small, Metadata)
df_small <- as.data.frame(df_small)
df_small <- df_small %>% select(Revised_Source_Niche)

rownames(df_small) <- rownames(df)

phylogenetic_tree <- gheatmap(
        p = p2,
        data = df_small,
        colnames_offset_y = -0.1,
        font.size = 1.5,
        hjust = 0,
        #colnames_position = "top",
        colnames = FALSE,
        #colnames_angle = 90,
        offset = 0.006,
        width = 0.05,
        color = 'grey'
) + 
        #ggplot2::ylim(NA, 700) #+
        theme(legend.position = "none") +
        #        legend.position = "bottom",
        #        legend.title = element_blank(),
        #        legend.text = element_text(size = 14),
        #        legend.box = "horizontal"
        #) #+
        #scale_fill_manual(
        #        aesthetics = c("colour", "fill"),
        #        values = colval,
        #        na.value = 'grey'
        #) +
#geom_tiplab2(
#        aes(image = Flag),
#        geom = "image",
#        size = 0.004,
#        align = TRUE,
#        linetype = NULL,
#        offset = 0.012
#) +
#geom_tiplab2(
#        aes(image = Revised_Source_Niche_img),
#        geom = "image",
#        size = 0.004,
#        align = TRUE,
#        linetype = NULL,
#        offset = 0.013
#) +
#geom_tiplab2(
#        aes(image = Pathogen_img),
#        geom = "image",
#        size = 0.004,
#        align = TRUE,
#        linetype = NULL,
#        offset = 0.014
#) +
geom_tiplab2(
        aes(label = pEC14_114),
        align = TRUE,
        linetype = NULL,
        offset = 0.017,
        size = 0.8
) +
geom_tiplab2(
        aes(label = ColV_percent_hit),
        align = TRUE,
        linetype = NULL,
        offset = 0.015,
        size = 0.8
        ) +
geom_tiplab2(
        aes(label = pMLST),#, colour = pMLST),
        align = TRUE,
        linetype = NULL,
        offset = 0.011,
        size = 0.8
) 

phylogenetic_tree

p <- ggtree(tree) %<+%
        Metadata +
        geom_tiplab(size = 0.5,
                    align = TRUE,
                    offset = 0.007,
                    linesize = 0.05#,
                    #aes(color = HC50_or_other)
        )# +
#geom_tippoint(size = 3, aes = (color = HC50_or_other))

genotypic_table <- gheatmap(
        p = p,
        data = df,
        colnames_offset_y = -0.1,
        font.size = 1.5,
        hjust = 0,
        colnames_position = "top",
        #colnames = FALSE,
        colnames_angle = 90,
        offset = 0.015,
        width = 2.5,
        color = 'grey'
        ) + 
        ggplot2::ylim(NA, 700) +
        theme(legend.position = "none") +
        #        legend.position = "bottom",
        #        legend.title = element_blank(),
        #        legend.text = element_text(size = 14),
        #        legend.box = "horizontal"
        #) +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = colval,
                na.value = 'grey'
        )
        #geom_tiplab(
        #        aes(image = Flag),
        #        geom = "image",
        #        size = 0.016,
        #        align = TRUE,
        #        linetype = NULL,
        #        offset = 6.5
        #) +
        #geom_tiplab(
        #        aes(image = Revised_Source_Niche_img),
        #        geom = "image",
        #        size = 0.013,
        #        align = TRUE,
        #        linetype = NULL,
        #        offset = 8.5
        #) +
        #geom_tiplab(
        #        aes(image = Pathogen_img),
        #        geom = "image",
        #        size = 0.013,
        #        align = TRUE,
        #        linetype = NULL,
        #        offset = 10
        #) +
        #geom_tiplab(
        #        aes(label = Revised_Collection_Year),
        #        align = TRUE,
        #        linetype = NULL,
        #        offset = 12,
        #        size = 1.7
        #) 

#pdf()

listy2 <- read_csv("delims/pCERC4_plasmid_coverage.csv")

listy_2 <- listy2[,2:ncol(listy2)]

rownames(listy_2) <- listy2$X1

plasmid_tree <- gheatmap(
        p = p,
        data = listy_2,
        #colnames_offset_y = -0.1,
        font.size = 1,
        hjust = 0,
        #colnames_position = "top",
        colnames = FALSE,
        colnames_angle = 90,
        offset = 0.02,
        width = 2.5,
        color = NULL,
        low = 'white',
        high = '#fb8072'
)+
        geom_tiplab(
                aes(image = Flag),
                geom = "image",
                size = 0.002,
                align = TRUE,
                linetype = NULL,
                offset = 0.012
        ) +
        geom_tiplab(
                aes(image = Revised_Source_Niche_img),
                geom = "image",
                size = 0.001,
                align = TRUE,
                linetype = NULL,
                offset = 0.0125
        ) +
        geom_tiplab(
                aes(image = Pathogen_img),
                geom = "image",
                size = 0.001,
                align = TRUE,
                linetype = NULL,
                offset = 0.013
        ) +
#        geom_tiplab(
#                aes(label = Revised_Collection_Year),
#                align = TRUE,
#                linetype = NULL,
#                offset = 19,
#                size = 0.5
#        ) +
        geom_tiplab(
                aes(label = ColV_percent_hit),
                align = TRUE,
                linetype = NULL,
                offset = 0.018,
                size = 0.5
        ) +
        geom_tiplab(
                aes(label = pMLST, color = pMLST),
                align = TRUE,
                linetype = NULL,
                offset = 0.014,
                size = 0.5
        ) +
        theme(legend.position = "none") #+
        #ggplot2::ylim(NA, 700)# +
#geom_image(x = 104.5, y = 48, image = "../Linear_Map.png", size = 0.630) 

Liu <- as.data.frame(df$`ColV (329/669)`)
Liu$working_name <- rownames(df)
ColV <- left_join(df6, Liu)

colnames(ColV)[3] <- "Liu"

ColV %>% filter(Liu == 0) %>% group_by(ColV_percent_hit) %>% summarise(counts = n()) %>% View()


ColV$ColV_70 <- ColV$ColV_percent_hit
ColV$ColV_70[ColV$ColV_70 < 70] <- 0
ColV$ColV_70[ColV$ColV_70 >= 70] <- 1
ColV$Liu[ColV$Liu == 4] <- 1

ColV <- ColV %>% select(working_name, ColV_percent_hit, ColV_70, Liu)

write_csv(ColV, path = "delims/ColV_groups.csv")
