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

path_to_tree <- "analysis/snippy/HC50_1106/fasttree/AVC171.clean.fullcore.tree"

path_to_abricate <- "analysis/abricate/HC50_1106/abricate.txt"

path_to_metadata <- "metadata/AVC171-HC50_curated.txt"

path_to_cgMLST <- "metadata/cgMLST.txt"

path_to_pointfinder <- "analysis/pointfinder/pointfinder_results.txt"

refname <- "AVC171"

#Changes working directory
setwd(path_to_repo)


#Read in the tree file
tree <-
        read.tree(file = path_to_tree)



#trim the names of the assemblies in the tree tip labels
tree$tip.label <- gsub("\\..*", "", tree$tip.label)

#reassign strain SG17-135's assembly barcode to its strain name (tree tip label)
tree$tip.label <- gsub("Reference", refname, tree$tip.label)

#Read in the abricate genotype data sheet (small number of rows for colname reassignment)
df <- read_delim(
        path_to_abricate,
        "\t",
        escape_double = FALSE,
        trim_ws = TRUE,
        n_max = 10
)

#Colname reassignment
colnames(df)[c(1, 10:11)] <-
        c("name", "perc_coverage", "perc_identity")
df_colnames <- colnames(df)

#Read in full abricate genotype data sheet
df <- read_delim(
        path_to_abricate,
        "\t",
        escape_double = FALSE,
        trim_ws = TRUE,
        col_names = FALSE,
        skip = 1
)

#Remove cases where there are multiple headers from concatenation of abricate reports
df <- df %>% filter(X2 != "SEQUENCE")

#Colname reassignment
colnames(df) <- df_colnames

#Convert percent coverage and identity to numeric type to allow filtering
df$perc_coverage <- as.numeric(df$perc_coverage)
df$perc_identity <- as.numeric(df$perc_identity)

#Filter to perc_coverage and perc_identity > 90%
df <-
        df %>% filter(perc_coverage > 90) %>% filter(perc_identity > 90)

#Trim excess characters the assembly names and reassign this to rownames
df$name <- gsub("\\..*", "", df$name)

df <- df %>% filter(name %in% tree$tip.label)

#Filter out non-virulence-associated genes
#df <-
#  df %>% filter(DATABASE %in% c("vfdb", "card", "plasmidfinder"))


################################################################################
####Remove columns that we dont need
################################################################################

#Change pap genes to be operons where they co-occur
df <- df[!grepl("pap(B|C|D|E|F|H|I|J)", df$GENE), ]
df$GENE <- gsub("papA", "papABCDEFHIJ", df$GENE)

#Change ybt genes to be operons where they co-occur
df <- df[!grepl("ybt(E|P|Q|S|U|X)", df$GENE), ]
df$GENE <- gsub("ybtA", "ybtAEPQSUX", df$GENE)

#Change yag genes to be operons where they co-occur
df <- df[!grepl("yag(W|Y|Z)", df$GENE), ]
df$GENE <- gsub("yagV", "yagVWYZ/ecpEDBA", df$GENE)

#Change iuc genes to be operons where they co-occur
df <- df[!grepl("iuc(B|C|D)", df$PRODUCT), ]
df$GENE <- gsub("iucA", "iucABCD", df$GENE)

#Change iro genes to be operons where they co-occur
df <- df[!grepl("iro(C|D|E|N)", df$GENE), ]
df$GENE <- gsub("iroB", "iroBCDEN", df$GENE)

#Change gsp genes to be operons where they co-occur
df <- df[!grepl("gsp(D|E|G|H|I|J|L|M)", df$GENE), ]
df$GENE <- gsub("gspC", "gspCDEGHIJLM", df$GENE)

#Change gsp genes to be operons where they co-occur
df <- df[!grepl("fim(B|C|E|F|G|H|I)", df$GENE), ]
df$GENE <- gsub("fimA", "fimABCEFGHI", df$GENE)

#Change chu genes to be operons where they co-occur
df <- df[!grepl("chu(S|T|U|V|W|X|Y)", df$GENE), ]
df$GENE <- gsub("chuA", "chuASTUVWXY", df$GENE)

#Change ent genes to be operons where they co-occur
df <- df[!grepl("ent(B|C|D|E|F|S)", df$GENE), ]
df$GENE <- gsub("entA", "entABCDEFS", df$GENE)

#Change ent genes to be operons where they co-occur
df <- df[!grepl("fep(B|C|D|G)", df$GENE), ]
df$GENE <- gsub("fepA", "fepABCDG", df$GENE)

#Change emr genes to be operons where they co-occur
df <- df[!grepl("emr(B|K|R|Y)", df$GENE), ]
df$GENE <- gsub("emrA", "emrABKRY", df$GENE)

#Change mdt genes to be operons where they co-occur
df <- df[!grepl("mdt(B|C|F|G|H|N|O|P)", df$GENE), ]
df$GENE <- gsub("mdtA", "mdtABCFGHNOP", df$GENE)

#Change mdt genes to be operons where they co-occur
df <- df[!grepl("eit(B|C|D)", df$GENE), ]
df$GENE <- gsub("eitA", "eitABCD", df$GENE)


################################################################################
####Clean gene names
################################################################################

#AMR genes
df$GENE <- gsub("^AAC", "aac", df$GENE)
df$GENE <- gsub("^ANT", "ant", df$GENE)
df$GENE <- gsub("^APH","aph", df$GENE)
df$GENE <- gsub("CMY-59", "blaCMY-59", df$GENE)
df$GENE <- gsub("CTX-M-55", "blaCTX-M-55", df$GENE)
df$GENE <- gsub("TEM-1", "blaTEM-1", df$GENE)
df$GENE <- gsub("FosA7", "fosA7", df$GENE)
df$GENE <- gsub("QnrS1", "qnrS1", df$GENE)

#Plasmid genes
df$GENE <- gsub("_Gamma_1", "(Gamma)", df$GENE)
df$GENE <- gsub("_1_Alpha", "(Alpha)", df$GENE)
df$GENE <- gsub("_1.*","",df$GENE)


################################################################################
####Remove unwanted AMR genes (not considered actual AMR genes by most)
################################################################################

df <- df %>% filter(GENE != "golS")
df <- df %>% filter(GENE != "sdiA")


################################################################################
####Replace strA/B names (more commonly used names for these genes)
################################################################################

df$GENE <- gsub("aph\\(3''\\)-Ib", "strA", df$GENE)
df$GENE <- gsub("aph\\(6\\)-Id", "strB", df$GENE)


#Prepend the database name to the gene name so we can identify specific databases later
df$GENE <- paste0(df$DATABASE, "_", df$GENE)

#trim db name to something smaller
df$GENE <- gsub("^vfdb", "v", df$GENE)
df$GENE <- gsub("^card", "r", df$GENE)
df$GENE <- gsub("^plasmidfinder", "p", df$GENE)
df$GENE <- gsub(".*EC_custom_Feb_2020_","c_", df$GENE)


#Read in pointfinder data on resistance snps
res_snps <- read_delim(path_to_pointfinder, 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

#fix colnames of res_snps for consistency
colnames(res_snps) <- gsub("name_file","name", colnames(res_snps))

#fix colnames of res_snps for consistency
colnames(res_snps) <- gsub("Mutation","GENE", colnames(res_snps))

#fix colnames of res_snps for consistency
res_snps$GENE <- gsub("^","r_", res_snps$GENE)

#fix sample names of res_snps for consistency
res_snps$name <- gsub("/.*.tsv","",res_snps$name)

#Cast the dataframe into a more usable format
res_snps2 <- dcast(data = res_snps, name ~ GENE, length, drop = FALSE)

#Select the columns "name" and "GENE"
df2 <- df %>% select(contains("name"), contains("GENE"))

#Cast the dataframe into a more usable format
df3 <- dcast(data = df2, name ~ GENE, length, drop = FALSE)

#Bind Res_SNP data to abricate data
df3 <- left_join(df3, res_snps2)

#Replace NAs (indicating a lack of a gene hit) with a zero
df3[is.na(df3)] <- 0

#Trim excess characters the assembly names and reassign this to rownames
rownames(df3) <- gsub("\\..*", "", df3$name)

#Remove the column that has assembly names
df3 <- df3[, 2:ncol(df3)]

#Change multiple hits to a binary hit (gene considered present or absent, regardless of copy number)
df3[df3 > 1] <- 1

#Sum the number of hits per gene
colSums(df3) -> sums

#reassign Colnames to include the number the format "gene_A 50% (n=41/82)"
colnames(df3) <-
        paste0(colnames(df3),
               " ",
               sums,
               "/",
               nrow(df3),
               " (",
               round(sums / nrow(df3) * 100),
               "%)")

#reassign strain SG17-135's assembly barcode to its strain name (data frame with genotypic data)
#rownames(df3) <- gsub("SAL_HC4750AA_AS", "SG17-135", rownames(df3))


#isolate the individual DBs so we can edit their contents for colorisation in the heatmap later
res <- df3 %>% select(starts_with("r"))
plas <- df3 %>% select(starts_with("p"))
vir <- df3 %>% select(starts_with("v_"))
cus <- df3 %>% select(starts_with("c_"))

#remove any remaining columns with no hits
res <- res[,colSums(res) > 0]
plas <- plas[,colSums(plas) > 0]
vir <- vir[,colSums(vir) > 0]
cus <- cus[,colSums(cus) > 0]


#DB specific cell content changing for colorisation in the heatmap later
res[res == 1] <- "R"
plas[plas == 1] <- "P"
vir[vir == 1] <- "V"
cus[cus == 1] <- "C"

#sort columnanmes alphabetically
res <- res %>% select(sort(names(.)))
plas <- plas %>% select(sort(names(.)))
vir <- vir %>% select(sort(names(.)))
cus <- cus %>% select(sort(names(.)))

#bind all the DBs again together
df4 <- cbind(res, plas, vir, cus)

#asign rownames again
#rownames(df4) <- df3$name

#replace 0's with Ns
df4[df4 == 0] <- "N"

#Read in metadata
Metadata <- read_delim(path_to_metadata,
                       "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

#Read in cgMLST data
cgMLST <- read_delim(path_to_cgMLST,
                     "\t",
                     escape_double = FALSE,
                     trim_ws = TRUE)

#Select HC columns from cgMLST
cgMLST <- cgMLST %>% select(Uberstrain, starts_with("HC"))

#bind Metadata and cgMLST tables
Metadata <- left_join(Metadata, cgMLST)

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

#Create a list of colnames for HCC columns from cgMLST table
cols <- Metadata %>% select(starts_with("HC")) %>% colnames()

#reassign the above columns as factors
Metadata[cols] <- lapply(Metadata[cols], factor)

##Designate strains as HC5-4181 or Other
#Metadata$HC5_or_other <- gsub("^4181$*", "HC5-4181", Metadata$HC5)
#Metadata$HC5_or_other <-
#        gsub("^[0-9].*", "Other", Metadata$HC5_or_other)

#Replace ND's in Source Niche with NA
Metadata$`Source Niche` <-
        gsub("ND", NA, Metadata$`Source Niche`)

#Original colour scheme for heatmap 
#Define colors for gene-type dependent coloring of gene hits
#colorgenotype <-
#  c(
#    "N" = "white",
#    "R" = "#bebada",
#    "V2" = "black",
#    "V" = "#fb8072",
#    "P" = "#80b1d3"
#  )

df4 <- df4 %>%
        #select(starts_with("v_"),starts_with("v2_"))
        select(starts_with("r_"),
               starts_with("p_"),
               starts_with("v_"),
               starts_with("c_"))
               
#rownames(df4) <- df3$name


#Remove the DB prepend from df4 column names
colnames(df4) <- gsub("^[^_]+_", "", colnames(df4))

#Remove spaces from column names for metadata - this usually causes issues
colnames(Metadata) <- gsub(" ", "_", colnames(Metadata))

#Provides clade numbers to colour by.
#If you dont know what the node labels are use the line below under "get node labels"
#You may need to change the sizing..
#tree2 <- groupClade(tree, c(87, 86))

Metadata$Flag <- Metadata$Revised_Country
Metadata$Flag <- gsub("Australia", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Australia.jpg", Metadata$Flag)
Metadata$Flag <- gsub("Denmark", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Denmark.jpg", Metadata$Flag)
Metadata$Flag <- gsub("United States", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/USA.png", Metadata$Flag)
Metadata$Flag <- gsub("^Ireland$", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Ireland.png", Metadata$Flag)
Metadata$Flag <- gsub("Vietnam", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Vietnam.png", Metadata$Flag)
Metadata$Flag <- gsub("United Kingdom", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/United_Kindgom.png", Metadata$Flag)
Metadata$Flag <- gsub("Mexico", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Mexico.png", Metadata$Flag)
Metadata$Flag <- gsub("Ghana", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Ghana.png", Metadata$Flag)
Metadata$Flag <- gsub("United Arab Emirates", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/UAE.png", Metadata$Flag)
Metadata$Flag <- gsub("Scotland", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Scotland.png", Metadata$Flag)
Metadata$Flag <- gsub("Germany", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Germany.png", Metadata$Flag)
Metadata$Flag <- gsub("Northern Ireland", "https://raw.githubusercontent.com/maxlcummins/SG17-135/resubmission/flags/Northern_Ireland.jpg", Metadata$Flag)
#Metadata$Flag[is.na(Metadata$Flag)] <- "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question.png"


Metadata$Classification_img <- Metadata$Classification
Metadata$Classification_img <- gsub("SEPEC", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/blood.png", Metadata$Classification_img)
Metadata$Classification_img[is.na(Metadata$Classification_img)] <- "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question.png"
Metadata$Classification_img <- gsub("Faecal", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/colon.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("RMAE", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/meat.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("APEC", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/poultry.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("UPEC", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/urine.png", Metadata$Classification_img)
Metadata$Classification_img <- gsub("NMEC", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/brain.png", Metadata$Classification_img)

Metadata$Pathogen <- Metadata$Classification
Metadata$Pathogen <- gsub("NMEC", "Systemic", Metadata$Pathogen)
Metadata$Pathogen <- gsub("APEC", "Systemic", Metadata$Pathogen) 
Metadata$Pathogen <- gsub("UPEC", "Urine", Metadata$Pathogen) 
Metadata$Pathogen <- gsub("SEPEC", "Systemic", Metadata$Pathogen) 
Metadata$Pathogen <- gsub("RMAE", "Flora", Metadata$Pathogen) 
Metadata$Pathogen <- gsub("Faecal", "Flora", Metadata$Pathogen)

Metadata$Pathogen_img <- Metadata$Pathogen
Metadata$Pathogen_img <- gsub("Flora", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/colon.png", Metadata$Pathogen_img)
Metadata$Pathogen_img <- gsub("Systemic", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/blood.png", Metadata$Pathogen_img)
Metadata$Pathogen_img <- gsub("Urine", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/urine.png", Metadata$Pathogen_img)
Metadata$Pathogen_img[is.na(Metadata$Pathogen_img)] <- "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question.png"

Metadata$Revised_Source_Niche_img <- Metadata$Revised_Source_Niche
Metadata$Revised_Source_Niche_img <- gsub("Poultry_or_human", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question_2.png", Metadata$Revised_Source_Niche_img)
Metadata$Revised_Source_Niche_img <- gsub("Canine", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/dog.png", Metadata$Revised_Source_Niche_img)
Metadata$Revised_Source_Niche_img <- gsub("Poultry", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/poultry.png", Metadata$Revised_Source_Niche_img)
Metadata$Revised_Source_Niche_img <- gsub("Human", "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/human.png", Metadata$Revised_Source_Niche_img)
Metadata$Revised_Source_Niche_img[is.na(Metadata$Revised_Source_Niche_img)] <- "https://raw.githubusercontent.com/maxlcummins/AVC171/max/images/Icons/question.png"


#Generate the tree
p <- ggtree(tree) %<+%
        Metadata +
        geom_tiplab(size = 2,
                    align = TRUE,
                    offset = 0.007#,
                    #aes(color = HC5_or_other)
        ) +
        geom_tiplab(aes(image = Flag), geom="image", size = 0.016, align = TRUE, linetype = NULL, offset = 0.0825
        ) +
        geom_tiplab(aes(image = Revised_Source_Niche_img), geom="image", size = 0.01, align = TRUE, linetype = NULL, offset = 0.1000
        ) +
        geom_tiplab(aes(image = Pathogen_img), geom="image", size = 0.01, align = TRUE, linetype = NULL, offset = 0.1125
        ) 

p2 <- p + geom_tiplab(aes(label=Revised_Collection_Year), align = TRUE, linetype = NULL, offset = 0.1265, size = 2)
         + geom_treescale(x = 0.1, y = -2.25, offset = -2.4, fontsize = 2)#+
#get node labels
#geom_text2(size = 2, aes(subset=!isTip, label=node), hjust=-.5)


colval <- c(#"black",# APEC
            "yellow",# C - Custom DB
            #"brown", # Faecal
            #"#609bce", # blue # Human
            "white", # N
            "#80b1d3", # P
            #"red", # Poultry
            "#bebada",#purple # R
            #"purple", # Retail Poultry Meat
            #"red", # SEPEC
            #"yellow", #UPEC
            "#fb8072" # red # V
            )
        
        
#        #'steelblue', #  'clade 1'                       √
#        #'firebrick', # 'middle clade?' 
#        #'steelblue', #  'clade 2'                       √
#        '#b6933f', #      'Source Niche - Food'           √
#        'black', #      'HC5_or_other (HC5-4181)        √
#        '#609bce', #      'Source Niche - Human'          √
#        'white', #      'Genes - 0 (N)                  √
#        '#80b1d3', #      'Plasmid Genes - 1 (P)          √       
#        '#bebada', #      'AMR genes - 1 - (R)            √ 
#        '#fb8072', #        'VFDB genes -1 - (V)'         √
#        '#cd2a18', #  'PAI genes - 1 (V2)                 √    
#        '#8d68ca' #  'Source Niche - Wild Animal'        √
#)

#pdf(file = "figures/alpha_tree.pdf", paper = "a4r")

#Generate the heatmap
a <- gheatmap(
        p = p2,
        data = df4,
        colnames_offset_y = -0.1,
        font.size = 1,
        hjust = 0,
        colnames_position = "top",
        #colnames = FALSE,
        colnames_angle = 90,
        offset = 0.15,
        width = 1.7,
        color = 'grey'
) +
        theme(legend.position = "none") +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = colval,
                na.value = 'grey'
        )# +
        #ggplot2::ylim(NA, 100)

a
#dev.off()

#Generate the heatmap
zz <- gheatmap(
        p = p2,
        data = listy_,
        #colnames_offset_y = -0.1,
        font.size = 1,
        hjust = 0,
        #colnames_position = "top",
        colnames = FALSE,
        colnames_angle = 90,
        offset = 0.15,
        width = 1.7,
        #color = 'grey',
        low = 'white',
        high = 'red'
) #+
  #      theme(legend.position = "none") +
  #      scale_fill_manual(
  #              aesthetics = c("colour", "fill"),
  #              values = colval,
  #              na.value = 'grey'
  #      )# +
#ggplot2::ylim(NA, 100)

