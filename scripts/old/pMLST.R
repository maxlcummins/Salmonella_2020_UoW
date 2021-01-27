library(ggtree)

#Load paths to files needed for the script
tree_path <-
        "analysis/snippy/AVC171_chromosome_rm_nometa/fasttree/AVC171_chromosome.clean.fullcore.tree"
abricate_path <- "analysis/abricate/ST95_all/abricate.txt"
pointfinder_path <-
        "analysis/pointfinder/all/ST95_all_pointfinder.txt"
ColV_data_path <- "analysis/abricate/ST95_all/colV_abricate.txt"
output_name <- "HC50_1106"
refname <- "AVC171"
cgMLST_path <- "metadata/curated_metadata_all.txt"

#Load in abricateR script to amalgamate data
source("scripts/abricateR.R")

#Run abricateR
abricateR(
        file = abricate_path,
        output = output_name,
        identity = 90,
        length = 90,
        writecsv = TRUE,
        pointfinder_data = pointfinder_path,
        ColV_Liu_data = ColV_data_path
)

pMLST_results <- read_delim("analysis/pMLST/pMLST_results.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

pMLST_results <- pMLST_results %>% filter(Allele != "No hit found")

pMLST_results$name_file <- gsub("\\/.*", "", pMLST_results$name_file)

simple_summary <- dcast(
        data = pMLST_results,
        name_file ~ Locus,
        value.var = 'Allele',
        drop = FALSE,
        fun.aggregate=function(x) paste(x, collapse = "/ ")
)
simple_summary <- apply(simple_summary, 2, function(y) gsub("/ [A-Z]+_", "/", y)) %>% as.data.frame()

simple_summary$FIA[simple_summary$FIA == ""] <- "A-NULL"
simple_summary$FIB[simple_summary$FIB == ""] <- "B-NULL"
simple_summary$FIC[simple_summary$FIC == ""] <- "C-NULL"
simple_summary$FIIK[simple_summary$FIIK == ""] <- "K-NULL"
simple_summary$FII[simple_summary$FII == ""] <- "F-NULL"

simple_summary <- simple_summary %>% select(name_file, FII, FIA, FIB, FIC, FIIK)

pMLST <- simple_summary[,1:2]

pMLST$FII <- paste(simple_summary$FII, simple_summary$FIA, simple_summary$FIB, simple_summary$FIC, simple_summary$FIIK, sep = ":")

pMLST$FII <- gsub("FII_","F",pMLST$FII)
pMLST$FII <- gsub("FIA_","A",pMLST$FII)
pMLST$FII <- gsub("FIB_","B",pMLST$FII)
pMLST$FII <- gsub("FIC_","C",pMLST$FII)
pMLST$FII <- gsub("FIIK_","K",pMLST$FII)

pMLST$FII <- gsub("NULL","",pMLST$FII)

colnames(pMLST) <- c("name", "pMLST")
