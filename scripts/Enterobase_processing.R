library(reshape2)
library(ggtree)
library(phytools)
library(tidyr)
library(readr)
library(magrittr)
library(dplyr)
library(stringr)

# Read in enterobase data and metadata
cgMLST_path <- "delims/all_enterobase/cgMLST_13_4_21.txt"
MLST_path <- "delims/all_enterobase/Salm_MLST_13_4_21.txt"
wgMLST_path <- "delims/all_enterobase/Salm_wgMLST_13_4_21.txt"
rMLST_path <- "delims/all_enterobase/Salm_rMLST_13_4_21.txt"
serotype_path <- "delims/all_enterobase/Salm_serotype_13_4_21.txt"
assembly_stats_path <- "delims/all_enterobase/Salm_assembly_stats_13_4_21.txt"

# Read in enterobase data streams
assembly_stats <- read_delim(assembly_stats_path, 
                             "\t", escape_double = FALSE, trim_ws = TRUE)
cgMLST <- read_delim(cgMLST_path, 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
MLST <- read_delim(MLST_path, 
                   "\t", escape_double = FALSE, trim_ws = TRUE)

wgMLST <- read_delim(wgMLST_path, 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
rMLST <- read_delim(rMLST_path, 
                    "\t", escape_double = FALSE, trim_ws = TRUE)
serotype <- read_delim(serotype_path, 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

#Change the ST column names in different ST schemes to prevent identical colnames
cgMLST <- cgMLST %>% rename(cgMLST = "ST")
wgMLST <- wgMLST %>% rename(wgMLST = "ST")


#As above for allelic differences
cgMLST <- cgMLST %>% rename(cgMLST_differences = "Differences")
wgMLST <- wgMLST %>% rename(wgMLST_differences = "Differences")
rMLST <- rMLST %>% rename(rST_differences = "Differences")
rMLST <- rMLST %>% rename("rST Predicted Serotype" = "Serotype (Predicted)")

#Remove shared columns between sheets
cgMLST <- cgMLST[,c(1,33:ncol(cgMLST))]
MLST <- MLST[,c(1,33:ncol(MLST))]
wgMLST <- wgMLST[,c(1,33:ncol(wgMLST))]
rMLST <- rMLST[,c(1,33:ncol(rMLST))]
serotype <- serotype[,c(1,33:ncol(serotype))]

#Combine sheets
enterobase_data <- assembly_stats %>%
        left_join(cgMLST, by = "Uberstrain") %>%
        left_join(MLST, by = "Uberstrain") %>%
        left_join(wgMLST, by = "Uberstrain") %>%
        left_join(rMLST, by = "Uberstrain") %>%
        left_join(serotype, by = "Uberstrain")

#Remove old tables from memory - they are big!
#rm(assembly_stats)
rm(cgMLST)
rm(MLST)
rm(rMLST)
rm(serotype)
rm(wgMLST)

#Remove underscores from colnames and replace with spaces
colnames(enterobase_data) <- gsub(" ", "_", colnames(enterobase_data))

#Remove strains with negative and NaN STs
enterobase_data <- enterobase_data %>% filter(ST > 0)

#Remove strains with negative and NaN cgMLSTs
enterobase_data <- enterobase_data %>% filter(cgMLST > 0)

#Remove strains with dodgey contig number
#Remove samples with 0 contigs > 200
enterobase_data <- enterobase_data %>% filter(`Contig_Number(>=200_bp)` > 0)
#Remove samples with 
enterobase_data <- enterobase_data %>% filter(`Contig_Number(>=200_bp)` <= 800)

Type_data <- enterobase_data %>% select(Uberstrain, `Data_Source(Accession_No.;Sequencing_Platform;Sequencing_Library;Insert_Size;Experiment;Bases;Average_Length;Status)`)

Type_data <- Type_data %>% rename(Data_Source = 'Data_Source(Accession_No.;Sequencing_Platform;Sequencing_Library;Insert_Size;Experiment;Bases;Average_Length;Status)')

Type_data$Data_Source <- gsub(";;",";_;",Type_data$Data_Source)

Data_source <- str_split_fixed(Type_data$Data_Source, ",", 96) %>% as.data.frame()

Data_source$Uberstrain <- enterobase_data$Uberstrain

Data_source <- melt(data = Data_source, id.vars = "Uberstrain", factorsAsStrings = FALSE)

Data_source_2 <- str_split_fixed(Data_source$value, ";", 8) %>% as.data.frame()

Data_source_3 <- cbind(Data_source, Data_source_2) %>% as.data.frame()

Data_source_3 <- Data_source_3 %>% select(-variable, -value)

table(Data_source_3$V6) %>% View()

colnames(Data_source_3) <- c("Uberstrain","Accession_No.","Sequencing_Platform","Sequencing_Library","Insert_Size","Experiment","Bases","Average_Length","Status")

Data_source_3 <- Data_source_3 %>% filter(Bases != "")

Data_source_meta <- left_join(Data_source_3, enterobase_data, by = "Uberstrain")


## seem to be missing a column? is it meant to end in "average length"? it seems to be... Cant seem to find any ending in assembly status...


HC0s <- enterobase_data %>% filter(Lab_Contact == "M. Cummins (UTS, Sydney)") %>% select(starts_with("HC0")) %>% unique()

HC0s <- enterobase_data %>% filter(`HC0_(indistinguishable)` %in% HC0s$`HC0_(indistinguishable)`)

#Generate a list of HC0 groups present in more than one niche
HC0_overlap <- HC0s %>% group_by(`HC0_(indistinguishable)`, Source_Niche) %>% summarise(counts = n()) %>% filter(!is.na(Source_Niche)) %>% group_by(`HC0_(indistinguishable)`) %>% summarise(number_of_sources = n()) %>% arrange(desc(number_of_sources)) %>% filter(number_of_sources > 1)



HC2s <- enterobase_data %>% filter(Lab_Contact == "M. Cummins (UTS, Sydney)") %>% select(starts_with("HC2")) %>% unique()

HC2s <- enterobase_data %>% filter(HC2 %in% HC2s$HC2)

#Generate a list of HC2 groups present in more than one niche
HC2_overlap <- HC2s %>% group_by(HC2, Source_Niche) %>% summarise(counts = n()) %>% filter(!is.na(Source_Niche)) %>% group_by(HC2) %>% summarise(number_of_sources = n()) %>% arrange(desc(number_of_sources)) %>% filter(number_of_sources > 1)


HC5s <- enterobase_data %>% filter(Lab_Contact == "M. Cummins (UTS, Sydney)") %>% select(starts_with("HC5")) %>% unique()

HC5s <- enterobase_data %>% filter(HC2 %in% HC5s$HC5)

#Generate a list of HC5 groups present in more than one niche
HC5_overlap <- HC5s %>% group_by(HC5, Source_Niche) %>% summarise(counts = n()) %>% filter(!is.na(Source_Niche)) %>% group_by(HC5) %>% summarise(number_of_sources = n()) %>% arrange(desc(number_of_sources)) %>% filter(number_of_sources > 1)



