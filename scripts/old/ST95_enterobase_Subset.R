library(readr)
library(dplyr)
library(magrittr)

#Metadata sheet for all strains identified as ST95 by enterobase.
#Note there is also a file at the path below that lists strains not included in
#this database as the release date for these strains was greater than the access
#date (Feb 19 2020)
#"/Users/maxcummins/Dropbox/Doctorate/Manuscripts/AVC171/AVC171/metadata/ST95_enterobase/ST95_enterobase_20-7-20.txt"

base_path <- "~/Dropbox/Doctorate/Manuscripts/AVC171/AVC171/"

setwd(base_path)

metadata_path <- "metadata/ST95_enterobase/ST95_enterobase_20-7-20.txt"

cgMLST_path <- "metadata/ST95_enterobase/ST95_enterobase_20-7-20_cgMLST.txt"

meta <-
        read_delim(
                metadata_path,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE
        )

#save original colnames so we can reassign them later.
#this can prevent headaches later - best to change colnames in script to let us
#do what we need to do and then change them back at the end
original_colnames <- colnames(meta)

#Remove underscores from column names
colnames(meta) <- gsub(
        " ",
        "_",
        colnames(meta)
)

#used to identify the niche/type/details combinations where details and niche
#were "NA".
#meta %>%
#        filter(is.na(Source_Details) == TRUE & is.na(Source_Niche) == TRUE) %>%
#        group_by(Source_Niche, Source_Type, Source_Details) %>%
#        summarise(count = n()) %>%
#        arrange(desc(count)) %>%
#        View()

#Remove strains where Source_Details and Source_Niche are NA
meta <- meta %>%
        filter(is.na(Source_Details) == FALSE & is.na(Source_Niche) == FALSE)

#Remove strains where Collection Year is NA
meta <- meta %>%
        filter(is.na(Collection_Year) == FALSE)

#Remove strains where Country is NA
meta <- meta %>%
        filter(is.na(Country) == FALSE)

#used to identify the niche/type/details combinations where details and niche
#were "NA".
#These may be manually curated later as they may be salvagable from other
#metadata associated with the samples in SRA.
#Otherwise they may just be included with country as Unknown
#meta %>%
#        filter(is.na(Continent) == TRUE) %>%
#        group_by(Source_Niche, Source_Type, Source_Details) %>%
#        summarise(count = n()) %>%
#        arrange(desc(count)) %>%
#        View()

#Used to identify strains that are not assembled and therefore cant be downloaded as assemblies
#None were identified as "Queued", only as "Assembled", therefore no steps
#needed to be taken here
#meta %>%
#        group_by(`Status(Assembly_stats)`) %>%
#        summarise(count = n()) %>%
#        arrange(desc(count)) %>%
#        View()

#change column name with illegal characters
colnames(meta) <-
        gsub("Contig_Number.*",
             "contig_num_greater_than_200_bp",
             colnames(meta))

#Remove strains with dodgy assembly stats
meta <- meta %>% filter(
        `Length` >= 4500000,
        `Length` <= 6500000,
        `Coverage` >= 20,
        contig_num_greater_than_200_bp <= 600,
        Low_Quality_Bases <= 50000
        #Low quality bases threshold informed by:
        #plot(meta$`Low Quality Bases`)
        #100 strains lost by implementing this cutoff, but these strains with high low quality base values are likely to mess up our trees
)

#reassign original colnames of meta
colnames(meta) <- original_colnames

#read in cgMLST data
cgMLST <-
        read_delim(
                cgMLST_path,
                "\t",
                escape_double = FALSE,
                trim_ws = TRUE
        )

#remove unwanted columns from cgMLST table
cgMLST <- cgMLST[,c(1,38:50)]

#Join metadata and cgMLST tables
meta <- left_join(meta, cgMLST)

#Remove strains where HC50 == 1106
meta <- meta %>%
        filter(HC50 != 1106)

#unhash the following line if a directory called delims doesn't exist
if(!dir.exists("delims")){
        dir.create("delims")
}
        
#Write our metadata sheet to file.
if(!file.exists("delims/metadata_subset.txt")){
        write_delim(x = meta,
                    path = "delims/metadata_subset.txt",
                    delim = "\t"
        )
}



