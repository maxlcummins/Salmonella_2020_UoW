library(ggtree)
library(phytools)
library(tidyr)
library(readr)
library(magrittr)
library(dplyr)

#### Variable setting ####
# Path to wd()
workdir <-
# Path to serovar data
serovar_data <-
# Path to kraken data
kraken_data <-

# Define our 'not in' command
'%nin%' <- Negate('%in%')

#### Serovar processing ####
# Read in serovar data
serovars <- read_csv(serovar_data)

# Create a column listing strain source as determined by strain name
serovars$strain_source <- serovars$genome

# Set hospital sourced names to source "Human"
serovars$strain_source <-
        gsub("SIML.*|ISLHD.*", "Human", serovars$strain_source)

# Set seagull sourced names to source "Gull"
serovars$strain_source <-
        gsub("SG.*|W_.*|Seagull.*", "Gull", serovars$strain_source)

# Generate a simple spreadsheet with serovar data
serovars_simple <-
        serovars %>% select(genome,
                            strain_source,
                            qc_status,
                            cgmlst_subspecies,
                            serovar_cgmlst)


#### Kraken Processing ####
# Read in kraken data
kraken <- read_delim(kraken_data,
                     "\t",
                     escape_double = FALSE,
                     trim_ws = TRUE)

# Clean names of samples so that their path isn't listed in their name
kraken$name_file <- gsub("\\/.*\\/", "", kraken$name_file)

# Remove dupliate headers
kraken <- kraken %>% filter(name_file != "name_file")

# Select only genus level hits
kraken2 <- kraken %>% filter(grepl("^G$", R))

# Select only the best hit for genus
simple_genus <- kraken2[!duplicated(kraken2$name_file),]

# Select only species level hits
kraken3 <- kraken %>% filter(grepl("^S", R))

# Select only the best hit for species
simple_species <- kraken3[!duplicated(kraken3$name_file),]

# Combine our genus and species hits
IDs <- left_join(simple_genus, simple_species, by = "name_file")


#### Combining our kraken and serovar data ####
# Select and rename our columns of interest
IDs <-
        IDs %>% select(name_file, root.x, root.y, `100.00.x`, `100.00.y`) %>%
        rename(
                "Genus" = root.x,
                "Species" = root.y,
                "Perc_frags_genus" = `100.00.x`,
                "Perc_frags_species" = `100.00.y`,
                "genome" = name_file)

# Join kraken and serovar data
IDs <- left_join(IDs, serovars_simple, by = "name_file")

# Pull out our Salmonella
salmonellae <-
        IDs %>% filter(!grepl("WARNING|FAIL", qc_status), Genus == "Salmonella")
