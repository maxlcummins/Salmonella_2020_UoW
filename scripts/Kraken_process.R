library(readr)
library(magrittr)
library(dplyr)

#### Variable setting ####
# Path to serovar data
serovar_data <- "output/UoW_Salmonella/sistr/serovars.csv"
# Path to kraken data
kraken_data <- "output/UoW_Salmonella/kraken2/kraken2_report.txt"

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
simple_genus <- kraken2[!duplicated(kraken2$name_file), ]

# Select only species level hits
kraken3 <- kraken %>% filter(grepl("^S", R))

# Select only the best hit for species
simple_species <- kraken3[!duplicated(kraken3$name_file), ]

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
                "genome" = name_file
        )

# Join kraken and serovar data
IDs <- left_join(IDs, serovars_simple, by = "genome")

# Produce simple

# Pull out our Salmonella which passed QC
salmonellae_pass <-
        IDs %>% filter(!grepl("WARNING|FAIL", qc_status), Genus == "Salmonella")

# Pull out our non-Salmonella
non_salmonellae <-
        IDs %>% filter(Genus != "Salmonella")

# Pull out our Salmonella which failed QC or had warnings
fail_salmonellae <-
        IDs %>% filter(grepl("WARNING|FAIL", qc_status), Genus == "Salmonella")

# Create a table listing serovar counts
serovar_table <-
        salmonellae_pass %>% group_by(serovar_cgmlst, strain_source) %>% summarise(counts = n()) %>% arrange(desc(counts))

# Create a list of files for us to move to salmonella only folder
salmonella_names <- paste0("output/UoW_Salmonella/shovill/final_assemblies/", salmonellae_pass$genome,".fasta")

#Move these Salmonella genomes over
message("Found ",length(salmonella_names)," Salmonella which met qc controls")
message("Copying them to output/UoW_Salmonella/shovill/salmonella")

file.copy(from = salmonella_names,
          to = "output/UoW_Salmonella/shovill/salmonella",
          recursive = FALSE,
          copy.mode = TRUE)

if(!dir.exists("delims")){
        dir.create("delims")
}

#Write CSV table
write_csv(serovar_table, "delims/serovar_table.csv")

# Write Salmonella pass qc table
write_csv(salmonellae_pass, "delims/salmonellae_pass.csv")

# Write Salmonella pass qc table
write_csv(fail_salmonellae, "delims/salmonellae_fail.csv")

# Write Salmonella pass qc table
write_csv(non_salmonellae, "delims/non_salmonellae.csv")