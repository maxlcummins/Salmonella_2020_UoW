Metadata1 <- read_delim("metadata/AVC171-HC50.txt",
                        "\t",
                        escape_double = FALSE,
                        trim_ws = TRUE)

cgmlst <- read_delim("metadata/cgMLST.txt",
                     "\t",
                     escape_double = FALSE,
                     trim_ws = TRUE)

cgmlst <- cgmlst %>% select(Uberstrain, starts_with("HC"))

Metadata1 <- left_join(Metadata1, cgmlst)

Metadata1$District <- as.character(Metadata1$District)

Metadata2 <- read_delim("metadata/ST95_enterobase/metadata_subset.txt",
                        "\t",
                        escape_double = FALSE,
                        trim_ws = TRUE)

full_metadata <- full_join(Metadata1, Metadata2)

write_delim(full_metadata,
            "metadata/all_metadata",
            "\t")
