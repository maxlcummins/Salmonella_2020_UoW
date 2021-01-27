ColV_vs_AMR <- Metadata  %>% select(class_res_counts, `ColV (328/668)`) %>% rename(ColV = `ColV (328/668)`)

ColV_vs_AMR$class_res_counts <- gsub("classes_res_","",ColV_vs_AMR$class_res_counts)

ColV_vs_AMR$ColV <- gsub("ColV_","ColV ",ColV_vs_AMR$ColV)
ColV_vs_AMR$ColV <- gsub("pos","positive",ColV_vs_AMR$ColV)
ColV_vs_AMR$ColV <- gsub("neg","negative",ColV_vs_AMR$ColV)

ColV_vs_AMR$class_res_counts <- as.numeric(ColV_vs_AMR$class_res_counts)


ggplot(data = ColV_vs_AMR,
       aes(x = ColV, y = class_res_counts)) +
        geom_boxplot(fill = "steelblue") +
        labs(y = "Number of classes of AMR determinants",
             x = "ColV Carriage",
             title = "Breadth of antimicrobial resistance versus ColV status")


ColV_vs_AMR$MDR <- ColV_vs_AMR$class_res_counts

ColV_vs_AMR$MDR[ColV_vs_AMR$MDR < 3] <- 0
ColV_vs_AMR$MDR[ColV_vs_AMR$MDR >= 3] <- 1


ColV_vs_AMR %>% group_by(ColV, MDR) %>% summarise(counts = n())

colv_mdr<- cbind(c(99,229), c(49,291))

fisher.test(colv_mdr)

ggplot(data = ColV_vs_AMR,
       aes(x = ColV, y = MDR)) +
        geom_boxplot(fill = "steelblue") +
        labs(y = "Number of classes of AMR determinants",
             x = "ColV Carriage",
             title = "Breadth of antimicrobial resistance versus ColV status")

