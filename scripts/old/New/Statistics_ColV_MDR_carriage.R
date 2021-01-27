ColV_vs_AMR <- Metadata  %>% select(class_res_counts, `ColV (328/668)`, pMLST, working_name, plas_counts, ESBL_ress, `intI1 (98/668)`) %>% rename(ColV = `ColV (328/668)`)

plasmids <- df_small2 %>% select(starts_with("p"), intI1)

plasmids$working_name <- rownames(df_small4)

ColV_vs_AMR <- left_join(ColV_vs_AMR, plasmids, by = "working_name")

ColV_vs_AMR$pU1_51_B10[ColV_vs_AMR$pU1_51_B10 < 75] <- 0
ColV_vs_AMR$pU1_51_B10[ColV_vs_AMR$pU1_51_B10 > 1] <- 1

ColV_vs_AMR$pUTI89[ColV_vs_AMR$pUTI89 < 75] <- 0
ColV_vs_AMR$pUTI89[ColV_vs_AMR$pUTI89 > 1] <- 1

ColV_vs_AMR$pAPEC_O2_ColV[ColV_vs_AMR$pAPEC_O2_ColV < 75] <- 0
ColV_vs_AMR$pAPEC_O2_ColV[ColV_vs_AMR$pAPEC_O2_ColV > 1] <- 1

ColV_vs_AMR$pBCE049_1[ColV_vs_AMR$pBCE049_1 < 75] <- 0
ColV_vs_AMR$pBCE049_1[ColV_vs_AMR$pBCE049_1 > 1] <- 1

ColV_vs_AMR$pSF_088_nores[ColV_vs_AMR$pSF_088_nores < 75] <- 0
ColV_vs_AMR$pSF_088_nores[ColV_vs_AMR$pSF_088_nores > 1] <- 1

ColV_vs_AMR$No_F_plasmid <- ColV_vs_AMR$pMLST

ColV_vs_AMR$No_F_plasmid <- gsub("[A-Z].*","0", ColV_vs_AMR$No_F_plasmid)
ColV_vs_AMR$No_F_plasmid <- gsub("^$","1", ColV_vs_AMR$No_F_plasmid)

ColV_vs_AMR <- ColV_vs_AMR %>% rename(plas_counts = "plas_counts.x")

ColV_vs_AMR$class_res_counts <- gsub("classes_res_","",ColV_vs_AMR$class_res_counts)
ColV_vs_AMR$plas_counts <- gsub("count_plas_","",ColV_vs_AMR$plas_counts)

ColV_vs_AMR$ColV <- gsub("ColV_","ColV ",ColV_vs_AMR$ColV)
ColV_vs_AMR$ColV <- gsub("pos","positive",ColV_vs_AMR$ColV)
ColV_vs_AMR$ColV <- gsub("neg","negative",ColV_vs_AMR$ColV)

ColV_vs_AMR$ESBL_ress <- gsub("ESBL_neg","0",ColV_vs_AMR$ESBL_ress)
ColV_vs_AMR$ESBL_ress <- gsub("ESBL_pos","1",ColV_vs_AMR$ESBL_ress)

ColV_vs_AMR$class_res_counts <- as.numeric(ColV_vs_AMR$class_res_counts)
ColV_vs_AMR$plas_counts <- as.numeric(ColV_vs_AMR$plas_counts)


ggplot(data = ColV_vs_AMR,
       aes(x = ColV, y = class_res_counts)) +
        geom_boxplot(fill = "steelblue") +
        labs(y = "Number of classes of AMR determinants",
             x = "ColV Carriage",
             title = "Breadth of antimicrobial resistance versus ColV status")


ColV_vs_AMR$MDR <- ColV_vs_AMR$class_res_counts

ColV_vs_AMR$MDR[ColV_vs_AMR$MDR < 3] <- 0
ColV_vs_AMR$MDR[ColV_vs_AMR$MDR >= 3] <- 1


table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pUTI89, MDR) %>% summarise(counts = n())

pUTI89_mdr<- cbind(c(table_trait_group[1,3]
                   ,table_trait_group[2,3]),
                 c(table_trait_group[3,3],
                   table_trait_group[4,3])
                 )

pUTI89_mdr <- matrix(unlist(pUTI89_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pSF_088_nores, MDR) %>% summarise(counts = n())

pSF_088_nores_mdr<- cbind(c(table_trait_group[1,3]
                   ,table_trait_group[2,3]),
                 c(table_trait_group[3,3],
                   table_trait_group[4,3])
)

pSF_088_nores_mdr <- matrix(unlist(pSF_088_nores_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pAPEC_O2_ColV, MDR) %>% summarise(counts = n())

pAPEC_O2_ColV_mdr<- cbind(c(table_trait_group[1,3]
                   ,table_trait_group[2,3]),
                 c(table_trait_group[3,3],
                   table_trait_group[4,3])
)

pAPEC_O2_ColV_mdr <- matrix(unlist(pAPEC_O2_ColV_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pBCE049_1, MDR) %>% summarise(counts = n())

pBCE049_mdr<- cbind(c(table_trait_group[1,3]
                   ,table_trait_group[2,3]),
                 c(table_trait_group[3,3],
                   table_trait_group[4,3])
)

pBCE049_mdr <- matrix(unlist(pBCE049_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pU1_51_B10, MDR) %>% summarise(counts = n())

pU1_51_B10_mdr <- cbind(c(table_trait_group[1,3]
                   ,table_trait_group[2,3]),
                 c(table_trait_group[3,3],
                   table_trait_group[4,3])
)

pU1_51_B10_mdr <- matrix(unlist(pU1_51_B10_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(intI1, MDR) %>% summarise(counts = n())

intI1_mdr<- cbind(c(table_trait_group[1,3]
                   ,table_trait_group[2,3]),
                 c(table_trait_group[3,3],
                   table_trait_group[4,3])
)

intI1_mdr <- matrix(unlist(intI1_mdr), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(ColV, MDR) %>% summarise(counts = n())

colv_mdr<- cbind(c(table_trait_group[1,3]
                    ,table_trait_group[2,3]),
                  c(table_trait_group[3,3],
                    table_trait_group[4,3]))

colv_mdr <- matrix(unlist(colv_mdr), 2)


fisher.test(pSF_088_nores_mdr)
fisher.test(pUTI89_mdr)
fisher.test(pU1_51_B10_mdr)
fisher.test(pAPEC_O2_ColV_mdr)
fisher.test(pBCE049_mdr)

fisher.test(intI1_mdr)
fisher.test(colv_mdr)

fisher.test(colv_mdr)


#### Check association between plasmid replicons and ESBL carriage ####

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pUTI89, ESBL_ress) %>% summarise(counts = n())

pUTI89_esbl<- cbind(c(table_trait_group[1,3]
                     ,table_trait_group[2,3]),
                   c(table_trait_group[3,3],
                     table_trait_group[4,3])
)

pUTI89_esbl <- matrix(unlist(pUTI89_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pSF_088_nores, ESBL_ress) %>% summarise(counts = n())

pSF_088_nores_esbl<- cbind(c(table_trait_group[1,3]
                            ,table_trait_group[2,3]),
                          c(table_trait_group[3,3],
                            table_trait_group[4,3])
)

pSF_088_nores_esbl <- matrix(unlist(pSF_088_nores_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pAPEC_O2_ColV, ESBL_ress) %>% summarise(counts = n())

pAPEC_O2_ColV_esbl<- cbind(c(table_trait_group[1,3]
                            ,table_trait_group[2,3]),
                          c(table_trait_group[3,3],
                            table_trait_group[4,3])
)

pAPEC_O2_ColV_esbl <- matrix(unlist(pAPEC_O2_ColV_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pBCE049_1, ESBL_ress) %>% summarise(counts = n())

pBCE049_esbl<- cbind(c(table_trait_group[1,3]
                      ,table_trait_group[2,3]),
                    c(table_trait_group[3,3],
                      table_trait_group[4,3])
)

pBCE049_esbl <- matrix(unlist(pBCE049_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pU1_51_B10, ESBL_ress) %>% summarise(counts = n())

pU1_51_B10_esbl <- cbind(c(table_trait_group[1,3]
                          ,table_trait_group[2,3]),
                        c(table_trait_group[3,3],
                          table_trait_group[4,3])
)

pU1_51_B10_esbl <- matrix(unlist(pU1_51_B10_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(intI1, ESBL_ress) %>% summarise(counts = n())

intI1_esbl<- cbind(c(table_trait_group[1,3]
                    ,table_trait_group[2,3]),
                  c(table_trait_group[3,3],
                    table_trait_group[4,3])
)

intI1_esbl <- matrix(unlist(intI1_esbl), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(ColV, ESBL_ress) %>% summarise(counts = n())

colv_esbl<- cbind(c(table_trait_group[1,3]
                   ,table_trait_group[2,3]),
                 c(table_trait_group[3,3],
                   table_trait_group[4,3]))

colv_esbl <- matrix(unlist(colv_esbl), 2)


fisher.test(pSF_088_nores_esbl)
fisher.test(pUTI89_esbl)
fisher.test(pU1_51_B10_esbl)
fisher.test(pAPEC_O2_ColV_esbl)
fisher.test(pBCE049_esbl)

fisher.test(intI1_esbl)
fisher.test(colv_esbl)

fisher.test(colv_esbl)


#### Check association between plasmid replicons and intI1 carriage ####

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pUTI89, intI1) %>% summarise(counts = n())

pUTI89_intI1<- cbind(c(table_trait_group[1,3]
                      ,table_trait_group[2,3]),
                    c(table_trait_group[3,3],
                      table_trait_group[4,3])
)

pUTI89_intI1 <- matrix(unlist(pUTI89_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pSF_088_nores, intI1) %>% summarise(counts = n())

pSF_088_nores_intI1<- cbind(c(table_trait_group[1,3]
                             ,table_trait_group[2,3]),
                           c(table_trait_group[3,3],
                             table_trait_group[4,3])
)

pSF_088_nores_intI1 <- matrix(unlist(pSF_088_nores_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pAPEC_O2_ColV, intI1) %>% summarise(counts = n())

pAPEC_O2_ColV_intI1<- cbind(c(table_trait_group[1,3]
                             ,table_trait_group[2,3]),
                           c(table_trait_group[3,3],
                             table_trait_group[4,3])
)

pAPEC_O2_ColV_intI1 <- matrix(unlist(pAPEC_O2_ColV_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pBCE049_1, intI1) %>% summarise(counts = n())

pBCE049_intI1<- cbind(c(table_trait_group[1,3]
                       ,table_trait_group[2,3]),
                     c(table_trait_group[3,3],
                       table_trait_group[4,3])
)

pBCE049_intI1 <- matrix(unlist(pBCE049_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(pU1_51_B10, intI1) %>% summarise(counts = n())

pU1_51_B10_intI1 <- cbind(c(table_trait_group[1,3]
                           ,table_trait_group[2,3]),
                         c(table_trait_group[3,3],
                           table_trait_group[4,3])
)

pU1_51_B10_intI1 <- matrix(unlist(pU1_51_B10_intI1), 2)

table_trait_group<- ColV_vs_AMR %>% filter(grepl("[A-Z]", pMLST)) %>% group_by(ColV, intI1) %>% summarise(counts = n())

colv_intI1<- cbind(c(table_trait_group[1,3]
                    ,table_trait_group[2,3]),
                  c(table_trait_group[3,3],
                    table_trait_group[4,3]))

colv_intI1 <- matrix(unlist(colv_intI1), 2)


fisher.test(pSF_088_nores_intI1)
fisher.test(pUTI89_intI1)
fisher.test(pU1_51_B10_intI1)
fisher.test(pAPEC_O2_ColV_intI1)
fisher.test(pBCE049_intI1)

fisher.test(colv_intI1)

fisher.test(colv_intI1)




