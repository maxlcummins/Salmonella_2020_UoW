library(ggtree)


#Provides clade numbers to colour by.
#If you dont know what the node labels are use the line below under "get node labels"
#You may need to change the sizing..
#tree2 <- groupClade(tree, c(87, 86))

#source("scripts/New/Genotype_Metadata_processing.R")

plas_names <- df_small %>% select(starts_with("p"),-plas_counts) %>% colnames()

plas_cuts <- c("0-50", "51-60", "61-70", "71-80", "81-90", "90-100")

plas_cells <- c()

for(i in plas_names){
        for(j in plas_cuts){
                plas_cells <- c(plas_cells, paste(i,j, sep = "_"))
}
}

plasmid1_col_gen <- colorRampPalette(c("white", "red"))
plasmid2_col_gen <- colorRampPalette(c("white", "black"))
plasmid3_col_gen <- colorRampPalette(c("white", "orange"))
plasmid4_col_gen <- colorRampPalette(c("white", "purple"))
plasmid5_col_gen <- colorRampPalette(c("white", "green"))


plasmid1_cols <- plasmid1_col_gen(6)
plasmid2_cols <- plasmid2_col_gen(6)
plasmid3_cols <- plasmid3_col_gen(6)
plasmid4_cols <- plasmid4_col_gen(6)
plasmid5_cols <- plasmid5_col_gen(6)

plot(x = 1:8, y = rep(1, 8), col = class_res_count_cols, pch = 19)

plas_cols <- c(plasmid1_cols, plasmid2_cols, plasmid3_cols, plasmid4_cols, plasmid5_cols)

col_els <- unique(c(Metadata$HC200_Other,
                    Metadata$class_res_counts,
                    Metadata$ESBL_ress,
                    df$`ColV (328/668)`,
                    df_small4$Pathogen,
                    plas_cells,
                    df_small4$Revised_Source_Niche,
                    df_small4$amr_counts,
                    df_small4$plas_counts,
                    df_small4$intI1
                    )) %>% 
        sort() %>%
        as.data.frame() %>%
        rename(variable = ".") 

col_els
                    

genotype_cols <- c(
        "white", # 0 - no gene present
        "black" # 1 - gene present
)

colV_cols <- c(
        "grey40", # ColV_neg - ColV present (Liu criteria)
        "red" # ColV_pos - ColV present (Liu criteria)
        )



amr_col_gen <- colorRampPalette(c("white", "blue", "red"))

amr_count_cols <- amr_col_gen(13)

class_res_count_cols <- amr_col_gen(8)

ESBL_res_cols <- c("white", # ESBL NEG
                   "black"  # ESBL POS
                   )

plas_count_cols <- c("gray100", #   plas_count - res_count_0
                  "gray75", #   plas_count - res_count_1
                  "gray50", #   plas_count - res_count_2
                  "gray25", #   plas_count - res_count_3
                  "gray1" #   plas_count - res_count_4
                          )

HC200_cols <- c(        "#60c757"	,	#	1104	#	HC200		
                       "#47037a"	,	#	1108	#	HC200		
                       "#726c00"	,	#	1592	#	HC200		
                       "#ff6ed0"	,	#	4252	#	HC200		
                       "#ff9954"	,	#	44	#	HC200		
                       "#525799"	,	#	55	#	HC200		
                       "#d62678"	,	#	6624	#	HC200		
                       "#d4abff"		#	8655	#	HC200
                       )

Pathogen_cols <- c("green",  # Pathogen status        # Environmental
                   "black", # Pathogen status         # Flora
                   "grey", # Pathogen status         # Other
                   "pink", # Pathogen status         # Raw Chicken
                   "red", # Pathogen status         # Systemic
                   "yellow") # Pathogen status         # Urine
                   

source_cols <- c(
        "brown"	,	#	Bovine	#	Source		
        "gold2"	,	#	Canine	#	Source		
        "springgreen4"	,	#	Environment	#	Source		
        "cyan2"	,	#	Human	#	Source		
        "black"	,	#	Other	#	Source		
        "mediumorchid1" #	Poultry	#	Source
)

col_els

col_els$colours <- c(Pathogen_cols,
                     genotype_cols,
                     HC200_cols,
                     class_res_count_cols,
                     colV_cols,
                     plas_count_cols,
                     ESBL_res_cols,
                     "black", #Other HC                    
                     plas_cols,
                     source_cols
)

cols_needed <- unique(c(df_small4$Pathogen,
                        Metadata$HC200_Other,
                    df$`ColV (328/668)`,
                    df_small4$plasmid_map,
                    df_small4$Revised_Source_Niche,
                    df_small4$class_res_counts,
                    df_small4$ESBL_ress,
                    df_small4$plas_counts,
                    df_small4$intI1
)) %>% 
        sort() %>%
        as.data.frame() %>%
        rename(variable = ".")

col_els2 <- left_join(cols_needed, col_els)




p2 <- ggtree(tree, layout = 'circular',
             open.angle = 0) %<+%
        Metadata +
        geom_tiplab(size = 0.5,
                    align = TRUE,
                    linesize = 0.15,
                    aes(color = as.factor(`ColV (328/668)`))
        ) +
        geom_tippoint(size = 0.2, aes(colour=as.factor(HC200_Other))
        )

phylogenetic_tree <- gheatmap(
        p = p2,
        data = df_small4,
        colnames_offset_y = -0.1,
        font.size = 1.5,
        hjust = 0,
        #colnames_position = "top",
        colnames = FALSE,
        #colnames_angle = 90,
        offset = 0.0025,
        width = 0.20,
        color = 'grey'
) + 
        #ggplot2::ylim(NA, 700) #+
        theme(legend.position = "bottom",
              legend.box = "vertical",
              legend.key.size = unit(1, "mm"),
              legend.title = element_blank(),
              legend.text = element_text(size = 8),
        ) +
        scale_fill_manual(
                aesthetics = c("colour", "fill"),
                values = col_els2$colours,
                na.value = 'grey'
        ) +
        geom_tiplab2(
                aes(label = pMLST),#, colour = pMLST),
                align = TRUE,
                linetype = NULL,
                offset = 0.0095,
                size = 0.8
        ) +
        geom_tiplab2(
                aes(label = H_type),#, colour = pMLST),
                align = TRUE,
                linetype = NULL,
                offset = 0.01225,
                size = 0.8
        ) +
        geom_tiplab2(
                aes(label = O_type),#, colour = pMLST),
                align = TRUE,
                linetype = NULL,
                offset = 0.01325,
                size = 0.8
        ) +
        geom_tiplab2(
                aes(label = fimH_type),#, colour = pMLST),
                align = TRUE,
                linetype = NULL,
                offset = 0.0150,
                size = 0.8
        )

require(cowplot)
leg1 <- get_legend(phylogenetic_tree)

plot(leg1)

phylogenetic_tree + theme(legend.position="none")


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
                values = colval2,
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

listy2 <- read_csv("delims/pEC14_114_plasmid_coverage.csv")

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
