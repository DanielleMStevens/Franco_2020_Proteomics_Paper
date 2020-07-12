#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 2/9/2020
# Script Purpose: Plotting Signifcant Go-terms from proteomics data to assess redundancy; plot final heatmap
# Inputs Necessary: Excel file with signficant go terms, protein id, child and parent go-terms, z-scores
# Outputs: Visualization of comparsion between go-terms called
#-----------------------------------------------------------------------------------------------

#library packages need to load
library(RColorBrewer)
library(extrafont)
font_import()
loadfonts(device = "win")
#windowsFonts()
library(sysfonts)
library(showtext)
library(readxl)
library(rlang)
library(vctrs)
library(tidyverse)
library(plyr)
library(magrittr)
library(dplyr)
library(xlsx)
library(cowplot)
library(grid)
library(forcats)
library(stringr)
library(caTools)
library(gplots)
library(remotes)
library(reshape2)
library(gridExtra)
library(circlize)
library(rafalib)
library(scales)
library(ggplot2)
library(Rcpp)
library(devtools)
library(ComplexHeatmap)

# NOTE: I have had issues sometimes loading the complex heatmap package (not sure why), 
# Try one of the many ways to download the package and if still running into troubles,
# consult google or contact me

#install_github("jokergoo/ComplexHeatmap")
#devtools::install_github("jokergoo/ComplexHeatmap")
#if (!requireNamespace("BiocManager", quietly=TRUE))
  #install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")


######################################################################

#preliminary plots to filter through redundent go-terms 

######################################################################

#raw go-terms file to process - load in file
file_to_open <- file.choose() #choose the file: v2_DEproteins_GOterms_09_20_2019.xlsx
upregulated_bio <- as.data.frame(read_excel(file_to_open, sheet=2, col_names = TRUE))
downregulated_bio <- as.data.frame(read_excel(file_to_open, sheet=3, col_names = TRUE))
upregulated_mol <- as.data.frame(read_excel(file_to_open, sheet=1, col_names = TRUE))
downregulated_mol <- as.data.frame(read_excel(file_to_open, sheet=4, col_names = TRUE))

######################################################################

#convert \t into underscores

######################################################################

for (i in 1:nrow(upregulated_bio)){
  upregulated_bio[i,4] <- gsub("\t","_", upregulated_bio[i,4])
}

for (i in 1:nrow(downregulated_bio)){
  downregulated_bio[i,4] <- gsub("\t","_", downregulated_bio[i,4])
}

for (i in 1:nrow(upregulated_mol)){
  upregulated_mol[i,4] <- gsub("\t","_", upregulated_mol[i,4])
}

for (i in 1:nrow(downregulated_mol)){
  downregulated_mol[i,4] <- gsub("\t","_", downregulated_mol[i,4])
}

######################################################################

# process and rewrite new file (excel) with go terms by number of hits by phytozome identifier

######################################################################

counts_upregulated_mol <- as.data.frame(table(upregulated_mol$`Protein Phytozome identifier`))
counts_upregulated_mol <- counts_upregulated_mol[order(table(upregulated_mol$`Protein Phytozome identifier`), decreasing = T),]
counts_upregulated_mol$Var1 <- as.character(counts_upregulated_mol$Var1)
total_names_up_mol <- upregulated_mol[order(match(upregulated_mol$`Protein Phytozome identifier`, counts_upregulated_mol$Var1)),]

counts_upregulated_bio <- as.data.frame(table(upregulated_bio$`Protein Phytozome identifier`))
counts_upregulated_bio <- counts_upregulated_bio[order(table(upregulated_bio$`Protein Phytozome identifier`), decreasing = T),]
counts_upregulated_bio$Var1 <- as.character(counts_upregulated_bio$Var1)
total_names_up_bio <- upregulated_bio[order(match(upregulated_bio$`Protein Phytozome identifier`, counts_upregulated_bio$Var1)),]

counts_downregulated_mol <- as.data.frame(table(downregulated_mol$`Protein Phytozome identifier`))
counts_downregulated_mol <- counts_downregulated_mol[order(table(downregulated_mol$`Protein Phytozome identifier`), decreasing = T),]
counts_downregulated_mol$Var1 <- as.character(counts_downregulated_mol$Var1)
total_names_down_mol <- downregulated_mol[order(match(downregulated_mol$`Protein Phytozome identifier`, counts_downregulated_mol$Var1)),]

counts_downregulated_bio <- as.data.frame(table(downregulated_bio$`Protein Phytozome identifier`))
counts_downregulated_bio <- counts_downregulated_bio[order(table(downregulated_bio$`Protein Phytozome identifier`), decreasing = T),]
counts_downregulated_bio$Var1 <- as.character(counts_downregulated_bio$Var1)
total_names_down_bio <- downregulated_bio[order(match(downregulated_bio$`Protein Phytozome identifier`, counts_downregulated_bio$Var1)),]

write.xlsx(total_names_up_mol, "GO_terms_decreasing_order_abundance.xlsx", sheetName = "up_mol", col.names = T, row.names = F)
write.xlsx(total_names_up_bio, "GO_terms_decreasing_order_abundance.xlsx", sheetName = "up_bio", col.names = T, row.names = F, append = T)
write.xlsx(total_names_down_mol, "GO_terms_decreasing_order_abundance.xlsx", sheetName = "down_mol", col.names = T, row.names = F, append = T)
write.xlsx(total_names_down_bio, "GO_terms_decreasing_order_abundance.xlsx", sheetName = "down_bio", col.names = T, row.names = F, append = T)

######################################################################

# plot total counts for each significant protein detected

######################################################################

#################### upregulated bio function go terms

fasta_names_upregulated_bio <- upregulated_bio[order(fct_infreq(upregulated_bio$`Protein Phytozome identifier`)),c(1,4)]
fasta_names_upregulated_bio <- dplyr::distinct(fasta_names_upregulated_bio)



upregulated_bio_parent <- ggplot(upregulated_bio, 
               aes(x = fct_infreq(as.factor(upregulated_bio$`Protein Phytozome identifier`)),
                   y = ..count.., fill = upregulated_bio$`Parent GO Term Name`)) +
  geom_bar (colour = "black") + 
  theme_bw() +
  scale_x_discrete(labels = fasta_names_upregulated_bio$`Fasta header`) +
  theme(axis.text.x = element_blank(), legend.key.size=unit(6,"point"),legend.position="bottom") +
  scale_y_continuous(breaks = seq(0,10, by = 1)) +
  guides(fill = guide_legend(title = "Parent GO Terms",nrow =4)) +
  xlab("")

upregulated_bio_child <- ggplot(upregulated_bio, 
                aes(x = fct_infreq(as.factor(upregulated_bio$`Protein Phytozome identifier`)), 
                    y = ..count.., fill = upregulated_bio$`Ontology Term Name`)) +
  geom_bar (colour = "black") + 
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 45, hjust = 1),legend.key.size=unit(6,"point"),legend.position="bottom") +
  scale_x_discrete(labels = fasta_names_upregulated_bio$`Fasta header`) +
  scale_y_continuous(breaks = seq(0,10, by = 1)) +
  guides(fill = guide_legend(title = "Child GO Terms", nrow = 4)) +
  xlab("Protein Name")


upregulated_bio_counts <- grid.arrange(upregulated_bio_parent, upregulated_bio_child, ncol=1,
             top = "Significantly Expressed Upregulated Proteins - GO Terms - Biological Function",
             heights = c(0.3,1))

ggsave(filename = "Upregulated_bio_counts.jpg", plot = upregulated_bio_counts, width = 11, height = 6, units = "in", dpi = 400)
dev.off()


#################################################################

# This plot is a subset of the data for the suplementary figure!!!

#################################################################

small_plot_for_supplementary_figure_fasta_name <- upregulated_bio[order(fct_infreq(upregulated_bio$`Protein Phytozome identifier`)),c(1,4)]
small_plot_for_supplementary_figure <- upregulated_bio[order(fct_infreq(upregulated_bio$`Protein Phytozome identifier`)),]

small_plot_for_supplementary_figure_fasta_name <- small_plot_for_supplementary_figure_fasta_name[1:19,]
small_plot_for_supplementary_figure <- small_plot_for_supplementary_figure[1:19,]
small_plot_for_supplementary_figure_fasta_name <- dplyr::distinct(small_plot_for_supplementary_figure_fasta_name)

ggplot(small_plot_for_supplementary_figure, 
       aes(x = fct_infreq(as.factor(small_plot_for_supplementary_figure$`Protein Phytozome identifier`)), 
           y = ..count.., fill = small_plot_for_supplementary_figure$`Ontology Term Name`)) +
  scale_x_discrete(label = small_plot_for_supplementary_figure_fasta_name$`Fasta header`) +
  geom_bar (colour = "black", width = 0.6, position = position_stack(vjust = 1), size = 0.3) + 
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 45, hjust = 1, size = 14, family = "Arial", color= "black"),
        axis.text.y = element_text(color = 'black', size = 14, family = "Arial"), 
          axis.title.y = element_text(color='black', size =16, family = "Arial"),
        axis.title.x = element_text(color = 'black', size = 16, family = "Arial"),
        aspect.ratio = 1.4,
        legend.key.size=unit(10,"point"), legend.text = element_text(size = 14, color = "black"), 
        legend.title = element_text(size = 16, color = "black"),
        legend.position = "right") +
  scale_y_continuous(breaks = seq(0,10, by = 1)) +
  guides(fill = guide_legend(title = "GO Terms", nrow = 4)) +
  xlab("\nProtein Annotation") 


#######################################################

#################### downregulated bio function go terms

fasta_names_downregulated_bio <- downregulated_bio[order(fct_infreq(downregulated_bio$`Protein Phytozome identifier`)),c(1,4)]
fasta_names_downregulated_bio <- distinct(fasta_names_downregulated_bio)

downregulated_bio_parent <- ggplot(downregulated_bio, 
                                 aes(x = fct_infreq(as.factor(downregulated_bio$`Protein Phytozome identifier`)), 
                                     y = ..count.., fill = downregulated_bio$`Parent GO Term Name`)) +
  geom_bar (colour = "black") + 
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.key.size=unit(6,"point"),legend.position="bottom") +
  scale_x_discrete(labels = fasta_names_downregulated_bio$`Fasta header`) +
  scale_y_continuous(breaks = seq(0,80, by = 10)) +
  guides(fill = guide_legend(title = "Parent GO Terms")) +
  xlab("")

downregulated_bio_child <- ggplot(downregulated_bio, 
                                aes(x = fct_infreq(as.factor(downregulated_bio$`Protein Phytozome identifier`)), 
                                    y = ..count.., fill = downregulated_bio$`Ontology Term Name`)) +
  geom_bar (colour = "black") + 
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 45, hjust = 1), legend.key.size=unit(6,"point"),legend.position="bottom") +
  scale_x_discrete(labels = fasta_names_downregulated_bio$`Fasta header`) +
  scale_y_continuous(breaks = seq(0,80, by = 10)) +
  guides(fill = guide_legend(title = "Child GO Terms")) +
  xlab("Protein Name")


downregulated_bio_counts <- grid.arrange(downregulated_bio_parent, downregulated_bio_child, ncol=1,
                                       top = "Significantly Expressed Downregulated Proteins - GO Terms - Biological Function",
                                       heights = c(0.4,1))

ggsave(filename = "Downregulated_bio_counts.jpg", plot = downregulated_bio_counts, width = 40, height = 26, units = "in", dpi = 400)

dev.off()

#########################################################

#################### upregulated mol function go terms

fasta_names_upregulated_mol <- upregulated_mol[order(fct_infreq(upregulated_mol$`Protein Phytozome identifier`)),c(1,4)]
fasta_names_upregulated_mol <- distinct(fasta_names_upregulated_mol)

upregulated_mol_parent <- ggplot(upregulated_mol, 
                                   aes(x = fct_infreq(as.factor(upregulated_mol$`Protein Phytozome identifier`)), 
                                       y = ..count.., fill = upregulated_mol$`Parent GO Term Name`)) +
  geom_bar (colour = "black") + 
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.key.size=unit(6,"point"),legend.position="bottom") +
  scale_x_discrete(labels = fasta_names_upregulated_mol$`Fasta header`) +
  scale_y_continuous(breaks = seq(0,10, by = 1)) +
  guides(fill = guide_legend(title = "Parent GO Terms")) +
  xlab("")

upregulated_mol_child <- ggplot(upregulated_mol, 
                                  aes(x = fct_infreq(as.factor(upregulated_mol$`Protein Phytozome identifier`)), 
                                      y = ..count.., fill = upregulated_mol$`Ontology Term Name`)) +
  geom_bar (colour = "black") + 
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 45, hjust = 1), legend.key.size=unit(6,"point"),legend.position="bottom") +
  scale_x_discrete(labels = fasta_names_upregulated_mol$`Fasta header`) +
  scale_y_continuous(breaks = seq(0,10, by = 1)) +
  guides(fill = guide_legend(title = "Child GO Terms")) +
  xlab("Protein Name")


upregulated_mol_counts <- grid.arrange(upregulated_mol_parent, upregulated_mol_child, ncol=1,
                                         top = "Significantly Expressed Upregulated Proteins - GO Terms - Molecular Function",
                                         heights = c(0.4,1))

ggsave(filename = "Upregulated_mol_counts.jpg", plot = upregulated_mol_counts, width = 32, height = 20, units = "in", dpi = 400)
dev.off()

######################################################

#################### upregulated mol function go terms

fasta_names_downregulated_mol <- downregulated_mol[order(fct_infreq(downregulated_mol$`Protein Phytozome identifier`)),c(1,4)]
fasta_names_downregulated_mol <- distinct(fasta_names_downregulated_mol)

downregulated_mol_parent <- ggplot(downregulated_mol, 
                                 aes(x = fct_infreq(as.factor(downregulated_mol$`Protein Phytozome identifier`)), 
                                     y = ..count.., fill = downregulated_mol$`Parent Ontology Term`)) +
  geom_bar (colour = "black") + 
  theme_bw() +
  theme(axis.text.x = element_blank(), legend.key.size=unit(6,"point"),legend.position="bottom") +
  scale_x_discrete(labels = fasta_names_downregulated_mol$`Fasta header`) +
  scale_y_continuous(breaks = seq(0,20, by = 2)) +
  guides(fill = guide_legend(title = "Parent GO Terms")) +
  xlab("")

downregulated_mol_child <- ggplot(downregulated_mol, 
                                aes(x = fct_infreq(as.factor(downregulated_mol$`Protein Phytozome identifier`)), 
                                    y = ..count.., fill = downregulated_mol$`Ontology Term Name`)) +
  geom_bar (colour = "black") + 
  theme_bw() +
  theme(axis.text.x =  element_text(angle = 45, hjust = 1), legend.key.size=unit(6,"point"),legend.position="bottom") +
  scale_x_discrete(labels = fasta_names_downregulated_mol$`Fasta header`) +
  scale_y_continuous(breaks = seq(0,20, by = 2)) +
  guides(fill = guide_legend(title = "Child GO Terms")) +
  xlab("Protein Name")


downregulated_mol_counts <- grid.arrange(downregulated_mol_parent, downregulated_mol_child, ncol=1,
                                       top = "Significantly Expressed Downregulated Proteins - GO Terms - Molecular Function",
                                       heights = c(0.4,1))

ggsave(filename = "Downregulated_mol_counts.jpg", plot = downregulated_mol_counts, width = 32, height = 20, units = "in", dpi = 400)
dev.off()

######################################################################

#filtered go terms

######################################################################

#choose go-terms file to process
file_to_open <- file.choose() #choose Filtered_GO_terms_to_plot.xlsx

upregulated_mol_filtered <- as.data.frame(read_excel(file_to_open, sheet = 1, col_names = TRUE))
upregulated_bio_filtered <- as.data.frame(read_excel(file_to_open, sheet =  2, col_names = TRUE))
downregulated_bio_filtered <- as.data.frame(read_excel(file_to_open, sheet = 4, col_names = TRUE))
downregulated_mol_filtered <- as.data.frame(read_excel(file_to_open, sheet = 3, col_names = TRUE))

#up mol filtered terms

all_filtered_up_mol_terms <- data.frame("Protein Phytozome identifier" = character(), "fasta header" = character(),
                                      "GO Term Relationship" = character(), "Ontology Term Name" = character(), 
                                      "Ontology Term ID" = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(upregulated_mol_filtered)){
  if(is.na(upregulated_mol_filtered[i,6]) == TRUE){
    temp_data_frame <- data.frame(upregulated_mol_filtered[i,1],upregulated_mol_filtered[i,4],"Parent GO",upregulated_mol_filtered[i,8],upregulated_mol_filtered[i,7])
    colnames(temp_data_frame) <- colnames(all_filtered_up_mol_terms)
    all_filtered_up_mol_terms <- rbind(all_filtered_up_mol_terms, temp_data_frame)
  }
  if(is.na(upregulated_mol_filtered[i,8]) == TRUE){
    temp_data_frame <- data.frame(upregulated_mol_filtered[i,1],upregulated_mol_filtered[i,4],"Child GO",upregulated_mol_filtered[i,6],upregulated_mol_filtered[i,5])
    colnames(temp_data_frame) <- colnames(all_filtered_up_mol_terms)
    all_filtered_up_mol_terms <- rbind(all_filtered_up_mol_terms, temp_data_frame)
  }
}
all_filtered_up_mol_terms$Protein.Phytozome.identifier <- as.character(all_filtered_up_mol_terms$Protein.Phytozome.identifier)

#up bio filtered terms

all_filtered_up_bio_terms <- data.frame("Protein Phytozome identifier" = character(), "fasta header" = character(),
                                        "GO Term Relationship" = character(), "Ontology Term Name" = character(), 
                                        "Ontology Term ID" = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(upregulated_bio_filtered)){
  if(is.na(upregulated_bio_filtered[i,6]) == TRUE){
    temp_data_frame <- data.frame(upregulated_bio_filtered[i,1],upregulated_bio_filtered[i,4],"Parent GO",upregulated_bio_filtered[i,8],upregulated_bio_filtered[i,7])
    colnames(temp_data_frame) <- colnames(all_filtered_up_bio_terms)
    all_filtered_up_bio_terms <- rbind(all_filtered_up_bio_terms, temp_data_frame)
  }
  if(is.na(upregulated_bio_filtered[i,8]) == TRUE){
    temp_data_frame <- data.frame(upregulated_bio_filtered[i,1],upregulated_bio_filtered[i,4],"Child GO",upregulated_bio_filtered[i,6],upregulated_bio_filtered[i,5])
    colnames(temp_data_frame) <- colnames(all_filtered_up_bio_terms)
    all_filtered_up_bio_terms <- rbind(all_filtered_up_bio_terms, temp_data_frame)
  }
}

#down mol filtered terms

all_filtered_down_mol_terms <- data.frame("Protein Phytozome identifier" = character(), "fasta header" = character(),
                                        "GO Term Relationship" = character(), "Ontology Term Name" = character(), 
                                        "Ontology Term ID" = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(downregulated_mol_filtered)){
  if(is.na(downregulated_mol_filtered[i,6]) == TRUE){
    temp_data_frame <- data.frame(downregulated_mol_filtered[i,1],downregulated_mol_filtered[i,4],"Parent GO",downregulated_mol_filtered[i,8],downregulated_mol_filtered[i,7])
    colnames(temp_data_frame) <- colnames(all_filtered_down_mol_terms)
    all_filtered_down_mol_terms <- rbind(all_filtered_down_mol_terms, temp_data_frame)
  }
  if(is.na(downregulated_mol_filtered[i,8]) == TRUE){
    temp_data_frame <- data.frame(downregulated_mol_filtered[i,1],downregulated_mol_filtered[i,4],"Child GO",downregulated_mol_filtered[i,6],downregulated_mol_filtered[i,5])
    colnames(temp_data_frame) <- colnames(all_filtered_up_mol_terms)
    all_filtered_down_mol_terms <- rbind(all_filtered_down_mol_terms, temp_data_frame)
  }
}

#down bio filtered terms

all_filtered_down_bio_terms <- data.frame("Protein Phytozome identifier" = character(), "fasta header" = character(),
                                        "GO Term Relationship" = character(), "Ontology Term Name" = character(), 
                                        "Ontology Term ID" = character(), stringsAsFactors = FALSE)
for (i in 1:nrow(downregulated_bio_filtered)){
  if(is.na(downregulated_bio_filtered[i,6]) == TRUE){
    temp_data_frame <- data.frame(downregulated_bio_filtered[i,1],downregulated_bio_filtered[i,4],"Parent GO",downregulated_bio_filtered[i,8],downregulated_bio_filtered[i,7])
    colnames(temp_data_frame) <- colnames(all_filtered_down_bio_terms)
    all_filtered_down_bio_terms <- rbind(all_filtered_down_bio_terms, temp_data_frame)
  }
  if(is.na(downregulated_bio_filtered[i,8]) == TRUE){
    temp_data_frame <- data.frame(downregulated_bio_filtered[i,1],downregulated_bio_filtered[i,4],"Child GO",downregulated_bio_filtered[i,6],downregulated_bio_filtered[i,5])
    colnames(temp_data_frame) <- colnames(all_filtered_down_bio_terms)
    all_filtered_down_bio_terms <- rbind(all_filtered_down_bio_terms, temp_data_frame)
  }
}


####################################################

#plotting molecular Function terms

####################################################

replot_mol_terms <- data.frame("Protein Identifier" = character(), "fasta header" = character(), "Regulation" = character(), "Ontology Term Name" = character(), "Ontology ID" = character())
for (i in 1:nrow(all_filtered_up_mol_terms)){
  temp_dataframe <- data.frame(all_filtered_up_mol_terms[i,1], all_filtered_up_mol_terms[i,2], "Upregulated", all_filtered_up_mol_terms[i,4], all_filtered_up_mol_terms[i,5])
  colnames(temp_dataframe) <- colnames(replot_mol_terms)
  replot_mol_terms <- rbind(replot_mol_terms, temp_dataframe)
}
for (i in 1:nrow(all_filtered_down_mol_terms)){
  temp_dataframe <- data.frame(all_filtered_down_mol_terms[i,1], all_filtered_down_mol_terms[i,2], "Downregulated", all_filtered_down_mol_terms[i,4], all_filtered_down_mol_terms[i,5])
  colnames(temp_dataframe) <- colnames(replot_mol_terms)
  replot_mol_terms <- rbind(replot_mol_terms, temp_dataframe)
}

#order_list_replot_mol <- replot_mol_terms %>%  dplyr::group_by(replot_mol_terms$Ontology.Term.Name, replot_mol_terms$Regulation) %>% 
 # dplyr::summarise(n=n()) %>% arrange(desc(n))

manual_ordering_based_on_order_list_and_personal_preferences_mol <- c("cysteine-type peptidase activity",
                                                                      "serine-type carboxypeptidase activity",
                                                                      "hydrolase activity",
                                                                      "oxidoreductase activity",
                                                                      "peptidase activity",
                                                                      "hydrolase activity, hydrolyzing O-glycosyl compounds",
                                                                      "catalytic activity",

                                                                      "serine-type endopeptidase activity",
                                                                      "threonine-type endopeptidase activity",
                                                                      "oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor",
                                                                      "phosphoglycerate mutase activity",                                                    
                                                                      "intramolecular transferase activity",
                                                                      "glutamate decarboxylase activity",                                                           
                                                                      "O-methyltransferase activity",   
                                                                      "structural molecule activity",
                                                                      "aspartic-type endopeptidase activity",
                                                                      "endopeptidase activity",
                                                                      "magnesium ion binding",
                                                                      "oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor",
                                                                      "glycolytic process",
                                                                      "ligase activity",
                                                                      "structural constituent of ribosome",
                                                                      
                                                                      "hydro-lyase activity",
                                                                      "carboxypeptidase activity",
                                                                      "hydrolase activity, acting on glycosyl bonds",
                                                                      "peptidase activity, acting on L-amino acid peptides",
                                                                      "peroxidase activity")

                                              

replot_mol_terms$Ontology.Term.Name <- factor(replot_mol_terms$Ontology.Term.Name, levels = manual_ordering_based_on_order_list_and_personal_preferences_mol)


sysfonts::font_add("Arial", "arial.ttf")

color_plot <- c("Upregulated" = "#BA3241", "Downregulated" = "#3C7AB2")

replot_mol_group1 <- ggplot(replot_mol_terms, 
                     aes(x = factor(replot_mol_terms$Ontology.Term.Name),
                         y=..count.., fill = replot_mol_terms$Regulation)) +
  geom_bar(colour = "black") + 
  theme_bw() +
  theme(axis.text.x =  element_text(hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), axis.title.x = element_text(color = "black"),
        text = element_text(size=9, family = "Arial",colour = "black"), legend.position="none") +
  xlab("") + 
  scale_fill_manual(values = color_plot) +
  ylab("Counts") +
  coord_flip(ylim = c(0,50))

replot_mol_group2 <- ggplot(replot_mol_terms, 
                            aes(x = factor(replot_mol_terms$Ontology.Term.Name),
                                y=..count.., fill = replot_mol_terms$Regulation)) +
  geom_bar(colour = "black") + 
  theme_bw() +
  theme(axis.text.y =  element_blank(), axis.text.x = element_text(color = "black"), axis.title.x = element_text(color = "black"),
        text = element_text(size=9, family = "Arial",colour = "black"),legend.key.size=unit(6,"point")) +
  xlab("") + 
  ylab("Counts") +
  scale_fill_manual(values = color_plot) +
  guides(fill = guide_legend(title = "GO Terms")) +
  coord_flip(ylim = c(75,150))


####################################################

#plotting biological process terms

####################################################

replot_bio_terms <- data.frame("Protein Identifier" = character(), "fasta header" = character(), "Regulation" = character(), "Ontology Term Name" = character(), "Ontology ID" = character())
for (i in 1:nrow(all_filtered_up_bio_terms)){
  temp_dataframe <- data.frame(all_filtered_up_bio_terms[i,1], all_filtered_up_bio_terms[i,2], "Upregulated", all_filtered_up_bio_terms[i,4], all_filtered_up_bio_terms[i,5])
  colnames(temp_dataframe) <- colnames(replot_bio_terms)
  replot_bio_terms <- rbind(replot_bio_terms, temp_dataframe)
}
for (i in 1:nrow(all_filtered_down_bio_terms)){
  temp_dataframe <- data.frame(all_filtered_down_bio_terms[i,1], all_filtered_down_bio_terms[i,2], "Downregulated", all_filtered_down_bio_terms[i,4], all_filtered_down_bio_terms[i,5])
  colnames(temp_dataframe) <- colnames(replot_bio_terms)
  replot_bio_terms <- rbind(replot_bio_terms, temp_dataframe)
}


#order_list_replot_bio <- replot_bio_terms %>% group_by(replot_bio_terms$Ontology.Term.Name, replot_bio_terms$Regulation) %>% 
 # dplyr::summarise(n=n()) %>% arrange(desc(n))

manual_ordering_based_on_order_list_and_personal_preferences_bio <- c("ubiquitin-dependent protein catabolic process",
                                                                      "glycolytic process",
                                                                      "carbohydrate metabolic process",
                                                                      "proteolysis",
                                                                      
                                                                      "cellular amino acid biosynthetic process",
                                                                      "amide biosynthetic process",
                                                                      "regulation of DNA replication",
                                                                      "'de novo' IMP biosynthetic process",
                                                                      "carbohydrate catabolic process",
                                                                      "serine-type carboxypeptidase activity",
                                                                      "cysteine-type peptidase activity",
                                                                      "hydrolase activity",
                                                                      "O-methyltransferase activity",
                                                                      "purine nucleotide biosynthetic process",
                                                                      "glucose catabolic process",
                                                                      "isocitrate metabolic process",
                                                                      "glutamate metabolic process",
                                                                      "cellulose biosynthetic process",
                                                                      "aspartic-type endopeptidase activity",
                                                                      "cellular amino acid metabolic process",
                                                                      "translation elongation factor activity",
                                                                      "nucleoside metabolic process",
                                                                      "methionine biosynthetic process",
                                                                      "small molecule metabolic process",
                                                                      "protein metabolic process",
                                                                      "tRNA aminoacylation for protein translation",
                                                                      "catalytic activity",
                                                                      "oxidation-reduction process",
                                                                      "metabolic process",
                                                                      "translation",
                                                                      
                                                                      "response to oxidative stress", 
                                                                      "peroxiredoxin activity",
                                                                      "glutathione peroxidase activity",
                                                                      "peroxidase activity")

replot_bio_terms$Ontology.Term.Name <- factor(replot_bio_terms$Ontology.Term.Name, levels = manual_ordering_based_on_order_list_and_personal_preferences_bio)


replot_bio_group1 <- ggplot(replot_bio_terms, 
                            aes(x = factor(replot_bio_terms$Ontology.Term.Name),
                                y=..count.., fill = replot_bio_terms$Regulation)) +
  geom_bar(colour = "black") + 
  theme_bw() +
  theme(axis.text.x =  element_text(hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), axis.title.x = element_text(color = "black"),
        text = element_text(size=9 ,family = "Arial", colour = "black"),legend.position="none") +
  xlab("") + 
  scale_fill_manual(values = color_plot) +
  ylab("Counts")+
  coord_flip(ylim = c(0,15)) 

replot_bio_group2 <- ggplot(replot_bio_terms, 
                            aes(x = factor(replot_bio_terms$Ontology.Term.Name),
                                y=..count.., fill = replot_bio_terms$Regulation)) +
  geom_bar(colour = "black") + 
  theme_bw() +
  theme(axis.text.y =  element_blank(), axis.text.x = element_text(color = "black"), axis.title.x = element_text(color = "black"),
        text = element_text(size= 9, family = "Arial", colour = "black"),legend.key.size=unit(6,"point")) +
  xlab("") + 
  scale_fill_manual(values = color_plot) +
  ylab("Counts") +
  guides(fill = guide_legend(title = "GO Terms")) +
  coord_flip(ylim = c(20,60))

#######################################################

#save plots to current directory
                           
######################################################

#note may need to reset directory (can use getwd() to determine current directory and setwd() to set a new directory path)

ggsave("MF_low.png", plot = replot_mol_group1, width = 5.8, height = 3.5, dpi = 300, units = "in")
ggsave("MF_high.png", plot = replot_mol_group2, width = 2.5, height = 3.5, dpi = 300, units = "in")


ggsave("BP_low.png", plot = replot_bio_group1, width = 3.5, height = 4, dpi = 300, units = "in")
ggsave("BP_high.png", plot = replot_bio_group2, width = 2.5, height = 4, dpi = 300, units = "in")

#######################################################

#heatmapploting

#######################################################


file_to_open <- file.choose() #choose heatmap Z scores
heatmap_values <- as.data.frame(read_excel(file_to_open,sheet =1, col_names = TRUE))
heatmap_values$`Fasta headers` <- as.factor(heatmap_values$`Fasta headers`)
heatmap_values <- heatmap_values[,c(1,2,4,8,9,3,5,6,7)]


heatmap_matrix <- as.matrix(heatmap_values[,2:9])
rownames_heatmap <- heatmap_values$`Fasta headers`
rownames_heatmap <- as.character(rownames_heatmap)
                                 

for (i in 1:length(rownames_heatmap)){
  rownames_heatmap[[i]] <- substring(rownames_heatmap[[i]],2,regexpr("[0-9]m", rownames_heatmap[[i]])+1)
  print(rownames_heatmap[[i]])
}

rownames_heatmap <- as.factor(rownames_heatmap)
rownames(heatmap_matrix) <- rownames_heatmap

all_go_terms <- data.frame("Upregulated Bio - GO" = character(length(rownames_heatmap)), 
                           "Downregulated Bio - GO" = character(length(rownames_heatmap)),
                           "Upregulated Mol - GO" = character(length(rownames_heatmap)), 
                           "Downregulated Mol - GO" = character(length(rownames_heatmap)), stringsAsFactors = FALSE)
all_go_terms <- cbind(rownames_heatmap, all_go_terms)
all_go_terms$rownames_heatmap <- as.character(all_go_terms$rownames_heatmap)

for (i in 1:nrow(all_go_terms)){
  #if the protein doesn't have a significant go term
  if(all_go_terms[i,1] %in% all_filtered_up_bio_terms$Protein.Phytozome.identifier == FALSE){
    all_go_terms[i,2] <- "NA"
  }
  if(all_go_terms[i,1] %in% all_filtered_down_bio_terms$Protein.Phytozome.identifier == FALSE){
     all_go_terms[i,3] <- "NA"
  }
  if(all_go_terms[i,1] %in% all_filtered_up_mol_terms$Protein.Phytozome.identifier == FALSE){
    all_go_terms[i,4] <- "NA"
  }
  if(all_go_terms[i,1] %in% all_filtered_down_mol_terms$Protein.Phytozome.identifier == FALSE){
    all_go_terms[i,5] <- "NA"
  }
  # if the protein has a go term
  if(all_go_terms[i,1] %in% all_filtered_up_bio_terms$Protein.Phytozome.identifier == TRUE){
    line <-all_filtered_up_bio_terms[all_filtered_up_bio_terms$Protein.Phytozome.identifier %in% all_go_terms[i,1],]
    all_go_terms[i,2] <- levels(line$Ontology.Term.Name)[line$Ontology.Term.Name]
  }
  if(all_go_terms[i,1] %in% all_filtered_down_bio_terms$Protein.Phytozome.identifier == TRUE){
    line <-all_filtered_down_bio_terms[all_filtered_down_bio_terms$Protein.Phytozome.identifier %in% all_go_terms[i,1],]
    all_go_terms[i,3] <- levels(line$Ontology.Term.Name)[line$Ontology.Term.Name]
  }
  if(all_go_terms[i,1] %in% all_filtered_up_mol_terms$Protein.Phytozome.identifier == TRUE){
    line <-all_filtered_up_mol_terms[all_filtered_up_mol_terms$Protein.Phytozome.identifier %in% all_go_terms[i,1],]
    all_go_terms[i,4] <- levels(line$Ontology.Term.Name)[line$Ontology.Term.Name]
  }
  if(all_go_terms[i,1] %in% all_filtered_down_mol_terms$Protein.Phytozome.identifier == TRUE){
    line <-all_filtered_down_mol_terms[all_filtered_down_mol_terms$Protein.Phytozome.identifier %in% all_go_terms[i,1],]
    all_go_terms[i,5] <- levels(line$Ontology.Term.Name)[line$Ontology.Term.Name]
  }
}


######################################################################

# Prep to plot heatmap via complex heatmap package

######################################################################


Infection_catagories <- c("Uninfected","Uninfected","Uninfected","Uninfected","Infected","Infected","Infected","Infected")
infection_dataframe <- data.frame(colnames(heatmap_matrix), Infection_catagories, stringsAsFactors = FALSE)
infection_dataframe$colnames.heatmap_matrix. <- as.character(infection_dataframe$colnames.heatmap_matrix.)


color_go_up_bio <- c(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(all_go_terms$Upregulated.Bio...GO))-1))
color_go_down_mol <- c(colorRampPalette(brewer.pal(8, "Set2"))(length(unique(all_go_terms$Downregulated.Mol...GO))-1))

color_go_up_mol <- c(colorRampPalette(brewer.pal(8, "Set3"))(length(unique(all_go_terms$Upregulated.Mol...GO))-1))
color_go_down_bio <- c(colorRampPalette(brewer.pal(8, "Set1"))(length(unique(all_go_terms$Downregulated.Bio...GO))-1))

go_term_annotation <- ComplexHeatmap::rowAnnotation("Upregulated Bio - GO" = all_go_terms$Upregulated.Bio...GO,
                                                    "Downregulated Bio - GO" = all_go_terms$Downregulated.Bio...GO,
                                                    "Upregulated Mol - GO" = all_go_terms$Upregulated.Mol...GO,
                                                    "Downregulated Mol - GO" = all_go_terms$Downregulated.Mol...GO,
                                                        col = list("Upregulated Bio - GO" =  c("NA" = "white",
                                                                                               "carbohydrate metabolic process" = "#1B9E77",
                                                                                               "proteolysis" = "#D95F02",
                                                                                               "peroxidase activity" = "#7570B3",                  
                                                                                               "ubiquitin-dependent protein catabolic process" = "#E7298A",
                                                                                               "response to oxidative stress" = "#66A61E",           
                                                                                               "glutathione peroxidase activity" = "#E6AB02",
                                                                                               "glycolytic process" = "#A6761D",                   
                                                                                               "peroxiredoxin activity" = "#666666"),
                                                                   "Downregulated Bio - GO" = c("NA" = "white",
                                                                                                "carbohydrate metabolic process" = "#E41A1C",              
                                                                                                "glycolytic process" =  "#BA3241",
                                                                                                "oxidation-reduction process" = "#904A67",             
                                                                                                "nucleoside metabolic process" = "#66628C",
                                                                                                "aspartic-type endopeptidase activity" = "#3C7AB2",     
                                                                                                "isocitrate metabolic process" = "#3B88A1",
                                                                                                "translation" =  "#409386",                      
                                                                                                "metabolic process" = "#469F6C",
                                                                                                "cellular amino acid metabolic process" = "#4BAB51",     
                                                                                                "ubiquitin-dependent protein catabolic process" = "#599E59",
                                                                                                "translation elongation factor activity" = "#6C866E",
                                                                                                "tRNA aminoacylation for protein translation" = "#7E6F84",
                                                                                                "protein metabolic process" = "#905899",        
                                                                                                "catalytic activity" = "#A6548C",
                                                                                                "small molecule metabolic process" = "#BF6065",             
                                                                                                "glucose catabolic process" = "#D76C3D",
                                                                                                "cellulose biosynthetic process" = "#F07816",  
                                                                                                "cellular amino acid biosynthetic process" = "#FF8C05",
                                                                                                "glutamate metabolic process" = "#FFAB11",
                                                                                                "cysteine-type peptidase activity" = "#FFCA1D",            
                                                                                                "methionine biosynthetic process" = "#FFE82A",         
                                                                                                "proteolysis" = "#F8F332",                           
                                                                                                "O-methyltransferase activity" = "#E3CA2F",              
                                                                                                "amide biosynthetic process" = "#CDA12C",    
                                                                                                "hydrolase activity" = "#B8782A",                  
                                                                                                "serine-type carboxypeptidase activity" = "#A8572D",      
                                                                                                "regulation of DNA replication" = "#BC6151",      
                                                                                                "purine nucleotide biosynthetic process" = "#CF6C76",       
                                                                                                "'de novo' IMP biosynthetic process" = "#E3769A",     
                                                                                                "carbohydrate catabolic process" = "#F781BF"),
                                                                   "Upregulated Mol - GO" = c("NA" = "white",
                                                                                              "peroxidase activity" = "#8DD3C7",
                                                                                              "catalytic activity" = "#D5EFBA",
                                                                                              "hydrolase activity, hydrolyzing O-glycosyl compounds" = "#EDECBD",
                                                                                              "hydrolase activity" = "#C3C0D6",
                                                                                              "serine-type carboxypeptidase activity" = "#DF9AA1",
                                                                                              "peptidase activity" = "#E48883",
                                                                                              "carboxypeptidase activity" =  "#96A8C1",                         
                                                                                              "cysteine-type peptidase activity" = "#B8B29F",
                                                                                              "peptidase activity, acting on L-amino acid peptides" = "#F6B762",
                                                                                              "oxidoreductase activity" = "#C7D267",
                                                                                              "hydrolase activity, acting on glycosyl bonds" = "#CDD796",    
                                                                                              "hydro-lyase activity" = "#FCCDE5"),
                                                                   "Downregulated Mol - GO" = c("NA" = "white", 
                                                                                                "hydrolase activity, hydrolyzing O-glycosyl compounds" = "#66C2A5",                                     
                                                                                                "glycolytic process" =  "#98B08E",                                                             
                                                                                                "oxidoreductase activity" =  "#CA9E78",                                                
                                                                                                "catalytic activity" = "#FC8D62",                                                   
                                                                                                "oxidoreductase activity, acting on the aldehyde or oxo group of donors, NAD or NADP as acceptor" = "#D79384",
                                                                                                "aspartic-type endopeptidase activity" = "#B299A7",                                                        
                                                                                                "structural constituent of ribosome" = "#8DA0CB",                                                
                                                                                                "endopeptidase activity" = "#AB98C8",                                                                  
                                                                                                "ligase activity" = "#C991C5",                                                                          
                                                                                                "peptidase activity" =  "#E78AC3",                                                                    
                                                                                                "magnesium ion binding" = "#D1A39E",                                                                
                                                                                                "phosphoglycerate mutase activity" = "#BBBD79",                                                            
                                                                                                "intramolecular transferase activity" = "#A6D854",                                                   
                                                                                                "glutamate decarboxylase activity" = "#C3D847",                                                
                                                                                                "cysteine-type peptidase activity" = "#E1D83B",                                                   
                                                                                                "O-methyltransferase activity" = "#FFD92F",                                                    
                                                                                                "oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor" = "#F6D250",   
                                                                                                "hydrolase activity" = "#EDCB72",                         
                                                                                                "structural molecule activity" = "#E5C494",                                      
                                                                                                "serine-type carboxypeptidase activity" = "#D4BE9E",                                         
                                                                                                "threonine-type endopeptidase activity" = "#C3B8A8",                                 
                                                                                                "serine-type endopeptidase activity" = "#B3B3B3")),
                                                    width = unit(6, "mm"))


z_score_color <-  circlize::colorRamp2(c(-2,0,2), c("#3C7AB2", "white","#BA3241"))
font_arial <- windowsFonts(Times=windowsFont("Arial"))



ComplexHeatmap::ht_opt(legend_border = "black", heatmap_border = TRUE, annotation_border = TRUE)

pdf("Heatmap Z scores_v2.pdf", width = 7, height = 11)

Heatmap_everything <- ComplexHeatmap::Heatmap(heatmap_matrix, 
                       border = TRUE,
                       col = z_score_color,
                       width = unit(4.5, "in"),
                       height = unit(9.5, "in"),


                       #annotate treatments 
                       top_annotation = HeatmapAnnotation(Treatment = anno_block(gp = gpar(fill = 3:2),
                                                                                 labels = c("Uninfected", "Infected"),
                                                                                 labels_gp = gpar(col = "black", fontsize = 14, fontfamily = "Arial")),
                                                          height = unit(0.8, "cm")),
                       #top_annotation = ComplexHeatmap::HeatmapAnnotation(Treatment = infection_dataframe$Infection_catagories, 
                        #                                                  col = list(Treatment =c("Uninfected" = "#599E59", "Infected" = "#FF8C05")),
                         #                                                 gp = gpar(fontsize = 14,fontfamily = "Arial")),

                                                                          
                                              
                       #column modifications
                       column_split = 2, 
                       show_column_names = FALSE,

                       #row modifications + dendrogram
                       show_row_names = FALSE,  
                       row_dend_width = unit(30, "mm"), 
                       #row_dend_gp = 

                       #details regarding modifying legend
                       heatmap_legend_param = list(col_fun = z_score_color, 
                                                   legend_width = unit(30, "mm"),
                                                   legend_height = unit(38, "mm"),
                                                   title = "Z-score", 
                                                   border = "black",
                                                   title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Arial"),
                                                   labels_gp = gpar(fontsize = 12,fontfamily = "Arial"),
                                                   #adjust_annotation_extension = TRUE,
                                                   grid_width = unit(0.5, "cm"),
                                                   legend_label_gp = gpar(col = "black",fontsize = 12, fontfamily = "Arial")))

ComplexHeatmap::draw(Heatmap_everything,    
                     heatmap_legend_side = "right", padding = unit(c(5,5,5,5), "mm"))

dev.off()

#annotation_legend_side = "right", merge_legend = TRUE,
##########################################################

#write to file go-terms supplemental table

#########################################################


#sup table
final_names_up_mol <- all_filtered_up_mol_terms[,c(1,2,4,5)]
final_names_up_bio <- all_filtered_up_bio_terms[,c(1,2,4,5)]
final_names_down_mol <- all_filtered_down_mol_terms[,c(1,2,4,5)]
final_names_down_bio <- all_filtered_down_bio_terms[,c(1,2,4,5)]

write.xlsx(final_names_up_mol, "Supplemental_Table_GO_terms.xlsx", sheetName = "Upregulated_Molecular_Function", col.names = T, row.names = F)
write.xlsx(final_names_up_bio, "Supplemental_Table_GO_terms.xlsx", sheetName = "Upregulated_Biological_Function", col.names = T, row.names = F, append = T)
write.xlsx(final_names_down_mol, "Supplemental_Table_GO_terms.xlsx", sheetName = "Downregulated_Molecular_Function", col.names = T, row.names = F, append = T)
write.xlsx(final_names_down_bio, "Supplemental_Table_GO_terms.xlsx", sheetName = "Downregulated_Biological_Function", col.names = T, row.names = F, append = T)



########################################################

# check glycoltic terms for reviewer during review process

#########################################################

# pull glytoitc process hits in processed file
filehold <- replot_bio_terms[grepl('glycolytic', replot_bio_terms$Ontology.Term.Name),]

# collct all bio terms from origin analysis before filtering + filter for hit for protein tags from above anaylsis
original_file <- rbind(total_names_up_bio, total_names_down_bio)
filehold_original_filtered <- colnames(original_file)
for (i in 1:nrow(filehold)){
  temp_df <- original_file[grepl(as.character(filehold[i,1]), original_file$`Protein Phytozome identifier`),]
  filehold_original_filtered <- rbind(filehold_original_filtered, temp_df)
}

write.xlsx(filehold, file = "Check_glycotylic_processes.xlsx", sheetName = "Filtered", col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(filehold_original_filtered, file = "Check_glycotylic_processes.xlsx", sheetName = "Original_Unfiltered", col.names = TRUE, row.names = TRUE, append = TRUE)

close(file("Check_glycotylic_processes.xlsx"))
