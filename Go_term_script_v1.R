#-----------------------------------------------------------------------------------------------
# Coaker Lab - Plant Pathology Department UC Davis
# Author: Danielle M. Stevens
# Last Updated: 10/1/19
# Script Purpose: Plotting Signifcant Go-terms from proteomics data to assess redundancy; plot final heatmap
# Inputs Necessary: Excel file with signficant go terms, protein id, child and parent go-terms, z-scores
# Outputs: Visualization of comparsion between go-terms called
#-----------------------------------------------------------------------------------------------

#library packages need to load
library(RColorBrewer)
library(readxl)
library(tidyversee)
library(plyr)
library(magrittr)
library(dplyr)
library(forcats)
library(stringr)
library(gplots)
library(remotes)
library(reshape2)
library(gridExtra)
library(circlize)
library(rafalib)
library(ggplot2)

#choose go-terms file to process
file_to_open <- file.choose()
upregulated_bio <- as.data.frame(read_excel(file_to_open, sheet=2, col_names = TRUE))
downregulated_bio <- as.data.frame(read_excel(file_to_open, sheet=3, col_names = TRUE))
upregulated_mol <- as.data.frame(read_excel(file_to_open, sheet=1, col_names = TRUE))
downregulated_mol <- as.data.frame(read_excel(file_to_open, sheet=4, col_names = TRUE))

#ignore for now - old ploting method for dendrogram
#hv <- as.matrix(heatmap_values[,2:9])
#d <- dist(scale(t(hv)))
#h.clust.d <- hclust(d)
#hv <- as.data.frame(hv)
#hv <- hv[,c(h.clust.d$order)]

#d2 <- dist(scale(hv))
#h.clust.d2 <- hclust(d2, method = "ward.D2")
#hv <- hv[c(h.clust.d2$order),]

#rownames_in_order <- heatmap_values$`Fasta headers`
#rownames_in_order <- rownames_in_order[h.clust.d2$order]
#hv <- cbind(rownames_in_order,hv)3testing <- melt(hv, id = c("rownames_in_order"))

#function to convert \t into underscores
remove_tabs <- function(n){
  for (i in 1:nrow(n)){
    n[i,4] <- gsub("\t","_",n[i,4])
    print(n[i,4])
  }
  return(n)
}

for (i in 1:nrow(upregulated_bio)){
  upregulated_bio[i,4] <- gsub("\t","_", upregulated_bio[i,4])
  print(upregulated_bio[i,4])
}

for (i in 1:nrow(downregulated_bio)){
  downregulated_bio[i,4] <- gsub("\t","_", downregulated_bio[i,4])
  print(downregulated_bio[i,4])
}
for (i in 1:nrow(upregulated_mol)){
  upregulated_mol[i,4] <- gsub("\t","_", upregulated_mol[i,4])
  print(upregulated_mol[i,4])
}
for (i in 1:nrow(downregulated_mol)){
  downregulated_mol[i,4] <- gsub("\t","_", downregulated_mol[i,4])
}

#broken function....
remove_tabs(upregulated_bio)
remove_tabs(downregulated_bio)
remove_tabs(upregulated_mol)
remove_tabs(downregulated_mol)

# process and rewrite new file (excel) with go terms by number of hits by phytozome identifier
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

#plot total counts for each significant protein detected
#################### upregulated bio function go terms
fasta_names_upregulated_bio <- upregulated_bio[order(fct_infreq(upregulated_bio$`Protein Phytozome identifier`)),c(1,4)]
fasta_names_upregulated_bio <- distinct(fasta_names_upregulated_bio)

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

######################--------------------filtered go terms











#############################--------------heatmapploting

file_to_open <- file.choose()
p_values <- as.data.frame(read_excel(file_to_open,sheet =1, col_names = TRUE))
#go_terms_citrus <- as.data.frame(read.csv(file_to_open, sep=",", header = TRUE))

file_to_open <- file.choose()
heatmap_values <- as.data.frame(read_excel(file_to_open,sheet =1, col_names = TRUE))
heatmap_values$`Fasta headers` <- as.factor(heatmap_values$`Fasta headers`)
heatmap_values <- heatmap_values[,c(1,2,4,8,9,3,5,6,7)]





heatmap_matrix <- as.matrix(heatmap_values[,2:9])
rownames_heatmap <- heatmap_values$`Fasta headers`
rownames(heatmap_matrix) <- rownames_heatmap


all_bio_go_terms <- rbind(upregulated_bio, downregulated_bio)
test_color <-  as.factor(all_bio_go_terms[1:1355,8])

color_bio_go <- colorRampPalette(brewer.pal(8, "Set2"))(length(unique(test_color)))
names(color_bio_go) <- unique(test_color)

heatmap(heatmap_matrix, col = greenred(300), labRow = FALSE)


col_fun = colorRamp2(c(-3, 0, 3), c("green", "black", "red"))
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5)
circos.initialize(rownames(heatmap_matrix), xlim = cbind(c(0,0), table(rownames(heatmap_matrix))))


#count_child_go_terms <- as.data.frame(table(upregulated_bio$`GO Term`))
#count_child_go_terms <- count_child_go_terms[c(order(count_child_go_terms$Freq,decreasing = T )),]

#count_parent_go_terms <- as.data.frame(table(upregulated_bio$`Parents GO ID`))
#count_parent_go_terms <- count_parent_go_terms[c(order(count_parent_go_terms$Freq, decreasing = T)),]





#alpha = upregulated_bio$`Ontology Term Name`
#  scale_fill_brewer(palette = "Spectral")
#scale_alpha_discrete()+
#  scale_y_continuous(limits = c(0, 10)) + 


ggplot(downregulated_bio, aes(x=downregulated_bio$`Fasta header`, y = ..count.., fill= downregulated_bio$`Parent GO Term Name`)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x =  element_text(angle = 45, hjust = 1)) 

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Ti/Ri plasmid colors
titypes<-c("TypeIa","TypeIb","TypeI", "TypeII", "TypeIII", "TypeIV", "TypeV", "RiTypeI","RiTypeII","RiTypeIII","TypeVI", "At", "chromid", "TypeIVa","TypeIVb","TypeIVc")
plasmidcolors<-c("#E41A1C", "#8B0000", brewer.pal(9,"Set1"), "darkturquoise", "seagreen1", "#6a3672", "#984EA3", "#c194c7")
names(plasmidcolors)<-titypes
#plasmidcolors
#barplot(rep(10,11), col=plasmidcolors)

#Species colors
specieslist<-c("G1","G4","G7","G7_Other","G8","G9","larrymoorei","rhizogenes","rubi","skierniewicense","vitis")
speciescolors<-gg_color_hue(11)
names(speciescolors)<-specieslist
#speciescolors
#barplot(rep(10,11), col=speciescolors)

#Biovar colors
biovarlist<-c("BiovarI","BiovarII","Brucella","BiovarIII","Rhizobium","Rhizobium_Other","Neorhizobium","Sinorhizobium-Ensifer","undicola","Mesorhizobium", "arsenijevicii", "G2", "blue", "black", "Ochrobactrum", "Shinella", "Aureimonas", "Martelella","rhizogenes-like")
biovarcolors<-c("#E41A1C","#4DAF4A","darkgreen","#377EB8","#04048C","cadetblue3","darkcyan","#984EA3","orchid","gray", "salmon", "purple", "blue", "black", "darkolivegreen1", "orange", "lightgoldenrod1", "pink", "dodgerblue3")
#biovarcolors<-c("#E41A1C","#4DAF4A","darkgreen","#377EB8","dodgerblue1","cadetblue3","darkcyan","#984EA3","orchid","gray", "salmon", "purple", "blue", "black", "green", "orange", "yellow", "pink")
names(biovarcolors)<-biovarlist
#biovarcolors

#Combined biovar/species colors
taxonomycolors<-c(speciescolors,biovarcolors)
#taxonomycolors
#barplot(rep(10,21), col=taxonomycolors)

