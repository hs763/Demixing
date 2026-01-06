library(Matrix)
matrix_dir = "/cephfs2/hannas/demixing/results/human_chimp_macaq/human_HCM_ref/outs/raw_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
human_mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(human_mat) = barcode.names$V1
rownames(human_mat) = feature.names$V1

matrix_dir = "/cephfs2/hannas/demixing/results/human_chimp_macaq/chimp_HCM_ref/outs/raw_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
chimp_mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(chimp_mat) = barcode.names$V1
rownames(chimp_mat) = feature.names$V1

matrix_dir = "/cephfs2/hannas/demixing/results/human_chimp_macaq/macaq_HCM_ref/outs/raw_feature_bc_matrix/"
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
macaq_mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(macaq_mat) = barcode.names$V1
rownames(macaq_mat) = feature.names$V1

##############################################################################################################################################
#preparing metadata
#examining dimentions
dim(human_mat)
[1]  45085 454945
> dim(chimp_mat)
[1]  45085 495782
> dim(macaq_mat)
[1]  45085 284410

#finding the rows with human genes and chimp genes
human_mat_row_h <- grep("^GRCh38", rownames(human_mat))
length(human_mat_row_h) #20070
human_mat_row_ch <- grep("^PanTro3", rownames(human_mat))
length(human_mat_row_ch) #23302
human_mat_row_m <- grep("^Mnem", rownames(human_mat))
length(human_mat_row_m) #1713

chimp_mat_row_h <- grep("^GRCh38", rownames(chimp_mat))
length(chimp_row_h) #20070
chimp_row_ch <- grep("^PanTro3", rownames(chimp_mat))
length(chimp_row_ch) #23302
chimp_row_m <- grep("^Mnem", rownames(chimp_mat))
length(chimp_row_m) #1713

chimp_mat_row_h <- grep("^GRCh38", rownames(macaq_mat))
length(chimp_row_h) #20070
chimp_row_ch <- grep("^PanTro3", rownames(macaq_mat))
length(chimp_row_ch) #23302
chimp_row_m <- grep("^Mnem", rownames(macaq_mat))
length(chimp_row_m) #1713

human_mat_h <- human_mat[1:20070,]
dim(human_mat_h) #[1]  20070 454945
h_hcm_counts_h <- colSums(human_mat_h)
sum(h_hcm_counts_h) #22851169

human_mat_ch <- human_mat[20071:43372,]
dim(human_mat_ch) #23302 454945
h_hcm_counts_ch <- colSums(human_mat_ch)
sum(h_hcm_counts_ch) #560619

human_mat_m <- human_mat[43373:45085,]
dim(human_mat_m) #1713 454945
h_hcm_counts_m <- colSums(human_mat_m)
sum(h_hcm_counts_m) #11553


chimp_mat_h <- chimp_mat[1:20070,]
dim(chimp_mat_h) #[1]  20070 495782
ch_hcm_counts_h <- colSums(chimp_mat_h)
sum(ch_hcm_counts_h) #6474102

chimp_mat_ch <- chimp_mat[20071:43372,]
dim(chimp_mat_ch) #23302 495782
ch_hcm_counts_ch <- colSums(chimp_mat_ch)
sum(ch_hcm_counts_ch) #21052271

chimp_mat_m <- chimp_mat[43373:45085,]
dim(chimp_mat_m) #1713 495782
ch_hcm_counts_m <- colSums(chimp_mat_m)
sum(ch_hcm_counts_m) #30124


macaq_mat_h <- macaq_mat[1:20070,]
dim(macaq_mat_h) #[1]  20070 284410
m_hcm_counts_h <- colSums(macaq_mat_h)
sum(m_hcm_counts_h) #1718257

macaq_mat_ch <- macaq_mat[20071:43372,]
dim(macaq_mat_ch) #23302 284410
m_hcm_counts_ch <- colSums(macaq_mat_ch)
sum(m_hcm_counts_ch) #1196237

macaq_mat_m <- macaq_mat[43373:45085,]
dim(macaq_mat_m) #1713 284410
m_hcm_counts_m <- colSums(macaq_mat_m)
sum(m_hcm_counts_m) #846869

metadata_h <- cbind(h_hcm_counts_h, h_hcm_counts_ch, h_hcm_counts_m)
metadata_ch <- cbind(ch_hcm_counts_h, ch_hcm_counts_ch, ch_hcm_counts_m)
metadata_m <- cbind(m_hcm_counts_h, m_hcm_counts_ch, m_hcm_counts_m)

empty_cells_h <- grep(0, metadata_h[,1])
length(empty_cells_h) #191106
metadata2_h <- metadata_h
metadata2_h <- metadata2_h[-c(empty_cells_h),]
write.csv(metadata2_h, "metadata_human_HCM_ref_QC.csv")
write.csv(metadata_h, "metadata_human_HCM_ref_full.csv")
dim(metadata_h)   #454945      3
dim(metadata2_h)  #263839      3

empty_cells_ch <- grep(0, metadata_ch[,2])
length(empty_cells_ch) #254013
metadata2_ch <- metadata_ch
metadata2_ch <- metadata2_ch[-c(empty_cells_ch),]
write.csv(metadata2_ch, "metadata_chimp_HCM_ref_QC.csv")
write.csv(metadata_ch, "metadata_chimp_HCM_ref_full.csv")
dim(metadata_ch)  #495782      3
dim(metadata2_ch) #241769      3

empty_cells_m <- grep(0, metadata_m[,3])
length(empty_cells_m) #217196
metadata2_m <- metadata_m
metadata2_m <- metadata2_m[-c(empty_cells_m),]
write.csv(metadata2_m, "metadata_macaq_HCM_ref_QC.csv")
write.csv(metadata_m, "metadata_macaq_HCM_ref_full.csv")
dim(metadata_m)  #284410    3
dim(metadata2_m) #67214     3

original <- rep("h", nrow(metadata2_h))  
metadata2_h <- cbind(metadata2_h, as.factor(original))
original <- rep("ch", nrow(metadata2_ch))  
metadata2_ch <- cbind(metadata2_ch, as.factor(original))
original <- rep("m", nrow(metadata2_m))  
metadata2_m <- cbind(metadata2_m, as.factor(original))

colnames(metadata2_h) <- c("X","h_counts","ch_counts","m_counts","original_species")
colnames(metadata2_ch) <- c("X","h_counts","ch_counts","m_counts","original_species")
colnames(metadata2_m) <- c("X","h_counts","ch_counts","m_counts","original_species")
combined_matrix <- rbind(metadata2_h, metadata2_ch, metadata2_m) #dim: 1044031       5
combined_matrix$total_counts <- combined_matrix$h_counts + combined_matrix$ch_counts + combined_matrix$m_counts
combined_matrix$h_counts_ratio <- combined_matrix$h_counts/combined_matrix$total_counts
combined_matrix$ch_counts_ratio <- combined_matrix$ch_counts/combined_matrix$total_counts
combined_matrix$m_counts_ratio <- combined_matrix$m_counts/combined_matrix$total_counts

##############################################################################################################################################
#predicting species:
max_ratio <- apply(combined_matrix[, c("h_counts_ratio", "ch_counts_ratio", "m_counts_ratio")], 1, max)
combined_matrix$predicted_species <- ifelse(max_ratio == combined_matrix$h_counts_ratio, "h",
                                            ifelse(max_ratio == combined_matrix$ch_counts_ratio, "ch", "m"))

write.csv(combined_matrix, "combined_matrix.csv")

mispredictions <- which(combined_matrix$predicted_species != combined_matrix$original_species)
length(mispredictions) #168495
dim(combined_matrix) #1044031 10
#% of cells that are mispredicted
168495/1044031 = 0.1613889

##############################################################################################################################################
#plotting all the histograms an the species predictorions 

library(ggplot2) 
#human mapping rates in all cells 
pdf("MappingRateDistribution_Human_CombinedMatrix.pdf")
ggplot(combined_matrix, aes(x = h_counts_ratio, fill = original_species)) + 
  geom_histogram(color = "#B19CD9", bins = 100) + 
  labs(x = "Human Transcript Mapping Rate", y = "Number of Cells", title = "Distribution of Mapping Rate to Human Genome", fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()

#chimp mapping rates in all cells 
pdf("MappingRateDistribution_Chimp_CombinedMatrix.pdf")
ggplot(combined_matrix, aes(x = ch_counts_ratio, fill = original_species)) + 
  geom_histogram(color="#FFB347", bins = 100) + 
  labs(x = "Chimpanzee Transcript Mapping Rate", y = "Number of Cells", title = "Distribution of Mapping Rate to Chimpanzee Genome", fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()

#ggplot(combined_matrix, aes(x=m_counts_ratio)) + geom_histogram(color="#87CEEB" , fill= "#87CEEB" , bins = 100) +labs(x ="Macaque transcript mapping rate", y= "Number of cells",title = "Mapping rate distribution for Macaque transcripts")
#ggplot(combined_matrix, aes(x=ch_counts_ratio)) + geom_histogram(color="#FFB347", fill= "#FFB347", bins = 100) +labs(x ="Chimpanzee transcript mapping rate", y= "Number of cells",title = "Mapping rate distribution for Chimpanzee transcripts")

pdf("MappingRateDistribution_Macaq_CombinedMatrix.pdf")
ggplot(combined_matrix, aes(x = m_counts_ratio, fill = original_species)) + 
  geom_histogram(color="#87CEEB", bins = 100) + 
  labs(x = "Macaque Transcript Mapping Rate", y = "Number of Cells", title = "Distribution of Mapping Rate to Macaque Genome", fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()

pdf("total_counts_distribution_combined_matrix.pdf")
ggplot(combined_matrix, aes(x = log(total_counts))) + 
  geom_density(color = "#B19CD9", fill = "#B19CD9") +  # Corrected color code for line
  labs(x = "Log(Total counts)", y = "Density", title = "Distribution of Log(Total Counts)")
dev.off()

##############################################################################################################################################
#applying a filter to remove the low quality cells 

dim(combined_matrix) # 1044031      11
dim(combined_matrix[combined_matrix$total_counts > 200,]) #10912    11
dim(combined_matrix[combined_matrix$total_counts > 400,]) #9472   11
summary(combined_matrix$total_counts)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
      0       0       1      47       2   41712 

cm <- combined_matrix
cm <_ cm[cm$total_counts >0,]
cm2 <- cm[cm$total_counts > 400,]
pdf("total_counts_distribution_combined_matrix_total_count>400 .pdf")
ggplot(cm2, aes(x = log(total_counts))) + 
  geom_density(color = "#B19CD9", fill = "#B19CD9") +  
  labs(x = "Log(Total counts)", y = "Density", title = "Distribution of Log(Total Counts)")
dev.off()

pdf("histogram_total_counts.pdf")
hp <- ggplot() +
  geom_histogram(data = cm, aes(x = log(total_counts)), color = "#B19CD9", fill = "#B19CD9", alpha = 0.5, bins=50) +
  geom_histogram(data = cm2, aes(x = log(total_counts)), color = "red", fill = "red", alpha = 0.5, bins=50) +  
  labs(x = "Log(Total counts)", y = "Density", title = "Distribution of Log(Total Counts)") +
  theme_minimal()
dev.off()

pdf("zoomed_in_histogram_total_counts.pdf")
zhp <- hp +
  xlim(5, NA) +  # Zoom in to x values between 4 and 9
  ylim(0, NA) +  # Keep the y-axis range unchanged
  labs(title = "Zoomed-in Section")
dev.off()

pdf("histogram_total_counts+zoom-in-section_500.pdf")
hp + annotation_custom(ggplotGrob(zhp), xmin = 5, xmax = Inf, ymin = 100000, ymax = 300000)
dev.off()

pdf("MappingRateDistribution_Human_CombinedMatrix_filter.pdf")
ggplot(cm2, aes(x = h_counts_ratio, fill = original_species)) + 
  geom_histogram(color = "#B19CD9", bins = 100) + 
  labs(x = "Human Transcript Mapping Rate", y = "Number of Cells", title = "Distribution of Mapping Rate to Human Genome", fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()

#chimp mapping rates in all cells 
pdf("MappingRateDistribution_Chimp_CombinedMatrix_filter.pdf")
ggplot(cm2, aes(x = ch_counts_ratio, fill = original_species)) + 
  geom_histogram(color="#FFB347", bins = 100) + 
  labs(x = "Chimpanzee Transcript Mapping Rate", y = "Number of Cells", title = "Distribution of Mapping Rate to Chimpanzee Genome", fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()

pdf("MappingRateDistribution_Macaq_CombinedMatrix_filter.pdf")
ggplot(cm2, aes(x = m_counts_ratio, fill = original_species)) + 
  geom_histogram(color="#87CEEB", bins = 100) + 
  labs(x = "Macaque Transcript Mapping Rate", y = "Number of Cells", title = "Distribution of Mapping Rate to Macaque Genome", fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()


##############################################################################################################################################
#mismatch befire filter 
mispredictions <- which(combined_matrix$predicted_species != combined_matrix$original_species)
length(mispredictions) #168495
dim(combined_matrix) #1044031 10
#% of cells that are mispredicted
168495/1044031 = 0.1613889

#mismatch after filter 
mispredictions <- which(cm2$predicted_species != cm2$original_species)
length(mispredictions) #1622
dim(cm2) #9472
#% of cells that are mispredicted
1622/9472 = 0.1712416

m <- cm2[cm2_m$m_counts_ratio > 0.1,]

mispredictions <- which(cm2_h$predicted_species != cm2_h$original_species) #1617/5712 = 0.2830882
mispredictions <- which(cm2_ch$predicted_species != cm2_ch$original_species) #5/3760 = 0.001329787
mispredictions <- which(cm2_m$predicted_species != cm2_m$original_species)

cm2$prediction_based_on_m <- ifelse(cm2$m_counts_ratio > 0.1, "m", "else")
cm2$prediction_based_on_ch <- ifelse(cm2$ch_counts_ratio > 0.5, "ch", "else")
cm2$prediction_based_on_h <- ifelse(cm2$h_counts_ratio > 0.8, "h", "else")

sum(rowSums(cm2[, c("prediction_based_on_m", "prediction_based_on_h", "prediction_based_on_ch")] == "else") < 2) 
#this is equla to 0 so there werenet any rows that had two species associated to them. 

sum(rowSums(cm2[, c("prediction_based_on_m", "prediction_based_on_h", "prediction_based_on_ch")] == "else") == 3)
#this is equla to 0 so there were no rows that had no species assigned to them. 

#final number of cells (after the filter) was about 10 000 out of 100 000. 
write.csv(cm2, "combined_matrix_after_filter_400.csv")

##############################################################################################################################################
cm2 <- cm[cm$total_counts >500,] #dim: 9081   11
#500 filter justification 
#<https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html>
#"The UMI counts per cell should generally be above 500, that is the low end of what we expect. 
#If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply."

mispredictions <- which(cm2$predicted_species != cm2$original_species)
length(mispredictions) #1553
dim(cm2) #9081
1553/9081=0.1710164

#the rest of the metrics are the same. 

dim(cm2[cm2$h_counts_ratio > 0.8,]) #3911
dim(cm2[cm2$ch_counts_ratio > 0.5,]) #3617 
dim(cm2[cm2$m_counts_ratio > 0.1,]) #1553

pdf("histogram_counts_filter_500.pdf")
ggplot() +
  geom_histogram(data = cm2, aes(x = log(total_counts)), color = "red", fill = "red", alpha = 0.5, bins=50) +  
  labs(x = "Log(Total counts)", y = "Density", title = "Distribution of Log(Total Counts)") +
  theme_minimal()
dev.off()




##############################################################################################################################################
cm <- read.csv("combined_matrix.csv")
cm2 <- read.csv("combined_matrix_after_filter_400.csv")
cm2 <- cm[cm$total_counts >500,]
write.csv(cm2, "combined_matrix_after_filter_500.csv")


pdf("a_MappingRateDistribution_Human_CombinedMatrix_no_color.pdf")
ggplot(cm, aes(x = h_counts_ratio)) + 
  geom_histogram(color = "#B19CD9", fill="#B19CD9", bins = 70) +
  labs(x = "Human Transcript Mapping Rate", y = "Number of Cells") +
  theme(aspect.ratio = 1)
dev.off()

#chimp mapping rates in all cells 
pdf("a_MappingRateDistribution_Chimp_CombinedMatrix_no_color.pdf")
ggplot(cm, aes(x = ch_counts_ratio)) + 
  geom_histogram(color="#FFB347", fill="#FFB347", bins = 70) + 
  labs(x = "Chimpanzee Transcript Mapping Rate", y = "Number of Cells") +
 theme(aspect.ratio = 1)
dev.off()

pdf("a_MappingRateDistribution_Macaq_CombinedMatrix_no_color.pdf")
ggplot(cm, aes(x = m_counts_ratio)) + 
  geom_histogram(color="#87CEEB",fill="#87CEEB", bins = 70) + 
  labs(x = "Macaque Transcript Mapping Rate", y = "Number of Cells") + 
theme(aspect.ratio = 1)
dev.off()




pdf("a_MappingRateDistribution_Human_CombinedMatrix_no_color_filter.pdf")
ggplot(cm2, aes(x = h_counts_ratio)) + 
  geom_histogram(color = "#B19CD9", fill="#B19CD9", bins = 70) +
  labs(x = "Human Transcript Mapping Rate", y = "Number of Cells") +
  theme(aspect.ratio = 1)
dev.off()

#chimp mapping rates in all cells 
pdf("a_MappingRateDistribution_Chimp_CombinedMatrix_no_color_filter.pdf")
ggplot(cm2, aes(x = ch_counts_ratio)) + 
  geom_histogram(color="#FFB347", fill="#FFB347", bins = 70) + 
  labs(x = "Chimpanzee Transcript Mapping Rate", y = "Number of Cells") +
 theme(aspect.ratio = 1)
dev.off()

pdf("a_MappingRateDistribution_Macaq_CombinedMatrix_no_color_filter.pdf")
ggplot(cm2, aes(x = m_counts_ratio)) + 
  geom_histogram(color="#87CEEB",fill="#87CEEB", bins = 70) + 
  labs(x = "Macaque Transcript Mapping Rate", y = "Number of Cells") + 
theme(aspect.ratio = 1)
dev.off()



pdf("a_MappingRateDistribution_Human_CombinedMatrix_colour_filter.pdf")
ggplot(cm2, aes(x = h_counts_ratio, color=original_species, fill = original_species)) + 
  geom_histogram(bins = 70) + 
  labs(x = "Human Transcript Mapping Rate", y = "Number of Cells", fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()

#chimp mapping rates in all cells 
pdf("a_MappingRateDistribution_Chimp_CombinedMatrix_colour_filter.pdf")
ggplot(cm2, aes(x = ch_counts_ratio, color=original_species, fill = original_species)) + 
  geom_histogram(bins = 70) + 
  labs(x = "Chimpanzee Transcript Mapping Rate", y = "Number of Cells", fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()

pdf("a_MappingRateDistribution_Macaq_CombinedMatrix_colour_filter.pdf")
ggplot(cm2, aes(x = m_counts_ratio, color=original_species, fill = original_species)) + 
  geom_histogram( bins = 70) + 
  labs(x = "Macaque Transcript Mapping Rate", y = "Number of Cells",fill="Original Species") +
  theme(aspect.ratio = 1)
dev.off()
