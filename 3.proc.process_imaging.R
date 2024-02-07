source('0.basics.R')


# Process H&E slide images for cell counts --------------------------

#inputs required
input_folder <- "EPICC/FFPE_samples/Imaging_data_raw/"

#main code
cc1 <- read_excel(paste0(input_folder, "CellStats_5_apr_22.xlsx"), sheet="Cell Counts")
names(cc1)[1] <- "sample"
cc2 <- read_excel(paste0(input_folder, "CellStats_26_apr_22_pt1.xlsx"), sheet="Cell Counts")
names(cc2)[1] <- "sample"
cc3 <- read_excel(paste0(input_folder, "CellStats_26_apr_22_pt2.xlsx"), sheet="Cell Counts")
names(cc3)[1] <- "sample"

cc4 <- cc3[grepl("norm", cc3$sample),] #only normal samples
cc <- rbind(cc1, cc2, cc4) #all tumour samples and normal samples from 1st and 2nd batch

# Process node samples
cc_node <- cc[grep("node|EC29_A25_H&E_1_LCM_EC29B25N1", cc$sample),]
cc_node$patient <- substr(cc_node$sample, 1, 4)
cc_node$slide_number <- substr(cc_node$sample, 1, 8)
all_slide_numbers <- unique(cc_node$slide_number)

#apply same rows as above and merge rows with multiple sub-regions per node
cc_node_final <-data.frame()

for (slide_no in all_slide_numbers){
  df1 <- cc_node[grep(slide_no, cc_node$slide_number),]
  if (nrow(df1)==1 | slide_no=="EC29_A22" | slide_no=="EC40_A20"){
    cc_node_final <- rbind(cc_node_final, df1)
  }else{
    for (j in 1:2){
      df2 <- df1[grep(paste0("node_",j), df1$sample),]
      if (nrow(df2)==1){
        cc_node_final <- rbind(cc_node_final, df2)
      }else if (nrow(df2)>1){
        num_sums <- df2 %>%
          summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
        df2$sample_simplified <- gsub("_a|_b|_c|_d|", "", df2$sample)
        
        df_for_bind <- data.frame(unique(df2$sample_simplified), num_sums, unique(df2$patient), unique(df2$slide_number))
        colnames(df_for_bind) <- colnames(cc_node)
        
        cc_node_final <- rbind(cc_node_final, df_for_bind)
      }
    }
  }
}

#manually specify below re merging rows
m1 <- cc_node_final[grep("EC07_A18|EC07_A19", cc_node_final$sample),]
num_sums1 <- m1 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
df_for_bind1 <-data.frame("EC07_A18A19_H&E_1_node_1", num_sums1, "EC07", "EC07_A18A19")
colnames(df_for_bind1) <- colnames(cc_node_final)
cc_node_final <- rbind(cc_node_final, df_for_bind1)
cc_node_final <- cc_node_final[!grepl("^EC07_A18$|^EC07_A19$", cc_node_final$slide_number),]

m2 <- cc_node_final[grep("EC12_A14|EC12_A15", cc_node_final$sample),]
num_sums2 <- m2 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
df_for_bind2 <-data.frame("EC12_A14A15_H&E_1_node_1", num_sums2, "EC12", "EC12_A14A15")
colnames(df_for_bind2) <- colnames(cc_node_final)
cc_node_final <- rbind(cc_node_final, df_for_bind2)
cc_node_final <- cc_node_final[!grepl("^EC12_A14$|^EC12_A15$", cc_node_final$slide_number),]


m3 <- cc_node_final[grep("EC29_A22", cc_node_final$sample),]
num_sums3 <- m3 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
df_for_bind3 <-data.frame("EC29_A22_H&E_1_node_1", num_sums3, "EC29", "EC29_A22")
colnames(df_for_bind3) <- colnames(cc_node_final)
cc_node_final <- rbind(cc_node_final, df_for_bind3)
cc_node_final <- cc_node_final[!grepl("EC29B22N", cc_node_final$sample),]

m4 <- cc_node_final[grep("EC40_A20", cc_node_final$sample),]
num_sums4 <- m4 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
df_for_bind4 <-data.frame("EC40_A20_H&E_1_node_1", num_sums4, "EC40", "EC40_A20")
colnames(df_for_bind4) <- colnames(cc_node_final)
cc_node_final <- rbind(cc_node_final, df_for_bind4)
cc_node_final <- cc_node_final[!grepl("EC40_A20.*_1_a|EC40_A20.*_1_b|EC40_A20.*_1_c", cc_node_final$sample),]

cc_node_final$sample_type <- "node"

# Process superficial, invasive margin and normal samples
cc_rest <- cc[grep("tumour|IM|norm-adj$|normal-adj$|normal-adjacent$|comp|norm-far$|normal-far$", cc$sample, ignore.case = TRUE),]
cc_rest$patient <- substr(cc_rest$sample, 1, 4)
cc_rest$slide_number <- substr(cc_rest$sample, 1, 11)
all_slide_numbers <- unique(cc_rest$slide_number)

n1 <- cc_rest[grep("EC10_block6_HE_4 - 2020-02-28 15.00.09........._EC10B06-norm-adj-comp", cc_rest$sample),]
num_sums5 <- n1 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
df_for_bind5 <-data.frame("EC10_A06_H&E_4_norm-adj", num_sums5, "EC10", "EC10_A06")
colnames(df_for_bind5) <- colnames(cc_rest)
cc_rest <- rbind(cc_rest, df_for_bind5)
cc_rest <- cc_rest[!grepl("EC10_block6_HE_4 - 2020-02-28 15.00.09........._EC10B06-norm-adj-comp", cc_rest$sample),]

n2 <- cc_rest[grep("EC14_A13_H&E_1_.*C14B13-norm-adj-comp", cc_rest$sample),]
num_sums6 <- n2 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
df_for_bind6 <-data.frame("EC14_A13_H&E_1_norm-adj", num_sums6, "EC14", "EC14_A13")
colnames(df_for_bind6) <- colnames(cc_rest)
cc_rest <- rbind(cc_rest, df_for_bind6)
cc_rest <- cc_rest[!grepl("EC14_A13_H&E_1_.*C14B13-norm-adj-comp", cc_rest$sample),]

n3 <- cc_rest[grep("EC07_A09_H&E_1_EC07B09-normal-adj-comp", cc_rest$sample),]
num_sums7 <- n3 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
df_for_bind7 <-data.frame("EC07_A09_H&E_1_norm-adj", num_sums7, "EC07", "EC07_A09")
colnames(df_for_bind7) <- colnames(cc_rest)
cc_rest <- rbind(cc_rest, df_for_bind7)
cc_rest <- cc_rest[!grepl("EC07_A09_H&E_1_EC07B09-normal-adj-comp", cc_rest$sample),]

cc_rest$sample_type <- NA

for (i in 1:nrow(cc_rest)){
  if (grepl("tumour", cc_rest$sample[i], ignore.case=TRUE)){
    cc_rest$sample_type[i] <- "tumour"
  }
  if (grepl("IM", cc_rest$sample[i], ignore.case=TRUE)){
    cc_rest$sample_type[i] <- "IM"
  }
  if (grepl("norm-adj$|normal-adj$|normal-adjacent$", cc_rest$sample[i], ignore.case=TRUE)){
    cc_rest$sample_type[i] <- "adjacent_normal_mucosa"
  }
  if (grepl("norm-far$|normal-far$", cc_rest$sample[i], ignore.case=TRUE)){
    cc_rest$sample_type[i] <- "normal_mucosa"
  }
}

all_cell_counts_epicc <- as.data.frame(rbind(cc_rest,cc_node_final))

#sum for total epithelial/ stromal
all_cell_counts_epicc$total_epithelial <- all_cell_counts_epicc$`Cancer Epithelial Cell` + all_cell_counts_epicc$`Normal Epithelial Cell` #total epi calculated. normal likely misclassification
all_cell_counts_epicc$total_stromal <- all_cell_counts_epicc$Fibroblast + all_cell_counts_epicc$Myocyte

#calculate number of immune cells per cancer cell
all_cell_counts_epicc$lym_to_all_epi <- all_cell_counts_epicc$Lymphocyte / all_cell_counts_epicc$total_epithelial 
all_cell_counts_epicc$neu_to_all_epi <- all_cell_counts_epicc$Neutrophil / all_cell_counts_epicc$total_epithelial
all_cell_counts_epicc$mac_to_all_epi <- all_cell_counts_epicc$Macrophage / all_cell_counts_epicc$total_epithelial 
all_cell_counts_epicc$endo_to_all_epi <- all_cell_counts_epicc$`Endothelial Cell` / all_cell_counts_epicc$total_epithelial 
all_cell_counts_epicc$stromal_to_all_epi <- all_cell_counts_epicc$total_stromal / all_cell_counts_epicc$total_epithelial 

#calculate fraction of immune cells
all_cell_counts_epicc$total_cells <- rowSums(all_cell_counts_epicc[2:9])
all_cell_counts_epicc$lym_fraction<- all_cell_counts_epicc$Lymphocyte/ all_cell_counts_epicc$total_cells
all_cell_counts_epicc$neu_fraction <- all_cell_counts_epicc$Neutrophil / all_cell_counts_epicc$total_cells
all_cell_counts_epicc$mac_fraction <- all_cell_counts_epicc$Macrophage / all_cell_counts_epicc$total_cells 
all_cell_counts_epicc$endo_fraction <- all_cell_counts_epicc$`Endothelial Cell` / all_cell_counts_epicc$total_cells 
all_cell_counts_epicc$stromal_fraction <- all_cell_counts_epicc$total_stromal / all_cell_counts_epicc$total_cells 
all_cell_counts_epicc$epi_fraction <- all_cell_counts_epicc$total_epithelial / all_cell_counts_epicc$total_cells
