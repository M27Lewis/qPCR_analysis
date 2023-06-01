library(tidyverse)
library(readxl)
library(writexl)
library(RColorBrewer)
library(broom)
library(DescTools)

### Custom ggplot Theme ###

theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5, margin = margin(b = 5)),
            text = element_text(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(0, 0, 0, 0),
            legend.title= element_blank(), 
            legend.background = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

############################

### Define pathway to the qPCR data
data <- "data/raw/2023-03-28_173637_First_qPCR_for_72h_transfection_PODXL_KD.xls" #path to your qPCR data file

### Name the experiment
experiment <- "PODXL_CasRx_72hr_Harry" #make a name for the experiment to create an output directory

### Create directories to save output files
plots_path <- file.path(paste0("plots/", Sys.Date(), "_", experiment)) #creates a variable to save the path to plots

if(!dir.exists(plots_path)){ #creates the new directory, if it doesn't already exist
  dir.create(plots_path)
}

tables_path <- file.path(paste0("tables/", Sys.Date(), "_", experiment)) #creates a variable to save the path to tables

if(!dir.exists(tables_path)){ #creates the new directory, if it doesn't already exist
  dir.create(tables_path)
}

### Assign name of control sample
ctrl <- "NT" #put the exact name of your control sample here. For example, NT or WT

### Assign names of housekeeping normalization gene(s)
ref <- c("GAPDH", "B-Actin") #put the exact names of the housekeeping normalization genes here so the script can assign them correctly

### Clean up data
qpcr_data <- read_excel(data, sheet = "Results", skip = 44, col_names = TRUE) |> #imports the data excel file and skips the first 44 rows, which don't have data
  slice(1:(n()-5)) |> #removes the last 5 rows, which don't have data
  dplyr::filter(!`Sample Name` %in% c("NTC", "No Template", "NT", "Blank")) |> #removes no template control data
  dplyr::filter(!`Sample Name` %in% c("-RT", "RT-", "NON RT")) |>   #removes the minus reverse transcriptase data
  separate(`Well Position`, into = c("row", "column"), sep = 1, convert = TRUE) |> #separates well position into row and column for plotting
  rename(primer = `Target Name`) #renames Target Name column to be primer

### Visualize plate layout
ggplot(qpcr_data, aes(x = column, y = row, fill = primer, label = `Sample Name`)) + #plots the layout of your qPCR plate. Confirm that it matches what you expect
  geom_tile(colour = "black") +
  geom_text() +
  scale_y_discrete(limits = c("P", "O", "N", "M", "L", "K", "J", "I", "H", "G", "F", "E", "D", "C", "B", "A")) +
  scale_x_continuous(breaks = 1:24)

### Summarize and clean up data further
summ_data <- qpcr_data |> 
  group_by(`Sample Name`, primer) |> 
  summarise(mean_Ct = mean(CT, na.rm = TRUE)) #summarizes the data to get mean Ct

summ_data <- summ_data |>
  separate(`Sample Name`, into = c("sample", "rep"), sep = " ", convert = TRUE) #splits the sample name so it is sample and replicate number. The sep = may need to be adjust if you don't use a space

ggplot(summ_data, aes(x = sample, y = mean_Ct, color = primer)) + #plots summarized Ct values for visual confirmation
  geom_point()

### Prepare reference data for analysis (this section prepares the housekeeping normalization data for use)
if (length(ref) > 1) {
  ref1_data <- summ_data |> 
    filter(primer == ref[1]) |> 
    rename("ref1_Ct" = "mean_Ct")
  
  ref2_data <- summ_data |> 
    filter(primer == ref[2]) |> 
    rename("ref2_Ct" = "mean_Ct")
  
  ref_data <- left_join(ref1_data, ref2_data, by = c("sample", "rep")) |> 
    select(!ends_with(".x") | ends_with(".y"))
  
  ref_data$ref_Ct <- rowMeans(subset(ref_data, select = c(ref1_Ct, ref2_Ct)), na.rm = TRUE)
} else {
  ref_data <- summ_data |> 
    filter(primer == ref) |> 
    rename("ref_Ct" = "mean_Ct") |> 
    select(!primer)
}

### Loop through data and analyze to create plots and perform statistical analysis of results
# this section splits the qPCR data by target gene, then loops through the list subset lists and analyzes the qPCR data creating a bar plot, table of results, and multiple tables of significance values based on the test used. The tests used are pairwise t-test, Dunnett's test, and one-way ANOVA followed by a Tukey test. I prefer the Dunnett's test for qPCR, so it is used in this loop to assign the significance indicators (asterisks) that are placed automatically on the plot. 
test_data <- summ_data |> 
  filter(!primer %in% ref) #remove housekeeping genes from the data

test_list <- split(test_data, test_data$primer) #split the data by target gene/primer

for (i in 1:length(test_list)){
  df_i <- test_list[[i]]
  combined_data <-  merge(df_i, ref_data, by = c("sample", "rep")) |>
    mutate(delta_Ct = mean_Ct - ref_Ct)
  
  treatment_summary <- combined_data |> 
    group_by(sample) |> 
    summarise(mean_delta_Ct = mean(delta_Ct))
  
  mean_ctrl <- filter(treatment_summary, sample == ctrl) |> 
    pull(mean_delta_Ct)
  
  final_data <- combined_data |> 
    mutate(delta_delta_Ct = delta_Ct - mean_ctrl)|> 
    mutate(rel_conc = 2^-delta_delta_Ct)
  
  final_summ_data <- final_data |> 
    group_by(sample, primer) |> 
    summarise(mean = mean(rel_conc), sd = sd(rel_conc), sem = sd(rel_conc)/sqrt(length(rel_conc)))
  
  
  final_summ_data <- rename(final_summ_data, rel_conc = mean)
  
  #Write fold change results to a table and export
  results_exp <- final_data |> 
    select(sample, rep, primer, delta_Ct, delta_delta_Ct, rel_conc) |> 
    rename(Fold_change = rel_conc)
  
  write_tsv(results_exp, paste(tables_path, "/", names(test_list)[i], "-", experiment, "_fold_change_results", ".txt", sep = ""))
  
  #Statistical Tests and export results
  #pairwise t-tests
  final_data$sample <- factor(final_data$sample)
  ttest_results <- pairwise.t.test(final_data$delta_Ct, final_data$sample, p.adjust.method = "BH", pool.sd = F)
  
  tt_pval <- cbind(" "=rownames(ttest_results$p.value), ttest_results$p.value)
  
  write_tsv(as.data.frame(tt_pval), paste(tables_path, "/", names(test_list)[i], "-", experiment, "_t-test_p-values_BH_correction", ".txt", sep = ""))
  
  #One-way ANOVA Test
  res.aov <- aov(delta_Ct ~ sample, data = final_data)
  summary(res.aov)
  
  #Tukey multiple pairwise comparison
  tuk_results <- tidy(TukeyHSD(res.aov))
  
  write_tsv(as.data.frame(tuk_results), paste(tables_path, "/", names(test_list)[i], "-", experiment, "_tukey_test_results", ".txt", sep = ""))
  
  #Dunnett's test as a post-hoc to compare to control
  final_data$sample <- relevel(final_data$sample, ctrl)
  DT <- DunnettTest(final_data$delta_Ct, final_data$sample)
  DTdf <- as.data.frame(DT[[ctrl]])
  dtpvs <- DTdf |> 
    rownames_to_column("sample") |> 
    select(c("sample", "pval")) |> 
    rename("p_val" = "pval") #creates dataframe with Dunnett's test results to assign significance asterisks to the plot
  dtpvs$sample <- substr(dtpvs$sample, 1, nchar(dtpvs$sample) - (nchar(ctrl) + 1)) #modifies the name of the samples to match the final_summ_data for joining
  
  write_tsv(as.data.frame(dtpvs), paste(tables_path, "/", names(test_list)[i], "-", experiment, "_Dunnetts_test_results", ".txt", sep = ""))
  
  #Join final data with pvals from Dunnett Test
  final_summ_data <- left_join(final_summ_data, dtpvs, by = "sample") #combines the final_summ_data and Dunnett's test results
  
  final_summ_data <- final_summ_data |> #creates a column that will assign the correct significance indicator based on Dunnett's test pvals
    mutate(
      label = case_when(
        p_val > 0.05 ~ "",
        p_val > 0.01 ~ "*",
        p_val > 0.001 ~ "**",
        !is.na(p_val) ~ "***",
        TRUE ~ NA_character_
      )
    )
  
  #Create bar plot of results and export
  p1 <- ggplot(final_data, aes(x = sample, y = rel_conc)) +
    geom_bar(final_summ_data, mapping = aes(x = factor(sample, level = c("NT", "3A", "20A", "111A", "1B", "SE66")), #set the factor level here to specify order of samples on the X-axis
                                            y = rel_conc, fill = sample), stat = "identity", color = "black", width = 0.5) +
    geom_point(final_data, mapping = aes(sample, rel_conc, fill = sample), size = 2, shape = 21, color = "black", alpha = 0.5, stroke = 0.5) + #add points to the bar plot
    geom_errorbar(final_summ_data, mapping = aes(sample, rel_conc, ymin = rel_conc - sem, ymax = rel_conc + sem), color = "black", width = 0.2) + #add error bars
    geom_text(final_summ_data, mapping = aes(label = label), nudge_y = 0.2) + #adds significance asterisks to plot, just above the top of the bars
    scale_y_continuous(expand = expansion(mult = c(0, .1))) + #this adjusts the X and Y axis to fit better
    labs(x = "", y = "Relative Expression", title = paste(names(test_list)[i], "\n Expression", sep = "")) + #plot labels
    theme_Publication() + #apply my personal ggplot theme, this can be skipped or modified
    scale_fill_brewer(palette = "Dark2") + #color fill palette I like. Can be manually set or use a different palette
    scale_color_brewer(palette = "Dark2") + #same as above, but for color outlines
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) #adjustments to remove legend and tilt the X-axis labels
  
  ggsave(paste(plots_path, "/", names(test_list)[i], "-", experiment, ".png", sep = ""), plot = p1, width = 3, height = 4, units = "in", dpi = 600) #save plot as .png
  ggsave(paste(plots_path, "/", names(test_list)[i], "-", experiment, ".pdf", sep = ""), plot = p1, width = 3, height = 4, units = "in", dpi = 600) #save plot as .pdf
  
}