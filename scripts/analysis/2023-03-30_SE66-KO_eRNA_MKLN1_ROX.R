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

# Define pathway to the qPCR data
data <- "data/raw/2023-03-30_132510_SE66-KO_eRNA_MKLN1_ROX.xls"

# Name the experiment
experiment <- "SE66-KO_eRNA_MKLN1_ROX"

# Create directories to save output files
plots_path <- file.path(paste0("plots/", Sys.Date(), "_", experiment))

if(!dir.exists(plots_path)){
  dir.create(plots_path)
}

tables_path <- file.path(paste0("tables/", Sys.Date(), "_", experiment))

if(!dir.exists(tables_path)){
  dir.create(tables_path)
}

# Assign name of control reference sample
ctrl <- "WT"

# Assign names of housekeeping normalization gene(s)
ref <- c("GAPDH", "B-Actin")

# Clean up data
qpcr_data <- read_excel(data, sheet = "Results", skip = 44, col_names = TRUE) |> 
  slice(1:(n()-5)) |> 
  dplyr::filter(!`Sample Name` %in% c("NTC", "No Template", "NT")) |> 
  dplyr::filter(!`Sample Name` %in% c("-RT", "RT-")) |>   
  separate(`Well Position`, into = c("row", "column"), sep = 1, convert = TRUE) |> 
  rename(primer = `Target Name`) |> 
  mutate(CT = as.numeric(na_if(CT, "Undetermined"))) #Converts "Undetermined" wells to "NA" and resets class to numeric

# Visualize plate layout
ggplot(qpcr_data, aes(x = column, y = row, fill = primer, label = `Sample Name`)) +
  geom_tile(colour = "black") +
  geom_text() +
  scale_y_discrete(limits = c("P", "O", "N", "M", "L", "K", "J", "I", "H", "G", "F", "E", "D", "C", "B", "A")) +
  scale_x_continuous(breaks = 1:24)

# Summarise and clean up data further
summ_data <- qpcr_data |> 
  group_by(`Sample Name`, primer) |> 
  summarise(mean_Ct = mean(CT, na.rm = TRUE)) 

summ_data <- summ_data |>
  separate(`Sample Name`, into = c("sample", "rep"), sep = " ", convert = TRUE)

ggplot(summ_data, aes(x = sample, y = mean_Ct, color = primer)) +
  geom_point()

# Prepare reference data for analysis
if (length(ref) > 1) {
  ref1_data <- summ_data |> 
    filter(primer == ref[1]) |> 
    rename("ref1_Ct" = "mean_Ct")
  
  ref2_data <- summ_data |> 
    filter(primer == ref[2]) |> 
    rename("ref2_Ct" = "mean_Ct")
  
  ref_data <- left_join(ref1_data, ref2_data, by = c("sample", "rep"))
  
  ref_data$ref_Ct <- rowMeans(subset(ref_data, select = c(ref1_Ct, ref2_Ct)), na.rm = TRUE)
} else {
  ref_data <- summ_data |> 
    filter(primer == ref) |> 
    rename("ref_Ct" = "mean_Ct")
}

# Loop through data and analyze to create plots and perform statistical analysis of results
test_data <- summ_data |> 
  filter(!primer %in% ref)

test_list <- split(test_data, test_data$primer)

for (i in 1:length(test_list)){
  df_i <- test_list[[i]]
  combined_data <- merge(df_i, ref_data, by = c("sample", "rep")) |>
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
  
  #Dunnett's test as a post-hoc to compare to control
  final_data$sample <- relevel(final_data$sample, ctrl)
  DT <- DunnettTest(final_data$delta_Ct, final_data$sample)
  DTdf <- as.data.frame(DT[[ctrl]])
  dtpvs <- DTdf |> 
    rownames_to_column("sample") |> 
    select(c("sample", "pval")) |> 
    rename("p_val" = "pval")
  dtpvs$sample <- substr(dtpvs$sample, 1, nchar(dtpvs$sample) - (nchar(ctrl) + 1))
  
  write_tsv(as.data.frame(dtpvs), paste(tables_path, "/", names(test_list)[i], "-", experiment, "_Dunnetts_test_results", ".txt", sep = ""))
  
  #One-way ANOVA Test
  res.aov <- aov(delta_Ct ~ sample, data = final_data)
  summary(res.aov)
  
  #Tukey multiple pairwise comparison
  tuk_results <- tidy(TukeyHSD(res.aov))
  
  write_tsv(as.data.frame(tuk_results), paste(tables_path, "/", names(test_list)[i], "-", experiment, "_tukey_test_results", ".txt", sep = ""))
  
  #Join final data with pvals from Dunnett Test
  final_summ_data <- left_join(final_summ_data, dtpvs, by = "sample")
  
  final_summ_data <- final_summ_data |> 
    mutate(
      label = case_when(
        p_val > 0.05 ~ "",
        p_val > 0.01 ~ "*",
        p_val > 0.001 ~ "**",
        !is.na(p_val) ~ "***",
        TRUE ~ NA_character_
      )
    )
  
  p1 <- ggplot(final_data, aes(x = sample, y = rel_conc)) +
    geom_bar(final_summ_data, mapping = aes(x = fct_reorder(sample, rel_conc, .desc = TRUE), y = rel_conc, fill = sample), stat = "identity", color = "black", width = 0.5) +
    geom_point(final_data, mapping = aes(sample, rel_conc, fill = sample), size = 2, shape = 21, color = "black", alpha = 0.5, stroke = 0.5) +
    geom_errorbar(final_summ_data, mapping = aes(sample, rel_conc, ymin = rel_conc - sem, ymax = rel_conc + sem), color = "black", width = 0.2) +
    geom_text(final_summ_data, mapping = aes(label = label), nudge_y = 0.2) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    labs(x = "", y = "Relative Expression", title = paste(names(test_list)[i], "\n Expression", sep = "")) +
    theme_Publication() +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste(plots_path, "/", names(test_list)[i], "-", experiment, ".png", sep = ""), plot = p1, width = 3, height = 4, units = "in", dpi = 600)
  ggsave(paste(plots_path, "/", names(test_list)[i], "-", experiment, ".pdf", sep = ""), plot = p1, width = 3, height = 4, units = "in", dpi = 600)
  
}

# Plot comparison of KO eRNA expression to NTC and -RT conditions
qpcr_data_PODe1 <- read_excel(data, sheet = "Results", skip = 44, col_names = TRUE) |> 
  slice(1:(n()-5)) |> 
  separate(`Well Position`, into = c("row", "column"), sep = 1, convert = TRUE) |> 
  rename(primer = `Target Name`) |> 
  mutate(CT = as.numeric(na_if(CT, "Undetermined"))) |>  #Converts "Undetermined" wells to "NA" and resets class to numeric
  filter(primer == "PODXLe1-P") |> 
  separate(`Sample Name`, into = c("sample", "rep"), sep = " ", convert = TRUE)

p2 <- ggplot(qpcr_data_PODe1, mapping = aes(x = fct_reorder(sample, CT, .desc = FALSE), y = CT, color = sample)) +
  geom_point() +
  theme_Publication() +
  labs(x = "") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave(paste(plots_path, "/KO_vs_NTC_RT-", experiment, ".png", sep = ""), plot = p2, width = 4, height = 4, units = "in", dpi = 600)
