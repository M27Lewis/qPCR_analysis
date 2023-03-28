library(tidyverse)
library(readxl)
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

data_path <- "data/raw/2023-03-28_ROX_vs_noROX_data.xls"

data <- read_excel(data_path, sheet = "Results", skip = 44, col_names = TRUE) |> 
  slice(1:(n()-5)) |> 
  dplyr::filter(!`Sample Name` %in% c("NTC", "No Template", "NT")) |> 
  dplyr::filter(!`Sample Name` %in% c("-RT", "RT-")) |>   
  separate(`Well Position`, into = c("row", "column"), sep = 1, convert = TRUE) |> 
  rename(primer = `Target Name`)

# Statistical testing of the results
tt_results <- t.test(data$`Ct SD` ~ data$Dye)
tt_results$p.value

# Create summarized dataframe and add significance indicator for plot
summ_data <- data |> 
  group_by(Dye) |> 
  summarise(mean = mean(`Ct SD`), sd = sd(`Ct SD`), sem = sd(`Ct SD`)/sqrt(length(`Ct SD`))) |> # Summarize the SDs from the data
  mutate(label = c("***", "")) # Add significance indicators based on results of the t-test

# Box plot of the variation in ROX vs noROX results
p1 <- ggplot(data, mapping = aes(x = Dye, y = `Ct SD`, fill = Dye)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(color = "black") +
  geom_text(data = summ_data, mapping = aes(label = label), nudge_y = 0.4, size = 6) +
  theme_Publication() +
  theme(legend.position = "none")

p1

# Create new folder to save output
newpath <- file.path(paste0("plots/", Sys.Date(), "_ROX_vs_noROX_comparison"))

if(!dir.exists(newpath)){
  dir.create(newpath)
}

ggsave(paste(newpath, "/ROX_vs_noROX_results_comparison.png", sep = ""), plot = p1, width = 3, height = 4, units = "in", dpi = 600) # Save plot


