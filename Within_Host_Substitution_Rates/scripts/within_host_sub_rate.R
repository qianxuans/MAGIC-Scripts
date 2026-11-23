# Libraries ----
required_packages <- c(
  "ggplot2",
  "openxlsx2",
  "dplyr",
  "patchwork"
)

# Stop with an install command if any package is missing
package_installed <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
missing_packages <- required_packages[!package_installed]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ", paste(missing_packages, collapse = ", "), "\n",
    "Install them with: install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "), "))",
    call. = FALSE
  )
}

# Load all packages
for (package_name in required_packages) {
  library(package_name, character.only = TRUE)
}
# Helper Functions ----

make_title <- function(slope, ci, r2, p, unit = "SNPs/day") {
  f <- function(x, d = 3) formatC(unname(x), format = "f", digits = d)
  bquote(atop(
    "slope = " * .(f(slope)) ~ .(unit) * " (95% CI " * .(f(ci[1])) * ", " * .(f(ci[2])) * ")",
    "adj. R"^2 ~ "=" ~ .(f(r2)) * "," ~ italic("p") ~ "=" ~ .(formatC(p, format = "g", digits = 2))
  ))
}

# Set the working directory to the folder containing this script
get_script_dir <- function() {
  # Run with source(), including the Source button in RStudio
  for (frame in rev(sys.frames())) {
    if (is.character(frame$ofile)) {
      return(dirname(normalizePath(frame$ofile)))
    }
  }
  
  # Run from the command line with Rscript
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    script_path <- gsub("~+~", " ", script_path, fixed = TRUE)  # Rscript encodes spaces as ~+~
    return(dirname(normalizePath(script_path)))
  }
  
  # Run line by line in RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(script_path)) {
      return(dirname(normalizePath(script_path)))
    }
  }
  
  stop("Could not find the script location. Save the script, then run it with Rscript, source(), or RStudio.")
}

script_dir <- get_script_dir()
parent_dir <- dirname(script_dir)
input_dir <- file.path(parent_dir,"input")
output_dir <- file.path(parent_dir,"output")
setwd(output_dir)

# Input ----
strains <- readRDS(file.path(input_dir,"MAGIC_thresher_plateau.RDS"))[["strains"]]

reference <- rbind(read.table(file.path(input_dir,"subject_strain_reference.txt"),
                              header = TRUE),
                   read.table(file.path(input_dir,"all_subject_strain_reference.txt"),
                              header = TRUE))

swab_dates <- read_xlsx(file.path(input_dir,"SwabDates.xlsx"))
snps_tab_dir <- file.path(input_dir,"snps_tab")

# Analysis ----

all_snp_tabs <- do.call(rbind,
                        lapply(seq_len(nrow(reference)),
                               function(row_entry){
                                 
                                 strain_genomes <- strains$genome[strains$new_strain_id == reference$strain[row_entry] &
                                                                    strains$subject == reference$subject[row_entry]]
                                 
                                 strain_genomes <- strain_genomes[strain_genomes != reference$reference[row_entry]]
                                 
                                 snps_genomes <- do.call(rbind,
                                                         lapply(strain_genomes,
                                                                function(genome){
                                                                  snps_tab <- read.delim(file.path(snps_tab_dir,
                                                                                                   paste0("strain",
                                                                                                          reference$strain[row_entry],
                                                                                                          "_subject",
                                                                                                          reference$subject[row_entry],
                                                                                                          "_",
                                                                                                          genome,
                                                                                                          "_snps.tab")))
                                                                  if(nrow(snps_tab) > 0){
                                                                    snps_tab <- cbind(strain = reference$strain[row_entry],
                                                                                      subject = reference$subject[row_entry],
                                                                                      genome = genome,
                                                                                      snps_tab)
                                                                    return(snps_tab)
                                                                  }
                                                                }))
                                 
                                 return(snps_genomes)
                                 
                               }))

strain_snp_tabs <- all_snp_tabs %>% 
  select(strain,
         subject,
         CHROM,
         POS,
         TYPE,
         EFFECT,
         LOCUS_TAG,
         GENE,
         PRODUCT) %>% 
  unique()


unique_snp_persistence <- do.call(rbind,
                                  lapply(seq_len(nrow(reference)),
                                         function(row_entry){
                                           
                                           unique_snps <- all_snp_tabs %>% 
                                             filter(strain == reference$strain[row_entry] &
                                                      subject == reference$subject[row_entry]) %>%
                                             select(CHROM,POS) %>%
                                             unique() %>%
                                             nrow()
                                           
                                           strain_subject_timepoint <- sort(unique(strains$timepoint[strains$new_strain_id == reference$strain[row_entry] &
                                                                                                       strains$subject == reference$subject[row_entry]]))
                                           
                                           # time point grouping based on consecutive differences
                                           timepoint_group <- setNames(strain_subject_timepoint,
                                                                       cumsum(c(1, diff(strain_subject_timepoint) != 1)))
                                           
                                           most_recent_date <- as.Date(swab_dates$`Date Collected`[swab_dates$`MAGIC ID` == reference$subject[row_entry] &
                                                                                                     swab_dates$`Time Point` == max(strain_subject_timepoint)])
                                           
                                           earliest_date <- as.Date(swab_dates$`Date Collected`[swab_dates$`MAGIC ID` == reference$subject[row_entry] &
                                                                                                  swab_dates$`Time Point` == min(strain_subject_timepoint)])
                                           
                                           strain_subject_culture_persistence <- 0
                                           
                                           for(group_entry in unique(names(timepoint_group))){
                                             
                                             group_most_recent_date <- as.Date(swab_dates$`Date Collected`[swab_dates$`MAGIC ID` == reference$subject[row_entry] &
                                                                                                             swab_dates$`Time Point` == max(timepoint_group[names(timepoint_group) == group_entry])])
                                             
                                             group_earliest_date <- as.Date(swab_dates$`Date Collected`[swab_dates$`MAGIC ID` == reference$subject[row_entry] &
                                                                                                          swab_dates$`Time Point` == min(timepoint_group[names(timepoint_group) == group_entry])])
                                             
                                             strain_subject_culture_persistence <- strain_subject_culture_persistence + as.integer(group_most_recent_date - group_earliest_date)
                                           }
                                           
                                           
                                           
                                           return(
                                             data.frame(
                                               strain = reference$strain[row_entry],
                                               subject = reference$subject[row_entry],
                                               culture_persistence = strain_subject_culture_persistence,
                                               phylogenomic_persistence = as.integer(most_recent_date - earliest_date),
                                               unique_snp = unique_snps
                                             )
                                           )
                                           
                                         }))

# Linear regression

# phylogenomic

phylogenomic_reg_model <- lm(unique_snp ~ phylogenomic_persistence,
                             data = unique_snp_persistence)

phylogenomic_r_squared <- summary(phylogenomic_reg_model)$adj.r.squared
phylogenomic_p_value <- summary(phylogenomic_reg_model)$coefficients["phylogenomic_persistence", "Pr(>|t|)"]
phylogenomic_slope <- coef(phylogenomic_reg_model)["phylogenomic_persistence"]
phylogenomic_slope_ci <- confint(phylogenomic_reg_model, level = 0.95)["phylogenomic_persistence", ]


phylogenomic_unique_snp_persistence_plot <- ggplot(unique_snp_persistence,
                                                   aes(x = phylogenomic_persistence,
                                                       y = unique_snp)) + 
  geom_point(color = "#011F5B",
             fill = "#29465B",
             size = 5,
             alpha = 0.35) + 
  geom_smooth(method = "lm",
              linewidth = 3,
              color = "#011F5B") + 
  scale_x_continuous(limits = c(0, 800)) +
  labs(x = "Phylogenomics Persistence\n/Strain/Subject (Days)",
       y = "Unique SNP",
       title = make_title(phylogenomic_slope, phylogenomic_slope_ci,
                          phylogenomic_r_squared, phylogenomic_p_value)) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(hjust = 0.5, 
                                  size = 20),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1,
                                    fill = "transparent"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        legend.text = element_text(size = 17.5),
        legend.title = element_text(size = 17.5),
        axis.line = element_line(colour = "black",
                                 linewidth = 1),
        axis.ticks.length=unit(0.25,"cm")
        
  )


culture_reg_model <- lm(unique_snp ~ culture_persistence,
                        data = unique_snp_persistence)

culture_r_squared <- summary(culture_reg_model)$adj.r.squared
culture_p_value <- summary(culture_reg_model)$coefficients["culture_persistence", "Pr(>|t|)"]
culture_slope <- coef(culture_reg_model)["culture_persistence"]
culture_slope_ci <- confint(culture_reg_model, level = 0.95)["culture_persistence", ]

culture_unique_snp_persistence_plot <- ggplot(unique_snp_persistence,
                                              aes(x = culture_persistence,
                                                  y = unique_snp)) + 
  geom_point(color = "#990000",
             size = 5,
             alpha = 0.35,
             fill = "#420D09") + 
  geom_smooth(method = "lm",
              linewidth = 3,
              color = "#990000") + 
  scale_x_continuous(limits = c(0, 800)) +
  labs(x = "Culture Persistence/Strain (Days)",
       y = "Unique SNP",
       title = make_title(culture_slope, culture_slope_ci,
                          culture_r_squared, culture_p_value)) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(hjust = 0.5, 
                                  size = 20),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(linewidth = 1,
                                    fill = "transparent"),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key = element_rect(colour = "transparent",
                                  fill = "transparent"),
        legend.text = element_text(size = 17.5),
        legend.title = element_text(size = 17.5),
        axis.line = element_line(colour = "black",
                                 linewidth = 1),
        axis.ticks.length=unit(0.25,"cm")
        
  )

combined_plot <- culture_unique_snp_persistence_plot + phylogenomic_unique_snp_persistence_plot +
  plot_layout(nrow=1)



pdf(file=paste0("Persistence_Unique_Snp_Comparison.pdf"),width = 12, height= 7.5)
print(combined_plot)
dev.off()
