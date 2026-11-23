# Libraries ----
required_packages <- c(
  "ggplot2",
  "openxlsx2",
  "patchwork",
  "ggExtra",
  "purrr",
  "ggrepel",
  "dplyr"
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

# Helper function to convert time point to days after birth
convert_tp_to_days_after_birth <- function(subject,
                                           timepoint,
                                           dob,
                                           swabs_dates){
  
  subject_dob <- as.Date(dob$DOB[dob$MAGIC_ID == subject])
  
  subject_timepoint_days_after_birth <- if(any(swabs_dates$`MAGIC ID` == subject & swabs_dates$`Time Point` == timepoint)){
    
    timepoint_date <- swabs_dates$`Date Collected`[swabs_dates$`MAGIC ID` == subject &
                                                     swabs_dates$`Time Point` == timepoint]
    
    timepoint_days_after_birth <- as.integer(timepoint_date - subject_dob)
    
    return(timepoint_days_after_birth)
  }else{
    
    timepoint_days_after_birth <- switch(
      as.character(timepoint),
      "1" = 0,
      "2" = 1 * 7,
      "3" = 2 * 7,
      "4" = 3 * 7,
      "5" = 4 * 7,
      "6" = 3 * 30,
      "7" = 6 * 30,
      "8" = 9 * 30,
      "9" = 12 * 30,
      "10" = 15 * 30,
      "11" = 18 * 30,
      "12" = 21 * 30,
      "13" = 24 * 30
    )
    
    return(timepoint_days_after_birth)
  }
}

# Input ----
abx_data <- read_xlsx(file.path(input_dir,"abx_data.xlsx"))
abx_class <- read_xlsx(file.path(input_dir,"abx_class.xlsx"))
positive_swabs <- read_xlsx(file.path(input_dir,"Genomes2Samples.xlsx"))
swabs_dates <- read_xlsx(file.path(input_dir,"SwabDates.xlsx"))
# Participant dates of birth are withheld as protected health information (PHI)
# dob <- read_xlsx(file.path(input_dir,"DOB.xlsx"))
colonization_timepoint <- readRDS("ersistence_sum.RDS")

# Function to visualize colonization and antibiotics usage ----

abx_visual <- function(abx_data,
                       abx_class,
                       swabs_dates,
                       dob,
                       colonization_timepoint,
                       plot_max_days,
                       plot_height,
                       plot_width){
  
  
  # Classify the antibiotics to S.aureus active or others
  
  active_abx <- c(abx_class$med_name[abx_class$simple_generic_name == "Mupirocin"],
                  as.character(na.omit(abx_class$med_name[abx_class$`Staph Acive` == "Yes"])))
  
  
  other_abx <- setdiff(abx_class$med_name,
                       active_abx)
  swabs_dates <- swabs_dates[!is.na(swabs_dates$`Date Collected`),]
  filtered_subjects <- sort(unique(swabs_dates$`MAGIC ID`))
  
  category_color_colors <- c("Culture" = "#990000",
                             "Phylogenomics" = "#011F5B")
  
  abx_fill_colors <- c("S.aureus Active" = "#B2C06E",
                       "Others" = "#C6BFBA")
  
  abx_color_colors <- c("S.aureus Active" = "#829726",
                        "Others" = "#766657")
  
  plot_list <- lapply(filtered_subjects, function(subject_entry){
    
    # Generate the data frame for the plot
    positive_timepoint <- colonization_timepoint %>%
      filter(subject == subject_entry) %>%
      select(subject, timepoint, category, line_connect) %>%
      mutate(days_after_birth = sapply(timepoint, function(timepoint_entry){
        convert_tp_to_days_after_birth(subject_entry, 
                                       timepoint_entry, 
                                       dob, 
                                       swabs_dates)})
      )
    
    # Negative swabs 
    swab_timepoint <- swabs_dates$`Time Point`[swabs_dates$`MAGIC ID` == subject_entry]
    
    phylogenomics_negative_timepoint <- setdiff(swab_timepoint, positive_timepoint$timepoint[positive_timepoint$category == "Phylogenomics"])
    
    culture_negative_timepoint <- setdiff(swab_timepoint, positive_timepoint$timepoint[positive_timepoint$category == "Culture"])
    # Phylogenomics Negative
    phylogenomics_negative_timepoint_df <- if(length(phylogenomics_negative_timepoint)>0){
      data.frame(
        subject = subject_entry,
        timepoint = phylogenomics_negative_timepoint,
        category = "Phylogenomics",
        days_after_birth = sapply(phylogenomics_negative_timepoint, 
                                  function(timepoint_entry) 
                                    convert_tp_to_days_after_birth(subject_entry, 
                                                                   timepoint_entry, 
                                                                   dob, 
                                                                   swabs_dates))
      )
    }else{
      data.frame()
    }
    # Clinical Negative
    culture_negative_timepoint_df <- if(length(culture_negative_timepoint)>0){
      data.frame(
        subject = subject_entry,
        timepoint = culture_negative_timepoint,
        category = "Culture",
        days_after_birth = sapply(culture_negative_timepoint, 
                                  function(timepoint_entry) 
                                    convert_tp_to_days_after_birth(subject_entry, 
                                                                   timepoint_entry, 
                                                                   dob, 
                                                                   swabs_dates))
      )
    }else{
      data.frame()
    }
    
    negative_timepoint <- rbind(
      phylogenomics_negative_timepoint_df,
      culture_negative_timepoint_df
    )
    
    # Find the antibiotics usage period
    
    subject_abx <- abx_data[abx_data$magic_id == subject_entry &
                              abx_data$abx_med_class == "Antibacterial",
                            c("magic_id",
                              "med_name",
                              "abx_med_class",
                              "age_at_start",
                              "age_at_stop")]
    
    
    if(nrow(subject_abx)>=1){
      subject_abx$age_at_stop <- sapply(subject_abx$age_at_stop,
                                        function(stop_day){
                                          if(stop_day >= plot_max_days | is.na(stop_day)) plot_max_days else stop_day})
    } 
    
    
    
    subject_abx$active_class <- sapply(subject_abx$med_name,
                                       function(med_name){
                                         if(med_name %in% active_abx){
                                           "S.aureus Active"
                                         }else if(med_name %in% other_abx){
                                           "Others"
                                         }
                                       })
    # Visualize the colonization and antibiotics usage
    
    # background color, if even grey, If odd transparent
    background_color <- if(which(filtered_subjects == subject_entry) %% 2 == 0) "lightgrey" else "transparent"
    
    if(subject_entry == last(filtered_subjects)){
      x_axis_line_config <- element_line()
      x_axis_ticks_config <- element_line()
      x_axis_text_config <- element_text(size = 20)
      x_axis_title_config <- element_text(size = 25,
                                          face = "bold")
    }else{
      x_axis_line_config <- element_blank()
      x_axis_ticks_config <- element_blank()
      x_axis_text_config <- element_blank()
      x_axis_title_config <- element_blank()
    }
    
    
    subject_abx_plot <- ggplot(data = positive_timepoint) + 
      # antibiotics block
      geom_rect(data = subject_abx,
                aes(xmin = age_at_start,
                    xmax = age_at_stop,
                    ymin = -Inf,
                    ymax = Inf,
                    fill = active_class,
                    color = active_class),
                alpha = 0.65,
                linewidth = 0.75) + 
      # positive point 
      geom_point(aes(x = days_after_birth,
                     y = interaction(category,
                                     subject),
                     color = category),
                 size = 2,
                 alpha = 1) +
      # phylogenetic positive line 
      geom_line(data = subset(positive_timepoint, category == "Phylogenomics"),
                aes(x = days_after_birth, 
                    y = interaction(category,
                                    subject),
                    group = interaction(subject,
                                        line_connect),
                    color = category),
                linewidth = 0.75,
                alpha = 1) +
      # clinical positive line
      geom_line(data = subset(positive_timepoint, category == "Culture"),
                aes(x = days_after_birth, 
                    y = interaction(category,
                                    subject),
                    group = interaction(subject,
                                        line_connect), 
                    color = category),
                linewidth = 0.75,
                alpha = 1) +
      # negative point
      {if(nrow(negative_timepoint)>0){
        geom_point(data = negative_timepoint,
                   aes(x = days_after_birth,
                       y = interaction(category,
                                       subject)),
                   shape = 13,
                   size = 2,
                   alpha = 1)}} +
      scale_color_manual(values = c(category_color_colors,
                                    abx_color_colors),
                         guide = "none") + 
      scale_fill_manual(values = abx_color_colors,
                        guide = "none") + 
      scale_x_continuous(name = "Days After Birth",
                         limits = c(0,plot_max_days)) + 
      scale_y_discrete(name = subject_entry) + 
      theme(axis.text.x = x_axis_text_config,
            axis.line.x = x_axis_line_config,
            axis.text.y = element_blank(),
            axis.title.y = element_text(size = 20,
                                        angle = 0,
                                        vjust = 0.5),
            axis.title.x = x_axis_title_config,
            axis.ticks.y = element_blank(),
            axis.ticks.x = x_axis_ticks_config,
            plot.background = element_blank(),
            panel.background = element_rect(fill = background_color),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_blank(),
            legend.title = element_blank(),
            legend.key = element_blank()
      )
    
    return(subject_abx_plot)
  })
  
  abx_plot <- reduce(plot_list, `/`) +
    plot_layout(heights = rep(1, length(plot_list))) +
    theme(plot.margin = margin(0, 0, 0, 0)
    )
  
  pdf(file="Antibitiocs_Visual.pdf",
      width=plot_width,
      height=plot_height)
  print(abx_plot)
  dev.off()
}

# Execute the function ----
abx_visual(abx_data,
           swabs_dates,
           dob,
           colonization_timepoint,
           plot_max_days = 810,
           plot_height = 60,
           plot_width = 17.5)
