# Libraries ----
required_packages <- c(
  "dplyr",
  "purrr",
  "openxlsx2",
  "ggplot2",
  "ggsignif",
  "parallel",
  "ggalluvial",
  "Polychrome",
  "scales"
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

# Input ----

genomes2samples <- read_xlsx(file.path(input_dir,"Genomes2Samples.xlsx"))
metadata <- readRDS(file.path(input_dir,"MAGIC_thresher_plateau.RDS"))[["strains"]]
swabs_dates <- read_xlsx(file.path(input_dir,"SwabDates.xlsx"))
persistence_list <- readRDS(file.path(input_dir,"persistence_list.RDS"))
persistence_sum <- readRDS(file.path(input_dir,"persistence_sum.RDS"))

# Analysis Functions ----


strain_phylogenomics_persistence <- lapply(persistence_list,
                                           function(list_entry){
                                             if(length(list_entry[["strain_phylogenomics_persistence"]]) == 0){
                                               return(list(
                                                 list(
                                                   subject = unique(list_entry$plot_df$subject),
                                                   strain = "Undetermined",
                                                   persistence = unique(list_entry$plot_df$timepoint[list_entry$plot_df$category == "Phylogenomics"])
                                                 )
                                               ))
                                             }else{
                                               return(list_entry[["strain_phylogenomics_persistence"]])
                                             }
                                           })

# Summarize Persistence pattern ----

sum_colonization_pattern <- function(persistence_sum,
                                     genomes2samples,
                                     metadata,
                                     swabs_dates,
                                     strain_phylogenomics_persistence){
  
  # all subjects
  # 3 subjects (440,442,444) have only stool samples which are not included in swabs data
  all_subjects <- c(setdiff(genomes2samples$subject,swabs_dates$`MAGIC ID`),
                    unique(swabs_dates$`MAGIC ID`))
  
  # culture positive subjects are those with positive swabs, including those with no genome due to QC or failed WGS
  culture_positive_subjects <- sort(unique(genomes2samples$subject))
  
  # culture 
  culture_persistent_subjects <- c()
  culture_intermittent_subjects <- c()
  culture_once_subjects <- c()
  culture_never_subjects <- setdiff(all_subjects,
                                    culture_positive_subjects)
  
  for(current_subject in culture_positive_subjects){
    
    subject_culture_timepoint <- sort(unique(genomes2samples$timepoint[genomes2samples$subject == current_subject]))
    
    if(length(subject_culture_timepoint) == 1){
      
      culture_once_subjects <- c(culture_once_subjects,
                                 current_subject)
      
    }else if(length(subject_culture_timepoint) > 1){
      
      rle <- rle(diff(subject_culture_timepoint))
      
      if(length(rle[rle$values == 1]) == 0){
        culture_intermittent_subjects <- c(culture_intermittent_subjects,
                                           current_subject)
      }else{
        culture_persistent_subjects <- c(culture_persistent_subjects,
                                         current_subject)
      }
    }
  }
  
  # Phylogenomics
  phylogenomics_persistent_subjects <- c()
  phylogenomics_intermittent_subjects <- c()
  phylogenomics_once_subjects <- c()
  # culture Never will always be Phylogenomics Never
  phylogenomics_never_subjects <- culture_never_subjects
  phylogenomics_undetermined_subjects <- c()
  
  
  for(current_subject in culture_positive_subjects) {
    
    subject_phylogenomics_persistence <- persistence_sum[persistence_sum$subject == current_subject & 
                                                           persistence_sum$category == "Phylogenomics",]
    
    if(nrow(subject_phylogenomics_persistence) == 0){
      # if there is no phylogenomics persistence for this subject
      
      if(# If all Failed WGS for this subject, the phylogenomics persistence pattern is phylogenomics undetermined
        all(unique(genomes2samples$genome[genomes2samples$subject == current_subject]) == "Failed WGS") | 
        # If there are both and Failed WGS and Not S.aureus and only Failed WGS, Not S.aureus
        # It's phylogenomics undetermined too
        all(unique(genomes2samples$genome[genomes2samples$subject == current_subject]) %in% c("Failed WGS","Not S.aureus")) &
        all(c("Failed WGS","Not S.aureus") %in% unique(genomes2samples$genome[genomes2samples$subject == current_subject]))){
        
        
        phylogenomics_undetermined_subjects <- c(phylogenomics_undetermined_subjects,
                                                 current_subject)
        
      }else if(all(unique(genomes2samples$genome[genomes2samples$subject == current_subject]) == "Not S.aureus")){
        
        
        # if there is culture positives but they turn out to be not S.aureus
        # and there are no phylogenomics positive time point left
        # this is a phylogenomics negative
        
        phylogenomics_never_subjects <- c(phylogenomics_never_subjects,
                                          current_subject)
        
      }
      
    }else{
      
      subject_entry <- which(sapply(strain_phylogenomics_persistence, function(x) x[[1]]$subject == current_subject))
      
      subject_strain <- strain_phylogenomics_persistence[[subject_entry]]
      
      if(length(subject_strain) == 1){
        
        if(length(subject_strain[[1]]$persistence) == 1){
          
          # If there is a culture presence that is missed by phylogenomics presence due to failed WGS
          # And this failed WGS could potential make it persistent
          # This is undetermined
          
          if(length(setdiff(genomes2samples$timepoint[genomes2samples$subject == current_subject & genomes2samples$genome == "Failed WGS"],
                            subject_strain[[1]]$persistence)) >= 1){
            
            phylogenomics_undetermined_subjects <- c(phylogenomics_undetermined_subjects, current_subject)
          }else{
            phylogenomics_once_subjects <- c(phylogenomics_once_subjects, current_subject)
          }
          
        } else {
          phylogenomics_persistent_subjects <- c(phylogenomics_persistent_subjects, current_subject)
        }
        
      } else {
        strain_timepoints <- sapply(subject_strain, function(strain) strain$persistence)
        all_once <- all(sapply(strain_timepoints, length) == 1)
        
        if(all_once) {
          if(length(unique(unlist(strain_timepoints))) == 1) {
            phylogenomics_once_subjects <- c(phylogenomics_once_subjects, current_subject)
          } else {
            phylogenomics_intermittent_subjects <- c(phylogenomics_intermittent_subjects, current_subject)
          }
        } else {
          phylogenomics_persistent_subjects <- c(phylogenomics_persistent_subjects, current_subject)
        }
      }
    }
  }
  
  return(
    list(
      culture_persistent_subjects = sort(culture_persistent_subjects),
      culture_intermittent_subjects = sort(culture_intermittent_subjects),
      culture_once_subjects = sort(culture_once_subjects),
      culture_never_subjects = sort(culture_never_subjects),
      phylogenomics_persistent_subjects = sort(phylogenomics_persistent_subjects),
      phylogenomics_intermittent_subjects = sort(phylogenomics_intermittent_subjects),
      phylogenomics_once_subjects = sort(phylogenomics_once_subjects),
      phylogenomics_never_subjects = sort(phylogenomics_never_subjects),
      phylogenomics_undetermined_subjects = sort(phylogenomics_undetermined_subjects)
    )
  )
}


# sankey plot ----

get_sankey_plot <- function(colonization_pattern_sum){
  
  # all subjects
  all_subjects <- sort(unique(unlist(sapply(colonization_pattern_sum,
                                            function(category_subjects){
                                              category_subjects
                                            }))))
  
  
  sankey_df <- data.frame(
    category = c(rep("Culture",length(all_subjects)),
                 rep("Phylogenomics",length(all_subjects))),
    subject = all_subjects,
    pattern = NA
  )
  
  sankey_df$pattern <- sapply(seq_len(nrow(sankey_df)),
                             function(row_entry){
                               if(sankey_df$category[row_entry] == "Culture"){
                                 if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$culture_persistent_subjects){
                                   "Persistent"
                                 }else if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$culture_intermittent_subjects){
                                   "Intermittent"
                                 }else if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$culture_once_subjects){
                                   "Once"
                                 }else if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$culture_never_subjects){
                                   "Never"
                                 }
                                 
                               }else{
                                 if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$phylogenomics_persistent_subjects){
                                   "Persistent"
                                 }else if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$phylogenomics_intermittent_subjects){
                                   "Intermittent"
                                 }else if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$phylogenomics_once_subjects){
                                   "Once"
                                 }else if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$phylogenomics_never_subjects){
                                   "Never"
                                 }else if(sankey_df$subject[row_entry] %in% colonization_pattern_sum$phylogenomics_undetermined_subjects){
                                   "Undetermined"
                                 }
                               }
                             })
  
  # Summarize the number of subjects with colonization patterns modified
  modified_subjects <- unlist(sapply(unique(sankey_df$subject),
                                     function(subject_entry){
                                       
                                       culture_pattern <- sankey_df$pattern[sankey_df$category == "Culture" &
                                                                             sankey_df$subject == subject_entry]
                                       phylogenomics_pattern <- sankey_df$pattern[sankey_df$category == "Phylogenomics" &
                                                                                   sankey_df$subject == subject_entry]
                                       
                                       if(culture_pattern != phylogenomics_pattern){
                                         return(subject_entry)
                                       }
                                     }))
  
  print(paste0("Modified Subject Count: ",length(modified_subjects)))
  
  pattern_colors <- setNames(c("#C72373",
                               "#70B4E2",
                               "#819825",
                               "#756757",
                               "#E3E0DC"),
                             c("Persistent",
                               "Intermittent",
                               "Once",
                               "Never",
                               "Undetermined"))
  
  sankey_df$category <- factor(sankey_df$category,
                              levels = c("Culture","Phylogenomics"))
  
  sankey_df$pattern <- factor(sankey_df$pattern,
                             levels = c("Persistent","Intermittent","Once","Never","Undetermined"))
  
  sankey_plot <- ggplot(sankey_df,
                       aes(x = category,
                           stratum = pattern,
                           alluvium = subject)) +
    geom_stratum(aes(fill = pattern),
                 color = "#00000000") +
    geom_text(stat = "stratum",
              aes(label = after_stat(count)),
              fontface = "bold",
              color = "white",
              size = 10) + 
    geom_flow(stat = "alluvium",
              lode.guidance = "frontback",
              aes(
                fill = pattern,
                color = "black"),
              linewidth = 0.5,
              alpha = 0.5) +
    scale_fill_manual(values = pattern_colors) +
    scale_color_manual(values = pattern_colors) +
    guides(fill=guide_legend(title="Colonization Pattern")) + 
    scale_size_identity() +
    labs(x = "Definition",
         y = "Subject Count") +
    scale_y_continuous(expand = c(0, 0),
                       limits = function(x) c(0, max(x))) +
    scale_x_discrete(expand = c(0.15, 0.15))+
    theme(axis.title.x =element_blank(),
          axis.text.y = element_text(angle = 0,
                                     size = 20),
          axis.title.y = element_text(size = 25,
                                      face = "bold"),
          axis.text.x = element_text(size = 20),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 15),
          legend.key.size = unit(1,"cm"),
          legend.position = "bottom",
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
    )
  
  pdf(file="Cohort_Colonization_Patterns.pdf",
      width=6.25,
      height=10)
  print(sankey_plot)
  dev.off()
}

# Persistence length comparison ----
get_persistence_plot <- function(persistence_sum,
                                 swabs_dates,
                                 colonization_pattern_sum){
  # Helper function to get the persistence length for subject
  get_persistence <- function(subject,
                              category,
                              persistence_sum,
                              swabs_dates){
    
    # First find how many colonization periods in this subject
    colonization_period <- unique(persistence_sum$line_connect[persistence_sum$subject == subject &
                                                                 persistence_sum$category == category])
    
    total_persistence <- 0
    
    # For each period, determine whether intermittent or persistent
    for(period in colonization_period){
      period_timepoint <-  unique(persistence_sum$timepoint[persistence_sum$subject == subject &
                                                              persistence_sum$category == category &
                                                              persistence_sum$line_connect == period])
      if(length(period_timepoint) > 1){
        # Only count persistence when the period is persistent
        min_date <- as.Date(swabs_dates$`Date Collected`[swabs_dates$`MAGIC ID` == subject &
                                                           swabs_dates$`Time Point` == min(period_timepoint)])
        
        max_date <- as.Date(swabs_dates$`Date Collected`[swabs_dates$`MAGIC ID` == subject &
                                                           swabs_dates$`Time Point` == max(period_timepoint)])
        
        persistence <- as.integer(max_date - min_date)
        
        total_persistence <- total_persistence + persistence
      }
    }
    return(
      data.frame(
        subject = subject,
        category = category,
        persistence = total_persistence
      )
    )
  }
  
  all_subjects <- sort(unique(unlist(sapply(colonization_pattern_sum,
                                            function(category_subjects){
                                              category_subjects
                                            }))))
  
  # Culture
  
  culture_persistence <- do.call(rbind,
                                 lapply(all_subjects,
                                        function(subject){
                                          get_persistence(subject,
                                                          category = "Culture",
                                                          persistence_sum,
                                                          swabs_dates)
                                        }))
  # Phylogenomics
  phylogenomics_persistence <- do.call(rbind,
                                       lapply(all_subjects,
                                              function(subject){
                                                get_persistence(subject,
                                                                category = "Phylogenomics",
                                                                persistence_sum,
                                                                swabs_dates)
                                              }))
  
  # Helper function to get stats for persistence
  get_persistence_stats <- function(persistence){
    persistence_day <- persistence$persistence[persistence$persistence != 0]
    category <- unique(persistence$category)
    stats_df <- data.frame(
      category = category,
      # how many subjects were persistently colonized using this approach
      persistent_subject_count = length(persistence_day),
      # mean of persistence using this approach
      average_persistence = mean(persistence_day),
      # median of persistence using this approach
      median_persistence = median(persistence_day),
      # SD of persistence of this approach
      sd_persistence = sd(persistence_day),
      # Q1 of persistence of this approach
      q1_persistence = quantile(persistence_day,0.25),
      # Q3 of persistence of this approach
      q3_persistence = quantile(persistence_day,0.75)
    )
    
    return(stats_df)
  }
  
  persistence_stats <- rbind(get_persistence_stats(culture_persistence),
                             get_persistence_stats(phylogenomics_persistence))
  
  
  # Visualization 
  plot_df <- cbind(culture_persistence[,c(1,3)],
                   phylogenomics_persistence[,3])
  
  colnames(plot_df)[c(2,3)] <- c("Culture",
                                 "Phylogenomics")
  
  plot_df$discrepancy <- plot_df$Phylogenomics - plot_df$Culture
  plot_df <- plot_df[!(plot_df$Culture == 0 &
                         plot_df$Phylogenomics == 0),]
  
  persistence_compare_boxplot <- ggplot(data = plot_df) +
    #culture
    geom_violin(aes(x = "Culture",
                    y = Culture),
                alpha = 0.5,
                trim = TRUE,
                fill = "#990000",
                color = "transparent") +
    geom_point(aes(x = "Culture",
                   y = Culture),
               color = "#420D09",
               alpha = 1,
               size = 1) + 
    #Phylogenetic
    geom_violin(aes(x = "Phylogenomics",
                    y = Phylogenomics),
                alpha = 0.5,
                trim = TRUE,
                fill = "#011F5B",
                color = "transparent") +
    geom_point(aes(x = "Phylogenomics",
                   y = Phylogenomics),
               color = "#29465B",
               alpha = 1,
               size = 1) + 
    geom_segment(aes(x = "Culture",
                     y = Culture,
                     xend = "Phylogenomics",
                     yend = Phylogenomics,
                     color = discrepancy),
                 alpha = 0.5,
                 linewidth = 0.75) + 
    scale_color_viridis_c(name = "Discrepancy (Days)",
                          breaks = pretty_breaks(n = 8)) +
    scale_y_continuous() + 
    labs(y = "Persistence/Subject (Days)") +
    theme_minimal() +
    scale_x_discrete(limits = c("Culture",
                                "Phylogenomics")) + 
    theme(axis.text.x = element_text(angle = 0,
                                     size = 25,
                                     face = "bold"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_text(angle = 0,
                                     size = 20),
          axis.title.y = element_text(size = 25,
                                      face = "bold"),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.key.size = unit(1.5,"cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA,
                                      linewidth = 1))
  
  pdf(file=paste0("Violin_Persistence_Comparison.pdf"),
      width=12.5,
      height=5.5)
  print(persistence_compare_boxplot)
  dev.off()
}


# Execute Functions ----


colonization_pattern_sum <- sum_colonization_pattern(persistence_sum,
                                                     genomes2samples,
                                                     metadata,
                                                     swabs_dates,
                                                     strain_phylogenomics_persistence)

get_sankey_plot(colonization_pattern_sum)


get_persistence_plot(persistence_sum,
                     swabs_dates,
                     colonization_pattern_sum)
