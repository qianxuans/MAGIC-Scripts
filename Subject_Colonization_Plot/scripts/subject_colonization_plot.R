#visualization of the culture persistence and phylogenomics persistence
# Libraries ----
required_packages <- c(
  "dplyr",
  "purrr",
  "ggExtra",
  "openxlsx2",
  "ggplot2",
  "colorspace",
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


# Helper functions ----
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

# Groups strain in the plot later
find_groups <- function(data_matrix) {
  # Get dimensions
  n_members <- nrow(data_matrix)
  n_times <- ncol(data_matrix)
  
  # Initialize adjacency matrix (who is connected to whom)
  adjacency <- matrix(0, n_members, n_members)
  
  # For each time point
  for(t in 1:n_times) {
    # Find members with value 1 at this time point
    active_members <- which(!is.na(data_matrix[,t]) & data_matrix[,t] == 1)
    
    # If there are multiple active members, they form connections
    if(length(active_members) > 1) {
      for(i in 1:(length(active_members)-1)) {
        for(j in (i+1):length(active_members)) {
          # Mark that these members are connected
          adjacency[active_members[i], active_members[j]] <- 1
          adjacency[active_members[j], active_members[i]] <- 1
        }
      }
    }
  }
  
  # Convert adjacency matrix to a distance matrix for clustering
  # (0 means connected, 1 means not connected)
  dist_matrix <- as.dist(1 - adjacency)
  
  # Perform single-linkage hierarchical clustering
  hc <- hclust(dist_matrix, method = "single")
  
  # Determine groups (clusters)
  # We need to cut the dendrogram where distance > 0
  # to separate disconnected components
  groups <- cutree(hc, h = 0.5)
  
  # Return result as a list of groups
  result <- list()
  for(g in unique(groups)) {
    result[[g]] <- which(groups == g)
  }
  
  return(result)
}

script_dir <- get_script_dir()
parent_dir <- dirname(script_dir)
input_dir <- file.path(parent_dir,"input")
output_dir <- file.path(parent_dir,"output")
setwd(output_dir)

# Input ----
metadata <- readRDS(file.path(input_dir,"MAGIC_thresher_plateau.RDS"))[["strains"]]
swabs_dates <- read_xlsx(file.path(input_dir,"SwabDates.xlsx"))
genomes2samples <- read_xlsx(file.path(input_dir,"Genomes2Samples.xlsx"))
# Persistence Plot for individual and all subjects  -----
persistence_plot <- function(metadata,
                             swabs_dates,
                             genomes2samples){
  
  # Positive subjects are those with positive swabs, including those with no S. aureus genome due to QC or failed WGS
  positive_subjects <- sort(unique(genomes2samples$subject))
  
  # all reasons to fail a time point
  fail_reasons <- c("Not S.aureus","Failed WGS","Contaminated S.aureus")
  
  persistence_list <- lapply(positive_subjects,
                             function(current_subject){
                               # timepoint with genomes 
                               metadata_subject <- metadata[metadata$subject == current_subject,]
                               # Add the timepoints failed WGS if any
                               if(any(fail_reasons %in% genomes2samples$genome[genomes2samples$subject == current_subject])){
                                 
                                 failed_wgs_subject <- data.frame(
                                   strain_id = "Undetermined",
                                   genome = "Undetermined",
                                   new_strain_id = "Undetermined",
                                   subject = current_subject,
                                   bodylocation = genomes2samples$bodysite[genomes2samples$subject == current_subject & genomes2samples$genome %in% fail_reasons],
                                   timepoint = genomes2samples$timepoint[genomes2samples$subject == current_subject & genomes2samples$genome %in% fail_reasons],
                                   collection_date = NA)
                                 
                                 for(failed_wgs_entry in seq_len(nrow(failed_wgs_subject))){
                                   
                                   current_tp <- failed_wgs_subject$timepoint[failed_wgs_entry]
                                   current_bodyloc <- failed_wgs_subject$bodylocation[failed_wgs_entry]
                                   failed_wgs_subject$new_strain_id[failed_wgs_entry] <- genomes2samples$genome[genomes2samples$subject == current_subject &
                                                                                                                  genomes2samples$bodysite == current_bodyloc &
                                                                                                                  genomes2samples$timepoint == current_tp]
                                   
                                   # collection date
                                   entry_collection_date <- swabs_dates$`Date Collected`[
                                     swabs_dates$`MAGIC ID` == current_subject & 
                                       swabs_dates$`Time Point` == current_tp
                                   ]
                                   
                                   entry_collection_date <- paste0("20", format(as.Date(entry_collection_date), "%y-%m-%d"))
                                   
                                   failed_wgs_subject$collection_date[failed_wgs_entry] <- entry_collection_date
                                   
                                 }
                                 
                                 
                                 
                                 metadata_subject <- rbind(metadata_subject,
                                                           failed_wgs_subject)
                               }
                               
                               metadata_subject$timepoint <- as.integer(metadata_subject$timepoint)
                               
                               # strain id
                               subject_strain_id <- unique(metadata_subject$new_strain_id[!metadata_subject$new_strain_id %in% fail_reasons])
                               #color of the strains
                               
                               strain_colors <- setNames(c(qualitative_hcl(length(subject_strain_id), palette = "Warm"),colorRampPalette(c("#E3E0DC", "#333333"))(length(fail_reasons))),
                                                         c(subject_strain_id,fail_reasons))
                               
                               # Data frame showing whether the samples were collected at specific time point
                               swabs_timepoints <- do.call(rbind,
                                                           lapply(seq_len(13),
                                                                  function(timepoint){
                                                                    return(
                                                                      data.frame(
                                                                        patient = current_subject,
                                                                        timepoint = timepoint,
                                                                        if_sample = timepoint %in% c(swabs_dates$`Time Point`[swabs_dates$`MAGIC ID` == current_subject],
                                                                                                     metadata_subject$timepoint[metadata_subject$subject == current_subject]),
                                                                        body_location = c("Axilla","Nares","Groin")
                                                                      )
                                                                    )
                                                                  }))
                               # Stool samples if applicable
                               swabs_timepoints <- if(any(metadata_subject$bodylocation == "Stool")){
                                 rbind(swabs_timepoints,
                                       data.frame(
                                         patient = current_subject,
                                         timepoint = metadata_subject$timepoint[metadata_subject$bodylocation == "Stool"],
                                         if_sample = TRUE,
                                         body_location = "Stool"
                                       ))
                               }else{
                                 swabs_timepoints
                               }
                               #check if there is overlapping, then use position_jitter()
                               overlap_timepoints <- do.call(rbind,
                                                             lapply(unique(metadata_subject$timepoint),
                                                                    function(timepoint_entry){
                                                                      meta_timepoint <- metadata_subject[metadata_subject$timepoint == timepoint_entry,]
                                                                      # At this time point, check all body locations to see if there is overlap
                                                                      location_overlap <- sapply(unique(swabs_timepoints$body_location),
                                                                                                 function(location){
                                                                                                   location_strains <- unique(meta_timepoint$strain[meta_timepoint$bodylocation == location])
                                                                                                   if(length(location_strains) > 1){
                                                                                                     return(TRUE)
                                                                                                   }else{
                                                                                                     return(FALSE)
                                                                                                   }
                                                                                                 })
                                                                      
                                                                      return(data.frame(timepoint = timepoint_entry,
                                                                                        body_location = unique(swabs_timepoints$body_location),
                                                                                        if_overlap = location_overlap) %>% filter(if_overlap))
                                                                    }))
                               
                               
                               # data frame for jitter dots 
                               jitter_dots <- if(nrow(overlap_timepoints) > 0){
                                 do.call(rbind,
                                         lapply(seq_len(nrow(overlap_timepoints)),
                                                function(row_entry){
                                                  overlap_timepoint <- overlap_timepoints[row_entry,]$timepoint
                                                  overlap_body_location <- overlap_timepoints[row_entry,]$body_location
                                                  
                                                  data.frame(
                                                    timepoint = overlap_timepoint,
                                                    body_location = overlap_body_location,
                                                    strain = metadata_subject$new_strain_id[metadata_subject$timepoint == overlap_timepoint &
                                                                                              metadata_subject$bodylocation == overlap_body_location]
                                                  )
                                                  
                                                }))
                               }else{
                                 data.frame()
                               }
                               
                               
                               # data frame for regular dots 
                               regular_dots <- do.call(rbind,
                                                       lapply(unique(metadata_subject$timepoint),
                                                              function(timepoint){
                                                                
                                                                timepoint_body_locations <- unique(metadata_subject$bodylocation[metadata_subject$timepoint == timepoint])
                                                                
                                                                regular_dots_tmp <- do.call(rbind,
                                                                                            lapply(timepoint_body_locations,
                                                                                                   function(location){
                                                                                                     
                                                                                                     strain <- if(nrow(overlap_timepoints) > 0){
                                                                                                       if(overlap_timepoints$timepoint == timepoint & overlap_timepoints$body_location == location){
                                                                                                         NA
                                                                                                       }else{
                                                                                                         unique(metadata_subject$new_strain_id[metadata_subject$bodylocation == location &
                                                                                                                                                 metadata_subject$timepoint == timepoint])
                                                                                                       }
                                                                                                     }else{
                                                                                                       unique(metadata_subject$new_strain_id[metadata_subject$bodylocation == location &
                                                                                                                                               metadata_subject$timepoint == timepoint])
                                                                                                     }
                                                                                                     data.frame(
                                                                                                       timepoint = timepoint,
                                                                                                       body_location = location,
                                                                                                       strain = strain
                                                                                                     )
                                                                                                   }))
                                                                
                                                              })) %>% 
                                 filter(!is.na(strain))
                               
                               # data frame for S.aureus negative samples 
                               negative_dots <- do.call(rbind,
                                                        lapply(unique(swabs_timepoints$timepoint[swabs_timepoints$if_sample]),
                                                               function(timepoint){
                                                                 negative_dots_tmp <- do.call(rbind,
                                                                                              lapply(unique(swabs_timepoints$body_location[swabs_timepoints$timepoint == timepoint]),
                                                                                                     function(location){
                                                                                                       if(any(metadata_subject$bodylocation == location & 
                                                                                                              metadata_subject$timepoint == timepoint)){
                                                                                                         return(
                                                                                                           data.frame(
                                                                                                             timepoint = 0,
                                                                                                             body_location = location,
                                                                                                             result = "Positive"
                                                                                                           )
                                                                                                         )
                                                                                                       }else{
                                                                                                         return(data.frame(
                                                                                                           timepoint = timepoint,
                                                                                                           body_location = location,
                                                                                                           result = "Negative"
                                                                                                         ))
                                                                                                       }
                                                                                                     }))
                                                                 return(negative_dots_tmp)
                                                               }))
                               
                               regular_dots$strain <- factor(regular_dots$strain,
                                                             levels = c(sort(unique(regular_dots$strain[!regular_dots$strain %in% fail_reasons])),
                                                                        fail_reasons))
                               
                               jitter_dots$strain <- factor(jitter_dots$strain,
                                                            levels = c(sort(unique(jitter_dots$strain[!jitter_dots$strain %in% fail_reasons])),
                                                                       fail_reasons))
                               
                               negative_dots$result <- factor(negative_dots$result,
                                                              levels = c("Positive","Negative"))
                               
                               isolates_plot <- ggplot() + 
                                 geom_point(data = regular_dots,
                                            aes(x = timepoint,
                                                y = body_location,
                                                colour = strain),
                                            size = 15,
                                            alpha = 0.75) +
                                 geom_point(data = negative_dots,
                                            aes(x = timepoint,
                                                y = body_location,
                                                shape = result),
                                            size = 15,
                                            alpha = 0.75) +
                                 scale_shape_manual(name = "Sampling Result",
                                                    values = setNames(c(1,13),c("Positive","Negative"))) +
                                 {if(nrow(jitter_dots)>0){
                                   geom_point(data = jitter_dots,
                                              aes(x = timepoint,
                                                  y = body_location,
                                                  color = strain),
                                              size = 15,
                                              alpha = 0.75,
                                              position = position_jitterdodge(jitter.width = 0,
                                                                              dodge.width = 0.5))
                                 }} +
                                 guides(shape = guide_legend(order = 1),
                                        color = guide_legend(order = 2)) + 
                                 geom_vline(xintercept = 0.5:13.5, 
                                            color = "grey80", 
                                            linewidth = 0.5) +
                                 geom_hline(yintercept = 0.5:3.5, 
                                            color = "grey80", 
                                            linewidth = 0.5) +
                                 scale_x_continuous(name = "Time Point",
                                                    breaks = seq(13),
                                                    limits = c(1,13)) + 
                                 scale_y_discrete(name = "Body Location",
                                                  expand = c(0,0.5),
                                                  limits = if("Stool" %in% metadata_subject$bodylocation) 
                                                    c("Stool","Groin","Axilla","Nares") 
                                                  else c("Groin","Axilla","Nares")) + 
                                 scale_color_manual(values = strain_colors, 
                                                    name = "Strain") + 
                                 theme(axis.title.x = element_blank(),
                                       axis.text.y = element_text(angle = 0,
                                                                  size = 20,
                                                                  face = "bold"),
                                       axis.title.y = element_blank(),
                                       axis.ticks = element_blank(),
                                       axis.text.x = element_blank(),
                                       plot.background = element_blank(),
                                       panel.background = element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       legend.text = element_text(size = 15),
                                       legend.title = element_text(size = 20),
                                       legend.key.size = unit(1,"cm"),
                                       panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
                                 ) + coord_fixed()
                               
                               # Persistence plot
                               
                               # Determine the phylogenomics persistence/presence time points for this patient
                               
                               phylogenomics_persistence_timepoint <- unique(as.integer(unlist(sapply(subject_strain_id,
                                                                                                      function(strain){
                                                                                                        strain_presence <- unique(metadata_subject$timepoint[metadata_subject$new_strain_id == strain])
                                                                                                        strain_persistence <- min(strain_presence):max(strain_presence)
                                                                                                        return(strain_persistence)
                                                                                                      }))))
                               # Add the time points for contaminated S.aureus genomes 
                               phylogenomics_persistence_timepoint <- sort(c(phylogenomics_persistence_timepoint,
                                                                             unique(metadata_subject$timepoint[metadata_subject$new_strain_id == "Contaminated S.aureus"])))
                               
                               # Determine the culture persistence/presence for this patient
                               persistence_df <- rbind(
                                 data.frame(
                                   timepoint = seq_len(13),
                                   category = "Culture",
                                   presence = sapply(seq_len(13),
                                                     function(timepoint){
                                                       if(timepoint %in% unique(metadata_subject$timepoint)) 
                                                         "Culture_Present"
                                                       else
                                                         "Absent"
                                                     })
                                 ),
                                 data.frame(
                                   timepoint = seq_len(13),
                                   category = "Phylogenomics",
                                   presence = sapply(seq_len(13),
                                                     function(timepoint){
                                                       if(timepoint %in% phylogenomics_persistence_timepoint) 
                                                         "Phylogenomics_Present"
                                                       else
                                                         "Absent"
                                                     })
                                 ))
                               
                               persistence_colors <- c("Culture_Present" = "#990000",
                                                       "Phylogenomics_Present" = "#011F5B",
                                                       "Absent" = "transparent")
                               
                               persistence_df$timepoint <- factor(as.character(persistence_df$timepoint),
                                                                  levels = as.character(seq_len(13)))
                               
                               persistence_plot <- ggplot() + 
                                 geom_tile(data = persistence_df,
                                           aes(x = timepoint,
                                               y = category,
                                               fill = presence),
                                           width = 1) + 
                                 scale_fill_manual(values = persistence_colors,
                                                   guide = "none") + 
                                 geom_hline(yintercept = 1.5, 
                                            color = "grey80", 
                                            linewidth = 0.5) +
                                 scale_x_discrete(name = "Time Point") + 
                                 scale_y_discrete(name = "Persistence",
                                                  limits = c("Phylogenomics","Culture")) + 
                                 theme(axis.text.x = element_text(size = 15),
                                       axis.text.y = element_text(angle = 0,
                                                                  size = 20,
                                                                  face = "bold"),
                                       # axis.title.y = element_text(size = 20,face = "bold"),
                                       # hide the y axis title
                                       axis.title.y = element_blank(),
                                       axis.title.x = element_text(size = 25),
                                       axis.ticks = element_blank(),
                                       plot.background = element_blank(),
                                       panel.background = element_blank(),
                                       panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(),
                                       panel.border = element_rect(colour = "black", fill = NA, linewidth = 2)
                                 ) + coord_fixed()
                               
                               combined_plot <- isolates_plot + persistence_plot +
                                 plot_layout(ncol=1) +
                                 plot_annotation(
                                   title = paste0("Subject",current_subject),
                                   theme = theme(
                                     plot.title = element_text(
                                       hjust = 0.5,
                                       size = 25,
                                       face = "bold"
                                     )
                                   )
                                 )
                               
                               # export the plot
                               pdf(file=paste0("Subject",
                                               current_subject,
                                               ".pdf"),
                                   width=15,
                                   height=10)
                               print(combined_plot)
                               dev.off()
                               
                               # return the persistence sum of current subject 
                               
                               subject_strain_phylogenomics_persistence <- lapply(subject_strain_id,
                                                                                  function(strain){
                                                                                    strain_presence <- unique(metadata_subject$timepoint[metadata_subject$new_strain_id == strain])
                                                                                    strain_persistence <- min(strain_presence):max(strain_presence)
                                                                                    return(
                                                                                      list(
                                                                                        subject = current_subject,
                                                                                        strain = strain,
                                                                                        persistence = strain_persistence
                                                                                      )
                                                                                    )
                                                                                  })
                               
                               persistence_df <- do.call(rbind, lapply(c("Culture", "Phylogenomics"),
                                                                       function(category) {
                                                                         data.frame(
                                                                           subject = current_subject,
                                                                           timepoint = seq_len(13),
                                                                           category = category,
                                                                           presence = as.character(persistence_df$presence[persistence_df$timepoint %in% seq_len(13) & 
                                                                                                                             persistence_df$category == category])
                                                                         )
                                                                       })) %>%
                                 filter(presence != "Absent")
                               
                               return(
                                 list(
                                   strain_phylogenomics_persistence = subject_strain_phylogenomics_persistence,
                                   plot_df = persistence_df,
                                   culture_phylogenomics_discrepancy = metadata_subject[metadata_subject$new_strain_id %in% fail_reasons,]
                                 )
                               )
                             })
  
  
  persistence_sum <- do.call(rbind,
                             lapply(persistence_list,function(list){
                               list[["plot_df"]]
                             }))
  
  persistence_sum <- persistence_sum %>%
    arrange(subject, category, timepoint) %>%
    group_by(subject, category) %>%
    #Only connect lines when consecutive for culture persistence
    mutate(line_connect = cumsum(c(1, diff(timepoint) > 1) * (category == "Culture")))
  
  # For phylogenomics persistence, not only when 2 time points are consecutive but also need to be the same strain 
  # Change the line_connect for phylogenomics in persistence_sum
  # This sapply will also return the subjects whose colonization pattern is modified 
  
  # Keep track of subjects in which patterns differ from culture and phylogenomics
  
  # First, add those subjects with no phylogenomics persistence but with culture persistence
  modified_subject <- as.character(unlist(sapply(unique(persistence_sum$subject),
                                                 function(subject_entry){
                                                   subject_presence <- persistence_sum$presence[persistence_sum$subject == subject_entry]
                                                   if(!"Phylogenomics_Present" %in% subject_presence){
                                                     return(subject_entry)
                                                   }else{
                                                     return(NULL)
                                                   }
                                                 })))
  
  
  for(subject in unique(persistence_sum$subject)){
    
    # Find the relevant list entry for this subject
    list_entry <- which(sapply(persistence_list,
                               function(list) {
                                 if(length(list[["strain_phylogenomics_persistence"]])>0){
                                   list[["strain_phylogenomics_persistence"]][[1]]$subject == subject
                                 }else{
                                   FALSE
                                 }
                               }))
    
    if(length(list_entry) == 1){
      
      subject_strain_index <- seq_along(persistence_list[[list_entry]][["strain_phylogenomics_persistence"]])
      
      phylogenomics_persistence <- matrix(nrow = 1,
                                          ncol = 13)
      culture_persistence <- matrix(nrow = 1,
                                    ncol = 13)
      overlapping_persistence <- matrix(nrow = length(subject_strain_index),
                                        ncol = 13)
      
      for (strain_entry in subject_strain_index) {
        overlapping_persistence[strain_entry,
                                persistence_list[[list_entry]][["strain_phylogenomics_persistence"]][[strain_entry]]$persistence] <- 1
        
        phylogenomics_persistence[1,persistence_list[[list_entry]][["strain_phylogenomics_persistence"]][[strain_entry]]$persistence] <- 1
        
        culture_persistence[1,
                            persistence_list[[list_entry]][["plot_df"]][persistence_list[[list_entry]][["plot_df"]]$category == "Culture","timepoint"]] <- 1
      }
      
      if (length(subject_strain_index) == 1) {
        # If there's only one strain for this subject
        persistence_timepoint <- persistence_list[[list_entry]][["strain_phylogenomics_persistence"]][[1]]$persistence
        
        # Directly modify persistence_sum
        persistence_sum$line_connect[persistence_sum$subject == subject &
                                       persistence_sum$category == "Phylogenomics" &
                                       persistence_sum$timepoint %in% persistence_timepoint] <- 1
      } else {
        # If there are multiple strains for this subject
        strains_group <- find_groups(overlapping_persistence)
        
        for (group_entry in seq_along(strains_group)) {
          group_persistence_timepoint <- as.integer(unique(unlist(sapply(strains_group[[group_entry]],
                                                                         function(strain_entry) {
                                                                           persistence_list[[list_entry]][["strain_phylogenomics_persistence"]][[strain_entry]]$persistence
                                                                         }))))
          # Directly modify persistence_sum
          persistence_sum$line_connect[persistence_sum$subject == subject &
                                         persistence_sum$category == "Phylogenomics" &
                                         persistence_sum$timepoint %in% group_persistence_timepoint] <- group_entry
        }
      }
      
      # Identify subjects where patterns differ from culture and phylogenomics
      
      if (!identical(phylogenomics_persistence, culture_persistence) |
          !identical(persistence_sum$line_connect[persistence_sum$subject == subject & persistence_sum$category == "Culture"],
                     persistence_sum$line_connect[persistence_sum$subject == subject & persistence_sum$category == "Phylogenomics"])) {
        
        modified_subject <- c(modified_subject, subject)
        
      }
    }
    
    
  }
  
  
  persistence_sum$subject <- factor(persistence_sum$subject,
                                    levels = sort(as.integer(unique(persistence_sum$subject))))
  
  get_persistence_plot <- function(persistence_sum,
                                   output_name,
                                   plot_height,
                                   plot_width){
    
    # Use persistence_sum to generate persistence summary plot
    category_colors <- c("Culture" = "#990000",
                         "Phylogenomics" = "#011F5B")
    
    persistence_plot_list <- lapply(sort(unique(persistence_sum$subject)),
                                    function(subject){
                                      
                                      subject_persistence <- persistence_sum[persistence_sum$subject == subject,]
                                      # background color, if even grey, If odd transparent
                                      background_color <- if(which(sort(unique(persistence_sum$subject)) == subject) %% 2 == 0) "lightgrey" else "transparent"
                                      
                                      if(subject == last(sort(unique(persistence_sum$subject)))){
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
                                      
                                      
                                      subject_persistence_plot <- ggplot(data = subject_persistence) + 
                                        geom_point(aes(x = timepoint,
                                                       y = interaction(category, subject),
                                                       color = category),
                                                   size = 5,
                                                   alpha = 1) +
                                        geom_line(data = subset(subject_persistence, category == "Phylogenomics" & line_connect != 0),
                                                  aes(x = timepoint, 
                                                      y = interaction(category, subject),
                                                      group = interaction(subject, line_connect),
                                                      color = category),
                                                  linewidth = 2.5,
                                                  alpha = 1) +
                                        geom_line(data = subset(subject_persistence, category == "Culture" & line_connect != 0),
                                                  aes(x = timepoint, 
                                                      y = interaction(category, subject),
                                                      group = interaction(subject, line_connect), 
                                                      color = category),
                                                  linewidth = 2.5,
                                                  alpha = 1)+
                                        scale_color_manual(values = category_colors,
                                                           guide = "none") + 
                                        scale_x_continuous(name = "Time Point",
                                                           breaks = seq(13),
                                                           limits = c(1,13)) +
                                        scale_y_discrete(name = subject) + 
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
                                    })
    
    persistence_plot <- reduce(persistence_plot_list, `/`) +
      plot_layout(heights = rep(1, length(persistence_plot_list))) +
      theme(plot.margin = margin(0, 0, 0, 0)
      )
    
    pdf(file=output_name,
        width=plot_width,
        height=plot_height)
    print(persistence_plot)
    dev.off()
  }
  
  # All subjects
  get_persistence_plot(persistence_sum = persistence_sum,
                       output_name = "All_Subjects_Persistence.pdf",
                       plot_width = 10,
                       plot_height = 50)
  saveRDS(persistence_sum,
          "persistence_sum.RDS")
  # Only subjects with persistence modified 
  modified_persistence_sum <- persistence_sum[persistence_sum$subject %in% modified_subject,]
  
  get_persistence_plot(persistence_sum = modified_persistence_sum,
                       output_name = "Modified_Subjects_Persistence.pdf",
                       plot_width = 10,
                       plot_height = 22.5)
  
  saveRDS(persistence_list,
          "persistence_list.RDS")
}


persistence_plot(metadata,
                 swabs_dates,
                 genomes2samples)
