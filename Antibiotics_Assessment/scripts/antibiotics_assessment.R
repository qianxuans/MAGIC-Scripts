# Libraries ----
required_packages <- c(
  "msm",
  "openxlsx2",
  "ggplot2",
  "parallel",
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

# Input ----
msm_input <- readRDS(file.path(input_dir,"hmm_abx_input.RDS"))
# Helper function to call msm ----
exec_abx_msm <- function(msm_input,
                         category){
  
  msm <- msm(pattern ~ days_after_birth,
             subject = subject,
             data = msm_input[[category]],
             cl = 0.95,
             covariates = ~ active_abx + others_abx,
             qmatrix = rbind(
               c(0.5,0.5),
               c(0.1,0.9)
             ),
             # underlying states as rows, and observed states as columns
             ematrix = rbind(
               c(0.8,0.2),
               c(0.2,0.8)
             ),
             censor = 3,
             censor.states = list(c(1,2)),
             control = list(fnscale=1000,
                            maxit=2000000000))
  
  # Transition probability matrix with antibiotics
  # The matrix of estimated transition probabilities 
  # Rows correspond to "from-state" and columns to "to-state".
  msm_both_abx_pmatrix <- pmatrix.msm(
    msm,
    covariates=list(active_abx = 1,
                    others_abx = 1),
    ci = "bootstrap",
    cl = 0.95,
    cores = 5,
    B = 100
  )
  
  msm_active_abx_pmatrix <- pmatrix.msm(
    msm,
    covariates=list(active_abx = 1,
                    others_abx = 0),
    ci = "bootstrap",
    cl = 0.95,
    cores = 5,
    B = 100
  )
  
  msm_others_abx_pmatrix <- pmatrix.msm(
    msm,
    covariates=list(active_abx = 0,
                    others_abx = 1),
    ci = "bootstrap",
    cl = 0.95,
    cores = 5,
    B = 100
  )
  
  msm_no_abx_pmatrix <- pmatrix.msm(
    msm,
    covariates=list(active_abx = 0,
                    others_abx = 0),
    ci = "bootstrap",
    cl = 0.95,
    cores = 5,
    B = 100
  )
  
  
  # Total length of stay
  # Mean sojourn times describe the average period in a single stay in a state
  
  msm_both_abx_length <- totlos.msm(msm,
                                    start = 1,
                                    end = NULL,
                                    fromt = 0,
                                    tot = 720,
                                    covariates=list(active_abx = 1,
                                                    others_abx = 1),
                                    ci = "bootstrap",
                                    cl = 0.95,
                                    cores = 5,
                                    B = 100)
  
  msm_active_abx_length <- totlos.msm(msm,
                                      start = 1,
                                      end = NULL,
                                      fromt = 0,
                                      tot = 720,
                                      covariates=list(active_abx = 1,
                                                      others_abx = 0),
                                      ci = "bootstrap",
                                      cl = 0.95,
                                      cores = 5,
                                      B = 100)
  
  msm_others_abx_length <- totlos.msm(msm,
                                      start = 1,
                                      end = NULL,
                                      fromt = 0,
                                      tot = 720,
                                      covariates=list(active_abx = 0,
                                                      others_abx = 1),
                                      ci = "bootstrap",
                                      cl = 0.95,
                                      cores = 5,
                                      B = 100)
  
  msm_no_abx_length <- totlos.msm(msm,
                                  start = 1,
                                  end = NULL,
                                  fromt = 0,
                                  tot = 720,
                                  covariates=list(active_abx = 0,
                                                  others_abx = 0),
                                  ci = "bootstrap",
                                  cl = 0.95,
                                  cores = 5,
                                  B = 100)
  
  
  # return the results
  list(
    full = list(
      category = category,
      msm = msm,
      # Transition Probability
      msm_both_abx_pmatrix = msm_both_abx_pmatrix,
      msm_active_abx_pmatrix = msm_active_abx_pmatrix,
      msm_others_abx_pmatrix = msm_others_abx_pmatrix,
      msm_no_abx_pmatrix = msm_no_abx_pmatrix,
      # Stay of Length
      msm_both_abx_length = msm_both_abx_length,
      msm_active_abx_length = msm_active_abx_length,
      msm_others_abx_length = msm_others_abx_length,
      msm_no_abx_length = msm_no_abx_length
    )
  )
}

## Get abx msm results ----
if(file.exists(file.path(output_dir,"hmm_abx_results.RDS"))){
  abx_msm_results <- readRDS(file.path(output_dir,"hmm_abx_results.RDS"))
}else{
  abx_msm_results <- lapply(c("culture_abx",
                              "phylogenomics_abx"),
                            function(category) {
                              exec_abx_msm(hmm_abx_input,
                                           category)
                            })
  
  saveRDS(abx_msm_results, file = "hmm_abx_results.RDS")
}

# Visualize the results ----

# Helper function to get the data frame for ggplot2 ----

get_plot_df_from_msm_results <- function(msm_result_list,
                                         abx_status,
                                         category){
  
  # Division:
  
  #abx_status is one of the following:
  # 1. both_abx
  # 2. active_abx
  # 3. others_abx
  # 4. no_abx
  
  # category is one of the following:
  # 1. pmatrix
  # 2. length
  
  entry <- paste0("msm_",
                  abx_status,
                  "_",
                  category)
  
  entry_index <- which(grepl(entry,names(msm_result_list$full)))
  
  if(category == "pmatrix"){
    result_matrix <- msm_result_list$full[[entry_index]]
    df <- data.frame(from = c("Positive","Positive","Negative","Negative"),
                     to = c("Positive","Negative","Positive","Negative"),
                     # Rows correspond to "from-state" and columns to "to-state".
                     probability = c(result_matrix$estimates[4]*100,
                                     result_matrix$estimates[2]*100,
                                     result_matrix$estimates[3]*100,
                                     result_matrix$estimates[1]*100),
                     low = c(result_matrix$L[4]*100,
                             result_matrix$L[2]*100,
                             result_matrix$L[3]*100,
                             result_matrix$L[1]*100),
                     up = c(result_matrix$U[4]*100,
                            result_matrix$U[2]*100,
                            result_matrix$U[3]*100,
                            result_matrix$U[1]*100))
    
    return(cbind(category = msm_result_list$full$category,
                 abx_status = abx_status,
                 result = paste0("From ",df$from,
                                 " To ",df$to),
                 df))
    
  }else if(category == "length"){
    
    result_matrix <- msm_result_list$full[[entry_index]]
    
    df <- data.frame(status = c("Positive","Negative"),
                     length = c(result_matrix[1,2],result_matrix[1,1]),
                     low = c(result_matrix[2,2],result_matrix[2,1]),
                     up = c(result_matrix[3,2],result_matrix[3,1]))
    return(
      cbind(category = msm_result_list$full$category,
            abx_status = abx_status,
            df)
    )
  }
}
# Function to make the plots ---- 
visualize_probability <- function(abx_msm_results){
  
  
  # culture transition
  culture_transition_df <- rbind(get_plot_df_from_msm_results(abx_msm_results[[1]],
                                                              abx_status = "no_abx",
                                                              category = "pmatrix"),
                                 get_plot_df_from_msm_results(abx_msm_results[[1]],
                                                              abx_status = "active_abx",
                                                              category = "pmatrix"),
                                 get_plot_df_from_msm_results(abx_msm_results[[1]],
                                                              abx_status = "others_abx",
                                                              category = "pmatrix"),
                                 get_plot_df_from_msm_results(abx_msm_results[[1]],
                                                              abx_status = "both_abx",
                                                              category = "pmatrix"))
  
  # phylogenomics transition
  phylogenomics_transition_df <- rbind(get_plot_df_from_msm_results(abx_msm_results[[2]],
                                                                    abx_status = "no_abx",
                                                                    category = "pmatrix"),
                                       get_plot_df_from_msm_results(abx_msm_results[[2]],
                                                                    abx_status = "active_abx",
                                                                    category = "pmatrix"),
                                       get_plot_df_from_msm_results(abx_msm_results[[2]],
                                                                    abx_status = "others_abx",
                                                                    category = "pmatrix"),
                                       get_plot_df_from_msm_results(abx_msm_results[[2]],
                                                                    abx_status = "both_abx",
                                                                    category = "pmatrix"))
  
  # culture Length
  culture_length_df <- rbind(get_plot_df_from_msm_results(abx_msm_results[[1]],
                                                          abx_status = "no_abx",
                                                          category = "length"),
                             get_plot_df_from_msm_results(abx_msm_results[[1]],
                                                          abx_status = "active_abx",
                                                          category = "length"),
                             get_plot_df_from_msm_results(abx_msm_results[[1]],
                                                          abx_status = "others_abx",
                                                          category = "length"),
                             get_plot_df_from_msm_results(abx_msm_results[[1]],
                                                          abx_status = "both_abx",
                                                          category = "length"))
  
  
  # phylogenomics Length
  phylogenomics_length_df <- rbind(get_plot_df_from_msm_results(abx_msm_results[[2]],
                                                                abx_status = "no_abx",
                                                                category = "length"),
                                   get_plot_df_from_msm_results(abx_msm_results[[2]],
                                                                abx_status = "active_abx",
                                                                category = "length"),
                                   get_plot_df_from_msm_results(abx_msm_results[[2]],
                                                                abx_status = "others_abx",
                                                                category = "length"),
                                   get_plot_df_from_msm_results(abx_msm_results[[2]],
                                                                abx_status = "both_abx",
                                                                category = "length"))
  
  # Category factor
  
  culture_transition_df$abx_status <- factor(
    dplyr::recode(
      culture_transition_df$abx_status,
      "no_abx" = "No Antibiotics",
      "active_abx" = "Active Antibiotics",
      "others_abx" = "Others Antibiotics",
      "both_abx" = "All Antibiotics"
    ),
    levels = c("No Antibiotics", "Active Antibiotics", "Others Antibiotics", "All Antibiotics")
  )
  
  phylogenomics_transition_df$abx_status <- factor(
    dplyr::recode(
      phylogenomics_transition_df$abx_status,
      "no_abx" = "No Antibiotics",
      "active_abx" = "Active Antibiotics",
      "others_abx" = "Others Antibiotics",
      "both_abx" = "All Antibiotics"
    ),
    levels = c("No Antibiotics", "Active Antibiotics", "Others Antibiotics", "All Antibiotics")
  )
  
  phylogenomics_length_df$abx_status <- factor(
    dplyr::recode(
      phylogenomics_length_df$abx_status,
      "no_abx" = "No Antibiotics",
      "active_abx" = "Active Antibiotics",
      "others_abx" = "Others Antibiotics",
      "both_abx" = "All Antibiotics"
    ),
    levels = c("No Antibiotics", "Active Antibiotics", "Others Antibiotics", "All Antibiotics")
  )
  
  culture_length_df$abx_status <- factor(
    dplyr::recode(
      culture_length_df$abx_status,
      "no_abx" = "No Antibiotics",
      "active_abx" = "Active Antibiotics",
      "others_abx" = "Others Antibiotics",
      "both_abx" = "All Antibiotics"
    ),
    levels = c("No Antibiotics", "Active Antibiotics", "Others Antibiotics", "All Antibiotics")
  )
  
  
  # color schemes
  bar_colors <- c("white",
                  "#B2C06E",
                  "#D779A9",
                  "#5A769D")
  
  names(bar_colors) <- c("No Antibiotics",
                         "Active Antibiotics",
                         "Others Antibioticss",
                         "All Antibiotics")
  
  text_colors <- c("black",
                   "#829726",
                   "#C72374",
                   "#315485")
  
  names(text_colors) <- c("No Antibiotics",
                          "Active Antibiotics",
                          "Others Antibioticss",
                          "All Antibiotics")
  
  ## Transition Plot ----
  culture_transition_df <- culture_transition_df[culture_transition_df$result %in% c("From Negative To Positive",
                                                                                     "From Positive To Negative"),]
  # Only Show from "From Positive to Negative" and "From Negative to Positive"
  culture_transition_plot <- ggplot(culture_transition_df,
                                    aes(x = abx_status,
                                        y = probability,
                                        fill = abx_status,
                                        color = abx_status)) +
    
    geom_bar(stat = "identity",
             position = position_dodge(),
             alpha = 0.5) +
    geom_text(aes(label = round(probability,
                                digits = 2),
                  color = abx_status),
              position = position_stack(vjust= 0.5),
              size = 7.5) + 
    scale_color_manual(values = text_colors) +
    scale_fill_manual(values = bar_colors) +
    geom_errorbar(aes(ymin = low, ymax = up),
                  width = 0.25,
                  position = position_dodge(0.9)) +
    facet_wrap(~ result,
               scales = "free_y",
               ncol = 2) +
    theme_minimal() +
    labs(title = "Culture Transition Probability",
         x = "Antibiotics Status", 
         y = "%Probability") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 25),
          strip.text = element_text(face = "bold",
                                    size = 20),
          axis.text.x = element_text(face = "bold",
                                     size = 15,
                                     angle = 45,
                                     hjust = 1,
                                     vjust = 1),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold",
                                     size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 25),
          axis.ticks.y = element_line(colour = "black",
                                      linewidth = 1),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 2))
  
  pdf(file= "Culture_Transition_Probability_Active_Others.pdf",
      width=12.5,
      height=5)
  print(culture_transition_plot)
  dev.off()
  
  phylogenomics_transition_df <- phylogenomics_transition_df[phylogenomics_transition_df$result %in% c("From Negative To Positive",
                                                                                                       "From Positive To Negative"),]
  phylogenomics_transition_plot <- ggplot(phylogenomics_transition_df,
                                          aes(x = abx_status,
                                              y = probability,
                                              fill = abx_status,
                                              color = abx_status)) +
    
    geom_bar(stat = "identity",
             position = position_dodge(),
             alpha = 0.5) +
    geom_text(aes(label = round(probability,
                                digits = 2),
                  color = abx_status),
              position = position_stack(vjust= 0.5),
              size = 7.5) + 
    scale_color_manual(values = text_colors) +
    scale_fill_manual(values = bar_colors) +
    geom_errorbar(aes(ymin = low, ymax = up),
                  width = 0.25,
                  position = position_dodge(0.9)) +
    facet_wrap(~ result,
               scales = "free_y",
               ncol = 2) +
    theme_minimal() +
    labs(title = "Phylogenomics Transition Probability",
         x = "Antibiotics Status", 
         y = "%Probability") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 25),
          strip.text = element_text(face = "bold",
                                    size = 20),
          axis.text.x = element_text(face = "bold",
                                     size = 15,
                                     angle = 45,
                                     hjust = 1,
                                     vjust = 1),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold",
                                     size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 25),
          axis.ticks.y = element_line(colour = "black",
                                      linewidth = 1),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 2))
  
  pdf(file= "Phylogenomics_Transition_Probability_Active_Others.pdf",
      width=12.5,
      height=5)
  print(phylogenomics_transition_plot)
  dev.off()
  
  
  
  culture_length_df <- culture_length_df[culture_length_df$status == "Positive",]
  culture_length_plot <- ggplot(culture_length_df,
                                aes(x = abx_status,
                                    y = length,
                                    fill = abx_status,
                                    color = abx_status)) +
    
    geom_bar(stat = "identity",
             position = position_dodge(),
             alpha = 0.5) +
    geom_text(aes(label = round(length,
                                digits = 2),
                  color = abx_status),
              position = position_stack(vjust= 0.5),
              size = 7.5) + 
    scale_color_manual(values = text_colors) +
    scale_fill_manual(values = bar_colors) +
    geom_errorbar(aes(ymin = low, ymax = up),
                  width = 0.25,
                  position = position_dodge(0.9)) +
    facet_wrap(~ status,
               scales = "free_y",
               ncol = 2) +
    theme_minimal() +
    labs(title = "Culture Persistence",
         x = "Antibiotics Status", 
         y = "Days") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 25),
          strip.text = element_blank(),
          axis.text.x = element_text(face = "bold",
                                     size = 15,
                                     angle = 45,
                                     hjust = 1,
                                     vjust = 1),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold",
                                     size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 25),
          axis.ticks.y = element_line(colour = "black",
                                      linewidth = 1),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 2))
  
  pdf(file= "Culture_Persistence_Active_Others.pdf",
      width=5,
      height=5)
  print(culture_length_plot)
  dev.off()
  
  phylogenomics_length_df <- phylogenomics_length_df[phylogenomics_length_df$status == "Positive",]
  
  phylogenomics_length_plot <- ggplot(phylogenomics_length_df,
                                      aes(x = abx_status,
                                          y = length,
                                          fill = abx_status,
                                          color = abx_status)) +
    
    geom_bar(stat = "identity",
             position = position_dodge(),
             alpha = 0.5) +
    geom_text(aes(label = round(length,
                                digits = 2),
                  color = abx_status),
              position = position_stack(vjust= 0.5),
              size = 7.5) + 
    scale_color_manual(values = text_colors) +
    scale_fill_manual(values = bar_colors) +
    geom_errorbar(aes(ymin = low, ymax = up),
                  width = 0.25,
                  position = position_dodge(0.9)) +
    facet_wrap(~ status,
               scales = "free_y",
               ncol = 2) +
    theme_minimal() +
    labs(title = "Phylogenomics Persistence",
         x = "Antibiotics Status", 
         y = "Days") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 25),
          strip.text = element_blank(),
          axis.text.x = element_text(face = "bold",
                                     size = 15,
                                     angle = 45,
                                     hjust = 1,
                                     vjust = 1),
          #axis.text.x = element_blank(),
          axis.text.y = element_text(face = "bold",
                                     size = 20),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face = "bold",
                                      size = 25),
          axis.ticks.y = element_line(colour = "black",
                                      linewidth = 1),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 2))
  
  pdf(file= "Phylogenomics_Persistence_Active_Others.pdf",
      width=5,
      height=5)
  print(phylogenomics_length_plot)
  dev.off()
  
  
}
# Execute the function ----
visualize_probability(abx_msm_results)