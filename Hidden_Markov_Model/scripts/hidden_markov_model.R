# Libraries ----
required_packages <- c(
  "msm",
  "dplyr",
  "openxlsx2",
  "ggplot2",
  "ggbreak"
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

msm_input <- readRDS(file.path(input_dir,"msm_input.RDS"))

# Analysis Functions ----
exec_msm <- function(msm_input,
                     category){
  # Use msm to perform the HMM 
  
  # A model is governed by a transition intensity matrix
  
  # along with a qmatrix whose non-zero entries represent the allowed transitions of the model
  
  # Fitting misclassification models with msm ----
  # To model observed states with misclassification, we define a matrix ematrix indicating the
  # states that can be misclassified.
  # Rows of this matrix correspond to true states, columns to observed states.
  
  # multi-state model with misclassification
  #set the seed
  set.seed(12050)
  
  msm <- msm(pattern ~ days_after_birth,
             subject = subject,
             data = msm_input[[category]],
             cl = 0.95,
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
  # Misclassification probability
  # Estimated misclassification probability matrix.
  # The rows correspond to true states, and columns observed states.
  msm_ematrix <- ematrix.msm(
    msm,
    covariates = 0,
    ci = "bootstrap",
    cl = 0.95,
    B = 100,
    cores = 8
  )
  
  
  # Total length of stay
  # Mean sojourn times describe the average period in a single stay in a state
  
  msm_length <- totlos.msm(msm,
                           start = 1,
                           end = NULL,
                           fromt = 0,
                           tot = 720,
                           cl = 0.95,
                           ci = "bootstrap",
                           B = 100,
                           cores = 8)
  
  # Transition probability matrix
  # The matrix of estimated transition probabilities 
  # Rows correspond to "from-state" and columns to "to-state".
  
  msm_pmatrix <- pmatrix.msm(
    msm,
    covariates = 0,
    ci = "bootstrap",
    cl = 0.95,
    cores = 8,
    B = 100
  )
  
  # return the results
  list(
    full = list(
      category = category,
      msm = msm,
      msm_pmatrix = msm_pmatrix,
      msm_ematrix = msm_ematrix,
      msm_length = msm_length
    ),
    transition_p = data.frame(from = c("Positive","Positive","Negative","Negative"),
                              to = c("Positive","Negative","Positive","Negative"),
                              # Rows correspond to "from-state" and columns to "to-state".
                              probability = c(msm_pmatrix$estimates[4]*100,
                                              msm_pmatrix$estimates[2]*100,
                                              msm_pmatrix$estimates[3]*100,
                                              msm_pmatrix$estimates[1]*100),
                              low = c(msm_pmatrix$L[4]*100,
                                      msm_pmatrix$L[2]*100,
                                      msm_pmatrix$L[3]*100,
                                      msm_pmatrix$L[1]*100),
                              up = c(msm_pmatrix$U[4]*100,
                                     msm_pmatrix$U[2]*100,
                                     msm_pmatrix$U[3]*100,
                                     msm_pmatrix$U[1]*100)),
    misclassification_p = data.frame(observed = c("Positive","Positive","Negative","Negative"),
                                     true = c("Positive","Negative","Positive","Negative"),
                                     result = c("True Positive","False Positive","False Negative","True Negative"),
                                     # The rows correspond to true states, and columns observed states.
                                     probability = c(msm_ematrix$estimates[4]*100,
                                                     msm_ematrix$estimates[3]*100,
                                                     msm_ematrix$estimates[2]*100,
                                                     msm_ematrix$estimates[1]*100),
                                     low = c(msm_ematrix$L[4]*100,
                                             msm_ematrix$L[3]*100,
                                             msm_ematrix$L[2]*100,
                                             msm_ematrix$L[1]*100),
                                     up = c(msm_ematrix$U[4]*100,
                                            msm_ematrix$U[3]*100,
                                            msm_ematrix$U[2]*100,
                                            msm_ematrix$U[1]*100)),
    stay_length = data.frame(status = c("Positive","Negative"),
                             length = c(msm_length[1,2],msm_length[1,1]),
                             low = c(msm_length[2,2],msm_length[2,1]),
                             up = c(msm_length[3,2],msm_length[3,1]))
  )
}


sort_emission_df <- function(msm_results){
  return(do.call(rbind,
                 lapply(msm_results,
                        function(result){
                          return(
                            cbind(category = result$full$category,
                                  result$misclassification_p)
                          )
                        })))
}
# Transition Probability
sort_transition_df <- function(msm_results){
  return(do.call(rbind,
                 lapply(msm_results,
                        function(result){
                          return(
                            cbind(category = result$full$category,
                                  result = paste0("From ",result$transition_p$from,
                                                  " To ",result$transition_p$to),
                                  result$transition_p)
                          )
                        })))
}

# Stay Length
sort_stay_length <- function(msm_results){
  return(do.call(rbind,
                 lapply(msm_results,
                        function(result){
                          return(
                            cbind(category = result$full$category,
                                  result$stay_length)
                          )
                        })))
}



# Helper to visualize the results 
visualize_probability <- function(total_transition_df,
                                  total_emission_df,
                                  total_stay_length){
  #misclassification 
  
  misclassification_df <- total_emission_df[total_emission_df$result %in% c("False Negative",
                                                                            "False Positive"),]
  
  misclassification_df$category <- gsub("phylogenomics", "Phylogenomics", 
                                        gsub("culture", "Culture",
                                             gsub("corrected_culture", "Corrected Culture",misclassification_df$category)))
  
  
  misclassification_df$category <- factor(misclassification_df$category, 
                                          levels = c("Culture","Corrected Culture","Phylogenomics"))
  
  #transition
  #use phylogenomics correct as phylogenomics
  total_transition_df <- total_transition_df[total_transition_df$result %in% c("From Negative To Positive",
                                                                               "From Positive To Negative"),]
  total_transition_df$category <- gsub("phylogenomics", "Phylogenomics", 
                                       gsub("culture", "Culture",
                                            gsub("corrected_culture", "Corrected Culture",total_transition_df$category)))
  
  total_transition_df$category <- factor(total_transition_df$category, 
                                         levels = c("Culture","Corrected Culture","Phylogenomics"))
  #stay length
  total_stay_length <- total_stay_length[total_stay_length$status == "Positive",]
  total_stay_length$category <- gsub("phylogenomics", "Phylogenomics", 
                                     gsub("culture", "Culture",
                                          gsub("corrected_culture", "Corrected Culture",total_stay_length$category)))
  
  total_stay_length$category <- factor(total_stay_length$category, 
                                       levels = c("Culture","Corrected Culture","Phylogenomics"))
  # color schemes
  bar_colors <- c("#990000",
                  "#7600bc",
                  "#011F5B")
  
  names(bar_colors) <- c("Culture",
                         "Corrected Culture",
                         "Phylogenomics")
  
  text_colors <- c("#420D09",
                   "#4c00b0",
                   "#29465B")
  
  names(text_colors) <- c("Culture",
                          "Corrected Culture",
                          "Phylogenomics")
  
  
  transition_p_plot <- ggplot(total_transition_df,
                              aes(x = category,
                                  y = probability,
                                  fill = category,
                                  color = category)) +
    
    geom_bar(stat = "identity",
             position = position_dodge(),
             alpha = 0.5) +
    geom_text(aes(label = round(probability,
                                digits = 2),
                  color = category),
              position = position_stack(vjust= 0.5),
              size = 7.5) + 
    scale_color_manual(values = bar_colors) +
    scale_fill_manual(values = text_colors) +
    geom_errorbar(aes(ymin = low, ymax = up),
                  width = 0.25,
                  position = position_dodge(0.9)) +
    facet_wrap(~ result,
               scales = "fixed",
               ncol = 2) +
    scale_y_break(c(20, 40),
                  scales = 0.25,
                  space = 0.075) + 
    theme_minimal() +
    labs(title = "Transition Probability",
         x = "Category", 
         y = "%Probability") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 25),
          strip.text = element_text(face = "bold",
                                    size = 20),
          axis.text.x = element_text(face = "bold",
                                     size = 20,
                                     angle = 45,
                                     hjust = 1),
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
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  pdf(file="Transition_Probability.pdf",
      width=10,
      height=5.5)
  print(transition_p_plot)
  dev.off()
  
  # misclassification
  misclassification_p_plot <- ggplot(misclassification_df,
                                     aes(x = category,
                                         y = probability,
                                         fill = category,
                                         color = category)) +
    
    geom_bar(stat = "identity",
             position = position_dodge(),
             alpha = 0.5) +
    geom_text(aes(label = round(probability,
                                digits = 2),
                  color = category),
              position = position_stack(vjust= 0.5),
              size = 7.5) + 
    scale_color_manual(values = bar_colors) +
    scale_fill_manual(values = text_colors) +
    geom_errorbar(aes(ymin = low, ymax = up),
                  width = 0.25,
                  position = position_dodge(0.9)) +
    facet_wrap(~ result,
               scales = "fixed",
               ncol = 2) +
    scale_y_break(c(30, 70),
                  scales = 0.25,
                  space = 0.075,
    ) + 
    theme_minimal() +
    labs(title = "Misclassification Probability",
         x = "Category", 
         y = "%Probability") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 25),
          strip.text = element_text(face = "bold",
                                    size = 17.5),
          axis.text.x = element_text(face = "bold",
                                     size = 20,
                                     angle = 45,
                                     hjust = 1),
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
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  pdf(file="Misclassification_Probability.pdf",
      width=10,
      height=5.5)
  print(misclassification_p_plot)
  dev.off()
  
  # stay length
  stay_length_plot <- ggplot(total_stay_length,
                             aes(x = category,
                                 y = length,
                                 fill = category,
                                 color = category)) +
    
    geom_bar(stat = "identity",
             position = position_dodge(),
             alpha = 0.5) +
    geom_text(aes(label = round(length,
                                digits = 2),
                  color = category),
              position = position_stack(vjust= 0.5),
              size = 7.5) + 
    scale_color_manual(values = bar_colors) +
    scale_fill_manual(values = text_colors) +
    geom_errorbar(aes(ymin = low, ymax = up),
                  width = 0.25,
                  position = position_dodge(0.9)) +
    facet_wrap(~ status,
               scales = "fixed",
               ncol = 1) +
    theme_minimal() +
    labs(title = "Persistence",
         x = "Category", 
         y = "Days") +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold",
                                    size = 25),
          strip.text = element_blank(),
          axis.text.x = element_text(face = "bold",
                                     size = 20,
                                     angle = 45,
                                     hjust = 1),
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
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  
  pdf(file=paste0("Length_Stay.pdf"),
      width=5,
      height=5.5)
  print(stay_length_plot)
  dev.off()
  
  
}


# Execute Analysis ----

if(file.exists(file.path(output_dir,"msm_results.RDS"))){
  msm_results <- readRDS(file.path(output_dir,"msm_results.RDS"))
}else{
  msm_results <- lapply(c("culture",
                          "corrected_culture",
                          "phylogenomics"),
                        function(category) {
                          exec_msm(msm_input,
                                   category)
                        })
  
  saveRDS(msm_results,
          file = "msm_results.RDS")
}

#misclassification 
total_emission_df <- sort_emission_df(msm_results)

#transition
total_transition_df <- sort_transition_df(msm_results)

# stay length
total_stay_length <- sort_stay_length(msm_results)

# visualize the results
visualize_probability(total_transition_df,
                      total_emission_df,
                      total_stay_length)