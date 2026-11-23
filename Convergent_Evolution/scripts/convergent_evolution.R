# Libraries ----
required_packages <- c(
  "dplyr",
  "openxlsx2",
  "ggplot2",
  "parallel"
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

# Settings ----
mutation_sets <- c("protein_altering", "synonymous", "all")
synonymous_terms <- c("synonymous_variant", "stop_retained_variant", "start_retained_variant")
min_lineages_parallel <- 2
n_null_simulations <- 100000
n_cores <- max(1, detectCores() - 4)

go_categories <- c("Molecular Function",
                   "Cellular Component",
                   "Biological Process")

# count each OG toward the parents of its annotated terms
propagate_go_ancestors <- TRUE 
# OGs per term. a single-OG term repeats the OG-level test
min_go_term_size <- 2
# drops near-root terms that span most of the genome
max_go_term_size <- 500
# reproducible random streams across mclapply workers
RNGkind("L'Ecuyer-CMRG")
set.seed(202611)

# Input ----
## Assembly-Scan results ----
assembly_scan_res <- readRDS(file.path(input_dir,"magic_assembly_scan_res.RDS"))
## Variants called by snippy ----
variant_df <- readRDS(file.path(input_dir,"strain_snps.RDS"))
subject_strain_reference <- readRDS(file.path(input_dir,"all_subject_strain_reference.RDS"))
## THRESHER strains ----
strains <- readRDS(file.path(input_dir,"MAGIC_thresher_plateau.RDS"))[["strains"]]
strains$lineage <- paste(strains$subject,strains$new_strain_id,sep = "|")
## Panaroo ortholog groups and gene sequences ----

gene_ortholog_group <- read.csv(file.path(input_dir,"gene_presence_absence.csv"),
                                check.names = FALSE)
gene_lengths_file <- file.path(input_dir,"gene_lengths.RDS")
## GO terms per ortholog group (from the GO-mapping script) ----
# columns used: ortholog_name, go_term, go_name, category
ortholog2go <- readRDS(file.path(input_dir,"ortholog2go.RDS"))
rownames(ortholog2go) <- NULL
## Results RDS path ----
# New file name for the whole-genome null, so results from the coding-length version aren't reloaded
all_results_list_file <- file.path(input_dir,"all_results_list_wg.RDS")

# Analysis ----

# Locus tag to ortholog group ----
genome_columns <- 4:ncol(gene_ortholog_group)

og_cells <- data.frame(
  ortholog_group = rep(gene_ortholog_group$Gene, times = length(genome_columns)),
  gene_ids = unlist(lapply(gene_ortholog_group[genome_columns], as.character), use.names = FALSE)
)
og_cells <- og_cells[!is.na(og_cells$gene_ids) & og_cells$gene_ids != "", ]

gene_id_list <- strsplit(og_cells$gene_ids, ";", fixed = TRUE)

# Paralogs are ';'-separated within a cell. Panaroo adds flag suffixes to some
# IDs (e.g. '_pseudo') and gives refound genes IDs that are not locus tags, so
locus_tag_og <- data.frame(
  ortholog_group = rep(og_cells$ortholog_group, lengths(gene_id_list)),
  gene_id = trimws(unlist(gene_id_list, use.names = FALSE))
)

# only the '<genome>_<5 digits>' locus tag is kept. 
locus_tag_og <- locus_tag_og[grepl("_\\d{5}", locus_tag_og$gene_id), ]
locus_tag_og$locus_tag <- sub("^(.+_\\d{5}).*$", "\\1", locus_tag_og$gene_id)
locus_tag_og$genome <- sub("_\\d{5}$", "", locus_tag_og$locus_tag)
locus_tag_og <- distinct(locus_tag_og, locus_tag, .keep_all = TRUE)

# Lineages and within-host mutations ----
# Use the strain in subject_strain_reference 
# This data frame records the strain with more than one isolates and 
# the isolates were not collected at the same date

subject_strain_reference$lineage  <- paste(subject_strain_reference$subject,
                                           subject_strain_reference$strain,
                                           sep = "|")

genic_variants <- variant_df %>%
  mutate(subject = as.character(subject),
         strain = as.character(strain),
         lineage = paste(subject, strain, sep = "|"),
         effect_term = sub(" .*$", "", EFFECT),
         # snpEff labels variants outside CDS (intergenic, RNA genes) with terms such as
         # intergenic_region or intragenic_variant, so classify by snippy's feature type
         mutation_class = case_when(
           is.na(FTYPE) | FTYPE != "CDS" ~ "non_coding",
           effect_term %in% synonymous_terms ~ "synonymous",
           TRUE ~ "protein_altering")) %>%
  filter(lineage %in% subject_strain_reference$lineage) %>%
  # a mutation is called in every later isolate that carries it; count it once
  distinct(lineage, CHROM, POS, REF, ALT, .keep_all = TRUE) 


# Create the og_mutations for downstream analysis 
og_mutations <- genic_variants

# Add the ortholog group to og_mutations

og_mutations$ortholog_group <- sapply(og_mutations$LOCUS_TAG,function(tag_entry){
  if(tag_entry != "" && !is.na(tag_entry)){
    group_entry <- locus_tag_og$ortholog_group[locus_tag_og$locus_tag == tag_entry]
    if(length(group_entry) == 0){
      NA_character_
    }else{
      group_entry
    }
  }else{
    NA_character_
  }
})

og_mutations$reference_genome <- sapply(og_mutations$lineage,function(lineage_entry){
  subject_strain_reference$reference[subject_strain_reference$lineage == lineage_entry]
})

# Coding target size and whole-genome size of each reference genome ----
reference_genomes <- unique(subject_strain_reference$reference)
reference_locus_tags <- locus_tag_og$locus_tag[locus_tag_og$genome %in% reference_genomes]
gene_lengths <- readRDS(gene_lengths_file)


reference_og_wg_length <- locus_tag_og %>%
  filter(genome %in% reference_genomes) %>%
  inner_join(gene_lengths, by = "locus_tag") %>%
  group_by(reference_genome = genome, ortholog_group) %>%
  summarise(og_length = sum(gene_length), n_copies = n(), .groups = "drop") %>%
  mutate(coding_length = sum(og_length), .by = reference_genome,
         wg_size = as.integer(assembly_scan_res$total_contig_length[match(reference_genome, assembly_scan_res$genome)]))

# 4. GO terms per ortholog group ----
og_map <- data.frame(target = unique(locus_tag_og$ortholog_group))
og_map$ortholog_group <- og_map$target

use_go_db <- requireNamespace("GO.db", quietly = TRUE)
go_category_labels <- c(MF = "Molecular Function", CC = "Cellular Component", BP = "Biological Process")

lookup_go_db <- function(go_terms) {
  valid_terms <- intersect(unique(go_terms), AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"))
  go_info <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = valid_terms,
                                                    columns = c("TERM", "ONTOLOGY"), keytype = "GOID"))
  data.frame(go_term = go_info$GOID,
             go_db_name = go_info$TERM,
             go_db_category = unname(go_category_labels[go_info$ONTOLOGY]))
}

og_go_direct <- ortholog2go %>%
  transmute(ortholog_group = as.character(ortholog_name),
            go_term = trimws(as.character(go_term)),
            go_name = as.character(go_name),
            category = as.character(category)) %>%
  filter(grepl("^GO:\\d{7}$", go_term), ortholog_group %in% og_map$ortholog_group) %>%
  distinct(ortholog_group, go_term, .keep_all = TRUE)

n_go_annotated_ogs <- n_distinct(og_go_direct$ortholog_group)
message(n_go_annotated_ogs, " of ", nrow(og_map), " ortholog groups have >= 1 GO term")

# terms missing from the QuickGO download: name and category from GO.db
if (use_go_db) {
  og_go_direct <- og_go_direct %>%
    left_join(lookup_go_db(og_go_direct$go_term), by = "go_term") %>%
    mutate(go_name = coalesce(go_name, go_db_name),
           category = coalesce(category, go_db_category)) %>%
    select(-go_db_name, -go_db_category)
}
n_terms_without_category <- n_distinct(og_go_direct$go_term[!og_go_direct$category %in% go_categories])
if (n_terms_without_category > 0) {
  message(n_terms_without_category, " GO term(s) without a category (obsolete or not found) are dropped")
}
og_go_direct <- filter(og_go_direct, category %in% go_categories)

og_go <- select(og_go_direct, ortholog_group, go_term, category)
if (propagate_go_ancestors) {
  if (!use_go_db) stop("propagate_go_ancestors = TRUE needs the GO.db package")
  ancestor_maps <- list("Biological Process" = GO.db::GOBPANCESTOR,
                        "Molecular Function" = GO.db::GOMFANCESTOR,
                        "Cellular Component" = GO.db::GOCCANCESTOR)
  og_go_ancestors <- bind_rows(lapply(go_categories, function(go_category) {
    category_terms <- filter(og_go, category == go_category)
    if (nrow(category_terms) == 0) return(NULL)
    # as.list() over the whole map instead of mget(): AnnotationDbi's mget()
    # method fails under some R / AnnotationDbi version combinations
    all_ancestors <- AnnotationDbi::as.list(ancestor_maps[[go_category]])
    ancestors <- all_ancestors[intersect(unique(category_terms$go_term), names(all_ancestors))]
    if (length(ancestors) == 0) return(NULL)
    term_ancestors <- data.frame(go_term = rep(names(ancestors), lengths(ancestors)),
                                 ancestor = unlist(ancestors, use.names = FALSE))
    category_terms %>%
      inner_join(term_ancestors, by = "go_term", relationship = "many-to-many") %>%
      filter(!is.na(ancestor), ancestor != "all") %>%
      transmute(ortholog_group, go_term = ancestor, category)
  }))
  og_go <- distinct(bind_rows(og_go, og_go_ancestors))
}

go_term_names <- distinct(filter(og_go_direct, !is.na(go_name)), go_term, go_name)
if (use_go_db) {
  go_term_names <- full_join(go_term_names, lookup_go_db(og_go$go_term), by = "go_term") %>%
    transmute(go_term, go_name = coalesce(go_name, go_db_name))
}
go_term_names <- distinct(go_term_names, go_term, .keep_all = TRUE)

go_term_size <- count(og_go, category, go_term, name = "term_size")
og_go_tested <- og_go %>%
  inner_join(go_term_size, by = c("category", "go_term")) %>%
  filter(term_size >= min_go_term_size, term_size <= max_go_term_size)
message(n_distinct(og_go_tested$go_term), " GO terms with ", min_go_term_size, "-",
        max_go_term_size, " ortholog groups are tested")

# Test functions ----
# P(X >= observed) for X = number of lineages hit, exact by convolution
poisson_binomial_upper_tail <- function(observed, hit_probabilities) {
  if (observed <= 0) return(1)
  count_distribution <- 1   # P(X = 0, 1, ..., n)
  for (p in hit_probabilities) {
    count_distribution <- c(count_distribution * (1 - p), 0) + c(0, count_distribution * p)
  }
  if (observed >= length(count_distribution)) return(0)
  sum(count_distribution[(observed + 1):length(count_distribution)])
}

# target_map: columns target, ortholog_group. A target is one OG or a set of
# OGs such as a GO term (an OG may belong to several targets). BH is applied
# across all targets in target_map.
# Whole-genome null: each lineage's mutations (lineage_mutations, counted genome-wide)
# land uniformly across its reference genome (wg_size), and a mutation that lands in a
# gene belongs to the tested class with probability class_fraction (1 for all mutations)
test_target_enrichment <- function(test_mutations, lineage_og_wg_length, target_map, class_fraction = 1) {
  null_hits <- lineage_og_wg_length %>%
    inner_join(target_map, by = "ortholog_group", relationship = "many-to-many") %>%
    group_by(target, lineage) %>%
    summarise(target_length = sum(og_length),
              n_copies = max(n_copies),
              wg_size = first(wg_size),
              lineage_mutations = first(lineage_mutations),
              .groups = "drop") %>%
    mutate(p_hit = -expm1(lineage_mutations * log1p(-class_fraction * target_length / wg_size))) %>%
    group_by(target) %>%
    summarise(n_lineages_at_risk = n(),
              expected_lineages = sum(p_hit),
              max_copies = max(n_copies),
              median_length = median(target_length),
              hit_probabilities = list(p_hit),
              .groups = "drop") %>%
    # targets carried only by lineages without mutations can't be hit, so they aren't tested
    filter(expected_lineages > 0)
  
  observed_hits <- test_mutations %>%
    inner_join(target_map, by = "ortholog_group", relationship = "many-to-many") %>%
    group_by(target) %>%
    summarise(n_lineages_mutated = n_distinct(lineage),
              n_subjects_mutated = n_distinct(subject),
              n_mutations = n(),
              max_mutations_one_lineage = max(table(lineage)),
              effects = paste(sort(unique(effect_term)), collapse = "; "),
              .groups = "drop")
  
  null_hits %>%
    left_join(observed_hits, by = "target") %>%
    mutate(across(c(n_lineages_mutated, n_subjects_mutated, n_mutations, max_mutations_one_lineage),
                  ~ coalesce(as.integer(.x), 0L)),
           obs_exp = n_lineages_mutated / expected_lineages,
           p_value = as.numeric(mapply(poisson_binomial_upper_tail, n_lineages_mutated, hit_probabilities)),
           q_value = p.adjust(p_value, method = "BH")) %>%
    select(-hit_probabilities) %>%
    relocate(target, n_lineages_mutated, expected_lineages, obs_exp, p_value, q_value) %>%
    arrange(p_value, desc(obs_exp))
}

# One GO category: the term-level test plus which OGs carry the mutations
test_go_category <- function(test_mutations, lineage_og_wg_length, go_category, class_fraction = 1) {
  go_map <- og_go_tested %>%
    filter(category == go_category) %>%
    transmute(target = go_term, ortholog_group)
  if (nrow(go_map) == 0) return(NULL)
  
  mutated_ogs <- test_mutations %>%
    inner_join(go_map, by = "ortholog_group", relationship = "many-to-many") %>%
    group_by(target, ortholog_group) %>%
    summarise(n_lineages = n_distinct(lineage), .groups = "drop") %>%
    arrange(target, desc(n_lineages), ortholog_group) %>%
    group_by(target) %>%
    summarise(mutated_ogs = paste0(ortholog_group, " (", n_lineages, ")", collapse = ", "),
              .groups = "drop")
  
  test_target_enrichment(test_mutations, lineage_og_wg_length, go_map, class_fraction) %>%
    left_join(mutated_ogs, by = "target") %>%
    rename(go_term = target) %>%
    left_join(go_term_names, by = "go_term") %>%
    left_join(go_term_size %>% filter(category == go_category) %>% select(go_term, term_size),
              by = "go_term") %>%
    mutate(category = go_category) %>%
    relocate(category, go_term, go_name, term_size, n_lineages_mutated, expected_lineages,
             obs_exp, p_value, q_value, mutated_ogs)
}

# Genome-wide parallelism: number of OGs mutated in >= min_lineages lineages
# when each lineage's ortholog-group mutations are re-placed at random within
# ortholog-group genes (coding-length null, conditioned on each lineage's gene mutations)
simulate_parallel_ogs <- function(lineage_og_length, n_simulations, min_lineages, n_cores) {
  og_levels <- unique(lineage_og_length$ortholog_group)
  lineage_targets <- lapply(split(lineage_og_length, lineage_og_length$lineage), function(x) {
    list(og_index = match(x$ortholog_group, og_levels),
         cumulative_length = cumsum(x$og_length),
         lineage_mutations = x$lineage_mutations[1])
  })
  
  null_counts <- mclapply(seq_len(n_simulations), function(i) {
    hit_og <- unlist(lapply(lineage_targets, function(target) {
      total_length <- target$cumulative_length[length(target$cumulative_length)]
      landing_site <- runif(target$lineage_mutations, 0, total_length)
      unique(target$og_index[findInterval(landing_site, target$cumulative_length) + 1])
    }), use.names = FALSE)
    sum(tabulate(hit_og, nbins = length(og_levels)) >= min_lineages)
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  
  unlist(null_counts)
}

# Run each mutation set ----
og_annotation <- data.frame(
  ortholog_group = gene_ortholog_group[[1]],   # Gene
  gene_name = gene_ortholog_group[[2]],        # Non-unique Gene name
  annotation = gene_ortholog_group[[3]]        # Annotation
)

# descriptive denominator: multi-genome lineages carrying the OG in any genome

og_lineage_presence <- locus_tag_og %>%
  distinct(ortholog_group, genome) %>%
  inner_join(select(strains, genome, lineage), by = "genome") %>%
  filter(lineage %in% subject_strain_reference$lineage) %>%
  group_by(ortholog_group) %>%
  summarise(n_lineages_with_og = n_distinct(lineage), .groups = "drop")

# Whole-genome null input, built once and used for every mutation set:
# all unique mutations per lineage anywhere in the genome (og_mutations holds every
# variant: CDS, RNA genes and intergenic), joined to each reference's gene lengths.
# Every lineage in subject_strain_reference is kept; lineages without mutations get 0
genome_wide_burden <- subject_strain_reference %>%
  select(lineage, reference_genome = reference) %>%
  left_join(count(og_mutations, lineage, name = "lineage_mutations"), by = "lineage") %>%
  mutate(lineage_mutations = coalesce(lineage_mutations, 0L))
lineage_og_wg_length <- inner_join(reference_og_wg_length, genome_wide_burden, by = "reference_genome")


if(file.exists(all_results_list_file)){
  all_results_list <- readRDS(all_results_list_file)
  og_results <- all_results_list$og_results
  go_results <- all_results_list$go_results
  genome_wide_results <- all_results_list$genome_wide_results
}else{
  og_results <- list()
  go_results <- list()
  genome_wide_results <- list()
  
  for (mutation_set in mutation_sets) {
    
    test_mutations <- if(mutation_set == "all"){
      og_mutations
    }else{
      filter(og_mutations, mutation_class == mutation_set)
    }
    
    # share of mutations in ortholog-group genes that belong to this set
    class_fraction <- if (mutation_set == "all") 1 else{
      mean(og_mutations$mutation_class[!is.na(og_mutations$ortholog_group)] == mutation_set)
    }
    
    # coding-based input for the genome-wide parallelism simulation only: this set's
    # mutations in ortholog-group genes (lineage_og_wg_length is not overwritten)
    lineage_burden <- count(filter(test_mutations, !is.na(ortholog_group)),
                            lineage, reference_genome, name = "lineage_mutations")
    
    lineage_og_length <- inner_join(reference_og_wg_length, lineage_burden, by = "reference_genome")
    
    ## Ortholog groups ----
    og_results[[mutation_set]] <- test_target_enrichment(test_mutations, lineage_og_wg_length, og_map, class_fraction) %>%
      rename(ortholog_group = target) %>%
      left_join(og_annotation, by = "ortholog_group") %>%
      left_join(og_lineage_presence, by = "ortholog_group") %>%
      mutate(lineages_mutated = sprintf("%d/%d (%.1f%%)", n_lineages_mutated, n_lineages_with_og,
                                        100 * n_lineages_mutated / n_lineages_with_og)) %>%
      relocate(ortholog_group, gene_name, annotation, lineages_mutated)
    
    ## GO terms, one category at a time ----
    # q_value is BH within a category; q_value_all_go is BH across all three
    go_set_results <- bind_rows(lapply(go_categories, function(go_category) {
      test_go_category(test_mutations, lineage_og_wg_length, go_category, class_fraction)
    }))
    if (nrow(go_set_results) > 0) {
      go_set_results$q_value_all_go <- p.adjust(go_set_results$p_value, method = "BH")
      go_set_results <- relocate(go_set_results, q_value_all_go, .after = q_value)
    }
    go_results[[mutation_set]] <- go_set_results
    
    ## Genome-wide parallelism ----
    null_parallel <- simulate_parallel_ogs(lineage_og_length, n_null_simulations,
                                           min_lineages_parallel, n_cores)
    observed_parallel <- sum(og_results[[mutation_set]]$n_lineages_mutated >= min_lineages_parallel)
    genome_wide_results[[mutation_set]] <- data.frame(
      mutation_set = mutation_set,
      n_lineages_with_mutations = nrow(lineage_burden),
      n_mutations = nrow(test_mutations),
      min_lineages = min_lineages_parallel,
      observed_ogs = observed_parallel,
      null_median_ogs = median(null_parallel),
      null_95_interval = paste(quantile(null_parallel, c(0.025, 0.975), names = FALSE), collapse = "-"),
      empirical_p = (1 + sum(null_parallel >= observed_parallel)) / (1 + n_null_simulations)
    )
  }
  
  all_results_list <- list(
    og_results = og_results,
    go_results = go_results,
    genome_wide_results = genome_wide_results
  )
  saveRDS(all_results_list, file = all_results_list_file)
}


# Output ----

# og_mutations holds every unique mutation per lineage, including non-coding ones
lineage_mutation_counts <- og_mutations %>%
  group_by(lineage) %>%
  summarise(n_mutations = n(),
            n_protein_altering = sum(mutation_class == "protein_altering"),
            n_synonymous = sum(mutation_class == "synonymous"),
            n_non_coding = sum(mutation_class == "non_coding"),
            n_in_ortholog_groups = sum(!is.na(ortholog_group)),
            .groups = "drop")

lineage_summary <- subject_strain_reference %>%
  select(lineage, subject, strain, reference_genome = reference) %>%
  left_join(lineage_mutation_counts, by = "lineage") %>%
  mutate(across(starts_with("n_"), ~ coalesce(.x, 0L))) %>%
  arrange(desc(n_mutations))


output_sheets <- c(og_results,
                   setNames(go_results, paste0("GO_", names(go_results))),
                   list(genome_wide = bind_rows(genome_wide_results),
                        lineages = lineage_summary))

write_xlsx(output_sheets, file = "SNP_Ortholog_Enrichment.xlsx")

## Observed vs expected plot, ortholog groups ----
# Every ortholog group at risk is shown; triangles have no mutated lineage.
# Triangles are drawn first so circles stay on top.
plot_data <- bind_rows(lapply(names(og_results), function(set_name) {
  og_results[[set_name]] %>%
    mutate(mutation_set = factor(set_name, levels = mutation_sets))
})) %>%
  mutate(point_class = factor(case_when(n_lineages_mutated == 0 ~ "Not mutated",
                                        q_value < 0.05 ~ "Mutated, q < 0.05",
                                        TRUE ~ "Mutated"),
                              levels = c("Not mutated", "Mutated", "Mutated, q < 0.05"))) %>%
  arrange(point_class)

plot_data$label <- ifelse(plot_data$point_class == "Mutated, q < 0.05", plot_data$ortholog_group, NA)

expected_exponent <- floor(log10(range(plot_data$expected_lineages)))
expected_range <- 10^(expected_exponent + c(-0.9, 0.9))
expectation_line <- data.frame(
  expected = 10^seq(log10(expected_range[1]), log10(expected_range[2]), length.out = 1000)
)

# Null 95% range: under random placement the number of mutated lineages is close to
# Poisson with mean = expected (the exact Poisson-binomial is slightly narrower).
# P(X <= k) = pgamma(expected, k + 1, lower.tail = FALSE) also holds for non-integer k,
# which gives smooth edges. A whole count sits above the upper edge exactly when its
# one-sided p is below 0.025, and below the lower edge when its lower-tail p is below 0.025.
continuous_poisson_quantile <- function(probability, expected) {
  vapply(expected, function(lambda) {
    uniroot(function(k) pgamma(lambda, k + 1, lower.tail = FALSE) - probability,
            interval = c(-1 + 1e-9, lambda + 100), tol = 1e-10)$root
  }, numeric(1))
}

expectation_band <- data.frame(
  expected = expectation_line$expected,
  lower = pmax(continuous_poisson_quantile(0.025, expectation_line$expected), -0.5),
  upper = continuous_poisson_quantile(0.975, expectation_line$expected) + 1
)

# Benjamini-Hochberg bar for the next discovery in each panel: with k genes already at
# q < 0.05 out of m tested, another gene needs p <= 0.05 * (k + 1) / m
# (0.05 / m for the first, 2 * 0.05 / m for the second, ...). The dotted line is the
# count at which a gene's Poisson p reaches that bar; genes above it have q < 0.05.
fdr_boundary <- bind_rows(lapply(names(og_results), function(set_name) {
  set_results <- og_results[[set_name]]
  next_discovery_p <- 0.05 * (sum(set_results$q_value < 0.05) + 1) / nrow(set_results)
  data.frame(mutation_set = factor(set_name, levels = mutation_sets),
             expected = expectation_line$expected,
             boundary = continuous_poisson_quantile(1 - next_discovery_p, expectation_line$expected) + 1)
}))

og_plot <- ggplot(plot_data, aes(x = expected_lineages, y = n_lineages_mutated)) +
  # band first so it sits behind the points and doesn't tint their colours
  geom_ribbon(data = expectation_band, aes(x = expected, ymin = lower, ymax = upper,linetype = "band"),
              inherit.aes = FALSE, fill = "#F2C94C", alpha = 0.25) +
  geom_point(aes(shape = point_class, fill = point_class, alpha = point_class),
             colour = "transparent",
             stroke = 0,
             size = 2) +
  geom_line(data = expectation_line, aes(x = expected, y = expected, linetype = "expected"),
            colour = "#B8860B", linewidth = 0.7, alpha = 0.5) +
  geom_line(data = fdr_boundary, aes(x = expected, y = boundary, linetype = "threshold"),
            colour = "#B22222", linewidth = 0.7, alpha = 0.5) +
  # gene names drawn last so they sit on top; each is moved one decade left
  # and an arrow points back to its gene
  ggrepel::geom_text_repel(aes(label = label), colour = "#EB5757", size = 5, fontface = "italic",
                           nudge_x = -1, min.segment.length = 0, point.padding = 0.2,
                           arrow = arrow(length = unit(1.5, "mm")),
                           seed = 1, na.rm = TRUE) +
  scale_x_log10() +
  # legend text only; the class names in plot_data stay the same
  scale_shape_manual(values = c("Not mutated" = 21, "Mutated" = 21, "Mutated, q < 0.05" = 21),
                     breaks = c("Not mutated", "Mutated", "Mutated, q < 0.05"),
                     labels = c("Not mutated", "Mutated", "Convergent Mutation"),
                     name = NULL) +
  # light blue (reassuring), light orange (caution), soft red (alarm)
  scale_fill_manual(values = c("Not mutated" = "#5DA9DD", "Mutated" = "#F2994A", "Mutated, q < 0.05" = "#EB5757"),
                    breaks = c("Not mutated", "Mutated", "Mutated, q < 0.05"),
                    labels = c("Not mutated", "Mutated", "Convergent Mutation"),
                    name = NULL) +
  # dense classes stay translucent; convergent points are opaque so the red reads clearly
  scale_alpha_manual(values = c("Not mutated" = 0.45, "Mutated" = 0.45, "Mutated, q < 0.05" = 1),
                     guide = "none") +
  # the reference curves are mapped to linetype only to give them a legend
  scale_linetype_manual(values = c(expected = "solid", band = "solid", threshold = "solid"),
                        breaks = c("expected", "band", "threshold"),
                        labels = c("Null Expectation", "95% Null Interval", "Convergence Threshold"),
                        name = NULL) +
  # shape and fill share an order so the point legend stays merged; the curve legend goes below
  guides(shape = guide_legend(order = 1, override.aes = list(size = 3)),
         fill = guide_legend(order = 1),
         linetype = guide_legend(order = 2)) +
  scale_y_continuous(breaks = scales::breaks_width(1)) +
  # panel titles only; the set names in the data stay the same
  facet_wrap(~ mutation_set, labeller = as_labeller(c(protein_altering = "Nonsynonymous",
                                                      synonymous = "Synonymous",
                                                      all = "All"))) +
  coord_cartesian(ylim = c(-0.5, max(plot_data$n_lineages_mutated) + 0.5)) +
  labs(x = "Expected Mutated Strain/Subject", y = "Observed Mutated Strain/Subject") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12.5),
    panel.grid.major.y = element_line(linewidth = 0.35,color = "grey70",linetype = "dotted"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

ggsave("SNP_Ortholog_Enrichment.pdf", og_plot, width = 12.5, height = 4.5)

## Top GO terms per category, protein-altering mutations ----
go_plot_data <- go_results[["protein_altering"]]
if (!is.null(go_plot_data) && any(go_plot_data$n_lineages_mutated > 0)) {
  go_plot_data <- go_plot_data %>%
    filter(n_lineages_mutated > 0) %>%
    group_by(category) %>%
    slice_min(p_value, n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(category = factor(category, levels = go_categories),
           term_label = reorder(paste0(go_name, " (", go_term, ")"), obs_exp),
           # same classes as the ortholog group plot; every term shown here is mutated
           point_class = factor(ifelse(q_value < 0.05, "Mutated, q < 0.05", "Mutated"),
                                levels = c("Mutated", "Mutated, q < 0.05")))
  
  go_plot <- ggplot(go_plot_data, aes(x = obs_exp, y = term_label)) +
    # obs/exp = 1 is the random-mutation expectation, styled like the expectation line in og_plot;
    # its own one-row data frame draws it once per panel (repeated lines would stack their alpha)
    geom_vline(data = data.frame(xintercept = 1), aes(xintercept = xintercept, linetype = "expected"),
               colour = "#B8860B", linewidth = 0.7, alpha = 0.5) +
    geom_point(aes(size = n_lineages_mutated, shape = point_class, fill = point_class, alpha = point_class),
               colour = "transparent",
               stroke = 0) +
    scale_x_log10() +
    # smallest bubble matches the point size in og_plot
    scale_size(range = c(2, 6)) +
    # legend text only; the class names in go_plot_data stay the same
    scale_shape_manual(values = c("Mutated" = 21, "Mutated, q < 0.05" = 21),
                       breaks = c("Mutated", "Mutated, q < 0.05"),
                       labels = c("Mutated", "Convergent Mutation"),
                       name = NULL) +
    scale_fill_manual(values = c("Mutated" = "#F2994A", "Mutated, q < 0.05" = "#EB5757"),
                      breaks = c("Mutated", "Mutated, q < 0.05"),
                      labels = c("Mutated", "Convergent Mutation"),
                      name = NULL) +
    # lighter orange keeps convergent terms in front; 0.6 rather than og_plot's 0.35
    # because the bubbles don't overlap here and small ones would almost disappear
    scale_alpha_manual(values = c("Mutated" = 0.6, "Mutated, q < 0.05" = 1),
                       guide = "none") +
    # the expectation line is mapped to linetype only to give it a legend
    scale_linetype_manual(values = c(expected = "solid"),
                          labels = c(expected = "Null Expectation"),
                          name = NULL) +
    # size keys need their own fill: the point fill is mapped, so the size legend would otherwise be blank
    guides(shape = guide_legend(order = 1, override.aes = list(size = 3)),
           fill = guide_legend(order = 1),
           size = guide_legend(order = 2, override.aes = list(shape = 21, fill = "#808285")),
           linetype = guide_legend(order = 3)) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Observed / Expected Mutated Strain/Subject", y = NULL, size = "Mutated Strain/Subject") +
    theme(
      strip.background = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.5, color = "grey70", linetype = "dotted"),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
  ggsave("SNP_GO_Enrichment.pdf", go_plot, width = 10, height = 1.5 + 0.25 * nrow(go_plot_data))
}
