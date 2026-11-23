# MAGIC *Staphylococcus aureus* Colonization Analysis Repository
Analysis scripts, inputs, and outputs for MAGIC (Microbiome, Antibiotics, and Growth Infant Cohort), a longitudinal study of *Staphylococcus aureus* colonization in 159 infants over the first two years of life. Colonization is resolved at the strain level with a phylogenomic definition that combines culture results with strain assignments from whole-genome sequencing.

## Repository Structure
Each analysis folder contains `scripts/`, `input/`, and `output/`. Each script sets its working directory to its own `output/` folder, reads from `input/`, and checks for the required R packages before running. 

| Folder | Analysis |
| --- | --- |
| `Subject_Colonization_Plot` | Per-subject swab results and persistence |
| `Colonization_Redefinition` | Colonization states and persistence under each definition |
| `Hidden_Markov_Model` | Misclassification and persistence by hidden Markov model |
| `Within_Host_Substitution_Rates` | SNP accumulation against persistence duration |
| `Convergent_Evolution` | Convergent within-host mutation |
| `Antibiotics_Assessment` | Antibiotic exposure and colonization transitions |
| `Antibiotics_Visual` | Colonization and antibiotic administration by subject |

## Genome Assembly and Quality Control
Raw sequencing reads were processed with Sunbeam (v4.3.7) for adapter trimming and host decontamination, and isolate genomes were assembled with SPAdes (v3.15.5).

Assembly quality was assessed using CheckM (v1.2.2), applying the following criteria:
- ≥95% CheckM completeness
- ≤5% CheckM contamination
- Lineage classification as "Staphylococcus (UID301)"

Species-level contamination was evaluated using Mash (v2.3). Of 313 culture-positive samples sequenced from 87 infants, 12 failed sequencing, 7 were contaminated, and 5 were not *S. aureus*, leaving 349 high-quality assemblies (181 nares, 37 axilla, 115 groin, 16 stool) from 289 samples of 78 infants for downstream analyses.

Assemblies were annotated with Bakta (v1.11.4, full database v6.0), and the pan-genome was constructed with Panaroo (v1.5.1).

## Strain Assignment
Strains were assigned with THRESHER (v0.4.1; https://github.com/microbialARC/THRESHER), which infers a phylothreshold, a SNP threshold corrected for phylogenetic topology, for each genetic background rather than applying a single fixed cutoff. The 349 assemblies were run with `--species sau --endpoint plateau --plateau_length 20`, which also returns MLST, clonal complex, and MRSA/MSSA status. This identified 82 strains across 78 infants.

The THRESHER results (`MAGIC_thresher_plateau.RDS`) are an input to the subject colonization, colonization redefinition, within-host substitution rate, and convergent evolution analyses.

## Colonization Redefinition
### Subject Colonization Plot
This script plots swab results, strain assignments, and colonization persistence for each of the 87 subjects with at least one culture-positive sample, and summarizes persistence across subjects.

#### Visualization:
- Upper panel: swab results at each body site (nares, axilla, groin, and stool where collected) across the 13 time points. Positive swabs are shown as circles colored by strain, and negative swabs as circles with a cross. Samples that failed WGS, were contaminated, or were not *S. aureus* are shown as separate categories.
- Lower panel: time points counted as colonized under the culture and phylogenomic definitions.

#### Output:
- `Subject<ID>.pdf`: one plot per subject
- `All_Subjects_Persistence.pdf`: culture (red) and phylogenomic (blue) persistence for all subjects
- `Modified_Subjects_Persistence.pdf`: subjects whose culture and phylogenomic persistence differ
- `persistence_list.RDS` and `persistence_sum.RDS`: per-subject results, used as input by `Colonization_Redefinition`

### Colonization Pattern and Persistence Comparison
This script assigns each subject a colonization state under the culture and phylogenomic definitions and compares persistence duration between the two.

| Colonization state | Culture | Phylogenomic |
| --- | --- | --- |
| Persistent | ≥2 consecutive positive time points | Same strain at ≥2 time points |
| Intermittent | >1 positive time point, none consecutive | >1 positive time point, none the same strain |
| Once | 1 positive time point | 1 positive time point |
| Never | No positive time points | No positive time points |

Under the phylogenomic definition, culture positives identified by WGS as another species are corrected to negative, and subjects are classified as undetermined when a time point that failed WGS could change their state. Persistence is the number of days between the first and last time point of each persistent period, summed per subject.

#### Visualization:
- `Cohort_Colonization_Patterns.pdf`: alluvial plot of subject counts in each colonization state under the culture and phylogenomic definitions
- `Violin_Persistence_Comparison.pdf`: persistence per subject under each definition, with lines connecting the same subject colored by the difference in days

## Hidden Markov Model Analysis and Visualization
This script fits a two-state (uncolonized, colonized) continuous-time hidden Markov model with misclassification using the msm R package. Observation times are days after birth, and a model is fitted separately for each colonization definition:

#### Visualization:
Estimates are shown with 95% confidence intervals from 100 bootstrap replicates.
- `Misclassification_Probability.pdf`: false-negative and false-positive probabilities
- `Length_Stay.pdf`: expected days colonized during the first 720 days of life, starting uncolonized
- `Transition_Probability.pdf`: one-day probabilities of acquisition and clearance

Fitted models and bootstrap estimates are saved in `msm_results.RDS`.

## Within-Host Evolution
### Within-Host Substitution Rates
This script regresses within-host SNP accumulation against colonization duration under each definition for 55 strain/subject timelines (54 strains in 52 subjects).
- SNPs were called with Snippy (v4.6.0) by mapping each later genome to the earliest genome of the same strain/subject (`input/snps_tab/strain<ID>_subject<ID>_<genome>_snps.tab`). SNP sites are pooled across the genomes of a strain/subject and counted once.
- Duration is calculated from swab collection dates. Under the culture definition, it is the summed span of consecutive positive time points; under the phylogenomic definition, it is the span from the first to the last detection of the strain.
- Unique SNPs are regressed on each duration by ordinary least squares, and the slope is the apparent substitution rate in SNPs/day.

#### Visualization:
- `Persistence_Unique_Snp_Comparison.pdf`: unique SNPs against culture (left) and phylogenomic (right) duration, with least-squares fits, 95% confidence bands, slope, adjusted R², and p value

### Convergent Evolution
This script tests whether within-host mutations arose in the same ortholog groups across independent strain/subjects more often than expected by chance.

#### Visualization:
- `SNP_Ortholog_Enrichment.pdf`: observed against expected numbers of mutated strain/subjects for each ortholog group, for nonsynonymous, synonymous, and all mutations, with the null expectation, 95% null interval, and convergence threshold
- `SNP_GO_Enrichment.pdf`: top GO terms in each category for nonsynonymous mutations
- `SNP_Ortholog_Enrichment.xlsx`: results for ortholog groups, GO terms, the genome-wide test, and mutation counts per strain/subject

Test results are saved in `input/all_results_list_wg.RDS`.

## Antibiotic Exposure and Colonization
### Hidden Markov Model with Antibiotic Exposure
This script adds antibiotic exposure to the two-state hidden Markov model to model how acquisition and clearance of colonization vary with exposure.

#### Visualization:
- `Phylogenomics_Transition_Probability_Active_Others.pdf`: daily probabilities of acquisition and clearance by exposure
- `Phylogenomics_Persistence_Active_Others.pdf`: expected days colonized by exposure
- `Culture_Transition_Probability_Active_Others.pdf` and `Culture_Persistence_Active_Others.pdf`: the same estimates under the culture definition

Fitted models and bootstrap estimates are saved in `hmm_abx_results.RDS`.

### Antibiotic Administration Plot
This script plots colonization status and antibiotic administration against days after birth for each subject with recorded swab collection dates.

#### Visualization:
- Each row is one subject, with culture (red) and phylogenomic (blue) colonization on separate tracks. Positive time points are dots connected by lines within each colonization period, and negative swabs are circles with a cross.
- Antibacterial courses (`abx_data.xlsx`) are shaded from start to stop day and classified with `abx_class.xlsx`: *S. aureus*-active agents, including mupirocin, in green and others in grey. Courses that continue past day 810 or have no recorded stop day extend to the end of the plot.
- Swab days after birth are calculated from collection dates and dates of birth. Time points without a recorded collection date are placed at their scheduled age (day 0, weekly through day 28, then every 90 days through day 720).

#### Output:
- `Antibiotics_Visual.pdf`: colonization and antibiotic administration for all plotted subjects

**Note: Participant dates of birth are withheld as protected health information (PHI).**
