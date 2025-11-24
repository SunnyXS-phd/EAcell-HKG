#---------#
#Author: Shiqi Huang
# Perform HouseKeeping Gene stability ranking
# using commonly available methods: ddCt, NormFinder, geNorm & BestKeeper, as a validation of RefFinder outputs
# the NormFinder source R script is from https://www.moma.dk/software/normfinder 
# R packages used: NormqPCR <geNorm, NormFinder> and ctrlGene <geNorm, BestKeeper>
#---------#

library(tidyverse)
library(ggplot2)
library(RColorBrewer)

#---------#
# deltaCt method is to calculate the Ct difference between a target gene and all other genes (pair by sample),
# then calculate the SD of each pairs of genes among all samples tested,
# lastly, take the mean SD of this target gene vs. all other genes, and rank from lowest(most stable)
#---------#
library(purrr)
# Load cleaned Cq data
# Transform Table
Cq.mean.dCt <- Cq.mean %>%
  pivot_wider(names_from = Detector, values_from = Cq.mean)

# Extract gene column names 
genes <- colnames(Cq.mean.dCt)[3:ncol(Cq.mean.dCt)]

# Generate all possible pairs of genes
gene_pairs <- combn(genes, 2, simplify = FALSE) #combn(): Generates all possible combinations of the gene columns.

## 1st Function to calculate dCt for each pair
calculate_dCt <- function(pair, data) {
  gene1 <- sym(pair[1])  # Convert string to symbol
  gene2 <- sym(pair[2])
  
  # Calculate dCt for each sample
  dCt_values <- data %>%
    mutate(Pair = paste(quo_name(gene1), "vs.", quo_name(gene2)),
           dCt = !!gene1 - !!gene2) %>%
    select(Sample, Pair, dCt)
  
  return(dCt_values)
}

# Apply the function to all gene pairs and combine results
dCt.results.all <- purrr::map_dfr(gene_pairs, ~calculate_dCt(.x, data = Cq.mean.dCt))

## 2nd, calculate Mean and SD for 
dCt.results.mean <- dCt.results.all %>%
  group_by(Pair) %>%
  summarise(
    Mean_dCt = mean(dCt, na.rm = TRUE),
    SD = sd(dCt, na.rm = TRUE)
  )

## Last, calculate mean SD across all pairs for any given target gene
# Initialize an empty data frame to store results
averageSD.results <- data.frame(Gene = character(), Average_SD = numeric())

for (target in genes) {
  # Filter out those contains the target gene
  tmp <- dCt.results.mean %>%
    filter(grepl(target, Pair))
  
  # Calculate average SD for the filtered pairs
  average_sd <- mean(tmp$SD, na.rm = TRUE)
  
  # Append results to average_sd_df
  averageSD.results <- rbind(averageSD.results, data.frame(Gene = target, Average_SD = average_sd))
}

# Plot
p.dCt <- ggplot(averageSD.results, aes(x = reorder(Gene, Average_SD), y = Average_SD)) +
  geom_col(fill = "plum") +
  theme_classic() +
  labs(x = "Gene", y = "Average standard deviation of deltaCt values")

#---------#
# NormFinder using source code
#---------#
# Retrieve the source code as a function:
source("r.NormOldStab5.txt")

# Import data with samples as columns (names in first row) and genes as rows (names in first column)

# Transform that into the correct format
meanCt.NormFinder <- meanCt.long |>
  select(Sample, Gene, Cq.mean) |>
  pivot_wider(names_from = Sample,
              values_from = Cq.mean
  )

Group_info <- meanCt.long |>
  select(Sample, Group) |>
  distinct() |>
  mutate(Group_code = case_when(
    Group == "D10+MG" ~ 1,
    Group == "D10-MG" ~ 2,
    Group == "D4+MG" ~ 3,
    Group == "D4-MG" ~ 4,
    TRUE ~ NA_integer_ # Handle unexpected cases
  )) |>
  select(Sample, Group_code) |> # otherwise will have additional column and rows
  pivot_wider(names_from = Sample, values_from = Group_code) |>
  mutate(Gene = "group") |>
  select(Gene, everything()) # To make the Gene column first

meanCt.NormFinder <- bind_rows(meanCt.NormFinder, Group_info) 

## Save this as the txt file required for the function
write.table(meanCt.NormFinder, file = "EA_meanCt_normFinder.txt", sep = "\t", row.names = FALSE)

# Apply function to get results
Result = Normfinder("EA_meanCt_normFinder.txt")

Result$Ordered 

# Plot results
res.NF1 <- as.data.frame(Result$Ordered) |> rownames_to_column(var = "Gene")

p.NF1 <- ggplot(res.NF1, aes(x = reorder(Gene, Stability), y = Stability)) +
  geom_col(fill = "skyblue") +
  theme_classic() +
  xlab("Gene")

#---------#
# NormqPCR package
#---------#

library(NormqPCR)
  
# Upload the tab-deliminated file for NormqPCR with ReadqPCR::read.qPCR() to creat qPCRBatch R object
EA.qPCRBatch = ReadqPCR::read.qPCR(filename = "CompiledMeanCq.txt")

Cq.mean <- read.table(file = "CompiledMeanCq.txt", sep = "\t", header = TRUE) 

## Create a phenoData/sample classification/grouping
sample.info <- data.frame(Sample = unique(Cq.mean$Sample),
                          Group = unique(Cq.mean$Group)) 

### Convert to annotatedDataDrame format
sample.info <- AnnotatedDataFrame(sample.info)
sampleNames(sample.info) <- sampleNames(EA.qPCRBatch)

EA.qPCRBatch <- ReadqPCR::read.qPCR(filename = "CompiledMeanCq.txt", phenoData = sample.info)

## Compute stability
Group <- pData(EA.qPCRBatch)[,"Group"]

### NormFinder results
res.NF <- selectHKs(EA.qPCRBatch, method = "NormFinder", Symbols = featureNames(EA.qPCRBatch), group = Group, log = FALSE, minNrHKs = 26)
#'minNrHKs' must be smaller than 'ncol(x)'
res.NF$ranking

# Should be the same as using:
res.NF2 <- stabMeasureRho(EA.qPCRBatch, group = Group, log = FALSE) |>
  sort()
print(res.NF2)

### geNorm results
res.gN <- selectHKs(EA.qPCRBatch, method = "geNorm", Symbols = featureNames(EA.qPCRBatch), log = FALSE, trace = FALSE)
res.gN$ranking

## Plot results
res.NF2 <- tibble::enframe(res.NF2, name = "Gene", value = "Rho")

p.NF2 <- ggplot(res.NF2, aes(x = reorder(Gene,Rho), y = Rho)) +
  geom_col(fill = "skyblue") +
  theme_classic() +
  labs(x = "Gene")

## Compiled the res.gN[[3]] meanM with ranking res.gN[[1]] mannually in excel
## Import this for plotting
res.gN.1 <- read.xlsx("EA_geNorm_results_NormqPCR.xlsx", sheet = "Sheet1")

p.gN2 <- ggplot(res.gN.1, aes(x = reorder(Gene,meanM), y = meanM)) +
  geom_col(fill = "blueviolet") +
  theme_classic() +
  labs(x = "Gene", y = "mean M")

# Plot V values 
v_data <- data.frame(
  comparison = factor(names(res.gN$variation), 
                      levels = rev(names(res.gN$variation))),  # Reserve order
  V_value = as.numeric(res.gN$variation)
)

ggplot(v_data, aes(x = comparison, y = V_value)) +
  geom_col(fill = "steelblue", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = sprintf("%.6f", V_value)), 
            vjust = -0.5, hjust = 0,
            size = 2.5, angle = 30, fontface = "bold") +
  labs(
    x = "Pairwise Comparison (n/n+1)",
    y = "Pairwise Variation (V value)",
    caption = "V < 0.15 indicates optimal number of reference genes is reached"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(face = "bold", size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


#---------#
# ctrlGene
#---------#
library(ctrlGene)

## Load data
Cq.mean.ctrlGene <- Cq.mean |>
  pivot_wider(names_from = Detector, values_from = Cq.mean) |>
  column_to_rownames(var = "Sample") |>
  dplyr::select(!Group)

res.BK.cG <-bestKeeper(Cq.mean.ctrlGene)
res.gN.cG <- geNorm(Cq.mean.ctrlGene, ctVal = TRUE)  

## Plot results
p.gN.cG <- ggplot(res.gN.cG, aes(x = reorder(Genes,Avg.M), y = Avg.M)) +
  geom_col(fill = "blueviolet") +
  theme_classic() +
  labs(x = "Gene", y = "Average M")


res.BK.cG.1 <- res.BK.cG[[1]] |>
  as.data.frame() |>
  t() |>
  as.data.frame() |>
  rownames_to_column(var = "Gene")
colnames(res.BK.cG.1)[colnames(res.BK.cG.1) == "SD[+/- CP]"] <- "SD_CP"

p.BK.cG <- ggplot(res.BK.cG.1, aes(x = reorder(Gene,SD_CP), y = SD_CP)) +
  geom_col(fill = "orchid") +
  theme_classic() +
  labs(x = "Gene", y = "CP Standard Deviation")

