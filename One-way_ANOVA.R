#---------#
#Data Processing & Analysis for Figure 3 using one-way ANOVA
#author: "Shiqi Huang"
#---------#

# Load necessary libraries
library(tidyverse)
library(openxlsx)
library(purrr)
library(outliers)
library(car)   # For Levene's test
library(DescTools)  # For post-hoc tests
library(PMCMRplus) # For post-hoc tests
library(ggpubr) # https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/ 
library(rstatix) # https://rpkgs.datanovia.com/rstatix/ 

## Import data

## Examine ΔCt
## Make a data frame with deltaCt to different combinations of HKGs
## Define DEGs
DEGs <- unique(Cq.mean$Detector)
DEGs <- DEGs[!DEGs%in%c("CAPZB", "FBXO7", "YWHAZ", "18S", "ACTB", "HPRT1")]

## Define all HKG combinations with labels
HKG_combinations <- list(
  sHKG3 = c("CAPZB", "FBXO7", "YWHAZ"),
  cHKG3 = c("18S", "ACTB", "HPRT1"),
  sHKG2 = c("CAPZB", "FBXO7"),
  cHKG2 = c("18S", "HPRT1"),
  CAPZB = c("CAPZB"),
  `18S` = c("18S")
)

## Function to calculate pseudo-housekeeping Cq for a given set of HKGs
calculate_pseudo_hk <- function(Cq_data, hkgs, label) {
  Cq_data %>%
    filter(Detector %in% hkgs) %>%
    group_by(Plate, Sample) %>%
    summarize(
      !!paste0("pseudoCq_", label) := mean(Cq.mean, na.rm = TRUE), 
      .groups = "drop"
    )
}

## Apply the function for each combination
pseudo_hk_results <- lapply(
  names(HKG_combinations),
  function(label) calculate_pseudo_hk(Cq.mean, HKG_combinations[[label]], label)
) %>%
  reduce(full_join, by = c("Plate", "Sample")) # Merge results side-by-side

## Join pseudo-housekeeping Cq with the original data
deltaCq.results <- Cq.mean %>%
  filter(Detector %in% DEGs) |>
  left_join(pseudo_hk_results, by = c("Sample", "Plate")) |>
  pivot_longer(
    cols = starts_with("pseudoCq_"), # Select pseudo-Cq columns
    names_to = "HKG_combination",    # Create column for HKG combination names
    values_to = "pseudoHKG_Cq"
  ) %>%
  mutate(
    ΔCq = Cq.mean - pseudoHKG_Cq,    # Calculate deltaCq (delta in Greek)
    ΔCq_label = paste0("Δ", Detector, "_", sub("pseudoCq_", "", HKG_combination))
  ) %>%
  select(Plate, Sample, Detector, Cq.mean, ΔCq_label, ΔCq, Days, MG)

# Make a data frame with ΔΔCt to different combinations of HKGs
deltadelta_Cq.full <- deltaCq.results %>%
  group_by(Plate, Detector, ΔCq_label) %>%
  mutate(
    RE = 2^(-ΔCq),
    control_deltaCq = max(ΔCq[Days == "D4" & MG == "+MG"], na.rm = TRUE),  # Set control
    ΔΔCq = ΔCq - control_deltaCq,  # Calculate ΔΔCq,
    FC = 2^(-ΔΔCq) # Calculate 2^-ΔΔCq
  ) %>%
  ungroup()

# Extract the Sample and Relative Expression data
dCq <- deltadelta_Cq.full |>
  select(Plate, Condition, Biorep, Detector, ΔCq_label, ΔCq, RE) 

# Convert Sample to factor and set Day4+MG as the baseline level, similar to RNAseq
dCq <- dCq %>%
  mutate(Condition = factor(Condition, levels = c("D4.+MG", "D10.+MG", "D10.-MG", "D4.-MG")))

# Assess Assumptions of ANOVA

## 1 Check for Outlier
### Add a unique row identifier for global indexing
dCq <- dCq %>%
  mutate(global_row = row_number())

# Pipeline to create df.outlier
dCq.outlier <- dCq %>%
  group_by(ΔCq_label) %>%
  group_modify(~ {
    current_data <- .x
    
    # 1. IQR-based outliers (grouped by factorA and factorB)
    iqr_outliers <- current_data %>%
      group_by(Condition) %>%
      rstatix::identify_outliers(RE) %>%
      dplyr::mutate(is.outlier_iqr = is.outlier,
                    is.extreme_iqr = is.extreme) %>%
      dplyr::select(global_row, is.outlier_iqr, is.extreme_iqr)
    
    # 2. Model-based outliers (car::outlierTest)
    model <- lm(RE ~ Condition, data = current_data)
    
    # Get Bonferroni p-values
    outlier_test <- tryCatch(
      car::outlierTest(model),
      error = function(e) NULL
    )
    
    if (!is.null(outlier_test)) {
      # Extract row indices relative to the current group
      outlier_indices <- as.numeric(names(outlier_test$rstudent))
      
      # Safely extract Bonferroni p-values (handles both vectors and matrices)
      bonferroni_p <- as.numeric(outlier_test$bonf.p)
      
      # Map to global_row and p-values
      bonf_p <- tibble(
        global_row = current_data$global_row[outlier_indices],
        bonferroni_p = bonferroni_p
      )
    } else {
      bonf_p <- tibble(global_row = numeric(), bonferroni_p = numeric())
    }
    
    # 3. Cook's distance for all rows in the current dataset
    cooks_d <- augment(model) %>%
      dplyr::mutate(global_row = current_data$global_row) %>%
      dplyr::select(global_row, cooks_d = .cooksd)
    
    # Combine all results
    current_data %>%
      left_join(iqr_outliers, by = "global_row") %>%
      left_join(bonf_p, by = "global_row") %>%
      left_join(cooks_d, by = "global_row")
  }) %>%
  ungroup() %>% 
  dplyr::mutate(
    is.outlier_iqr = replace_na(is.outlier_iqr, FALSE),
    is.extreme_iqr = replace_na(is.extreme_iqr, FALSE),
    bonferroni_p = replace_na(bonferroni_p, 1) # Replace NA with 1 (non-significant, since this column is num)
  )

### View results
outlier <- dCq.outlier %>%
  filter(is.outlier_iqr | is.extreme_iqr | (bonferroni_p <= 0.05) | (cooks_d >= 0.5)) %>%
  select(ΔCq_label, Condition, Biorep, RE, 
         is.outlier_iqr, is.extreme_iqr, bonferroni_p, cooks_d)

### Remove outliers: Tier 1: extreme : ≥ 3xIQR or cooks > 1
outlier %>% dplyr::filter(is.extreme_iqr | cooks_d > 1) %>% print(n = Inf)

dCq.outlierRemoved <- dCq.outlier %>%
  dplyr::filter(!(is.extreme_iqr | cooks_d > 1))

### Remove outliers: Tier 2: bonf_p < 0.05 && cooks > 0.5

outlier %>% dplyr::filter(bonferroni_p <= 0.05 & cooks_d >= 0.5) %>% print(n = Inf) 

dCq.outlierRemoved %>% dplyr::filter(bonferroni_p <= 0.05 & cooks_d >= 0.5) %>% print(n = Inf)

dCq.outlierRemoved <- dCq.outlierRemoved %>%
  dplyr::filter(!(bonferroni_p <= 0.05 & cooks_d >= 0.5)) 

### See sample size per group to ensure n≥3
dCq.outlierRemoved %>% dplyr::count(ΔCq_label, Condition) %>% dplyr::filter(n < 3) %>% print(n = Inf)

## 2. Test Normality & 3. Equal Variance
res.assumption.1 <- dCq.outlierRemoved %>%
  group_by(ΔCq_label) %>%
  summarise(
    # Perform ANOVA
    ANOVA = list(aov(RE ~ Condition)),
    # Check Assumptions
    ShapiroP = list(shapiro.test(residuals(aov(RE ~ Condition)))$p.value),
    LeveneP = list(leveneTest(RE ~ Condition)$`Pr(>F)`[1])
  )

# Highlight significant results. ie. not normal and/or not equal variance
res.assumption.1 |>
  dplyr::filter(ShapiroP < 0.05 | LeveneP < 0.05) |>
  print(n = Inf) 

res.assumption.1 |>
  dplyr::filter(ShapiroP < 0.05) |>
  print(n = Inf) 

res.assumption.1 |>
  dplyr::filter(ShapiroP >= 0.05 & LeveneP < 0.05) |>
  print(n = Inf) 

## 4. Transform the data if any

# Apply log transformation on original data
dCq.transformed <- dCq %>%
  filter(Detector %in% "conditions based on res.assumption.1") %>%
  mutate(LogRE = log10(RE)) 

## 5. Re-test
# Re-test outlier
dCq.transformed.outlier <- dCq.transformed %>%
  group_by(ΔCq_label) %>%
  group_modify(~ {
    current_data <- .x
    
    # 1. IQR-based outliers (grouped by factorA and factorB)
    iqr_outliers <- current_data %>%
      group_by(Condition) %>%
      rstatix::identify_outliers(LogRE) %>%
      dplyr::mutate(is.outlier_iqr = is.outlier,
                    is.extreme_iqr = is.extreme) %>%
      dplyr::select(global_row, is.outlier_iqr, is.extreme_iqr)
    
    # 2. Model-based outliers (car::outlierTest)
    model <- lm(LogRE ~ Condition, data = current_data)
    
    # Get Bonferroni p-values
    outlier_test <- tryCatch(
      car::outlierTest(model),
      error = function(e) NULL
    )
    
    if (!is.null(outlier_test)) {
      # Extract row indices relative to the current group
      outlier_indices <- as.numeric(names(outlier_test$rstudent))
      
      # Safely extract Bonferroni p-values (handles both vectors and matrices)
      bonferroni_p <- as.numeric(outlier_test$bonf.p)
      
      # Map to global_row and p-values
      bonf_p <- tibble(
        global_row = current_data$global_row[outlier_indices],
        bonferroni_p = bonferroni_p
      )
    } else {
      bonf_p <- tibble(global_row = numeric(), bonferroni_p = numeric())
    }
    
    # 3. Cook's distance for all rows in the current dataset
    cooks_d <- augment(model) %>%
      dplyr::mutate(global_row = current_data$global_row) %>%
      dplyr::select(global_row, cooks_d = .cooksd)
    
    # Combine all results
    current_data %>%
      left_join(iqr_outliers, by = "global_row") %>%
      left_join(bonf_p, by = "global_row") %>%
      left_join(cooks_d, by = "global_row")
  }) %>%
  ungroup() %>% 
  dplyr::mutate(
    is.outlier_iqr = replace_na(is.outlier_iqr, FALSE),
    is.extreme_iqr = replace_na(is.extreme_iqr, FALSE),
    bonferroni_p = replace_na(bonferroni_p, 1) # Replace NA with 1 (non-significant, since this column is num)
  )



### View results
outlier_Transformed <- dCq.transformed.outlier %>%
  filter(is.outlier_iqr | is.extreme_iqr | (bonferroni_p <= 0.05) | (cooks_d >= 0.5)) %>%
  select(ΔCq_label, Condition, Biorep, RE, 
         is.outlier_iqr, is.extreme_iqr, bonferroni_p, cooks_d)

### Remove outliers: Tier 1: extreme : ≥ 3xIQR or cooks > 1:
outlier_Transformed %>% dplyr::filter(is.extreme_iqr | cooks_d > 1) %>% print(n = Inf)

dCq.transformed.outlierRemoved <- dCq.transformed.outlier %>%
  dplyr::filter(!(is.extreme_iqr | cooks_d > 1)) # from 1956 to 1950


### See sample size per group to ensure n≥3
dCq.transformed.outlierRemoved %>% 
  dplyr::count(ΔCq_label, Condition) %>% dplyr::filter(n < 3) %>% print(n = Inf)


## Retest assumptions
res.assumption.2 <- dCq.transformed.outlierRemoved %>%
  group_by(ΔCq_label) %>%
  summarise(
    # Perform ANOVA
    ANOVA = list(aov(LogRE ~ Condition)),
    # Check Assumptions
    ShapiroP = list(shapiro.test(residuals(aov(LogRE ~ Condition)))$p.value),
    LeveneP = list(leveneTest(LogRE ~ Condition)$`Pr(>F)`[1])
  )
# Highlight significant results. ie. not normal and/or not equal variance
res.assumption.2 |>
  dplyr::filter(ShapiroP < 0.05 | LeveneP < 0.05) |>
  print(n = Inf)

res.assumption.2 |>
  dplyr::filter(ShapiroP < 0.05) |>
  print(n = Inf)

res.assumption.2 |>
  dplyr::filter(ShapiroP >= 0.05 & LeveneP < 0.05) |>
  print(n = Inf) 

# Conduct ANOVA & Post hoc tests

## normal 1way-ANOVA

### Function to perform ANOVA analysis 
anova_auto <- function(data, value) {
  # Fit ANOVA model
  anova_model <- aov(as.formula(paste(value, "~ Condition")), data = data) #as.formula(paste(value, "~ Sample")) to convert the "value" string to formula suitable for aov & levene setting
  
  # Check Assumptions
  shapiro_p <- shapiro.test(residuals(anova_model))$p.value
  levene_p <- leveneTest(as.formula(paste(value, "~ Condition")), data = data)$`Pr(>F)`[1]
  
  # Perform Post-Hoc Test
  tukey <- TukeyHSD(anova_model)
  duncan <- tryCatch(PMCMRplus::duncanTest(anova_model), error = function(e) NULL)  # tryCatch() to handle potential errors that may break the flow
  
  # Return all results
  list(
    ANOVA = summary(anova_model),
    ShapiroP = shapiro_p,
    LeveneP = levene_p,
    TukeyR = tukey,
    DuncanR = duncan
  )
}

#1 For those meet assumptions with raw data
dCq_anovaRE <- dCq.outlierRemoved %>%
  filter(ΔCq_label %in% "conditions")

res.anova.RE <- dCq_anovaRE %>%
  group_split(ΔCq_label) %>%
  set_names(dCq_anovaRE %>% group_by(ΔCq_label) %>% group_keys() %>% pull()) %>%
  map(~anova_auto(.x, "RE"))

#2 For those meet assumptions with LogRE
dCq_anovaLogRE <- dCq.transformed.outlierRemoved %>%
  dplyr::filter("conditions")

res.anova.LogRE <- dCq_anovaLogRE %>%
  group_split(ΔCq_label) %>%
  set_names(dCq_anovaLogRE %>% group_by(ΔCq_label) %>% group_keys() %>% pull()) %>%
  map(~anova_auto(.x, "LogRE"))


### Tidy results
# Extract essential info from the results
tidy_anova_results <- function(result, label) {
  # Extract Tukey p-values for pairwise comparisons
  tukey_results <- if (!is.null(result$TukeyR) && is.matrix(result$TukeyR[[1]])) {
    tukey_df <- as.data.frame(result$TukeyR[[1]])
    comparisons <- rownames(tukey_df)
    # Convert to a named vector, with comparisons as column names
    tukey_df <- setNames(tukey_df[["p adj"]], paste0("TukeyP_", comparisons))
    # Convert to a single-row data frame
    as.data.frame(as.list(tukey_df))
  } else {
    NULL
  }
  
  # Extract Duncan p-values for pairwise comparisons
  duncan_results <- if (!is.null(result$DuncanR) && is.matrix(result$DuncanR$p.value)) {
    duncan_df <- as.data.frame(result$DuncanR$p.value)
    comparisons <- expand.grid(rownames(duncan_df), colnames(duncan_df)) %>%
      apply(1, function(x) paste(x, collapse = "-"))
    # Convert to a named vector, with comparisons as column names
    duncan_vec <- setNames(as.vector(as.matrix(duncan_df)), paste0("DuncanP_", comparisons))
    # Convert to a single-row data frame
    as.data.frame(as.list(duncan_vec))
  } else {
    NULL
  }
  
  # Combine all results into a single tibble row
  base_row <- tibble::tibble(
    Label = label,
    anovaF = result$ANOVA[[1]][["F value"]][1],
    anovaP = result$ANOVA[[1]][["Pr(>F)"]][1],
    ShapiroP = result$ShapiroP,
    LeveneP = result$LeveneP
  )
  # Append Tukey results if they exist
  base_row <- bind_cols(base_row, as_tibble(tukey_results), as_tibble(duncan_results))
}

#### Apply the function
res.anovaRE.export <- purrr::imap_dfr(res.anova.RE, tidy_anova_results)

res.anovaLogRE.export <- purrr::imap_dfr(res.anova.LogRE, tidy_anova_results)


### Welch's ANOVA for non-equal variance:

# If normality met, but not equal variance, use Welch ANOVA
welch_auto <- function(data, value) {
  # Fit model
  welch_model <- oneway.test(as.formula(paste(value, "~ Condition")), data = data, var.equal = FALSE)
  
  # Check Assumptions
  shapiro_p <- shapiro.test(residuals(lm(as.formula(paste(value, "~ Condition")), data = data)))$p.value # oneway.test() donesn't return a model object with residuals 
  levene_p <- leveneTest(as.formula(paste(value, "~ Condition")), data = data)$`Pr(>F)`[1]
  
  # Perform Post-Hoc Test
  dunnett <- PMCMRplus::dunnettT3Test(as.formula(paste(value, "~ Condition")), data = data)
  
  # Return all results
  list(
    Welch = welch_model,
    ShapiroP = shapiro_p,
    LeveneP = levene_p,
    DunnettT3R = dunnett
  )
}

#### For those normal on row data
dCq_welch.RE <- dCq.outlierRemoved %>%
  filter(ΔCq_label %in% "conditions")

res.welch.RE <- dCq_welch.RE %>%
  group_split(ΔCq_label) %>%
  set_names(dCq_welch.RE %>% group_by(ΔCq_label) %>% group_keys() %>% pull()) %>%
  map(~welch_auto(.x, "RE"))     



#### Tidy results

# Tidy up the results to extract essential info
tidy_welch_results <- function(result, label) {
  # Extract Dunnett's T3 p-values for pairwise comparisons
  dunnett_pvalues <- if (!is.null(result$DunnettT3R) && is.matrix(result$DunnettT3R$p.value)) {
    tmp <- result$DunnettT3R$p.value
    # Convert to a named vector, with comparisons as column names
    comparisons <- outer(rownames(tmp), colnames(tmp), paste, sep = "-")
    comparisons <- comparisons[!is.na(tmp)]  # Filter out NA comparisons
    setNames(as.vector(tmp[!is.na(tmp)]), paste0("DunnettT3P_", comparisons))
  } else {
    NULL
  }
  
  # Combine all results into a single row
  base_row <- tibble::tibble(
    Label = label,
    WelchF = result$Welch$statistic["F"],
    WelchP = result$Welch$p.value,
    ShapiroP = result$ShapiroP,
    LeveneP = result$LeveneP
  )
  # Append Dunnett results if they exist
  base_row <- bind_cols(base_row, as_tibble(as.list(dunnett_pvalues)))
}

# Apply the function
res.welchRE.export <- purrr::imap_dfr(res.welch.RE, tidy_welch_results) 

### Krustal.test on original data
# Extract the data:
dCq_KW <- dCq.OutlierRemoved |>
  filter(ΔCq_label == "conditions")

# Automate Kruskal-Wallis and Post-Hoc Tests
res.kw <- dCq_KW %>%
  group_by(ΔCq_label) %>% 
  summarise(
    # Perform Kruskal-Wallis test
    Kruskal = list(kruskal.test(RE ~ Condition)),
    
    # Extract results from Kruskal-Wallis test
    kwChi_squared = kruskal.test(RE ~ Condition)$statistic,
    kwDF = kruskal.test(RE ~ Condition)$parameter,
    kwP = kruskal.test(RE ~ Condition)$p.value,
    
    # Perform Dunn's post-hoc test if Kruskal is significant, store full result
    Dunn = list(if (kruskal.test(RE ~ Condition)$p.value < 0.05) {
      PMCMRplus::kwAllPairsDunnTest(RE ~ Condition, p.adjust.method = "BH")
    } else {
      NULL
    }),
    
    # Perform pairwise Wilcoxon test if Kruskal is significant, store full result
    Wilcoxon = list(if (kruskal.test(RE ~ Condition)$p.value < 0.05) {
      pairwise.wilcox.test(RE, Condition, p.adjust.method = "BH")
    } else {
      NULL
    })
  )

#### Tidy
# Extract and reshape Dunn's test results to wide format
dunn_results <- res.kw %>%
  filter(!map_lgl(Dunn, is.null)) %>%  # Filter groups with Dunn's test results
  mutate(
    Dunn_Table = map(Dunn, ~ {
      # Extract the p-value matrix
      pval_matrix <- .x$p.value
      
      # Convert to long format
      as.data.frame(as.table(pval_matrix)) %>%
        rename(Group1 = Var1, Group2 = Var2, PValue = Freq) %>%  # Rename columns
        filter(!is.na(PValue))  # Remove NA entries
    })
  ) %>%
  select(ΔCq_label, Dunn_Table) %>%  # Keep relevant columns
  unnest(Dunn_Table) %>%  # Flatten into a single data frame
  mutate(Comparison = paste(Group1, Group2, sep = " vs ")) %>%  # Create a comparison column
  select(ΔCq_label, Comparison, PValue) %>%  # Keep only necessary columns
  pivot_wider(
    names_from = Comparison,  # Use pairwise comparisons as column names
    values_from = PValue      # Fill columns with p-values
  )

wilcoxon_results <- res.kw %>%
  filter(!map_lgl(Wilcoxon, is.null)) %>%  # Filter groups with Wilcoxon test results
  mutate(
    Wilcoxon_Table = map(Wilcoxon, ~ {
      # Extract the p-value matrix
      pval_matrix <- .x$p.value
      
      # Convert to long format
      as.data.frame(as.table(pval_matrix)) %>%
        rename(Group1 = Var1, Group2 = Var2, PValue = Freq) %>%  # Rename columns
        filter(!is.na(PValue))  # Remove NA entries
    })
  ) %>%
  select(ΔCq_label, Wilcoxon_Table) %>%  # Keep relevant columns
  unnest(Wilcoxon_Table)  %>%  # Flatten into a single data frame
  mutate(Comparison = paste(Group1, Group2, sep = " vs ")) %>%  # Create a comparison column
  select(ΔCq_label, Comparison, PValue) %>%  # Keep only necessary columns
  pivot_wider(
    names_from = Comparison,  # Use pairwise comparisons as column names
    values_from = PValue      # Fill columns with p-values
  )

# Comine post hoc results for kw
kw_posthoc <- full_join(dunn_results, wilcoxon_results, by = "ΔCq_label", suffix = c(".dunn", ".wilcoxon"))

# Extract Kruskal-Wallis results
kw_results <- res.kw %>%
  select(ΔCq_label, kwChi_squared, kwDF, kwP)



## Plot selected significant genes
### Compile all significance p value into one single sheet with the same format first (manually now)

# Import compiled p-value file
stats <- read.xlsx("combinedVcStats_results.xlsx")

# Change format to suit stat_pvalue_manual()
sigGenes <- c("CD36", "CDH5", "HIF1A", "PCNA", "PECAM1", "VWF")
gene_pattern <- paste(sigGenes, collapse = "|")

stats <- stats %>%
  filter(str_detect(Label, gene_pattern)) %>% # str_detect() does not inherently handle a vector of patterns (sigGenes) for matching against another vector (label). Instead, it expects a single pattern (or regex) as the second argument
  pivot_longer(
    cols = starts_with("P_"),
    names_to = "comparison", # Temporary column for names
    values_to = "p.adj" # Values go here
  ) %>%
  mutate(
    group1 = str_extract(comparison, "(?<=P_)(D(4|10)\\.[+-]MG)"), # Extract group1, ?<=... look for behind this
    group2 = str_extract(comparison, "(?<=-)(D(4|10)\\.[+-]MG)")  # Extract group2
  ) %>%
  select(Label, model_P, group1, group2, p.adj) # Keep only desired columns

# Add Significant sign
stats <- stats %>%
  mutate(p.adj.signif = case_when(p.adj > 0.1 ~ NA_character_,
                                  p.adj <= 0.1 & p.adj > 0.05 ~ sprintf("%.4f", p.adj),
                                  p.adj <= 0.05 & p.adj > 0.01 ~ "*",
                                  p.adj <= 0.01 & p.adj > 0.001 ~ "**",
                                  p.adj <= 0.001 & p.adj > 0.0001 ~ "***",
                                  p.adj <= 0.0001 ~ "****"),
         ΔCq_label = Label)

# Add Y position value
y_positions <- dCq.outlierRemoved %>%
  filter(Detector %in% sigGenes) %>%
  group_by(ΔCq_label, Condition) %>% # start with the highest among this panel, not only for a Sample
  summarise(max_RE = quantile(RE, 0.75) + IQR(RE)*1.5) %>%
  ungroup() %>%
  group_by(ΔCq_label) %>%
  mutate(y.position = max(max_RE)*1.01) %>% # Add a small value to set the y.position slightly higher
  ungroup() %>%
  dplyr::select(ΔCq_label, y.position) %>%
  distinct()

## Merge y.position with the original data
stats <- stats %>%
  dplyr::filter(!is.na(p.adj.signif)) %>%
  left_join(y_positions, by = "ΔCq_label") %>%
  group_by(ΔCq_label) %>%
  mutate(
    y.position = ifelse(duplicated(y.position) | duplicated(y.position, fromLast = TRUE),
                        y.position + seq_along(y.position) * (0.1*y.position),  # Increment for duplicates by 5% (if seq_along(y.position) * 0.01 is add a fix 0.01 instead)
                        y.position)) %>%
  ungroup() %>%
  dplyr::select(ΔCq_label, model_P, group1, group2, p.adj, p.adj.signif, y.position) 

# Plot with significance annotation
dCq.outlierRemoved |>
  filter(Detector %in% sigGenes) |>
  ggboxplot(x = "Condition", y = "RE", xlab = "Growth States", ylab = "Relative Expression",
            color = "Condition", legend = "right",
            add = "jitter", shape = "Condition",
            facet.by = "ΔCq_label", ncol = 6, scales = "free_y") +
  stat_pvalue_manual(stats, label = "p.adj.signif", label.size = 3, tip.length = 0.01, step.increase = 0.02, step.group.by = "ΔCq_label") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

