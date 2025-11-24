#---------#
#Data Processing & Analysis for Figure 4 using two-way ANOVA
#author: "Shiqi Huang"
#---------#

# Load necessary libraries
library(openxlsx) 
library(tidyverse)
library(purrr)
library(ggpubr) # https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/ 
library(rstatix) # https://rpkgs.datanovia.com/rstatix/ 
library(outliers)
library(car)   # For Levene's test
library(DescTools)  # For post-hoc tests
library(PMCMRplus) # For post-hoc tests
library(broom) # tidy stats objects into tibbles
library(emmeans) # 2-way ANOVA post-hoc
library(nlme) # gls for non-equal variance ANOVA
library(ARTool) # art for non-parametric ANOVA (Rank based)
library(WRS2) # t2way for robust 2-way ANOVA (Handling outliers and unequal variances without assuming normality.)

## Import data
## Create a ΔCt data frame as in One-way_ANOVA.R, but with additional HKG combinations:
HKG_combinations <- list(
  sHKG3 = c("CAPZB", "FBXO7", "YWHAZ"),
  cHKG3 = c("18S", "ACTB", "HPRT1"),
  `s2_18S` = c("18S", "CAPZB", "FBXO7"),
  sHKG2 = c("CAPZB", "FBXO7"),
  cHKG2 = c("18S", "HPRT1"),
  `Cap_18S` = c("18S", "CAPZB"),
  CAPZB = c("CAPZB"),
  `18S` = c("18S")
)

# Assess Assumptions of ANOVA

## 1 Check for Outlier
### auto function:
detect_outliers <- function(data, group_var, response_var, factor1, factor2 = NULL, global_row_col = "global_row") {
  # Check required packages
  require(dplyr, quietly = TRUE)
  require(tidyr, quietly = TRUE)
  require(rstatix, quietly = TRUE)
  require(broom, quietly = TRUE)
  require(car, quietly = TRUE)
  
  tryCatch({
    if (!global_row_col %in% colnames(data)) {
      data <- data %>% 
        mutate(!!global_row_col := row_number())
    }
    
    data %>%
      group_by(across({{ group_var }})) %>%
      group_modify(~ {
        current_data <- .x
        
        # Convert factors to factor type
        factors_to_convert <- factor1
        if (!is.null(factor2)) factors_to_convert <- c(factors_to_convert, factor2)
        
        current_data <- current_data %>%
          dplyr::mutate(dplyr::across(tidyselect::all_of(factors_to_convert), ~ factor(.x)))
        
        # Determine IQR grouping variables
        iqr_group_vars <- factor1
        if (!is.null(factor2)) iqr_group_vars <- c(iqr_group_vars, factor2)
        
        # Build model formula
        if (is.null(factor2)) {
          model_formula <- stats::as.formula(paste(response_var, "~", factor1))
        } else {
          model_formula <- stats::as.formula(paste(response_var, "~", factor1, "*", factor2))
        }
        
        # 1. IQR-based outlier detection
        iqr_outliers <- tryCatch({
          current_data %>%
            dplyr::group_by(dplyr::across(tidyselect::all_of(iqr_group_vars))) %>%
            rstatix::identify_outliers(!!rlang::sym(response_var)) %>%
            dplyr::mutate(
              is.outlier_iqr = is.outlier,
              is.extreme_iqr = is.extreme
            ) %>%
            dplyr::select(!!rlang::sym(global_row_col), is.outlier_iqr, is.extreme_iqr)
        }, error = function(e) {
          message("IQR detection failed: ", e$message)
          tibble::tibble(!!rlang::sym(global_row_col) := integer(), 
                         is.outlier_iqr = logical(), 
                         is.extreme_iqr = logical())
        })
        
        # 2. Model-based outlier detection
        model <- tryCatch(
          stats::lm(model_formula, data = current_data),
          error = function(e) NULL
        )
        
        if (is.null(model)) {
          bonf_p <- tibble::tibble(!!rlang::sym(global_row_col) := numeric(), 
                                   bonferroni_p = numeric())
          cooks_d <- tibble::tibble(!!rlang::sym(global_row_col) := numeric(), 
                                    cooks_d = numeric())
        } else {
          # Bonferroni p-values
          outlier_test <- tryCatch(
            car::outlierTest(model),
            error = function(e) NULL
          )
          
          bonf_p <- if (!is.null(outlier_test)) {
            tibble::tibble(
              !!rlang::sym(global_row_col) := current_data[[global_row_col]][as.numeric(names(outlier_test$rstudent))],
              bonferroni_p = as.numeric(outlier_test$bonf.p)
            )
          } else {
            tibble::tibble(!!rlang::sym(global_row_col) := numeric(), 
                           bonferroni_p = numeric())
          }
          
          # Cook's distance
          cooks_d <- tryCatch({
            broom::augment(model) %>%
              dplyr::mutate(!!rlang::sym(global_row_col) := current_data[[global_row_col]]) %>%
              dplyr::select(!!rlang::sym(global_row_col), cooks_d = .cooksd)
          }, error = function(e) {
            message("Cook's distance failed: ", e$message)
            tibble::tibble(!!rlang::sym(global_row_col) := numeric(), 
                           cooks_d = numeric())
          })
        }
        
        # Combine all results
        current_data %>%
          dplyr::left_join(iqr_outliers, by = global_row_col) %>%
          dplyr::left_join(bonf_p, by = global_row_col) %>%
          dplyr::left_join(cooks_d, by = global_row_col) %>%
          dplyr::mutate(
            is.outlier_iqr = tidyr::replace_na(is.outlier_iqr, FALSE),
            is.extreme_iqr = tidyr::replace_na(is.extreme_iqr, FALSE),
            bonferroni_p = tidyr::replace_na(bonferroni_p, 1)
          )
        
      }, error = function(e) {
        message("Outlier detection failed: ", e$message)
        return(data %>% dplyr::mutate(error = paste("Error:", e$message)))
      })
  })
}

### Detect
ddCq.outlier.RE <- detect_outliers(
  data = ddCq,
  group_var = "ΔCq_label",
  response_var = "RE",
  factor1 = "Condition",
  factor2 = "Sample")

### View
outlier.RE <- ddCq.outlier.RE %>%
  filter(is.outlier_iqr | is.extreme_iqr | (bonferroni_p <= 0.05) | (cooks_d >= 0.5)) %>%
  select(ΔCq_label, Condition, Sample, RE, 
         is.outlier_iqr, is.extreme_iqr, bonferroni_p, cooks_d) # 267

### Remove outliers: 
#### Tier 1: extreme : ≥ 3xIQR or cooks ≥ 1:
outlier.RE %>% dplyr::filter(is.extreme_iqr | cooks_d >= 1) %>% print(n = Inf) # 0
#### Tier 2: bonf_p ≤ 0.05 && cooks ≥ 0.5:
outlier.RE %>% dplyr::filter(bonferroni_p <= 0.05 & cooks_d >= 0.5) %>% print(n = Inf) # 26
# Remove those but ensure sample size per group n≥3

## 2. Test Normality & Equal V
### General helper function
get_anova_p <- function(model) {
  tidy_df <- broom::tidy(model) %>%
    filter(term != "Residuals")
  p_values <- tidy_df$p.value
  names(p_values) <- tidy_df$term
  p_values
}
### Test
res.assumption.RE <- ddCq.REoutlierRemoved %>%
  group_by(ΔCq_label) %>%
  summarise(
    # Perform ANOVA and Store the named vector of p-values in a list-column
    ANOVA_p = list(get_anova_p(aov(RE ~ Condition * Sample))),
    # Check Assumptions
    ShapiroP = list(shapiro.test(residuals(aov(RE ~ Condition * Sample)))$p.value),
    LeveneP = list(leveneTest(RE ~ Condition * Sample)$`Pr(>F)`[1])
  ) %>%
  # Spread the list-column into separate columns.
  unnest_wider(ANOVA_p, names_sep = "_") %>%
  unnest(cols = c(ShapiroP, LeveneP))

### Highlight significant results. ie. not normal and/or not equal variance
res.assumption.RE |>
  dplyr::filter(ShapiroP < 0.05) |>
  dplyr::select(ΔCq_label, ShapiroP, LeveneP) %>%
  print(n = Inf) # n=73

res.assumption.RE |>
  dplyr::filter(ShapiroP < 0.01) |>
  dplyr::select(ΔCq_label, ShapiroP, LeveneP) %>%
  print(n = Inf) # n=69

res.assumption.RE |>
  dplyr::filter(LeveneP < 0.05) |>
  dplyr::select(ΔCq_label, ShapiroP, LeveneP) %>%
  print(n = Inf) # n=34 

## 3. Transform the data
### Apply log transformation
dCq.transformed.RE <- ddCq %>%
  filter(Detector %in% "conditions based on res.assumption.RE") %>%
  mutate(LogRE = log10(RE)) # 3264

### Re-test Outlier
dCqTrans.outlier.RE <- detect_outliers(
  data = dCq.transformed.RE,
  group_var = "ΔCq_label",
  response_var = "LogRE",
  factor1 = "Condition",
  factor2 = "Sample")

### Repeat the View and Remove process accordingly

### Re-test Normality & euqal variance
res.assumption.RE.1 <- dCq.transformed.RE %>%
  group_by(ΔCq_label) %>%
  summarise(
    # Perform ANOVA and Store the named vector of p-values in a list-column
    ANOVA_p = list(get_anova_p(aov(LogRE ~ Condition * Sample))),
    # Check Assumptions
    ShapiroP = list(shapiro.test(residuals(aov(LogRE ~ Condition * Sample)))$p.value),
    LeveneP = list(leveneTest(LogRE ~ Condition * Sample)$`Pr(>F)`[1])
  ) %>%
  # Spread the list-column into separate columns.
  unnest_wider(ANOVA_p, names_sep = "_") %>%
  unnest(cols = c(ShapiroP, LeveneP))

### Highlight significant results. ie. not normal and/or not equal variance
res.assumption.RE.1 |>
  dplyr::filter(ShapiroP < 0.05) |>
  dplyr::select(ΔCq_label, ShapiroP, LeveneP) %>%
  print(n = Inf) # n=10

res.assumption.RE.1 |>
  dplyr::filter(ShapiroP < 0.01) |>
  dplyr::select(ΔCq_label, ShapiroP, LeveneP) %>%
  print(n = Inf) # n=3

res.assumption.RE.1 |>
  dplyr::filter(LeveneP < 0.05) |>
  dplyr::select(ΔCq_label, ShapiroP, LeveneP) %>%
  print(n = Inf) # n=2 

## 4. Run ANOVA
### 2-way ANOVA
#### wrapper functions
##### Function to perform 2-way ANOVA + posthoc analysis 
anova_auto_2way <- function(data, value, factor1, factor2, interaction = TRUE, ss_type = 3, white_adjust = FALSE, adjust = "tukey") {
  # Convert to factors with safe conversion
  data <- data %>%
    dplyr::mutate(
      !!factor1 := factor(!!sym(factor1)),
      !!factor2 := factor(!!sym(factor2))
    )
  
  # Build formula safely
  formula_str <- if(interaction) {
    paste(value, "~", factor1, "*", factor2)
  } else {
    paste(value, "~", factor1, "+", factor2)
  }
  
  # Fit model with error handling
  model_lm <- try(lm(as.formula(formula_str), data = data))
  if(inherits(model_lm, "try-error")) return("Model fitting failed")
  
  # ANOVA analysis
  anova_result <- try(car::Anova(model_lm, type = ss_type, white.adjust = white_adjust))
  if(inherits(anova_result, "try-error")) return("ANOVA failed")
  
  # Post-hoc analysis
  posthoc <- list()
  
  # Main effects
  for(f in c(factor1, factor2)) {
    emm <- emmeans(model_lm, specs = f)
    posthoc[[paste0("main_", f)]] <- pairs(emm, adjust = adjust)
  }
  
  # Interaction decomposition
  if(interaction) {
    # Get full interaction means (using * instead of :)
    emm_int <- emmeans::emmeans(model_lm, 
                                specs = reformulate(paste(factor1, factor2, sep = "*")))
    
    # Within-factor comparisons with dynamic names
    posthoc[[paste0("Within_", factor1)]] <- emmeans::contrast(
      emm_int,
      method = "pairwise",
      by = factor1,
      adjust = adjust
    )
    
    posthoc[[paste0("Within_", factor2)]] <- emmeans::contrast(
      emm_int,
      method = "pairwise",
      by = factor2,
      adjust = adjust
    )
  }
  
  list(
    anova = anova_result,
    posthoc = posthoc,
    assumptions = list(
      shapiro = shapiro.test(residuals(model_lm)),
      levene = car::leveneTest(model_lm)
    ),
    factors = list(f1 = factor1, f2 = factor2)  # Store factor names for tidy function
  )
}

##### Tidy result to wide data fame with function
tidy_anova_export <- function(anova_list, group_name) {
  # Initialize base tibble with group name
  base_tibble <- tibble::tibble(Group = group_name)
  
  # Early error handling
  if (!is.list(anova_list)) return(tibble::tibble(Group = group_name, Error = "Invalid result"))
  if (!"factors" %in% names(anova_list)) return(tibble::tibble(Group = group_name, Error = "No factor info"))
  
  # Extract factors
  factor1 <- anova_list$factors$f1
  factor2 <- anova_list$factors$f2
  
  # Enhanced extraction function for different effect types
  extract_ph <- function(ph_obj, effect_type, factor_name) {
    if (is.null(ph_obj)) return(tibble::tibble())
    
    tryCatch({
      ph_df <- ph_obj %>% 
        broom::tidy() %>%
        dplyr::mutate(contrast = gsub(" - ", "vs", contrast) %>% gsub("[()]", "", .))
      
      # Get factor levels directly from emmeans object
      num_levels <- tryCatch({
        if (inherits(ph_obj, "emmGrid")) {
          length(ph_obj@levels$contrast)
        } else {
          NA
        }
      }, error = function(e) NA)
      
      # Select appropriate p-value column
      p_value_col <- if (is.na(num_levels)) {
        # Fallback to column existence check
        if ("adj.p.value" %in% names(ph_df)) "adj.p.value" else "p.value"
      } else {
        if (num_levels > 2) "adj.p.value" else "p.value"
      }
      
      # Handle interaction layers first
      if (any(c(factor1, factor2) %in% names(ph_df))) {
        by_var <- intersect(names(ph_df), c(factor1, factor2))
        ph_df <- ph_df %>%
          tidyr::unite(col = "contrast",
                       contrast, !!by_var,
                       sep = "@")  # Unique separator
      }
      
      # Choose p_value_col (if this step before interaction one, will cause error)
      ph_df <- ph_df %>%
        dplyr::select(contrast, !!rlang::sym(p_value_col)) %>%
        dplyr::rename(p_value = !!rlang::sym(p_value_col))
      
      ph_df %>%
        dplyr::select(contrast, p_value) %>%
        tidyr::pivot_wider(
          names_from = contrast,
          values_from = p_value,
          names_prefix = paste0(effect_type, "_", factor_name, "_")
        )
    }, error = function(e) {
      tibble::tibble(!!paste0(effect_type, "_", factor_name, "_error") := "Posthoc failed")
    })
  }
  
  # ANOVA table processing
  anova_p <- tryCatch({
    anova_data <- anova_list$anova
    
    # If it's a car::Anova object or other supported objects:
    if (inherits(anova_data, "Anova.mlm") || inherits(anova_data, "anova")) {
      broom::tidy(anova_data) %>%
        dplyr::filter(term != "Residuals") %>%
        dplyr::select(term, p.value) %>%
        tidyr::pivot_wider(
          names_from = term,
          values_from = p.value,
          names_prefix = "ANOVA_"
        )
    } else if (inherits(anova_data, "data.frame")) {
      # For GLS anova output: treat it directly as data.frame
      anova_data %>%
        tibble::rownames_to_column("term") %>%
        dplyr::filter(term != "Residuals") %>%
        dplyr::select(term, `p-value`) %>%
        tidyr::pivot_wider(
          names_from = term,
          values_from = `p-value`,
          names_prefix = "ANOVA_"
        )
    } else {
      tibble::tibble(ANOVA_error = "Unsupported anova object")
    }
  }, error = function(e) {
    tibble::tibble(ANOVA_error = "ANOVA processing failed")
  })
  
  # Process all effects
  main1 <- extract_ph(anova_list$posthoc[[paste0("main_", factor1)]], "Main", factor1)
  main2 <- extract_ph(anova_list$posthoc[[paste0("main_", factor2)]], "Main", factor2)
  int1 <- extract_ph(anova_list$posthoc[[paste0("Within_", factor1)]], "Int", factor1)
  int2 <- extract_ph(anova_list$posthoc[[paste0("Within_", factor2)]], "Int", factor2)
  
  # Combine components
  dplyr::bind_cols(
    base_tibble,
    anova_p,
    main1,
    main2,
    int1,
    int2
  )
}

#### Apply to those normal (or mild unnormal with ShapiroP >0.01) and equal V
dCq_anovaRE <- ddCq.REoutlierRemoved %>%
  filter(ΔCq_label %in% "conditions match")

res.anova.RE <- dCq_anovaRE %>%
  group_split(ΔCq_label) %>%
  set_names(dCq_anovaRE %>% group_by(ΔCq_label) %>% group_keys() %>% pull()) %>%
  imap(~ tryCatch(
    anova_auto_2way(.x, "RE", "Condition", "Sample", interaction = TRUE), # use LogRE for those fit with transformed
    error = function(e) {
      # Return a custom error message that includes the group name (provided by .y)
      paste("Error in group", .y, ":", e$message)
    }
  )
  )

res.anova.RE.export <- purrr::imap_dfr(res.anova.RE, ~ {
  if (is.character(.x)) {
    tibble::tibble(Group = .y, Error = .x)
  } else {
    tidy_anova_export(.x, .y)
  }
}) %>%
  dplyr::select(Group, dplyr::everything())

#### Apply to those normal (or mild unnormal with ShapiroP ≥ 0.01) but mild unequal V (0.05 > LeveneP ≥ 0.01): with white.adjust
dCq_anovaRE.white <- ddCq.REoutlierRemoved %>%
  dplyr::filter(ΔCq_label %in% "conditions match")

res.anova_white.RE <- dCq_anovaRE.white %>%
  group_split(ΔCq_label) %>%
  set_names(dCq_anovaRE.white %>% group_by(ΔCq_label) %>% group_keys() %>% pull()) %>%
  imap(~ tryCatch(
    anova_auto_2way(.x, "RE", "Condition", "Sample", interaction = TRUE, white_adjust = TRUE), # use LogRE for those fit with transformed
    error = function(e) {
      # Return a custom error message that includes the group name (provided by .y)
      paste("Error in group", .y, ":", e$message)
    }
  )
  )

res.anova_white.RE.export <- purrr::imap_dfr(res.anova_white.RE, ~ {
  if (is.character(.x)) {
    tibble::tibble(Group = .y, Error = .x)
  } else {
    tidy_anova_export(.x, .y)
  }
}) %>%
  dplyr::select(Group, dplyr::everything())

### GLS model for those normal but severe unequal V
#### wrapper functions
gls_model_selector <- function(df, response, factor1, factor2, verbose = TRUE) {
  require(nlme)
  require(MuMIn)
  
  # Input validation
  if (!all(c(factor1, factor2, response) %in% names(df))) {
    stop("One or more variables not found in dataframe")
  }
  
  # Observation check
  obs_check <- table(df[[factor1]], df[[factor2]])
  if(any(obs_check < 2)) {
    stop("Insufficient observations for variance estimation")
  }
  
  # Create interaction term as a new column
  interaction_name <- "interaction_term"  # Explicitly define the column name
  df[[interaction_name]] <- interaction(df[[factor1]], df[[factor2]], drop = TRUE)
  df[[interaction_name]] <- as.factor(df[[interaction_name]])  # Use [[ with string
  
  ## Create the main model formula
  form <- stats::as.formula(paste(response, "~", factor1, "*", factor2))
  
  ## Variance structure formulas
  var_form_factor1 <- stats::as.formula(paste("~1 |", factor1))
  var_form_factor2 <- stats::as.formula(paste("~1 |", factor2))
  var_form_interaction <- stats::as.formula(paste("~1 |", interaction_name))
  
  ## Fit models
  mod_homogeneous <- nlme::gls(form, data = df)
  mod_factor1 <- nlme::gls(form, data = df, 
                           weights = nlme::varIdent(form = var_form_factor1))
  mod_factor2 <- nlme::gls(form, data = df, 
                           weights = nlme::varIdent(form = var_form_factor2))
  mod_interaction <- tryCatch(
    nlme::gls(form, data = df, 
              weights = nlme::varIdent(form = var_form_interaction)),
    error = function(e) {
      message("Interaction model error: ", e$message)
      NULL
    }  # Gracefully handle failed interaction model
  )
  
  model_list <- list(homogeneous = mod_homogeneous, 
                     factor1 = mod_factor1, 
                     factor2 = mod_factor2, 
                     interaction = mod_interaction)
  # Remove NULL entries
  model_list <- Filter(Negate(is.null), model_list)  
  # Check if any models remain
  if (length(model_list) == 0) {
    stop("All models failed to fit.")
  }
  # Model selection
  aicc <- sapply(model_list, MuMIn::AICc)
  best <- model_list[[which.min(aicc)]]
  
  if(verbose) {
    cat("Best model:", names(which.min(aicc)), "\n")
    print(aicc)
    print(stats::anova(best))
  }
  
  # Diagnostic plots (if verbose)
  if (verbose) {
    par(mfrow = c(1, 2))
    plot(fitted(best), resid(best, type = "pearson"),
         main = "Residuals vs Fitted", xlab = "Fitted values", ylab = "Residuals")
    abline(h = 0, col = "red")
    qqnorm(resid(best), main = "Q-Q Plot")
    qqline(resid(best), col = "red")
  }
}
#**NOTE** emmeans() with gls models often fails with dynamic specs for interactions, so build gls model and emm_int outside master functions, then pass it in

gls_posthoc <- function(model, emm_int, df, factor1, factor2, adjust = "tukey") {
  require(emmeans)
  # ANOVA
  anova_res <- stats::anova(model)
  
  # Post-hoc analysis
  posthoc <- list() 
  
  # Pairwise comparisons
  ## Main effects
  for(f in c(factor1, factor2)) {
    emm <- try(emmeans::emmeans(model, specs = f, data = df), silent = TRUE)
    if(!inherits(emm, "try-error")) {
      posthoc[[paste0("main_", f)]] <- pairs(emm, adjust = adjust)
    }
  }
  
  ## Within-factor comparisons with dynamic names
  posthoc[[paste0("Within_", factor1)]] <- emmeans::contrast(
    emm_int,
    method = "pairwise",
    by = factor1,
    adjust = adjust
  )
  
  posthoc[[paste0("Within_", factor2)]] <- emmeans::contrast(
    emm_int,
    method = "pairwise",
    by = factor2,
    adjust = adjust
  )
  
  resd <- residuals(model)
  return(list(
    anova = anova_res,
    posthoc = posthoc,
    assumptions = list(
      shapiro = shapiro.test(resd),
      levene = car::leveneTest(resd ~ df[[factor1]] * df[[factor2]])
    ),
    factors = list(f1 = factor1, f2 = factor2)
  ))
}

#### Apply to data (in our case, on transformed data only)
dCq_gls <- dCq.transformed.RE %>%
  dplyr::filter(ΔCq_label == "condition match")
##### Test model
gls_model_selector(dCq_glsLogRE, "LogRE", "Condition", "Sample") #best model is factor1, so "Condition"

##### Apply model
res.glsLogRE <- dCq_glsLogRE %>%
  group_split(ΔCq_label) %>%
  set_names(dCq_glsLogRE %>% group_by(ΔCq_label) %>% group_keys() %>% pull()) %>%
  imap(~ tryCatch({
    best_model <- nlme::gls(LogRE ~ Condition * Sample, data = .x, weights = varIdent(form = ~1|Condition))
    emm_int <- emmeans::emmeans(best_model, specs = ~ Condition * Sample, data = .x) #gls model don't store data argument cleanly as in lm model, so must specify if call emmeans inside another function
    gls_posthoc(best_model, emm_int, df = .x, factor1 = "Condition", factor2 = "Sample")
  }, error = function(e) {
    paste("Error in group", .y, ":", e$message)
  }))

##### Tidy
res.glsLogRE.export <- purrr::imap_dfr(res.glsLogRE, ~ {
  if (is.character(.x)) {
    tibble::tibble(Group = .y, Error = .x)
  } else {
    tidy_anova_export(.x, .y)  # this function will now handle gls-based results
  }
}) %>%
  dplyr::select(Group, dplyr::everything())

### non-parametric
#### ART
art_auto_2way <- function(data, value, factor1, factor2, adjust = "tukey", inter = FALSE) {
  require(ARTool)
  require(emmeans)
  
  data <- data %>%
    dplyr::mutate(
      !!factor1 := factor(!!sym(factor1)),
      !!factor2 := factor(!!sym(factor2))
    )
  
  formula_str <- as.formula(paste(value, "~", factor1, "*", factor2))
  
  # Fit ART model
  art_model <- ARTool::art(formula_str, data = data)
  
  # ANOVA table
  anova_res <- anova(art_model)
  
  # Post-hoc
  posthoc <- list()
  
  for(f in c(factor1, factor2)) {
    posthoc[[paste0("main_", f)]] <- art.con(art_model, f, adjust = adjust)
  }
  
  posthoc[["interaction"]] <- art.con(art_model, 
                                      paste(factor1, factor2, sep = ":"), 
                                      adjust = adjust,
                                      interaction = inter) # For interaction contrasts
  
  list(
    anova = anova_res,
    posthoc = posthoc,
    assumptions = list(
      shapiro = shapiro.test(art_model$residuals),
      levene = car::leveneTest(art_model$residuals ~ data[[factor1]] * data[[factor2]])
    ),
    factors = list(f1 = factor1, f2 = factor2)
  )
}
#**Note** this post-hoc has all possible pairwise comparisons (39 more), inflated number of tests and may affect type I error rate

tidy_art_export <- function(art_list, group_name) {
  # Initialize base tibble with group name
  base_tibble <- tibble::tibble(Group = group_name)
  
  # Early error handling
  if (!is.list(art_list)) return(tibble::tibble(Group = group_name, Error = "Invalid result"))
  if (!"factors" %in% names(art_list)) return(tibble::tibble(Group = group_name, Error = "No factor info"))
  
  # Extract factors
  factor1 <- art_list$factors$f1
  factor2 <- art_list$factors$f2
  
  # Custom ART ANOVA processor
  process_art_anova <- function(anova_data) {
    tryCatch({
      anova_data %>%
        tibble::rownames_to_column("term") %>%
        dplyr::filter(term != "Residuals") %>%
        dplyr::select(
          term, 
          p.value = dplyr::any_of(c("Pr(>F)", "p.value"))
        ) %>%
        tidyr::pivot_wider(
          names_from = term,
          values_from = p.value,
          names_prefix = "ANOVA_"
        )
    }, error = function(e) {
      tibble::tibble(ANOVA_error = "ANOVA processing failed")
    })
  }
  
  # Generic main effect processor (updated)
  process_main_effect <- function(ph_obj, effect_type, factor_name) {
    if (is.null(ph_obj)) return(tibble::tibble())
    
    tryCatch({
      ph_obj %>%
        as.data.frame()  %>% #broom::tidy() don't support art objects
        dplyr::mutate(
          contrast = contrast %>%
            as.character() %>%
            gsub(" - ", "vs", .) %>%
            gsub("[()]", "", .)
        ) %>%
        dplyr::select(contrast, p.value) %>%
        tidyr::pivot_wider(
          names_from = contrast,
          values_from = p.value,
          names_prefix = paste0(effect_type, "_", factor_name, "_")
        )
    }, error = function(e) {
      tibble::tibble(!!paste0(effect_type, "_", factor_name, "_error") := "Posthoc processing failed")
    })
  }
  
  # Enhanced interaction processor
  process_art_interaction <- function(ph_obj) {
    if (is.null(ph_obj)) return(tibble::tibble())
    
    tryCatch({
      ph_df <- ph_obj %>%
        as.data.frame()  %>%
        dplyr::mutate(contrast = contrast %>%
                        as.character() %>%
                        gsub(" - ", "vs", .) %>%  # Replace contrasts
                        gsub("[()]", "", .))      # Remove parentheses
      
      # Flexible splitting using regular expressions
      processed <- ph_df %>%
        tidyr::extract(
          col = contrast,
          into = c("f1_left", "f2_left", "f1_right", "f2_right"),
          regex = "^([^,]+),([^vs]+)vs([^,]+),(.+)$",
          remove = FALSE
        )
      
      # Dynamic contrast identification
      processed %>%
        dplyr::mutate(
          # Determine which factor is constant
          constant_factor = dplyr::case_when(
            f1_left == f1_right ~ factor1,
            f2_left == f2_right ~ factor2,
            TRUE ~ "mixed"
          ),
          contrast_type = dplyr::case_when(
            constant_factor == factor1 ~ paste0("Int_", factor1),
            constant_factor == factor2 ~ paste0("Int_", factor2),
            TRUE ~ "Complex"
          ),
          contrast_values = dplyr::case_when(
            constant_factor == factor1 ~ paste0(f2_left, "vs", f2_right, "@", f1_left),
            constant_factor == factor2 ~ paste0(f1_left, "vs", f1_right, "@", f2_left),
            TRUE ~ paste0(f1_left, ",", f2_left, "vs", f1_right, ",", f2_right)
          )
        ) %>%
        dplyr::filter(constant_factor != "mixed") %>%
        dplyr::mutate(
          full_name = paste0(contrast_type, "_", contrast_values)
        ) %>%
        dplyr::distinct(full_name, p.value) %>%
        tidyr::pivot_wider(
          names_from = full_name,
          values_from = p.value
        )
    }, error = function(e) {
      tibble::tibble(Int_error = "Interaction processing failed")
    })
  }
  
  # Process components
  anova_p <- process_art_anova(art_list$anova)
  main1 <- process_main_effect(art_list$posthoc[[paste0("main_", factor1)]], "Main", factor1)
  main2 <- process_main_effect(art_list$posthoc[[paste0("main_", factor2)]], "Main", factor2)
  int_effects <- process_art_interaction(art_list$posthoc$interaction)
  
  # Combine all components
  dplyr::bind_cols(
    base_tibble,
    anova_p,
    main1,
    main2,
    int_effects
  )
}

#### WRS2
wrs2_auto_2way <- function(data, value, factor1, factor2, tr = 0.2) {
  require(WRS2)
  tryCatch({
    data <- data %>%
      dplyr::mutate(
        !!factor1 := factor(!!sym(factor1)),
        !!factor2 := factor(!!sym(factor2))
      )
    
    # Generate all pairwise comparisons
    f1_levels <- levels(data[[factor1]])
    f2_levels <- levels(data[[factor2]])
    
    formula_str <- as.formula(paste(value, "~", factor1, "*", factor2))
    
    # Robust two-way ANOVA
    ## t2way(): a simple and efficient robust ANOVA method for data with potential outliers using trimmed means
    ## pbad2way(): a more robust approach that can handle a wider range of potential deviations from normality, use M-estimators and medians, post-hoc with mcp2a()
    anova_res <- WRS2::t2way(formula_str, data = data, tr = tr)
    
    # Official post-hoc using mcp2atm: 
    posthoc <- tryCatch({WRS2::mcp2atm(formula_str, data = data, tr = tr) # difference of difference, can't accept contast matrix argument
    }, error = function(e) {
      # Check for DF error specifically
      if (grepl("degrees of freedom must be greater than or equal to 2", e$message)) {
        warning("\n! DEGREES OF FREEDOM WARNING!",
                "\nTrimming level (tr = ", tr, ") may be too aggressive",
                "\nOriginal error: ", e$message, "\n")
      }
      return(conditionMessage(e))  # Return error message as character
    })
    
    list(
      anova = anova_res,
      posthoc = posthoc,
      factors = list(f1 = factor1, f2 = factor2)
    )
  }, error = function(e) {
    return(conditionMessage(e))  # Return other errors as character
  })
}

tidy_wrs2_export <- function(wrs2_list, group_name) {
  base_tibble <- tibble::tibble(Group = group_name)
  
  if (!is.list(wrs2_list)) return(tibble::tibble(Group = group_name, Error = "Invalid result"))
  if (!"factors" %in% names(wrs2_list)) return(tibble::tibble(Group = group_name, Error = "No factor info"))
  
  factor1 <- wrs2_list$factors$f1
  factor2 <- wrs2_list$factors$f2
  
  #Extract ANOVA results
  anova_p <- tibble::tibble(
    "ANOVA_{factor1}" := wrs2_list$anova$A.p.value,
    "ANOVA_{factor2}" := wrs2_list$anova$B.p.value,
    "ANOVA_{factor1}:{factor2}" := wrs2_list$anova$AB.p.value
  )
  
  #Post hoc results
  if ("posthoc" %in% names(wrs2_list) && 
      is.list(wrs2_list$posthoc) &&
      "effects" %in% names(wrs2_list$posthoc) &&
      "contrasts" %in% names(wrs2_list$posthoc)) {
    ph <- wrs2_list$posthoc
    
    # Helper to parse contrast names
    parse_contrast <- function(col_name, contrast_df) {
      # Get factor components from row names
      row_parts <- strsplit(rownames(contrast_df), "_")
      f1_levels <- sapply(row_parts, `[`, 1)
      f2_levels <- sapply(row_parts, `[`, 2)
      
      # Find compared groups
      pos <- which(contrast_df[[col_name]] == 1)
      neg <- which(contrast_df[[col_name]] == -1)
      
      if (length(pos) == 0 || length(neg) == 0) return(NA_character_)
      
      # Main effect comparisons
      if (grepl(paste0("^", factor1, "[0-9]+$"), col_name)) {
        lvls <- unique(f1_levels[c(pos, neg)])
        return(paste0("Main_", factor1, "_", lvls[1], "vs", lvls[2]))
      }
      
      if (grepl(paste0("^", factor2, "[0-9]+$"), col_name)) {
        lvls <- unique(f2_levels[c(pos, neg)])
        return(paste0("Main_", factor2, "_", lvls[1], "vs", lvls[2]))
      }
      
      # Interaction comparisons
      if (grepl(paste0(factor1, "[0-9]+:", factor2, "[0-9]+"), col_name)) {
        parts <- strsplit(col_name, ":")[[1]]
        f1_comp <- parse_contrast(parts[1], contrast_df) %>% gsub("^Main_", "", .)
        f2_comp <- parse_contrast(parts[2], contrast_df) %>% gsub("^Main_", "", .)
        return(paste0("Interaction_", f1_comp, "_x_", f2_comp))
      }
      
      return(NA_character_)
    }
    
    # Process posthoc effects
    ph_results <- purrr::map(names(ph$effects), function(effect_type) {
      pvals <- ph$effects[[effect_type]]$p.value
      # Dynamic regex based on effect type
      pattern <- case_when(
        effect_type == factor1 ~ paste0("^", factor1, "[0-9]+$"),
        effect_type == factor2 ~ paste0("^", factor2, "[0-9]+$"),
        effect_type == paste0(factor1, ":", factor2) ~ paste0("^", factor1, "[0-9]+:", factor2, "[0-9]+"),
        TRUE ~ "^$" # No match
      )
      
      contrast_cols <- grep(pattern, names(ph$contrasts), value = TRUE)
      comp_names <- purrr::map_chr(contrast_cols, ~parse_contrast(.x, ph$contrasts))
      
      # Ensure equal lengths before setting names
      if (length(pvals) != length(comp_names)) {
        stop(glue::glue("Mismatch in {effect_type}: {length(pvals)} p-values vs {length(comp_names)} names"))
      }
      
      stats::setNames(as.list(pvals), comp_names)
    }) |> purrr::flatten()
    
    ph_results <- tibble::as_tibble(ph_results)
  } else {
    ph_results <- tibble::tibble(Posthoc_Error = "Invalid/missing posthoc structure")
  }
  
  # Combine all results
  final_df <- dplyr::bind_cols(
    base_tibble,
    anova_p,
    ph_results
  )
  
  return(final_df)
}


####Application
dCq_nonparamRE <- ddCq.REoutlierRemoved %>%
  dplyr::filter(ΔCq_label %in% "conditions match")

##### ART model
res.artRE <- dCq_nonparamRE %>%
  group_split(ΔCq_label) %>%
  set_names(dCq_nonparamRE %>% group_by(ΔCq_label) %>% group_keys() %>% pull()) %>%
  imap(~ tryCatch({
    art_auto_2way(.x, value = "RE", factor1 = "Condition", factor2 = "Sample")
  }, error = function(e) {
    paste("Error in group", .y, ":", e$message)
  })) 

res.artRE.export <- purrr::imap_dfr(res.artRE, ~ {
  if (is.character(.x)) {
    tibble::tibble(Group = .y, Error = .x)
  } else {
    tidy_art_export(.x, .y)  # this function will now handle gls-based results
  }
}) %>%
  dplyr::select(Group, tidyselect::everything())

##### ART model
res.wrsRE <- dCq_nonparamRE %>%
  group_split(ΔCq_label) %>%
  set_names(dCq_nonparamRE %>% group_by(ΔCq_label) %>% group_keys() %>% pull()) %>%
  imap(~ tryCatch({
    wrs2_auto_2way(.x, value = "RE", factor1 = "Condition", factor2 = "Sample", tr = 0.2)
  }, error = function(e) {
    paste("Error in group", .y, ":", e$message)
  }))

res.wrsRE.export <- purrr::imap_dfr(res.wrsRE, ~ {
  if (is.character(.x)) {
    tibble::tibble(Group = .y, Error = .x)
  } else {
    tidy_wrs2_export(.x, .y)  # this function will now handle gls-based results
  }
}) %>%
  dplyr::select(Group, tidyselect::everything())

## Export
### List of data frames with sheet names
res.list <- list("ANOVA_RE" = res.anova.RE.export, 
                 "ANOVA_logRE" = res.anovaLogRE.export,
                 "ANOVA_white_RE" = res.anova_white.RE.export,
                 "ANOVA_white_logRE" = res.anova_white.LogRE.export,
                 "GLS" = res.glsLogRE.export,
                 "ART" = res.artRE.export,
                 "WRS2" = res.wrsRE.export)

# Create a workbook
wb <- createWorkbook()

# Add each data frame to the workbook
for (sheet_name in names(res.list)) {
  addWorksheet(wb, sheet_name)                          # Add a worksheet
  writeData(wb, sheet_name, res.list[[sheet_name]])  # Write data to the worksheet
}

# Save the workbook as an Excel file
saveWorkbook(wb, "combinedDHA_StatsResults.xlsx", overwrite = TRUE)

## Plot
### Import compiled p-value file
stats.RE <- read.xlsx("combinedDHA_StatsResults.xlsx")

### Change format to suit stat_pvalue_manual(): group1, group2, label, y.position
stats.RE <- stats.RE %>%
  mutate(
    ΔCq_label = Group,
    factor1_sig = ANOVA_Condition < 0.05,
    factor2_sig = ANOVA_Sample < 0.05,
    interaction_sig = `ANOVA_Condition:Sample` < 0.1  # Use 0.1 for interaction
  ) %>%
  pivot_longer(
    cols = matches("Main_|Int_"),
    names_to = "comparison", # Temporary column for names
    values_to = "p.adj" # Values go here
  ) %>%
  # Parse comparison type
  mutate(
    comparison_type = case_when(
      str_detect(comparison, "Main_Condition_") ~ "main_Condition",
      str_detect(comparison, "Main_Sample_") ~ "main_Sample",
      str_detect(comparison, "Int_Condition_") ~ "interaction_Condition",
      str_detect(comparison, "Int_Sample_") ~ "interaction_Sample"
    ),
    p.adj.signif = case_when(p.adj > 0.1 ~ NA_character_,
                             p.adj <= 0.1 & p.adj > 0.05 ~ sprintf("%.4f", p.adj),
                             p.adj <= 0.05 & p.adj > 0.01 ~ "*",
                             p.adj <= 0.01 & p.adj > 0.001 ~ "**",
                             p.adj <= 0.001 & p.adj > 0.0001 ~ "***",
                             p.adj <= 0.0001 ~ "****"),
  ) 

### For those with Main effects only
statRE.main1 <- stats.RE %>%
  dplyr::filter(interaction_sig == FALSE & factor1_sig == TRUE) %>%
  mutate(
    group1 = str_extract(comparison, "(?<=_)(D(4|10)\\.[+-]MG)(?=vs)"), # Extract group1, ?<=... look for behind this
    group2 = str_extract(comparison, "(?<=vs)(D(4|10)\\.[+-]MG)")  # Extract group2
  ) %>%
  dplyr::filter(comparison_type == "main_Condition") %>%
  dplyr::select(ΔCq_label, group1, group2, p.adj, p.adj.signif)


statRE.main2 <- stats.RE %>%
  dplyr::filter(interaction_sig == FALSE & factor2_sig == TRUE) %>%
  mutate(
    group1 = str_extract(comparison, "(?<=_)[A-Z0-9]+(?=vs)"), 
    group2 = str_extract(comparison, "(?<=vs)[A-Z0-9]+")
  ) %>%
  dplyr::filter(comparison_type == "main_Sample") %>%
  dplyr::select(ΔCq_label, group1, group2, p.adj, p.adj.signif)

### For Interaction effects
statRE.Int.C <- stats.RE %>%
  dplyr::filter(interaction_sig == TRUE & comparison_type == "interaction_Condition") %>%
  mutate(
    Condition = str_extract(comparison, "(?<=@)(D(4|10)\\.[+-]MG)"), # Extract within factor
    group1 = str_extract(comparison, "(?<=_)[A-Z0-9]+(?=vs)"), 
    group2 = str_extract(comparison, "(?<=vs)[A-Z0-9]+(?=@)")  
  ) %>%
  dplyr::select(ΔCq_label, Condition, group1, group2, p.adj, p.adj.signif)

statRE.Int.S <- stats.RE %>%
  dplyr::filter(interaction_sig == TRUE & comparison_type == "interaction_Sample") %>%
  mutate(
    Sample = str_extract(comparison, "(?<=@)[A-Z0-9]+"),
    group1 = str_extract(comparison, "(?<=_)(D(4|10)\\.[+-]MG)(?=vs)"), 
    group2 = str_extract(comparison, "(?<=vs)(D(4|10)\\.[+-]MG)(?=@)")  
  ) %>%
  dplyr::select(ΔCq_label, Sample, group1, group2, p.adj, p.adj.signif)

### Add xy position for interaction terms
#### Add X position values
x_position.C.RE <- ddCq.REoutlierRemoved %>%
  dplyr::group_by(ΔCq_label, Condition) %>%
  t_test(RE ~ Sample, p.adjust.method = "bonferroni") %>%
  add_x_position(x = "Condition") %>%
  dplyr::select(ΔCq_label, Condition, group1, group2, xmin, xmax)

x_position.S.RE <- ddCq.REoutlierRemoved %>%
  dplyr::group_by(ΔCq_label, Sample) %>%
  t_test(RE ~ Condition, p.adjust.method = "bonferroni") %>%
  add_x_position(x = "Condition", group = "Sample") %>%
  dplyr::select(ΔCq_label, Sample, group1, group2, xmin, xmax)

#### Add Y position values
y_positions.C.RE <- ddCq.REoutlierRemoved %>%
  group_by(ΔCq_label, Condition) %>% # start with the highest among this panel, not only for a Sample
  summarise(max_RE = max(RE, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(y.position = max_RE*1.01)

y_positions.S.RE <- ddCq.REoutlierRemoved %>%
  group_by(ΔCq_label) %>% # start with the highest among this panel, not only for a Sample
  summarise(max_RE = max(RE, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(y.position = max_RE*1.05)

#### Merge y.position with p value & x.position
statRE.Int.C <- statRE.Int.C %>%
  dplyr::filter(!is.na(p.adj.signif)) %>%
  left_join(y_positions.C.RE, by = c("ΔCq_label", "Condition")) %>%
  group_by(ΔCq_label, Condition) %>%
  mutate(
    y.position = ifelse(duplicated(y.position) | duplicated(y.position, fromLast = TRUE),
                        y.position + seq_along(y.position) * (0.09*y.position),  # Increment for duplicates by 10% (if seq_along(y.position) * 0.01 is add a fix 0.01 instead)
                        y.position)) %>%
  ungroup() %>%
  dplyr::select(ΔCq_label, Condition, group1, group2, p.adj, p.adj.signif, y.position) %>%
  dplyr::left_join(x_position.C.RE, by = c("ΔCq_label", "Condition", "group1", "group2")) 

statRE.Int.S <- statRE.Int.S %>%
  dplyr::filter(!is.na(p.adj.signif)) %>%
  left_join(y_positions.S.RE, by = "ΔCq_label") %>%
  group_by(ΔCq_label) %>%
  mutate(
    y.position = ifelse(duplicated(y.position) | duplicated(y.position, fromLast = TRUE),
                        y.position + seq_along(y.position) * (0.09*y.position),  # Increment for duplicates by 5% (if seq_along(y.position) * 0.01 is add a fix 0.01 instead)
                        y.position)) %>%
  ungroup() %>%
  dplyr::select(ΔCq_label, Sample, group1, group2, p.adj, p.adj.signif, y.position) %>%
  dplyr::left_join(x_position.S.RE, by = c("ΔCq_label", "Sample", "group1", "group2")) 

### Add xy position for main terms
x_position.main1.RE <- ddCq.REoutlierRemoved %>%
  dplyr::group_by(ΔCq_label) %>%
  t_test(RE ~ Condition, p.adjust.method = "bonferroni") %>%
  add_x_position(x = "Condition") %>%
  dplyr::select(ΔCq_label, group1, group2, xmin, xmax)

y_positions.main.RE <- ddCq.REoutlierRemoved %>%
  group_by(ΔCq_label) %>% # start with the highest among this panel, not only for a Sample
  summarise(max_RE = max(RE, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(y.position = max_RE*1.01) 

statRE.main1 <- statRE.main1 %>%
  dplyr::filter(!is.na(p.adj.signif)) %>%
  left_join(y_positions.main.RE, by = "ΔCq_label") %>%
  group_by(ΔCq_label) %>%
  mutate(
    y.position = ifelse(duplicated(y.position) | duplicated(y.position, fromLast = TRUE),
                        y.position + seq_along(y.position) * (0.05*y.position),  # Increment for duplicates by 5% (if seq_along(y.position) * 0.01 is add a fix 0.01 instead)
                        y.position)) %>%
  ungroup() %>%
  dplyr::select(ΔCq_label, group1, group2, p.adj, p.adj.signif, y.position) %>%
  dplyr::left_join(x_position.main1.RE, by = c("ΔCq_label", "group1", "group2")) 

x_position.main2.RE <- ddCq.REoutlierRemoved %>%
  dplyr::group_by(ΔCq_label) %>%
  t_test(RE ~ Sample, p.adjust.method = "bonferroni") %>%
  add_x_position(x = "Sample") %>%
  dplyr::select(ΔCq_label, group1, group2, xmin, xmax)

statRE.main2 <- statRE.main2 %>%
  dplyr::filter(!is.na(p.adj.signif)) %>%
  left_join(y_positions.main.RE, by = "ΔCq_label") %>%
  group_by(ΔCq_label) %>%
  mutate(
    y.position = ifelse(duplicated(y.position) | duplicated(y.position, fromLast = TRUE),
                        y.position + seq_along(y.position) * (0.05*y.position),  # Increment for duplicates by 5% (if seq_along(y.position) * 0.01 is add a fix 0.01 instead)
                        y.position)) %>%
  ungroup() %>%
  dplyr::select(ΔCq_label, group1, group2, p.adj, p.adj.signif, y.position) %>%
  dplyr::left_join(x_position.main2.RE, by = c("ΔCq_label", "group1", "group2")) 

### Add text to plot
anova_text.RE <- stats.RE %>%
  select(ΔCq_label, ANOVA_Condition, ANOVA_Sample, `ANOVA_Condition:Sample`) %>%
  dplyr::distinct(ΔCq_label, .keep_all = TRUE) %>%
  pivot_longer(
    cols = -ΔCq_label,
    names_to = "Effect",
    values_to = "p.value"
  ) %>%
  mutate(
    Effect = case_when(
      Effect == "ANOVA_Condition" ~ "Condition",
      Effect == "ANOVA_Sample" ~ "Sample",
      TRUE ~ "Condition \u00D7 Sample" #\u00D7 is the Unicode for multiply symbol
    ),
    label = paste0(Effect, ": ", ifelse(p.value < 0.001, "p < 0.001", sprintf("p = %.3f", p.value)))
  ) %>%
  # Calculate y positions (adjust based on your data range)
  left_join(y_positions.main.RE %>% select(ΔCq_label, max_RE), by = "ΔCq_label") %>%
  group_by(ΔCq_label) %>%
  arrange(
    match(Effect, c("Condition", "Sample", "Condition \u00D7 Sample")),
    .by_group = TRUE
  ) %>%
  mutate(
    # all labels sit above the max point, staggered by 10% of the panel height:
    y.position = Inf,
    vjust = 1.5 + (row_number() - 1) * 1.25, 
    # pick an x (e.g. at the left or center of the x-axis, Inf for right most):
    x.position = 0.5
  ) %>%
  ungroup()

### Plot by individual gene
#### wrapper functions
plot_gene_byf1 <- function(gene_data, y_value, stats1_list, stats2_list, statsmain_list, anova_text) {
  gene_name <- unique(gene_data$Detector)
  
  # Subset the p-value data frames for the current gene
  stats1 <- stats1_list %>% filter(str_detect(ΔCq_label, gene_name))
  stats2 <- stats2_list %>% filter(str_detect(ΔCq_label, gene_name))
  stats.main <- statsmain_list %>% filter(str_detect(ΔCq_label, gene_name))
  anno_text <- anova_text %>% filter(str_detect(ΔCq_label, gene_name))
  
  # Create the ggboxplot
  p <- ggboxplot(gene_data, x = "Condition", y = y_value, ylab = y_value,
                 color = "Sample", legend = "right",
                 add = "jitter", shape = "Sample",
                 facet.by = "ΔCq_label", ncol = 8, scales = "free_y") +
    stat_pvalue_manual(stats1, label = "p.adj.signif", label.size = 2, tip.length = 0.01)
  
  p<- p + 
    stat_pvalue_manual(stats2, label = "p.adj.signif", label.size = 2, tip.length = 0.03, step.increase = 0.02, step.group.by = "ΔCq_label") +
    stat_pvalue_manual(stats.main, label = "p.adj.signif", label.size = 2, tip.length = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  if (nrow(anno_text) > 0) {
    p <- p + geom_text(
      data = anno_text,
      aes(x = x.position, y = y.position, label = label, vjust = vjust),
      hjust = 0, size = 2.5) # geom_text with a data‑frame per facet works beautifully :contentReference[oaicite:1]{index=1}
  }
  
  return(p)
}

plot_gene_byf2 <- function(gene_data, y_value, stats_list, anova_text) {
  gene_name <- unique(gene_data$Detector)
  
  # Subset the p-value data frames for the current gene
  stats <- stats_list %>% filter(str_detect(ΔCq_label, gene_name))
  anno_text <- anova_text %>% filter(str_detect(ΔCq_label, gene_name))
  
  # Create the ggboxplot
  p <- ggboxplot(gene_data, x = "Sample", y = y_value, ylab = y_value,
                 color = "Condition", legend = "right",
                 add = "jitter", shape = "Condition", 
                 facet.by = "ΔCq_label", ncol = 8, scales = "free_y") +
    stat_pvalue_manual(stats, label = "p.adj.signif", label.size = 2, tip.length = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
  
  if (nrow(anno_text) > 0) {
    p <- p + geom_text(
      data  = anno_text,
      aes(x = x.position, y = y.position, label = label, vjust = vjust),
      hjust = 0, size = 2.5)
  }
  
  return(p)
}

save_plots <- function(plot_list, 
                       prefix = "myplot", 
                       date = Sys.Date(), 
                       path = NULL,
                       width = NA,
                       height = NA,
                       units = c("in", "cm", "mm", "px"),
                       dpi = 300,
                       limitsize = FALSE
) {
  # Get the names of the plots from the list
  plot_names <- names(plot_list)
  
  # Loop through each plot and save
  for (i in seq_along(plot_list)) {
    # Extract the name of the plot (e.g., "xxx" from plot_list$xxx)
    name <- plot_names[i]
    
    # Create the filename: "prefix_xxx_date.png"
    filename <- sprintf(
      "%s_%s_%s.png", 
      prefix, 
      name, 
      format(date, "%Y%m%d") # Date format: YYYYMMDD (adjust as needed)
    )
    
    # Save the plot using ggsave
    ggsave(
      filename = filename,
      plot = plot_list[[i]],
      path = path,
      width = width,
      height = height,
      units = units,
      dpi = dpi
    )
  }
}

#### Application:
# Split the data by Detector to plot
plotRE.gene.list <- split(ddCq.REoutlierRemoved,  ddCq.REoutlierRemoved$Detector)

# Apply the plotting function to each gene's data
plot_listRE1 <- lapply(plotRE.gene.list, plot_gene_byf1, y_value = "RE",
                      stats1_list = statRE.Int.C, stats2_list = statRE.Int.S, statsmain_list = statRE.main1, anova_text = anova_text.RE)

# for main effect by factor 2
plot_listRE2 <- lapply(plotRE.gene.list, plot_gene_byf2, y_value = "RE",
                       stats_list = statRE.main2, anova_text = anova_text.RE)
# Save
save_plots(
  plot_list = plot_listRE1,
  prefix = "RE_8refC", 
  path = file.path(DataDir, "RE_plots"),
  width = 80, 
  height = 10,
  units = "cm",
  dpi = 300
)

save_plots(
  plot_list = plot_listRE2,
  prefix = "RE_8refS", 
  path = file.path(DataDir, "RE_plots"),
  width = 80, 
  height = 10,
  units = "cm",
  dpi = 300
)