library(forestplot)

credible_snip_names <- c("credible1", "credible2", "credible3")

#mycolors <- c(male = "#619CFF", female = "#F8766D", pre = "#B79F00", post = "#00BA38") # ggplot
#mycolors <- c(male = "#0072B2", female = "#D55E00", pre = "#E69F00", post = "#009E73") # Okabe-Ito (default palette in R4)
mycolors <- c(male = "#2297E6", female = "#DF536B", pre = "#F5C710", post = "#61D04F") # R4


show_package_versions <- function(packages = NULL) {
  all <- inner_join(as_tibble(available.packages()), as_tibble(installed.packages()), 
                    by="Package", suffix=c("_available", "_installed"))
  all <- all %>% select(package=Package, installed=Version_installed, available=Version_available)
  if (is.null(packages)) {
    all
  } else {
    all %>% filter(package %in% packages)
  }
}

# Returns also a list (with class "mywarnings") of warnings 
capture_warnings <- function (expr, classes = "warning") 
{
  start <- lubridate::now()
  warnings <- list()
  result <- withCallingHandlers(expr, 
                                warning = function(w) {
                                  #cat("Got warning!\n")
                                  key <- w$message
                                  if (key %in% names(warnings)) {
                                    warnings[[key]]$count <<- warnings[[key]]$count + 1
                                    warnings[[key]]$lst   <<- c(warnings[[key]]$lst, list(w))
                                    if (inherits(w, classes)) {
                                      tryInvokeRestart("muffleWarning")
                                    }
                                    
                                  } else {
                                    warnings[[key]] <<- list(count=1, lst=list(w))   # Store the warning to the well-know global variable
                                  }
                                })
  stop <- lubridate::now()
  attr(warnings, "start") <- start
  attr(warnings, "stop")  <- stop
  attr(warnings, "class")  <- c("mywarnings", "list")
  return(list(result = result, warnings = warnings))
}

# Methods print and as_tibble of the mywarnings class:
print.mywarnings <- function(warnings) {
  my_format <- function(t, include_date = TRUE) {
    format <- if (include_date) "%Y-%m-%d %H:%M:%S" else "%H:%M:%S"
    #print(format)
    strftime(t, format = format, tz = 'GMT')
  }
  start <- attr(res$warnings, "start")
  stop <- attr(res$warnings, "stop")
  same_date <- as_date(start) == as_date(stop)
  diff <- stop-start
  # if (same_date) {
  #   start <- hms::as_hms(start)
  #   stop  <- hms::as_hms(stop)
  # }
  sprintf("Computating started at %s and stopped at %s\n", my_format(start), my_format(stop, include_date = ! same_date)) %>% cat()
  sprintf("Elapsed time was %.2f %s\n", diff, attr(diff, "units")) %>% cat()
  for (key in names(warnings)) {
    sprintf("%i copies of %s", warnings[[key]]$count, key) %>% cat()
  }
}
as_tibble.mywarnings <- function(warnings) {
  tibble(
    message = names(warnings),
    count = warnings %>% map_int("count") %>% unname()
  )
}


# Add missing combinations to the input dataframe so that the plots have each column of equal width.
# The dots give the column names whose combination you are interested in.
add_missing_combinations <- function(df, ...) {
  #variables <- rlang::list2(...)     # without this !!! in arguments won't work
  variables <- ensyms(...)          # without this symbols as arguments won't work, only strings
  tmp <- df %>% select(!!!variables) %>% as.list() %>% map(unique)
  #print(tmp)
  to_add <- expand_grid(!!!tmp) %>% anti_join(df %>% select(!!!variables))
  #print(class(to_add))
  #str(to_add)
  bind_rows(df, to_add)
}
add_missing_combinations2 <- function(df, ...) {
  variables <- rlang::list2(...)     # without this !!! in arguments won't work
  #variables <- ensyms(...)          # without this symbols as arguments won't work, only strings
  tmp <- df %>% select(!!!variables) %>% as.list() %>% map(unique)
  #print(tmp)
  to_add <- expand_grid(!!!tmp) %>% anti_join(df %>% select(!!!variables))
  #print(class(to_add))
  #str(to_add)
  bind_rows(df, to_add)
}

my_correlation <- function(df, levels=NULL) {
  tmp <- df %>% cor() %>% as_tibble(pearson, rownames="x") %>% pivot_longer(cols=-x, names_to="y")
  if (is.null(levels))
    tmp %>% mutate(across(c(x, y), function(x) factor(x)))
  else tmp %>%
    mutate(across(c(x, y), function(x) factor(x, levels=levels)))
}


#plot_correlation <- function(df, square=TRUE, size = GeomText$default_aes$size) {
# The fontsize can now be tweaked with theme(geom = element_geom(fontsize = 15))
plot_correlation <- function(df, square=TRUE) {
    df %>%
    mutate(value = if (square) value**2 else value) %>%
    mutate(y=fct_rev(y)) %>%
    ggplot(aes(x, y, fill=value)) + 
    geom_raster() +
    geom_text(aes(label=sprintf("%.2f", value)), color="white") + #, size=size) +
    scale_fill_continuous(breaks=seq(0, 1, length.out = 5), limits=c(0,1)) +
    scale_x_discrete(expand=expansion(add=0)) +
    scale_y_discrete(expand=expansion(add=0)) +
    labs(x="SNP 1", y="SNP 2", title="Linkage disequilibrium using 'cor'", fill=if (square) bquote(R^2) else bquote(R))
}

compute_cis_zi <- function(fit, exponentiate = FALSE, model="count") {
  cis <- confint(fit) %>% as_tibble(rownames = "term")
  if (model != "full")
    cis <- cis %>% filter(str_starts(term, model))
  df <- cis %>%
    mutate(estimate=coef(fit, model=model)) %>% 
    #as.data.frame() %>% 
    #rownames_to_column("term") %>%
    #as_tibble(rownames = "term") %>%
    rename(conf.low = "2.5 %", conf.high = "97.5 %") %>%
    filter(!str_detect(term, "(Intercept)")) %>%
    mutate(term = str_remove(term, sprintf("^%s_", model)))
  if (exponentiate)
    df <- df %>% mutate(across(-term, exp))
  df
}

compute_cis <- function(fit, model="count") {
  if ("zeroinfl" %in% class(fit)) {
    cis <- compute_cis_zi(fit, exponentiate = TRUE, model=model)
  } else {
    cis <- fit %>% 
      tidy(conf.int = TRUE, exponentiate = TRUE) %>%
      filter(term != "(Intercept)")
  }
  cis <- cis %>% mutate(term = str_remove_all(term, "`"))   # Remove ` characters. Interaction terms are surrounded by these since the colon character
  cis
}

# Following function computes confidence intervals for sum of predictors.
# The 'variables' column of dataframe 'df' lists all the variable names to be summed.
# The function adds columns estimate, conf.low, and conf.high to the input dataframe.
get_odds_ratio_sums <- function(fit, df, exponentiate = FALSE) {
  helper <- function(fit, variables, exponentiate) {
    n <- length(variables)
    cm <- vcov(fit)
    sub <- tryCatch(
      cm[variables, variables],
      error = function(e) { NULL }
    )
    if (is.null(sub))         # Some variables were not found in the model
      return(tibble(estimate=NA, conf.low = NA, conf.high = NA))
    # Used instructions from the following page
    # https://www.quora.com/How-would-you-calculate-a-confidence-interval-if-you-wanted-to-do-it-for-the-summation-of-the-data-set
    v <- sum(sub)   # variance of the sum of variables
    #se <- sqrt(v*n) # not completely sure about this line! Yes, this is incorrect!
    se <- sqrt(v)
    estimate  <- sum(fit$coefficients[variables])
    conf.low  <- estimate - 1.96 * se
    conf.high <- estimate + 1.96 * se
    if (exponentiate) {
      tibble(estimate = exp(estimate), conf.low = exp(conf.low), conf.high = exp(conf.high))
    } else {
      tibble(estimate = estimate, conf.low = conf.low, conf.high = conf.high)
    }
  }
  result <- df %>% #head(5) %>% 
    rowwise() %>% 
    mutate(tmp = helper(fit, variables, exponentiate = exponentiate)) %>%
    ungroup() %>% unnest(tmp)
  result
}


pretty_conv2 <- tibble::tribble(
  ~phenocode,                ~phenotype,
  "age", "Age",
  "BL_AGE", "Baseline age",
  "diagnosis_age", "Diagnosis age",
  "blood_donor", "Blood donor",
  "donations_in_five_years", "Donations in five previous years",
  "venesections", "Venesections",
  "hemochromatosis", "Hemochromatosis",
  "sex2_pre", "Premenopausal female",
  "sex2_post", "Postmenopausal female",
  "sex2pre", "Premenopausal female",
  "sex2post", "Postmenopausal female",
  "sex_female", "Female",
  "sexfemale", "Female",
  "sex", "Sex",
  "M13_ARTHROSIS",          "Arthrosis (1 year)",
  "K11_CHOLELITH",          "Cholelithiasis (5 years)",
  "I9_HYPTENSESS",          "Hypertension (1 year)",
  "E4_HYTHYNAS",            "Hypothyroidism (5 years)",
  "C3_CANCER",              "Malignant neoplasm (5 years)",
  "credible1", "Causal1 (rs181949568)",
  "credible2", "Causal2 (rs1284018747)",
  "credible3", "Causal3 (rs1300021441)",
  "rs181949568", "Causal1 (rs181949568)",
  "rs114179634", "HLA1 (rs114179634)",
  "rs9261354",   "HLA4 (rs9261354)",   # too high VIF
  "rs3094093",   "HLA2 (rs3094093)",
  "supersnip", "HLA3 (rs9272324)",
  "H1_2",      "H1-2",
  "H63D",      "H63D",
  "S65C",      "S65C",
  "C282Y",     "C282Y",
  "C282Y_X1",     "C282Y heterozygote",
  "C282Y_X2",     "C282Y homozygote",
  "C282Y1",     "C282Y heterozygote",
  "C282Y2",     "C282Y homozygote",
  "H63D_X1",      "H63D heterozygote",
  "H63D_X2",      "H63D homozygote",
  "H63D1",      "H63D heterozygote",
  "H63D2",      "H63D homozygote",
  "H1_2_X1",      "H1-2 heterozygote",
  "H1_2_X2",      "H1-2 homozygote",
) %>% deframe()
# Form pretty interaction names
create_interactions <- function(x, y, input_sep="_x_", output_sep=" x ") {
  interactions <- bind_rows(   # Make interactions in both ways
    expand_grid(a=x, b=y),
    expand_grid(a=y, b=x)
  )
  interactions %>% 
    mutate(interaction = paste(a, b, sep=input_sep), pretty=paste(pretty_conv2[a], pretty_conv2[b], sep=output_sep)) %>% 
    select(interaction, pretty) %>%
    deframe()
}
interactions1 <- create_interactions(c("C282Y_X1", "C282Y_X2"), 
                                    c("H63D_X1", "H63D_X2", "S65C", "supersnip", "donations_in_five_years", "sex2_pre", "sex2_post", "sex_female", 
                                        "rs114179634", "rs3094093", credible_snip_names))
interactions2 <- create_interactions(c("C282Y1", "C282Y2"), 
                                     c("H63D1", "H63D2", "S65C", "supersnip", "donations_in_five_years", "sex2pre", "sex2post", "sexfemale", 
                                       "rs114179634", "rs3094093", credible_snip_names),
                                     input_sep = ":")
#mutate(interaction=paste(a, b, sep="_x_"), pretty=paste(pretty_conv[a], pretty_conv[b], sep=" x ")) %>% select(interaction, pretty)
pretty_conv <- c(pretty_conv2, interactions1, interactions2)

prettify <- function(v) coalesce(pretty_conv[v], v)


# Variable order
non_genetic_variables <- c("age", "donations_in_five_years", "sex_female", "sex2pre", "sex2post", "sex2_pre", "sex2_post",   "sexfemale",
                           "M13_ARTHROSISTRUE", "K11_CHOLELITHTRUE", "I9_HYPTENSESSTRUE", "E4_HYTHYNASTRUE", "C3_CANCERTRUE"
) 
genetic_variables <- c("S65C", "rs114179634", "rs3094093", "supersnip", "credible1", "credible2", "credible3", 
                       "H63D1", "H63D2", "C282Y1", "C282Y2", 
                       "donations_in_five_years_x_C282Y_X1", 
                       "donations_in_five_years_x_C282Y_X2", 
                       "donations_in_five_years:C282Y1",
                       "donations_in_five_years:C282Y2", 
                       "sex_female_x_C282Y_X1", "sex_female_x_C282Y_X2",  "sexfemale:C282Y1", "sexfemale:C282Y2",
                       "sex2pre:C282Y1", 
                       "sex2_pre_x_C282Y_X1", 
                       "sex2pre:C282Y2", 
                       "sex2_pre_x_C282Y_X2", 
                       "sex2post:C282Y1", 
                       "sex2_post_x_C282Y_X1", 
                       "sex2post:C282Y2",
                       "sex2_post_x_C282Y_X2", 
                       "S65C:C282Y1", 
                       "supersnip_x_C282Y_X1",
                       "supersnip:C282Y1", 
                       "supersnip:C282Y2",
                       "supersnip_x_C282Y_X2",
                       "H63D1:C282Y1", 
                       "H63D_X1", "H63D_X2",  "C282Y_X1", "C282Y_X2", "H63D_X1_x_C282Y_X1", "S65C_x_C282Y_X1")

# Extract N, prevalence and R2 from a model
model_stats2 <- function(fit) {
  if (length(intersect(class(fit), c("glm", "lm"))) > 0) {   
    fit2 <- fit
    y <- fit2[["y"]]
    N <- length(y)
    cases <- sum(y==1)
    prevalence <- mean(y==1)
    R2 <- fit2 %>% fmsb::NagelkerkeR2() %>% pluck("R2")
  } else if ("coxph" %in% class(fit)) {
    m <- fit[["y"]] %>% unclass()
    N <- nrow(m)
    cases <- sum(m[, 2])
    prevalence <- mean(m[, 2])
    R2 <- NA
  } else {
    fit2 <- fit %>% extract_fit_engine()
    y <- fit2[["data"]] %>% pull(..y)
    N <- length(y)
    cases <- sum(y == "case")
    prevalence <- mean(y == "case")
    R2 <- fit2 %>% fmsb::NagelkerkeR2() %>% pluck("R2")
  }
  tibble(N=N, cases=cases, prevalence=prevalence, R2=R2)
}
model_stats <- function(fit) {
  df <- model_stats2(fit)
  sprintf("N=%i, cases=%i, prevalence=%.2f%%, R2=%.2f", df$N, df$cases, 100*df$prevalence, df$R2)
}

# prettify <- function(v, pretty_conv) coalesce(pretty_conv[v], v)  # If element does not exist in pretty_conv, keep it as it is

# Single color
make_prettify_function <- function(pretty_conv) { function(v) coalesce(pretty_conv[v], v) } # Convert label if found in pretty_conv
plot_cis <- function(cis, title, subtitle=NULL, lasso=FALSE, conversion = pretty_conv) {
  cis %>% 
    mutate(term = term %>% 
             str_remove("TRUE") %>% 
             fct_relevel(intersect(names(conversion), .)) %>%    # This could be useless now  
             fct_relabel(make_prettify_function(conversion)) %>%  # Convert to pretty labels 
             fct_rev(),                                          # Order labels from top to bottom
           hollow_color = if_else((conf.low < 1 & conf.high > 1) | (lasso & estimate == 1), NA, "black")) %>%
    ggplot(aes(x=estimate, xmin=conf.low, xmax=conf.high, y=term, fill=I(hollow_color))) + 
    geom_vline(xintercept = 1, color="gray") +
    geom_linerange() + 
    geom_point(shape=21, size=3) +
    scale_fill_discrete(na.value=NA, guide="none") +  
    scale_x_log10() +
    labs(title=title, subtitle=subtitle, x="Odds ratio", y="Variable")
}

# Multi color
plot_cis2 <- function(cis, title, subtitle=NULL, lasso = FALSE, na.rm = FALSE, conversion = pretty_conv) {
  cis %>% 
    mutate(term = term %>% 
             #str_remove("TRUE") %>% 
             fct_relabel(function(v) str_remove(v, "TRUE")) %>%  # Removing it this way, does not mess with the order of levels
             #fct_relevel(intersect(names(pretty_conv), .)) %>% 
             fct_relabel(make_prettify_function(conversion)) %>% 
             fct_rev(),
           hollow_color = if_else((conf.low < 1 & conf.high > 1) | (lasso & estimate == 1), NA, model)) %>%
    ggplot(aes(x=estimate, xmin=conf.low, xmax=conf.high, y=term, color=model, fill=hollow_color)) + 
    geom_vline(xintercept = 1, color="gray") +
    geom_linerange(position = position_dodge(width = 0.5, preserve = "total"), na.rm=na.rm) + 
    geom_point(position = position_dodge(width = 0.5, preserve = "total"), shape=21, size=3, na.rm=na.rm, show.legend = TRUE) +
    guides(color = guide_legend(override.aes = list(shape = 16) ) ) +
    scale_fill_discrete(na.value=NA, guide="none", drop = FALSE) +  
    #scale_fill_discrete(na.value=NA, drop = FALSE) +  
    scale_x_log10() +
    labs(title=title, subtitle=subtitle, x="Odds ratio", y="Variable")
}


# Helpers for named vectors.


# This sets unnamed elements to the element values and keeps existing names.
# set_names() would overwrite all names if even one name is missing
# v is a named vector
fill_names <- function(v) { 
  old_names <- names(v) %>% na_if("")
  new_names <- coalesce(old_names, unname(v))
  set_names(v, new_names)
}
# Move names to values, and values to names
# v is a named vector
swap_names_and_values <- function(v) {
  set_names(names(v), unname(v))
}
is_whole_number <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

# Show the fixed predictors (either mean (for age) or mode (for others)
# important_predictors is a (named) vector
get_fixed_predictors <- function(newdata, important_predictors) {
  helper <- function (value, name) {
    value <- unique(value)
    glue("{name}={if_else(is_whole_number(value) | is.logical(value), as.character(value), sprintf('%.1f', value))}")
  }
  newdata %>% 
    select(-all_of(important_predictors)) %>% 
    #distinct() %>% 
    as.list() %>% 
    imap(helper) %>% 
    paste(collapse=" ")
}


# Returns a tibble with confidence intervals for risks (probability of hc), covariates and train counts
compute_risk_cis <- function(lf, newdata, important_predictors, count_data = train, conf.level = 0.95) {
  variables <- important_predictors %>% unname()
  fit <- lf %>% extract_fit_engine()
  # Compute confidence intervals for predictions
  # Note: We use confidence interval instead of prediction intervals, since for logistic regression
  # prediction intervals could only be [0,0], [1,1] or [0,1].
  # Same would apply for Poisson regression as well.  
  # Confidence intervals don't have this problem as they are for the mean, not for a single observation.
  prediction_cis <- lf %>%
    extract_recipe() %>%
    bake(newdata) %>% 
    ciTools::add_ci(fit, names = c("lower", "upper"), yhatName = "mean", alpha = 1 - conf.level) %>%
    mutate(type = "parametric")
  count_variables <- setdiff(variables, "donations_in_five_years")
  counts <- count_data %>% count(across(all_of(count_variables)))
  result <- 
    bind_cols(
      newdata %>% select(all_of(variables)),
      prediction_cis %>% select(mean, lower, upper)) %>% 
    left_join(counts, by=count_variables) %>%
    replace_na(list(n=0))
  result
}

# Pivoted version of the risk table
# model is a tidymodels fitted object
# newdata is a dataframe, containing covariates we want to use for prediction
# important_predictors is a named character vector of the predictors we want to show in the risk table
# count_data is a dataframe from which the count of variable combinations is computed
compute_risk_table <- function(model, newdata, important_predictors, count_data, col="sex2", conf.level = 0.95) {
  cis <- compute_risk_cis(model, newdata, important_predictors, count_data = count_data, conf.level = conf.level)
  grouping_variables <- setdiff(unname(important_predictors), c("sex2", "sex", "donations_in_five_years"))
  tmp2 <- cis %>% 
    group_by(across(all_of(grouping_variables))) %>% 
    mutate(max_mean = max(mean), min_n = min(n), max_n = max(n)) %>%     # Compute maximums over demographic groups (sexes) and numbers of donations
    group_by(max_mean, min_n, max_n, .add=TRUE) %>% nest() %>% ungroup() %>%
    #filter(max_n >= 10) %>%
    filter(min_n >= 5) %>%
    select(-min_n) %>%
    slice_max(order_by=max_mean, n=1000) %>%   # select 32 combinations that have the highest probability (maximum over demographic groups)
    unnest(data)
  # Create string 'value' which combines the probability and count
  if (col == "sex") {
    tmp2 <- tmp2 %>%
      mutate(value = sprintf("%.3f (%i)", mean, n))
  } else {  # col == "sex2", the count could be less than 5
    tmp2 <- tmp2 %>%
      mutate(n = if_else(donations_in_five_years == 0 | TRUE, n, NA),
             value = sprintf("%.3f (%s)", mean, n))
  }
  
  labels <- tmp2 %>% 
    select(-c(max_mean, max_n, mean, lower, upper, n)) %>% 
    pivot_wider(names_from=starts_with("sex"), values_from=value)   # either sex or sex2
  
  if ("donations_in_five_years" %in% important_predictors) {
    result <- tmp2 %>% inner_join(labels, by=union(grouping_variables, "donations_in_five_years"))
  } else { # col == "sex"
    result <- tmp2 %>% inner_join(labels, by=grouping_variables)
  }
  result %>% mutate(conf.level = conf.level)
}


# This is the final version of the risk table plotting function.


# conf.level is a numeric vector of length 1 or 2 with values in range [0,1]
plot_risk_table <- function(result, important_predictors, title = NULL, col="sex2",
                            conf.level,
                            number_of_blocks = 6,
                            max_x = NULL,
                            line.margin = unit(4, 'mm')    # Space between confidence intervals
                            #conf.level = 0.95, conf.level2 = 0.80
) 
{
  # Select of top 'number_of_blocks' combinations
  if ("donations_in_five_years" %in% colnames(result)) {
    number_of_blood_donation_levels <- result$donations_in_five_years %>% n_distinct()
  } else {
    number_of_blood_donation_levels <- NA
  }
  rows_in_block = if (col == "sex2") number_of_blood_donation_levels * 3 else 1 * 2    # number of donation levels times number of sex levels
  number_of_rows <- number_of_blocks * rows_in_block
  result <- result %>% group_by(conf.level) %>% slice_head(n = number_of_rows) %>% ungroup()

  if (is.null(max_x)) {
    max_x <- result %>% pull(upper) %>% max()
    max_x <- as.integer(max_x/0.05 + 1)*0.05   # Round up to the next multiple of 0.05
  }
  # Which column to pivot by
  if (col == "sex2") {
    #important_predictors2 <- c(important_predictors, Male="male", Pre="pre", Post="post") %>% discard_at("Sex2")
    legend <- c(Male="male", Premenopausal="pre", Postmenopausal="post")
    important_predictors2 <- c(important_predictors, legend) %>% discard_at("Sex2")
    #colors <- c("blue", "darkred", "green")
    #colors <- c(mycolors["male"], mycolors["pre"], mycolors["post"])
    colors <- mycolors[legend]
    if (number_of_blood_donation_levels == 3)
      stripes <- c("#ffffff", "#EFEFEF", "#DFDFDF")
    else
      stripes <- c("#ffffff", "#DFDFDF")
  } else {  # col == "sex"
    legend <- c(Male="male", Female="female")
    important_predictors2 <- c(important_predictors, legend) %>% discard_at("Sex")
    #colors <- c("blue", "darkred")
    colors <- mycolors[legend]
    stripes <- c("#ffffff", "#DFDFDF")
  }
  variables <- important_predictors2 %>% unname()
  header <- important_predictors2 %>% swap_names_and_values() # forestplot wants these in stupid order
  
  if (col == "sex2")
    colgap <- grid::convertX(unit(2, "mm"), "npc")    # default is 6 mm, but npc units must be used
  else
    colgap <- grid::convertX(unit(6, "mm"), "npc")    # default is 6 mm, but npc units must be used
  xlab   <- "Hemochromatosis probability"
  #xticks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.30, 0.35, 0.4)  # forestplot::forestplot does not accept xlim. I use this to
  xticks <- seq(0, max_x, 0.05)  # forestplot::forestplot does not accept xlim. I use this to
  # set the bounds of the x axis
  confidence_line_color <- "black"
  txt_gp <- fpTxtGp(
    label = gpar(cex = 0.6),
    xlab = gpar(
      cex = 0.7  # double the x axis font size
      #fontsize = 7
    ),
    ticks = gpar(cex = 0.6)     # default is 0.5
  )
  txt_gp2 <- fpTxtGp(
    xlab = gpar(
      cex = 0.7,  # double the x axis font size
      alpha = 0
      #fontsize = 7
    ),
    label = gpar(alpha = 0, cex = 0.6),    # Make invisible
    legend = gpar(alpha = 0),
    #ticks = gpar(alpha = 0),
    ticks = gpar(alpha = 0, cex = 0.6)     # default is 0.5
  )
  #boxsize = 0.2
  boxsize = 0.3
  legend_names <- names(legend)
  #line.margin <- unit(4, 'mm')
  # Two different confidence levels are achieved by overlaying two plots with different confidence levels.
  # The second plot only draws the confidence lines and boxes visibly.
  rlang::inject(
    p1 <- result %>%
      filter(conf.level == .env$conf.level[[1]]) %>%
      group_by(across(all_of(col))) %>%     # This seems counter-intuitive to normal use of group_by
      forestplot(labeltext = !!variables,
                 legend = !!legend_names,
                 title = title,
                 new_page = FALSE,
                 colgap = colgap,
                 xlab = xlab,
                 #xlim = c(0, 0.4),   # This does not work
                 xticks = xticks,
                 line.margin = line.margin,
                 txt_gp = txt_gp,
                 col = fpColors(lines = colors), 
                 lwd.ci = 1.5,    # line width
                 boxsize = boxsize) %>%
      fp_set_style(box = colors %>% map(function(x) gpar(fill = x, col = x)),
                   default = gpar(lineend = "square", vertices = TRUE)) %>%
      fp_add_header(!!!header) %>% 
      fp_set_zebra_style(!!!stripes)
  )
  #par(new=TRUE)
  rlang::inject(
    p2 <- result %>%
      filter(conf.level == .env$conf.level[[2]]) %>%
      group_by(across(all_of(col))) %>%     # This seems counter-intuitive to normal use of group_by
      forestplot(labeltext = !!variables,
                 legend = !!legend_names,
                 title = "",
                 new_page = FALSE,
                 colgap = colgap,
                 zero = NA,
                 xlab = " ",
                 #xlim = c(0, 0.4),   " This does not work
                 xticks = xticks,
                 line.margin = line.margin,
                 txt_gp = txt_gp2, 
                 col = fpColors(zero = NA, lines = colors), 
                 #legend_args = fpLegend(gp = gpar(alpha = 0)),   # odd error message, cannot hide legend
#                 lwd.ci = 6,    # line width
                 lwd.ci = 4.5,    # line width
                 boxsize = boxsize) %>%
      fp_set_style(box = colors %>% map(function(x) gpar(fill = x, col = x)),
                   default = gpar(lineend = "square", vertices = TRUE)) %>%
      fp_add_header(!!!header) #%>% 
    # fp_set_zebra_style(!!!stripes)
  )
  list(p1, p2)
}


# Computes the number of exact zeros in the fitted values, and the number of times of the probability of extra zeros
# to being exactly 0 or 1.
exact_zeros_or_ones <- function(model) {
  # The fitted value is exactly 0
  n <- model %>% fitted() %>% near(0) %>% sum()
  if ("zeroinfl" %in% class(model)) {
    probabilities <- predict(model, type="zero")   # Probability of an extra zero
    zero_probability <- near(probabilities, 0) %>% sum()
    one_probability <- near(probabilities, 1) %>% sum()
  } else if ("glm" %in% class(model) && model$family$family == "binomial") {
    probabilities <- predict(model, type="response")   # Probability is exactly zero or one
    zero_probability <- near(probabilities, 0) %>% sum()
    one_probability <- near(probabilities, 1) %>% sum()
  } else {
    zero_probability <- NA
    one_probability <- NA
  }
  tibble("Zero count" = n, "Probability zero"=zero_probability, "Probability one"=one_probability)
}

constant_information_transformation <- function(fit) {
  if (! "family" %in% names(fit))
    return (function(x) x)    # Use identity function
  
  family <- fit$family$family
  # Constant-information scale transformation of the linear predictor.
  # See page 308 of Generalized linear models with examples in R.
  transformation <- if (family == "binomial") {
    function(x) asin(sqrt(x))
  } else if (family == "poisson") {
    function(x) sqrt(x)
  } else if (family == "gamma") {
    function(x) log(x)
  } else if (str_detect(family, "Negative Binomial\\(.*\\)")) {
    function(x) x
  } else {
    function(x) { stop("not implemented")}
  }
}

my_ggdiagnostics <- function(fit) {
  transformation <- constant_information_transformation(fit)
  #par(mfrow = c(2, 2))
  qr <- qresid(fit)
  #qqnorm(qr, las=1); qqline(qr)
  use_deviance_residuals <- TRUE
  if (use_deviance_residuals) {
    label <- "Standardized\ndeviance residuals"
    standardized_residuals <- rstandard(fit)
  } else {
    label <- "Standardized\nquantile residuals"
    standardized_residuals <- qresid(fit) / sqrt(1 - hatvalues(fit))
  }
  df <- tibble(qr = qr, 
               standardized_residuals = standardized_residuals,
               cit_fitted = transformation(fitted(fit)),
               linear_predictor = fit$linear.predictors,
               cook = cooks.distance(fit),
               hat  = hatvalues(fit),
               z = resid(fit, type="working") + fit$linear.predictors,   # Working responses
               index = 1:length(qr))
  g1 <- df %>% ggplot(aes(sample=qr)) + 
    geom_abline(slope=1, intercept = 0, color="red") +
    geom_qq() + labs(x="Theoretical quantiles", y="Sample quantiles")
  #plot(qr ~ asin(sqrt(fitted(fit))), las=1)   # sin^-1(sqrt(mu)) is the constant-information scale transformation for Binomial distribution    
  g2 <- df %>%
    ggplot(aes(x=cit_fitted, y=qr)) +
    geom_point() + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    labs(x="Constant information fitted", y="Quantile residual")
  g3 <- df %>% 
    ggplot(aes(x=index, y=cook, xend=index, yend=0)) +
    geom_segment() + labs(x="Index", y="Cook's distance")
  g4 <- df %>% 
    ggplot(aes(x=index, y=hat, xend=index, yend=0)) +
    geom_segment() + labs(x="Index", y="Hatvalue")
  
  # Why does the fitted curve look different from the one created by base graphics?
  g5 <- df %>%
    ggplot(aes(x=cit_fitted, y=standardized_residuals)) +
    geom_point() + 
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    #geom_smooth(method = "loess", span=2/3, method.args=list(span = 2/3, degree = 1, family="gaussian"), se=FALSE) +
    labs(y=label, x = "Constant information fitted")
  # scatter.smooth(standardized_residuals ~ transformation(fitted(fit)), las=1,
  #                ylab=label,
  #                xlab="Constant information fitted")
  
  g6 <- df %>%
    ggplot(aes(x=linear_predictor, y=z)) +
    geom_point() + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
    labs(x="Working responses", y="Linear predictor")
  # plot(z ~ fit$linear.predictors, las=1, xlab="Working responses",
  #      ylab="Linear predictor")
  # abline(0, 1)
  #plot(cooks.distance(fit), type="h", las=1)
  #plot(hatvalues(fit), type="h", las=1)
  (g1 + g2) / (g3 + g4) / (g5 + g6)
}

variable_plots <- function(fit, variables="age") {
  df1 <- fit$data %>% select(all_of(variables)) %>%
    mutate(y = rstandard(fit)) %>% 
    pivot_longer(cols = all_of(variables), values_to="x", names_to = "variable")
  # df1 <- tibble(standardized_residuals = rstandard(fit),
  #               value = fit$data$age)
  # g1 <- df1 %>% ggplot(aes(x=value, y=standardized_residuals)) + geom_point() +
  #   labs(x="age", y="Standardized deviance residuals")
  # Partial residual plot
  # This probably does not make sense when model contains interactions
  tmp <- termplot(fit, data=fit$data, 
                  partial.resid = TRUE, terms = variables, las=1, plot = FALSE)
  df2 <- bind_rows(tmp, .id="variable") %>% as_tibble()
  df <- bind_rows(full=df1, partial=df2, .id="type")
  #g2 <- df2 %>% ggplot(aes(x, y)) + geom_point() + labs(x="age")
  g <- df %>% ggplot(aes(x, y)) + geom_point() +
    labs(x="value", y="residual") +
    ggh4x::facet_grid2(rows=vars(variable), cols=vars(type),
                       #variable ~ type, 
                       scales="free", independent=TRUE) 
  #g1 + g2
  g
}


# A copy of the code that cut function uses
breaks_default <- function(x, number_of_breaks) 
{
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  if (length(number_of_breaks) != 1L)
    stop("Length of breaks must be one")
  if (is.na(number_of_breaks) || number_of_breaks < 2L) 
    stop("invalid number of intervals")
  nb <- as.integer(number_of_breaks + 1)
  dx <- diff(rx <- range(x, na.rm = TRUE))
  if (dx == 0) {
    dx <- if (rx[1L] != 0) 
      abs(rx[1L])
    else 1
    breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000, 
                      length.out = nb)
  }
  else {
    breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
    breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] + 
                             dx/1000)
  }
  breaks
}
# This solution unites neighbouring bins until they are large enough.
# Note, that division into bins is slightly different, than what geom_histogram does by default.
# If df is grouped, then the breaks are common in all groups.
cut_helper_variable_width <- function(df, column="a", number_of_bins=30, limit=5, debug=FALSE) {
  if (is.vector(df) || is.Date(df)) {    # Convert vector to a one-column dataframe
    df <- tibble(a = df)
  }
  if (is.numeric(df[[column]])) {
    breaks <- breaks_default(df[[column]], number_of_bins)
  } else {  # Assume the object type is date
    breaks <- breaks_default(unclass(df[[column]]), number_of_bins)
    breaks <- as_date(breaks)
  }
  #x <- cut(x, breaks=breaks)
  helper <- function(x, breaks) {
    number_of_bins <- length(breaks) - 1
    cut(x, breaks = breaks) %>% 
      table() %>% as_tibble() %>% rename(bin=".") %>% 
      mutate(low = breaks[1:number_of_bins], high = breaks[2:(number_of_bins+1)])
  }
  combine_bins <- function(df, first, second) {
    df[[first, "n"]] <- df[[first, "n"]] + df[[second, "n"]]   # Add the counts in the bins
    df[[first, "high"]] <- df[[second, "high"]]                # Extend the first bin
    df <- df[-second,]                                         # Delete the second bin
    df
  }
  my_min <- function(v) {
    v <- as.double(v)
    v[v==0] <- Inf
    min(v)
  }
  #df <- helper(x, breaks)
  
  grouping_variables <- group_vars(df)
  binned <- df %>%
    #group_by(operation) %>%
    group_map(function(df, key)  helper(df[[column]], breaks) %>% bind_cols(key)) %>%
    bind_rows()
  min_count <- binned %>% group_by(low, high) %>% summarise(n = my_min(n), .groups = "drop")
  #min_count[min_count==0] <- Inf
  if (length(grouping_variables) > 0)
    binned <- binned %>% group_by(across(all_of(grouping_variables)))
  # Start uniting the smallest bin with its neighbour
  i <- which.min(min_count$n)  # Unite the smallest bin
  #print(i)
  while (min_count[[i, "n"]] < limit) {
    #cat(sprintf("i=%i, class(i)=%s\n", i, class(i)))
    n <- nrow(min_count)
    i2 <-               # with this neighbouring bin
      if (i==1) 2 else {
        if (i==n) n-1 else {
          if (min_count[[i-1, "n"]] < min_count[[i+1, "n"]]) { 
            i-1         # Unite with left bin
          } else {
            i+1         # Unite with right bin
          }
        }
      }
    first <- min(i, i2)
    second <- max(i, i2)
    if (debug) cat(sprintf("Combining rows %i and %i from a tibble with %i rows\n", first, second, n))
    min_count[min_count==Inf] <- 0
    binned <- binned %>% group_map(function(df, key) combine_bins(df, first, second) %>% bind_cols(key)) %>% bind_rows()
    min_count <- binned %>% group_by(low, high) %>% summarise(n = my_min(n), .groups = "drop")
    #min_count[min_count==0] <- Inf
    if (length(grouping_variables) > 0)
      binned <- binned %>% group_by(across(all_of(grouping_variables)))
    i <- which.min(min_count$n)  # Unite the smallest bin
    #print(i)
    #df <- unite_small_bins(df, limit)
  }
  breaks <- c(min_count$low, min_count$high[nrow(min_count)])
  return(list(breaks=breaks, df=binned))
}

# Gets the size of the smallest bin. Returns a tibble, possibly with many rows.
min_bin_size <- function(g, limit=5) {
  df <- layer_data(g)
  m <- min(setdiff(df$count, 0))
  number_of_bins <- df %>% select(xmin, xmax) %>% distinct() %>% nrow()
  msg <- sprintf("Size of the smallest bin is %i, number of bins is %i", m, number_of_bins)
  if (m < limit) {
    warning(msg)
    df %>% filter(count == m)
  } else {
    message(msg)
    invisible(NULL)
  }
  
}
small_bins <- function(g, limit=5) {
  df <- layer_data(g)
  df %>% filter(count < limit)
}




