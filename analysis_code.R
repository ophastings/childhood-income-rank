# The file contains all of the analysis in the paper:
# "Growing up Different(ly than Last Time We Asked): Social Status and Changing Reports of Childhood Income Rank"
#
# Amber Obermaier and Orestes P. Hastings
#
# Last updated: Aug 11, 2025
#
# Uses the version of the gssr package that was updated on Nov 11, 2024 
# (data types of some variables differs from previous versions, so extra data wrangling could be required)

library(tidyverse)
library(huxtable)
library(haven)
library(janitor)

library(sandwich) 
library(lmtest)

library(fixest) 
library(lme4) 
library(lmerTest) 

library(skimr)

# install.packages("devtools") 
# devtools::install_github("davidsjoberg/ggsankey", force = T)
library(ggsankey)

# Install 'gssr' from 'ropensci' universe
#install.packages('gssr', repos =
#                   c('https://kjhealy.r-universe.dev', 'https://cloud.r-project.org'))

# Also recommended: install 'gssrdoc' as well
#install.packages('gssrdoc', repos =
#                 c('https://kjhealy.r-universe.dev', 'https://cloud.r-project.org'))

library(gssr, gssrdoc)

#######################
# Change this to the directory where you want your result to go
setwd("results")
#######################

#### Load the data ##### 

data(gss_panel10_long)
data(gss_panel08_long)
data(gss_panel06_long)


# a uniqe id is created by adding a 6 digit number to each id
gss_panel06_long <- gss_panel06_long %>% 
  mutate(id = as.numeric(paste0(firstid)) + 100000) %>%
  mutate(firstyear = 2006 ) %>%
  relocate(firstyear, id)

gss_panel08_long <- gss_panel08_long %>% 
  mutate(id = as.numeric(paste0(firstid)) + 200000) %>%
  mutate(firstyear = 2008) %>%
  relocate(firstyear, id)

gss_panel10_long <- gss_panel10_long %>% 
  mutate(id = as.numeric(paste0(firstid)) + 300000) %>%
  mutate(firstyear = 2010)

# these two variables (cshutyp06 and version) were giving some trouble before combing rows, and are not needed...
gss_panel06_long <- gss_panel06_long %>%
  select(-version)

gss_panel08_long <- gss_panel08_long %>%
  select(-cshutyp06) %>%
  select(-version)
  
gss_panel10_long <- gss_panel10_long %>%
  select(-cshutyp06) %>%
  select(-version)


#### Merge the data and clean ####

gss_panel_long <- bind_rows(gss_panel06_long, gss_panel08_long, gss_panel10_long) |>
  relocate(firstyear, id)


# Drop if no 2nd wave weight, which means the respondent did not appear at least twice
gss_panel_long <- gss_panel_long %>%
  filter(!is.na(wtpannr12))


# standardize incom16
gss_panel_long <- gss_panel_long %>%
  mutate(zincom16 =as.numeric(scale(as.numeric(incom16))))

# create income measures
gss_panel_long <- gss_panel_long %>%
  mutate(realinc = as.numeric(as.character(realinc)),
         log_income = log(realinc*2.78),
         zinc = as.numeric(scale(log_income)),
         inc1000 = realinc*2.78/1000,
         incequiv = realinc/sqrt(as.numeric(hompop)))

# reverse satin so higher is better and standardize
gss_panel_long <- gss_panel_long %>%
  mutate(reversesatfin = 3 - as.numeric(satfin)) %>%
  mutate(zsatfin =as.numeric(scale(reversesatfin)))

# standardize finrela
gss_panel_long <- gss_panel_long %>%
  mutate(zfinrela =as.numeric(scale(as.numeric(finrela))))

gss_panel_long <- gss_panel_long %>%
  mutate(year_f = factor(year))

gss_panel_long <- gss_panel_long %>%
  mutate(marital_f = as_factor(marital),
         marital_f = factor(tolower(marital_f)),
         marital_f = relevel(marital_f, ref = "married"))
  
gss_panel_long <- gss_panel_long %>%
  mutate(sex = as_factor(sex),
         sex = factor(tolower(sex)),
         sex = relevel(sex, ref = "male"))
  
gss_panel_long <- gss_panel_long %>%
  mutate(age = as.numeric(ifelse(age == "89 or Older", "89", as.character(age))),
         educ = as.numeric(as.character(educ)),
         childs = as.numeric(as.character(childs)),
         attend = as.numeric(attend) - 1) # because it was coded from 0 (but doesn't actually matter)

# create a race/ethnicity variable
gss_panel_long <- gss_panel_long %>%
  mutate(raceeth = case_when(
    as.numeric(hispanic) == 1 & as.numeric(race) == 1 ~ "nh white",
    as.numeric(hispanic) == 1 & as.numeric(race) == 2 ~ "nh black",
    as.numeric(hispanic) > 1 & as.numeric(hispanic) <= 50 ~ "hispanic",
    as.numeric(hispanic)==1 & as.numeric(race)==3~"other", TRUE~NA)) %>%
  mutate(raceeth = factor(raceeth, levels = c("nh white", "nh black", "hispanic", "other"))) 

# keep variables still to be used in the analysis
gss_panel_sum <- gss_panel_long %>%
  dplyr::select(id, wave, zincom16, incom16, log_income, inc1000, incequiv, realrinc, realinc, sex, raceeth, finrela, satfin, finalter, age, born, wtpannr123, wtpannr12, log_income, zsatfin, zfinrela, zinc, marital_f,attend,childs, educ,income, realrinc, hompop, year, year_f, firstyear)

# Create wave 1 time-invariant measures
gss_panel_wide <- gss_panel_sum %>%
  pivot_wider(names_from = wave, values_from = -c(wave, id, wtpannr12, wtpannr123, firstyear))

# drop any respondents not born in the U.S.
gss_panel_wide <- gss_panel_wide |>
  filter(born_1 != "No" & born_2 != "No" & born_3 != "No")

# check dropped in not born in US
gss_panel_wide |>
  tabyl(born_1)

# flag for an older than 70 at any point across waves
gss_panel_wide <- gss_panel_wide %>%
  mutate(oldie = case_when(
    age_1 > 70 | age_2 > 70 | age_3 > 70 ~ TRUE,
    TRUE ~ FALSE))

# flag if sex changes
gss_panel_wide <- gss_panel_wide %>%
  mutate(sex_change = case_when(
    # Check for changes between sex_1 and sex_2, sex_2 and sex_3, or sex_1 and sex_3, ignoring NA
    sex_1 != sex_2 & !is.na(sex_1) & !is.na(sex_2) ~ TRUE,
    sex_2 != sex_3 & !is.na(sex_2) & !is.na(sex_3) ~ TRUE,
    sex_1 != sex_3 & !is.na(sex_1) & !is.na(sex_3) ~ TRUE,
    TRUE ~ FALSE ))

# flag if raceeth changes
gss_panel_wide <- gss_panel_wide %>%
  mutate(raceeth_change = case_when(
    # Check for changes between raceeth_1 and raceeth_2, raceeth_2 and raceeth_3, or raceeth_1 and raceeth_3, ignoring NA
    raceeth_1 != raceeth_2 & !is.na(raceeth_1) & !is.na(raceeth_2) ~ TRUE,
    raceeth_2 != raceeth_3 & !is.na(raceeth_2) & !is.na(raceeth_3) ~ TRUE,
    raceeth_1 != raceeth_3 & !is.na(raceeth_1) & !is.na(raceeth_3) ~ TRUE,
    TRUE ~ FALSE ))

# flag if incom16 changes
gss_panel_wide <- gss_panel_wide %>%
  mutate(incom16_change = case_when(
    incom16_1 != incom16_2 & !is.na(incom16_1) & !is.na(incom16_2) ~ TRUE,
    incom16_2 != incom16_3 & !is.na(incom16_2) & !is.na(incom16_3) ~ TRUE,
    incom16_1 != incom16_3 & !is.na(incom16_1) & !is.na(incom16_3) ~ TRUE,
    TRUE ~ FALSE ))

# How many people change?
gss_panel_wide %>%
  tabyl(incom16_change)

# incom16_change  n  percent
# FALSE         1513 0.435521
# TRUE          1961 0.564479

gss_panel_wide <- gss_panel_wide |>
  relocate(incom16_1, incom16_2, incom16_3, incom16_change)

gss_panel_wide %>%
  tabyl(sex_change)

gss_panel_wide %>%
  tabyl(raceeth_change)

# create for variable for race and sex at wave 1
gss_panel_wide <- gss_panel_wide %>%
  mutate(Isex = sex_1)

gss_panel_wide <- gss_panel_wide %>%
  mutate(Irace = raceeth_1)   


# reshape to long form
gss_panel_long2 <- gss_panel_wide %>%
  pivot_longer(
    cols = -c(id, wtpannr12, wtpannr123, Isex, Irace, sex_change, raceeth_change, incom16_change, firstyear, oldie),
    names_to = c(".value", "wave"),  
    names_pattern = "(.*)_(\\d+)")

# For now I'm dropping the oldies
gss_panel_wide <- gss_panel_wide |>
  filter(oldie == FALSE)

gss_panel_long2 <- gss_panel_long2 |>
  filter(oldie == FALSE)



#### Create the Sankey Diagram #####

df2 <- gss_panel_wide |>
  dplyr::select(incom16_1, incom16_2, incom16_3) 

df2 <- df2 |> mutate("wave 1"= factor(tolower(as_factor(incom16_1))),
                     "wave 2" = factor(tolower(as_factor(incom16_2))),
                     "wave 3" = factor(tolower(as_factor(incom16_3)))
) |>
  drop_na("wave 1", "wave 2", "wave 3")


dflong <- df2 |> 
  make_long("wave 1", "wave 2", "wave 3")


dflong$node <- factor(dflong$node,levels = c("far below average", "below average", "average", "above average", "far above average"))
dflong$next_node <- factor(dflong$next_node,levels = c("far below average", "below average", "average", "above average", "far above average"))


pl <- ggplot(dflong, aes(x = x,                        
                         next_x = next_x,                                     
                         node = node,
                         next_node = next_node,        
                         fill = node))

pl <- pl +geom_sankey(flow.alpha = 0.5,          #This Creates the transparency of your node 
                      node.color = "black",     # This is your node color        
                      show.legend = TRUE,  # This determines if you want your legend to show
                      width = .2)        #

# pl + labs(x = "waves are two years apart",
#           fill = "childhood income\nrank at age 16") +
#   theme_sankey(base_size = 12) +
#   theme(legend.position = "bottom") +   # Move legend to bottom
#   scale_fill_manual(breaks = c("far above average", "above average", "average", "below average", "far below average"), 
#                     values = c("#66C2A5" ,"#FC8D62" ,"#8DA0CB", "#E78AC3", "#A6D854"),
#                     guide = guide_legend(nrow = 3))

pl + labs(x = "waves are two years apart",
          fill = "childhood income\nrank at age 16") +
  theme_sankey(base_size = 12) +
  theme(legend.position = "right") +   # Move legend to bottom
  scale_fill_manual(breaks = c("far above average", "above average", "average", "below average", "far below average"), 
                    values = c("#66C2A5" ,"#FC8D62" ,"#8DA0CB", "#E78AC3", "#A6D854"),
                    )

ggsave("sankey.pdf", width = 8, height = 4, units = "in")

# ggsave("sankey2.pdf", width = 6, height = 5.5, units = "in")


#### Descriptives #### 

descriptives <- na.omit(gss_panel_long2[, c("incom16", "finrela", "satfin", "log_income", 
                                            "zincom16", "zinc", "zfinrela", "zsatfin",
                                            "age", "educ", "marital_f", "childs", "attend", "year_f",
                                            "raceeth", "sex"
                                            )]) %>%
  mutate(incom16 = as.numeric(incom16),
         finrela = as.numeric(finrela),
         satfin = 3 - as.numeric(satfin))

d <- skim(descriptives)

# Extract mean and sd for numeric variables from 'd'
numeric_summary <- d %>%
  filter(skim_type == "numeric") %>%
  dplyr::select(skim_variable, numeric.mean, numeric.sd) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) # round

# Process factor variables to get counts and proportions
factor_vars <- names(descriptives)[sapply(descriptives, is.factor)]

# For each factor variable, calculate counts and proportions
factor_summary <- descriptives %>%
  select(all_of(factor_vars)) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "level") %>%
  group_by(variable, level) %>%
  summarize(n = n(), .groups = 'drop') %>%
  group_by(variable) %>%
  mutate(proportion = n / sum(n)) |>
  mutate(across(where(is.numeric), ~ round(.x, 2))) # round

# Display the summaries
print("Numeric Variables Summary:")
print(numeric_summary)

print("Factor Variables Summary:")
print(factor_summary)

write_csv(numeric_summary, "descriptives_numeric_new.csv")
write_csv(factor_summary, "descriptives_factor_new.csv")

# By hand combine these into one nice table in word



#### Data for regression models #### 

# Subset the dataframe to include only the variables used in the regression model
gss_panel_long2_complete <- na.omit(gss_panel_long2[, c("zincom16", "zinc", "zfinrela", "zsatfin", 
                                                        "age", "educ", "marital_f", "childs", "attend", "year_f", "wtpannr12", "id", 
                                                        "raceeth", "sex", "Irace", "Isex", "raceeth_change", "sex_change"
                                                        )])

## , "inc1000", "incequiv" (Add these two back into the list above is checking for robustness)


#### Pooled Models #### 

# Define the formulas for the models
formulas <- list(
  zincom16 ~ zinc + age + educ + raceeth + sex + marital_f + childs + attend + year_f ,
  zincom16 ~ zsatfin + age + educ + raceeth + sex + marital_f + childs + attend + year_f ,
  zincom16 ~ zfinrela + age + educ + raceeth + sex + marital_f + childs + attend + year_f ,
  zincom16 ~ zinc + zfinrela + zsatfin + age + educ + raceeth + sex + marital_f + childs + attend + year_f 
)

# huxreg(models[[1]], models[[2]], models[[3]], models[[4]], coef = c("log_income" = "income (logged)"))

# Initialize an empty list to store results
results <- list()
models <- list()

# Fit models and store results
for (i in 1:length(formulas)) {
  model_name <- paste("Model", i ,sep = "")
  temp_model <- lm(formulas[[i]], data = gss_panel_long2_complete, weights = wtpannr12)
  cluster <- coeftest(temp_model, vcov = vcovCL, cluster = ~ id, save = T) # this gives the clustered standard errors
  models[[i]] <- cluster 
  coefficients <- coef(cluster)
  conf_intervals <- confint(cluster)
  r <- data.frame(coef = coefficients, lower = conf_intervals[, 1], upper = conf_intervals[, 2])
  results[[i]] <- r |> mutate(model = model_name,
                              variable = rownames(r))
}

huxreg(models[[1]])

result_pooled <- huxreg(models[[1]], models[[2]], models[[3]], models[[4]], 
       coefs = c("income (logged, standardized)" = "zinc",
                 "financial satisfaction (standardized)" = "zsatfin",
               "perceived relative income (standardized)" = "zfinrela" , 
               "female" = "sexfemale",
               "non-Hispanic black" = "raceethnh black",
               "Hispanic" = "raceethhispanic",
               "other race/ethnicity" = "raceethother",  
               "age" = "age",
               "years of education" = "educ",
               "number of children" = "childs",
               "religious service attendance" = "attend",
               "marital status: widowed" = "marital_fwidowed",
               "marital status: divorced" = "marital_fdivorced",
               "marital status: separated" = "marital_fseparated",
               "marital status: never married" = "marital_fnever married"
               ),
       statistics = c(
         "N (Observations)" = "nobs",
         "R-squared" = "r.squared"
       )
)

quick_docx(result_pooled, file = "pooled.docx")

combined_df <- NULL
for (i in 1:length(formulas)) {
  combined_df <- rbind(combined_df, results[[i]])
}


combined_df_keep <- combined_df |>
  filter(variable == "zinc" | variable == "zsatfin" | variable == "zfinrela" ) |>
  mutate(model = str_replace_all(model, "Model1", "Separate"),
         model = str_replace_all(model, "Model2", "Separate"),
         model = str_replace_all(model, "Model3", "Separate"),
         model = str_replace_all(model, "Model4", "Together"),
         variable = str_replace_all(variable, "zinc", "income (logged)"),
         variable = str_replace_all(variable, "zfinrela", "perceived relative income"),
         variable = str_replace_all(variable, "zsatfin", "financial satisfaction"),
         model = factor(model, levels = c("Together", "Separate")), # Reorder model to display Separate first
  )


ggplot(combined_df_keep, aes(x = fct_rev(fct_relevel(variable, "income (logged)")), y = coef, color = model, shape = model)) +
  geom_hline(yintercept=0, lwd=1, colour="grey50") +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), width = 0.2, linewidth = .8) +
  theme_minimal() +
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Variable", y = "Coefficient (Standardized)") +    # title = "Coefficient Plots from Pooled Models (no Fixed Effects)"
  coord_flip() +
  scale_color_manual(values = c("Separate" = "#1F77B4", "Together" = "#FF7F0E"), breaks = c("Separate", "Together")) + 
  scale_shape_manual(values = c("Separate" = 16, "Together" = 17), breaks = c("Separate", "Together")) 


ggsave("pooled.pdf", width = 8, height = 4, units = "in")


#### Functions for Fixed Effects Regressions ####

# This function is for adding the number of respondents (level 2) to the model output
add_ngrps_to_model <- function(model, fe_variable = "id") {
  # Extract the number of groups (respondents) from the model's fixed effects
  fixef_sizes <- summary(model)$fixef_sizes
  
  # Check if the fixed effects variable is present
  if (fe_variable %in% names(fixef_sizes)) {
    ngrps_value <- fixef_sizes[fe_variable]
  } else {
    warning(paste("Fixed effect variable", fe_variable, "not found in the model's fixed effects."))
    ngrps_value <- NA
  }
  
  # Use tidy_override to add 'ngrps' to the glance statistics
  model2 <- tidy_override(model, 
                          glance = list(
                            ngrps = as.integer(ngrps_value)
                          ), 
                          extend = TRUE
  )
  return(model2)
}

# Function to extract coefficients and standard errors for the figures
extract_coefficients <- function(model) {
  coef <- coef(model)
  se <- se(model)
  coef_se_df <- as.data.frame(cbind(coef, se))
  
  coef_se_df <- coef_se_df |>
    mutate(variable = rownames(coef_se_df)) |>
    mutate(lower = coef - 1.96*se) |>
    mutate(upper = coef + 1.96*se)
  
  return(coef_se_df)
}

#### Check on the Variance

vars <- c("zinc", "zfinrela", "zsatfin")

decomp_tbl <- gss_panel_long2_complete %>%
  select(id, all_of(vars)) %>%    # keep id and vars (weights not used here)
  pivot_longer(-id, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    var_total   = var(value, na.rm = TRUE),
    # between = variance of the id-means
    var_between = var( tapply(value, id, mean, na.rm = TRUE), na.rm = TRUE ),
    # within = variance of (value − person_mean)
    var_within  = var( value - ave(value, id, FUN = function(z) mean(z, na.rm = TRUE)),
                       na.rm = TRUE ),
    pct_between = var_between / var_total,
    pct_within  = var_within  / var_total
  ) %>%
  arrange(desc(pct_within))

decomp_tbl |>
  mutate(across(where(is.numeric), ~ round(., 2))) |>
  write_csv("variance_decomposition_rounded.csv")

write_csv(decomp_tbl, "variance_decomposition.csv")

#### Main Fixed Effects Models ####

m1afe<-feols(zincom16 ~ zinc + marital_f+childs+attend + year_f|id, data = gss_panel_long2_complete , weights=~wtpannr12)
m1bfe<-feols(zincom16 ~ zfinrela + marital_f+childs+attend + year_f|id, data = gss_panel_long2_complete, weights=~wtpannr12)
m1cfe<-feols(zincom16 ~  zsatfin + marital_f+childs+attend + year_f|id, data = gss_panel_long2_complete, weights=~wtpannr12)
m1allfe<-feols(zincom16 ~ zinc + zfinrela + zsatfin + marital_f+childs+attend + year_f|id, data = gss_panel_long2_complete, weights=~wtpannr12)

huxreg(m1afe, m1bfe, m1cfe, m1allfe) 

summary(m1allfe)$r2["within"]

# check that equivaled and logged income produce similar null results (they do)
#m1allfe_alt1<-feols(zincom16 ~ incequiv + zfinrela + zsatfin + marital_f+childs+attend + year_f|id, data = gss_panel_long2_complete, weights=~wtpannr12)
#m1allfe_alt2<-feols(zincom16 ~ inc1000 + zfinrela + zsatfin + marital_f+childs+attend + year_f|id, data = gss_panel_long2_complete, weights=~wtpannr12)
#huxreg(m1allfe, m1allfe_alt1, m1allfe_alt2, statistics = c("AIC", "BIC"))


m1afe2 <- add_ngrps_to_model(m1afe)
m1cfe2 <- add_ngrps_to_model(m1cfe)
m1bfe2 <- add_ngrps_to_model(m1bfe)
m1allfe2 <- add_ngrps_to_model(m1allfe)


resultfe <- huxreg(m1afe2, m1cfe2, m1bfe2, m1allfe2, 
                 coefs = c("income (logged, standardized)" = "zinc",
                           "financial satisfaction (standardized)" = "zsatfin",
                           "perceived relative income (standardized)" = "zfinrela",
                           "number of children" = "childs",
                           "religious service attendance" = "attend",
                           "marital status: widowed" = "marital_fwidowed",
                           "marital status: divorced" = "marital_fdivorced",
                           "marital status: separated" = "marital_fseparated",
                           "marital status: never married" = "marital_fnever married"
                 ),
                 statistics = c(
                   "N (Observations)" = "nobs",
                   "N (Respondents)" = "ngrps"
                 )
)

quick_docx(resultfe, file = "fe_main_new.docx")



# Extract coefficients and standard errors for each model separately
coef_m1afe <- extract_coefficients(m1afe)
coef_m1bfe <- extract_coefficients(m1bfe)
coef_m1cfe <- extract_coefficients(m1cfe)
coef_m1allfe <- extract_coefficients(m1allfe)


# Add a column to indicate the model
coef_m1afe$model <- "Separate Models"
coef_m1bfe$model <- "Separate Models"
coef_m1cfe$model <- "Separate Models"
coef_m1allfe$model <- "Combined Model"


# Combine the results
combined_coef_se <- rbind(coef_m1afe, coef_m1bfe, coef_m1cfe, coef_m1allfe) |>
  filter(variable == "zinc" | variable == "zsatfin" | variable == "zfinrela" ) |>
  mutate(
    variable = str_replace_all(variable, "zinc", "income (logged)"),
    variable = str_replace_all(variable, "zfinrela", "perceived relative income"),
    variable = str_replace_all(variable, "zsatfin", "financial satisfaction"),
    model = factor(model, levels = c("Combined Model", "Separate Models")), # Reorder model to display Separate first
  )

ggplot(combined_coef_se, aes(x = fct_rev(fct_relevel(variable, "income (logged)")), y = coef, color = model, shape = model)) +
  geom_hline(yintercept=0, lwd=1, colour="grey50") +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), width = 0.2, linewidth = .8) +
  theme_minimal() +
  labs(x = "Variable", y = "Coefficient (Standardized)") +   #, title = "Fixed Effects Coefficients") +
  coord_flip() +
  scale_color_manual(values = c("Separate Models" = "#1F77B4", "Combined Model" = "#FF7F0E"), breaks = c("Separate Models", "Combined Model")) + 
  scale_shape_manual(values = c("Separate Models" = 16, "Combined Model" = 17), breaks = c("Separate Models", "Combined Model")) 


ggsave("fe_main.pdf", width = 8, height = 4, units = "in")


#### Models by Race ####

gss_panel_long2_raceeth <- gss_panel_long2_complete |>
  filter(raceeth_change == FALSE)

gss_nhwhite <- gss_panel_long2_raceeth |>
  filter(Irace == "nh white")

mfe_white <- feols(zincom16 ~ zinc + zfinrela + zsatfin + marital_f+childs+attend + year_f|id, data = gss_nhwhite, weights=~wtpannr12)


gss_nhblack <- gss_panel_long2_raceeth |>
  filter(Irace == "nh black")

mfe_black <- feols(zincom16 ~ zinc + zfinrela + zsatfin + marital_f+childs+attend + year_f|id, data = gss_nhblack, weights=~wtpannr12)


gss_hispanic <- gss_panel_long2_raceeth |>
  filter(Irace == "hispanic")

mfe_hispanic <- feols(zincom16 ~ zinc + zfinrela + zsatfin + marital_f+childs+attend + year_f|id, data = gss_hispanic, weights=~wtpannr12)

# interaction term model
mfe_raceall <- feols(zincom16 ~ Irace*(zinc + zfinrela + zsatfin + marital_f+childs+attend + year_f) |id, data = gss_panel_long2_raceeth, weights=~wtpannr12)

huxreg(mfe_raceall)

# add the number of respondents (level 2) to each model
mfe_white2 <- add_ngrps_to_model(mfe_white)
mfe_black2 <- add_ngrps_to_model(mfe_black)
mfe_hispanic2 <- add_ngrps_to_model(mfe_hispanic)


resultfe_race <- huxreg("NH White" = mfe_white2, "NH Black" = mfe_black2, "Hispanic" = mfe_hispanic2, 
                   coefs = c("income (logged, standardized)" = "zinc",
                             "financial satisfaction (standardized)" = "zsatfin",
                             "perceived relative income (standardized)" = "zfinrela",
                             "number of children" = "childs",
                             "religious service attendance" = "attend",
                             "marital status: widowed" = "marital_fwidowed",
                             "marital status: divorced" = "marital_fdivorced",
                             "marital status: separated" = "marital_fseparated",
                             "marital status: never married" = "marital_fnever married"
                   ),
                   statistics = c(
                     "N (Observations)" = "nobs",
                     "N (Respondents)" = "ngrps"
                   )
)

quick_docx(resultfe_race, file = "fe_race.docx")

## interaction model
mfe_raceall2 <- add_ngrps_to_model(mfe_raceall)

# Get the exact names of the coefficients in the model
#coef_names <- names(coef(mfe_raceall))
#print(coef_names)

# Create the huxreg table for mfe_raceall
resultfe_raceall <- huxreg(
  mfe_raceall2,
  coefs = c(
    # Main effects
    "Income (standardized)" = "zinc",
    "Perceived Relative Income (standardized)" = "zfinrela",
    "Financial Satisfaction (standardized)" = "zsatfin",
    
    # Interactions with NH Black
    "Income × NH Black" = "Iracenh black:zinc",
    "Perceived Relative Income × NH Black" = "Iracenh black:zfinrela",
    "Financial Satisfaction × NH Black" = "Iracenh black:zsatfin",
    
    # Interactions with Hispanic
    "Income × Hispanic" = "Iracehispanic:zinc",
    "Perceived Relative Income × Hispanic" = "Iracehispanic:zfinrela",
    "Financial Satisfaction × Hispanic" = "Iracehispanic:zsatfin",
    
    # Interactions with Other Race
    "Income × Other Race" = "Iraceother:zinc",
    "Perceived Relative Income × Other Race" = "Iraceother:zfinrela",
    "Financial Satisfaction × Other Race" = "Iraceother:zsatfin"
  ),
  statistics = c(
    "N (Observations)" = "nobs",
    "N (Respondents)" = "ngrps"
  )
)

# Export the table to a Word document
quick_docx(resultfe_raceall, file = "fe_raceall.docx")
###

# Extract coefficients and standard errors for each model separately
coef_white <- extract_coefficients(mfe_white)
coef_black <- extract_coefficients(mfe_black)
coef_hispanic <- extract_coefficients(mfe_hispanic)


# Add a column to indicate the model
coef_white$model <- "NH White"
coef_black$model <- "NH Black"
coef_hispanic$model <- "Hispanic"


# Combine the results
combined_coef_se_race <- rbind(coef_white, coef_black, coef_hispanic) |>
  filter(variable == "zinc" | variable == "zsatfin" | variable == "zfinrela" ) |>
  mutate(
    variable = str_replace_all(variable, "zinc", "income (logged)"),
    variable = str_replace_all(variable, "zfinrela", "perceived relative income"),
    variable = str_replace_all(variable, "zsatfin", "financial satisfaction"),
  )


plot_race <- ggplot(combined_coef_se_race, aes(x = fct_rev(fct_relevel(variable, "income (logged)")), y = coef, color = model, shape = model)) +
  geom_hline(yintercept=0, lwd=1, colour="grey50") +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), width = 0.2, linewidth = .8) +
  theme_minimal() +
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Variable", y = "Coefficient (Standardized)") + # , title = "Fixed Effects Coefficients by Race/Ethnicity") +
  coord_flip() +
  scale_color_manual(values = c("NH White" = "#1F77B4", "NH Black" = "#FF7F0E", "Hispanic" = "#2CA02C"), breaks = c("NH White", "NH Black", "Hispanic")) + 
  scale_shape_manual(values = c("NH White" = 16, "NH Black" = 17, "Hispanic" = 15), breaks = c("NH White", "NH Black", "Hispanic")) 

plot_race

ggsave("fe_race.pdf", width = 8, height = 4, units = "in")


#### Models by Sex #### 
# By sex
gss_panel_long2_sex <- gss_panel_long2_complete |>
  filter(sex_change == FALSE)


gss_male <- gss_panel_long2_sex |>
  filter(Isex == "male")

mfe_male <- feols(zincom16 ~ zinc + zfinrela + zsatfin + marital_f+childs+attend + year_f|id, data = gss_male, weights=~wtpannr12)

gss_female <- gss_panel_long2_sex |>
  filter(Isex == "female")

mfe_female <- feols(zincom16 ~ zinc + zfinrela + zsatfin + marital_f+childs+attend + year_f|id, data = gss_female, weights=~wtpannr12)

huxreg(mfe_male, mfe_female)

mfe_allsex <- feols(zincom16 ~ (zinc + zfinrela + zsatfin + marital_f+childs+attend + year_f)*Isex|id, data = gss_panel_long2_sex, weights=~wtpannr12)
huxreg(mfe_allsex)


# add the number of respondents (level 2) to each model
mfe_male2 <- add_ngrps_to_model(mfe_male)
mfe_female2 <- add_ngrps_to_model(mfe_female)
mfe_allsex2 <- add_ngrps_to_model(mfe_allsex)


resultfe_sex <- huxreg("Male" = mfe_male2, "Female" = mfe_female2,
                        coefs = c("income (logged, standardized)" = "zinc",
                                  "financial satisfaction (standardized)" = "zsatfin",
                                  "perceived relative income (standardized)" = "zfinrela",
                                  "number of children" = "childs",
                                  "religious service attendance" = "attend",
                                  "marital status: widowed" = "marital_fwidowed",
                                  "marital status: divorced" = "marital_fdivorced",
                                  "marital status: separated" = "marital_fseparated",
                                  "marital status: never married" = "marital_fnever married"
                        ),
                       statistics = c(
                         "N (Observations)" = "nobs",
                         "N (Respondents)" = "ngrps"
                       )
)

quick_docx(resultfe_sex, file = "fe_sex.docx")


## interaction model
# Get the exact names of the coefficients in the model
# coef_names <- names(coef(mfe_allsex2))

# Create the huxreg table for mfe_raceall
resultfe_sexall <- huxreg(
  mfe_allsex2,
  coefs = c(
    # Main effects
    "Income (standardized)" = "zinc",
    "Perceived Relative Income (standardized)" = "zfinrela",
    "Financial Satisfaction (standardized)" = "zsatfin",
    
    # Interactions with Female
    "Income × Female" = "zinc:Isexfemale",
    "Financial Satisfaction × Female" = "zsatfin:Isexfemale", 
    "Perceived Relative Income × Female" = "zfinrela:Isexfemale"
    
  ),
  statistics = c(
    "N (Observations)" = "nobs",
    "N (Respondents)" = "ngrps"
  )
)

# Export the table to a Word document
quick_docx(resultfe_sexall, file = "fe_sexall.docx")


# Extract coefficients and standard errors for each model separately
coef_male <- extract_coefficients(mfe_male)
coef_female <- extract_coefficients(mfe_female)


# Add a column to indicate the model
coef_male$model <- "Male"
coef_female$model <- "Female"


# Combine the results
combined_coef_se_sex <- rbind(coef_female, coef_male) |>
  filter(variable == "zinc" | variable == "zsatfin" | variable == "zfinrela" ) |>
  mutate(
    model = factor(model, levels = c("Male", "Female")), # Reorder model to display Female first
    variable = str_replace_all(variable, "zinc", "income (logged)"),
    variable = str_replace_all(variable, "zfinrela", "perceived relative income"),
    variable = str_replace_all(variable, "zsatfin", "financial satisfaction"),
  )

plot_sex <- ggplot(combined_coef_se_sex, aes(x = fct_rev(fct_relevel(variable, "income (logged)")), y = coef, color = model, shape = model)) +
  geom_hline(yintercept=0, lwd=1, colour="grey50") +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5), width = 0.2, linewidth = .8) +
  theme_minimal() +
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Variable", y = "Coefficient (Standardized)") + #, title = "Fixed Effects Coefficients by Sex") +
  coord_flip() +
  scale_color_manual(values = c("Female" = "#1F77B4", "Male" = "#FF7F0E"), breaks = c("Female", "Male")) + 
  scale_shape_manual(values = c("Female" = 16, "Male" = 17), breaks = c("Female", "Male"))

plot_sex
ggsave("fe_sex.pdf", plot_sex, width = 8, height = 4, units = "in")

# If combinging sex and race into one plot
# plot_sex / plot_race

# combined_plot <- plot_sex / plot_spacer() / plot_race + 
#   plot_layout(heights = c(1, 0.2, 1)) # Adjust spacer height here
# 
# combined_plot

# ggsave("fe_het.pdf", width = 6, height = 6, units = "in")

