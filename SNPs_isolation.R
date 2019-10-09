##### Install required libraries

GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("data.table", "dplyr"))

##### Load in the genotype matrices #####

# Testing the SNP loading
SNP_list_test <- fread(input = "~/Documents/twins_ML_project/plink/temp", sep = " ")
SNP_list_test_noHET <- SNP_list_test[, grep("HET", colnames(SNP_list_test)):=NULL]
SNP_Fcell_test <- fread(input = "~/Twins/FID_ID_FCellFACS.csv")

SNP_merged_test <- merge(x = SNP_list_test_noHET, y = SNP_Fcell_test, by = "IID")
setnafill(SNP_merged_test, fill = 0)

SNP_merged_test <- SNP_merged_test[, grep("AT", colnames(SNP_merged_test)):=NULL]
SNP_merged_test <- SNP_merged_test[, grep("SEX", colnames(SNP_merged_test)):=NULL]
SNP_merged_test <- SNP_merged_test[, grep("PHENOTYPE", colnames(SNP_merged_test)):=NULL]

SNP_merged_test_Twin1 <- SNP_merged_test[!duplicated(SNP_merged_test[, 2])]
SNP_merged_test_Twin2 <- SNP_merged_test[duplicated(SNP_merged_test[, 2])]

SNP_merged_test_Twin1 <- SNP_merged_test_Twin1[, grep("ID", colnames(SNP_merged_test_Twin1)):=NULL]
SNP_merged_test_Twin2 <- SNP_merged_test_Twin2[, grep("ID", colnames(SNP_merged_test_Twin2)):=NULL]

class(SNP_merged_test_Twin1)

SNP_merged_test_Twin1[features]

m2 <- randomForest::tuneRF(
  x = SNP_merged_test_Twin1,
  y = SNP_merged_test_Twin1$geom_mean_FCFACS,
  ntreeTry = 500,
  mtryStart = 10,
  stepFactor = 1.5,
  improve = 0.01,
  trace = T  # show real-time progress 
)

m2[1:19]/5005

# Split the glycomics data into training and test sets
Twins1_SNPs_fcell_split <- rsample::initial_split(SNP_merged_test_Twin1, prop = 0.8)
Twins1_SNPs_fcell_train <- rsample::training(Twins1_SNPs_fcell_split)
Twins1_SNPs_fcell_test  <- rsample::testing(Twins1_SNPs_fcell_split)

# Using h2o for grid optimisation and random forest generation, setup a cluster with 32Gb of RAM
h2o.init(max_mem_size = "16G")

# Turn training set into h2o object
train.h2o <- as.h2o(Twins1_SNPs_fcell_train)

length(colnames(Twins1_SNPs_fcell_test))

ncol(SNP_merged_test_Twin1)

# Hyperparameter grid
hyper_grid.h2o <- list(
  ntrees = seq(501, 801, by = 100),
  mtries = seq((ncol(SNP_merged_test_Twin1)/5), (ncol(SNP_merged_test_Twin1)/5)+5000, by = 1000),
  sample_rate = c(.55, .6, .65, .7)
)

1645+10+40+30
1750-900
sqrt(5000)

# Random grid search criteria
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.005,
  stopping_rounds = 10,
  max_runtime_secs = 30*60
)

# Set values for x and y
y <- "geom_mean_FCFACS"
x <- setdiff(names(SNP_merged_test_Twin1), y)

# Build the grid search 
grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid",
  x = x, 
  y = y, 
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

# Collect the results and sort by our model performance metric of choice
best_grid <- h2o.getGrid(
  grid_id = "rf_grid", 
  sort_by = "mse",
  decreasing = F
)

# Print results
print(best_grid)

# Try to get some plots of the optimisation metrics
plot(best_grid@summary_table[["mtries"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["ntrees"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["sample_rate"]], best_grid@summary_table[["mse"]])

# Grab the model_id for the top model, chosen by validation error
best_model_id <- best_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Now let’s evaluate the model performance on a test set
test.h2o <- as.h2o(glycomics_fcell_test)
best_model_performance <- h2o.performance(model = best_model, newdata = test.h2o)

# RMSE of best model
h2o.mse(best_model_performance) %>% sqrt()
best_model_performance

# Make a prediction from the h2o random forest on the test set
pred_h2o <- predict(best_model, test.h2o)
head(pred_h2o)

# Plot 25 most important variables
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 76, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 25 variables which most reduce the OOB RMSE") +
  theme_minimal()

h2o.shutdown(prompt = F)  # Shuts down the Java cluster# Split the glycomics data into training and test sets
glycomics_fcell_split <- rsample::initial_split(twins_fcell_glycomics, prop = 0.8)
glycomics_fcell_train <- rsample::training(glycomics_fcell_split)
glycomics_fcell_test  <- rsample::testing(glycomics_fcell_split)

# Using h2o for grid optimisation and random forest generation, setup a cluster with 32Gb of RAM
h2o.init(max_mem_size = "16G")

# Turn training set into h2o object
train.h2o <- as.h2o(glycomics_fcell_train)

# Hyperparameter grid
hyper_grid.h2o <- list(
  ntrees = seq(200, 600, by = 100),
  mtries = seq(2, 5, by = 1),
  sample_rate = c(.55, .6, .65, .7)
)

# Random grid search criteria
search_criteria <- list(
  strategy = "RandomDiscrete",
  stopping_metric = "mse",
  stopping_tolerance = 0.005,
  stopping_rounds = 10,
  max_runtime_secs = 30*60
)

# Set values for x and y
y <- "geom_mean_FCFACS"
x <- setdiff(names(twins_fcell_glycomics), y)

# Build the grid search 
grid <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid",
  x = x, 
  y = y, 
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

# Collect the results and sort by our model performance metric of choice
best_grid <- h2o.getGrid(
  grid_id = "rf_grid", 
  sort_by = "mse",
  decreasing = F
)

# Print results
print(best_grid)

# Try to get some plots of the optimisation metrics
plot(best_grid@summary_table[["mtries"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["ntrees"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["sample_rate"]], best_grid@summary_table[["mse"]])

# Grab the model_id for the top model, chosen by validation error
best_model_id <- best_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Now let’s evaluate the model performance on a test set
test.h2o <- as.h2o(glycomics_fcell_test)
best_model_performance <- h2o.performance(model = best_model, newdata = test.h2o)

# RMSE of best model
h2o.mse(best_model_performance) %>% sqrt()
best_model_performance

# Make a prediction from the h2o random forest on the test set
pred_h2o <- predict(best_model, test.h2o)
head(pred_h2o)

# Plot 25 most important variables
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 76, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 25 variables which most reduce the OOB RMSE") +
  theme_minimal()

h2o.shutdown(prompt = F)  # Shuts down the Java cluster