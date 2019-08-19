##### Load the packages and the data #####

# Load/install the required packages
# for tidyverse needed libxml2-dev, libssl-dev, libcurl4-openssl-dev in the OS
install_github("kassambara/easyGgplot2", force = T)

GetPackages <- function(required.packages) {
  packages.not.installed <- 
    required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){
    install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "reshape2", "devtools", "PerformanceAnalytics", 
              "easyGgplot2", "ggplot2", "car", "ggpubr", "lubridate", "bootstrap", "randomForestExplainer", "corrplot", "ggraph"))

# Import the data
twins_data_bmi <- read.delim("~/Downloads/E957-_Data_transfer_15%2f02%2f2019_%28Height_and_Weight%29/bmi.csv", sep = ",")
twins_data_demographics <- read.delim("~/Downloads/E957-_Data_transfer_15%2f02%2f2019_%28Height_and_Weight%29/demographics.csv", sep = ",")
twins_data_glycomics <- read.delim("~/Downloads/Glycomics/glycans.igg.global.combat.scale.processed.txt")
twins_data_metabolomics <- read.delim("~/Downloads/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse.txt")

##### Assessing the overlap in the datasets #####

count(
  merge(
    as.data.frame(sort(unique(twins_data_glycomics$PublicID))),
    as.data.frame(sort(unique(twins_data_metabolomics$PublicID))),
    by.x = 'sort(unique(twins_data_glycomics$PublicID))', 
    by.y = 'sort(unique(twins_data_metabolomics$PublicID))'
  )
)

##### Adding additional information to the datasets #####

# Add month/year DOB information to new data frame
twins_dob <- data.frame(matrix(vector(), dim(twins_data_metabolomics)[1], 2, dimnames=list(
  c(), c("PublicID", "birth_year"))), stringsAsFactors=F)
twins_dob$PublicID <- twins_data_metabolomics$PublicID
twins_dob$birth_year <- twins_data_metabolomics$month_year_birth
twins_dob$birth_year <- dmy(twins_dob$birth_year)
twins_dob$birth_year <- format(as.Date(twins_dob[["birth_year"]]), "%Y")

# Remove duplicates from the dob list
twins_dob <- twins_dob[!duplicated(twins_dob$PublicID), ]

# Merge information into glycomics
twins_data_gly_dob <- merge(twins_data_glycomics, twins_dob, by = "PublicID")
twins_data_gly_dob$date <- as.Date(twins_data_gly_dob$date, format = "%d/%m/%Y")
twins_data_gly_dob$age_at_test <- as.numeric(format(as.Date(twins_data_gly_dob[["date"]]), "%Y")) - as.numeric(twins_data_gly_dob$birth_year)
drops <- c("PublicID", "batch", "plate", "date", "Year.Reported", "Month.Reported", "month_year_birth", "birth_year", "ID_visityear")
twins_data_gly_dob <- twins_data_gly_dob[ , !(names(twins_data_gly_dob) %in% drops)]

##### Glycomics Regression #####

# Melt the data for glycomics
gly <- melt(twins_data_glycomics[c(5:80)])

##### Log Rank #####

# names(twins_data_gly_dob)[1:76]
# noquote(colnames(twins_data_gly_dob))
# colnames(twins_data_gly_dob)[1:76]

# Multiple Linear Regression Example 
fit <- lm(age_at_test ~ GP1 + GP2 + GP4 + GP5 + GP6 + GP7 + GP8 + GP9 + GP10 + GP11 + GP12 + GP13 + GP14 + GP15 + 
            GP16 + GP17 + GP18 + GP19 + GP22 + GP23 + GP24 + GP2021 + IGP24 + IGP25 + IGP26 + IGP27 + IGP28 + 
            IGP29 + IGP30 + IGP31 + IGP32 + IGP33 + IGP34 + IGP35 + IGP36 + IGP37 + IGP38 + IGP39 + IGP40 + 
            IGP41 + IGP42 + IGP43 + IGP44 + IGP45 + IGP46 + IGP47 + IGP48 + IGP49 + IGP50 + IGP51 + IGP52 + 
            IGP53 + IGP54 + IGP55 + IGP56 + IGP57 + IGP58 + IGP59 + IGP60 + IGP61 + IGP62 + IGP63 + IGP64 + 
            IGP65 + IGP66 + IGP67 + IGP68 + IGP69 + IGP70 + IGP71 + IGP72 + IGP73 + IGP74 + IGP75 + IGP76 + 
            IGP77, data=twins_data_gly_dob)

summary(fit) # show results

# Make some summary plots
unique(twins_data_glycomics$PublicID)
chart.Correlation(twins_data_gly_dob[70:77]) # subset of the data for purposes of visualisation

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
p.mat <- cor.mtest(twins_data_gly_dob)$p # p-values for the regressions
corrplot(cor(twins_data_gly_dob[60:77]), method = "color", col = col(200), # subset of the data for purposes of visualisation
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = FALSE)

# Plotting a couple of predictors for closer inspection
plot(age_at_test ~ IGP33, data = twins_data_gly_dob)
fit_IG33 <- lm(age_at_test ~ IGP33, data = twins_data_gly_dob)
abline(fit_IG33)
summary(fit_IG33)

plot(age_at_test ~ IGP55, data = twins_data_gly_dob)
fit_IG55 <- lm(age_at_test ~ IGP55, data = twins_data_gly_dob)
abline(fit_IG55)
summary(fit_IG55)

# Other useful functions 
coefficients(fit) # model coefficients
confint(fit, level = 0.95) # CIs for model parameters 
fitted(fit) # predicted values
residuals(fit) # residuals
anova(fit) # anova table 
vcov(fit) # covariance matrix for model parameters 
influence(fit) # regression diagnostics

# diagnostic plots 
layout(matrix(c(1:4), 2, 2)) # optional 4 graphs/page
plot(fit)

# Assessing R2 shrinkage using 10-Fold Cross-Validation

# matrix of predictors
x <- as.matrix(twins_data_gly_dob[c(
  "GP1", "GP2", "GP4", "GP5", "GP6", "GP7", "GP8", "GP9", "GP10", "GP11",       
  "GP12", "GP13", "GP14", "GP15", "GP16", "GP17", "GP18", "GP19", "GP22", "GP23",       
  "GP24", "GP2021", "IGP24", "IGP25", "IGP26", "IGP27", "IGP28", "IGP29", "IGP30", "IGP31", 
  "IGP32", "IGP33", "IGP34", "IGP35", "IGP36", "IGP37", "IGP38", "IGP39", "IGP40", "IGP41", 
  "IGP42", "IGP43", "IGP44", "IGP45", "IGP46", "IGP47", "IGP48", "IGP49", "IGP50", "IGP51", 
  "IGP52", "IGP53", "IGP54", "IGP55", "IGP56", "IGP57", "IGP58", "IGP59", "IGP60", "IGP61", 
  "IGP62", "IGP63", "IGP64", "IGP65", "IGP66", "IGP67", "IGP68", "IGP69", "IGP70", "IGP71", 
  "IGP72", "IGP73", "IGP74", "IGP75", "IGP76", "IGP77")])

# vector of predicted values
y <- as.matrix(twins_data_gly_dob[c("age_at_test")]) 

# define functions 
theta.fit <- function(x, y){lsfit(x, y)}
theta.predict <- function(fit, x){cbind(1, x) %*% fit$coef}

results <- crossval(x, y, theta.fit, theta.predict, ngroup=10)
cor(y, fit$fitted.values)**2 # raw R^2 
cor(y, results$cv.fit)**2 # cross-validated R2

##### Random Forest Glycomics #####

# Line histograms method
ggplot(gly) + 
  geom_line(aes(colour = variable, x = value), stat = "density", size = 0.25) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(
          size = .1, color = "black"),
        legend.position="none") +
  ylim(0, 0.6) +
  xlim(-5, 5)

# Plot a qq graph for further visual inspection of the data
ggqqplot(gly$value)
qqPlot(gly$value)
hist(gly$value) # Right skew

# Test for normality
ks.test(gly$value, y = pnorm)

# Testing "normality" of variables being used as predictors in a regression model
# before the fit is unwarranted. Does make sense to test the normality of the
# residuals since that is what is assumed in the modeling theory.

# Overlapping density histograms methods

# Glycomics
ggplot(gly, aes(x = value, fill = variable)) + theme_minimal() + 
  geom_density(alpha = 0.1, binwidth = 0.05) +
  guides(fill = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  xlim(-5, 5)

##### random forest modellling using caret #####

model_rf <- caret::train(YEAR_BIRTH ~ .,
                         data = twins_data2,
                         method = "rf",
                         preProcess = c("scale", "center"),
                         trControl = trainControl(method = "repeatedcv", 
                                                  number = 10, 
                                                  repeats = 10, 
                                                  savePredictions = TRUE, 
                                                  verboseIter = FALSE))

model_rf_smol_gly <- caret::train(age_at_test ~ .,
                                  data = twins_data_gly_dob,
                                  method = "rf",
                                  preProcess = c("scale", "center"),
                                  trControl = trainControl(method = "repeatedcv", 
                                                           number = 10, 
                                                           repeats = 10, 
                                                           savePredictions = TRUE, 
                                                           verboseIter = FALSE))

model_rf_smol_gly_ntree10 <- caret::train(as.factor(age_at_test) ~ .,
                                          data = twins_data_gly_dob,
                                          method = "rf",
                                          nodesize = 100,
                                          preProcess = c("scale", "center"),
                                          trControl = trainControl(method = "repeatedcv", 
                                                                   number = 10, 
                                                                   repeats = 10,
                                                                   savePredictions = TRUE, 
                                                                   verboseIter = FALSE))

# Plotting the random forest trees

tree_func <- function(final_model, tree_num) {
  # get tree by index
  tree <- randomForest::getTree(final_model, 
                                k = tree_num, 
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
  
  # prepare data frame for graph
  graph_frame <- data.frame(
    from = rep(tree$rowname, 2),
    to = c(tree$`left daughter`, tree$`right daughter`)
    )
  
  # convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")
  
  # set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- as.character(tree$prediction)
  V(graph)$split <- as.character(round(tree$`split point`, digits = 2))
  
  # plot
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), na.rm = TRUE, repel = TRUE) +
    geom_node_label(aes(label = split), vjust = 2.5, na.rm = TRUE, fill = "white") +
    geom_node_label(aes(label = leaf_label, fill = leaf_label), na.rm = TRUE, 
                    repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))
  
  # print plot
  print(plot)
}

tree_num <- which(model_rf$finalModel$forest$ndbigtree == min(model_rf$finalModel$forest$ndbigtree))
tree_func(final_model = model_rf$finalModel, tree_num)
getTree(model_rf$results, k = 1, labelVar=TRUE)

# Using just the glycomics data (min number of nodes)
tree_num <- which(model_rf_smol_gly$finalModel$forest$ndbigtree == min(model_rf_smol_gly$finalModel$forest$ndbigtree))
tree_func(final_model = model_rf_smol_gly$finalModel, tree_num)

# Using just the glycomics data (min number of nodes, 100 minimum per node)
tree_num <- which(model_rf_smol_gly_ntree10$finalModel$forest$ndbigtree == min(model_rf_smol_gly_ntree10$finalModel$forest$ndbigtree))[1]
tree_func(final_model = model_rf_smol_gly_ntree10$finalModel, tree_num)

# tree_num <- which(model_rf$finalModel$forest$ndbigtree == min(model_rf$finalModel$forest$ndbigtree))
# tree_func(final_model = model_rf$finalModel, tree_num)

##### Summarising the random forest results #####

# Age as a continuous numeric variable
rf_gly_1 <- randomForest(age_at_test ~ ., data = twins_data_gly_dob, importance = TRUE)
rf_met_1 <- randomForest(age_at_test ~ ., data = twins_data_met_dob, importance = TRUE)

# Age as a discrete factor
rf_gly_2 <- randomForest(as.factor(age_at_test) ~ ., data = twins_data_gly_dob, importance = TRUE, localImp = TRUE)
rf_met_2 <- randomForest(as.factor(age_at_test) ~ ., data = twins_data_met_dob, importance = TRUE) # n.a values replace with 0

# Make the html reports
explain_forest(rf_gly_1, interactions = TRUE, data = twins_data_gly_dob)
explain_forest(rf_gly_2, interactions = TRUE, data = twins_data_gly_dob)
explain_forest(rf_met_1, interactions = TRUE, data = twins_data_met_dob)
explain_forest(rf_met_2, interactions = TRUE, data = twins_data_met_dob)

randomForestExplainer::important_variables(rf_gly_2, k = 5)
