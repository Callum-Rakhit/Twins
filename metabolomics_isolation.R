##### Metabolimic Individual Analysis #####

# Load/install the required packages for tidyverse needed libxml2-dev, libssl-dev, libcurl4-openssl-dev in the OS
GetPackages <- function(required.packages) {
  packages.not.installed <- 
    required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){
    install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = TRUE))
}

GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "reshape2", "devtools",
              "PerformanceAnalytics", "ggplot2", "car", "ggpubr", "lubridate", "bootstrap", 
              "randomForestExplainer", "corrplot", "ggraph"))

install_github("kassambara/easyGgplot2", force = T)  # install_github is a devtools package, so load this first
library(easyGgplot2)

# Import the various datasets
# twins_data_bmi <- read.delim("~/Downloads/E957-_Data_transfer_15%2f02%2f2019_%28Height_and_Weight%29/bmi.csv", sep = ",")
twins_data_demographics <- read.delim("~/Downloads/E957-_Data_transfer_15%2f02%2f2019_%28Height_and_Weight%29/demographics.csv", sep = ",")
twins_data_metabolomics <- read.delim("~/Downloads/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse.txt")

##### Adding additional information to the datasets #####

# Extract year of birth information to new data frame
twins_data_demographics <- data.frame(matrix(vector(), dim(twins_data_metabolomics)[1], 2, dimnames = list(
  c(), c("PublicID", "birth_year"))), stringsAsFactors = F)
twins_data_demographics$PublicID <- twins_data_metabolomics$PublicID
twins_data_demographics$birth_year <- twins_data_metabolomics$month_year_birth
twins_data_demographics$birth_year <- dmy(twins_data_demographics$birth_year)
twins_data_demographics$birth_year <- format(as.Date(twins_data_demographics[["birth_year"]]), "%Y")

# Remove duplicates from the dob list
twins_data_demographics <- twins_data_demographics[!duplicated(twins_data_demographics$PublicID), ]

# Merge information into metabolomics
twins_data_metabolomics <- merge(twins_data_metabolomics, twins_data_demographics, by = "PublicID")
twins_data_metabolomics$age_at_test <- as.numeric(twins_data_metabolomics$year_hli) - as.numeric(twins_data_metabolomics$birth_year)
drops <- c("PublicID", "month_year_birth.x", "age_hli", "sex", "zygosity", "bmi_hli", "date_hli", "year_hli", 
           "visit_hli", "batch_hli", "Year.Reported", "Month.Reported", "month_year_birth", "birth_year")
twins_data_metabolomics <- twins_data_metabolomics[ , !(names(twins_data_metabolomics) %in% drops)]

##### Metabolomics Regression #####

# Multiple Linear Regression Example
metabolites <- paste(colnames(twins_data_metabolomics)[-length(twins_data_metabolomics)], collapse = " + ")

fit_met <- lm(paste("age_at_test ~ ", metabolites, sep = ""), 
        data = twins_data_metabolomics,
        na.action = na.exclude)

summary(fit_met) # show results

##### Playing around with ggplots for lm models #####

library(data.table)
devtools::install_github("cardiomoon/ggiraphExtra")
library(ggiraphExtra)

multiple_regression <- function(number_of_)

as.data.frame(a[["coefficients"]])
  
a <- summary(fit_met)
a <- as.data.frame(a$coefficients)
a <- setDT(a[order(a$`t value`),][1:4,], keep.rownames = TRUE)[]

postitives
negitives
a <- setDT(a[order(a$`Pr(>|t|)`),][1:4,], keep.rownames = TRUE)[]
b <- paste(a$rn[2:4], collapse = " + ")
a <- lm(paste("age_at_test ~ ", b, sep = ""), 
        data = twins_data_metabolomics,
        na.action = na.exclude)

par(mfrow = c(2,2))
plot(a)

ggPredict(a, se = T,  interactive = T)

d <- twins_data_metabolomics

# Obtain predicted and residual values
d$predicted <- predict(a)
d$residuals <- residuals(a)

# Create plot
d %>% 
  gather(key = "iv", value = "x", -age_at_test, -predicted, -residuals) %>%
  ggplot(aes(x = x, y = age_at_test)) +
  geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
  geom_point(aes(color = residuals)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  guides(color = FALSE) +
  geom_point(aes(y = predicted), shape = 1) +
  facet_grid(~ iv, scales = "free_x") +
  theme_minimal()

#####################################################

# Plotting some individual graphs
plot(age_at_test ~ m52682, data = twins_data_metabolomics)
fit1 <- lm(age_at_test ~ m52682, data = twins_data_metabolomics)
abline(fit1)
plot(age_at_test ~ m34400, data = twins_data_metabolomics)
fit1 <- lm(age_at_test ~ m34400, data = twins_data_metabolomics)
abline(fit1)

# Other useful functions 
coefficients(fit_met) # model coefficients
confint(fit_met, level = 0.95) # CIs for model parameters 
fitted(fit_met) # predicted values
residuals(fit_met) # residuals
anova(fit_met) # anova table 
vcov(fit_met) # covariance matrix for model parameters 
influence(fit_met) # regression diagnostics

# diagnostic plots 
layout(matrix(c(1:4), 2, 2)) # optional 4 graphs/page
plot(fit_met) # abline plots a straight line, so only takes two variables into account

##### Assessing R2 shrinkage using 10-Fold Cross-Validation  #####

# matrix of predictors
x <- as.matrix(twins_data_met_dob[c(names(twins_data_met_dob)[2:757])])

# vector of predicted values
y <- as.matrix(twins_data_met_dob[c("age_at_test")]) 

# define functions 
theta.fit.met <- function(x, y){lsfit(x, y)}
theta.predict.met <- function(fit_met, x){cbind(1, x) %*% fit_met$coef} 

# Doesn't work?
results_met <- crossval(x, y, theta.fit.met, theta.predict.met, ngroup = 10)
results_met <- crossval(x, y, theta.fit.met, theta.predict.met, ngroup = 700)

cor(y, fit_met$fitted.values)**2 # raw R^2 
cor(y, results_met$cv.fit)**2 # cross-validated R2

##### Random Forest Metabolomics #####

# Melt the data  for metabolomics
met <- melt(twins_data_metabolomics[c(11:766)])
# Can remove certain values if required (metnozero <- met[ which( met$value > 0 | met$value < 0) , ])

# Line histograms method
ggplot(met) + 
  geom_line(aes(colour = variable, x = value), stat = "density", size = 0.25) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(
          size = .1, color = "black"),
        legend.position = "none")  +
  ylim(0, 0.6) +
  xlim(-5, 5)

# Plot a qq graph for further visual inspection of the data
ggqqplot(met$value)
qqPlot(met$value)
hist(met$value) # Looks normal

# Overlapping density histograms methods
ggplot(met, aes(x = value, fill = variable)) + theme_minimal() + 
  geom_density(alpha = 0.1, binwidth = 0.05) +
  guides(fill = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  ylim(0, 0.6) +
  xlim(-5, 5)

# Testing for normality
ks.test(met$value, y = pnorm)