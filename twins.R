GetPackages <- function(required.packages) {
  packages.not.installed <- 
    required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){
    install.packages(packages.not.installed)}
  suppressMessages(lapply(required.packages, require, character.only = TRUE))
}

GetPackages(c("readxl", "ggplot2"))

library(ggplot2)

data <- read.delim("~/Twins_data/Glycomics/glycans.igg.global.combat.scale.processed.txt")

# Histogram of BMI

twins1 <- read_xlsx("~/Twins_data/E957_150219_PID.xlsx", col_names = T, sheet = 1)
names(twins1) = c("PublicID", "Study", "Birth_Year", "Sex", "Zygosity", "Ethnicity", "Anomaly")
twins2 <- read_xlsx("~/Twins_data/E957_150219_PID.xlsx", col_names = T, sheet = 2)
names(twins2) = c("PublicID", "Visit_Date", "Height_cm", "Weight_cm", "BMI")

data_merged <- merge(data, twins1, by = "PublicID")
data_merged <- merge(data_merged, twins2, by = "PublicID")
data_merged$age <- as.Date(data_merged$date) - as.Date(data_merged$Birth_Year)

data_merged$Birth_Year

data_merged$Test_Year <- format(as.Date(data_merged$date, format="%d/%m/%Y"),"%Y")

data_merged$Age_at_test <- as.numeric(data_merged$Test_Year) - data_merged$Birth_Year

class(data_merged$Birth_Year)

gglot(data_merged$Age_at_test, data_merged$GP1)
plot(data_merged$Age_at_test, data_merged$GP2)

ggplot(data_merged, aes(x = Age_at_test, y = GP4)) +
  geom_point() +  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(size = .1, color = "black")) +
  labs(title = "Correlation between age and glycome",
       y = "GP4 (glycome)", x = "Age") +
  stat_smooth(method = "lm", col = "red")

summary(lm(data_merged$Age_at_test ~ data_merged$GP4))

filtered <- twins[twins$BMI < 100,]
filtered <- filtered[filtered$BMI != 0,]

mean(filtered$BMI)
max(filtered$BMI)
hist(filtered$BMI, breaks = 20)