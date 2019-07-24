##### Metabolimic Individual Analysis #####

# Load/install the required packages for tidyverse needed libxml2-dev, 
# libssl-dev, libcurl4-openssl-dev in the OS
install_github("kassambara/easyGgplot2", force = T)

GetPackages <- function(required.packages) {
  packages.not.installed <- 
    required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){
    install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = TRUE))
}

GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "reshape2", "devtools", "PerformanceAnalytics", 
              "easyGgplot2", "ggplot2", "car", "ggpubr", "lubridate", "bootstrap", "randomForestExplainer", "corrplot", "ggraph"))

# Import the data
twins_data_bmi <- read.delim("~/Downloads/E957-_Data_transfer_15%2f02%2f2019_%28Height_and_Weight%29/bmi.csv", sep = ",")
twins_data_demographics <- read.delim("~/Downloads/E957-_Data_transfer_15%2f02%2f2019_%28Height_and_Weight%29/demographics.csv", sep = ",")
twins_data_metabolomics <- read.delim("~/Downloads/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse/metabolon_2015_scaleRunDayMedian_imputeRunDayMin_normInverse.txt")

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

# Merge information into metabolomics
twins_data_met_dob <- merge(twins_data_metabolomics, twins_dob, by = "PublicID")
twins_data_met_dob$age_at_test <- as.numeric(twins_data_met_dob$year_hli) - as.numeric(twins_data_met_dob$birth_year)
drops <- c("PublicID", "month_year_birth.x", "age_hli", "sex", "zygosity", "bmi_hli", "date_hli", "year_hli", 
           "visit_hli", "batch_hli", "Year.Reported", "Month.Reported", "month_year_birth", "birth_year")
twins_data_met_dob <- twins_data_met_dob[ , !(names(twins_data_met_dob) %in% drops)]

##### Metabolomics Regression #####

##### Log Rank #####
met <- melt(twins_data_metabolomics[c(11:766)])

# Multiple Linear Regression Example 

# Trying to get this to work
# fit_met <- lm(age_at_test ~ c(names(twins_data_met_dob)[10:751]), data=twins_data_met_dob)
# (needs to include the + symbol)

fit_met <- lm(age_at_test ~ 
                m52701 + m44621 + m52689 + m52673 + m52630 + m52672 + m52682 + m52677 + m52478 + m52477 + m52716 + m52474 + m39270 + m52613 + m52475 + m52704 + m52476 + m52612 + m52614 + m52702
              + m39271 + m48762 + m52603 + m18985 + m19130 + m34404 + m32391 + m20675 + m34400 + m33971 + m33972 + m32497 + m37536 + m37752 + m38768 + m38168 + m39609 + m38296 + m52698 + m46325
              + m33228 + m35186 + m34214 + m34397 + m45456 + m33821 + m44630 + m48341 + m33871 + m35153 + m33822 + m44633 + m37231 + m45675 + m44563 + m47888 + m44560 + m34393 + m45951 + m52710
              + m27447 + m52690 + m34419 + m36600 + m36594 + m52500 + m52499 + m34391 + m44682 + m37419 + m30460 + m32350 + m27665 + m34395 + m34389 + m52707 + m19258 + m35625 + m45453 + m47087
              + m52705 + m52697 + m52711 + m46798 + m52453 + m46799 + m21184 + m48258 + m35628 + m36602 + m52683 + m52632 + m52431 + m33230 + m52706 + m52462 + m52464 + m52467 + m52454 + m52610
              + m52465 + m52463 + m52684 + m52633 + m42446 + m42449 + m52450 + m52737 + m52461 + m19263 + m52669 + m52470 + m52616 + m52634 + m21127 + m33955 + m35631 + m35305 + m20458 + m52497
              + m52498 + m37418 + m42450 + m52447 + m52449 + m52629 + m52611 + m52466 + m52699 + m52700 + m52452 + m52446 + m52468 + m52738 + m52438 + m42448 + m33961 + m42398 + m19324 + m52703
              + m38276 + m46115 + m6146  + m42374 + m43761 + m43343 + m43266 + m19266 + m33387 + m32815 + m52602 + m46203 + m36746 + m52281 + m42489 + m52294 + m18281 + m22036 + m35675 + m1432 
              + m17945 + m15667 + m32506 + m35257 + m36593 + m45095 + m31928 + m52282 + m21232 + m48259 + m45455 + m47118 + m33419 + m35253 + m45452 + m43400 + m31675 + m41220 + m39223 + m35635
              + m32197 + m34399 + m1566  + m31787 + m32397 + m531   + m542   + m22053 + m39600 + m1549  + m32457 + m22001 + m48448 + m31943 + m27672 + m48763 + m12017 + m46165 + m44526 + m15676
              + m36749 + m38667 + m46547 + m46548 + m15677 + m32445 + m15749 + m3155  + m1558  + m37181 + m37209 + m37202 + m37203 + m37211 + m37210 + m36099 + m15681 + m48441 + m35527 + m541  
              + m1669  + m32349 + m22116 + m46146 + m36098 + m1418  + m34424 + m37190 + m37198 + m37200 + m37196 + m33968 + m46297 + m31938 + m437   + m15685 + m1419  + m35136 + m1494  + m43231
              + m36776 + m35114 + m34390 + m32198 + m554   + m32980 + m1126  + m1107  + m41530 + m22132 + m46537 + m528   + m1561  + m31591 + m4970  + m48885 + m48255 + m1118  + m1110  + m1638 
              + m512   + m443   + m22175 + m18362 + m48492 + m15778 + m55    + m12129 + m3141  + m38100 + m32586 + m43807 + m2137  + m32412 + m569   + m1642  + m32489 + m32492 + m15500 + m35320
              + m48782 + m1563  + m22842 + m63    + m15506 + m34396 + m38637 + m38178 + m52304 + m47076 + m1564  + m2132  + m1712  + m1769  + m27718 + m513   + m37104 + m15705 + m1868  + m22176
              + m37443 + m56    + m514   + m33941 + m32425 + m36747 + m1114  + m40700 + m17805 + m35718 + m601   + m36808 + m5086  + m32415 + m44675 + m32504 + m32388 + m48407 + m31548 + m18467
              + m33587 + m33973 + m37459 + m1552  + m20699 + m42420 + m15765 + m47112 + m48195 + m1643  + m27719 + m44876 + m37092 + m37063 + m33934 + m36738 + m2730  + m33949 + m18245 + m34456
              + m18369 + m44872 + m33422 + m33364 + m2734  + m43829 + m52473 + m18280 + m587   + m48152 + m15443 + m57    + m53    + m396   + m44664 + m1572  + m15122 + m43847 + m37455 + m47155
              + m15990 + m58    + m32346 + m18476 + m32599 + m18477 + m32620 + m39379 + m21029 + m33954 + m18357 + m43802 + m32446 + m1573  + m46957 + m1644  + m35678 + m32328 + m35436 + m15753
              + m59    + m41912 + m22137 + m22138 + m40473 + m35322 + m43264 + m590   + m3127  + m15716 + m40730 + m27513 + m18349 + m32405 + m1123  + m33441 + m35437 + m1125  + m40049 + m44656
              + m34407 + m35107 + m1417  + m15140 + m527   + m1645  + m34534 + m60    + m40010 + m36756 + m1105  + m34035 + m46223 + m1301  + m1303  + m20676 + m15872 + m46142 + m48153 + m1121 
              + m1302  + m44878 + m18374 + m46144 + m1584  + m48429 + m15745 + m1124  + m1365  + m32418 + m48182 + m33952 + m31536 + m40469 + m15650 + m35137 + m1498  + m36752 + m35157 + m36713
              + m1585  + m33942 + m22185 + m37432 + m43488 + m48434 + m15720 + m33943 + m27710 + m33946 + m33967 + m1587  + m1589  + m32377 + m33950 + m37496 + m37076 + m48187 + m33959 + m32390
              + m1591  + m43249 + m2829  + m48433 + m594   + m37431 + m1356  + m42092 + m44877 + m36754 + m33936 + m52285 + m38102 + m35160 + m46111 + m1493  + m1505  + m35172 + m45413 + m20694
              + m1336  + m33447 + m52434 + m38165 + m37506 + m44681 + m1508  + m18254 + m36103 + m48841 + m32553 + m15958 + m35126 + m64    + m40192 + m38150 + m40016 + m41377 + m22130 + m566  
              + m42109 + m52719 + m52720 + m15704 + m1444  + m33935 + m32619 + m40708 + m32562 + m38170 + m35127 + m1898  + m40703 + m40731 + m32452 + m31932 + m33442 + m31555 + m46225 + m22194
              + m31522 + m48428 + m48990 + m18335 + m1899  + m1806  + m15772 + m27731 + m1515  + m32398 + m1648  + m2342  + m42049 + m39592 + m52605 + m42463 + m52433 + m52615 + m37529 + m48490
              + m48491 + m52495 + m47153 + m47154 + m42459 + m52435 + m52437 + m17747 + m34445 + m34384 + m1358  + m33969 + m19503 + m34409 + m15730 + m1437  + m41888 + m37058 + m46960 + m15336
              + m20693 + m2125  + m18494 + m32807 + m12261 + m36850 + m35669 + m22206 + m18392 + m18394 + m27738 + m1284  + m604   + m36095 + m2761  + m35428 + m32306 + m607   + m32401 + m40406
              + m54    + m37097 + m1299  + m42395 + m1604  + m1670  + m606   + m1649  + m42072 + m1567  + m3147  + m46500 + m46497 + m46355 + m48159 + m47478 + m32578 + m46592 + m46593 + m46390
              + m49613 + m46594 + m32698 + m46597 + m46598 + m46599 + m46360 + m46443 + m46460 + m46601 + m46602 + m47651 + m46914 + m46603 + m47594 + m46604 + m32857 + m46260 + m33132 + m33140
              + m46509 + m46466 + m46507 + m47641 + m46608 + m46367 + m47905 + m46588 + m47687 + m47642 + m47702 + m47909 + m46612 + m46358 + m46613 + m46614 + m47709 + m46616 + m46618 + m46619
              + m46620 + m46628 + m46458 + m47201 + m47715 + m47716 + m47718 + m48217 + m46894 + m46622 + m46624 + m46997 + m46471 + m46633 + m46634 + m46364 + m46636 + m47788 + m46643 + m47441
              + m46646 + m46417 + m46632 + m47225 + m46649 + m47787 + m46657 + m46661 + m46640 + m46266 + m46662 + m46473 + m46977 + m47800 + m47801 + m46668 + m46644 + m46697 + m47013 + m47804
              + m46681 + m46511 + m46683 + m46517 + m47993 + m46689 + m46691 + m46692 + m46693 + m46406 + m46698 + m46490 + m46412 + m40097 + m46672 + m47000 + m47821 + m48001 + m46366 + m47670
              + m46690 + m46695 + m46368 + m46701 + m46259 + m46295 + m46363 + m46384 + m46388 + m46398 + m46410 + m46486 + m46489 + m46493 + m46694 + m46696 + m46740 + m46898 + m46904 + m46905
              + m46972 + m46987 + m47026 + m47072 + m47417 + m47640 + m48067 + m48076 + m48081 + m48160 + m48168 + m48537 + m48992 + m49704 + m49228 + m49463 + m49469 + m49515 + m49521 + m49583
              + m49663 + m49681 + m49883 + m52050 + m52090 + m52094 + m52126 + m52154 + m52278 + m52280 + m52284 + m52413 + m52297 + m52504 + m52655 + m48312, data=twins_data_met_dob)

summary(fit_met) # show results

# Plotting some individual graphs
plot(age_at_test ~ m52682, data=twins_data_met_dob)
fit1 <- lm(age_at_test ~ m52682, data=twins_data_met_dob)
abline(fit1)
plot(age_at_test ~ m34400, data=twins_data_met_dob)
fit1 <- lm(age_at_test ~ m34400, data=twins_data_met_dob)
abline(fit1)

# Other useful functions 
coefficients(fit_met) # model coefficients
confint(fit_met, level=0.95) # CIs for model parameters 
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
results_met <- crossval(x, y, theta.fit.met, theta.predict.met, ngroup=10)
results_met <- crossval(x, y, theta.fit.met, theta.predict.met, ngroup=700)

cor(y, fit_met$fitted.values)**2 # raw R^2 
cor(y, results_met$cv.fit)**2 # cross-validated R2

##### Random Forest Metabolomics #####

# Melt the data  for metabolomics
met <- melt(twins_data_metabolomics[c(11:766)])
# Can remove certain values if required (metnozero <- met[ which( met$value > 0 | met$value < 0) , ])

# Line histograms method
ggplot(met) + 
  geom_line(aes(colour= variable, x=value), stat="density", size=0.25) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(
          size = .1, color = "black"),
        legend.position="none")  +
  ylim(0, 0.6) +
  xlim(-5, 5)

# Plot a qq graph for further visual inspection of the data
ggqqplot(met$value)
qqPlot(met$value)
hist(met$value) # Looks normal

# Overlapping density histograms methods
ggplot(met, aes(x=value, fill=variable)) + theme_minimal() + 
  geom_density(alpha = 0.1, binwidth=0.05) +
  guides(fill=FALSE) +
  scale_color_brewer(palette="Dark2") +
  ylim(0, 0.6) +
  xlim(-5, 5)

# Testing for normality
ks.test(met$value, y = pnorm)