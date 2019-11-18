##### Install required libraries

GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "rsample",
              "reshape2", "devtools", "PerformanceAnalytics", "ggplot2", "lubridate",
              "bootstrap", "corrplot", "ggraph", "doParallel", "ranger", "data.table", "h2o",
              "sparsio"))

twin1_1 <- as.data.frame(Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_train_Twin1*"))
twin1_1$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_train_Twin1*")` <- as.character(twin1_1$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_train_Twin1*")`)
colnames(twin1_1) <- c("Twin1")
twin1_2 <- as.data.frame(Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_train_Twin1*"))
twin1_2$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_train_Twin1*")` <- as.character(twin1_2$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_train_Twin1*")`)
colnames(twin1_2) <- c("Twin1")
twin1_3 <- as.data.frame(Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_train_Twin1*"))
twin1_3$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_train_Twin1*")` <- as.character(twin1_3$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_train_Twin1*")`)
colnames(twin1_3) <- c("Twin1")
twin_1 <- as.data.frame(rbind(twin1_1, twin1_2, twin1_3))

twin2_1 <- as.data.frame(Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_test_Twin2*"))
twin2_1$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_test_Twin2*")` <- as.character(twin2_1$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*1-4_SNP_merged_test_Twin2*")`)
colnames(twin2_1) <- c("Twin2")
twin2_2 <- as.data.frame(Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_test_Twin2*"))
twin2_2$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_test_Twin2*")` <- as.character(twin2_2$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*5-10_SNP_merged_test_Twin2*")`)
colnames(twin2_2) <- c("Twin2")
twin2_3 <- as.data.frame(Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_test_Twin2*"))
twin2_3$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_test_Twin2*")` <- as.character(twin2_3$`Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*11-22_SNP_merged_test_Twin2*")`)
colnames(twin2_3) <- c("Twin2")
twin_2 <- as.data.frame(rbind(twin2_1, twin2_2, twin2_3))

files <- as.data.frame(Sys.glob(paths = "~/Documents/twins_ML_project/plink/ML_project/*top100*"))
files$Twin1 <- twin_1
files$Twin2 <- twin_2
colnames(files) <- c("top100", "Twin1", "Twin2")
files$top100 <- as.character(files$top100)
rm(list = ls(pattern = "^twin")) 

# rs2071348

SNPselector <- function(files){
  top_100 <- read.csv(file = files[,1][1], header = T, stringsAsFactors = F)
  top_100[nrow(top_100) + 1,] <- c("geom_mean_FCFACS")
  SNPs_and_FACS_Twins1 <- readRDS(file = files[,2][1,])
  SNPs_and_FACS_Twins2 <- readRDS(file = files[,3][1,])
  SNPs_and_FACS_train_Twin1_filtered <<- subset(SNPs_and_FACS_Twins1, select = top_100[,1])
  SNPs_and_FACS_test_Twin2_filtered <<- subset(SNPs_and_FACS_Twins2, select = top_100[,1])
  for(i in 2:9){
    top_100 <- read.csv(file = files[,1][i], header = T, stringsAsFactors = F)
    SNPs_and_FACS_Twins1 <- readRDS(file = files[,2][i,])
    SNPs_and_FACS_Twins2 <- readRDS(file = files[,3][i,])
    SNPs_and_FACS_train_Twin1_filtered_i <- subset(SNPs_and_FACS_Twins1, select = top_100[,1])
    SNPs_and_FACS_train_Twin1_filtered <<- cbind(SNPs_and_FACS_train_Twin1_filtered, SNPs_and_FACS_train_Twin1_filtered_i)
    SNPs_and_FACS_test_Twin2_filtered_i <- subset(SNPs_and_FACS_Twins2, select = top_100[,1])
    SNPs_and_FACS_test_Twin2_filtered <<- cbind(SNPs_and_FACS_test_Twin2_filtered, SNPs_and_FACS_test_Twin2_filtered_i)
  }
}

FACS_data <- rbind(SNPs_and_FACS_train_Twin1_filtered, SNPs_and_FACS_test_Twin2_filtered)

output <- ggplot(FACS_data, aes(x = geom_mean_FCFACS)) +
  geom_histogram(color = "darkblue", fill = "lightblue", binwidth = 1) +
  theme(
    # Centre title
    plot.title = element_text(hjust = 0.5),
    # Lengends to the top
    legend.position = "top",
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank())

ggsave("~/Dropbox/STP/MSc Bioinformatics/MSc Project/FACS_Fcell_histogram.pdf", output, width = 16*0.75, height = 9*0.75)

SNPselector(files)

rm(list = c("SNPs_and_FACS_Twins1", "SNPs_and_FACS_Twins2", "top_100", "files"))
gc()

saveRDS(SNPs_and_FACS_train_Twin1_filtered, file = "~/Documents/twins_ML_project/plink/ML_project/SNPs_and_FACS_train_Twin1_filtered.rds")  # Save locally as compact rds file
saveRDS(SNPs_and_FACS_test_Twin2_filtered, file = "~/Documents/twins_ML_project/plink/ML_project/SNPs_and_FACS_test_Twin2_filtered.rds")  # Save locally as compact rds file

SNPs_and_FACS_train_Twin1_filtered.rds <- readRDS(file = "~/Documents/twins_ML_project/plink/ML_project/SNPs_and_FACS_train_Twin1_filtered.rds")
SNPs_and_FACS_test_Twin2_filtered.rds <- readRDS(file = "~/Documents/twins_ML_project/plink/ML_project/SNPs_and_FACS_test_Twin2_filtered.rds")

m2 <- randomForest::tuneRF(
  x = SNPs_and_FACS_train_Twin1_filtered,
  y = SNPs_and_FACS_train_Twin1_filtered$geom_mean_FCFACS,
  ntreeTry = 500,
  mtryStart = 100,
  stepFactor = 1.25,
  improve = 0.01,
  trace = T  # show real-time progress
)

saveRDS(m2, file = "~/Documents/twins_ML_project/plink/ML_project/mtry_testing_900SNPs.rds")

# Split the glycomics data into training and test sets
# Twins1_SNPs_fcell_split <- rsample::initial_split(SNP_merged_test_Twin1, prop = 0.8)
# Twins1_SNPs_fcell_train <- rsample::training(Twins1_SNPs_fcell_split)
# Twins1_SNPs_fcell_test  <- rsample::testing(Twins1_SNPs_fcell_split)

Twins1_SNPs_fcell_train <- SNPs_and_FACS_train_Twin1_filtered
Twins1_SNPs_fcell_test <- SNPs_and_FACS_test_Twin2_filtered

rm(list = c("SNPs_and_FACS_train_Twin1_filtered", "SNPs_and_FACS_test_Twin2_filtered", "Twins1_SNPs_fcell_split", "m2"))
gc()

# Using h2o for grid optimisation and random forest generation, setup a cluster with 16Gb of RAM and all cores
h2o.init(max_mem_size = "16G", nthreads = -1)

# Use all cores
# registerDoParallel(parallel::detectCores())

# Turn training set into h2o object (slow)
# train.h2o <- as.h2o(Twins1_SNPs_fcell_train)
# test.h2o <- as.h2o(Twins1_SNPs_fcell_test)

# as.h2o very slo, use h2o.importFile instead (better to write to disk then load into h2o)
# data clipped? Maximum of ~265,000 columns?
data.table::fwrite(x = Twins1_SNPs_fcell_train, file = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins1_900_SNPs_fcell_train.csv")
data.table::fwrite(x = Twins1_SNPs_fcell_test, file = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins2_900_SNPs_fcell_test.csv")

train.h2o <- h2o.importFile(path = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins1_900_SNPs_fcell_train.csv")
test.h2o <- h2o.importFile(path = "/home/callumrakhit/Documents/twins_ML_project/plink/ML_project/Twins2_900_SNPs_fcell_test.csv")

rm(Twins1_SNPs_fcell_test)
rm(Twins1_SNPs_fcell_train)
gc()

# head(names(test.h2o))
# tail(names(test.h2o))  # Check if geom_mean_FCFACS is the last column

# Hyperparameter grid
hyper_grid.h2o <- list(
  ntrees = seq(101, 401, by = 50),
  mtries = seq(250, 500, by = 50),
  sample_rate = c(.55, .6, .65, .7)
)

# Crappy hyperparameter grid
hyper_grid.h2o <- list(
  ntrees = 601,
  mtries = 50000,
  sample_rate = c(0.6)
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
x <- setdiff(names(train.h2o), y)

# Need to enlarge the stack size in bash (ulimit -s 16384) - enlarge stack limit to 16Mb
ulimit::memory_limit(size = 16384)

# Identify stack limit
Cstack_info()["size"]

# Build the grid search 
grid_900 <- h2o.grid(
  algorithm = "randomForest",
  grid_id = "rf_grid_900",
  x = x, 
  y = y,  
  training_frame = train.h2o,
  hyper_params = hyper_grid.h2o,
  search_criteria = search_criteria
)

# Save the model locally
saveRDS(object = grid_900, file = "~/Documents/twins_ML_project/plink/ML_project/rf_900.rds")

# Load model
grid <- readRDS("~/Documents/twins_ML_project/plink/ML_project/rf_900.rds")

# Collect the results and sort by our model performance metric of choice
best_grid <- h2o.getGrid(
  grid_id = "rf_grid_900", 
  sort_by = "mse",
  decreasing = F
)

# Print results
print(best_grid)

# Try to get some plots of the optimisation metrics
plot(best_grid@summary_table[["mtries"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["ntrees"]], best_grid@summary_table[["mse"]])
plot(best_grid@summary_table[["sample_rate"]], best_grid@summary_table[["mse"]])

# Get the model_id for the top model, chosen by validation error
best_model_id <- best_grid@model_ids[[1]]
best_model <- h2o.getModel(best_model_id)

# Evaluate the model performance on a test set
best_model_performance <- h2o.performance(model = best_model, newdata = test.h2o)

# RMSE of best model
h2o.mse(best_model_performance) %>% sqrt()
best_model_performance

# Make a prediction from the h2o random forest on the test set
pred_h2o <- predict(best_model, test.h2o)
head(pred_h2o)

best_model@model[["variable_importances"]]

# Plot 100 most important variables
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 50, wt = relative_importance) %>%
  ggplot(aes(reorder(variable, relative_importance), relative_importance)) +
  geom_col() +
  coord_flip() +
  ggtitle("The 50 variables which most reduce the OOB RMSE") +
  theme_minimal()

# Pick out the top 10,000 rsIDs (ranked by relative importance)
as.data.frame(best_model@model[["variable_importances"]]) %>%
  dplyr::arrange(desc(relative_importance)) %>%
  dplyr::top_n(n = 100, wt = relative_importance) %>% 
  {. ->> rf_1_top_100_predictors }

best_model@model[["variable_importances"]]

rf_1_top_100_predictors <- rf_1_top_100_predictors[,1]

View(rf_1_top_100_predictors)

rf_1_top_100_predictors == "rs2071348"

write_csv(x = as.data.frame(rf_1_top_100_predictors), path = "~/Documents/twins_ML_project/plink/ML_project/rf_900_top100.csv")

h2o.shutdown(prompt = F)  # Shuts down the Java cluster

install.packages("randomForestExplainer")
library(randomForestExplainer)

??randomForestExplainer

plot_multi_way_importance(best_model@model[["variable_importances"]][["variable"]], 
                          x_measure = best_model@model[["variable_importances"]][["scaled_importance"]], 
                          y_measure = best_model@model[["variable_importances"]][["relative_importance"]], 
                          size_measure = best_model@model[["variable_importances"]][["percentage"]], 
                          no_of_labels = 5)

install.packages("rsnps")

library(rsnps)
annotations(snp = 'rs766432', output = 'all')
ncbi_snp_query2(SNPs = "rs332")
ncbi_snp_query2("rs766432")
ncbi_snp_summary("rs420358")

rs766432
rs9399137

SNPs <- as.data.frame(gsub('.{2}$', '', noquote(best_model@model[["variable_importances"]][["variable"]])))
SNPs
View(SNPs)

install.packages("BiocManager")
BiocManager::install("biomaRt")
library(biomaRt)

snp_mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

snp_ids <- c("rs16828074", "rs17232800")
snp_ids <- SNPs$`gsub(".{2}$", "", noquote(best_model@model[["variable_importances"]][["variable"]]))`

snp_attributes <- c("refsnp_id", "chr_name", "chrom_start")

snp_locations <- getBM(attributes = snp_attributes, filters = "snp_filter", values = snp_ids, mart = snp_mart)

manhatten <- as.data.frame(snp_locations)

SNPs <- as.data.frame(gsub('.{2}$', '', noquote(best_model@model[["variable_importances"]][["variable"]])))
SNPs$variable.importance <- best_model@model[["variable_importances"]][["percentage"]]
SNPs$`gsub(".{2}$", "", noquote(best_model@model[["variable_importances"]][["variable"]]))` <- as.character(SNPs$`gsub(".{2}$", "", noquote(best_model@model[["variable_importances"]][["variable"]]))`)
colnames(SNPs) <- c("refsnp_id", "variable.importance")

SNPs <- merge(manhatten, y = SNPs, by = "refsnp_id")
SNPs <- SNPs[!is.na(as.numeric(as.character(SNPs$chr_name))),]
SNPs$chr_name <- as.numeric(SNPs$chr_name)
SNPs$chrom_start <- as.numeric(SNPs$chrom_start)

options(scipen = 999)

View(reorder(x = SNPs$chrom_start, X = SNPs$chr_name))

ggplot(SNPs, aes(x = reorder(x = SNPs$chrom_start, X = SNPs$chr_name), y = SNPs$variable.importance)) +
  
  # Show all points
  geom_point(aes(color = as.factor(SNPs$chr_name)), alpha = 0.8, size = 1.3) +
  # scale_color_manual(values = rep(c("grey", "skyblue"), 24)) +
  
  # custom X axis:
  # scale_x_continuous(label = SNPs$chr_name, breaks= SNPs$chr_name) +
  # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  # ylim(0, 9) +
  
  # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
)

class(SNPs$chr_name)
class(SNPs$chrom_start)
class(SNPs$variable.importance)

ggplot(SNPs) + 
  geom_point(aes(
    reorder(x = SNPs$chrom_start, X = SNPs$chr_name), 
    y = variable.importance, 
    color = as.factor(SNPs$chr_name)
    )) +
  # scale_x_continuous(label = SNPs$chr_name, breaks= SNPs$chr_name) +
  theme(
    # Lengends to the top
    legend.position = "none",
    # Remove the y-axis
    # axis.title.y = element_blank(),
    # Remove panel border
    panel.border = element_blank(),
    # Remove panel grid lines
    panel.grid.major.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size = .25, color = "black"),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank())

class(SNPs$chrom_start)

manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level = NA, 
                         annotate = NULL,
                         ann.default = list(),
                         should.thin = T, 
                         thin.pos.places = 2, 
                         thin.logp.places = 2, 
                         xlab = "Chromosome", 
                         # ylab = expression(-log[10](p-value)),
                         ylab = expression("Relative Importance in Random Forest Model"),
                         col = c("gray","darkgray"), 
                         # col = c("blue", "lightblue"),
                         panel.extra = NULL, pch = 20, cex = 0.8,...) {
  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  
  # make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  
  # make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  
  # calculate absolute genomic position from relative chromosomal positions
  posmin <- tapply(pos, chr, min);
  posmax <- tapply(pos, chr, max);
  posshift <- head(c(0, cumsum(posmax)), - 1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos <- function(cchr, cpos) {
    p <- posshift[as.character(cchr)] + cpos
    return(p)
  }
  
  # parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default <- list(x = "peak", y = "peak", adj = NULL, pos = 3, offset = 0.5, 
                      col = NULL, fontface = NULL, fontsize = NULL, show = F)
  parse.label <- function(rawval, groupname) {
    r <- list(text = groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times = length(pvalue)))
  }
  
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) + 1 ], 
                              fill=lfills[(i-2) %% length(lfills) + 1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  
  # reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      # logp = round(-log10(pvalue), thin.logp.places),
      logp = pvalue,
      pos = round(genpos, thin.pos.places), 
      chr = chr,
      grp = grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    # logp <- -log10(pvalue)
    logp <- pvalue
  }
  rm(pos, pvalue)
  gc()
  
  # custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  
  # make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A <- list();
    # maxy <- ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0))) + 0.5;
    maxy <- ceiling(max(y, ifelse(!is.na(sig.level), sig.level, 0))) - 0.85;
    A$ylim = c(0, maxy);
    A;
  }
  
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis = axis.chr, ann.settings=ann.settings, 
         prepanel = prepanel.chr, scales=list(axs="i"),
         panel = function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             # add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           # allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
  }

ann <- rep(1, length(SNPs$variable.importance))
ann[with(SNPs, chr_name == 2)] <- 2
ann[with(SNPs, chr_name == 6)] <- 3
ann[with(SNPs, chr_name == 11)] <- 4
ann <- factor(ann, levels = 1:4, labels=c("", "Chromsome 2", "Chromsome 6", "Chromsome 11"))

output <- manhattan.plot(factor(SNPs$chr_name, levels=c(1:22)), SNPs$chrom_start, SNPs$variable.importance, annotate = ann)
ggsave("~/Desktop/Manhatten.pdf", output, width = 16*1, height = 9*1)


View(ann)

