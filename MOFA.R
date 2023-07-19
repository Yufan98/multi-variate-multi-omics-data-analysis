#-------------------------------Initialization----------------------------------
library(data.table)
library(MOFA2)
library(reticulate)
library(tidyr)
library(dplyr)
library(ggplot2)

#-------------------------------Import data-------------------------------------
## Read data generated from MOFA_Preprocess.py
proteomics <- read.csv("./protein4mofa.csv")
phosphoproteomics <- read.csv("./phosphos4mofa.csv")
transcriptomics <- read.csv("./genes4mofa.csv")
metadata <- read.csv("./datasets/metadata.csv")

## Rearrange
new_col <- c('sample', 'group', 'feature', 'view', 'value')
process_df <- function(df) {
  df <- df %>%
    gather(-c(feature,view), key=sample, value = value) %>%
    separate(col=sample, sep = "_", into=c("group","sample")) %>%
    select(new_col)
  return(df)
}

pro_single <- process_df(proteomics)
phos_single <- process_df(phosphoproteomics)
trans_single <- process_df(transcriptomics)

##joint
multi_all <-rbind.data.frame(pro_single, phos_single, trans_single)

#------------------------------Create MOFA--------------------------------------
MOFAobject <- create_mofa_from_df(multi_all,extract_metadata = TRUE)
print(MOFAobject)
metadata$Concentration <-as.factor(metadata$Concentration)
samples_metadata(MOFAobject) <- metadata
plot_data_overview(MOFAobject)
# Train model
## Prepare options
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 8
head(model_opts)
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
head(train_opts)
## Set options
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#------------------------------Train MOFA--------------------------------------
outfile = file.path(getwd(),"caffeine_model.hdf5")
## Train model
model <- run_mofa(MOFAobject, outfile)
#---------------------------------Results---------------------------------------
# Retrieve the data 
outfile = file.path(getwd(),"caffeine_model.hdf5")
model <- load_model(outfile)

plot_data_overview(model)
# Overall results
Nsamples = sum(model@dimensions$N)
colSums(model@cache[["variance_explained"]][["r2_per_factor"]][["T10"]])
colSums(model@cache[["variance_explained"]][["r2_per_factor"]][["T6"]])
colSums(model@cache[["variance_explained"]][["r2_per_factor"]][["T24"]])

model@cache[["variance_explained"]][["r2_per_factor"]][["T10"]]
model@cache[["variance_explained"]][["r2_per_factor"]][["T6"]]
model@cache[["variance_explained"]][["r2_per_factor"]][["T24"]]

plot_variance_explained(model, x="view", y="factor")
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
cor <- plot_factor_cor(model)
# Plot
plot_factor(model, groups = "all",
            factor = 1, color_by = "Concentration", shape_by = "group",
            scale = FALSE,dot_size = 4
)
plot_factor(model,
            factor = 2,color_by = "Concentration", shape_by = "group",
            scale = FALSE,dot_size = 4
)
plot_factor(model,
            factor = 3,color_by = "Concentration", shape_by = "group",
            scale = FALSE,dot_size = 4
)
plot_factor(model,
            factor = 4,color_by = "Concentration", shape_by = "group",
            scale = FALSE,dot_size = 4
)
plot_factor(model,
            factor = 5,color_by = "Concentration", shape_by = "group",
            scale = FALSE,dot_size = 4
)
plot_factor(model,
            factor = 6,color_by = "Concentration", shape_by = "group",
            scale = FALSE,dot_size = 4
)
plot_factor(model,
            factor = 7,color_by = "Concentration", shape_by = "group",
            scale = FALSE,dot_size = 4
)

# Exact data
factors <- get_factors(model, factors = "all")
weights <- get_weights(model, views = "all", factors = "all")
data <- get_data(model)

# Function to calculate p-values
calculate_p_values <- function(data, factors, view_name, time_point) {
  view_data <- as.data.frame(cbind(data[[view_name]][[time_point]], data[[view_name]]$T24))
  view_data <- view_data[sort(rownames(view_data)),]
  factors_data <- as.data.frame(rbind(factors[[time_point]], factors$T24))
  p_values <- as.data.frame(weights[[view_name]])
  p_values <- p_values[sort(rownames(p_values)),]
  
  for (j in 1:dim(factors_data)[2]) {
    for (i in 1:dim(view_data)[1]) {
      tryCatch({
        p_values[i, j] <- cor.test(unlist(view_data[i, ]), factors_data[, j])$p.value
      }, error=function(e){})
    }
  }
  
  p_values$corr <- 'na'
  
  for (k in 1:dim(p_values)[1]) {
    if (min(unlist(p_values[k, -8])) < 0.05) {
      p_values[k, 8] <- names(which.min(p_values[k, c(1:7)]))
    } else {
      p_values[k, 8] <- 'Not relative'
    }
  }
  
  return(p_values)
}

# Calculate p-values for each data type
p_values_phos <- calculate_p_values(data, factors, "Phosphoproteomics", "T16")
p_values_pro <- calculate_p_values(data, factors, "Proteomics", "T10")
p_values_trans <- calculate_p_values(data, factors, "Transcriptomics", "T16")

# Write results to CSV files
write.csv(p_values_phos, file = "mappiing_results_phos.csv")
write.csv(p_values_pro, file = "mappiing_results_pro.csv")
write.csv(p_values_trans, file = "mappiing_results_trans.csv")
