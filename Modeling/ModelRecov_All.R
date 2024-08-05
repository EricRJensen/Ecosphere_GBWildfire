library(tidyverse)
library(randomForest)
library(rfUtilities)
library(caret)
library(ggpubr)
library(raster)

# Function to clean dataframes
drop_cols = function(df){
  return (df %>% 
            dplyr::select(-c(`.geo`, Year, Fire_Name, `system:index`)) %>%
            mutate(ID = `paste(___)`) %>%
            dplyr::select(-`paste(___)`))}

# ---------------------- Read in CSVs ---------------------------------
# Topographic
setwd('C:/Users/erjensen/Documents/Thesis/ThesisAIM/data/Recov_CSVs/')
topo_vars2 <- read_csv('Biophysical_Predictors_Mode.csv') %>%
  mutate(landform = as.factor(mode)) %>%
  drop_cols() %>%
  dplyr::select(-mode)

topo_vars <- read_csv('Biophysical_Predictors_Mean.csv') %>%
  mutate(ID = `paste(___)`) %>%
  dplyr::select(-c(`.geo`, `system:index`,`paste(___)`)) %>% 
  left_join(topo_vars2, by = 'ID') %>%
  mutate(Fire_Name = as.factor(Fire_Name),
         Year = as.factor(Year))
remove(topo_vars2)

# Seeding
seed_vars <- read_csv('Seeding_fires.csv') %>%
  drop_cols() %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(Aerial = as.factor(Aerial),
         Drill = as.factor(Drill))

# Soils
soil_vars <- read_csv('Soil_Predictors_Fires.csv') %>%
  drop_cols()

# Species Richness
spec_vars <- read_csv('SpecRich_fires.csv') %>%
  drop_cols() %>%
  dplyr::select(-mean)

# Fire
fire_vars <- read_csv('Severity_fires.csv') %>%
  drop_cols()

# Climate averages                     
clia_files <- list.files('ClimateAverages/', full.names = TRUE, pattern = 'csv')
clia_vars = read_csv(clia_files[1], na=c("","NA")) %>% dplyr::select(-'mean')
for(i in clia_files){
  var_str <- str_split(i, "/")[[1]][2] %>% str_remove('_fires.csv')
  print(var_str)
  name = var_str
  indv_cli_df <- read_csv(i,na=c("","NA")) %>% dplyr::select(c(`paste(___)`, contains('mean')))
  indv_cli_df[, paste(name, '_avg', sep = '')] <- indv_cli_df$mean
  indv_cli_df <- dplyr::select(indv_cli_df,-mean)
  clia_vars <- left_join(clia_vars, indv_cli_df, by = 'paste(___)') %>% unique() }
remove(indv_cli_df, i, var_str, name, clia_files)
clia_vars = clia_vars %>% drop_cols()

# Prefire RAP
RAP_vars = read.csv('Recov_Response.csv')%>%
  mutate(ID = `paste.___.`) %>%
  dplyr::select(-`paste.___.`) %>%
  dplyr::select(c(ID, NDPDI_pre, AFGC_pre, SHR_pre, PFGC_pre, BG_pre))

# Post-fire climate
clip_files <- list.files('ClimatePostfire/', full.names = TRUE, pattern = 'csv')
clip_vars = read_csv(clip_files[1], na=c("","NA")) %>% drop_cols()
for(i in clip_files[2:length(clip_files)]){
  indv_cli_df <- read_csv(i,na=c("","NA")) %>% drop_cols() 
  clip_vars = rbind(clip_vars, indv_cli_df)}
var_names = names(clip_vars)[1:22] %>% word(2, sep = '_') %>% paste('_pst', sep = '')
var_names = c(var_names, 'ID')
names(clip_vars) = var_names
remove(i, clip_files, indv_cli_df,var_names)

# Merge all tables
predictors = topo_vars %>%
  left_join(seed_vars, by = 'ID') %>%
  left_join(soil_vars, by = 'ID') %>%
  left_join(clia_vars, by = 'ID') %>%
  left_join(clip_vars, by = 'ID') %>%
  left_join(fire_vars, by = 'ID') %>%
  left_join(spec_vars, by = 'ID') %>%
  left_join(RAP_vars, by = 'ID') %>%
  dplyr::select(ID,everything()) %>%
  drop_na()
names(predictors)[2:ncol(predictors)] <-paste('p.', names(predictors)[2:ncol(predictors)], sep = '')

# Read in response variables
responses = read.csv('Recov_Response.csv') %>%
  mutate(ID = `paste.___.`) %>%
  dplyr::select(c(ID, NDPDI_post))
names(responses)[2:2] <-paste('r.', names(responses)[2:2], sep = '')
remove(clia_vars, seed_vars, soil_vars, spec_vars, topo_vars, clip_vars, fire_vars, drop_cols, RAP_vars)


# Join the responses and predictors
# ID strings to separate into training and validation ID sets
valid_ids <- responses$ID %>% sample(size = round(nrow(responses)*.2))
train_ids <- responses$ID[!(responses$ID%in% valid_ids)]

# Functions for producing modeling and validation sets
train_df <- function(col){ 
  return(responses %>%
    left_join(predictors, by = 'ID') %>%
    filter(ID %in% train_ids) %>%
    dplyr::select(c(col, starts_with('p.'))) %>%
    drop_na()) }

valid_df <- function(col){ 
  return(responses %>%
    left_join(predictors, by = 'ID') %>%
    filter(ID %in% valid_ids) %>%
      dplyr::select(c(col, starts_with('p.'))) %>%
    drop_na()) }


# # Years to recovery dataframes
# RecovCo_train <- train_df('r.RecovCont')
# RecovCo_train$r.RecovCont[RecovCo_train$r.RecovCont > 100] <- 100
# RecovCo_valid <- valid_df('r.RecovCont')
# RecovCo_valid$r.RecovCont[RecovCo_valid$r.RecovCont > 100] <- 100
# 
# # Binary will-recover/won't-recover datafranes
# RecovBi_train <- train_df('r.RecovBi')
# RecovBi_valid <- valid_df('r.RecovBi')

# Contemporary condition
Post_train <- train_df('r.NDPDI_post')
Post_valid <- valid_df('r.NDPDI_post')

# # Postfire - prefire 
# Post_diff_train <- train_df('r.NDPDI_diff')
# Post_diff_valid <- valid_df('r.NDPDI_diff')

remove(responses,predictors,train_ids, valid_ids, train_df, valid_df)




# --------------------------------------------------------------------------------------
# -------------------------- Apply random forest models --------------------------------
# ---------------------- Contemporary condition ----------------------------
### Generate random number to set seed
# sample(1:1000000,1)
model_seed = 517953
set.seed(model_seed)

# # # Remove multicollinear predictors -- first remove  factors to avoid error
MCtable = Post_train[,-which(sapply(Post_train, class) == "factor")]
vars_mc <- multi.collinear(MCtable, leave.out = TRUE, p=0.05)
remove(MCtable)

# Remove collinear variables
Post_train <- Post_train %>% dplyr::select(-c(p.latitude, p.longitude, p.Fire_Name, vars_mc))
Post_valid <- Post_valid %>% dplyr::select(-c(p.latitude, p.longitude, p.Fire_Name, vars_mc))

# Perform Random Forest variable selection
Post_sel = rf.modelSel(xdata = Post_train %>% dplyr::select(starts_with('p.')), ydata = Post_train$r.NDPDI_post, r = seq(.5,.9, .1))

# Create modeling dataframe based on variable selection output
Post_train <- Post_train %>% dplyr::select(c(r.NDPDI_post, Post_sel$parameters[[4]])) %>% dplyr::select(-p.Year)

# Train Random forest model
Post_model <- randomForest(r.NDPDI_post ~ ., data = Post_train, importance = TRUE, ntree = 500)
Post_model
View(Post_model$importance)
varImpPlot(Post_model)

# Predict and plot predictions
# Predict on validation set
rf_preds <- predict(object = Post_model, newdata = Post_valid %>% dplyr::select(-c(r.NDPDI_post)))

# Bind prediction to actual data
rf_val_preds <- tibble(Predicted = rf_preds,
                       Response = Post_valid$r.NDPDI_post) %>%
  mutate(Error = abs(Predicted - Response)) %>%
  mutate(SqError = (Predicted - Response)^2)

# Calculate validation stats
postResample(rf_val_preds$Predicted, rf_val_preds$Response)
rf_rmse <- postResample(rf_val_preds$Predicted, rf_val_preds$Response)[[1]] %>% round(4)
rf_rsq <- postResample(rf_val_preds$Predicted, rf_val_preds$Response)[[2]] %>% round(4)

# Plot random forest
rf_rsq_grob <- text_grob(paste("R-squared =", as.character(rf_rsq)), x = .14, y = .95)
rf_rmse_grob <- text_grob(paste("RMSE =", as.character(rf_rmse)), x = .13, y = .9)

lm(Predicted~Response, rf_val_preds)

ggplot(rf_val_preds, aes(x=Response, y = Predicted)) +
  geom_point(alpha=.2, color = 'blue',  shape = 16) +
  geom_smooth(method='lm',formula=y~x, color = 'red') +
  geom_abline(intercept = 0, slope = 1, cex = 1) +
  labs(x = "Contemporary NDPDI from RAP", y = "Modeled NDPDI", title = "Modeling Contemporary NDPDI from RAP") +
  annotation_custom(rf_rsq_grob) +
  annotation_custom(rf_rmse_grob) +
  scale_x_continuous(limits = c(0,2), breaks = seq(0,2,.4)) +
  scale_y_continuous(limits = c(0,2), breaks = seq(0,2,.4)) +
  theme_minimal()

ggplot(Post_train)+
  geom_point(mapping = aes(x = p.NDPDI_pre, y = r.NDPDI_post, color = p.Fire_Name))+
  labs(x = 'Pre-fire species richness', y = 'Pre-fire NDPDI - post-fire NDPDI')+
  theme_minimal()
#theme(legend.position = "none")



#------------------------ Apply predictive model to Holloway Fire --------------------------------
# Apply predictive model for 15 year post-fire condition
tifs <- list.files('C:/Users/erjensen/Documents/Thesis/ThesisAIM/data/Recov_rasters/', pattern='tif', full.names = TRUE)
ModelRasters = stack(tifs)
names(ModelRasters) <- paste('p.', names(ModelRasters), sep = '')
LongDraw = predict(ModelRasters, Post_model)
plot(LongDraw)

# Assess post-fire LongDraw visually
LongDrawPost <- raster('C:/Users/erjensen/Desktop/NDPDI_post.tif')
plot(LongDrawPost)

writeRaster(LongDraw, 'C:/Users/erjensen/Desktop/NDPDI_Predpost.tif')
