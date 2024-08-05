library(tidyverse)
library(randomForest)
library(rfUtilities)
library(caret)
library(ggpubr)

# Function to clean dataframes
drop_cols = function(df){
  return (df %>% 
            dplyr::select(-c(`.geo`, Year, Fire_Name, `system:index`)) %>%
            mutate(ID = `paste(___)`) %>%
            select(-`paste(___)`))}


# ---------------------- Read in Predictors CSVs ---------------------------------
# Topographic
setwd('C:/Users/erjensen/Documents/Thesis/ThesisAIM/data/Recov_CSVs/')
topo_vars2 <- read_csv('Biophysical_Predictors_Mode.csv') %>%
  mutate(landform = as.factor(mode)) %>%
  drop_cols() %>%
  select(-mode)

topo_vars <- read_csv('Biophysical_Predictors_Mean.csv') %>%
  mutate(ID = `paste(___)`) %>%
  select(-c(`.geo`, `system:index`,`paste(___)`)) %>% 
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
  select(-mean)

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
  select(-`paste.___.`) %>%
  select(c(ID, NDPDI_pre))

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
  select(ID,everything()) %>%
  drop_na()
names(predictors)[2:ncol(predictors)] <-paste('p.', names(predictors)[2:ncol(predictors)], sep = '')

# Read in response variables
responses = read.csv('Recov_Response.csv') %>%
  mutate(ID = `paste.___.`) %>%
  select(c(ID, NDPDI_diff, NDPDI_post, RecovBi, RecovCont)) %>%
  mutate(RecovBi = as.factor(RecovBi))
names(responses)[2:5] <-paste('r.', names(responses)[2:5], sep = '')
remove(clia_vars, seed_vars, soil_vars, spec_vars, topo_vars, clip_vars, fire_vars, drop_cols, RAP_vars)


# -------------- Functions for producing modeling and validation sets -----------------

train_df <- function(col){ 
  return(responses %>%
    left_join(predictors, by = 'ID') %>%
    select(c(col, starts_with('p.'))) %>%
    select(-c(p.longitude, p.latitude)) %>%
    dplyr::rename(response = col) %>%
    drop_na()) }

# Years to recovery dataframes
RecovCo_df <- train_df('r.RecovCont')
RecovCo_df$r.RecovCont[RecovCo_df$r.RecovCont > 100] <- 100

# Binary will-recover/won't-recover datafranes
RecovBi_df <- train_df('r.RecovBi') 

# Contemporary condition
Post_df <- train_df('r.NDPDI_post')

# Postfire - prefire 
Post_diff_df <- train_df('r.NDPDI_diff')

remove(responses,predictors, train_df)


# ------------------ Modeling post-fire condition / recovery -------------------

# Function for leave-one-out variable selection and validation for random forest for regression
LeaveOneOut_reg <- function(df, varSelThres){
  
  # Generate random number to set seed
  model_seed = 517953
  set.seed(model_seed)
  
  # Remove multicollinear predictors -- first remove  factors to avoid error
  MCtable = df[,-which(sapply(df, class) == "factor")]
  vars_mc <- multi.collinear(MCtable, leave.out = TRUE, p=0.05)
  df <- df %>% select(-vars_mc)
  
  # Perform Random Forest variable selection
  varSel = rf.modelSel(xdata = df %>% dplyr::select(starts_with('p.')), ydata = df$response, r = varSelThres)
  
  # Keep select columns for modeling from variable selection
  df_sel <- df %>% select(c(response, varSel$parameters[[2]])) 
  
  # Create tibble to store validation values
  val_tb <- tibble()
  
  #List of fires to eventually map nested function over
  fire_list = as.character(df$p.Fire_Name %>%unique())
  
  # Function to train and apply random forest model using leave-one-out validation
  apply_model <- function(df, fire){
    
    # Filter into training and validation sets
    df_train <- filter(df, p.Fire_Name != fire) %>%
      select(-p.Fire_Name)
    df_valid <- filter(df, p.Fire_Name == fire) %>%
      select(-p.Fire_Name)
    
    # Model the response variable
    # Run tuning function
    set.seed(model_seed)
    tuned <-tuneRF(y = df_train$response, x = df_train %>% dplyr::select(starts_with('p.')), ntreeTry=200)
    
    # Get the mtry value that minimizes error
    mtry_opt <- tuned[, 'mtry'][which.min(tuned[,'OOBError'])]
    
    # Train the model
    model = randomForest(response ~ ., data = df_train, importance = TRUE, ntree = 500, mtry = mtry_opt)
    
    # Predict on validation set
    rf_preds <- predict(object = model, newdata = df_valid %>% dplyr::select(starts_with('p.')))
    
    # Bind prediction to actual data
    rf_val_preds <- tibble(Predicted = rf_preds,
                           Response = df_valid$response,
                           Fire = fire) %>%
      mutate(Error = abs(Predicted - Response)) %>%
      mutate(SqError = (Predicted - Response)^2)
    
    # Calculate validation stats
    postResample(rf_val_preds$Predicted, rf_val_preds$Response)
    rf_rmse <- postResample(rf_val_preds$Predicted, rf_val_preds$Response)[[1]] %>% round(4)
    rf_rsq <- postResample(rf_val_preds$Predicted, rf_val_preds$Response)[[2]] %>% round(4)
    
    val_one <- tibble(Fire = fire,
                      rsq = rf_rsq,
                      rmse = rf_rmse)
    
    val_tb <- bind_rows(val_tb, val_one)
    
    return(val_tb) }
  
  # Map function over list of fires
  val_full_df <- map(fire_list, apply_model, df = df_sel)
  return(val_full_df)
}

# Contemporary condition
Post_valid = LeaveOneOut_reg(Post_df, .8) %>% bind_rows()
varImpPlot(Post_valid[[2]])

Post_valid_indv <- Post_valid[[2]]
Post_valid_indv$importance 

mean(Post_valid$rmse)
mean(Post_valid$rsq)


# -------------- Works after returning rf_val_preds -----------------------

# Plot the outputs
plot<-ggplot(Post_valid_8)+
  geom_point(mapping = aes(x = Response, y = Predicted, color = as.factor(Fire)), alpha = .3)+
  lims(x = c(0,2), y = c(0,2))+
  facet_wrap(~Fire)+
  geom_abline(slope = 1, intercept = 0)+
  theme_minimal() +
  theme(legend.position = 'none') 

png('C:/Users/erjensen/Desktop/Plots_ContempModel.png', type = "cairo", units ="in", width = 8, height = 8, pointsize = 16, res = 400)
plot
dev.off()


# ---------------------------------------------------------------------

# Difference between post-fire and pre-fire
Post_diff_valid_8 = LeaveOneOut_reg(Post_diff_df, .8) %>% bind_rows()


# ------------------ Modeling binary classification -------------------

# Function for leave-one-out variable selection and validation for random forest for binary classification
LeaveOneOut_bin <- function(df, varSelThres){
  
  ### Generate random number to set seed
  # sample(1:1000000,1)
  model_seed = 517953
  set.seed(model_seed)
  
  
  # Perform Random Forest variable selection
  varSel = rf.modelSel(xdata = df %>% dplyr::select(starts_with('p.')), ydata = df$response, r = varSelThres)
  
  # Keep select columns for modeling from variable selection
  df_sel <- df %>% select(c(response, p.Fire_Name, varSel$parameters[[2]])) 
  
  # Create tibble to store validation values
  val_tb <- tibble()
  
  #List of fires to eventually map nested function over
  fire_list = as.character(df$p.Fire_Name %>%unique())
  
  # Function to train and apply random forest model using leave-one-out validation
  apply_model <- function(df, fire){
    
    # Filter into training and validation sets
    df_train <- dplyr::filter(df, p.Fire_Name != fire) %>%
      select(-p.Fire_Name)
    df_valid <- dplyr::filter(df, p.Fire_Name == fire) %>%
      select(-p.Fire_Name)
    
    # Model the response variable
    # Run tuning function
    set.seed(model_seed)
    tuned <-tuneRF(y = df_train$response, x = df_train %>% dplyr::select(starts_with('p.')), ntreeTry=200)
    
    # Get the mtry value that minimizes error
    mtry_opt <- tuned[, 'mtry'][which.min(tuned[,'OOBError'])]
    
    # Train the model
    model = randomForest(response ~ ., data = df_train, importance = TRUE, ntree = 500, mtry = mtry_opt)
    
    # Predict on validation set
    rf_preds <- predict(object = model, newdata = df_valid %>% dplyr::select(starts_with('p.')))
    
    # Bind prediction to actual data
    rf_val_preds <- tibble(Predicted = rf_preds,
                           Response = df_valid$response,
                           Fire = fire) %>%
      mutate(Error = abs(Predicted - Response)) %>%
      mutate(SqError = (Predicted - Response)^2)
    
    # Calculate validation stats
    postResample(rf_val_preds$Predicted, rf_val_preds$Response)
    rf_rmse <- postResample(rf_val_preds$Predicted, rf_val_preds$Response)[[1]] %>% round(4)
    rf_rsq <- postResample(rf_val_preds$Predicted, rf_val_preds$Response)[[2]] %>% round(4)
    
    val_one <- tibble(Fire = fire,
                      rsq = rf_rsq,
                      rmse = rf_rmse)
    
    val_tb <- bind_rows(val_tb, val_one)
    
    return(rf_val_preds)}
  
  # Map function over list of fires
  val_full_df <- map(fire_list, apply_model, df = df_sel)
  return(val_full_df)
}

RecovBi_8 <- LeaveOneOut_bin(RecovBi_df, .8)
RecovBi_8_df <- RecovBi_8 %>% bind_rows()

sensitivity(RecovBi_8_df$Predicted, RecovBi_8_df$Response)
specificity(RecovBi_8_df$Predicted, RecovBi_8_df$Response)


# --------------------------------------------------------------------------------------
# -------------------------- Apply random forest models --------------------------------

# ---------------------------- Years to Recovery ------------------------------
### Generate random number to set seed
# sample(1:1000000,1)
model_seed = 517953
set.seed(model_seed)

# # Remove multicollinear predictors -- first remove  factors to avoid error
MCtable = RecovCo_train[,-which(sapply(RecovCo_train, class) == "factor")]
vars_mc <- multi.collinear(MCtable, leave.out = TRUE, p=0.05)
remove(MCtable)

# Remove collinear variables
RecovCo_train <- RecovCo_train %>% dplyr::select(-c(vars_mc, p.latitude, p.longitude))
RecovCo_valid <- RecovCo_valid %>% dplyr::select(-c(vars_mc, p.latitude, p.longitude))

# Create modeling dataframe based on variable selection output
RecovCo_train <- RecovCo_train %>% select(c(r.RecovCont, RecovCo_sel$parameters[[4]]))

# Train Random forest model
RecovCo_model <- randomForest(r.RecovCont ~ ., data = RecovCo_train, importance = TRUE, ntree = 500)
RecovCo_model
plot(RecovCo_model)

# Predict and plot predictions
# Predict on validation set
rf_preds <- predict(object = RecovCo_model, newdata = RecovCo_valid %>% dplyr::select(-c(r.RecovCont)))

# Bind prediction to actual data
rf_val_preds <- tibble(Predicted = rf_preds,
                       Response = RecovCo_valid$r.RecovCont) %>%
  mutate(Error = abs(Predicted - Response)) %>%
  mutate(SqError = (Predicted - Response)^2)

# Calculate validation stats
postResample(rf_val_preds$Predicted, rf_val_preds$Response)
rf_rmse <- postResample(rf_val_preds$Predicted, rf_val_preds$Response)[[1]] %>% round(4)
rf_rsq <- postResample(rf_val_preds$Predicted, rf_val_preds$Response)[[2]] %>% round(4)

# Plot random forest
rf_rsq_grob <- text_grob(paste("R-squared =", as.character(rf_rsq)), x = .14, y = .95)
rf_rmse_grob <- text_grob(paste("RMSE =", as.character(rf_rmse)), x = .13, y = .9)

# Generate linear regression
lm(Predicted~Response, rf_val_preds)

# Visualize the results
ggplot(rf_val_preds, aes(x=Response, y = Predicted)) +
  geom_point(alpha=.2, color = 'blue',  shape = 16) +
  geom_smooth(method='lm',formula=y~x, color = 'red') +
  geom_abline(intercept = 0, slope = 1, cex = 1) +
  labs(x = "RAP Dataset Years to Recovery", y = "Modeled Years to Recovery", title = "Modeling Years to Recovery of Pre-fire NDPDI") +
  annotation_custom(rf_rsq_grob) +
  annotation_custom(rf_rmse_grob) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20)) +
  theme_minimal()


# ---------------------------- Binary recovery ------------------------------

### Generate random number to set seed
set.seed(model_seed)

# Remove collinear variables
RecovBi_train <- RecovBi_train %>% dplyr::select(-c(vars_mc, p.latitude, p.longitude))
RecovBi_valid <- RecovBi_valid %>% dplyr::select(-c(vars_mc, p.latitude, p.longitude))

# Create modeling dataframe based on variable selection output
RecovBi_train <- RecovBi_train %>% select(c(r.RecovBi, RecovBi_sel$parameters[[4]]))

# Train Random forest model
RecovBi_model <- randomForest(r.RecovBi ~ ., data = RecovBi_train, importance = TRUE, ntree = 500)
RecovBi_model

# Predict and plot predictions
# Predict on validation set
rf_preds <- predict(object = RecovBi_model, newdata = RecovBi_valid %>% dplyr::select(-c(r.RecovBi)))

# Bind prediction to actual data
rf_val_preds <- tibble(Predicted = rf_preds,
                       Response = RecovBi_valid$r.RecovBi) 
sensitivity(rf_val_preds$Predicted, rf_val_preds$Response)
specificity(rf_val_preds$Predicted, rf_val_preds$Response)


# ---------------------- Contemporary condition ----------------------------
### Generate random number to set seed
set.seed(model_seed)

# Remove collinear variables
Post_train <- Post_train %>% dplyr::select(-c(p.latitude, p.longitude, p.Fire_Name))
Post_valid <- Post_valid %>% dplyr::select(-c(p.latitude, p.longitude, p.Fire_Name))

# Perform Random Forest variable selection
Post_sel = rf.modelSel(xdata = Post_train %>% dplyr::select(starts_with('p.')), ydata = Post_train$r.NDPDI_post, r = seq(.5,.9, .1))

#Assess the output from variable selection
plot(Post_sel)
Post_sel

# Create modeling dataframe based on variable selection output
Post_train <- Post_train %>% select(c(r.NDPDI_post, Post_sel$parameters[[4]]))

# Train Random forest model
Post_model <- randomForest(r.NDPDI_post ~ ., data = Post_train, importance = TRUE, ntree = 500)
Post_model
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

# ---------------------- Postfire-prefire differenced ----------------------------
### Generate random number to set seed
set.seed(model_seed)

# Remove collinear variables
Post_diff_train <- Post_diff_train %>% dplyr::select(-c(vars_mc, p.latitude, p.longitude, p.Fire_Name))
Post_diff_valid <- Post_diff_valid %>% dplyr::select(-c(vars_mc, p.latitude, p.longitude, p.Fire_Name))

# Perform Random Forest variable selection
# Post_diff_sel = rf.modelSel(xdata = Post_diff_train %>% dplyr::select(starts_with('p.')), ydata = Post_diff_train$r.NDPDI_diff, r = seq(.5,.9, .1))
#Assess the output from variable selection
# plot(Post_diff_sel)
# Post_diff_sel


# ------------------------ RUN THIS -------------------------------------

# Create modeling dataframe based on variable selection output
Post_diff_train <- Post_diff_train %>% select(c(r.NDPDI_diff, Post_diff_sel$parameters[[4]]))

# Train Random forest model
Post_diff_model <- randomForest(r.NDPDI_diff ~ ., data = Post_diff_train, importance = TRUE, ntree = 500)
Post_diff_model
varImpPlot(Post_diff_model)

# Predict and plot predictions
# Predict on validation set
rf_preds <- predict(object = Post_diff_model, newdata = Post_diff_valid %>% dplyr::select(-c(r.NDPDI_diff)))

# Bind prediction to actual data
rf_val_preds <- tibble(Predicted = rf_preds,
                       Response = Post_diff_valid$r.NDPDI_diff) %>%
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
  scale_x_continuous(limits = c(-2,2), breaks = seq(-2,2,.4)) +
  scale_y_continuous(limits = c(-2,2), breaks = seq(-2,2,.4)) +
  theme_minimal()


