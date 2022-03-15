source("heart_explore.R")

#Load libraries
library(kknn)

'
To do:
1.remove known predictors that are questionable
2. change to num_b
3. plot confusion matrix and comment on reuslts compared to literature. 
4.Determine the K value to report in manuscript
'
#Remove ambiguous predictors 

#change num to a factor 
heart_train <- heart_train %>% 
  mutate(num = ifelse(num > 0,"disease","healthy"))

heart_train <- heart_train %>% 
  select(-c(thalach, sex))


#Predictor columns
names <- colnames(heart_train %>% select(-num))

#Tibble to store accuracies:
accuracies <- tibble(size = integer(),
                     model_string = character(),
                     accuracy = numeric())

#create a model spec
knn_spec <- nearest_neighbor(weight_func = "rectangular",
                             neighbors = tune()) %>%
  set_engine("kknn") %>% 
  set_mode("classification")

#Create a 5-fold cross-validation object
heart_vfold <- vfold_cv(heart_train,
                        v = 5,
                        strata = num)

#Store the total number of predictors 
n_total <- length(names)

#Store selected predictors
selected <- c()

for (i in 1:n_total) {
  accs <- list()
  models <- list()
  for (j in 1:length(names)){
    preds_new <- c(selected, names[[j]])
    model_string <- paste("num", "~", paste(preds_new, collapse = '+'))
    
    #Create recipe from model string
    heart_recipe <- recipe(as.formula(model_string),
                            data = heart_train) %>% 
      step_scale(all_predictors()) %>% 
      step_center(all_predictors())
    
    #Tune KNN classifier with these predictors
    #collect accuracy for best k
    
    acc <- workflow() %>% 
      add_recipe(heart_recipe) %>% 
      add_model(knn_spec) %>% 
      tune_grid(resamples = heart_vfold, grid = 10) %>% 
      collect_metrics() %>% 
      filter(.metric == "accuracy") %>% 
      summarise(mx = max(mean))
    acc <- acc$mx %>% unlist()
    
    #add result to dataframe
    accs[[j]] <- acc
    models[[j]] <- model_string
  }
  jstar <- which.max(unlist(accs))
  accuracies <- accuracies |> 
    add_row(size = i, 
            model_string = models[[jstar]], 
            accuracy = accs[[jstar]])
  selected <- c(selected, names[[jstar]])
  names <- names[-jstar]
}

write.table(accuracies, file = "heart_forward2.csv", row.names = FALSE)

# heart_train_explore %>% 
#   ggpairs()
#Notes from TA: 
"
Binarize/factor the nums to disease & healthy,look at stacked histograms and 
compare the distribution.
For a given value if the means between the correlation histograms are different
enough then that predictor can be used to distinguish b/w disease states, that 
predictor can be used. If the means/distributions are too similar/overlap,
then the algorithim will have diffculty classifying a value. 

"
