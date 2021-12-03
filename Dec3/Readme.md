### getting started .. ##

Run below code in your 'terminal'
```
cd gpfs/scratch/[your PSU ID]/
mkdir -p DAWG2021F-3rd
cd DAWG2021F-2rd

wget https://github.com/DAWG-psu/DAWG2021F/raw/main/Dec3/DAWG_2021F_3rd.R
wget https://github.com/DAWG-psu/DAWG2021F/raw/main/Dec3/DAWG_2021F_3rd.RData
```


# DAWG 2021F 3rd workshop 
## Machine learning application on microbiome dataset

## Table of contents
1. What is machine learning
2. Why it could be useful for microbiome data analysis
3. Random forest
4. Machine learning procedure
5. R Script


## What is machine learning?

Let machine learn from the data, and predict or decide what will happen next 

<img width="724" alt="Screen Shot 2021-11-30 at 10 49 53 PM" src="https://user-images.githubusercontent.com/77017866/144168784-d4d11e06-02f5-4e56-a20f-8be1ba4cdb25.png">


Machine learning algorithms build a model based on sample data (trainind data), in order to make predictions or decisions without being explicitly programmed to do so. 

<img width="545" alt="image" src="https://user-images.githubusercontent.com/77017866/144168951-0d9a64a5-b579-41f1-91cc-97e202a8e87c.png">


From the above figure, you can see that the subset of data can be used as "Training Dataset" which can be used for building a model. Then, "Testing Dataset"(Or Validation dataset) can be used to evaluate the model. Model can be used to predict the new datapoint. 

Doesn't it sound familiar?

![image](https://user-images.githubusercontent.com/77017866/144310725-c2d9388a-27d4-4f29-8bd6-20c0d88980b8.png)

Yes, linear regression! Using data point to build a model (equation), and new data point can be used for prediction. Linear regression is one of the most well known algorithms in machine learning for less complicated data.


## Why it could be useful for microbiome data analysis

![image](https://user-images.githubusercontent.com/77017866/144312193-9f9dc239-688e-4a35-9d2b-2f90f3086158.png)
Cammarota et al., 2020

In most cases of the microbiome analysis, we have i) taxonomic count table (either ASVs or OTUs), ii) taxonomic informatio (taxonomy table), and iii) metadata (or sample data). It is obviously true that microbiome data is somewhat too complicated and too many variables (ASVs or OTUs count) in them. Thus, it is NOT simple to find out the significant relationship bewteen each taxon to the sample characteristics (from sample data). Thus, we hope that machine learning can explain the relationship between sample characteristics with microbiome composition.

Please note that, machine learning is not for getting a perfect answer (always true), but most plausible answer (most likely true). In most cases of the research, machine learning is used for 'hypothesis generating' which requires further confirmation in either *in vivo* or *in vitro*.

## Random Forest

Random forest is the one of the most common machine learning algorithms which many micorbiome researchers alwasys used. When I first use RF for my research, there were only a few paper were already published about microbiome data analysis using random forest, however, now it is very common to find papers applying RF.

Random forest is based on decision tree.

<img width="830" alt="Screen Shot 2021-12-01 at 3 23 44 PM" src="https://user-images.githubusercontent.com/77017866/144308363-fd509df7-1bfe-45b8-a3ea-2049097c4dbb.png">


<img width="838" alt="Screen Shot 2021-12-01 at 3 24 35 PM" src="https://user-images.githubusercontent.com/77017866/144308412-1ba28ecd-7d88-4e56-9807-00c7dc3c1415.png">

## Machine learning procedure

### Tuning parameters
#### empirical tuning (grid search vs random search)

![image](https://user-images.githubusercontent.com/77017866/144506862-f60da1c2-d7e0-42fa-9f2e-e76bce660270.png)

![image](https://user-images.githubusercontent.com/77017866/144521135-1e62b6e0-c458-4786-970d-9e725550b4d9.png)


**mtry**: Number of variables randomly sampled as candidates at each split

**ntree**: Number of trees to grow

**nodesize**: Minimum number of sample to stop splitting

### Model evaluation

#### Area under the curve
![image](https://user-images.githubusercontent.com/77017866/144507018-0e66dad9-2467-48be-ac41-44bd80816362.png)

ROC curve is a performance measurement for the classification problem based on the model sensitivity and specificity. It is a probability curve and its area under the curve (AUC) represents the degree of separability. The ROC curve is plotted with TPR against the FPR where TPR is on the y-axis and FPR is on the x-axis.

#### Kappa statistics

Kappa is a chance-corrected measure of agreement (reliability of the classification). A mean for evaluating the prediction performance of classifiers

![image](https://user-images.githubusercontent.com/77017866/144539573-19c0042b-4217-4864-bfcd-e45e7b819ecc.png)

### Cross validation VS hold-out validation


#### Cross validation
![image](https://user-images.githubusercontent.com/77017866/144508953-4bd49313-fe0a-440e-a789-82e38fadd910.png)

Cross-validation or ‘k-fold cross-validation’ is when the dataset is randomly split up into ‘k’ groups. One of the groups is used as the test set and the rest are used as the training set. The model is trained on the training set and scored on the test set. Then the process is repeated until each unique group as been used as the test set.

For example, for 4-fold cross validation, the dataset would be split into 5 groups, and the model would be trained and tested 5 separate times so each group would get a chance to be the test set. This can be seen in the graph below.

#### Hold-out validation

Hold-out is when you split up your dataset into a ‘train’ and ‘test’ set. The training set is what the model is trained on, and the test set is used to see how well that model performs on unseen data. A common split when using the hold-out method is using 80% of data for training and the remaining 20% of the data for testing.

#### Hold-out VS Cross-validation

Hold-out has its advantage of using completely independent dataset to validate the model. Also, it is computationally less intensive because it did not require multiple validations. However, Hold-out might subject to higher variance given the samller size of the data. Also, it is important to split training and testing data in similar ratio (but completely random)

Cross validation is using entire dataset to evaluate the model. Also, it requires less number of sample than hold-out method. However, it is computationally intense and still have chance to have over-fitting problem.

### Variable importance measures

Varialbe importance is indicating how much a given model uses the variable to make accurate predictions. In other words, how important the variable is for classification. Random forest R package has two default variable importance measures (i.e., Mean decrease accuracy, Mean decrease Gini), Mean decrease AUC is also widely used to calculate the variable importance. 

#### Mean decrease accuracy

The Mean Decrease Accuracy  expresses how much accuracy the model losses by excluding each variable. The more the accuracy suffers, the more important the variable is for the successful classification.

#### Mean decrease Gini index

The Mean Decrease in Gini coefficient is a measure of how each variable contributes to the homogeneity of the nodes and leaves in the resulting random forest. 

![image](https://user-images.githubusercontent.com/77017866/144508150-f376f9a7-e324-442e-8104-9d52e435265f.png)

#### Mean decrease AUC

The Mean Decrease AUC is a measure of how much AUC the model losses by excluding each variable. 

## Rscript and annotation

### Processing data (phyloseq object to table)

1. Convert to desired taxonomy level (e.g, genus or family)
```
phyloseq_genus <- taxa_level(phyloseq, "Genus")
```
2. Generating relative abundance table of ASV table
```
phyloseq_relative_abundance <- transform_sample_counts(phyloseq_genus, function(x) x/sum(x))
taxa_are_rows(phyloseq_relative_abundance) ##Taxa should be in column
phyloseq_relative_abundance <- load("phyloseq_object.rds")
asv_relative_abundance <- otu_table(phyloseq_relative_abundance)
metadata <- sample_data(phyloseq_relative_abundance) ## metadata
metadata$salmonella 
```

3. Creating one table containing all the information we need for the classification
```
model_salmonella <- as.data.frame(cbind(asv_relative_abundance, "salmonella" = metadata$salmonella))
model_salmonella$salmonella <- as.factor(model_salmonella$salmonella) # Change the feature as 'factor'
```

4. Filter out unnecessary characters in taxon name
```
names(model_salmonella) = gsub(pattern = " ", replace= "_", x = names(model_salmonella))
names(model_salmonella)[names(model_salmonella) == 'f__'] <- 'f__unclassified'
names(model_salmonella) = gsub(pattern = "f__", replace= "", x = names(model_salmonella))
names(model_salmonella)[names(model_salmonella) == 'g__'] <- 'g__unclassified'
names(model_salmonella) = gsub(pattern = "g__", replace= "", x = names(model_salmonella))
names(model_salmonella) = gsub(pattern = '\\[', replace= "", x = names(model_salmonella))
names(model_salmonella) = gsub(pattern = '\\]', replace= "", x = names(model_salmonella))
names(model_salmonella) = gsub(pattern = '\\/' , replace= "OR", x = names(model_salmonella))
names(model_salmonella) = gsub(pattern = '-' , replace = "_", x=names(model_salmonella))
names(model_salmonella) = gsub(pattern = "\\(", replace = "_", x = names(model_salmonella))
names(model_salmonella) = gsub(pattern = "\\)", replace = "_", x = names(model_salmonella))
model_salmonella[sapply(model_salmonella, is.character)] <- lapply(model_salmonella[sapply(model_salmonella, is.character)], 
                                                             as.factor)

```

### Tuning parameters of random forest and making model using MLR package in R 

Parameters will be tuned in R based on the models' accuracy (AUC) and reliability (Kappa)

```
require(mlr)
require(party)
require(randomForest)
require(ggplot2)

```

1. Making "task" (makeClassifTask)
```
trainTask<- makeClassifTask(data = model_salmonella, # dataset
                            target = "salmonella", # name of the column in the data rame with Salmonella presence/absence data in it
                            positive = "1") # how positive Salmonella samples are marked in that column)
```

2. Setting "parameters" (makeParamSet)

* You can find a full lists of the learning methods and their parameters integrated in mlr here : https://mlr.mlr-org.com/articles/tutorial/integrated_learners.html

```
parameters <- makeParamSet(
    makeDiscreteParam("ntree",values=c(501)), # always use an odd number of trees
    makeDiscreteParam("mtry", values=seq(5,ncol(otu_table(phyloseq_relative_abundance)), 50)), # this says test every 5th number - adapt range based on number of features in dataset
    makeDiscreteParam("nodesize", values=seq(5,40,10)) # sets max possible size of each tree - adapt range based on number of samples in dataset
    )
print(parameters)
```

3. Setting "evaluation and optimization algorithm
```
set_cv <- makeResampleDesc("RepCV", reps = 3, folds=10) # 10-fold cross validation repeated 3 times
control <- makeTuneControlGrid() # grid search of the optimal parameter sets
```

4. Make learner
```
learner <- makeLearner("classif.randomForest", predict.type = 'prob') # set to prob so the outcome is the predicted probability of Salmonella being present
```

5. Start the tuning process

```
tune.RF <- tuneParams(learner = learner, 
                      resampling = set_cv, 
                      task = trainTask, 
                      par.set = parameters, 
                      control = control, 
                      measures = list(auc, kappa)) # tune based on auc and kappa
tune.RF
```

6. Set final parameters for actual model
```
parameter_final <- setHyperPars(learner, par.vals = tune.RF$x)
```

7. Making randomforest model
```
trained.model.RF <- mlr::train(parameter_final, trainTask) # train the final model, see https://mlr.mlr-org.com/articles/tutorial/train.html
learner.RF<-(getLearnerModel(trained.model.RF)) # get the trained forest in a format you can use with other packages
learner.RF
```

* Alternatively, you can use RandomForest package directly after tuning - more flexibility
```
library(randomForest)
randomforest_model <- randomForest(salmonella~., data = model_salmonella, 
                                   mtry = tune.RF$x$mtry, ntree = tune.RF$x$ntree, nodesize = tune.RF$x$nodesize, 
                                   importance=TRUE)

```

### Variable importance measures

1. Extract mean decrease accuracy and mean decrease gini

```
varImp_RF <- as.data.frame(importance(randomforest_model, scale=T))
varImpPlot(randomforest_model)
varImp_gini <- varImp_RF[order(varImp_RF$MeanDecreaseGini, decreasing = T),]
varImp_acc <- varImp_RF[order(varImp_RF$MeanDecreaseAccuracy, decreasing = T),]
```

2. Create bar plot with top 15 most informative genus
```
varimp_plotting <- function(varImp_object){
  varImp_object <- varImp_object[1:15,]
  varImp_object[,5] <- varImp_object[,4]/sum(varImp_object[,4])
  p <- ggplot(varImp_object, aes(x = reorder(row.names(varImp_object),V5), y = V5)) + 
    geom_bar(stat = "identity",color='black', fill='gray') + 
    coord_flip()+ 
    xlab("") +
    ylab("") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text=element_text(size=10), axis.title=element_text(size=10),
          text=element_text(size=10, color="black"))
  return(p)
}

varimp_plotting(varImp_gini)
varimp_plotting(varImp_acc)
```

![image](https://user-images.githubusercontent.com/77017866/144528820-93f99f75-b1b2-4e25-902e-3bd6d4dbf071.png)

This looks good. We identified a few genus that is associated with presence of Salmonella in the sample. However, it is not the end. We have to check a few things before we decided to make a conclusion

### Check the results

1. Correlation between genus

```
varimp_15 <- varImp_gini[1:15,]
genus_rf <- model_salmonella[, rownames(varimp_15)]
corr <- cor(genus_rf)
library(corrplot)
corrplot(corr, method="color", order = "AOE", type="lower", tl.col = "black", cl.ratio =0.2, tl.srt =45, tl.cex =0.7)
```
![image](https://user-images.githubusercontent.com/77017866/144529368-e135b605-3473-44f7-bc12-6ca7ad160b80.png)

Based on the plot, we can find a huge correlation between a number of genus in the sample, this type of data structure needs to be either addressed or mitigated. One possible solution is to choose different machine learning algorithm, which can take care of the collinearity problem (i.e., conditional random forest)

2. taxa prevalence

```
boxplot(genus_rf, col=rainbow(20), las=2)
```
![image](https://user-images.githubusercontent.com/77017866/144530920-36a03ad6-22ca-4f19-ae0c-a7290d7e9436.png)




