##making a function called 'taxa_level'. Very similar as tax_glom from phyloseq package
taxa_level <- function(physeq,which_level){
  #enforce orientation
  OTU <- otu_table(physeq)
  SAM <- sample_data(physeq)
  OTU_taxonomy <- tax_table(physeq)
  new_abund_table<-NULL
  if(which_level=="Otus"){
    OTU_tree <- phy_tree(physeq)
    new_abund_table<-OTU
  } else {
    list<-na.omit(unique(OTU_taxonomy[,which_level]))
    new_abund_table<-NULL
    for(i in list){
      rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i])
      tmp<-data.frame(rowSums(OTU[,rt]))
      if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
      if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
    }
  }
  OTU<-as.data.frame(as(new_abund_table,"matrix"))
  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(SAM)
  #reconstruct the phyloseq object
  physeq<-NULL
  if(which_level=="Otus"){
    physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
  } else {
    physeq<-merge_phyloseq(phyloseq(OTU),SAM)
  }
  return(physeq)
}

#Convert to desired taxonomy level (i.e., genus)
phyloseq_genus <- taxa_level(phyloseq, "Genus")

#Generating relative abundance table
phyloseq_relative_abundance <- transform_sample_counts(phyloseq_genus, function(x) x/sum(x))
taxa_are_rows(phyloseq_relative_abundance) ##Taxa should be in column
#phyloseq_relative_abundance <- load("phyloseq_object.rds")
asv_relative_abundance <- otu_table(phyloseq_relative_abundance)
metadata <- sample_data(phyloseq_relative_abundance) ## metadata
metadata$salmonella 
#Creating one table containing all the information we need for the classification
model_salmonella <- as.data.frame(cbind(asv_relative_abundance, "salmonella" = metadata$salmonella))
model_salmonella$salmonella <- as.factor(model_salmonella$salmonella) # Change the feature as 'factor'

#Filter out unnecessary characters in taxon name
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

#Tuning hyperparameter using mlr package
install.packages("mlr")
install.packages("party")
install.packages("randomForest")
install.packages("ggplot2")

require(mlr)
require(party)
require(randomForest)
require(ggplot2)

#task - providing dataset
trainTask<- makeClassifTask(data = model_salmonella, # dataset
                            target = "salmonella", # name of the column in the data rame with Salmonella presence/absence data in it
                            positive = "1") # how positive Salmonella samples are marked in that column)

#define parameters
#parameters <- makeParamSet(
#  makeDiscreteParam("ntree",values=c(501)), # always use an odd number of trees
#  makeDiscreteParam("mtry", values=seq(5,ncol(otu_table(phyloseq_genus)), 50)), # this says test every 5th number - adapt range based on number of features in dataset
#  makeDiscreteParam("nodesize", values=seq(5,40,10)) # sets max possible size of each tree - adapt range based on number of samples in dataset
#)
#print(parameters)

parameters <- makeParamSet(
  makeDiscreteParam("ntree",values=c(501)), # always use an odd number of trees
  makeDiscreteParam("mtry", values=seq(5,ncol(otu_table(phyloseq_relative_abundance)), 200)), # this says test every 5th number - adapt range based on number of features in dataset
  makeDiscreteParam("nodesize", values=c(5)) # sets max possible size of each tree - adapt range based on number of samples in dataset
)
print(parameters)

#define evaluation and optimization algoritm
set_cv <- makeResampleDesc("RepCV", reps = 3, folds=10) # 10-fold cross validation repeated 3 times
control <- makeTuneControlGrid() # grid search of the optimal parameter sets

# Make learner
learner <- makeLearner("classif.randomForest", predict.type = 'prob')

#Start the tuning process (it might take a while)
tune.RF <- tuneParams(learner = learner, 
                      resampling = set_cv, 
                      task = trainTask, 
                      par.set = parameters, 
                      control = control, 
                      measures = list(auc, kappa)) # tune based on auc and kappa
tune.RF
#After tuning, set final parameters for actual model development
parameter_final <- setHyperPars(learner, par.vals = tune.RF$x)

#Making RF model and extract the model object
trained.model.RF <- mlr::train(parameter_final, trainTask) # train the final model, see https://mlr.mlr-org.com/articles/tutorial/train.html
learner.RF<-(getLearnerModel(trained.model.RF)) # get the trained forest in a format you can use with other packages
learner.RF
## Variable importance ##
varImp <- as.data.frame(importance(learner.RF, scale=T, type = 2))
varImpPlot(learner.RF)

## OR we can use random forest package directly after tuning - more flexibility
library(randomForest)
randomforest_model <- randomForest(salmonella~., data = model_salmonella, 
                                   mtry = tune.RF$x$mtry, ntree = tune.RF$x$ntree, nodesize = tune.RF$x$nodesize, 
                                   importance=TRUE)
varImp_RF <- as.data.frame(importance(randomforest_model, scale=T))
varImpPlot(randomforest_model)
varImp_gini <- varImp_RF[order(varImp_RF$MeanDecreaseGini, decreasing = T),]
varImp_acc <- varImp_RF[order(varImp_RF$MeanDecreaseAccuracy, decreasing = T),]

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


## Check the correlation between variables
varimp_15 <- varImp_gini[1:15,]
genus_rf <- model_salmonella[, rownames(varimp_15)]
corr <- cor(genus_rf)
library(corrplot)
corrplot(corr, method="color", order = "AOE", type="lower", tl.col = "black", cl.ratio =0.2, tl.srt =45, tl.cex =0.7)

## Check the relative abundance of each variables
boxplot(genus_rf, col=rainbow(20), las=2)

