breastcancer=breast_cancer_wisconsin[,-1]
colnames(breastcancer) <- c("ClumpThickness","CellSizeUniformity",
                            "CellShapeUnformity","MarginalAdhesion",
                            "SingleEpithelialCellSize","BareNuclei",
                            "BlandChromatin","NormalNucleoli","Mitoses","Class")

breastcancer[, 6] <- sapply(breastcancer[, 6], as.integer)
#breastcancer[, 1] <- sapply(breastcancer[, 1], as.factor)
breastcancer[, 10] <- breastcancer[, 10]/2
breastcancer[, 10] <- breastcancer[, 10]-1
breastcancer=breastcancer[!(is.na(breastcancer$BareNuclei) | breastcancer$BareNuclei==""), ]
View(breastcancer)


##########################EDA code

rbind(
  apply(cancer.train[,], 2, min),
  apply(cancer.train[,], 2, median),
  apply(cancer.train[,], 2, max),
  apply(cancer.train[,], 2, mean),
  apply(cancer.train[,], 2, sd)
)
##boxplot
boxplot(cancer.train[,-10], notch=T, col = c("red","blue","green","yellow","white",
                                             "pink","grey","brown","orange"))

#clumpthickness plot
plot(density(cancer.train$ClumpThickness))
hist(cancer.train$ClumpThickness, prob=T, col="green", main="Histogram and Density of Clump Thickness", xlim=c(0,10), xlab="Clump Thickness")
lines(density(cancer.train$ClumpThickness), col="red", lwd=2)

abline(v=mean(cancer.train$ClumpThickness), col="blue", lty=2, lwd=1.5)

#Bland Chromatin Plot
plot(density(cancer.train$BlandChromatin))
hist(cancer.train$BlandChromatin, prob=T, col="yellow", main="Histogram and Density of Bland Chromatin", xlim=c(0,10), xlab="Bland Chromatin")
lines(density(cancer.train$BlandChromatin), col="red", lwd=2)

abline(v=mean(cancer.train$BlandChromatin), col="blue", lty=2, lwd=1.5)


## Cell size uniformity plot
plot(density(cancer.train$CellSizeUniformity))
hist(cancer.train$CellSizeUniformity, prob=T, col="orange", main="Histogram and Density of Cell Size Uniformity", xlim=c(0,10), xlab="Cell Size Uniformity")
lines(density(cancer.train$CellSizeUniformity), col="red", lwd=2)
abline(v=mean(cancer.train$CellSizeUniformity), col="blue", lty=2, lwd=1.5)

##Cell shape unformity plot
plot(density(cancer.train$CellShapeUnformity))
hist(cancer.train$CellShapeUnformity, prob=T, col="pink", main="Histogram and Density of Cell Shape Unformity", xlim=c(0,10), xlab="Cell Shape Uniformity")
lines(density(cancer.train$CellShapeUnformity), col="red", lwd=2)
abline(v=mean(cancer.train$CellShapeUnformity), col="blue", lty=2, lwd=1.5)


##scatterplot

pairs(cancer.train[,1:5])

pairs(cancer.train[,5:9])

plot(cancer.train$CellSizeUniformity, cancer.train$CellShapeUnformity, col= "blue")
abline(lm(cancer.train$CellSizeUniformity~cancer.train$CellShapeUnformity), col="red")

plot(cancer.train$ClumpThickness, cancer.train$MarginalAdhesion, col= "blue")
abline(lm(cancer.train$ClumpThickness~cancer.train$MarginalAdhesion), col="red")





################################ Losgistic regression and clasiffication tree

## random sample training and testing 
index<- sample(nrow(breastcancer), 0.8*nrow(breastcancer))
breastcancer.train<- breastcancer[index,]
breastcancer.test<- breastcancer[-index, ]

### build model
breastcancer.glm0<- glm(Class~., family=binomial, data=breastcancer.train)
## more details
sum.glm0<- summary(breastcancer.glm0)
sum.glm0$deviance
AIC(breastcancer.glm0)
BIC(breastcancer.glm0)

#### Prediction
pred.glm0<- predict(breastcancer.glm0, newdata = breastcancer.train, type = "response")
hist(pred.glm0)
boxplot(pred.glm0)


### classification
# simple way to choose cut-off probability -- sample proportion of "1"
pcut1<- mean(breastcancer.train$Class)
class.glm0<- (pred.glm0>pcut1)*1

table(class.glm0)

### confusion matrix
table(breastcancer.train$Class, class.glm0, dnn=c("Truth","Predicted"))

### error rate
FP=sum((breastcancer.train$Class==0 & class.glm0==1))
FN=sum((breastcancer.train$Class==1 & class.glm0==0))
FPR=sum((breastcancer.train$Class==0 & class.glm0==1))/sum(breastcancer.train$Class==0)
FNR=sum((breastcancer.train$Class==1 & class.glm0==0))/sum(breastcancer.train$Class==1)

# misclassification rate
(FP+FN)/nrow(breastcancer.train)
# another way, much easier!
mean(breastcancer.train$Class!=class.glm0)


##############################3

### write a cost function
cost<- function(obs, pred.p, pcut){
  class.p<- (pred.p>pcut)*1
  
  # define weight (managerial)
  weight.FP<- 1
  weight.FN<- 5
  
  #FP<- sum((obs==0 & class.p==1))
  #FN<- sum((obs==1 & class.p==0))
  #cost<- weight.FP*FP+weight.FN*FN
  
  FP<- (obs==0) & (class.p==1)
  FN<- (obs==1) & (class.p==0)
  cost<- mean(weight.FP*FP+weight.FN*FN)
  return(cost)
}

# apply the function (input values)
cost(obs=breastcancer.train$Class, pred.p=pred.glm0, pcut=pcut1) #this is the cost assoiated to using the pcut from symmetric cost


### choose optimal cutoff
pcut.seq<- seq(0, 1, length.out = 100)
cost.seq<- rep(0, length(pcut.seq))
for(i in 1:length(pcut.seq)){
  cost.seq[i]<- cost(breastcancer.train$Class, pred.glm0, pcut.seq[i])
}
# plot
plot(pcut.seq, cost.seq)
# find the optimal cut-off
optimal.pcut.glm0 = pcut.seq[which(cost.seq==min(cost.seq))]

# use this optimal cut-off to get confusion matrix and cost (weighted MR)
class.glm0.opt<- (pred.glm0>optimal.pcut.glm0)*1
table(breastcancer.train$Class, class.glm0.opt, dnn = c("Truth","Predicted"))
MR1=12/(350+183+12)
# cost (apply the cost function)
cost(breastcancer.train$Class, pred.glm0, optimal.pcut.glm0)

#so as we can see by using the optimal pcut, the cost that we get is lower than when we use the pcut from symetric cost

######################
## ROC
install.packages('verification')
library(verification)
roc.plot(breastcancer.train$Class == '1', pred.glm0)$roc.vol

install.packages('ROCR')
library(ROCR)
pred <- prediction(pred.glm0, breastcancer.train$Class)
perf <- performance(pred, "tpr", "fpr")
plot(perf, colorize=TRUE)
unlist(slot(performance(pred, "auc"), "y.values"))



############################
## out-of-sample performance for full model
pred.glm0.test<- predict(breastcancer.glm0, newdata = breastcancer.test, type = "response")


## use the same optimal cut-off from training sample
# confusion matrix
class.glm0.test<- (pred.glm0.test>optimal.pcut.glm0)*1
table(breastcancer.test$Class, class.glm0.test, dnn = c("Truth","Predicted"))
MR2=5/(78+54+5)
# asymmetric cost (weighted MR)
cost(breastcancer.test$Class, pred.glm0.test, optimal.pcut.glm0)

## out of sample roc

#Roc curve
pred2 <- prediction(pred.glm0.test,breastcancer.test$Class)
perf2 <- performance(pred2, "tpr", "fpr")
plot(perf2, colorize=TRUE)
#Get the AUC (area under curve)
unlist(slot(performance(pred2, "auc"), "y.values"))
# the ROC curve is very important in the report


######
#CV







######## best model

#stepwise
breastcancer.glm.step <- step(breastcancer.glm0) 
sum.glm1<- summary(breastcancer.glm.step)
sum.glm1$deviance
AIC(breastcancer.glm.step)
BIC(breastcancer.glm.step)
## out-of-sample performance
pred.glm1<- predict(breastcancer.glm.step, newdata = breastcancer.train, type = "response")
# confusion matrix
pcut3=optimal.pcut.glm0 
class.glm2<- (pred.glm1>pcut3)*1
table(breastcancer.train$Class, class.glm2, dnn = c("Truth","Predicted"))

#ROC
pred <- prediction(pred.glm1, breastcancer.train$Class)
perf <- performance(pred, "tpr", "fpr")
plot(perf, colorize=TRUE)
#Get the AUC (area under curve)
unlist(slot(performance(pred, "auc"), "y.values"))
# the ROC curve is very important in the report

cost(breastcancer.train$Class, pred.glm1, optimal.pcut.glm0)


######out of sample performance
## out-of-sample performance for full model
pred.glm0.test1<- predict(breastcancer.glm.step, newdata = breastcancer.test, type = "response")

## use the same optimal cut-off from training sample
# confusion matrix
class.glm0.test1<- (pred.glm0.test1>optimal.pcut.glm0)*1
table(breastcancer.test$Class, class.glm0.test1, dnn = c("Truth","Predicted"))
MR2=5/(78+54+5)
# asymmetric cost (weighted MR)
cost(breastcancer.test$Class, pred.glm0.test, optimal.pcut.glm0)

#Roc curve
pred3 <- prediction(pred.glm0.test1,breastcancer.test$Class)
perf3 <- performance(pred3, "tpr", "fpr")
plot(perf3, colorize=TRUE)
#Get the AUC (area under curve)
unlist(slot(performance(pred3, "auc"), "y.values"))
# the ROC curve is very important in the report






#########

install.packages("rpart")
library(rpart)
#clasification tree

### classification tree
breastcancer.tree<- rpart(Class~., data = breastcancer.train, method = "class" )
breastcancer.tree
install.packages('rpart.plot')
library(rpart.plot)
prp(breastcancer.tree, digits = 4, extra = 1)

pred0<- predict(breastcancer.tree, type="class")
table(breastcancer.train$Class, pred0, dnn = c("True", "Pred"))

### asymmetric loss
breastcancer.tree<- rpart(Class~. , data = breastcancer.train, method = "class", 
                    parms = list(loss=matrix(c(0,5,1,0), nrow = 2)))
breastcancer.tree
prp(breastcancer.tree, digits = 4, extra = 1)

sum.glm2<- summary(breastcancer.tree)

### prediction
breastcancer.tree.pred<- predict(breastcancer.tree, type = "class")
table(breastcancer.train$Class, breastcancer.tree.pred, dnn = c("True", "Pred"))
MR4=22/(341+182)

## exercise
breastcancer.tree.pred0<- predict(breastcancer.tree, newdata=breastcancer.train, type = "prob")
unique(breastcancer.tree.pred0[,2])
## proportion of the 1 in terminal nodes


#### ROC curve out of sample
library(ROCR)
pred = prediction(breastcancer.tree.pred0[,2], breastcancer.train$Class)
perf = performance(pred, "tpr", "fpr")
plot(perf, colorize=TRUE)
slot(performance(pred, "auc"), "y.values")[[1]]

### cost
cost <- function(r, pi){
  weight1 = 5
  weight0 = 1
  c1 = (r==1)&(pi==0) #logical vector - true if actual 1 but predict 0
  c0 = (r==0)&(pi==1) #logical vector - true if actual 0 but predict 1
  return(mean(weight1*c1+weight0*c0))
}

cost(breastcancer.train$Class, breastcancer.tree.pred)




###############out of sample
# testing error (out of sample)
breastcancer.tree.pred.test<- predict(breastcancer.tree, newdata=breastcancer.test, type = "class")
table(breastcancer.test$Class, breastcancer.tree.pred.test, dnn = c("True", "Pred"))
MR5=6/(77+6+54)
## exercise
breastcancer.tree.pred.test0<- predict(breastcancer.tree, newdata=breastcancer.test, type = "prob")
unique(breastcancer.tree.pred.test0[,2])
## proportion of the 1 in terminal nodes

#### ROC curve out of sample
library(ROCR)
pred = prediction(breastcancer.tree.pred.test0[,2], breastcancer.test$Class)
perf = performance(pred, "tpr", "fpr")
plot(perf, colorize=TRUE)
slot(performance(pred, "auc"), "y.values")[[1]]


### cost
cost <- function(r, pi){
  weight1 = 5
  weight0 = 1
  c1 = (r==1)&(pi==0) #logical vector - true if actual 1 but predict 0
  c0 = (r==0)&(pi==1) #logical vector - true if actual 0 but predict 1
  return(mean(weight1*c1+weight0*c0))
}

cost(breastcancer.test$Class, breastcancer.tree.pred.test)




######## Pruned tree



breastcancer.largetree<- rpart(Class~. , data = breastcancer.train, cp=0.001)
prp(breastcancer.largetree, digits = 4, extra = 1)
plotcp(breastcancer.largetree)
breastcancer.largetree$cptable

pruned.tree<- prune(breastcancer.largetree, cp=0.069)
prp(pruned.tree, digits = 4, extra = 1)

###########

pruned.tree<- rpart(Class~., data = breastcancer.train,method = "class", cp=0.069)
prp(pruned.tree, digits = 4, extra = 1)







################### 
### prediction

pruned.tree.pred<- predict(pruned.tree, type = "class")
table(breastcancer.train$Class, pruned.tree.pred, dnn = c("True", "Pred"))
MR6=(35)/(352+158+35)
# testing error
pruned.tree.pred.test<- predict(pruned.tree, newdata=breastcancer.test, type = "class")
table(breastcancer.test$Class, pruned.tree.pred.test, dnn = c("True", "Pred"))
MR5=13/(80+44+12)
## exercise
pruned.tree.pred.test0<- predict(pruned.tree, newdata=breastcancer.test, type = "prob")
unique(pruned.tree.pred.test0[,2])
## proportion of the 1 in terminal nodes

## exercise
pruned.tree.pred.test1<- predict(pruned.tree, newdata=breastcancer.train, type = "prob")
unique(pruned.tree.pred.test1[,2])
## proportion of the 1 in terminal nodes
#### ROC curve
library(ROCR)
pred = prediction(pruned.tree.pred.test1[,2], breastcancer.train$Class)
perf = performance(pred, "tpr", "fpr")
plot(perf, colorize=TRUE)
slot(performance(pred, "auc"), "y.values")[[1]]


### cost
cost <- function(r, pi){
  weight1 = 5
  weight0 = 1
  c1 = (r==1)&(pi==0) #logical vector - true if actual 1 but predict 0
  c0 = (r==0)&(pi==1) #logical vector - true if actual 0 but predict 1
  return(mean(weight1*c1+weight0*c0))
}

cost(breastcancer.train$Class, pruned.tree.pred)





#### ROC curve  test
library(ROCR)
pred = prediction(pruned.tree.pred.test0[,2], breastcancer.test$Class)
perf = performance(pred, "tpr", "fpr")
plot(perf, colorize=TRUE)
slot(performance(pred, "auc"), "y.values")[[1]]


### cost
cost <- function(r, pi){
  weight1 = 5
  weight0 = 1
  c1 = (r==1)&(pi==0) #logical vector - true if actual 1 but predict 0
  c0 = (r==0)&(pi==1) #logical vector - true if actual 0 but predict 1
  return(mean(weight1*c1+weight0*c0))
}

cost(breastcancer.test$Class, pruned.tree.pred.test)




######################################
## Cross-validation
#  rewrite the cost function that pcut cannot be an input parameter
#  this is because cv.glm requires that the cost function can only have two vector inputs
#  need pre-specify pcut
pcut = optimal.pcut.glm0
cost<- function(obs, pred.p){
  class.p<- (pred.p>pcut)*1
  
  # define weight (managerial)
  weight.FP<- 1
  weight.FN<- 5
  
  FP<- (obs==0) & (class.p==1)
  FN<- (obs==1) & (class.p==0)
  cost<- mean(weight.FP*FP+weight.FN*FN)
  return(cost)
}

# perform cross-validation using FULL DATA
library(boot)
breastcancer.glm1<-glm(Class~.,family=binomial,breastcancer.train);  
cv.result = cv.glm(data=breastcancer.train, glmfit=breastcancer.glm1, cost=costfunc, K=10) 
cv.result$delta[2]



####################################  Advanced tree models


#Bagging
#install.packages("ipred")
library(ipred)

ntree<- c(1, 3, 5, seq(10, 200, 20))
MSE.test<- rep(0, length(ntree))
for(i in 1:length(ntree)){
  cancer.bag1<- bagging(Class~.-ID, data = cancer.train, nbagg=ntree[i])
  cancer.bag.pred1<- predict(cancer.bag1, newdata = cancer.test)
  MSE.test[i]<- mean((cancer.test$Class-cancer.bag.pred1)^2)
}
plot(ntree, MSE.test, type = 'l', col=2, lwd=2)

#Prediction error on bagged tree
cancer.bag<- bagging(Class~.-ID, data = cancer.train, nbagg=100)
#cancer.bag

cancer.bag.pred<- predict(cancer.bag, newdata = cancer.test)
mean((cancer.test$Class-cancer.bag.pred)^2)

#Compared to one tree
library(rpart)
cancer.tree<- rpart(Class~.-ID, data = cancer.train)
cancer.tree.pred<- predict(cancer.tree, newdata = cancer.test)
mean((cancer.test$Class-cancer.tree.pred)^2)

#####
##OOB
cancer.bag.oob<- bagging(Class~.-ID, data = cancer.train, coob=T, nbagg=length(breastcancer))
cancer.bag.oob

#Comparison with 150 tree
cancer.bag.oob1<- bagging(Class~.-ID, data = cancer.train, coob=T, nbagg=150)
cancer.bag.oob1

#####
##Random Forests

breastcancer$Class <- as.factor(breastcancer$Class)
cancer.train$Class <- as.factor(cancer.train$Class)
cancer.test$Class <- as.factor(cancer.test$Class)
#install.packages("randomForest")
library(randomForest)

cancer.rf<- randomForest(Class~.-ID, data = cancer.train, type="classification")
cancer.rf

##plot
plot(cancer.rf)
legend("right", legend = c("OOB Error", "FPR", "FNR"),
       lty = c(1,2,3), col = c("black", "red", "green"))

#cost function
cancer.rf.pred<- predict(cancer.rf, type = "prob")[,2]
costfunc = function(obs, pred.p, pcut){
  weight1 = 5   # define the weight for "true=1 but pred=0" (FN)
  weight0 = 1    # define the weight for "true=0 but pred=1" (FP)
  c1 = (obs==1)&(pred.p<pcut)    # count for "true=1 but pred=0"   (FN)
  c0 = (obs==0)&(pred.p>=pcut)   # count for "true=0 but pred=1"   (FP)
  cost = mean(weight1*c1 + weight0*c0)  # misclassification with weight
  return(cost) # you have to return to a value when you write R functions
} 
p.seq = seq(0.01, 0.5, 0.01)
cost = rep(0, length(p.seq))  
for(i in 1:length(p.seq)){ 
  cost[i] = costfunc(obs = cancer.train$Class, pred.p = cancer.rf.pred, pcut = p.seq[i])  
}
plot(p.seq, cost)

optimal.pcut = min(p.seq[which(cost==min(cost))])

#Confusion Matrix on Training Sample
cancer.rf.class<- (cancer.rf.pred>optimal.pcut)*1
table(cancer.train$Class, cancer.rf.class, dnn = c("True", "Pred"))

#Confusion Matrix on Testing Sample
## out-of-sample
cancer.rf.pred.test<- predict(cancer.rf, newdata=cancer.test, type = "prob")[,2]
cancer.rf.class.test<- (cancer.rf.pred.test>optimal.pcut)*1
table(cancer.test$Class, cancer.rf.class.test, dnn = c("True", "Pred"))


#####Error count on Training Sample
TFP = sum((cancer.train$Class==0 & cancer.rf.class==1))
TFN = sum((cancer.train$Class==1 & cancer.rf.class==0))

###Misclassification Rates
sum(TFP)/sum(nrow(cancer.train))
sum(TFN)/sum(nrow(cancer.train))
#misclassification rate
sum((TFP)+(TFN))/sum(nrow(cancer.train))


#####Error count on Testing Sample
FP = sum((cancer.test$Class==0 & cancer.rf.class.test==1))
FN = sum((cancer.test$Class==1 & cancer.rf.class.test==0))

#Misclassification Rates
sum(FP)/sum(nrow(cancer.test))
sum(FN)/sum(nrow(cancer.test))
#misclassification rate
sum((FP)+(FN))/sum(nrow(cancer.test))



##testing different inputs
#default settings were found to be optimal
cancer.rf2<- randomForest(Class~.-ID, data = cancer.train, type="classification", ntree=500, mtry=3)
cancer.rf2

plot(cancer.rf2)
legend("right", legend = c("OOB Error", "FPR", "FNR"),
       lty = c(1,2,3), col = c("black", "red", "green"))

cancer.glm<-glm(Class~.-ID,family=binomial,cancer.train)
summary(cancer.glm)














