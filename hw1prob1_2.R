data <- read.csv("ozone.csv")
require(MASS)
######
#Problem 1
X <- as.matrix(data[,-1]);X.orig <- X; X <- scale(X,center=T,scale=T)
Y <- as.numeric(data[,1]); Y.mean <- mean(Y); Y <- Y - Y.mean

##Bound data proportionaly
n = nrow(data)
bound_train <- floor((nrow(data)/5)*3)
bound_valid <- floor((nrow(data)/5)*4)
bound_test <- floor((nrow(data)/5)*5)


#Randomly split data 60/20/20
df <- cbind(Y,X)[sample(1:n, replace=FALSE),]

df.train = as.matrix(df[1:bound_train,])
df.valid = as.matrix(df[(bound_train+1):bound_valid,])
df.test <- as.matrix(df[(bound_valid+1):bound_test,])
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]

(beta.OLS <- solve(t(X.tr)%*%X.tr)%*%t(X.tr)%*%Y.tr) #Yes, it matches lm with scaled data!!! good.
OLS.pred <- df.valid[,-1]%*%beta.OLS
(MSE.OLS <- sum((df.valid[,1] - OLS.pred)^2)/nrow(df.valid))

####Implement Ridge
penalty_vec <- seq(.01,100, by= .1 )

MSE.R_vec <- c()

for (lam in penalty_vec) {
  beta.R <- solve(t(X.tr)%*%X.tr + diag(rep(lam,8)))%*%t(X.tr)%*%Y.tr
  resid.R <- df.valid[,1] -  (df.valid[,-1] %*% beta.R)
  MSE.R <- sum(resid.R^2)/nrow(df.valid)
  MSE.R_vec <- c(MSE.R_vec, MSE.R)
}

min(MSE.R_vec)
#Graphics
require(ggplot2)
prob1_df <- data.frame(lambda = penalty_vec, MSE.Ridge = MSE.R_vec, MSE.LS = rep(MSE.R_vec[1],times = length(MSE.R_vec)))
ggplot(prob1_df) + geom_line(aes(x = lambda, y = MSE.Ridge, color = 'MSE.Ridge')) +
  geom_line(aes(x = lambda, y = MSE.LS, color = 'MSE.LS')) + ylab("MSE") +
  ggtitle("MSE Existence Theorem on O-zone data")+ labs(color='Normalized MSE')

#####
#Problem 2


##### 
#MANUAL BUILD OF KERNEL FUNCTION!!!! Matches with package results quite closely, but validation on more than 1 parameter is too slow.
#Building out a Kernel Ridge Regression system
#####

kern.funct <- function(z1 = "pvec1",z2 = "pvec2", type = c("IP", "Quad", "Poly", "Exp"), d = 3, c = 1,sigma = 3) {
  if(type == "IP") {
    return(z1 %*% z2)
  }
  else if (type == "Quad") {
    return((z1 %*% z2)^2)
  }
  else if (type == "Poly") {
    return((c + z1 %*% z2)^d)
  }
  else if (type == "Exp") {
    w <- z1 - z2
    return(exp(-(w %*% w)^2/sigma^2))
  }
}


kern.matrix <- function(data = c(), kern.type = "IP", data.new = NA, sigma = 1) {
  
  if(is.matrix(data.new) == F) {
    n <- nrow(data)
    mat <- matrix(nrow=n, ncol = n) #Holder matrix
    for (i in 1:n) {
      for(j in 1:n){
        mat[i,j] <- kern.funct(data[i,], data[j,], kern.type, sigma = sigma) #kernel function on p-vectors
      }
    }
    return(mat)
  }
  
  else{
    rows = nrow(data.new)
    cols = nrow(data)
    
    mat <- matrix(nrow=rows, ncol = cols) #Holder matrix
    for (i in 1:rows) {
      for(j in 1:cols){
        mat[i,j] <- kern.funct(data.new[i,], data[j,], kern.type, sigma = sigma) #kernel function on p-vectors
      }
    }
    return(mat)
  }
  
}

kern.reg <- function(data = df.train, data.pred = df.valid, type = "IP", penalty = penalty_vec, sigma = 1) {
  
  #Fit the model with train, check the kernels with valid and report the best model.
  #RIDGE
  X <- as.matrix(data[,-1])
  Y <- as.numeric(data[,1])
  K <- kern.matrix(data = X, kern.type = type, sigma = sigma)
  n <- nrow(X)
  MSE.pred_vec <- c()
  for(lambda in penalty) {
    alpha <- solve(K + diag(rep(lambda, times = nrow(K))))%*%Y
    K.pred <- kern.matrix(data = X, kern.type = type, data.new = data.pred[,-1], sigma = sigma)
    Y.pred <- K.pred%*%alpha #prediction on validation set
    MSE.pred <- sum((data.pred[,1]-Y.pred)^2)/nrow(data.pred)
    MSE.pred_vec <- c(MSE.pred_vec,MSE.pred)
  }
  return(list(MSE = MSE.pred_vec, Y.predict = Y.pred, lambda = penalty))
}

#####
#Normalize data
X <- as.matrix(data[,-1]); X = scale(X,center=T,scale=T)
Y <- as.numeric(data[,1]); Y.mean = mean(Y); Y = Y - Y.mean

##Bound data proportionaly
n = nrow(data)
bound_train <- floor((nrow(data)/5)*3)
bound_valid <- floor((nrow(data)/5)*4)
bound_test <- floor((nrow(data)/5)*5)


require(glmnet)
penalty_vec <- seq(.01,500, by= 10 )

#Hold pred MSE's
ridge.t_errors <- c() 
quad.t_errors <- c()
poly.t_errors <- c()
exp.t_errors <- c()

results <- data.frame(ridge = ridge.t_errors, quad = quad.t_errors, poly = poly.t_errors, exp = exp.t_errors)
colnames(results) <- c("Inner Product", "Simple Poly", "Full Poly", "Gaussian/RBF")
colMeans(results)


for (iter in 1:10) {
  
  #Randomly split data 60/20/20
  df <- cbind(Y,X)[sample(1:n, replace=FALSE),]
  
  df.train = as.matrix(df[1:bound_train,])
  df.valid = as.matrix(df[(bound_train+1):bound_valid,])
  df.test = as.matrix(df[(bound_valid+1):bound_test,])
  
  ######
  
#####  
  #Ridge/Inner product Kernel
  fit.ridge <- kern.reg() 
  fit.ridge$MSE
  ridge.lambda <- penalty_vec[which(fit.ridge$MSE == min(fit.ridge$MSE))] #Choose lambda from validation set
  (ridge.t_error <- kern.reg(penalty = ridge.lambda, data.pred = df.test)$MSE) #Gets the MSE on the unseen test data
  ridge.t_errors <- c(ridge.t_errors, ridge.t_error)
  
  #Simple/Naive quadratic kernel
  fit.quad <- kern.reg(type = "Quad") 
  fit.quad$MSE
  quad.lambda <- penalty_vec[which(fit.quad$MSE == min(fit.quad$MSE))]
  (quad.t_error <- kern.reg(penalty = quad.lambda, data.pred = df.test)$MSE)
  quad.t_errors <- c(quad.t_errors, quad.t_error)
  
  
  #Multipolynomial
  fit.poly <- kern.reg(type = "Poly") 
  fit.poly$MSE
  poly.lambda <- penalty_vec[which(fit.poly$MSE == min(fit.poly$MSE))]
  (poly.t_error <- kern.reg(penalty = poly.lambda, data.pred = df.test)$MSE)
  poly.t_errors <- c(poly.t_errors, poly.t_error)
  #Gaussian/RBF
  fit.exp <- kern.reg(type = "Exp")
  #RBF kernel: Matrix tuning
  sigmas <- seq(.1, 8, by = .5)
  tune.RBF <- matrix(nrow = length(sigmas),ncol = length(penalty_vec))
  for(i in sigmas) {
    tune.RBF[which(sigmas == i),] <- kern.reg(type = "Exp", sigma= i)$MSE
  }
  
  tune.index <- which(tune.RBF == min(tune.RBF), arr.ind = T) #index of best parameter set (sigma/lambda)
  tune.param <- c(sigmas[tune.index[1]], penalty_vec[tune.index[2]]) #best parameter set (sigma/lambda)
  
  (exp.t_error <- kern.reg(penalty = tune.param[2], data.pred = df.test, sigma = tune.param[1])$MSE)
  exp.t_errors <- c(exp.t_errors, exp.t_error)
}


#####
#Problem 2b

#####
#Normalize data
data.normalize <- function() {
  X <<- as.matrix(data[,-1]);X.orig <<- X; X <<- scale(X,center=T,scale=T)
  Y <<- as.numeric(data[,1]); Y.mean <<- mean(Y); Y <<- Y - Y.mean
  
  ##Bound data proportionaly
  n = nrow(data)
  bound_train <- floor((nrow(data)/5)*3)
  bound_valid <- floor((nrow(data)/5)*4)
  bound_test <- floor((nrow(data)/5)*5)
  
  
  #Randomly split data 60/20/20
  df <- cbind(Y,X)[sample(1:n, replace=FALSE),]
  
  df.train = as.matrix(df[1:bound_train,])
  df.valid = as.matrix(df[(bound_train+1):bound_valid,])
  df.train <<- rbind(df.train,df.valid) #We'll be using CV for part 2B, so combine to get 80% training/validation data
  df.test <<- as.matrix(df[(bound_valid+1):bound_test,])
}

X <- as.matrix(data[,-1]);X.orig = X; X = scale(X,center=T,scale=T)
Y <- as.numeric(data[,1]); Y.mean = mean(Y); Y = Y - Y.mean

##Bound data proportionaly
n = nrow(data)
bound_train <- floor((nrow(data)/5)*3)
bound_valid <- floor((nrow(data)/5)*4)
bound_test <- floor((nrow(data)/5)*5)


#Randomly split data 60/20/20
df <- cbind(Y,X)[sample(1:n, replace=FALSE),]

df.train = as.matrix(df[1:bound_train,])
df.valid = as.matrix(df[(bound_train+1):bound_valid,])
df.train = rbind(df.train,df.valid) #We'll be using CV for part 2B, so combine to get 80% training/validation data
df.test = as.matrix(df[(bound_valid+1):bound_test,])

###Hold Results
LS.results <- c()
ridge.results <- c()
subsets.results <- c()
lasso.results <- c()
enet.results <- c()
alasso.results <- c()
mc.results <- c()

results.trad_linear <- data.frame(LS = LS.results, Ridge = ridge.results, BestSubsets = subsets.results,
                                  Lasso = lasso.results, ElasticNet = enet.results, AdaptiveLasso = alasso.results, MCplus = mc.results)
write.csv(results.trad_linear,"trad_linear_pred_results2.csv")
colMeans(results.trad_linear)
#Least Squares
for(i in 1:10){
data.normalize()
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]

(beta.OLS <- solve(t(X.tr)%*%X.tr)%*%t(X.tr)%*%Y.tr) #Yes, it matches lm with scaled data!!! good.
OLS.pred <- df.test[,-1]%*%beta.OLS

(MSE.OLS <- sum((df.test[,1] - OLS.pred)^2)/nrow(df.test))
LS.results <- c(LS.results, MSE.OLS)
}
#Ridge
library(glmnet)
for(i in 1:10){
data.normalize()
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]

ridge.lambda <- cv.glmnet(x = X.tr, y = Y.tr, alpha = 0)$lambda.min
ridge.fit <- glmnet(x = X.tr, y = Y.tr, alpha = 0, lambda = ridge.lambda)
ridge.pred <- predict(ridge.fit, newx = df.test[,-1])
(MSE.ridges <- sum((df.test[,1] - ridge.pred)^2)/nrow(df.test))
ridge.results <- c(ridge.results, MSE.ridges)

coef(ridge.fit)

ridge.plot <- glmnet(x = X.tr, y = Y.tr, alpha = 0)
plot(ridge.plot)
}
#Best Subsets
library(leaps)
for(i in 1:10){
data.normalize()
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]

subsets <- leaps(X.tr, Y.tr) #iterate through all subsets and calcualting mallow Cp
subsets.index <- subsets$which[which(subsets$Cp == min(subsets$Cp)),] #Select best subset based on cp
X.BS <- X.tr[,subsets.index] #Model matrix selected by Best Subsets
(beta.BS <- solve(t(X.BS)%*%X.BS)%*%t(X.BS)%*%Y.tr)
X.BS_test <- df.test[,-1]; X.BS_test <- X.BS_test[,subsets.index]
BS.pred <- X.BS_test%*%beta.BS
(MSE.BS <- sum((df.test[,1] - BS.pred)^2)/nrow(df.test))
subsets.results <- c(subsets.results, MSE.BS)
}
#Lasso
for(i in 1:10){

data.normalize()
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]

lasso.lambda <- cv.glmnet(x = X.tr, y = Y.tr, alpha = 1)$lambda.min
lasso.fit <- glmnet(x = X.tr, y = Y.tr, alpha = 1, lambda = lasso.lambda)
lasso.pred <- predict(lasso.fit, newx = df.test[,-1])
(MSE.lasso <- sum((df.test[,1] - lasso.pred)^2)/nrow(df.test))

coef(lasso.fit)

lasso.results <- c(lasso.results, MSE.lasso)
testt <- glmnet(x = X.tr, y = Y.tr, alpha = 1)


plot(testt)
}
#Elastic Net
require(caret)
for(i in 1:10){
data.normalize()
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]

lambda.grid <- 10^seq(2, -2, length= 40)
alpha.grid <- seq(0,1,length = 10)

#CV method for train function
trnCtrl = trainControl(method = "repeatedCV", number = 10, repeats = 5)

#Set up search grid
srchGrd = expand.grid(.alpha = alpha.grid, .lambda = lambda.grid)

#Tune parameters
netFit <- train(x = X.tr, y = Y.tr, method = 'glmnet', trControl = trnCtrl, tuneGrid = srchGrd)
net.alpha <- netFit$bestTune[1]; net.lambda <- netFit$bestTune[2]

#Get predictions
enet.fit <- glmnet(x = X.tr, y = Y.tr, alpha = net.alpha, lambda = net.lambda)
enet.pred <- predict(enet.fit, newx = df.test[,-1])
(MSE.net <- sum((df.test[,1] - enet.pred)^2)/nrow(df.test))
enet.results <- c(enet.results, MSE.net)

enet.plot <- glmnet(x = X.tr, y = Y.tr, alpha = net.alpha)
plot(enet.plot)

}

#Adaptive Lasso
require(parcor)
for(i in 1:10) {
data.normalize()
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]

alasso.coefs <- adalasso(X.tr, Y.tr, k = 10, use.Gram = TRUE,both=TRUE)$coefficients.adalasso
alasso.index <- (which(alasso.coefs != 0))

Xtr.alasso <- X.tr[,alasso.index] #Model matrix selected by alasso
(beta.alasso <- solve(t(Xtr.alasso)%*%Xtr.alasso)%*%t(Xtr.alasso)%*%Y.tr)
X.test_alasso <- df.test[,-1]; X.test_alasso <- X.test_alasso[,alasso.index]
alasso.pred <- X.test_alasso%*%beta.alasso
(MSE.alasso <- sum((df.test[,1] - alasso.pred)^2)/nrow(df.test))

alasso.results <- c(alasso.results, MSE.alasso)
}
#MC+
require(ncvreg)
for(i in 1:10){
data.normalize()
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]

mcp.lambda <- cv.ncvreg(X.tr,Y.tr,family="gaussian",penalty="MCP")$lambda.min
mcp.fit <- ncvreg(X.tr,Y.tr,family="gaussian",penalty="MCP",lambda= mcp.lambda)
mcp.pred <- predict(mcp.fit, X = df.test[,-1])
(MSE.mcp <- sum((df.test[,1] - mcp.pred)^2)/nrow(df.test))
mc.results <- c(mc.results, MSE.mcp)

coef(mcp.fit)

mcp.plot <-ncvreg(X.tr,Y.tr,family="gaussian",penalty="MCP")
plot(mcp.plot)
}

#Kernel ridge regression: Package implementation
library(CVST)

MSE.KRR_RBF <- c()
MSE.KRR_FullP <- c()
MSE.KRR_P <- c()
MSE.KRR_Inner <- c()

results.kernel <- data.frame(Gaussian = MSE.KRR_RBF, FullPoly = MSE.KRR_FullP , SimplePoly = MSE.KRR_P[-1], InnerProduct= MSE.KRR_Inner )
write.csv(results.kernel, "kernel_results_CLVT2.csv")
colMeans(results.kernel) #FullPoly is by far the best! d ~ 2 / ~ 36

##1) Gaussian Kernel
for(i in 1:10){
  
data.normalize()
X.tr <- df.train[,-1]; Y.tr <- df.train[,1]


d = constructData(x=X.tr, y=Y.tr) ## Structure data in CVST format
krr_learner = constructKRRLearner()   ## Build the base learner

#Create search grid 
param = constructParams(kernel="rbfdot", sigma=10^seq(2, -2, length= 10), lambda=4^seq(2, -2, length= 10))

#run CV
cv.results <- fastCV(d, krr_learner, params= param,constructCVSTModel())
kern.rbf_sigma = as.numeric(cv.results[[1]][2]) ; kern.rbf_lambda = as.numeric(cv.results[[1]][3])



cv.params = list(kernel="rbfdot", sigma= kern.rbf_sigma, lambda= kern.rbf_lambda)
krr_trained = krr_learner$learn(d, cv.params)

##To predict, format test data, 'dTest', the same way as in 'd'
dTest <- constructData(x=df.test[,-1], y=df.test[,1])
pred = krr_learner$predict(krr_trained, dTest)

MSE.KRR <- sum((df.test[,1] - pred)^2)/nrow(df.test)
MSE.KRR_RBF <- c(MSE.KRR_RBF, MSE.KRR)

}
##2) Full Poly Kernel
for(i in 1:10){
  
data.normalize()
d = constructData(x=X.tr, y=Y.tr) ## Structure data in CVST format
krr_learner = constructKRRLearner()   ## Build the base learner

#Create search grid [scale = gamma here, but we keep this as 1 to follow lecture notes as well as make life easier..]
param = constructParams(kernel="polydot", degree=seq(1, 5, length = 10), scale = 1, offset= 10^seq(2, -2, length= 10), lambda=10^seq(2, -2, length= 10))

#run CV
cv.results <- fastCV(d, krr_learner, params= param,constructCVSTModel())
kern.poly_degree = as.numeric(cv.results[[1]][2]) ; kern.poly_scale = as.numeric(cv.results[[1]][3]); kern.poly_offset = as.numeric(cv.results[[1]][4]);kern.poly_lambda = as.numeric(cv.results[[1]][5])

cv.params = list(kernel="polydot", degree = kern.poly_degree,scale = kern.poly_scale,offset = kern.poly_offset,lambda = kern.poly_lambda)
krr_trained = krr_learner$learn(d, cv.params)


##To predict, format test data, 'dTest', the same way as in 'd'
dTest <- constructData(x=df.test[,-1], y=df.test[,1])
pred = krr_learner$predict(krr_trained, dTest)

(MSE.KRR_poly <- sum((df.test[,1] - pred)^2)/nrow(df.test))
MSE.KRR_FullP <- c(MSE.KRR_FullP, MSE.KRR_poly)

}
##3) Simple/Naive Poly Kernel
for(i in 1:10){
  
data.normalize()
d = constructData(x=X.tr, y=Y.tr) ## Structure data in CVST format
krr_learner = constructKRRLearner()   ## Build the base learner

#Create search grid [scale = gamma here, but we keep this as 1 to follow lecture notes as well as make life easier..]
param = constructParams(kernel="polydot", degree=seq(1, 5, length = 20), scale = 1, offset= 0, lambda=10^seq(2, -2, length= 20))

#run CV
cv.results <- fastCV(d, krr_learner, params= param,constructCVSTModel())
kern.poly_degree = as.numeric(cv.results[[1]][2]) ; kern.poly_scale = as.numeric(cv.results[[1]][3]); kern.poly_offset = as.numeric(cv.results[[1]][4]);kern.poly_lambda = as.numeric(cv.results[[1]][5])

cv.params = list(kernel="polydot", degree = kern.poly_degree,scale = kern.poly_scale,offset = kern.poly_offset,lambda = kern.poly_lambda)
krr_trained = krr_learner$learn(d, cv.params)


##To predict, format test data, 'dTest', the same way as in 'd'
dTest <- constructData(x=df.test[,-1], y=df.test[,1])
pred = krr_learner$predict(krr_trained, dTest)

(MSE.KRR_poly <- sum((df.test[,1] - pred)^2)/nrow(df.test))
MSE.KRR_P <- c(MSE.KRR_P, MSE.KRR_poly)

}

##4) Inner Product Poly Kernel
for(i in 1:10)
{

data.normalize()
d = constructData(x=X.tr, y=Y.tr) ## Structure data in CVST format
krr_learner = constructKRRLearner()   ## Build the base learner

#Create search grid [scale = gamma here, but we keep this as 1 to follow lecture notes as well as make life easier..]
param = constructParams(kernel="polydot", degree=1, scale = 1, offset= 0, lambda=10^seq(2, -2, length= 100))

#run CV
cv.results <- fastCV(d, krr_learner, params= param,constructCVSTModel())
kern.poly_degree = as.numeric(cv.results[[1]][2]) ; kern.poly_scale = as.numeric(cv.results[[1]][3]); kern.poly_offset = as.numeric(cv.results[[1]][4]);kern.poly_lambda = as.numeric(cv.results[[1]][5])

cv.params = list(kernel="polydot", degree = kern.poly_degree,scale = kern.poly_scale,offset = kern.poly_offset,lambda = kern.poly_lambda)
krr_trained = krr_learner$learn(d, cv.params)


##To predict, format test data, 'dTest', the same way as in 'd'
dTest <- constructData(x=df.test[,-1], y=df.test[,1])
pred = krr_learner$predict(krr_trained, dTest)

(MSE.KRR_poly <- sum((df.test[,1] - pred)^2)/nrow(df.test))
MSE.KRR_Inner <- c(MSE.KRR_Inner, MSE.KRR_poly)
}


######Problem 2c: Visualisation of prediction comparisons

#load results data that we saved from above code.
trad_reg.results <- read.csv("trad_linear_pred_results3.csv")[,-1]
kern_reg.results <- read.csv("kernel_results_CLVT2.csv")[,-1] 

colMeans(kern_reg.results) #Full Poly is best kernel
colMeans(trad_reg.results) #Adaptive Lasso did the best out of these

results.df <- as.data.frame(cbind(trad_reg.results, kern_reg.results)) #Data frame of moving average
results.df <-  data.frame(cumsum(results.df)/(1:10), Trial = 1:10)
##Visualize
library(ggplot2)
library("reshape2")

#Plot comparison MSE's over time
results.mdf <- melt(results.df, id.vars="Trial", value.name="Avg_Prediction_MSE", variable.name="Method")
ggplot(results.mdf) + geom_line(aes(x = Trial, y = Avg_Prediction_MSE, factor = Method, color = Method)) +
  scale_x_continuous(breaks = 1:10) + ylab("Rolling Average Prediction MSE") + ggtitle("Predictive performance of regression methods")

#####
#Plot the fit of the best model (Kernel: Full Poly)
#saved: MSE= 10.11187 /lambda 100/ offset 35.93814 / d = 2.77778/ scale 1
best.preds <- pred #Just a saved instance of the best model
test <- df.test[,1]

best.result.df <- data.frame(pred = best.preds, test = test, num = 1:66)
best.result.mdf <- melt(best.result.df, id.vars="num", value.name="Pred", variable.name="Method")
ggplot(best.result.mdf) + geom_point(aes(x = num, y = Pred, factor = Method, color = Method)) + geom_line(aes(x = num, y = Pred, factor = Method, color = Method)) + 
  ylab("O-zone") + ggtitle("Illustrative fit of Best Model (Kernel Full Polynomial)") + labs(color="Result from:")+
  scale_color_manual(labels = c("Model Prediction", "Test Set"),values = c("blue", "red")) + xlab("Index of test observation")


######
#Problem 3
gene.X <- read.csv("X_SRBCT.csv"); gene.Y <-read.csv("Y_SRBCT.csv")
X <- gene.X; Y<- gene.Y

#Normalize data


