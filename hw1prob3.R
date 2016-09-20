#Stat 444 HW1 Problem 3
######
#Problem 3
library("ggplot2")
library("reshape2")
library("gtools")
library(glmnet)
library(ncvreg)
library(MASS)
#Plotter function
key_predictor_plot <- function(x, ..., num_predictors=10, alpha_min=0.4){
  stopifnot(require("ggplot2"))
  stopifnot(require("reshape2"))
  stopifnot(require("gtools"))
  
  if(inherits(x, "glmnet")){
    beta <- x$beta # Direct access of beta rather than coef() avoids intercept
  } else if(inherits(x, "ncvreg")){
    ## ncvreg includes an intercept term (unlike glmnet) so we omit it
    beta <- if (length(x$penalty.factor) == NROW(x$beta)){
      coef(x)
    } else {
      coef(x)[-1, , drop = FALSE] ## (intercept)
    }
  }
  lambda <- x$lambda
  
  active = glmnet::nonzeroCoef(beta)
  nactive = length(active)
  beta <- as.matrix(beta[active, , drop=FALSE])
  index = log(lambda)
  xlab=expression(log(lambda))
  
  ## Identify the first num_predictors to enter the model
  ## (or fewer, in the case when we add two variables simultaneously)
  key_predictors <- names(which(beta[,max(which(colSums(beta != 0) <= num_predictors))] != 0))
  if((length(key_predictors) > num_predictors) | (length(key_predictors) < 2)){
    ## If we have lots of predictors (e.g., ridge)
    ## just pull out the predictors with the largest |\hat{\beta}_j|
    ## and carry on
    key_predictors <- rownames(beta)[rank(abs(beta[,ncol(beta)]))][1:num_predictors]
  }
  
  
  ## Now we rearrange our data in a 'long' form preferred by ggplot2
  mb <- melt(beta)
  colnames(mb) <- c("Variable", "Step", "Beta")
  mb$Variable <- as.character(mb$Variable)
  mb <- cbind(mb, key_predictor=mb$Variable %in% key_predictors)
  mb <- cbind(mb, LogLambda=rep(index, each=nrow(beta)))
  mb <- cbind(mb, plot_color=ifelse(mb$key_predictor, mb$Variable, "Other"))
  
  ## And add some extra columns to use as plot parameters.
  g <- ggplot(mb, aes(x=LogLambda, y=Beta, group=Variable))
  if(!all(mb$key_predictor)){
    g <- g + geom_line(aes(col=plot_color, alpha=key_predictor))
    g <- g + scale_alpha_discrete(aes(alpha=key_predictor), guide="none", range=c(alpha_min, 1))
  } else {
    g <- g + geom_line(aes(col=plot_color))
  }
  
  g <- g + xlab(expression(log(lambda))) + ylab(expression(hat(beta[j])))
  g <- g + scale_color_discrete(guide=guide_legend(title="Predictor"), 
                                breaks=c(mixedsort(unique(mb$Variable)), "Other"))
  
  ## Conventionally, these are plotted with decreasing LogLambda
  g <- g + scale_x_reverse()
  
  g 
}
#####
#Test function
p <- 500
n <- 100
s <- 10

Sigma <- toeplitz(0.9^(0:(p-1)))
X <- mvrnorm(n, mu=rep(0, p), Sigma)

beta <- rep(0, p); 
beta[1:s] <- runif(s, 2, 3)

y <- X %*% beta + rnorm(n)

g1 <- glmnet(X, y)

g <- key_predictor_plot(g1) + theme_bw() + ggtitle("Lasso Regression Example")
print(g)

#####
#Start problem
gene.X <- read.csv("X_SRBCT.csv"); gene.Y <-read.csv("Y_SRBCT.csv")
X<- gene.X ; Y <- gene.Y
colnames(Y) = "p53"
#Normalize data
X <- as.matrix(X);X.orig <- X; X <- scale(X,center=T,scale=T)
Y <- as.numeric(as.matrix(Y)); Y.mean <- mean(Y); Y <- Y - Y.mean
#Elastic Net
require(caret)

enet.fit <- glmnet(x = X, y = Y, alpha = .5)
g1 <- key_predictor_plot(enet.fit) + theme_bw() + ggtitle("Elastic Net Regression Example")
print(g1)

#Lasso
lasso.fit <- glmnet(x = X, y = Y, alpha = 1)

g2 <- key_predictor_plot(lasso.fit) + theme_bw() + ggtitle("Lasso Regression Example")
print(g2)

#MC+

mcp.fit <- ncvreg(X,Y,family="gaussian",penalty="MCP")

g3 <- key_predictor_plot(mcp.fit) + theme_bw() + ggtitle("MC+ Regression Example")
print(g3)

