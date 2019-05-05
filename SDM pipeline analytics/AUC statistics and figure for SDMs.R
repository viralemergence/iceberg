
setwd('D:/ICEBERG/RawENMs/PPM/Statistics')
file.list <- list.files()

train <- c()
test <- c()

for(file in file.list) {
  data <- read.csv(file)
  train <- c(train,data$cv.mean.train.auc[1])
  test <- c(test,data$cv.mean.test.auc[1])
}

mean(train)
mean(test)
par(mfrow=c(2,1))
hist(train, main='Model validation on training data', xlab='AUC',
     xlim=c(0,1))
hist(test, main='Model validation on test data', xlab='AUC',
     xlim=c(0,1))
