library(matrix)
library(ROCR)

# WHERE DO HOSTADJ AND ALLSUMS COME FROM 

true.names <- rownames(HostAdj)
pred.names <- rownames(AllSums)

goodnames <- intersect(true.names, pred.names)

true <- HostAdj[rownames(HostAdj) %in% goodnames,
                colnames(HostAdj) %in% goodnames]

pred <- AllSums[rownames(AllSums) %in% goodnames,
                colnames(AllSums) %in% goodnames]

df <- data.frame(as.vector(true[lower.tri(true)]), as.vector(pred[lower.tri(pred)]))
colnames(df) <- c('observed','predicted')

df$observed[df$observed>1] = 1

pred <- prediction(df$predicted, df$observed)
performance(pred,"auc")
plot(performance(pred, "tpr", "fpr"),main='Receiver operator curve (AUC = 0.782)')
abline(0,1,col='red')

perf.tss <- performance(pred,"sens","spec")
tss.list <- (perf.tss@x.values[[1]] + perf.tss@y.values[[1]] - 1)
tss.df <- data.frame(alpha=perf.tss@alpha.values[[1]],tss=tss.list)
plot(tss.df,type='l')

thresh <- tss.df$alpha[which(tss.df$tss==max(tss.df$tss))]
perf.tss@y.values[[1]][which(perf.tss@alpha.values[[1]]==thresh)] # Sensitivity
perf.tss@x.values[[1]][which(perf.tss@alpha.values[[1]]==thresh)] # Specificity
1-perf.tss@y.values[[1]][which(perf.tss@alpha.values[[1]]==thresh)] # Type I error rate
1-perf.tss@x.values[[1]][which(perf.tss@alpha.values[[1]]==thresh)] # Type II error rate
