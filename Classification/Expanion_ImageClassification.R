rm(list=ls())

library(MASS)
library(glmnet)
library(RCurl)
library(ROCR)
library(fBasics) 
library(pROC)
library(randomForest)

Fname="Features_Exp"
setwd("R:/Beck Lab/_USERS/Octavian/Expansion Microscopy/RESULTS/_MAX projection SAMPLES/Analysis_2nd_Submission/Expansion")

f <- read.csv(paste(Fname,".csv",sep=""))
# Remove features of those images which are recommend to delete by pathologists (Images have bad quality and border line cases)
f <- f[f[,4]!=1,]
f <- f[f[,4]!=2,]
f <- f[f[,4]!=3,]
f <- f[f[,4]!=4,]

# Extract image labels from Image name
labels = as.character(f[,1])
labels = gsub("-", "_", labels)
y = unlist(lapply(strsplit(labels,split="_",fixed=T),function(x)(x[1])))
table(y)

x=f[,-c(1,2,3,4)]

# Remove those features which have NA or NAN or NULL or Zero values
sds=apply(x,2,sd)
x=x[,!is.na(sds)]
sds=apply(x,2,sd)
x=x[,!(sds==0)]

# Expansion Factor Normalization. Divide each features with image expansion factor
x2<-x
ef<-as.numeric(f[,3])
for (i in 1:nrow(x)){
  for (j in 1:ncol(x)){
    x2[i,j] = x[i,j]/ef[i]
  }
}
x2$label = y
write.csv(x2, 'Data_7Summary_4Classes.csv',row.names=FALSE)

# For Selecting 2 Summary Features only (Mean & Standard Deviation)
Features =  colnames(x)
Features.Cat = unlist(lapply(strsplit(Features,split="_",fixed=T),function(x)(x[1])))
x.Mean = x2[,(Features.Cat == "Mean")]
x.STD = x2[,(Features.Cat == "STD")]

x2 = cbind(x.Mean,x.STD)
x2$label = y
write.csv(x2, 'Data_2Summary_4Classes.csv', row.names=FALSE)

case = 'Data_7Summary_4Classes'
f=read.csv(paste(case,".csv",sep=""))

x = f[,-ncol(f)]
y = f[,ncol(f)]
table(y)

SelectedFolds = 6

##########################################################################################################
########################################   UDH vs ADH
##########################################################################################################
c1 = "UDH"
c2 = "ADH"
nfolds = SelectedFolds

Selected <- (y == c1 | y== c2)
xx = as.matrix(x[Selected,])
yy = factor(y[Selected])
table(yy)

xx.2 = x[Selected,]
xx.2$label = yy
write.csv(xx.2,paste(case,'_',c1,'_',c2,'.csv',sep=''),row.names=F)

### 6-CV
alpha=1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nfolds,alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_",nfolds,"F-CV_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

### LOO-CV
nfolds = "LOO"
alpha=1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nrow(xx),alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_LOO_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

fit=glmnet(xx,yy,family="bin",lambda=cv$lambda.min,alpha=alpha)
b=data.matrix(fit$beta[as.numeric(fit$beta)!=0,])
b=b[order(abs(b[,1]),decreasing=T),]
b
write.table(b,paste("_CoefOnNuclear_",c1,"_",c2,"_LOO_.txt",sep=""),sep="\t",col.names=NA)

##########################################################################################################
########################################   Normal vs UDH
##########################################################################################################
c1 = "Normal"
c2 = "UDH"
nfolds = SelectedFolds

Selected <- (y == c1 | y== c2)
xx = as.matrix(x[Selected,])
yy = factor(y[Selected])
table(yy)

xx.2 = x[Selected,]
xx.2$label = yy
write.csv(xx.2,paste(case,'_',c1,'_',c2,'.csv',sep=''),row.names=F)

### 6-CV
alpha = 1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nfolds,alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_",nfolds,"F-CV_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

### LOO-CV
nfolds = "LOO"
alpha = 1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nrow(xx),alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_LOO_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

fit=glmnet(xx,yy,family="bin",lambda=cv$lambda.min,alpha=alpha)
b=data.matrix(fit$beta[as.numeric(fit$beta)!=0,])
b=b[order(abs(b[,1]),decreasing=T),]
b
write.table(b,paste("_CoefOnNuclear_",c1,"_",c2,"_LOO_.txt",sep=""),sep="\t",col.names=NA)

##########################################################################################################
########################################   Normal vs ADH
##########################################################################################################
c1 = "Normal"
c2 = "ADH"
nfolds = SelectedFolds

Selected <- (y == c1 | y== c2)
xx = as.matrix(x[Selected,])
yy = factor(y[Selected])
table(yy)

xx.2 = x[Selected,]
xx.2$label = yy
write.csv(xx.2,paste(case,'_',c1,'_',c2,'.csv',sep=''),row.names=F)

### 6-CV
alpha = 1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nfolds,alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_",nfolds,"F-CV_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

### LOO-CV
nfolds = "LOO"
set.seed(12345)
alpha = 1
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nrow(xx),alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_LOO_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

fit=glmnet(xx,yy,family="bin",lambda=cv$lambda.min,alpha=alpha)
b=data.matrix(fit$beta[as.numeric(fit$beta)!=0,])
b=b[order(abs(b[,1]),decreasing=T),]
b
write.table(b,paste("_CoefOnNuclear_",c1,"_",c2,"_LOO_.txt",sep=""),sep="\t",col.names=NA)

##########################################################################################################
########################################   UDH vs DCIS
##########################################################################################################
c1 = "UDH"
c2 = "DCIS"
nfolds = SelectedFolds

Selected <- (y == c1 | y== c2)
xx = as.matrix(x[Selected,])
yy = factor(y[Selected])
table(yy)

xx.2 = x[Selected,]
xx.2$label = yy
write.csv(xx.2,paste(case,'_',c1,'_',c2,'.csv',sep=''),row.names=F)

### 6-CV
alpha = 1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nfolds,alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_",nfolds,"F-CV_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

### Loo-CV
nfolds = "LOO"
alpha=1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nrow(xx),alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_LOO_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

fit=glmnet(xx,yy,family="bin",lambda=cv$lambda.min,alpha=alpha)
b=data.matrix(fit$beta[as.numeric(fit$beta)!=0,])
b=b[order(abs(b[,1]),decreasing=T),]
b
write.table(b,paste("_CoefOnNuclear_",c1,"_",c2,"_LOO_.txt",sep=""),sep="\t",col.names=NA)


##########################################################################################################
########################################   ADH vs DCIS
##########################################################################################################
c1 = "ADH"
c2 = "DCIS"
nfolds = SelectedFolds

Selected <- (y == c1 | y== c2)
xx = as.matrix(x[Selected,])
yy = factor(y[Selected])
table(yy)

xx.2 = x[Selected,]
xx.2$label = yy
write.csv(xx.2,paste(case,'_',c1,'_',c2,'.csv',sep=''),row.names=F)

### 6-CV
alpha=1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nfolds,alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_",nfolds,"F-CV_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

### LOO-CV
nfolds = "LOO"
alpha=1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nrow(xx),alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_LOO_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

fit=glmnet(xx,yy,family="bin",lambda=cv$lambda.min,alpha=alpha)
b=data.matrix(fit$beta[as.numeric(fit$beta)!=0,])
b=b[order(abs(b[,1]),decreasing=T),]
b
write.table(b,paste("_CoefOnNuclear_",c1,"_",c2,"_LOO_.txt",sep=""),sep="\t",col.names=NA)

##########################################################################################################
########################################   Normal vs DCIS
##########################################################################################################
# Selecting 2 Summary Features for this category of comparison
f= x2
x = f[,-ncol(f)]
y = f[,ncol(f)]

c1 = "Normal"
c2 = "DCIS"
nfolds = SelectedFolds

Selected <- (y == c1 | y== c2)
xx = as.matrix(x[Selected,])
yy = factor(y[Selected])
table(yy)

xx.2 = x[Selected,]
xx.2$label = yy
write.csv(xx.2,paste(case,'_',c1,'_',c2,'.csv',sep=''),row.names=F)

### 6-CV
alpha=1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nfolds,alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_",nfolds,"F-CV_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

### LOO-CV
nfolds = "LOO"
alpha=1
set.seed(12345)
cv=cv.glmnet(xx,yy,family="bin",type="class",nfolds=nrow(xx),alpha=alpha,keep=T,grouped=T)
cv.preds=cv$fit.preval[,cv$lambda==cv$lambda.min]
if(length(cv.preds) == length(yy)){
  pred <- prediction(cv.preds, yy)
}else{
  pred <- prediction(cv.preds[,1], yy)  
}
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
auc=round(performance(pred,"auc")@y.values[[1]],2)
auc
tiff(paste("_AUC_",c1,"-",c2,"_LOO_alpha-",alpha,".tiff",sep=""),width=1000,height=500)
par(mfrow=c(1,2))
plot(cv) 
plot(perf,col="red",main=paste(c1," vs ",c2," ",nfolds,"-Fold CV, Alpha=",alpha," AUC=",auc,sep=""))
text(x = 0.5,y=0.5,paste("AUC = ",auc,sep=""))
dev.off()

fit=glmnet(xx,yy,family="bin",lambda=cv$lambda.min,alpha=alpha)
b=data.matrix(fit$beta[as.numeric(fit$beta)!=0,])
b=b[order(abs(b[,1]),decreasing=T),]
b
write.table(b,paste("_CoefOnNuclear_",c1,"_",c2,"_LOO_.txt",sep=""),sep="\t",col.names=NA)

