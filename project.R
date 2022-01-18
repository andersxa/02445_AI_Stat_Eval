library(nnet) #For classification
library(cvTools) #For CV
library(beeswarm)
library(lawstat)
library(caret)
set.seed(1234) # for reproducibility
load("armdata.RData")
#PROCESS DATA
{
  armdata.array = array(armdata,dim=c(16,10,10,100,3))
  summed.data = data.frame()
  distances = c(15,22.5,30,37.5,45.0,"N")
  obstacles = c("S","M","T")
  for(ex in 1:16){
    for(per in 1:10){
      for(ite in 1:10){
        x = armdata[[ex]][[per]][[ite]][,1]
        y = armdata[[ex]][[per]][[ite]][,2]
        z = armdata[[ex]][[per]][[ite]][,3]
        distance = distances[ceiling(ex / 3)]
        obstacle = obstacles[(ex-1) %% 3+1]
        
        center.mass = mean(sqrt((x-x[1])^2+(y-y[1])^2+(z-z[1])^2),na.rm=TRUE)
        
        Lx = lm(x ~ polym(1:length(x),degree=3,raw=TRUE))
        Ly = lm(y ~ polym(1:length(y),degree=3,raw=TRUE))
        Lz = lm(z ~ polym(1:length(z),degree=3,raw=TRUE))
        
        sum.row = data.frame(
          row.names=c(paste0(ex,".",per,".",ite)),
          experiment=ex,
          person=per,
          repetition=ite,
          distance=distance,
          obstacle=obstacle,
          center.mass=center.mass,
          px.I = coef(Lx)[1],
          px.1 = coef(Lx)[2],
          px.2 = coef(Lx)[3],
          px.3 = coef(Lx)[4],
          py.I = coef(Ly)[1],
          py.1 = coef(Ly)[2],
          py.2 = coef(Ly)[3],
          py.3 = coef(Ly)[4],
          pz.I = coef(Lz)[1],
          pz.1 = coef(Lz)[2],
          pz.2 = coef(Lz)[3],
          pz.3 = coef(Lz)[4]
        )
        summed.data = rbind(summed.data,sum.row)
      }
    }
  }
  
  summed.data$obstacle[summed.data$experiment == 16] = "N"
  
  summed.data$experiment = factor(summed.data$experiment,ordered=FALSE)
  summed.data$person = factor(summed.data$person,ordered=TRUE)
  summed.data$repetition = factor(summed.data$repetition,ordered=FALSE)
  summed.data$distance = factor(summed.data$distance,ordered=TRUE,levels=c("N",15,22.5,30,37.5,45.0))
  summed.data$obstacle = factor(summed.data$obstacle,ordered=TRUE,levels=c("N",obstacles))
  
  summed.data.preNA = summed.data
  summed.data=na.omit(summed.data)
}
#PLOT
{
  pdf("OBSTACLE_DIST.pdf",12,6)
  par(mfrow = c(2,3))
  par(mar = c(par("mar")[1], par("mar")[2], 0.1, 0.1))
  boxplot(mean.x ~ distance, summed.data, xlab="Distance", ylab="Mean x", col=2)
  boxplot(mean.y ~ distance, summed.data, xlab="Distance", ylab="Mean y", col=3)
  boxplot(mean.z ~ distance, summed.data, xlab="Distance", ylab="Mean z", col=4)
  boxplot(mean.x ~ obstacle, summed.data, xlab="Height", ylab="Mean x", col=2)
  boxplot(mean.y ~ obstacle, summed.data, xlab="Height", ylab="Mean y", col=3)
  boxplot(mean.z ~ obstacle, summed.data, xlab="Height", ylab="Mean z", col=4)
  dev.off()
  
  pdf("OBSTACLE_DIST_CM.pdf",12,4)
  par(mfrow = c(1,3))
  par(mar = c(par("mar")[1], par("mar")[2], 0.1, 0.1))
  boxplot(center.mass ~ distance, summed.data, xlab="Distance", ylab="Curve Center of Mass", col=2)
  boxplot(center.mass ~ obstacle, summed.data, xlab="Height", ylab="Curve Center of Mass", col=3)
  boxplot(center.mass ~ person, summed.data, xlab="Person", ylab="Curve Center of Mass", col=4)
  dev.off()
  
  pdf("OBSTACLE_DIST_BEESWARM.pdf",12,6)
  par(mfrow = c(2,3))
  par(mar = c(par("mar")[1], par("mar")[2], 0.1, 0.1))
  beeswarm(mean.x ~ distance, summed.data, xlab="Distance", ylab="Mean x", col=2,corral = "random", priority="density")
  beeswarm(mean.y ~ distance, summed.data, xlab="Distance", ylab="Mean y", col=3,corral = "random", priority="density")
  beeswarm(mean.z ~ distance, summed.data, xlab="Distance", ylab="Mean z", col=4,corral = "random", priority="density")
  beeswarm(mean.x ~ obstacle, summed.data, xlab="Height", ylab="Mean x", col=2,corral = "random", priority="density")
  beeswarm(mean.y ~ obstacle, summed.data, xlab="Height", ylab="Mean y", col=3,corral = "random", priority="density")
  beeswarm(mean.z ~ obstacle, summed.data, xlab="Height", ylab="Mean z", col=4,corral = "random", priority="density")
  dev.off()
  
  colors = c(
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf")
  
  pdf("INTERACTION_OBST_DIST.pdf",12,6)
  par(mfrow = c(1,2))
  par(mar = c(par("mar")[1], par("mar")[2], 0.1, 0.1))
  Person = summed.data$person
  interaction.plot(summed.data$distance,Person,summed.data$center.mass, xlab="Distance", ylab="Mean center of mass", col=colors, leg.bty="o", fixed=TRUE,lwd=2)
  interaction.plot(summed.data$obstacle,Person,summed.data$center.mass, xlab="Height", ylab="Mean center of mass", col=colors, leg.bty="o", fixed=TRUE,lwd=2)
  dev.off()
  
  pdf("INTERACTION_EXPERIMENT.pdf",6,6)
  par(mfrow = c(1,1))
  par(mar = c(par("mar")[1], par("mar")[2], 0.1, 0.1))
  Person = summed.data$person
  interaction.plot(summed.data$experiment,Person,summed.data$center.mass, xlab="Experiment", ylab="Mean center of mass", col=colors, leg.bty="o", fixed=TRUE,lwd=2)
  dev.off()
}

#OVERALL NORMALITY
pdf("QQ_PLOT_CENTER_MASS_NORMALITY_HOMOGENUITY.pdf",12,6)
par(mfrow=c(1,2))
qqnorm(summed.data$center.mass);qqline(summed.data$center.mass)
shapiro.test(summed.data$center.mass)
#TEST FOR HOMOGENUITY - I.E. EQUAL VARIANCE IN GROUPS
group.tests = lapply(split(summed.data, summed.data$person), function(x) shapiro.test(x$center.mass))
group.means.p.values = array(unlist(lapply(group.tests, function(x) x$p.value)))
plot(sort(group.means.p.values),xlab="Index of person sorted by p-value",ylab="p value",main="Group Normality", pch=order(group.means.p.values));abline(h=0.05,col=2,lwd=2)
legend("topleft", legend=1:10, pch=1:10)
dev.off()
print(group.means.p.values)
kruskal.test(center.mass ~ person,data=summed.data)


#####
#PERSON EFFECT ON EXPERIMENT TESTS
#Test for equal variance in groups
levene.test(summed.data$center.mass,interaction(summed.data$person,summed.data$experiment))
#Test for person interaktion med experiment
L = lm(center.mass ~ person*experiment, data=summed.data)
L2 = lm(center.mass ~ experiment, data=summed.data)
#PLOT
pdf("EXPERIMENT_REGRESSION.pdf",12,6)
par(mfrow=c(1,2));plot(L,1);plot(L,2)
dev.off()
#TEST FOR PER GROUP NORMALITY / HOMOSCEDASTICITY
shapiro.test(residuals(L))

#Test of main effect
kruskal.test(center.mass ~ experiment,data=summed.data)

#####
#PERSON EFFECT ON OBSTACLE TESTS
#Test for equal variance in groups
levene.test(summed.data$center.mass,interaction(summed.data$person,summed.data$obstacle))
#Test for person interaktion med experiment
L = lm(center.mass ~ person*obstacle, data=summed.data)
#PLOT
pdf("OBSTACLE_REGRESSION.pdf",12,6)
par(mfrow=c(1,2));plot(L,1);plot(L,2)
dev.off()
#Test for interaction
kruskal.test(center.mass ~ interaction(person,obstacle),data=summed.data)
#Test for main effect
kruskal.test(center.mass ~ obstacle,data=summed.data)

#####
#PERSON EFFECT ON DISTANCE TESTS
#Test for equal variance in groups
levene.test(summed.data$center.mass,interaction(summed.data$person,summed.data$distance))
#Test for person interaktion med experiment
L = lm(center.mass ~ person*distance, data=summed.data)
#PLOT
pdf("DISTANCE_REGRESSION.pdf",12,6)
par(mfrow=c(1,2));plot(L,1);plot(L,2)
dev.off()
#Test for interaction
kruskal.test(center.mass ~ interaction(person,distance),data=summed.data)
#Test for mai neffect
kruskal.test(center.mass ~ distance,data=summed.data)

#Classification distance ALL POLYNOMIAL COEFFICIENTS
cv.classify.poly = function(obstacle){
  px.I <- summed.data$px.I
  py.I <- summed.data$py.I
  pz.I <- summed.data$pz.I
  px.1 <- summed.data$px.1
  py.1 <- summed.data$py.1
  pz.1 <- summed.data$pz.1
  px.2 <- summed.data$px.2
  py.2 <- summed.data$py.2
  pz.2 <- summed.data$pz.2
  px.3 <- summed.data$px.3
  py.3 <- summed.data$py.3
  pz.3 <- summed.data$pz.3
  
  N = length(px.I)
  
  K = 10;
  set.seed(1234) # for reproducibility
  random.order = array(dim=c(16*10*10))
  for(i in 1:(16*10)){
    random.order[((i-1)*10+1):(i*10)] = sample(10)
  }
  
  # Create k-fold crossvalidation partition
  # Create K-fold split
  CV <- list(n=N,K=K,R=1,subsets=1:N,which=random.order[complete.cases(summed.data.preNA)]) #cvFolds(length(obstacle), K=K)
  max.length = 0
  for(k in 1:K){
    max.length = max(length(CV$subsets[CV$which==k]),max.length)
  }
  
  #Array to store predictions
  test.true = matrix(rep(NA, times=K*max.length), nrow=K)
  xyz.test.pred = matrix(rep(NA, times=K*max.length), nrow=K)
  xy.test.pred = matrix(rep(NA, times=K*max.length), nrow=K)
  yz.test.pred = matrix(rep(NA, times=K*max.length), nrow=K)
  xz.test.pred = matrix(rep(NA, times=K*max.length), nrow=K)
  x.test.pred = matrix(rep(NA, times=K*max.length), nrow=K)
  y.test.pred = matrix(rep(NA, times=K*max.length), nrow=K)
  z.test.pred = matrix(rep(NA, times=K*max.length), nrow=K)
  
  for(k in 1:K){
    print(paste('Crossvalidation fold ', k, '/', K, sep=''))
    not.fold = CV$which!=k
    #Extract train data 
    px.I_train <- px.I[not.fold];
    px.1_train <- px.1[not.fold];
    px.2_train <- px.2[not.fold];
    px.3_train <- px.3[not.fold];
    py.I_train <- py.I[not.fold];
    py.1_train <- py.1[not.fold];
    py.2_train <- py.2[not.fold];
    py.3_train <- py.3[not.fold];
    pz.I_train <- pz.I[not.fold];
    pz.1_train <- pz.1[not.fold];
    pz.2_train <- pz.2[not.fold];
    pz.3_train <- pz.3[not.fold];
    obstacle_train <- obstacle[not.fold];
    data_frame_train <- data.frame(px.I=px.I_train,
                                   px.1=px.1_train,
                                   px.2=px.2_train,
                                   px.3=px.3_train,
                                   py.I=py.I_train,
                                   py.1=py.1_train,
                                   py.2=py.2_train,
                                   py.3=py.3_train,
                                   pz.I=pz.I_train,
                                   pz.1=pz.1_train,
                                   pz.2=pz.2_train,
                                   pz.3=pz.3_train,
                                   obstacle=obstacle_train)
    
    fold = CV$which==k
    #Extract test data 
    px.I_test <- px.I[fold];
    px.1_test <- px.1[fold];
    px.2_test <- px.2[fold];
    px.3_test <- px.3[fold];
    py.I_test <- py.I[fold];
    py.1_test <- py.1[fold];
    py.2_test <- py.2[fold];
    py.3_test <- py.3[fold];
    pz.I_test <- pz.I[fold];
    pz.1_test <- pz.1[fold];
    pz.2_test <- pz.2[fold];
    pz.3_test <- pz.3[fold];
    obstacle_test <- obstacle[fold];
    test_df <- data.frame(px.I=px.I_test,
                          px.1=px.1_test,
                          px.2=px.2_test,
                          px.3=px.3_test,
                          py.I=py.I_test,
                          py.1=py.1_test,
                          py.2=py.2_test,
                          py.3=py.3_test,
                          pz.I=pz.I_test,
                          pz.1=pz.1_test,
                          pz.2=pz.2_test,
                          pz.3=pz.3_test,
                          obstacle=obstacle_test)
    
    xyz.multinom.model <- multinom(obstacle ~ px.I + px.1 + px.2 + px.3 + py.I + py.1 + py.2 + py.3 + pz.I + pz.1 + pz.2 + pz.3, data = data_frame_train, trace=FALSE)
    xy.multinom.model <- multinom(obstacle ~ px.I + px.1 + px.2 + px.3 + py.I + py.1 + py.2 + py.3, data = data_frame_train, trace=FALSE)
    yz.multinom.model <- multinom(obstacle ~ py.I + py.1 + py.2 + py.3 + pz.I + pz.1 + pz.2 + pz.3, data = data_frame_train, trace=FALSE)
    xz.multinom.model <- multinom(obstacle ~ px.I + px.1 + px.2 + px.3 + pz.I + pz.1 + pz.2 + pz.3, data = data_frame_train, trace=FALSE)
    x.multinom.model <- multinom(obstacle ~ px.I + px.1 + px.2 + px.3, data = data_frame_train, trace=FALSE)
    y.multinom.model <- multinom(obstacle ~ py.I + py.1 + py.2 + py.3, data = data_frame_train, trace=FALSE)
    z.multinom.model <- multinom(obstacle ~ pz.I + pz.1 + pz.2 + pz.3, data = data_frame_train, trace=FALSE)
    
    xyz.pred = predict(xyz.multinom.model,test_df,"class")
    xy.pred = predict(xy.multinom.model,test_df,"class")
    yz.pred = predict(yz.multinom.model,test_df,"class")
    xz.pred = predict(xz.multinom.model,test_df,"class")
    x.pred = predict(x.multinom.model,test_df,"class")
    y.pred = predict(y.multinom.model,test_df,"class")
    z.pred = predict(z.multinom.model,test_df,"class")
    
    xyz.test.pred[k,1:length(xyz.pred)] = xyz.pred
    xy.test.pred[k,1:length(xy.pred)] = xy.pred
    yz.test.pred[k,1:length(yz.pred)] = yz.pred
    xz.test.pred[k,1:length(xz.pred)] = xz.pred
    x.test.pred[k,1:length(x.pred)] = x.pred
    y.test.pred[k,1:length(y.pred)] = y.pred
    z.test.pred[k,1:length(z.pred)] = z.pred
    
    test.true[k,1:length(obstacle_test)] = obstacle_test
  }
  
  #Array to store predictions
  xyz.pred.flat = xyz.test.pred;dim(xyz.pred.flat) = 16*10*10;xyz.pred.flat=factor(xyz.pred.flat)
  xy.pred.flat = xy.test.pred;dim(xy.pred.flat) = 16*10*10;xy.pred.flat=factor(xy.pred.flat)
  yz.pred.flat = yz.test.pred;dim(yz.pred.flat) = 16*10*10;yz.pred.flat=factor(yz.pred.flat)
  xz.pred.flat = xz.test.pred;dim(xz.pred.flat) = 16*10*10;xz.pred.flat=factor(xz.pred.flat)
  x.pred.flat = x.test.pred;dim(x.pred.flat) = 16*10*10;x.pred.flat=factor(x.pred.flat)
  y.pred.flat = y.test.pred;dim(y.pred.flat) = 16*10*10;y.pred.flat=factor(y.pred.flat)
  z.pred.flat = z.test.pred;dim(z.pred.flat) = 16*10*10;z.pred.flat=factor(z.pred.flat)
  
  test.true.flat = test.true;dim(test.true.flat) = 16*10*10;test.true.flat=factor(test.true.flat)
  
  #CM = confusionMatrix(xyz.pred.flat,test.true.flat)
  
  classifier.mcnemar.test = function(pred1,pred2,true){
    N11 = sum((pred1 == true) & (pred2 == true),na.rm=TRUE)
    N21 = sum((pred1 != true) & (pred2 == true),na.rm=TRUE)
    N12 = sum((pred1 == true) & (pred2 != true),na.rm=TRUE)
    N22 = sum((pred1 != true) & (pred2 != true),na.rm=TRUE)
    cont.mat = matrix(c(N11,N12,N21,N22),nrow = 2,dimnames = list("Classifier 1" = c("True", "False"), "Classifier 2" = c("True", "False")))
    print(cont.mat)
    ci.p <- prop.test(max(cont.mat[1,2],cont.mat[2,1]),cont.mat[2,1]+cont.mat[1,2])$conf
    print("mcnemar test conf:")
    print(ci.p)
    return(list(test=mcnemar.test(cont.mat),conf.lb=ci.p[1],conf.ub=ci.p[2]))
  }
  
  classifiers = list(xyz.test.pred,xy.test.pred,yz.test.pred,xz.test.pred,x.test.pred,y.test.pred,z.test.pred)
  classifier.names = c("xyz","xy","yz","xz","x","y","z")
  C = length(classifiers)
  d.class.p.values = matrix(rep(NA, times=C*C), nrow=C, dimnames=list(classifier.names,classifier.names))
  test.conf.lbs = matrix(rep(NA, times=C*C), nrow=C, dimnames=list(classifier.names,classifier.names))
  test.conf.ubs = matrix(rep(NA, times=C*C), nrow=C, dimnames=list(classifier.names,classifier.names))
  for(i in 1:C){
    for(j in 1:C){
      if(i>=j) next
      print(paste("CL1:",classifier.names[i],"vs.","CL2:",classifier.names[j]))
      test.out = classifier.mcnemar.test(classifiers[[i]],classifiers[[j]],test.true)
      d.class.p.values[i,j] = test.out$test$p.value
      test.conf.lbs[i,j] = test.out$conf.lb
      test.conf.ubs[i,j] = test.out$conf.ub
      
    }
  }
  print("McNemar's test proportion lower bound")
  print(test.conf.lbs)
  print("McNemar's test proportion upper bound")
  print(test.conf.ubs)
  d.class.p.values[!is.na(d.class.p.values)] = p.adjust(d.class.p.values[!is.na(d.class.p.values)],method="holm")
  print(d.class.p.values < 0.05)
  
  for(i in 1:C){
    print(paste("Classifier:", classifier.names[i], "Accuracy:", mean(classifiers[[i]] == test.true,na.rm=TRUE)))
  }
  return(d.class.p.values)
}
cv.classify.poly(summed.data$distance)
cv.classify.poly(summed.data$obstacle)

#full.scatter
xyzdata = (data.frame(center.mass=summed.data$center.mass,px.I=summed.data$px.I,px.1=summed.data$px.1,px.2=summed.data$px.2,px.3=summed.data$px.3,py.I=summed.data$py.I,py.1=summed.data$py.1,py.2=summed.data$py.2,py.3=summed.data$py.3,pz.I=summed.data$pz.I,pz.1=summed.data$pz.1,pz.2=summed.data$pz.2,pz.3=summed.data$pz.3))
cor(xyzdata)


