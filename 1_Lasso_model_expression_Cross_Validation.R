############### initialization #####################
setwd("/Users/xuhaifeng/Documents/PhD_Project_1/scripts")
source("my_function.r")
initialization()

### import the data ##########
setwd("/Users/xuhaifeng/Documents/PhD_Project_1/data")
expr_data = readRDS("expr_data_clean.rds")

# Run this line if you want to randomlize the labels 
#expr_data[,1] = sample(expr_data[,1])
expr_data = expr_data[order(rownames(expr_data)), ]

### folds creation ########
# divide the data to 10 folds and check if they@re duplicated
#test_index <- createDataPartition(randomized_data$label, p = 0.1, times = 10)

#Folds with subtypes creation
LumB = expr_data[expr_data$subtype == "LumB",]
LumA = expr_data[expr_data$subtype == "LumA",]
Her2 = expr_data[expr_data$subtype == "Her2",]
Basal = expr_data[expr_data$subtype == "Basal",]
Normal = expr_data[expr_data$subtype == "Normal",]

test_index_LumB = createFolds(y = LumB$label, k = 10, list = TRUE, returnTrain = FALSE)
test_index_LumA = createFolds(y = LumA$label, k = 10, list = TRUE, returnTrain = FALSE)
test_index_Her2 = createFolds(y = Her2$label, k = 10, list = TRUE, returnTrain = FALSE)
test_index_Basal = createFolds(y = Basal$label, k = 10, list = TRUE, returnTrain = FALSE)
test_index_Normal = createFolds(y = Normal$label, k = 10, list = TRUE, returnTrain = FALSE)

# initialization of the nested cross-validation
results = matrix(0,5,10)
mse_results = matrix(0,1,10)
GeneList = list()
i = 1
par(mfrow=c(3,4))
tmp = matrix(as.factor(0), 1, 3)
my_pal <- seecol(pal_unikn_light)

#rando = sample.int(100,1)
rando = 100

#### nested cross-validation ####
for (i in 1:10) {
  
  set.seed(rando)
  test_data <- rbind(LumB[as.numeric(test_index_LumB[[i]]), ], 
                     LumA[as.numeric(test_index_LumA[[i]]), ],
                     Her2[as.numeric(test_index_Her2[[i]]), ],
                     Basal[as.numeric(test_index_Basal[[i]]), ],
                     Normal[as.numeric(test_index_Normal[[i]]), ]
  )
  train_data <- rbind(LumB[as.numeric(-test_index_LumB[[i]]), ], 
                      LumA[as.numeric(-test_index_LumA[[i]]), ],
                      Her2[as.numeric(-test_index_Her2[[i]]), ],
                      Basal[as.numeric(-test_index_Basal[[i]]), ],
                      Normal[as.numeric(-test_index_Normal[[i]]), ]
  )
  
  #remove the labels for x.train and x.test using column_3-column_N
  x.test = test_data[, 3:ncol(expr_data)]
  y.test = test_data[, 1]
  x.train = train_data[, 3:ncol(expr_data)]
  y.train = train_data[, 1]
  
  y.train = as.character(y.train)
  y.train[y.train == "IDC"] = 0
  y.train[y.train == "DCIS"] = 1
  y.train = as.numeric(y.train)
  y.train = factor(y.train, levels = c(0,1), ordered = TRUE)
  
  y.test = as.character(y.test)
  y.test[y.test == "IDC"] = 0
  y.test[y.test == "DCIS"] = 1
  y.test = as.numeric(y.test)
  
  
  pr_label = vector(mode = "logical", length = length(y.test))
  for (pr_i in 1:length(y.test)) {
    if(y.test[pr_i] == 1){
      pr_label[pr_i] = TRUE
    }
  }
  
  y.test = factor(y.test, levels = c(0,1), ordered = TRUE)
  cv_output <- cv.glmnet(x = as.matrix(x.train), y = as.factor(y.train), 
                         alpha = 1, family = "binomial", type.measure="auc")
  #plot(cv_output)
  model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
                alpha = 1, family = "binomial", type.measure="auc")
  
  #get the common genes
  betaGene = model$beta[,model$lambda==cv_output$lambda.min]
  betaGene = betaGene[betaGene != 0]
  betaGene = as.data.frame(betaGene)
  GeneList[[i]] = rownames(betaGene)
  
  
  result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)
  label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)
  result = as.numeric(result)
  
  #calculate average mse for this fold
  mse = cbind(as.data.frame(result), as.data.frame(as.numeric(as.character(y.test))))
  temp_c = as.data.frame(matrix(0,nrow(mse),1))
  mse = cbind(mse, temp_c)
  for (j in 1:nrow(mse)) {
    mse[j,3] = (mse[j,1] - mse[j,2])^2
  }
  mse_results[1,i] = mean(mse[j,3])
  
  roc = roc(response = y.test, predictor = result, quiet = TRUE, direction = "<")
  plot(roc)
  auc = auc(roc)
  
  #PR curve
  wpr<-pr.curve(result, weights.class0 = pr_label, curve = TRUE,
                max.compute = T, min.compute = T, rand.compute = T)
  if (i == 1){
    plot(wpr, add = FALSE, color = my_pal[i], yaxt = "n",cex.axis=1.5, cex.lab = 1.5, main = "")
  }else{
    plot(wpr, add = FALSE, color = my_pal[i], yaxt = "n",cex.axis=1.5, cex.lab = 1.5, main = "")
  }
  
  au_prc = wpr$auc.integral
  no_skill = length(y.test[y.test==1]) / length(y.test)
  lines(x= seq(from=0, to=1, by=0.05), y = rep_len(no_skill, 21),  type = "c", col = "red", lwd = 3)
  text(0.5, 0.07, round(no_skill,4))
  axis(side = 2, labels = FALSE)
  
  ## Draw the y-axis.
  axis(side = 2,
       ## Rotate the labels.
       las = 2,
       ## Adjust the label position.
       mgp = c(3, 0.75, 0), cex.axis=1.5, cex.lab = 1.5)
  
  #balanced accuracy at cutoff 0.5
  CM = confusionMatrix(as.factor(label_result), as.factor(as.character(y.test)))
  
  results[1,i] = auc
  results[2,i] = CM$byClass[11]
  results[3,i] = au_prc
  results[4,i] = CM$byClass[1]
  results[5,i] = CM$byClass[2]
}
#tmp = tmp[2:nrow(tmp),]

rownames(results) = c("AUC", "balanced_accuracy","auprc","spec", "sens(DCIS)")

setwd("/Users/xuhaifeng/Documents/phd1/results/basic")
#saveRDS(results, "lasso_expr_new.rds")

final_auc = mean(results[1,])
final_auc
sd_auc = sd(results[1,])
sd_auc

final_ba = mean(results[2,]) 
final_ba
sd_ba = sd(results[2,])
sd_ba

final_prc = mean(results[3,]) 
final_prc
sd_prc = sd(results[3,])
sd_prc

final_mse = mean(mse_results)
final_mse
sd_mse = sd(mse_results)
sd_mse

# calculate the feature number
num = 0
nums = matrix(0,1,10)
for (i in 1:10) {
  nums[i] = length(GeneList[[i]])
}
num = mean(nums)
num
sd_num = sd(nums)
sd_num

#find the common gene
common_gene = GeneList[[1]]
for (j in 2:10) {
  common_gene = intersect(common_gene, GeneList[[i]])
}
common_gene
