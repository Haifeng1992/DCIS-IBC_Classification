############### initialization #####################
setwd("/Users/xuhaifeng/Documents/PhD_project_1/scripts")
source("my_function.r")
initialization()

############### import data #####################
setwd("/Users/xuhaifeng/Documents/PhD_project_1/data/validation_data")
expr_data = readRDS("train_data_10_genes.rds")
vali_data = readRDS("vali_data_10_genes_new.rds")
vali_all = readRDS("validation_processed_new.rds")


setwd("/Users/xuhaifeng/Documents/PhD_project_1/data")

rando = sample.int(100,1)
set.seed(rando)

expr_data[,2:ncol(expr_data)] = apply(expr_data[,2:ncol(expr_data)], 2, my_zsore) 
vali_data[,2:ncol(vali_data)] = apply(vali_data[,2:ncol(vali_data)], 2, as.character)
vali_data[,2:ncol(vali_data)] = apply(vali_data[,2:ncol(vali_data)], 2, as.numeric)
vali_data[,2:ncol(vali_data)] = apply(vali_data[,2:ncol(vali_data)], 2, my_zsore)

expr_data$label = as.character(expr_data$label)
vali_data$label = as.character(vali_data$label)

vali_micro = vali_data[vali_data$label == "MI-DCIS",]
vali_rest = vali_data[vali_data$label != "MI-DCIS",]
vali_rest = vali_rest[vali_rest$label != "MI-DCIS (PIK)",]

test_data <- vali_rest
train_data <- expr_data

x.test = test_data[, 2:ncol(test_data)]
y.test = test_data[, 1]
x.train = train_data[, 2:ncol(train_data)]
y.train = train_data[, 1]

y.train = as.character(y.train)
y.train[y.train == "IBC"] = 0
y.train[y.train == "DCIS"] = 1
y.train = as.numeric(y.train)
y.train = factor(y.train, levels = c(0,1), ordered = TRUE)

y.test = as.character(y.test)
y.test[y.test == "IBC"] = 0
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
plot(cv_output)

model= glmnet(x = as.matrix(x.train), y = as.matrix(as.factor(y.train)), 
              alpha = 1, family = "binomial", type.measure="auc")

betaGene = model$beta[,model$lambda==cv_output$lambda.min]
betaGene = betaGene[betaGene != 0]
betaGene = as.data.frame(betaGene)
GeneList = rownames(betaGene)

#predidct on training
result = predict(model, as.matrix(x.train), type = "response", s=cv_output$lambda.min)
result = as.numeric(result)
y.train = as.numeric(as.character(y.train))
k = roc(y.train, result)
par(mfrow=c(1,1))
plot(k, col = "black",  yaxt = "n",cex.axis=1.25, cex.lab = 1.25, lwd = 2)
auc(k)


#predict on validation
options(scipen=999)
result = predict(model, as.matrix(x.test), type = "response", s=cv_output$lambda.min)



y.test = as.numeric(as.character(y.test))
result = as.numeric(result)


k = roc(y.test, result, direction = "<")
plot(k, col = "red", add = TRUE, lwd = 2)
#title("result of validation set")
auc(k)

# This line uses the default 0.5 as the cutoff but it turns to no optimal results
#label_result = predict(model, as.matrix(x.test),type = "class",  s=cv_output$lambda.min)

# Use smaller cutoff, reported in the paper. 
#This is likely due to the differences between 
#the microarray gene expression data (training cohort) 
#and RNA-sequencing data (validation cohort)

label_result = result
the_cutoff = 0.03
label_result[label_result>the_cutoff] = 1
label_result[label_result<=the_cutoff] = 0
confusionMatrix(as.factor(as.character(label_result)), as.factor(as.character(y.test)))
