library(Boruta)
library(glmnet)
# Specify the file paths of the CSV files
geneloc <- "C:/Users/matth/anaconda3/envs/data1/data1/lihc/data/gene.csv"
methylloc <- "C:/Users/matth/anaconda3/envs/data1/data1/lihc/data/methyl.csv"
miRNAloc <- "C:/Users/matth/anaconda3/envs/data1/data1/lihc/data/miRNA.csv"
labelloc <- "C:/Users/matth/anaconda3/envs/data1/data1/lihc/label.csv"

# Read the CSV files into separate data frames
Gene <- read.csv(geneloc)
Methyl <- read.csv(methylloc)
miRNA <- read.csv(miRNAloc)
label <- read.csv(labelloc)


#data cleaning Gene dataset
Gene <- data.frame(t(Gene))
colnames(Gene) <-Gene[1,]
Gene <- Gene[-1,]
Gene <- cbind(label$label1,Gene)
colnames(Gene)[1]<- "target"
Gene <- sapply(Gene,as.numeric)
Gene <- data.frame(Gene)


#data cleaning Methyl dataset
Methyl <- data.frame(t(Methyl))
colnames(Methyl) <-Methyl[1,]
Methyl <- Methyl[-1,]
Methyl <- cbind(label$label1,Methyl)
colnames(Methyl)[1]<- "target"
Methyl <- sapply(Methyl,as.numeric)
Methyl <- data.frame(Methyl)


#data cleaning miRNA dataset
miRNA <- data.frame(t(miRNA))
colnames(miRNA) <-miRNA[1,]
miRNA <- miRNA[-1,]
miRNA <- cbind(label$label1,miRNA)
colnames(miRNA)[1]<- "target"
miRNA <- sapply(miRNA, as.numeric)
miRNA <- data.frame(miRNA)


###  feature selection 
## Gene
Boruta_Gene <- Boruta(Gene[,-1],Gene$target)
Gene_Boruta_features <- getSelectedAttributes(Boruta_Gene)
write.csv(Gene_Boruta_features, "C:/Users/matth/anaconda3/envs/data1/data1/Lasso_Boruta_Selection/lihc/Boruta_gene.csv", row.names = FALSE)


##lasso

Gene_lasso <-cv.glmnet(x= as.matrix(Gene[,-1]),y = Gene$target,alpha=1,nfolds=5,family="binomial")
# Looks at Gene_lasso$lambda.min and Gene_lasso$lambda.1se
lambda_Gene <- Gene_lasso$lambda.min
Gene_Coefficients <- coef(Gene_lasso,s=lambda_Gene)
important_features_Gene <- which(Gene_Coefficients != 0)
feature_names_Gene <- colnames(Gene[, -1])
important_feature_names_Gene <- feature_names_Gene[important_features_Gene]
write.csv(important_feature_names_Gene, "C:/Users/matth/anaconda3/envs/data1/data1/Lasso_Boruta_Selection/lihc/Lasso_gene.csv", row.names = FALSE)

## Methyl

#Boruta
Boruta_Methyl <- Boruta(Methyl[-1],Methyl$target)
Methyl_Boruta_features <- getSelectedAttributes(Boruta_Methyl)
write.csv(Methyl_Boruta_features, "C:/Users/matth/anaconda3/envs/data1/data1/Lasso_Boruta_Selection/lihc/Boruta_methyl.csv", row.names = FALSE)

#lasso
Methyl_lasso <-cv.glmnet(x=as.matrix(Methyl[,-1]),y=Methyl$target,alpha=1,nfolds=5,family = "binomial")
#  Looks at  Methyl_lasso$lambda.min and Methyl_lasso$lambda.1se
lambda_Methyl <- Methyl_lasso$lambda.min
Methyl_Coefficients <- coef(Methyl_lasso,s=lambda_Methyl)
important_features_Methyl <- which(Methyl_Coefficients != 0)
feature_names_Methyl <- colnames(Methyl[, -1])
important_feature_names_Methyl <- feature_names_Methyl[important_features_Methyl]
write.csv(important_feature_names_Methyl, "C:/Users/matth/anaconda3/envs/data1/data1/Lasso_Boruta_Selection/lihc/Lasso_methyl.csv", row.names = FALSE)
## miRNA

#Boruta
Boruta_miRNA <- Boruta(miRNA[,-1],miRNA$target)
miRNA_Boruta_features <- getSelectedAttributes(Boruta_miRNA)
write.csv(miRNA_Boruta_features, "C:/Users/matth/anaconda3/envs/data1/data1/Lasso_Boruta_Selection/lihc/Boruta_miRNA.csv", row.names = FALSE)

#lasso
miRNA_lasso <-cv.glmnet(x= as.matrix(miRNA[,-1]),y = miRNA$target,alpha=1,nfolds=5,family="binomial")
##  Looks at  miRNA_lasso$lambda.min and miRNA_lasso$lambda.1se
lambda_miRNA <- miRNA_lasso$lambda.min
miRNA_Coefficients <- coef(miRNA_lasso,s=lambda_miRNA)
important_features_miRNA <- which(miRNA_Coefficients != 0)
feature_names_miRNA <- colnames(miRNA[, -1])
important_feature_names_miRNA <- feature_names_miRNA[important_features_miRNA]
write.csv(important_feature_names_miRNA, "C:/Users/matth/anaconda3/envs/data1/data1/Lasso_Boruta_Selection/lihc/Lasso_miRNA.csv", row.names = FALSE)

