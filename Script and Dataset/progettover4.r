###UTIILITA'###
write.table(Dataset_finale, file="Dataset_finale_tipolec.csv", quote = F, sep = ";", dec = ",", na = "", row.names = T)
save(Dataset_finale, file="Dataset_finale.RData")
remove(distance_matrix)
str(Dataset_finale, list.len=ncol(df))
is.numeric(Dataset_finale)
remove(Dataset_finale)
colnames(Dataset_finale)
OutVals <- boxplot(Prova)$out
colnames(Dataset_finale)
plot(Dataset$ALL_GSM1637066,  xlim=c(0, 28),  main="Without Outliers", xlab="frequency", pch="*", col="red", cex=2)

##########INSERIMENTO TIPO LEUCEMIA IN COLONNE##########

###TABELLA CORRISPONDENZA DATASET_LEUKEMIA_DEFINITIVO###
Dataset_Leukemia <- cbind(Dataset_Leukemia_definitivo[1], Dataset_Leukemia_definitivo[10])
Dataset_Leukemia <- as.matrix(Dataset_Leukemia)

###SOSTITUZIONE NOMI COLONNE DATA FINALE CON TIPO LEUCEMIA###
for(i in 1:606) {
  for(j in 1:606) {
    if(colnames(Dataset_finale[i]) == Dataset_Leukemia [j,1])
      colnames(Dataset_finale)[colnames(Dataset_finale)==colnames(Dataset_finale[i])] <-  
        if (Dataset_Leukemia[j,2] == "acute_myeloid_leukemia")  {
          paste("AML",colnames(Dataset_finale[i]), sep="_", collapse = NULL) }
    else
      paste("ALL",colnames(Dataset_finale[i]), sep="_", collapse = NULL) 
  }
}

##########PCA##########

###RICOSTRUZIONE DATASET (NA)###
for(i in 1:19340) {
  
  for(j in 1:606) {
    
    if(is.na(Dataset_finale[i,j])){
      
      Dataset_finale[i,j] = (sum(Dataset_finale[i,1:606], na.rm=TRUE)/606)
      
    }
  }
}

write.table(Dataset_finale, file="Dataset_finale_ricostruito.csv", quote = F, sep = ";", dec = ",", na = "", row.names = T)
save(Dataset_finale, file="Dataset_finale_ricostruito.RData")

###N. GENI 19340###

###NORMALIZZAZIONE PRIMA DI PCA###
Dataset_matrix <- as.matrix(Dataset_finale)
normalized_matrix <- normalizeQuantiles(Dataset_matrix, ties=TRUE)
Dataset_normalized <- as.data.frame(normalized_matrix)
Dataset_normalized_T <- t(Dataset_normalized)
Dataset_normalized_T  <- as.data.frame(Dataset_normalized_T)

###ELABORAZIONE DATI PCA###
Dataset.pca <- prcomp(Dataset_normalized_T,center = TRUE)
pred <- predict(Dataset.pca, Dataset_normalized_T_na)
summary(Dataset.pca)
str(Dataset.pca)

###VISUALIZZAZIONE GRAFICA PC1 PC2####
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
typedata <- c(rep("AML",222),rep("ALL",180), rep("AML", 10), rep("ALL", 143), rep("AML",22), rep("ALL",29))
ggbiplot(Dataset.pca, ellipse=TRUE, circle = TRUE, obs.scale = 1, var.scale = 1,var.axes=FALSE, groups=typedata)

##########PRE-PROCESSING##########

#########RICERCA OUTLIER##########
string <-c(rep(0,19340))

for(i in 1:606) {
  Outlier_removal <- Dataset_finale[1:19340,1:2]
  Outlier_removal[1:nrow(Dataset_finale),1] <- rownames(Dataset_finale)
  outlier_values <- boxplot(Dataset_finale[1:nrow(Dataset_finale),i], plot = FALSE, range = 1.1)$out
  
  if(length(outlier_values) != 0 ){
    
    Outlier_removal[1:nrow(Dataset_finale),2] <- Dataset_finale[1:nrow(Dataset_finale),i]
    
    for(j in 1:19340){
      for(m in 1:length(outlier_values)){
        
        if(Outlier_removal[j,2] == outlier_values[m]) {
          
          string[j] <- rownames(Dataset_finale[which(rownames(Dataset_finale) %in% Outlier_removal[j,1]),])
          
        }
      }
    }
    
  }
  
}

###ELIMINAZIONE OUTLIER###
Dataset_finale<- Dataset_finale[-which(rownames(Dataset_finale) %in% string),]

###N. GENI 16836###

###VISUALIZZAZIONE GRAFICA OUTLIER###
plot(Dataset_finale$ALL_GSM1637066,  xlim=c(0, 28),  main="Without Outliers", xlab="frequency", pch="*", col="red", cex=2)
OutVals <- boxplot(Dataset_finale$ALL_GSM1637066)$out

library(limma)

###CALCOLO MATRICE DELLE DISTANZE###
Dataset_finale_T <- t(Dataset_finale)
Dataset_finale_T <- as.data.frame(Dataset_finale_T)
distance_matrix <- dist(Dataset_finale_T)
hc <- hclust(distance_matrix)
plot(hc)

##########NORMALIZZAZIONE##########
Dataset_matrix <- as.matrix(Dataset_finale)
normalized_matrix <- normalizeQuantiles(Dataset_matrix, ties=TRUE)
Dataset_normalized <- as.data.frame(normalized_matrix)
Dataset_normalized_T <- t(Dataset_normalized)
Dataset_normalized_T <- as.data.frame(Dataset_normalized_T)

###CALCOLO MATRICE DELLE DISTANZE POST-NORMALIZZAZIONE###
distance_matrix_normal <- dist(Dataset_normalized_T)
class(distance_matrix_normal)
hc_norm <- hclust(distance_matrix_normal)
plot(hc_norm)

###COSTRUZIONE MATRICE DEI CONTRASTI###
labels <- c(rep("AML",222),rep("ALL",180), rep("AML", 10), rep("ALL", 143), rep("AML",22), rep("ALL",29))
pdata <- as.data.frame(labels)
pdata <- as.factor(pdata$labels)
model_limma <- model.matrix(~-1+labels, data=pdata)
write.table(model_limma, file="model_limma.csv", quote = F, sep = ";", dec = ",", na = "", row.names = T)
colnames(model_limma) <- c("ALL", "AML")
fit <- lmFit(Dataset_normalized, model_limma)
contrasts_names <- "ALL-AML"
contrasts <- makeContrasts(ALL-AML, levels= model_limma)
contrasts_fit <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(contrasts_fit)
genes_names <- nrow(Dataset_normalized)
p_table <- topTable(fit2, number = genes_names, coef = 1)
p_table$score <- sign(p_table$logFC)*-log10(p_table$P.Value)
p_table <- p_table[order(p_table$score, decreasing = TRUE),]
save(p_table, file="p_table.RData")
write.table(p_table, file="p_table.csv", quote = F, sep = ";", dec = ",", na = "", row.names = T)
score_mean <- mean(p_table$score)

##########ESTRAZIONE INDICI GENI CON SCORE > -0.359091823584378############
indexes_g <- which(p_table[7] > -0.359091823584378)
p_table2 <- p_table[indexes_g,]

###CONFRONTO GENI INDIVIDUATI###
idx <- which(rownames(Dataset_finale) %in% rownames(p_table2))
Dataset_finale= Dataset_finale[idx,]
save(Dataset_finale, file="dataset_dopomedia.RData")

###N. GENI 8769###

############BORUTA############
library(Boruta)
set.seed(111)
y=c(rep(0,222),rep(1,180), rep(0, 10), rep(1, 143), rep(0,22), rep(1,29))
boruta.train <- Boruta(y ~ .,data = Dataset_finale_T, doTrace = 2, maxRuns = 100)
print(boruta.train)

dat<-c(rep(0,8769))

for(i in 1:8769) {
  if(boruta.train$finalDecision[i]=="Rejected")
    dat[i] <-  i
}

############CONFRONTO BOURTA GSEA############
data <- Dataset_finale[-dat,]
idx <- which(rownames(data) %in% rownames(Dataset_finale))
Dataset_GSEA_Boruta <- Dataset_finale[idx,]

###N. GENI 92###

########MLP WITH KERAS########

Dataset_finale <- Dataset_GSEA_Boruta

###CARICAMENTO DATASET FINALE###
library(keras)
Dataset_finale_T <- as.data.frame(t(Dataset_finale))
labels <- c(rep("AML",222),rep("ALL",180), rep("AML", 10), rep("ALL", 143), rep("AML",22), rep("ALL",29))
binary_labels <- c(rep(0,222),rep(1,180), rep(0, 10), rep(1, 143), rep(0,22), rep(1,29))
genes <- rownames((Dataset_finale))

###NORMALIZE FUNCTION (NORMALIZZAZIONE PER COLONNA)###
normalize_df <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

###NORMALIZATION (MIN-MAX) PER GENE###
Dataset_finale_norm_per_gene <- as.data.frame(lapply(Dataset_finale_T, normalize_df))
rownames(Dataset_finale_norm_per_gene) <- rownames(Dataset_finale_T)

write.table(Dataset_finale_norm_per_gene, file="Dataset_finale_norm_per_gene.csv", quote = F, sep = ";", dec = ",", na = "", row.names = T)
save(Dataset_finale_norm_per_gene, file = "Dataset_finale_norm_per_gene.RData")

Dataset_finale_norm_per_gene <- cbind(Dataset_finale_norm_per_gene, binary_labels)

###NORMALIZATION (MIN-MAX) PER CAMPIONE###
Dataset_finale_norm_per_campione <- as.data.frame(lapply(Dataset_finale, normalize_df))
rownames(Dataset_finale_norm_per_campione) <- genes
Dataset_finale_norm_per_campione <- as.data.frame(t(Dataset_finale_norm_per_campione))

write.table(Dataset_finale_norm_per_campione, file="Dataset_finale_norm_per_campione.csv", quote = F, sep = ";", dec = ",", na = "", row.names = T)
save(Dataset_finale_norm_per_campione, file = "Dataset_finale_norm_per_campione.RData")

Dataset_finale_norm_per_campione <- cbind(Dataset_finale_norm_per_campione, binary_labels)


##################### DA QUI IN POI IL DATASET È NORMALIZZATO ##############################

###SPLITTING DATASETINTO TRAINING/TESTING WITH KTABLE PER GENE###
library(ade4)
blocks <- c(121,121,121,121,122)
Dataset_finale_norm_per_gene  <- as.data.frame(t(Dataset_finale_norm_per_gene))

Dataset_finale_norm_per_gene <- Dataset_finale_norm_per_gene[,sample(ncol(Dataset_finale_norm_per_gene))]
anal <- ktab.data.frame(Dataset_finale_norm_per_gene, blocks = blocks)
anal[6:12] <- NULL

train_x_1 <- cbind(anal$Ana1, anal$Ana2, anal$Ana3, anal$Ana4)
train_y_1 <-train_x_1[93,]
train_x_1 <- train_x_1[-93,]

train_y_1 <-t(train_y_1)
train_x_1 <- t(train_x_1)

test_y_1 <- anal$Ana5[93,] 
test_x_1 <- anal$Ana5[-93,] 

test_y_1 <-t(test_y_1)
test_x_1 <- t(test_x_1)

train_x_2 <- cbind(anal$Ana1, anal$Ana2, anal$Ana3, anal$Ana5)
train_y_2 <-train_x_2[93,]
train_x_2 <- train_x_2[-93,]

train_y_2 <-t(train_y_2)
train_x_2 <- t(train_x_2)

test_y_2 <- anal$Ana4[93,] 
test_x_2 <- anal$Ana4[-93,] 

test_y_2 <-t(test_y_2)
test_x_2 <- t(test_x_2)

train_x_3 <- cbind(anal$Ana1, anal$Ana2, anal$Ana4, anal$Ana5)
train_y_3 <-train_x_3[93,]
train_x_3 <- train_x_3[-93,]

train_y_3 <-t(train_y_3)
train_x_3 <- t(train_x_3)

test_y_3 <- anal$Ana3[93,] 
test_x_3 <- anal$Ana3[-93,] 

test_y_3 <-t(test_y_3)
test_x_3 <- t(test_x_3)

train_x_4 <- cbind(anal$Ana1, anal$Ana3, anal$Ana4, anal$Ana5)
train_y_4 <-train_x_4[93,]
train_x_4 <- train_x_4[-93,]

train_y_4 <-t(train_y_4)
train_x_4 <- t(train_x_4)

test_y_4 <- anal$Ana2[93,] 
test_x_4 <- anal$Ana2[-93,] 

test_y_4 <-t(test_y_4)
test_x_4 <- t(test_x_4)

train_x_5 <- cbind(anal$Ana2, anal$Ana3, anal$Ana4, anal$Ana5)
train_y_5 <-train_x_5[93,]
train_x_5 <- train_x_5[-93,]

train_y_5 <-t(train_y_5)
train_x_5 <- t(train_x_5)

test_y_5 <- anal$Ana1[93,] 
test_x_5 <- anal$Ana1[-93,] 

test_y_5 <-t(test_y_5)
test_x_5 <- t(test_x_5)

############ A) DATASET NORMALIZZATO PER GENE ##################


############### TEST 1 ###############

#DEFINING THE MODEL
model_norm_per_gene_1 <- keras_model_sequential()

model_norm_per_gene_1 %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_gene_1)

model_norm_per_gene_1 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_gene_1 <- model_norm_per_gene_1 %>% fit(train_x_1, train_y_1, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_gene_1)

score_per_gene_1 <- model_norm_per_gene_1 %>% evaluate(test_x_1, test_y_1)
score_per_gene_1

model_norm_per_gene_1 %>% predict_classes(test_x_1)


############### TEST 2 ###############

#DEFINING THE MODEL
model_norm_per_gene_2 <- keras_model_sequential()

model_norm_per_gene_2 %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_gene_2)

model_norm_per_gene_2 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_gene_2 <- model_norm_per_gene_2 %>% fit(train_x_2, train_y_2, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_gene_2)

score_per_gene_2 <- model_norm_per_gene_2 %>% evaluate(test_x_2, test_y_2)
score_per_gene_2

model_norm_per_gene_2 %>% predict_classes(test_x_2)


############### TEST 3 ###############

#DEFINING THE MODEL
model_norm_per_gene_3 <- keras_model_sequential()

model_norm_per_gene_3  %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_gene_3)

model_norm_per_gene_3 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_gene_3 <- model_norm_per_gene_3 %>% fit(train_x_3, train_y_3, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_gene_3)

score_per_gene_3 <- model_norm_per_gene_3 %>% evaluate(test_x_3, test_y_3)
score_per_gene_3

model_norm_per_gene_3 %>% predict_classes(test_x_3)


############### TEST 4 ###############

#DEFINING THE MODEL
model_norm_per_gene_4 <- keras_model_sequential()

model_norm_per_gene_4 %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_gene_4)

model_norm_per_gene_4 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_gene_4 <- model_norm_per_gene_4 %>% fit(train_x_4, train_y_4, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_gene_4)

score_per_gene_4 <- model_norm_per_gene_4 %>% evaluate(test_x_4, test_y_4)
score_per_gene_4

model_norm_per_gene_4 %>% predict_classes(test_x_4)


############### TEST 5 ###############

#DEFINING THE MODEL
model_norm_per_gene_5 <- keras_model_sequential()

model_norm_per_gene_5 %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_gene_5)

model_norm_per_gene_5 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_gene_5 <- model_norm_per_gene_5 %>% fit(train_x_5, train_y_5, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_gene_5)

score_per_gene_5 <- model_norm_per_gene_5 %>% evaluate(test_x_5, test_y_5)
score_per_gene_5

model_norm_per_gene_5 %>% predict_classes(test_x_5)

######MEDIA KERAS GENE DA STAMPARE
MEDIA_LOSS_GENE = (score_per_gene_1$loss +  score_per_gene_2$loss + score_per_gene_3$loss + score_per_gene_4$loss + score_per_gene_5$loss)/5
MEDIA_ACC_GENE = (score_per_gene_1$acc +  score_per_gene_2$acc + score_per_gene_3$acc + score_per_gene_4$acc + score_per_gene_5$acc)/5

MEDIA_LOSS_GENE
MEDIA_ACC_GENE

###SPLITTING DATASETINTO TRAINING/TESTING WITH KTABLE PER CAMPIONE###
library(ade4)
blocks <- c(121,121,121,121,122)
Dataset_finale_norm_per_campione  <- as.data.frame(t(Dataset_finale_norm_per_campione))

Dataset_finale_norm_per_campione <- Dataset_finale_norm_per_campione[,sample(ncol(Dataset_finale_norm_per_campione))]
anal <- ktab.data.frame(Dataset_finale_norm_per_campione, blocks = blocks)
anal[6:12] <- NULL

train_x_1 <- cbind(anal$Ana1, anal$Ana2, anal$Ana3, anal$Ana4)
train_y_1 <-train_x_1[93,]
train_x_1 <- train_x_1[-93,]

train_y_1 <-t(train_y_1)
train_x_1 <- t(train_x_1)

test_y_1 <- anal$Ana5[93,] 
test_x_1 <- anal$Ana5[-93,] 

test_y_1 <-t(test_y_1)
test_x_1 <- t(test_x_1)

train_x_2 <- cbind(anal$Ana1, anal$Ana2, anal$Ana3, anal$Ana5)
train_y_2 <-train_x_2[93,]
train_x_2 <- train_x_2[-93,]

train_y_2 <-t(train_y_2)
train_x_2 <- t(train_x_2)

test_y_2 <- anal$Ana4[93,] 
test_x_2 <- anal$Ana4[-93,] 

test_y_2 <-t(test_y_2)
test_x_2 <- t(test_x_2)

train_x_3 <- cbind(anal$Ana1, anal$Ana2, anal$Ana4, anal$Ana5)
train_y_3 <-train_x_3[93,]
train_x_3 <- train_x_3[-93,]

train_y_3 <-t(train_y_3)
train_x_3 <- t(train_x_3)

test_y_3 <- anal$Ana3[93,] 
test_x_3 <- anal$Ana3[-93,] 

test_y_3 <-t(test_y_3)
test_x_3 <- t(test_x_3)

train_x_4 <- cbind(anal$Ana1, anal$Ana3, anal$Ana4, anal$Ana5)
train_y_4 <-train_x_4[93,]
train_x_4 <- train_x_4[-93,]

train_y_4 <-t(train_y_4)
train_x_4 <- t(train_x_4)

test_y_4 <- anal$Ana2[93,] 
test_x_4 <- anal$Ana2[-93,] 

test_y_4 <-t(test_y_4)
test_x_4 <- t(test_x_4)

train_x_5 <- cbind(anal$Ana2, anal$Ana3, anal$Ana4, anal$Ana5)
train_y_5 <-train_x_5[93,]
train_x_5 <- train_x_5[-93,]

train_y_5 <-t(train_y_5)
train_x_5 <- t(train_x_5)

test_y_5 <- anal$Ana1[93,] 
test_x_5 <- anal$Ana1[-93,] 

test_y_5 <-t(test_y_5)
test_x_5 <- t(test_x_5)
############ A) DATASET NORMALIZZATO PER  ##################

############### TEST 1 ###############

#DEFINING THE MODEL
model_norm_per_campione_1 <- keras_model_sequential()

model_norm_per_campione_1 %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_campione_1)

model_norm_per_campione_1 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_campione_1 <- model_norm_per_campione_1 %>% fit(train_x_1, train_y_1, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_campione_1)

score_per_campione_1 <- model_norm_per_campione_1 %>% evaluate(test_x_1, test_y_1)
score_per_campione_1

model_norm_per_campione_1 %>% predict_classes(test_x_1)

############### TEST 2 ###############

#DEFINING THE MODEL
model_norm_per_campione_2 <- keras_model_sequential()

model_norm_per_campione_2 %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_campione_2)

model_norm_per_campione_2 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_campione_2 <- model_norm_per_campione_2 %>% fit(train_x_2, train_y_2, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_campione_2)

score_per_campione_2 <- model_norm_per_campione_2 %>% evaluate(test_x_2, test_y_2)
score_per_campione_2

model_norm_per_campione_2 %>% predict_classes(test_x_2)

############### TEST 3 ###############

#DEFINING THE MODEL
model_norm_per_campione_3 <- keras_model_sequential()

model_norm_per_campione_3  %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_campione_3)

model_norm_per_campione_3 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_campione_3 <- model_norm_per_campione_3 %>% fit(train_x_3, train_y_3, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_campione_3)

score_per_campione_3 <- model_norm_per_campione_3 %>% evaluate(test_x_3, test_y_3)
score_per_campione_3

model_norm_per_campione_3 %>% predict_classes(test_x_3)

############### TEST 4 ###############

#DEFINING THE MODEL
model_norm_per_campione_4 <- keras_model_sequential()

model_norm_per_campione_4 %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_campione_4)

model_norm_per_campione_4 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_campione_4 <- model_norm_per_campione_4 %>% fit(train_x_4, train_y_4, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_campione_4)

score_per_campione_4 <- model_norm_per_campione_4 %>% evaluate(test_x_4, test_y_4)
score_per_campione_4

model_norm_per_campione_4 %>% predict_classes(test_x_4)

############### TEST 5 ###############

#DEFINING THE MODEL
model_norm_per_campione_5 <- keras_model_sequential()

model_norm_per_campione_5 %>%
  layer_dense(units=27, activation="relu", input_shape = c(92)) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=8, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units=2, activation="relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 1, activation="sigmoid")

summary(model_norm_per_campione_5)

model_norm_per_campione_5 %>% compile( 
  loss = 'binary_crossentropy',
  optimizer = optimizer_rmsprop(lr = 0.0001),
  metrics = c('accuracy') 
)

history_per_campione_5 <- model_norm_per_campione_5 %>% fit(train_x_5, train_y_5, epochs=150, validation_split = 0.3, batch_size = 32, shuffle = F)
plot(history_per_campione_5)

score_per_campione_5 <- model_norm_per_campione_5 %>% evaluate(test_x_5, test_y_5)
score_per_campione_5

model_norm_per_campione_5 %>% predict_classes(test_x_5)

######MEDIA KERAS CAMPIONE DA STAMPARE
MEDIA_LOSS_CAMPIONE = (score_per_campione_1$loss +  score_per_campione_2$loss + score_per_campione_3$loss + score_per_campione_4$loss + score_per_campione_5$loss)/5
MEDIA_ACC_CAMPIONE = (score_per_campione_1$acc +  score_per_campione_2$acc + score_per_campione_3$acc + score_per_campione_4$acc + score_per_campione_5$acc)/5

MEDIA_LOSS_CAMPIONE
MEDIA_ACC_CAMPIONE

################## SVM ##################
Dataset_finale <- Dataset_GSEA_Boruta
library(e1071)

### PREPARING THE DATA
Dataset_finale_T <- as.data.frame(t(Dataset_finale))
labels <- c(rep("AML",222),rep("ALL",180), rep("AML", 10), rep("ALL", 143), rep("AML",22), rep("ALL",29))

Dataset_finale_T <- cbind(Dataset_finale_T, labels)
Dataset_finale_T<- Dataset_finale_T[sample(nrow(Dataset_finale_T)),]
Data_Name <- Dataset_finale_T[1:606,93]

Dataset_finale_T[,93] <- NULL
Dataset_finale <- as.data.frame(t(Dataset_finale_T))
Data_Name <- as.data.frame(t(Data_Name))

#TRAINING DATA
library(ade4)
blocks <- c(121,121,121,121,122)

anal <- ktab.data.frame(Dataset_finale, blocks = blocks)
anal[6:12] <- NULL

names <- ktab.data.frame(Data_Name, blocks = blocks)
names[6:12] <- NULL

train_svm_x_1 <- as.data.frame(t(cbind(anal$Ana1, anal$Ana2, anal$Ana3, anal$Ana4)))
train_svm_y_1 <- as.factor(t(cbind(names$Ana1, names$Ana2, names$Ana3, names$Ana4)))

test_svm_y_1 <- as.factor(t(names$Ana5))
test_svm_x_1 <- as.data.frame(t(anal$Ana5))


train_svm_x_2 <- as.data.frame(t(cbind(anal$Ana1, anal$Ana2, anal$Ana3, anal$Ana5)))
train_svm_y_2 <- as.factor(t(cbind(names$Ana1, names$Ana2, names$Ana3, names$Ana5)))

test_svm_y_2 <- as.factor(t(names$Ana4))
test_svm_x_2 <- as.data.frame(t(anal$Ana4))


train_svm_x_3 <- as.data.frame(t(cbind(anal$Ana1, anal$Ana2, anal$Ana4, anal$Ana5)))
train_svm_y_3 <- as.factor(t(cbind(names$Ana1, names$Ana2, names$Ana4, names$Ana5)))

test_svm_y_3 <- as.factor(t(names$Ana3))
test_svm_x_3 <- as.data.frame(t(anal$Ana3))


train_svm_x_4 <- as.data.frame(t(cbind(anal$Ana1, anal$Ana3, anal$Ana4, anal$Ana5)))
train_svm_y_4 <- as.factor(t(cbind(names$Ana1, names$Ana3, names$Ana4, names$Ana5)))

test_svm_y_4 <- as.factor(t(names$Ana2))
test_svm_x_4 <- as.data.frame(t(anal$Ana2))


train_svm_x_5 <- as.data.frame(t(cbind(anal$Ana2, anal$Ana3, anal$Ana4, anal$Ana5)))
train_svm_y_5 <- as.factor(t(cbind(names$Ana2, names$Ana3, names$Ana4, names$Ana5)))

test_svm_y_5 <- as.factor(t(names$Ana1))
test_svm_x_5 <- as.data.frame(t(anal$Ana1))

### TEST 1 ###
model_svm_1 <- svm(train_svm_x_1, train_svm_y_1)
summary(model_svm_1)
pred_test_1 <- predict(model_svm_1, test_svm_x_1)
table(pred_test_1, test_svm_y_1)

### TEST 2 ###
model_svm_2 <- svm(train_svm_x_2, train_svm_y_2)
summary(model_svm_2)
pred_test_2 <- predict(model_svm_2, test_svm_x_2)
table(pred_test_2, test_svm_y_2)

### TEST 3 ###
model_svm_3 <- svm(train_svm_x_3, train_svm_y_3)
summary(model_svm_3)
pred_test_3 <- predict(model_svm_3, test_svm_x_3)
table(pred_test_3, test_svm_y_3)

### TEST 4 ###
model_svm_4 <- svm(train_svm_x_4, train_svm_y_4)
summary(model_svm_4)
pred_test_4 <- predict(model_svm_4, test_svm_x_4)
table(pred_test_4, test_svm_y_4)

### TEST 5 ###
model_svm_5 <- svm(train_svm_x_5, train_svm_y_5)
summary(model_svm_5)
pred_test_5 <- predict(model_svm_5, test_svm_x_5)
table(pred_test_5, test_svm_y_5)

######TOTALE SVM DA STAMPARE
test_table = as.factor(c(test_svm_y_1, test_svm_y_2,test_svm_y_3,test_svm_y_4,test_svm_y_5))
pred_table = as.factor(c(pred_test_1, pred_test_2,pred_test_3,pred_test_4,pred_test_5))

table(pred_table, test_table) ###1 = ALL  ###2 = AML

