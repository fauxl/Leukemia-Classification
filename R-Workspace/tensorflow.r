library(keras)

############### TEST 1 ###############

#DEFINING THE MODEL
model_norm_per_gene_1 <- keras_model_sequential()

model_norm_per_gene_1 %>%
  layer_dense(units=216, activation="relu", input_shape = c(720)) %>%
  layer_dropout(rate = FLAGS$dropout1) %>%
  layer_dense(units=64, activation="relu") %>%
  layer_dropout(rate = FLAGS$dropout2) %>%
  layer_dense(units=19, activation="relu") %>%
  layer_dropout(rate = FLAGS$dropout3) %>%
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

