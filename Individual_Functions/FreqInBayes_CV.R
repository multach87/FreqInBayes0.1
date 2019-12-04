FreqInBayes.CV <- function(data , freq_K , K , Diffusion , 
                           family_brms = c("gaussian","student","skew_normal")) {
       #print current status
       cat("n = " , data[[4]]$n , "eta.x = " , data[[4]]$eta.x , ", eta.y = " , data[[4]]$eta.y , 
           ", iter = " , data[[4]]$iter , "\n") 
       #make dataframe
       #\\\Only takes data with X and Y pieces currently, 
       ##///needs to be able to take flexibel data structure
       X <- data$X
       Y <- data$Y
       K <- K
       subset.ind <- kfold_subsetter.vec(n = length(X) , k = K)
       cat("subset.ind = " , subset.ind , "\n")
       temp <- data.frame(X , Y , subset = subset.ind)
       out <- list()
       for(i in 1 : (K + 1)) {
              if(i == 1) {
                     out[[i]] <- list(subset = subset.ind)
                     cat("i1 = " , i , "\n")
              } else {
                     #status updater
                     cat("n = " , data[[4]]$n , "eta.x = " , data[[4]]$eta.x , 
                         ", eta.y = " , data[[4]]$eta.y , 
                         ", iter = " ,  data[[4]]$iter , "\n") 
                     cat("i = " , i-1 , "\n")
                     train <- temp[temp$subset != (i-1) , c("X" , "Y")]
                     test <- temp[temp$subset == (i-1) , c("X" , "Y")]
                     #RUN here
                     ###Frequentist Cross-validated regression
                     #NOTE requires package 'caret'
                     data_ctrl <- trainControl(method = "cv", number = freq_K)
                     model_caret <- train(Y ~ X , data = train , trControl = data_ctrl ,
                                          method = "lm" , na.action = na.pass)
                     model_cv_final <- model_caret$finalModel
                     model_cv_final <- summary(model_cv_final)
                     sigma <- model_cv_final$sigma
                     model_coeff <- model_cv_final$coefficients
                     intercept_b <- model_coeff[1 , 1]
                     intercept_SE <- model_coeff[1 , 2]
                     X_b <- model_coeff[2 , 1]
                     X_SE <- model_coeff[2 , 2]
                     
                     
                     ###Diffusion parameter
                     #\\\Diffusion = 0.1 #this should be changeable in the function (0 to 1)
                     Diffusion <- Diffusion*100
                     if (Diffusion == 0){
                            Diffusion <- 1
                     }
                     SD_intercept <- intercept_SE*sqrt(Diffusion)
                     SD_b <- abs(intercept_b*sqrt(Diffusion))
                     sigma <- sigma * sqrt(Diffusion)
                     
                     ###Bayesian section
                     #NOTE: requires package 'brms'
                     #\\\family_brms = "gaussian" #this should be changable (gaussian,student,skew_normal)
                     
                     jj_intercept <- paste("normal(",intercept_b,',',SD_intercept, ")")
                     jj_b <- paste("normal(",X_b,',',SD_b, ")")
                     jj_sigma <- paste("student_t(", sigma, ',', 0, ',', 5, ")")
                     
                     
                     m1priors <- c(
                            prior_string(jj_intercept, class = "Intercept"),
                            prior_string(jj_b, class = "b"),
                            prior_string(jj_sigma, class = "sigma")
                     )
                     
                     
                     #freqin model
                     m1 <- brm(
                            Y~X,
                            data = test,
                            prior = m1priors,
                            family = family_brms,
                            seed = 1,
                            iter = 4000,
                            chains = 4
                     )
                     summ.m1 <- summary(m1)
                     #freqinbayes R2
                     {
                            y.m1 <- get_y(m1) 
                            ypred.m1 <- posterior_linpred(m1, transform = TRUE) 
                            e.m1 <- -1 * sweep(ypred.m1, 2, y.m1) 
                            var_ypred.m1 <- apply(ypred.m1 , 1 , var) 
                            var_e.m1 <- apply(e.m1, 1, var) 
                            Bayes_R2.m1 <- var_ypred.m1 / (var_ypred.m1 + var_e.m1)
                     }
                     
                     out[[i]] <- list(Bayes_R2.m1 = median(Bayes_R2.m1) , 
                                      Summary = summ.m1)
              }
              
       }
       #uninformed bayesian model
       m0 <- brm(Y ~ X, data = temp, 
                 prior = c(prior(normal(0, 1), class = "Intercept"), 
                           prior(normal(0, 1), class = "b", coef = "X"), 
                           prior(student_t(4, 0, 1), class = "sigma")), 
                 seed = 2302
       )
       summ.m0 <- summary(m0)$fixed
       #uninformed R2
       {
              y.m0 <- get_y(m0) 
              ypred.m0 <- posterior_linpred(m0, transform = TRUE) 
              e.m0 <- -1 * sweep(ypred.m0, 2, y.m0) 
              var_ypred.m0 <- apply(ypred.m0 , 1 , var) 
              var_e.m0 <- apply(e.m0, 1, var) 
              Bayes_R2.m0 <- var_ypred.m0 / (var_ypred.m0 + var_e.m0)
              Bayes_R2.m0 <- median(Bayes_R2.m0)
       }
       #frequentist model
       {
              mf <- summary(lm(Y ~ X , data = temp))
              R2.mf <- mf$adj.r.squared
              summ.mf <- mf$coefficients
       }
       FreqInBayes <- extract_info(data = out , n.iter = K)
       
       return(list(Conditions = data[[4]],
                   FreqInBayes = FreqInBayes , 
                   UnInBayes = list(Bayes_R2 = Bayes_R2.m0 , 
                                    model.Int = summ.m0[1 , ] , 
                                    model.B = summ.m0[2 , ]) , 
                   FreqLm = list(Adj_R2 = R2.mf , 
                                 model.Int = summ.mf[1 , ] , 
                                 model.B = summ.mf[2 , ])
       )
       )
}