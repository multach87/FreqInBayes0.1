extract_info <- function(data , n.iter) {
       R2 <- numeric(n.iter)
       model.Int <- data.frame(matrix(nrow = n.iter , ncol = 7))
       model.B <- data.frame(matrix(nrow = n.iter , ncol = 7))
       for(j in 1:5) {
              R2[j] <- data[[(j + 1)]]$Bayes_R2.m1
              model.Int[j , ] <- data[[(j + 1)]]$Summary$fixed[1 , ]
              model.B[j , ] <- data[[(j + 1)]]$Summary$fixed[2 , ]
              
       }
       colnames(model.Int) <- colnames(model.B) <- colnames(data[[(j + 1)]]$Summary$fixed)
       return(list(Bayes_R2.m1 = R2 , 
                   model.Int = model.Int , 
                   model.B = model.B))
}