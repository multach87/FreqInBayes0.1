#Balanced kfold subsetting
kfold_subsetter.vec <- function(n , k , seed = 7) {
       #determine number of larger subsets (when unequal subsets)
       nsams.large <- n %% k
       
       #determine number of smaller subsets (total number when equal subsets)
       nsams.small <- k - nsams.large
       
       #determine sample size of larger subsets (when unequal subsets)
       samsize.large <- ceiling(n / k) * (nsams.large != 0)
       
       #determine sample size of smaller subsets (all subset size when equal subsets)
       samsize.small <- floor(n / k)
       
       #indicator for which subset
       subset.indicator <- c(rep((1 : k) , floor(n / k)) ,
                             rep((1 : (nsams.large) ) , (1 * (nsams.large != 0)) ))
       
       #fix random assignment process
       if(seed) {
              set.seed(seed)
       }
       
       #combine subset indicator with original data  
       kfold.subset <- return(sample(subset.indicator))
}