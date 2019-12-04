create_sim_repped <- function(sim.structure , num.conds = nrow(sim.structure) , 
                              num.pars ,  num.iters) {
       sim.structure.repped <- as.data.frame(matrix(ncol = (num.pars + 2) , 
                                                    nrow = (num.conds * num.iters))) 
       colnames(sim.structure.repped) <- c("beta", "n" , "eta.x" , "eta.y" , 
                                           "seed" , "iter")
       for(i in 1:nrow(sim.structure)) {
              sim.structure.repped[((num.iters * (i - 1)) + 1) : (num.iters * i) , (1 : num.pars)] <- 
                     purrr::map_dfr(seq_len(num.iters) , ~sim.structure[i , ])
       }
       sim.structure.repped[ , "seed"] <- rnorm((num.conds * num.iters))
       sim.structure.repped[ , "iter"] <- rep(1:5)
       return(sim.structure.repped)
}