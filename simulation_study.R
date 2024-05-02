
# Simulation study of the p-value defined in Rosenberg et al. 2022

# This code reproduces the simulation described in Section B of the Supplementary Information
# and produces Figures 1 and S1.

pv_simulation <- function(B, K, Nsim = 200, N = 2000, Nmax = 200, bootstrap = TRUE){
  
  pv_list <- replicate(Nsim,{
  
  nb_list <- c()
  pv_boot <- c()
  
  if(bootstrap){
    xs <- runif(N) # Original samples
    ys <- runif(N)
  }
  
  for(b in 1:B){
    
    if(bootstrap){
      x <- sample(xs, Nmax, replace = TRUE) # b-th bootstrapped samples
      y <- sample(ys, Nmax, replace = TRUE)
    }else{
      x <- stats::runif(Nmax)
      y <- stats::runif(Nmax)
    }
    
    xy <- c(x,y)
    t_b <- c()
    t_sample <- stats::wilcox.test(x,y)$statistic
    
    for(j in 1:K){
      
      xy_j <- sample(xy, 2*Nmax)
      t_bk <- stats::wilcox.test(xy_j[1:Nmax], xy_j[(Nmax+1):(2*Nmax)])$statistic
      t_b <- c(t_b, t_bk)
      
    }
    
    pb <- (sum(t_b > t_sample) + 1)/(K+1)
    pv_boot <- c(pv_boot, pb)
    nb_list <- c(nb_list, sum(t_b > t_sample))
    
  }

  pv_ros <- (sum(nb_list) + 1)/(B*K+1)
  pv_ros}) # p-value defined in Rosenberg et al. 2022

  return(pv_list)}

# Simulate the distribution of p-values for different values of (B, K)

# Independence setting (no bootstrap)

data_01 <- data.frame('pv' = pv_simulation(B = 25, K = 100, bootstrap = FALSE), 'B' = 'B = 25 (indep)', 'K' = 'K = 100')
data_02 <- data.frame('pv' = pv_simulation(B = 25, K = 200, bootstrap = FALSE), 'B' = 'B = 25 (indep)', 'K' = 'K = 200')
data_03 <- data.frame('pv' = pv_simulation(B = 25, K = 300, bootstrap = FALSE), 'B' = 'B = 25 (indep)', 'K' = 'K = 300')
data_sim <- rbind(data_01, data_02, data_03)

# Bootstrap setting

for(b in c(25, 50, 100)){
  for(k in c(100, 200, 300)){
   
    data_bk <- data.frame('pv' = pv_simulation(B = b, K = k, bootstrap = TRUE), 'B' = paste('B =',b), 'K' = paste('K =',k))
    data_sim <- rbind(data_sim, data_bk)
    
  }
}

data_sim$t <- (data_sim$pv - 0.5)*sqrt(12*data_sim$B)
data_sim$B <- ordered(data_sim$B, levels = c('B = 25 (indep)','B = 25', 'B = 50', 'B = 100'))

# Figure S1 in the Supplementary Information

library(latex2exp)
library(ggplot2)

theme_set(theme_bw())
ggplot(data_sim, aes(x = t))+
  geom_histogram(aes(y = ..density..), bins = 30, fill = 'lightblue', col = 'black')+
  geom_density(col = 'darkblue', size = 0.75)+
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), col = 'red', size = 0.75)+
  facet_grid(rows = vars(B), cols = vars(K))+
  labs(x = TeX("$\\sqrt{12B}(\\p_{(c,c'),X)}-\\0.5)$"), y = 'Density')

# Figure 1 in the main text

p1 <- ggplot(data_sim[which(data_sim$B == 'B = 25' & data_sim$K == 'K = 200'),], aes(x = 0.5 + t/sqrt(12*25)))+
  geom_histogram(aes(y = ..density..), bins = 40, fill = 'lightblue', col = 'black')+
  geom_density(col = 'darkblue', size = 1)+
  xlim(c(0,1))+
  stat_function(fun = dnorm, args = list(mean = 0.5, sd = 1/sqrt(12*25)), col = 'red', size = 1)+
  stat_function(fun = dunif, args = list(min = 0, max = 1), col = 'darkgreen', size = 1)+
  labs(x = TeX("$\\p_{(c,c'),X}$"), y = 'Density')+
  theme(text = element_text(size = 20), legend.position = 'none')

p2 <- ggplot(data_sim[which(data_sim$B == 'B = 25' & data_sim$K == 'K = 200'),], aes(x = 0.5 + t/sqrt(12*25)))+
  stat_ecdf(aes(col = 'Empirical'), size = 1)+
  xlim(c(0,1))+
  stat_function(fun = pnorm, args = list(mean = 0.5, sd = 1/sqrt(12*25)), aes(col = 'Gaussian'), size = 1)+
  stat_function(fun = punif, args = list(min = 0, max = 1), aes(col = 'Uniform'), size = 1)+
  labs(x = TeX("$\\p_{(c,c'),X}$"), y = 'CDF')+
  theme(text = element_text(size = 20))+
  scale_colour_manual("Distribution", values = c("darkblue", "red", "darkgreen"))

# Produce two-panel figure:
ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = 'bottom', labels = c('(a)', '(b)'))


