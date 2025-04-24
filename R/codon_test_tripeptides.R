rm(list = ls())
library(readr)
library(dplyr)

# Uncomment to install torustest
#devtools::install_github("https://github.com/gonzalez-delgado/torustest")

# Load replication data downloaded from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/5P81D4
data_path <- '/path_to_data'
data_precs <- readr::read_csv(paste(c(data_path, "data-precs.csv"), collapse = '/'))

# Get tripeptide information
source("get_tripeptides.R")
data_trip <- get_tripeptides(data_precs, N_cores = 5)

# Amino-acid and tripeptide list
aa_list <- unique(c(data_trip$res_left, data_trip$res, data_trip$res_right))
trip_list <- expand.grid(aa_list, aa_list, aa_list)

# Merge != H,E structures into one category
data_trip$secondary[which(! data_trip$secondary %in% c('H', 'E'))] <- 'Other'

# Angle-based secondary structure classification
data_trip$secondary_angles <- 'Other'
data_trip$secondary_angles[which(data_trip$phi > -180 & data_trip$phi <= 0 & data_trip$psi > -120 & data_trip$psi <= 50)] <- 'H'
data_trip$secondary_angles[which(data_trip$phi > -180 & data_trip$phi <= 0 & data_trip$psi > 50 & data_trip$psi <= 240)] <- 'E'
do_angles <- FALSE # If the angle-based classification has to be considered, set to TRUE

# Remove ambiguous codon assignments
data_trip <- data_trip[-which(data_trip$codon_score < 1),]

# Aggregate redundant data points
data_trip <- data_trip %>% group_by(unp_id, res_id, res, codon, secondary) %>% summarize(
  phi = atan2(sum(sin(phi)), sum(cos(phi)))*180/pi,
  psi = atan2(sum(sin(psi)), sum(cos(psi)))*180/pi)

# Simulate null statistics
N_sim = 2000
sim_null <- torustest::sim.null.stat(N_sim, NC = 4)

# Empty data.frame to fill it with results
data_codon <- data.frame('L' = NA, 'C' = NA, 'R' = NA, 'C1' = NA, 'C2' = NA, 'pv' = NA, 'sec' = NA, 'n1' = NA, 'n2' = NA) 

# Perform tests

for (sec in c('H', 'E', 'Other')){
  
  cat(paste0('Performing tests for ',sec,' structures...\n'))
  
  for (i in 1:nrow(trip_list)){
    
    cat(paste0('Performing goodness-of-fit geodesic test for ',as.vector(trip_list[i,1]),'-',as.vector(trip_list[i,2]),'-',as.vector(trip_list[i,3]),' tripeptide...\n'))
  
    if(!do_angles){ # DSSP-based secondary structure classification
      
      data_i <- data_trip[which(data_trip$res_left == as.vector(trip_list[i,1]) & data_trip$res == as.vector(trip_list[i,2]) & data_trip$res_right == as.vector(trip_list[i,3]) & data_trip$secondary == sec), c('phi', 'psi', 'codon')] # Conformations in "sec" (DSSP) structure corresponding to "aa" amino-acid
    
    }else{ # Angle-based secondary structure classification
      
       data_i <- data_trip[which(data_trip$res_left == as.vector(trip_list[i,1]) & data_trip$res == as.vector(trip_list[i,2]) & data_trip$res_right == as.vector(trip_list[i,3]) & data_trip$secondary_angles == sec), c('phi', 'psi', 'codon')] # Conformations in "sec" (angles-based definition) structure corresponding to "aa" amino-acid
    
    }
       
    if(nrow(data_i) == 0){next}
    
    data_i <- data_i[!(is.na(data_i$phi) | is.na(data_i$psi)), ] # Remove NAs
    if(nrow(data_i) == 0){next}
    
    data_i$phi <- (data_i$phi + 180)/360 # Rescale Ramachandran map to [0,1] x [0,1]
    data_i$psi <- (data_i$psi + 180)/360
    
    codons <- unique(data_i$codon) # List of synonymous codons for "aa" amino-acid
    if('---' %in% codons){codons <- setdiff(codons, '---')}
    if(length(codons) <= 1){next} # Only 1 codon
    
    comb_cod <- t(combn(codons, 2)) # List of combinations of synonymous codons
    
    for (k in 1:nrow(comb_cod)){
      
      cat(paste0('Codon pair ',k,' out of ',nrow(comb_cod),'...\n'))
     
      data_c1 <- data_i[which(data_i$codon == comb_cod[k,1]), c('phi','psi')]
      data_c2 <- data_i[which(data_i$codon == comb_cod[k,2]), c('phi','psi')]
      
      if(nrow(data_c1) < 10 | nrow(data_c2) < 10){next} # Keep only samples with >= 10 points
      
      #p-value for the goodness-of-fit test
      pv_k <- torustest::twosample.geodesic.torus.test(data_c1, data_c2, n_geodesics = 2, NC_geodesic = 2, sim_null = sim_null)
        
      data_i_k <- c(as.vector(trip_list[i,1]), as.vector(trip_list[i,2]), as.vector(trip_list[i,3]), comb_cod[k,1], comb_cod[k,2], pv_k, 
                      sec, nrow(data_c1), nrow(data_c2))

      data_codon <- rbind(data_codon, data_i_k)
        
    }
  }
}


# Results formatting

data_codon <- data_codon[!is.na(data_codon$pv), ] # Remove NA

# Format numeric variables
data_codon$pv <- as.numeric(as.vector(data_codon$pv))
data_codon$n1 <- as.numeric(as.vector(data_codon$n1))
data_codon$n2 <- as.numeric(as.vector(data_codon$n2))

# Set lower-bound for p-values (as computed after MC simulation)
data_codon$pv[which(data_codon$pv < 1/N_sim)] <- 1/N_sim

# BH correction for multiplicity
data_codon$pv_corrected <- NA
for(s in c('H', 'E', 'Other')){data_codon$pv_corrected[which(data_codon$sec == s)] <- p.adjust(data_codon$pv[which(data_codon$sec == s)], method = 'BH')}

# Compute rejection proportions
rej_E <- mean(data_codon$pv_corrected[which(data_codon$secondary == 'E')] < 0.05)
rej_H <- mean(data_codon$pv_corrected[which(data_codon$secondary == 'H')] < 0.05)

# Plot results (Figure S5)
library(ggplot2)
theme_set(theme_bw(base_size = 15))
ggplot(data_codon, aes(x = pv_corrected, col = secondary))+
  stat_ecdf(size = 1.2)+
  labs(x='BH p-value',y='ECDF',col='Secondary structure')+
  theme(legend.position = 'bottom', text = element_text(size = 20),legend.box="vertical")+
  geom_vline(xintercept = 0.05, alpha = 0.5, col = 'darkblue', linetype = 'dashed')+
  annotate("label", x = 0.4, y = 0.88, size = 5, label = paste0('Rejection proportion:\n H: ',round(rej_H,2),'    E: ',round(rej_E,2))) 




