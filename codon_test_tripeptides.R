rm(list = ls())
library(readr)

# Load functions from wgof_torus (https://github.com/gonzalez-delgado/wgof_torus)
wgof_torus_path <- '/path_to_downloaded_wgof_torus_scripts' # Folder containing (only) the R scripts of wgof_torus
for(Rfile in list.files(wgof_torus_path)){if(grepl("\\.R$", Rfile)){source(paste(c(wgof_torus_path, Rfile), collapse = '/'))}}

# Load replication data downloaded from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/5P81D4
data_path <- '/path_to_data'
data_precs <- read_csv(paste(c(data_path, "data-precs.csv"), collapse = '/'))

# Get tripeptide information
source("get_tripeptides.R")
data_trip <- get_tripeptides(data_precs, N_cores = 5)

# Amino-acid and tripeptide list
aa_list <- unique(c(data_trip$res_left, data_trip$res, data_trip$res_right))
trip_list <- expand.grid(aa_list, aa_list, aa_list)

# Merge != H,E structures into one category
data_trip$secondary[which(! data_trip$secondary %in% c('H', 'E'))] <- 'Other'

# Define secondary structure using angles
data_trip$secondary_angles <- 'Other'
data_trip$secondary_angles[which(data_trip$phi > -180 & data_trip$phi <= 0 & data_trip$psi > -120 & data_trip$psi <= 50)] <- 'H'
data_trip$secondary_angles[which(data_trip$phi > -180 & data_trip$phi <= 0 & data_trip$psi > 50 & data_trip$psi <= 240)] <- 'E'

# Simulate null statistics
N_sim = 2000
sim_null <- sim.null.stat(N_sim, NC = 4)

# Empty data.frame to fill it with results
data_codon <- data.frame('L' = NA, 'C' = NA, 'R' = NA, 'C1' = NA, 'C2' = NA, 'pv' = NA, 'sec' = NA, 'n1' = NA, 'n2' = NA) 

# Perform tests

for (sec in c('H', 'E', 'Other')){
  
  cat(paste0('Performing tests for ',sec,' structures...\n'))
  
  for (i in 1:nrow(trip_list)){
    
    cat(paste0('Performing goodness-of-fit geodesic test for ',as.vector(trip_list[i,1]),'-',as.vector(trip_list[i,2]),'-',as.vector(trip_list[i,3]),' tripeptide...\n'))

    data_i <- data_trip[which(data_trip$res_left == as.vector(trip_list[i,1]) & data_trip$res == as.vector(trip_list[i,2]) & data_trip$res_right == as.vector(trip_list[i,3]) & data_trip$secondary == sec), c('phi', 'psi', 'codon')] # Conformations in "sec" (DSSP) structure corresponding to "aa" amino-acid
    #data_i <- data_trip[which(data_trip$res_left == as.vector(trip_list[i,1]) & data_trip$res == as.vector(trip_list[i,2]) & data_trip$res_right == as.vector(trip_list[i,3]) & data_trip$secondary_angles == sec), c('phi', 'psi', 'codon')] # Conformations in "sec" (angles-based definition) structure corresponding to "aa" amino-acid
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
      
      if(nrow(data_c1) < 30 | nrow(data_c2) < 30){next} # Keep only samples with >= 30 points
      
      #p-value for the goodness-of-fit test
      pv_k <- twosample.geodesic.torus.test(data_c1, data_c2, n_geodesics = 2, NC_geodesic = 2, sim_null = sim_null)
        
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

# Save results

save_path <- '/path_to_save_results'
saveRDS(data_codon, paste(c(save_path, 'codon_test_tripeptide_results.Rda'), collapse = '/'))





