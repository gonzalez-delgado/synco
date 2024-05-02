rm(list = ls())
library(readr)

# Uncomment to install torustest
#devtools::install_github("https://github.com/gonzalez-delgado/torustest")

# Experimental data downloaded from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/5P81D4
data_path <- '/path_to_data'

# AlphaFold2 data in the "nt_structure_2024.txt" file (don't needed if AF2 data is not used)
AF2_path <- '/path_to_AF2_data' # Path to the file "nt_structure_2024.txt".

# Load experimental data
data_precs <- readr::read_csv(paste(c(data_path, "data-precs.csv"), collapse = '/'))
aa_list <- unique(data_precs$name) # List of amino-acids

do_AF2 <- FALSE # To run the analysis on the AlphaFold2 database, set to TRUE

if(do_AF2){
  
  # Load AF2 database
  AF2_data <- readr::read_delim(paste(c(AF2_path,"nt_structure_2024.txt"), collapse = '/'), delim = "\t", escape_double = FALSE, trim_ws = TRUE)

  AF2_data <- AF2_data[which(AF2_data$UniProt %in% unique(data_precs$unp_id)), ] # Keep only proteins present in the experimental database
  AF2_data <- AF2_data[-which(AF2_data$pLDDT24 < 90),] # Filter by pLDDT score
  AF2_data <- AF2_data[-which(AF2_data$Omega < 20 & AF2_data$Omega > -20), ] # Remove CIS conformations
  AF2_data <- AF2_data[, c('UniProt', 'Position', 'Aa', 'Codon', 'Structure', 'Phi', 'Psi')]
  colnames(AF2_data) <- colnames(data_precs)
  data_precs <- AF2_data # Rename database

}

# Merge != H,E structures into one category
data_precs$secondary[which(! data_precs$secondary %in% c('H', 'E'))] <- 'Other'

# Secondary structure angle-based definition
data_precs$secondary_angles <- 'Other'
data_precs$secondary_angles[which(data_precs$phi > -180 & data_precs$phi <= 0 & data_precs$psi > -120 & data_precs$psi <= 50)] <- 'H'
data_precs$secondary_angles[which(data_precs$phi > -180 & data_precs$phi <= 0 & data_precs$psi > 50 & data_precs$psi <= 240)] <- 'E'
do_angles <- FALSE # If the angle-based classification has to be considered, set to TRUE

# Remove ambiguous codon assignments
data_precs <- data_precs[-which(data_precs$codon_score < 1),]

# Aggregate redundant data points
data_precs <- data_precs %>% group_by(unp_id, res_id, name, codon, secondary) %>% summarize(
  phi = atan2(sum(sin(phi)), sum(cos(phi)))*180/pi,
  psi = atan2(sum(sin(psi)), sum(cos(psi)))*180/pi)

# Simulate null statistics
N_sim = 2000
sim_null <- torustest::sim.null.stat(N_sim, NC = 4)

# Empty data.frame to fill it with results
data_codon <- data.frame('AA' = NA, 'C1' = NA, 'C2' = NA, 'pv' = NA, 'sec' = NA, 'n1' = NA, 'n2' = NA) 

# Perform tests

for (sec in c('H', 'E', 'Other')){
  
  cat(paste0('Performing tests for ',sec,' structures...\n'))
  
  for (aa in aa_list){
    
    cat(paste0('Performing tests for ',aa,' amino-acid...\n'))
    
    if(!do_angles){ # DSSP-based secondary structure classification
      
      data_aa <- data_precs[which(data_precs$name == aa & data_precs$secondary == sec), c('phi', 'psi','codon')] # Conformations in "sec" (DSSP) structure corresponding to "aa" amino-acid
    
    }else{ # Angle-based secondary structure classification
      
      data_aa <- data_precs[which(data_precs$name == aa & data_precs$secondary_angles == sec), c('phi', 'psi','codon')] # Conformations in "sec" (angles-based definition) structure corresponding to "aa" amino-acid
    
    }
    
    if(nrow(data_aa) == 0){next}
    
    data_aa <- data_aa[!(is.na(data_aa$phi) | is.na(data_aa$psi)), ] # Remove NAs
    if(nrow(data_aa) == 0){next}
    
    data_aa$phi <- (data_aa$phi + 180)/360 # Rescale Ramachandran map to [0,1] x [0,1]
    data_aa$psi <- (data_aa$psi + 180)/360
    
    codons <- unique(data_aa$codon) # List of synonymous codons for "aa" amino-acid
    if('---' %in% codons){codons <- setdiff(codons, '---')}
    if(length(codons) <= 1){next} # Only 1 codon
    
    comb_cod <- t(combn(codons, 2)) # List of combinations of synonymous codons

    for (k in 1:nrow(comb_cod)){
      
      cat(paste0('Codon pair ',k,' out of ',nrow(comb_cod),'...\n'))
     
      data_c1 <- data_aa[which(data_aa$codon == comb_cod[k,1]), c('phi','psi')]
      data_c2 <- data_aa[which(data_aa$codon == comb_cod[k,2]), c('phi','psi')]
      
      if(nrow(data_c1) < 30 | nrow(data_c2) < 30){next} # Keep only samples with >= 30 points
    
      #p-value for the goodness-of-fit test
      pv_k <- torustest::twosample.geodesic.torus.test(data_c1, data_c2, n_geodesics = 2, NC_geodesic = 2, sim_null = sim_null)
     
      data_aa_k <- c(aa, comb_cod[k,1], comb_cod[k,2], pv_k, sec, nrow(data_c1), nrow(data_c2))
      data_codon <- rbind(data_codon, data_aa_k)
        
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
saveRDS(data_codon, paste(c(save_path, 'codon_test_results.Rda'), collapse = '/'))

