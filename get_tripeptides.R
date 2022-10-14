
#' Get tripeptide information from the dataset used in [1]
#' 
#' From the dataset used in [1], available at https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/5P81D4,
#' the function uses the provided sequence information to extract the identities of neighboring amino-acids for each point.
#'
#' @param data_precs The dataset used in [1], downloaded from https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/5P81D4 and
#' loaded using readr::read_csv function.
#' @param N_cores The number of cores for parallel computation.
#' 
#' @return A dataframe containing as columns:
#' \item{res_left}{The identity of the left neighbor.}
#' \item{res}{The identity of the amino-acid.}
#' \item{res_right}{The identity of the right neighbor.}
#' \item{codon}{The identity of the translated codon.}
#' \item{phi}{The value (in degrees) of the phi angle.}
#' \item{psi}{The value (in degrees) of the psi angle.}
#' \item{secondary}{The secondary structure associated to the conformation.}
#' 
#' @references [1] Rosenberg, A.A., Marx, A. & Bronstein, A.M. Codon-specific Ramachandran plots show amino acid backbone conformation depends on identity of the translated codon. Nat Commun 13, 2815 (2022). https://doi.org/10.1038/s41467-022-30390-9
#' 
#' @export

library(dplyr)

get_tripeptides <- function(data_precs, N_cores){
  
  # List of proteins in database
  prot_list <- unique(data_precs$unp_id)
  
  doParallel::registerDoParallel(N_cores)
  
  trip_data <- foreach::foreach(prot = prot_list, .combine = rbind)%dopar%{
    
    prot_data <- data_precs[which(data_precs$unp_id == prot), ]
    
    data_prot <- foreach::foreach(res = unique(prot_data$res_id), .combine = rbind)%do%{
      
      res_central <- which(prot_data$res_id == res) # Central residue identity
      res_left <- which(prot_data$res_id == res - 1) # Left residue identity
      res_right <- which(prot_data$res_id == res + 1) # Right residue identity
      
      if((length(res_left) == 0 | length(res_right) == 0) | (length(res_left) != length(res_right))){data_res <- rep(NA, 7)}else{
        
        # Tripeptide information
        data_res <- cbind(prot_data$name[res_left], prot_data$name[res_central], prot_data$name[res_right], 
                          prot_data$codon[res_central], prot_data$phi[res_central], prot_data$psi[res_central], 
                          prot_data$secondary[res_central]) 
      } 
      
      data_res # Residue data
    }
    
    data_prot # Protein data
    
  }
  
  # Format tripeptide data
  trip_data <- as.data.frame(trip_data)
  colnames(trip_data) <- c('res_left', 'res', 'res_right', 'codon', 'phi', 'psi', 'secondary')
  
  # Format numeric variables
  trip_data$phi <- as.numeric(as.vector(trip_data$phi))
  trip_data$psi <- as.numeric(as.vector(trip_data$psi))
  
  # Remove NA
  trip_data <- trip_data[!(is.na(trip_data$codon) | is.na(trip_data$phi) | is.na(trip_data$psi)), ]
  
  return(trip_data)
  
}

