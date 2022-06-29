#reorder
Reorder <- function(res = all_results, mod_det = mod_det){
  #split vectors for reordering
  #host states
  h_states <-colnames(res[, grep("^[^slnbtKm]*$", colnames(res))])
  #uninfected
  zero_endo <- h_states[grep("00", h_states)]
  #only A
  A <- rep(NA, mod_det$endo_no_per_sp)
  for (i in 1:mod_det$endo_no_per_sp){
    A[i] <- h_states[grep(paste0(i, "0"), h_states)]
  }

  #double infected
  if(mod_det$endo_species==2){
  B <- rep(NA, mod_det$endo_no_per_sp)
  for (i in 1:mod_det$endo_no_per_sp){
    B[i] <- h_states[grep(paste0("0", i), h_states)]
    dubs <- paste0("N", c(1:mod_det$endo_no_per_sp) * 11)
    #non_zero co-infections
    to_remove <- c(dubs, B, A, zero_endo)
    h_states_coinf <- h_states[!h_states %in% to_remove]
    h_states_coinf1 <- h_states_coinf[1:length(h_states_coinf) / 2]
    h_states_coinf2 <- tail(h_states_coinf, length(h_states_coinf) / 2)
    col <-c(zero_endo, A, h_states_coinf1, dubs, h_states_coinf2, B)
  }
  return(col)
}

  #concatenate vectors

  if(mod_det$endo_species==1){
    col <-c(zero_endo, A)
  }
  return(col)
}
