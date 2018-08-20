sangerseq_function <-
function(input_file, ratio_value, make_chromatogram, outp1){
  #input_file <- "./test/HG_A01_015.ab1"
  input_name <- sub('.*\\/', '', input_file)
  #allele1 <- paste(c("\n> ", input_name, "_1"), collapse = "")
  #allele2 <- paste(c("\n> ", input_name, "_2"), collapse = "")
  outp2 <- list()
  input_read <- readsangerseq(input_file)
  input_base_calls <- makeBaseCalls(input_read, ratio = ratio_value)
  outp2[[1]] <- input_base_calls
  input_prim<- primarySeq(input_base_calls, string = T)
  input_sec<- secondarySeq(input_base_calls, string = T)
  
  prim_list <- unlist(strsplit(input_prim, ""))
  sec_list <- unlist(strsplit(input_sec, ""))
  
  both_alleles <- c()
  
  for (i in 1:(length(prim_list))) {
    if (prim_list[i] == sec_list[i]) both_alleles = c(both_alleles, prim_list[i])
    else {
      if ((prim_list[i] == "A"||sec_list[i] == "A"||prim_list[i] == "W"||sec_list[i] == "W") & 
          (prim_list[i] == "T"||sec_list[i] == "T"||prim_list[i] == "W"||sec_list[i] == "W")) both_alleles = c(both_alleles, "W")
      else if ((prim_list[i] == "C"||sec_list[i] == "C"||prim_list[i] == "Y"||sec_list[i] == "Y") & 
               (prim_list[i] == "T"||sec_list[i] == "T"||prim_list[i] == "Y"||sec_list[i] == "Y")) both_alleles = c(both_alleles, "Y")
      else if ((prim_list[i] == "G"||sec_list[i] == "G"||prim_list[i] == "K"||sec_list[i] == "K") &
               (prim_list[i] == "T"||sec_list[i] == "T"||prim_list[i] == "K"||sec_list[i] == "K")) both_alleles = c(both_alleles, "K")
      else if ((prim_list[i] == "G"||sec_list[i] == "G"||prim_list[i] == "S"||sec_list[i] == "S") &
               (prim_list[i] == "C"||sec_list[i] == "C"||prim_list[i] == "S"||sec_list[i] == "S")) both_alleles = c(both_alleles, "S")
      else if ((prim_list[i] == "A"||sec_list[i] == "A"||prim_list[i] == "M"||sec_list[i] == "M") &
               (prim_list[i] == "C"||sec_list[i] == "C"||prim_list[i] == "M"||sec_list[i] == "M")) both_alleles = c(both_alleles, "M")
      else if ((prim_list[i] == "A"||sec_list[i] == "A"||prim_list[i] == "R"||sec_list[i] == "R") &
               (prim_list[i] == "G"||sec_list[i] == "G"||prim_list[i] == "R"||sec_list[i] == "R")) both_alleles = c(both_alleles, "R")
      else both_alleles = c(both_alleles, "N")
    }
  }
  
  out_str <- paste(input_name, "\nPrimary Base Calls", input_prim, "\nSecondary Base Calls", input_sec, 
                   "\nCombined", paste(both_alleles, collapse =""), sep = "\n")
  
  #param_string1 <- paste0("Signal/Noise Ratio for Peak Detection: ", ratio_value)

  write(out_str, file=outp1[1], append=TRUE)
  
  
  if (make_chromatogram & (length(both_alleles)> 100)) {
    chrom_file <- paste(c(input_file, ".pdf"), collapse = "")
    chromatogram(input_base_calls, trim5 = 30, trim3 = 0, showcalls = "both", 
                 width = 100, height = 2, filename = chrom_file, showhets = TRUE)
  }
  
  outp2[[2]] <- both_alleles
  
  return(outp2)
}
