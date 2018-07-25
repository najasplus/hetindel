
# create new directory and unzip files there, create file names for outputs

files_fun <- function(archive, prefix_file, submit_count) {
  warning1 <- c()
  output_prefix_fun <- c()
 if (submit_count) {
  
  if (is.null(archive)) {warning1 <- "Please upload a zip file"
  }
    else {
    extract_dir <- paste0("./", prefix_file)
    dir.create(extract_dir)
    unzip(archive[[4]], exdir = extract_dir)

    sequences_outputs <- paste(c(extract_dir, "/" ,prefix_file, "_sequence_outputs.txt"), collapse = "") 
    match_file <- paste(c(extract_dir, "/", prefix_file, "_match.txt"), collapse = "")
    output_prefix_fun <- c(sequences_outputs, match_file, extract_dir)
    }
 

  output_prefix_fun <- c(output_prefix_fun, warning1)
  
  return (output_prefix_fun)
  }
}

sangerseq_function <- function(input_file, ratio_value, make_chromatogram, outp1){
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

phaseShift_function <- function(both_alleles, input_reference, beginning_start, homo_match_len, offset_3p, outp1, input_file, sangerobj){
  input_name <- sub('.*\\/', '', input_file)
  homo_mismatch <- as.integer(homo_match_len/10)
  both_alleles_str <- paste(both_alleles, collapse ="")
  reference <- DNAString(input_reference)
  rev_ref <- reverseComplement(reference)
  print(input_name)
  # sink(outp1[2], append = T)
  # 
  # print(c("alleles: ", alleles), quote = F)
  # sink()
  
  write(c("\n", input_name), file = outp1[2], append = T)
  
  if (length(both_alleles)<(beginning_start + homo_match_len + offset_3p)) {
    sink(outp1[2], append = T)
    print("Sequence is too short", quote = F)
    sink()
  }
  else {
  match_fw <- matchPattern(DNAString(both_alleles_str, start = beginning_start, nchar = homo_match_len), 
                           reference, fixed = F, max.mismatch = homo_mismatch, with.indels = F)
  match_rev <- matchPattern(DNAString(both_alleles_str, start = beginning_start, nchar = homo_match_len), 
                            rev_ref, fixed = F, max.mismatch = homo_mismatch, with.indels = F)
  
  # check if the beginning of the sequence is a single match to the reference
  
  
  
  if (length(match_fw) == 0 && length(match_rev) == 0) {
    write("\nSequence doesn't match the reference", file = outp1[2], append = T)
  } else if (length(match_fw) > 1 || length(match_rev) > 1) {
    write("\nMatch to the reference is not unambiguous", file = outp1[2], append = T)
    
  } else {       
    #check the orientation of the sequence and pick right orientation
    if (length(match_fw) == 1) {
      refseq <- DNAStringSet(c(Ref = as.character(reference)))
      beginning_match <- match_fw
    }   
    else if (length(match_rev) == 1) {
      refseq <- DNAStringSet(c(Ref_rv = as.character(rev_ref)))
      beginning_match <- match_rev
    }
    
    
    reference_beginning_match <- start(beginning_match) # this position on the reference corresponds to beginning_start position on the tested sequence
    
    start_pos <- nchar(both_alleles_str) - offset_3p # start search x positions away from the 3' end of the sequence 
    match_str_len <- 35 #length of a fragment on the end of the sequence used for matching
    
    #check if the 3' of the sequence is homozygous
    
    num_hom <- sum(both_alleles[start_pos:(start_pos+match_str_len)] %in% c("A", "T", "G", "C"))/match_str_len 
    
    if(num_hom > 0.8) {
      aligned_str <- pairwiseAlignment(refseq, DNAString(both_alleles_str), type="local-global") #local-global doesn't like end gaps
      sink(outp1[2], append = T)
      writePairwiseAlignments(aligned_str)
      sink()
    } else {
      # if the sequence is not homozygous, try to find phase shift
      # initialize end_matches (matches at the end)
      end_matches <- matchPattern(DNAString(both_alleles_str, start = start_pos, match_str_len), refseq[[1]], fixed = F, 
                                  max.mismatch = 3, with.indels = F)
      
      while (match_str_len > 15) {
        if (length(end_matches) == 2) {
          tested_length <- start_pos - beginning_start #length of the analyzed sequence between matches
          ref_length1 <- start(end_matches[1]) - reference_beginning_match # length of a corresponding sequence on a reference sequence for match1 (allele1) and 2
          ref_length2 <- start(end_matches[2]) - reference_beginning_match
          sink(outp1[2], append = T)

          print(c("phase shift:", end(end_matches[2]) - end(end_matches[1])), quote = F)

          alleles <- c(tested_length -ref_length1, tested_length - ref_length2)
          print(c("alleles: ", alleles), quote = F)
          sink()
          

          break
        }
        else {
          start_pos <- start_pos + 3
          match_str_len <- match_str_len -2
          end_matches <- matchPattern(DNAString(both_alleles_str, start = start_pos, match_str_len), refseq[[1]], fixed = F, 
                                      max.mismatch = 3, with.indels = F)
          #print(input_file)
          #end_matches
        }
      }
      
      # if search didn't return proper match
      if (length(end_matches) != 2) {
        write("\nCouldn't solve phase shift", file = outp1[2], append = T)

      }

      both_str_indel <- substr(both_alleles_str, start = 0, stop = nchar(both_alleles_str))

      maxshift <- max(abs(alleles)) + 5
      deconvolved <- calculate.indel(both_str_indel, maxshift)
      seq_vector <- c(refseq, DNAStringSet(unlist(deconvolved)))
      aln <- msa(seq_vector, method = "ClustalOmega", order = "input" )
                   
      
      sink(outp1[2], append = T)
      print(aln, show="complete", showConsensus=FALSE)
      sink()
      
      allele1 <- paste(c("\n> ", input_name, "_1"), collapse = "")
      allele2 <- paste(c("\n> ", input_name, "_2"), collapse = "")
      
      deconv_str <- paste(allele1, deconvolved[1], allele2, deconvolved[2], "\n", sep = "\n")
      
      write(deconv_str, file=outp1[1], append=TRUE)
      
    }
      }
  }
}
