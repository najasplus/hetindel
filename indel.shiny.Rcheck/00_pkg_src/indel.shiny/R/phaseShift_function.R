phaseShift_function <-
function(both_alleles, input_reference, beginning_start, homo_match_len, offset_3p, outp1, input_file, sangerobj){
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
