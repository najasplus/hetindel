files_fun <-
function(archive, prefix_file, submit_count) {
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
