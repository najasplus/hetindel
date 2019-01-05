#' Shiny app server function
#'
#'

EXAMPLE.DATA.PATH <- "data-raw/example.zip"

shinyAppServer <- function(input, output, session) {
  
  output$downloadTestdata <- downloadHandler(
      filename = function(){
        paste("test_data","zip", sep=".")
      },
      content = function(file){
        file.copy(EXAMPLE.DATA.PATH, file)
      },
      contentType = "application/zip"
    )
  observeEvent(input$submit, {
    # con <- file(paste0(input$prefix_file, ".log"))
    # sink(con, append=TRUE)
    # sink(con, append=TRUE, type="message")
    # print()
    
    outp1 <- files_fun(archive = input$archive, input$prefix_file, submit_count = input$submit)
    
    param_str <- paste0("Signal/Noise Ratio for Peak Detection: ", input$ratio_value, "\n", 
                            "Offset for matching the 5' homozygous part of the sequence to the reference: ", input$beginning_start,"\n",
                            "Length of the homozygous part of the sequence to match the reference: ", input$homo_match_len, "\n",
                            "Offset at 3' end of the sequence to find matches for heterozygous part: ", input$offset_3p, "\n")
    
    write(param_str, file=outp1[2], append=TRUE)
    
    # sink(outp1[2], append = T)
    #     print(param_str, quote = F)
    # sink()
    # 
    
    ref_file <- list.files(path = outp1[3], pattern = "reference.txt", full.names = T, recursive = T)
    if (length(ref_file) != 1){
      output$error_txt <- renderText({outp1[[1]]})
      output$match_txt <- renderText({paste0(" ERROR: Reference file not found. Check that it has correct name.")
      })
      #stop("reference not found")
      #tryCatch(stop("reference not found"))
    } else {
      input_reference <- gsub("[\r\n]", "", (read_file(ref_file)))
      input_reference <- gsub("[\n]", "", input_reference)
      
      seq_list <- list.files(path = outp1[3], pattern = ".ab1$", full.names = T, recursive = T)
      for (j in seq_list) {
        input_file <- j
        withProgress(message = "Analyzing chromatograms", min = 0, max = 1, value = 1, {
          outp2 <- sangerseq_function(input_file, ratio_value = input$ratio_value, 
                                      make_chromatogram = input$make_chromatogram,
                                      outp1 = outp1) })
        #if (length(outp2[[2]])<(input$beginning_start + input$homo_match_len + input$offset_3p)) {next}
        
        withProgress(message = "Calculating allele shift", min = 0, max = 1, value = 1, {
          try(phaseShift_function(both_alleles = outp2[[2]], input_reference, beginning_start = input$beginning_start,
                                  homo_match_len = input$homo_match_len,
                                  offset_3p = input$offset_3p, 
                                  outp1 = outp1, input_file = input_file, sangerobj = outp2 [[1]]) ) })
        
      }
      
      out_zip_name <- paste0(Sys.Date(), "_", input$prefix_file, ".zip")
      
      withProgress(message = "Archiving...", min = 0, max = 1, value = 1, {
        out_archive <- zip(paste0(outp1[3], "/", out_zip_name), list.files(outp1[3], full.names = T)) })
      output$download <- downloadHandler(filename = out_zip_name, 
                                         content = function(file) {
                                           file.copy(paste0(outp1[3], "/", out_zip_name), file)
                                         },
                                         contentType = "application/zip")
      
      
      # observeEvent(input$reset, {
      #   unlink(outp1[3], recursive = T)
      #   unlink(out_zip_name)
      # })
      
      
      
      output$sequences_txt <- renderText({
        paste(out_zip_name)
      })
      
      output$match_txt <- renderText({
        "Now you can download the output files as a .zip archive "
      })  
    }
    
    session$onSessionEnded(function() { unlink(outp1[3], recursive = TRUE) } )
    
  }, once = F) 
  
  
  
}