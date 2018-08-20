#library(shiny)
#library(knitr)
#library(readr)

# options(repos = BiocInstaller::biocinstallRepos())
# getOption("repos")
#library(sangerseqR)
#library(Biostrings)
#library(BiocGenerics)
#library("msa")
#library("readr")
#library("crayon")

#source("20180724_shiny_sangerseq_CRISPR_fun_het_al.R")
#source("indelForShiny.R")


ui <- fluidPage(
  
  titlePanel("Analyze Heterozygous Indels"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      fileInput(inputId = "archive",
                label = "Upload your .ab1 sequences together with the reference as a .zip archive:"),
      textInput(inputId = "prefix_file",
                label = "Prefix for output files",
                value = "test"),
      sliderInput(inputId = "ratio_value",
                  label = "Signal/Noise Ratio for Peak Detection:",
                  value = 0.35,
                  min = 0.05, max = 0.75, step = 0.05),
      checkboxInput(inputId = "make_chromatogram",
                    label = "Make Chromatograms",
                    value = TRUE),
      numericInput(inputId = "beginning_start",
                   label = "Offset for matching the 5' homozygous part of the sequence to the reference",
                   value =60),
      numericInput(inputId = "homo_match_len",
                   label = "Length of the homozygous part of the sequence to match the reference",
                   value = 50),
      numericInput(inputId = "offset_3p",
                   label = "Offset at 3' end of the sequence to find matches for heterozygous part",
                   value = 80, min = 40),
      
      actionButton(inputId = "submit", label = "Submit!"),
      
      width = 4
    ),
    
    mainPanel(
      
      
      downloadButton(outputId = "download", label = "Download results"),
      
      span(textOutput(outputId = "error_txt"), style="color:red"),
      textOutput(outputId = "sequences_txt"),
      textOutput(outputId = "match_txt"),
      
      
      doc <- tags$html(
        tags$head(
          tags$title('Analyze Heterozygous Indels')
        ),
        tags$body(
          h2('About'),
          div(id='about', class='simpleDiv',
              'This app is designed to analyze heterozygous and homozygous indels 
              using', 
              strong('.ab1'), 
              'Sanger sequnce files and reference sequence in fasta format as a',
              strong('.txt'), 'file '),
          
          h2('Instructions'),
          h3('Input'),
          div(id='instructions', class='simpleDiv',
              'Create a .zip file containing your .ab1 chromatograms and corresponding reference sequence in
              fasta format as an individual .txt file with a name ending with "reference.txt"'
          ),
          tags$ul(tags$li("The expected indel should occur 70 or more bases from the 5\' and 3\' ends of the sequence"),
                  tags$li('sequences may be in forward and reverse orientation to the reference'), 
                  tags$li('the reference sequence should correspond to the sequenced PCR product'),
                  tags$li('No non-IUPAC characters are allowed in the reference sequence')),
          div(id='instructions', class='simpleDiv',
              'Upload your .zip archive to the app'),
          h3('Parameters'),
          tags$ul(tags$li('you can change the chromatogram peak detection with Signal/Noise Ratio setting. 
                          Reduce it when primary and secondary peaks are different in heigth. Increase when the noise level is high'),
                  tags$li('check "Make Chromatograms" if you want to get the sequence chromatograms as pdf 
                          to analyze peak detection (recommended)'), 
                  tags$li('Choose 5\' offset and sequence length for matching the homozygous part of sequence to the reference. 
                          Use this setting to avoid the areas rich in polymorphysms or bad sequence'),
                  tags$li('Choose 3\' offset and sequence length for matching the homozygous part of sequence to the reference'),
                  tags$li('For more reliable allele determination it is recommeded that tested homozygous and heterozygous 
                          parts of the sequence would be in the vicinity of the indel')),
          h3('Output'),
          div(id='output', class='simpleDiv',
              'As an output you receive a .zip file containing your original files, the chromatograms as .pdf (optionally),
              and two txt files. The fils *_sequences.txt contains the sequence data as a IUPAC codes, and predicted deconvolved alleles' ),
          div(id='output', class='simpleDiv',
              'File *_match.txt contains the information about the parameters used for the sequences processing, and for each individual sequence phase shift between the alleles and alignment of predicted alleles to the reference. '),
          
          div(id='output', class='simpleDiv',
              'If a sequence is homozygous or heterozygous with one of the alleles wild-type, you will receive a sequence 
              alignment to the reference sequence.'),

          div(id='output', class='simpleDiv',
              'Last update: 24 July 2018'),
          
          br()
          ),
        
        
        width = 8
                  ) 
        )
    )
)

server <- function(input, output, session) {
  
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

shinyApp(ui=ui, server = server)
