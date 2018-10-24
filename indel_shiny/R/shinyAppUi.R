#' Shiny app server object
#'
#' @import shiny

shinyAppUi <- fluidPage(
  
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