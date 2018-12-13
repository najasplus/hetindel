library("optparse")
library("readr")
library("crayon")
library("seqGen")

#make command line options 
option_list <- list(
	make_option(c("-i","--input-file"), type="character",default=NULL,
		    help="file path to sequence data. If set -s is ignored."),
	make_option(c("-s","--string"), type="character",default=NULL,
		    help="sequence string"),
	make_option(c("-n","--number"), type="integer",default=1,help="number of sequence in file. 
		[default %default"),
	make_option(c("-m","--max-shift"), type="integer",default=15,
		help="The maximum potential phase shift size to be considered during analysis; 
		    can be set up to 1/2 of the length of the analyzed fragment. 
		    For fragments with multiple phase shifts, 
		    maximum phase shift of 1/10 of the length of the analyzed fragment is recommended.
		    [default %default]",metavar="number"),
	make_option(c("-p","--penalty"),type="integer",default=2,
		help="The cost of transition from one phase shift to another 
		(equivalent of the gap opening cost in alignments)[default %default].",metavar="number"),
	make_option(c("-f","--fix-shifts"),type="character",default=NULL,
		help="Restricts analysis to phase shifts of selected magnitudes 
		    (usually for re-analysis). Separate multiple values by a comma.",
		metavar="list of numbers"),
	make_option(c("-r","--right"), action="store_true",default=TRUE,
		help="Display gaps aligned to the right. [default]"),
	make_option(c("-l","--left"),action="store_false", dest="right",
		help="Display gaps aligned to the left."),
	make_option(c("-a","--align"),action="store_true",default=TRUE, 
		help="Option for displaying reconstructed allelic sequences. [default]"),
	make_option(c("-A","--not-align"),action="store_false",dest="align",
		help="Option for not displaying reconstructed allelic sequences."),
	make_option(c("-g","--long-indel"),action="store_true",default=FALSE,
		help="Display alternative, longer variants of reconstructed insertions.
		[default \"%default\"]"),
	make_option(c("-G","--generate"),type="list",default=NULL,
		help="Generate a random BioString with introduced shifts in IUPAC Notation. 
		Input as: [length, exact, occurence] length: length of sequence; one of 150, 300 ,700.
		exact: whether or not the shift is exact or has a variance in occurence; one of TRUE, FALSE.
		occurence: occurence of the shift, one of b (begin), m (middle), q1 (first quater), q3 (third quater)." )
	)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


options(error=traceback)
countCharOccurrences <- function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
}



alignLeft <- function(matrix,seq.length,seq.matrix){
	z1 <- 0
	for(i in seq.length:2){
		matrix[4,i] <- 0
		if((seq.length - i-1 - matrix[2,i]) >0){
			j <- matrix[2,i+1]
			x <- matrix[2,i]
			if(x<j){
				matrix[4,i] <- 1
				if(any(matrix[1,i]==c(0,1))&&any(matrix[1,i+j]==c(0,1))&&
					seq.matrix[1,i-1]==seq.matrix[2,i+j-1]){
					##########################################
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(2,0))&&
					seq.matrix[2,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(1,0))&&
					seq.matrix[2,i-1]==seq.matrix[2,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(2,0))&&
					seq.matrix[1,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j
				}else{
					matrix[4,i] <- 0
				}
			}else if(x>0&&j<0){
				matrix[4,i] <- 1
				if(i<j){
					matrix[4,i] <- 0
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i-j]==c(1,0))&&
					seq.matrix[1,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i-j]==c(2,0))&&
					seq.matrix[2,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i-j]==c(1,0))&&
					seq.matrix[1,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i-j]==c(2,0))&&
					seq.matrix[2,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <-j
				}else{
					matrix[4,i] <- 0
				}
			}else if( x>j){
				matrix[4,i] <- 1
				if(i<j){
					matrix[4,i] <- 0
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(1,0))&&
					any(matrix[1,i-j]==c(1,0))&&seq.matrix[1,i-1]==seq.matrix[2,i+j-1]&&
					seq.matrix[1,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(1,0))&&
					any(matrix[1,i-j]==c(2,0))&&seq.matrix[1,i-1]==seq.matrix[2,i+j-1]&&
					seq.matrix[2,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(2,0))&&
					any(matrix[1,i-j]==c(1,0))&&seq.matrix[2,i-1]==seq.matrix[1,i+j-1]&&
					seq.matrix[1,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(2,0))&&
					any(matrix[1,i-j]==c(2,0))&&seq.matrix[2,i-1]==seq.matrix[1,i+j-1]&&
					seq.matrix[2,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(1,0))&&
					any(matrix[1,i-j]==c(1,0))&&seq.matrix[2,i-1]==seq.matrix[2,i+j-1]&&
					seq.matrix[1,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(1,0))&&
					any(matrix[1,i-j]==c(2,0))&&seq.matrix[2,i-1]==seq.matrix[2,i+j-1]&&
					seq.matrix[2,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(2,0))&&
					any(matrix[1,i-j]==c(1,0))&&seq.matrix[1,i-1]==seq.matrix[1,i+j-1]&&
					seq.matrix[1,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(2,0))&&
					any(matrix[1,i-j]==c(2,0))&&seq.matrix[1,i-1]==seq.matrix[1,i+j-1]&&
					seq.matrix[2,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else{
					matrix[4,i] <- 0
				}
			}
		}
		if(matrix[4,i] > 0){
			z1 <- z1+1
		}
	}
	return(list(matr=matrix,val=z1))
}
alignRight <- function(matrix,seq.length,seq.matrix){
	z1 <- 0
	for(i in 2:seq.length+1){
		matrix[4,i] <- 0
		if((i- matrix[2,i]) >0){
			j <- matrix[2,i-1]
			x <- matrix[2,i]
			if(x<j){
				matrix[4,i] <- 1
				if(i<j){
					matrix[4,i] <- 0
				}else  if(any(matrix[1,i-j]==c(1,0))&&any(matrix[1,i]==c(1,0))&&
					seq.matrix[1,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else  if(any(matrix[1,i-j]==c(2,0))&&any(matrix[1,i]==c(2,0))&&
					seq.matrix[2,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else  if(any(matrix[1,i-j]==c(1,0))&&any(matrix[1,i]==c(2,0))&&
					seq.matrix[1,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else  if(any(matrix[1,i-j]==c(2,0))&&any(matrix[1,i]==c(1,0))&&
					seq.matrix[2,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else{
					matrix[4,i] <- 0
				}
			}else if(j>0 && x<0){
				matrix[4,i] <- 1
				if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(1,0))&&
					seq.matrix[1,i-1]==seq.matrix[2,i+j-1]){
					matrix[2,i] <- j 
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(2,0))&&
					seq.matrix[1,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j 
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(1,0))&&
					seq.matrix[2,i-1]==seq.matrix[2,i+j-1]){
					matrix[2,i] <- j 
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(2,0))&&
					seq.matrix[1,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j 
				}else{
					matrix[4,i] <- 0
				}
			}else if( x>j){
				matrix[4,i] <- 1
				if(i<j){
					matrix[4,i] <- 0
				}else if(any(matrix[1,i-j]==c(1,0))&&any(matrix[1,i]==c(1,0))&&
					any(matrix[1,i+j]==c(1,0))&&seq.matrix[1,i-j-1]==seq.matrix[2,i-1]&&
					seq.matrix[1,i-1]==seq.matrix[2,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i-j]==c(1,0))&&any(matrix[1,i]==c(1,0))&&
					any(matrix[1,i+j]==c(2,0))&&seq.matrix[1,i-j-1]==seq.matrix[2,i-1]&&
					seq.matrix[1,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i-j]==c(2,0))&&any(matrix[1,i]==c(2,0))&&
					any(matrix[1,i+j]==c(1,0))&&seq.matrix[2,i-j-1]==seq.matrix[1,i-1]&&
					seq.matrix[2,i-1]==seq.matrix[2,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i-j]==c(2,0))&&any(matrix[1,i]==c(2,0))&&
					any(matrix[1,i+j]==c(2,0))&&seq.matrix[2,i-j-1]==seq.matrix[1,i-1]&&
					seq.matrix[2,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i-j]==c(1,0))&&any(matrix[1,i]==c(2,0))&&
					any(matrix[1,i+j]==c(1,0))&&seq.matrix[1,i-j-1]==seq.matrix[1,i-1]&&
					seq.matrix[2,i-1]==seq.matrix[2,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i-j]==c(1,0))&&any(matrix[1,i]==c(2,0))&&
					any(matrix[1,i+j]==c(2,0))&&seq.matrix[1,i-j-1]==seq.matrix[1,i-1]&&
					seq.matrix[2,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i-j]==c(2,0))&&any(matrix[1,i]==c(1,0))&&
					any(matrix[1,i+j]==c(1,0))&&seq.matrix[2,i-j-1]==seq.matrix[2,i-1]&&
					seq.matrix[1,i-1]==seq.matrix[2,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i-j]==c(2,0))&&any(matrix[1,i]==c(1,0))&&
					any(matrix[1,i+j]==c(2,0))&&seq.matrix[2,i-j-1]==seq.matrix[2,i-1]&&
					seq.matrix[1,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j
				}else{
					matrix[4,i] <- 0
				}
			}
		}
		if(matrix[4,i] >0 ){
			z1 <- z1+1
		}
	}
	for(i in 1:seq.length){
		if(matrix[4,i]>0){
			if(matrix[4,i-1]>0){
				matrix[4,i] <- matrix[4,i-1]
			}else{
				for(j in i:seq.length){
					if(matrix[4,j]==0){
						break
					}
				}
				matrix[4,i] <-j-i
			}
		}
	}
	return(list(matr=matrix,value=z1))
}


#print(opt$input)
#print(opt$string)
#print(opt$max)
#print(opt$penalty)
#print(opt$fix)

#print(opt$right)
#print(opt$left)
#print(opt$align)
#print(opt$not)
#print(opt$long)


seq.number <- opt$number

#check if input or string is set
if(!is.null(opt$input)){
	#parse input file
	sequence <- trimws(readLines(opt$input))

	#get combined sequence
	sequence <- sequence[grep("/.",sequence)+1]

	if(length(sequence)<seq.number){
		print_help(opt_parser)
		stop("File does not contain ",seq.number," sequences, but only ",
		length(sequence),".",call.=FALSE)
	}else{
		sequence <- sequence[seq.number]
	}
}else if (!is.null(opt$string)){
	sequence <- opt$string
}else if(length(opt$generate)>0){
	params <- unlist(strsplit(opt$generate,","))
	if(length(params) != 3){
		print_help(opt_parser)
		stop("Please provide exactly 3 parameters for generate in the form of: %length, %exact, %occurence",call=FALSE)
	}
	generate.length <- as.numeric(params[1])
	if(!any(generate.length==c(150,300,700))){
		print_help(opt_parser)
		stop("Please only choose length of 150, 300, or 700 for generated Sequences.",call=FALSE)
	}
	generate.exact <- as.logical(params[2])
	generate.occurence <- params[3]
	if(!any(generate.occurence==c("m","b","q1","q3"))){
		print_help(opt_parser)
		stop("Please only choose occurence of \"b\", \"m\", \"q1\", \"q3\".",call=FALSE)
	}
	sequence <- generateRandBioString(generate.length,generate.exact,generate.occurence)
}else{
	print_help(opt_parser)
	stop("Please provide an input file or string.\n",call.=FALSE)
}	

#set variables
max.shift <- opt$max
penalty <- opt$penalty
if(!is.null(opt$fix)){
	fixed.shifts <- as.numeric(unlist(strsplit(opt$fix, split=",")))
	}else{fixed.shifts <- NULL}
is.right <- opt$right
is.align <- opt$align
is.longindel <- opt$long

#convert sequence to all uppercase
sequence <- toupper(sequence)

seq.asvector=strsplit(sequence,"")[[1]]

ambig <- 0

seq.length <- nchar(sequence)
#maximum shift is limited to half of sequence length
if(max.shift*2>seq.length){
	max.shift=floor(nchar(sequence)/2)
}
if (seq.length!=0){ #if sequence exists
	

	#set scores as in indelligent
	match.score=1
	mismatch.score=0
	

	#initialize matrix seq.matrix
	seq.matrix <- matrix(,nrow=2,ncol=seq.length)
	#fill seq.matrix
	pos.of=seq.asvector=="R"
	seq.matrix[1,pos.of]="A"
	seq.matrix[2,pos.of]="G"
	pos.of=seq.asvector=="Y"
	seq.matrix[1,pos.of]="C"
	seq.matrix[2,pos.of]="T"
	pos.of=seq.asvector=="S"
	seq.matrix[1,pos.of]="G"
	seq.matrix[2,pos.of]="C"
	pos.of=seq.asvector=="W"
	seq.matrix[1,pos.of]="A"
	seq.matrix[2,pos.of]="T"
	pos.of=seq.asvector=="M"
	seq.matrix[1,pos.of]="A"
	seq.matrix[2,pos.of]="C"
	pos.of=seq.asvector=="K"
	seq.matrix[1,pos.of]="G"
	seq.matrix[2,pos.of]="T"
	#fill remaining characters
	pos.of=is.na(seq.matrix[1,])
	seq.matrix[1,pos.of]=seq.asvector[pos.of]
	seq.matrix[2,pos.of]=seq.asvector[pos.of]
	#check for 3-fold degenerate bases
	if (!(seq.matrix=="A" || seq.matrix=="T" || seq.matrix=="G" || seq.matrix=="C")){
		bd <- TRUE
		ambig1 <- sum(!(seq.matrix=="A" | seq.matrix=="T" | seq.matrix=="G" | seq.matrix=="C"))
	}else{
		bd <- FALSE
	}


	
#	print(seq.matrix)




	#init calculation matrices
	vorward.score <- array(0,dim=c(seq.length+2,max.shift+1,2))
	scores <- array(0,dim=c(seq.length+2,max.shift+1,2))
	
	#init conditions for vorward.score and scores (first 1 to max.shift+1 as the corresponding index 
	#in vorward.score, last max.shift to seq.length+1 as falling cooresponding index e.g. 10 9 8 ... 0)
	for ( i in 2:(max.shift+1)){
		vorward.score[i,i:(max.shift+1),] <- (i-1)
		scores[(seq.length-i+3),i:(max.shift+1),] <- (i-1)
	}
		
#	print(vorward.score)
#	print(scores)	
	#calculate vorward.score in forward direction
	for ( i in 1:seq.length){
		jMax <- max.shift
		if ((i-1) < max.shift){	#hier evtl fehler?
			jMax <- i-1
		}
		for (j in 1:(jMax+1)){
			for (z in 1:2){
				SCMax <- -111111
				for (x in 1:(jMax+1)){
					if (j==1 && x==1 && seq.matrix[1,i] != seq.matrix[2,i]){
						
						SC <- max(vorward.score[i,x,1],vorward.score[i,x,2]) - mismatch.score
					}else if (j==1 && seq.matrix[1,i]==seq.matrix[2,i]){
						SC <- max(vorward.score[i,x,1],vorward.score[i,x,2]) + match.score
					}else if (j != x){
						SC <- max(vorward.score[i,x,1],vorward.score[i,x,2]) + match.score
					}else if (((seq.matrix[3-z,i]==seq.matrix[1,i-(j-1)]) && 
					(vorward.score[i-(j-2),j,1] >= vorward.score[i-(j-2),j,2])) || ((seq.matrix[3-z,i]==seq.matrix[2,i-(j-1)]) &&
					( vorward.score[i-(j-2),j,1] <= vorward.score[i-(j-2),j,2]))){
						SC  <- max(vorward.score[i,j,1],vorward.score[i,j,2])+match.score
					} else {
						SC  <-  max(vorward.score[i,j,1],vorward.score[i,j,2])- mismatch.score
					}
					if (x!=j){
						SC <- SC - penalty - abs(x-j)
					}
					if (SC > SCMax){
						SCMax  <-  SC
					}
				}
				vorward.score[i+1,j,z]  <- SCMax
			}
		}
	}

#	print(vorward.score)
	#calculate scores in reverse direction
	for (i in seq.length:1){
		jMax <- max.shift
		if((seq.length-i) < max.shift){
			jMax <- seq.length - i
		}
		for ( j in 1:(jMax+1)){
			for (z in 1:2){
				SCMax <- -111111
				for (x in 1:(jMax+1)){
					if(j==1 && x == 1 && seq.matrix[1,i] != seq.matrix[2,i]){
						SC <- max(scores[i+2,x,1],scores[i+2,x,2]) -mismatch.score
					}else if(j==1 && seq.matrix[1,i]==seq.matrix[2,i]){
						SC <- max(scores[i+2,x,1],scores[i+2,x,2])+match.score
					}else if(j!=x){
						SC=max(scores[i+2,x,1],scores[i+2,x,2]) + match.score
					}else if(((seq.matrix[z,i]==seq.matrix[2,i+(j-1)]) && (scores[i+j,j,1] >= scores[i+j,j,2])) || 
						((seq.matrix[z,i]==seq.matrix[1,i+(j-1)]) && (scores[i+j,j,1] <=scores[i+j,j,2]))){
						SC <- max(scores[i+2,j,1],scores[i+2,j,2]) +match.score
					}else{
						SC <- max(scores[i+2,j,1],scores[i+2,j,2]) - mismatch.score
					}
					if (x!=j){
						SC <- SC-penalty-abs(x-j)
					}
					if(SC > SCMax){
						SCMax <- SC
					}
				}
				scores[i+1,j,z]=SCMax
			}
		}
	}
#	print(scores)
	#add vorward.score and scores store result in scores
	scores <- vorward.score+scores

	
	#make list of possible indels
	x <- c()
	for (i in 2:(seq.length+1)){
		SC <- -111111
		SCc <- 0
		for (j in 0:max.shift){
			if(is.null(fixed.shifts) || any(fixed.shifts==j)){
				if (scores[i,j+1,1]> SC){
					SC <- scores[i,j+1,1]
					SCc <- j
				}
				if (scores[i,j+1,2]>SC){
					SC <- scores[i,j+1,2]
					SCc <- j
				}
			}
		}
		if(!any(x==SCc)){
			x <- c(x,SCc)
		}
	}
#	print(x)
	#apply scores from scores and resolve positions in seq.matrix. Results are stored in phase.shift.matrix.
	#phase.shift.matrix[1,i] the value indicating which base from seq.matrix goes into the first string. 
	#If phase.shift.matrix[1,i]==0 then the position is ambiguous.
	#phase.shift.matrix[2,i] the value indicating the phase shift which resolves the position.

	SCd <- 0
	phase.shift.matrix <- matrix(0,nrow=4,ncol=(seq.length+1))
	repeat{
		for(i in 1:seq.length){
			SC <- -111111
			SCa <- 0	#possibility to solve in position 1
			SCb <- 0	#possibility to solve in position 2
			SCc <- 0	#shift size
			
			for(j in 0:max.shift){
				if((is.null(fixed.shifts) && SCd==0) || (any(fixed.shifts==j) && SCd==0) || 
					(any(x==j) && SCd==1)){
					if(scores[i+1,j+1,1] > SC){
						SC <- scores[i+1,j+1,1]
						SCa <- 0
						SCb <- 0
						SCc <- j
					}
					if(scores[i+1,j+1,2] > SC){
						SC <- scores[i+1,j+1,2]
						SCa <- 0
						SCb <- 0
						SCc <- j
					}
					if(((scores[i+1,j+1,1]==SC) || (scores[i+1,j+1,2]==SC)) && (0+j)==phase.shift.matrix[2,i]){
						SCc <- j
					}
					if((scores[i+1,j+1,1]==SC) && (any(x==j))){
						SCa <- 1
					}
					if((scores[i+1,j+1,2]==SC) && (any(x==j))){
						SCb <- 1
					}
				}
			}
		
			if(seq.matrix[1,i]==seq.matrix[2,i]){
				phase.shift.matrix[1,i+1] <- 1
			}else if (SCa==1 && SCb==0){
				phase.shift.matrix[1,i+1] <- 1
			}else if (SCb==1 && SCa ==0){
				phase.shift.matrix[1,i+1] <- 2
			}else{
				phase.shift.matrix[1,i+1] <- 0
			}
			phase.shift.matrix[2,i+1] <- SCc
		}
		#remove phase shifts recovered at the number of consecutive positions smaller 
		#than the phase shift magnitude

		
		u <- c()
		SCa <- 0

		for(i in 1:seq.length){
			if(phase.shift.matrix[2,i+1] != phase.shift.matrix[2,i]){
				for(z in 1:phase.shift.matrix[2,i+1]){
					if(i+z > seq.length){
						break
					}
					if(phase.shift.matrix[2,i+z+1]!=phase.shift.matrix[2,i+1]){
						SCa <- 1
						break
					}
					if (z==phase.shift.matrix[2,i+1]){
						z <- z+1
					}
				}
			}
			if(phase.shift.matrix[2,i+1] != phase.shift.matrix[2,2]){
				SCc <- 1
			}
			if(z>phase.shift.matrix[2,i+1] && !any(u==phase.shift.matrix[2,i+1])){
				u <- c(u,phase.shift.matrix[2,i+1])
			}
		}
		if(u != x && SCd==0){
			   SCa <- 1
		}
		x <- u
		if(SCa==1 || SCd==1){
			SCd <- SCd+1
		}
		if(SCd!=1){		#recalculate phase.shift.matrix if short indel was removed
			break
		}
	}
	#print(phase.shift.matrix)

	#resolve the remaining ambiguities.
	#resolved is used to determine the possibility of resolving an ambiguous position with the 
	#same phase shift as the preceding positions or with the same shift as 
	#the following positions
	resolved <- matrix(,nrow=4,ncol=2)

	#mark and rotate positions for long indels 
	if (is.longindel){
		SCa = phase.shift.matrix[2,2]
		for(i in 1:seq.length){
			if(phase.shift.matrix[2,i+1]>SCa){
				for(j in 0:(SCa-1)){
					phase.shift.matrix[2,i+j+1] <- SCa
					phase.shift.matrix[1,i+j+1] <- 0
				}
				i <- i+SCa
				SCa <- phase.shift.matrix[2,i+1]
			}else if(phase.shift.matrix[2,i+1]<SCa){
				SCa <- phase.shift.matrix[2,i+1]
				for(j in 1:SCa){
					phase.shift.matrix[1,i-j+1] <- 0
					phase.shift.matrix[2,i-j+1] <- SCa
				}
			}
		}
		SCa <- phase.shift.matrix[2,2]
		SCd <- 0
		for(i in 1:seq.length){
			if(phase.shift.matrix[2,i+1] != SCa){
				for(j in max(1,i-10):min(seq.length,i+10)){
					if(scores[j,SCa,2]==scores[j,phase.shift.matrix[2,i+1],2] && scores[j,SCa,2] >= scores[j,SCa,1] && 
						scores[j,phase.shift.matrix[2,i+1],2]>= scores[j,phase.shift.matrix[2,i+1],1]){
						phase.shift.matrix[1,j+1] <- 0
					}else if(scores[j,SCa,1] == scores[j,phase.shift.matrix[2,i+1],1] && 
						scores[j,SCa,2] <= scores[j,SCa,1] && 
						scores[j,phase.shift.matrix[2,i+1],2] <=scores[j,phase.shift.matrix[2,i+1],1]){
						phase.shift.matrix[1,j] <- 0
					}
					if(seq.matrix[1,j] == seq.matrix[2,j]){
						phase.shift.matrix[1,j+1] <- 1
					}
				}
				if(SCa != 0){
					SCd <- abs(SCd-1)
				}
				SCa <- phase.shift.matrix[2,i+1]
			}
			if(SCd==1){
				phase.shift.matrix[2,i+1] <- -phase.shift.matrix[2,i+1]
				if(phase.shift.matrix[1,i+1] >0 && seq.matrix[1,i] != seq.matrix[2,i]){
					phase.shift.matrix[1,i+1]  <-  3-phase.shift.matrix[1,i+1]
				}
			}
		}
	}


	if(SCc==1){
		repeat{
			phase.shift.matrix <- alignLeft(phase.shift.matrix,seq.length,seq.matrix)$matr
			phase.shift.matrix <- alignRight(phase.shift.matrix,seq.length,seq.matrix)$matr

			#print(phase.shift.matrix)
			#calculate and mark the ambiguities that could potentially be resolved




			for(i in 1:seq.length){
				if(phase.shift.matrix[1, i+1]==0){
					j <- abs(phase.shift.matrix[2,i+1])


					#print(j*SCc+i+1)
					#print(i)
					if(i<=j){
						phase.shift.matrix[3,i+1] <- 1	#1 could be resolved, 0- cannot be resolved
					}else if(j==0){
						phase.shift.matrix[3,i+1] <- 1
					}else if(phase.shift.matrix[3,i-j+1] == 0 && abs(phase.shift.matrix[2,i-j+1])==j && phase.shift.matrix[1,i-j+1]==0){
						phase.shift.matrix[3,i+1] <- 0
					}else if(phase.shift.matrix[3,i-j+1]==1 && abs(phase.shift.matrix[2,i-j+1])==j && phase.shift.matrix[1,i-j+1]==0){
						phase.shift.matrix[3,i+1] <- 1
					}else if(j*SCc+i+1 <=seq.length){
						SCc <- 1
						while(phase.shift.matrix[1,j*SCc+i+1]==0 && abs(phase.shift.matrix[2,j*SCc+1+i])==j && 
							j*SCc+i+1<=seq.length){
							SCc <- SCc+1
						}
						#cat("erste Stelle",3-phase.shift.matrix[1,j*SCc+i+1],"\n")
						#cat("zweite Stelle",j*SCc+i,"\n")
						#cat("Wert",seq.matrix[3-phase.shift.matrix[1,j*SCc+i+1],j*SCc+i],"\n")


						if(abs(phase.shift.matrix[2,i-j+1]) < j && phase.shift.matrix[2,i+1] >0){
							phase.shift.matrix[3,i+1] <- 1
						}else if((j*SCc +i) > seq.length){ # evtl +1
							phase.shift.matrix[3, i+1] <- 1
						}else if(abs(phase.shift.matrix[2,j*SCc+i+1]) != j){
							phase.shift.matrix[3,i+1] <- 1
						}else if(abs(phase.shift.matrix[2,j*SCc+i+1+phase.shift.matrix[4,i+1]]) != j){
							phase.shift.matrix[3,i+1] <- 1


						}else if(phase.shift.matrix[1,i-j+1]>0 && (3-phase.shift.matrix[1,j*SCc+i+1]<3 )&& 
							(3-phase.shift.matrix[1,j*SCc+i+1])>0){

							if((phase.shift.matrix[2,i+1] > 0) && ((SCc / 2 ) != floor(SCc/2)) && 
							((seq.matrix[phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] && 
									seq.matrix[1,i] == seq.matrix[3-phase.shift.matrix[1,j*SCc+i+1],j*SCc+i]) || 
								(seq.matrix[phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i] && 
									seq.matrix[2,i]==seq.matrix[3-phase.shift.matrix[1,j*SCc+i+1],j*SCc+i]))){
								

								phase.shift.matrix[3,i+1] <- 1
							}
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*SCc+i+1])<3 && 
							3-phase.shift.matrix[1,j*SCc+i+1]>0){

							if(phase.shift.matrix[2,i+1] >0 && (SCc/2) == floor(SCc/2) && 
								((seq.matrix[phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i] && 
									seq.matrix[2,i] == seq.matrix[3-phase.shift.matrix[1,j*SCc+i+1],j*SCc+i])||
								(seq.matrix[phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[1,i] && 
									seq.matrix[1,i] == seq.matrix[3- phase.shift.matrix[1,j*SCc+i+1],j*SCc+i]))){
								phase.shift.matrix[3,i+1] <- 1
							}
						}else if(abs(phase.shift.matrix[2,i-j+1]) != j && phase.shift.matrix[2,i+1] < 0){
							phase.shift.matrix[3,i+1] <- 1
						}else if(phase.shift.matrix[1,i-j+1] ==0){
							phase.shift.matrix[3,i+1] <- 0
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*SCc+i+1])<3 && 
							3-phase.shift.matrix[1,j*SCc+i+1]>0){

							if(phase.shift.matrix[2,i+1] < 0 && (SCc/2) != floor(SCc/2) && 
								((seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[1,i] && 
									seq.matrix[2,i]== seq.matrix[phase.shift.matrix[1,j*SCc+i+1],j*SCc+i]) || 
								(seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] && 
									seq.matrix[1,i]==seq.matrix[phase.shift.matrix[1,j*SCc+i+1],j*SCc+i]))){
								phase.shift.matrix[3,i+1] <- 1

						}
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*SCc+i+1])<3 && 
							3-phase.shift.matrix[1,j*SCc+i+1]>0){

							if(phase.shift.matrix[2,i+1] < 0 && (SCc/2) == floor(SCc/2) && 
								(( seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i-j] && 
									seq.matrix[1,i] == seq.matrix[phase.shift.matrix[1,j*SCc+i+1],j*SCc+i]) || 
								(seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] &&
									seq.matrix[2,i]==seq.matrix[phase.shift.matrix[1,j*SCc+i+1],j*SCc+i]))){
								phase.shift.matrix[3,i+1] <- 1
							}
						}else{

							phase.shift.matrix[3,i+1] <- 0
						}
					}
				}
			}

			j <- 0

			for(i in 1:seq.length){
				

				if(phase.shift.matrix[1,i+1]==0 && phase.shift.matrix[3,i+1]==1){

					for(z in 1:(max.shift+phase.shift.matrix[4,1])){
						if(z>=i || z> seq.length -i){
							break
						}
						if(phase.shift.matrix[2,i-z+1] != phase.shift.matrix[2,i+z+1]){
							break
						}
					}

					if(z<=i && z-phase.shift.matrix[3,i+1] <= max.shift+1 && z <= seq.length -(i+1)){

						ind1 <- abs(phase.shift.matrix[2,i-z+1])
						ind2 <- abs(phase.shift.matrix[2,i+z+1])
						if(!is.longindel){
							if((i> ind1) && (i>ind2) && (seq.length -i)>ind1 && 
								(seq.length-i)>ind2){
								if(((seq.matrix[2,i]==seq.matrix[1,i-ind1] && 
									scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2]) || 
								(seq.matrix[2,i] == seq.matrix[2,i-ind1] && 
									scores[i-ind1+1,ind1+1,1] <= scores[i-ind1+1,ind1+1,2])) && 
								((seq.matrix[1,i] == seq.matrix[2,i+ind1] && 
									scores[i+ ind1+1,ind1+1,1] >= scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[1,i] == seq.matrix[1,i+ind1] && 
									scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
									resolved[1,1] <- 1
								}else{
									resolved[1,1] <- 0
								}
								if(((seq.matrix[1,i]==seq.matrix[1,i-ind1] && 
									scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
								(seq.matrix[1,i] == seq.matrix[2,i-1] &&
									scores[i-ind1+1,ind1+1,1] <= scores[i-ind1+1,ind1+1,2]))&&
								((seq.matrix[2,i]==seq.matrix[2,i+ind1] && 
									scores[i+ind1+1,ind1+1,1]>=scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[1,i+ind1] && 
									scores[i+ind1+1,ind1+1,1]<=scores[i+ind1+1,ind1+1,2]))){
									resolved[1,2] <- 1
								}else{
									resolved[1,2] <- 0
								}


								if(((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]>=scores[i-ind2+1,ind2+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]<= scores[i-ind2+1,ind2+1,2]))&&
								((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
								(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[2,1] <- 1
								}else{
									resolved[2,1] <- 0
								}


								if(((seq.matrix[1,i]==seq.matrix[1,i-ind2] && 
									scores[i-ind2+1,ind2+1,1]>=scores[i-ind2+1,ind2+1,2])||
								(seq.matrix[1,i]==seq.matrix[2,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]<=scores[i-ind2+1,ind2+1,2]))&&
								((seq.matrix[2,i]==seq.matrix[2,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
								(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[2,2] <- 1
								}else{
									resolved[2,2] <- 0
								}
								if((ind1 < ind2) &&
									((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&
										scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
									(seq.matrix[2,i]==seq.matrix[2,i-ind1]&&
										scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2]))&&
									((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&
										scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
									(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&
										scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[4,1] <- 1
								}else{
									resolved[4,1] <- 0
								}
								if((ind1<ind2)&&
									((seq.matrix[1,i]==seq.matrix[1,i-ind1]&&
										scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
									(seq.matrix[1,i]==seq.matrix[2,i-ind1]&&
										scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2]))&&
									((seq.matrix[2,i]==seq.matrix[2,i+ind2]&&
										scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
									(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&
										scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[4,2] <- 1
								}else{
									resolved[4,1] <- 0
								}
							}else if(i<=ind1){
								if(((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&
									scores[i+ind1+1,ind1+1,1]>=scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&
									scores[i+ind1+1,ind1+1,1]<=scores[i+ind1+1,ind1+1,2]))){
									resolved[1,1] <- 1
								}else{
									resolved[1,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&
									scores[i+ind1+1,ind1+1,1]>=scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&
									scores[i+ind1+1,ind1+1,1]<=scores[i+ind1+1,ind1+1,2]))){
									resolved[1,2] <- 1
								}else{
									resolved[1,2] <- 0
								}
								resolved[2,] <- 0
							}else if(seq.length-i< ind2){
								resolved[1,] <- 0
								if(((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]>=scores[i-ind2+1,ind2+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]<=scores[i.ind2+1,ind2+1,2]))){
									resolved[2,1] <- 1
								}else{
									resolved[2,1] <- 0
								}
								if(((seq.matrix[1,i]==seq.matrix[1,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]>=scores[i-ind2+1,ind2+1,2])||
								(seq.matrix[1,i]==seq.matrix[2,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]<=scores[i-ind2+1,ind2+1,2]))){
									resolved[2,2] <- 1
								}else{
									resolved[2,2] <- 0
								}
							}							
							if(phase.shift.matrix[2,i-z+1]<phase.shift.matrix[2,i+z+1]||i<=ind1){
								
								if(((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
								(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[3,1] <- 1
								}else{
									resolved[3,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
								(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[3,2] <- 1
								}else{
									resolved[3,2] <- 0
								}
							}else{

								if(ind1==0||((seq.matrix[2,i]==seq.matrix[1,i-ind1]&&
									scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind1]&&
									scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2]))){
									resolved[3,1] <- 1
								}else{
									resolved[3,1] <- 0
								}
								if(ind1==0||((seq.matrix[1,i]==seq.matrix[1,i-ind1]&&
									scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind1]&&
									scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2]))){
									resolved[3,2] <- 1
								}else{
									resolved[3,2] <- 0
								}
							}
							if(ind1==0||ind2==0){
								resolved[4,] <- 0
							}	

						}else{  #long indel
							if(phase.shift.matrix[2,i-z+1]>0){
								if(((seq.matrix[2,i]==seq.matrix[1,i-ind1]&&
									phase.shift.matrix[1,i-ind1+1] !=2)||
								(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind1+1]!=1))){
									resolved[3,1] <- 1
								}else{
									resolved[3,1] <- 0
								}
								if(((seq.matrix[1,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[1,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1]!=1))){
									resolved[3,2] <- 1
								}else{
									resolved[3,2] <- 0
								}
							}else{
								if(((seq.matrix[1,i]==seq.matrix[2,i-ind1] && phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[1,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))){
									resolved[4,1] <- 1
								}else{
									resolved[4,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[2,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))){
									resolved[4,2] <- 1
								}else{
									resolved[4,2] <- 0
								}
							}
							if(phase.shift.matrix[2,i+z+1]>0){
								if(((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
									(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+2]!=1))){
									resolved[4,1] <- 1
								}else{
									resolved[4,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
									(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
									resolved[4,2] <- 1
								}else{
									resolved[4,2] <- 0
								}
							}else{
								if(((seq.matrix[2,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
									(seq.matrix[2,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
									resolved[4,1] <- 1
								}else{
									resolved[4,1] <- 0
								}
								if(((seq.matrix[1,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
									(seq.matrix[1,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
									resolved[4,2] <- 1
								}else{
									resolved[4,2] <- 0
								}
							}
							if(i >= ind1 && i>= ind2 && (seq.length-i)>ind1&&
								(seq.length-i)<ind2){
								if(phase.shift.matrix[2,i-z+1]>0){
									if(((seq.matrix[2,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))&&
									((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[1,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1]!=1))&&
									((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}else{
									if(((seq.matrix[1,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))&&
									((seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))&&
									((seq.matrix[1,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}
								if(phase.shift.matrix[2,i+z+1]>0){
									if(((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=1))&&
									((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
										resolved[2,1] <- 1
									}else{
										resolved[2,1] <- 0
									}
									if(((seq.matrix[1,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=1))&&
									((seq.matrix[2,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
										resolved[2,2] <- 1
									}else{
										resolved[2,2] <- 0
									}
								}else{
									if(((seq.matrix[1,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=1))&&
									((seq.matrix[2,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2]!=2)||
										(seq.matrix[2,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
										resolved[2,1] <- 1
									}else{
										resolved[2,1] <- 0
									}
									if(((seq.matrix[2,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind2]!=1))&&
									((seq.matrix[1,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
										resolved[2,2] <- 1
									}else{
										resolved[2,2] <- 0
									}
								}
							}else if(i<=ind1){
									if(((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
									resolved[2,] <- 0
							}else if((seq.length-i)<ind2){
									resolved[1,] <- 0
									if(((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind1]!=2)||
										(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind1]!=1))){
										resolved[2,1] <- 1
									}else{
										resolved[2,1] <- 0
									}
									if(((seq.matrix[1,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind1]!=1))){
										resolved[2,2] <- 1
									}else{
										resolved[2,2] <- 0
									}
							}
						}

						if(ind1==0){
							resolved[1,] <- 0
						}
						if(ind2==0){
							resolved[2,] <- 0
						}
						if(!is.longindel){
							if(ind1>ind2){
								if(phase.shift.matrix[2,i]==ind1){
									if(z<=ind1){
										resolved[1,] <- 0
										if(z<=ind2){
											resolved[3,] <- 0
										}
									}

									if((phase.shift.matrix[4,i+1]-z)<0){
										resolved[2,] <- 0
									}
								}else{
									resolved[1,] <- 0
									resolved[3,] <- 0
								}
							}else{ #ind1<ind2
								if(phase.shift.matrix[2,i+1]==ind1){
									if(phase.shift.matrix[3,i+1]<=ind2){
										resolved[2,] <- 0
										if(ind1 != 0 && z > phase.shift.matrix[4,i+1]){
											resolved[3,] <- 0
										}
									}
									if(z>ind2){
										resolved[2,] <- 0
									}
								}else{
									if(z+phase.shift.matrix[4,i+1]<=ind2){
										resolved[2,] <- 0
									}
									if(z>ind1){
										resolved[4,] <- 0
									}
									resolved[1,] <- 0
								}
							}
						}else{
							if(ind1>ind2){
								if(abs(phase.shift.matrix[2,i+1])==ind1){
									if(z<=ind1){
										resolved[1,] <- 0
									}
									if((phase.shift.matrix[4,i+1]-z)<0){
										resolved[c(2,4),] <- 0
									}
								}else{
									resolved[c(1,3),] <- 0
								}
							}else{
								if(phase.shift.matrix[2,i+1] <= ind2){
									if(phase.shift.matrix[4,i+1]<=ind2){
										resolved[2,] <- 0
									}
									if(ind1!=0 &&z > phase.shift.matrix[4,i+1]){
										resolved[4,] <- 0
									}
									if(z>ind2){
										resolved[2,] <- 0
									}
								}else{
									if((z+phase.shift.matrix[4,i+1])<=ind2){
										resolved[2,] <- 0
									}
									resolved[c(1,3),] <- 0
								}
							}
						}
						if(!is.longindel){
							if(all(resolved[1,]==1)){
							}else if(all(resolved[2,]==1)){
							}else if(all(resolved[c(1,2),1]==1)&&all(resolved[c(1,2),2]==0)){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[c(1,2),1]==0)&&all(resolved[c(1,2),2]==1)){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}else if((resolved[1,1]==1&&resolved[2,2]==1)||(resolved[1,2]==1&&resolved[2,1]==1)){
							}else if(all(resolved[c(1,3),1]==1)&&all(resolved[c(1,3),2]==0)){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[c(1,3),1]==0)&&all(resolved[c(1,3),2]==1)){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}else if(((resolved[1,1]==1&&resolved[3,2]==1)||(resolved[1,2]==1&&resolved[3,1]==1))
								&&z< ind1 &&phase.shift.matrix[4,i+1]>ind1){
							}else if(((resolved[2,1]==1&&resolved[3,2]==1)||(resolved[2,2]==1&&resolved[3,1]==1))
								&&z<=ind1 &&phase.shift.matrix[4,i+1]>ind2 &&ind2 >ind1){
							}else if(((resolved[2,1]==1&&resolved[3,2]==1)||(resolved[2,2]==1&&resolved[3,1]==1))
								&&phase.shift.matrix[2,i+1]==ind2&&phase.shift.matrix[1,i+ind2+1]==0){
							}else if(((resolved[2,1]==1&&resolved[3,2]==1)||(resolved[2,2]==1&&resolved[3,1]==1))
								&&phase.shift.matrix[2,i+1]==ind1&&phase.shift.matrix[4,i+1]>=z){
							}else if(((resolved[1,1]==1&&resolved[3,2]==1)||(resolved[1,2]==1&&resolved[3,1]==1))
								&&phase.shift.matrix[2,i+1]==ind1&&ind2>ind1&&phase.shift.matrix[4,i+1]>=(z+ind1)){
							}else if(xor(resolved[1,1]==1,resolved[1,2]==1)){
								if(resolved[1,1]==1){
									phase.shift.matrix[1,i+1] <- 1
									scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
								}else{
									phase.shift.matrix[1,i+1] <- 2
									scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
								}
							}else if(xor(resolved[2,1]==1,resolved[2,2]==1)){
								if(resolved[2,1]==1){
									phase.shift.matrix[1,i+1] <- 1
									scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
								}else{
									phase.shift.matrix[1,i+1] <- 2
									scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
								}
							}else if(all(resolved[3,]==1)&&all(resolved[4,]==0)){
							}else if(all(resolved[3,]==1)&&all(resolved[4,]==c(1,0))){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[3,]==1)&&all(resolved[4,]==c(0,1))){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}else if(resolved[3,1]==1){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(resolved[3,2]==1){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}
						}else{

							if(all(resolved[1,]==1)){
							}else if(all(resolved[2,]==1)){
							}else if(((resolved[1,1]==1&&resolved[3,2]==1)||(resolved[1,2]==1&&resolved[3,1]==1))){
							}else if(((resolved[2,1]==1&&resolved[4,2]==1)||(resolved[2,2]==1&&resolved[4,1]==1))){
							}else if(((resolved[1,1]==1&&resolved[4,2]==1)||(resolved[1,2]==1&&resolved[4,1]==1))){
							}else if(((resolved[2,1]==1&&resolved[3,2]==1)||(resolved[2,2]==1&&resolved[3,1]==1))){
							}else if(((resolved[4,1]==1&&resolved[3,2]==1)||(resolved[4,2]==1&&resolved[3,1]==1))){
							}else if(xor(resolved[1,1]==1,resolved[1,2]==1)){
								if(resolved[1,1]==1){
									phase.shift.matrix[1,i+1] <- 1
									scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
								}else{
									phase.shift.matrix[1,i+1] <- 2
									scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
								}
							}else if(xor(resolved[2,1]==1,resolved[2,2]==1)){
								if(resolved[2,1]==1){
									phase.shift.matrix[1,i+1] <- 1
									scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
								}else{
									phase.shift.matrix[1,i+1] <- 2
									scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
								}
							}else if(all(resolved[3,]==1)){
							}else if(all(resolved[4,]==1)){
							}else if(all(resolved[4,]==c(1,0))){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[4,]==c(0,1))){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[3,]==c(1,0))){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[3,]==c(0,1))){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}
						}
					}else{
						if(i==1||i==seq.length||phase.shift.matrix[2,i]==phase.shift.matrix[2,i+2]){
							if(i!=seq.length){
								ind1 <- abs(phase.shift.matrix[2,i+2])
								ind2 <- phase.shift.matrix[2,i+2]
							}else{
								ind1 <- abs(phase.shift.matrix[2,i])
								ind2 <- phase.shift.matrix[2,i]
							}
							if(i<=ind1){
								if(((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&
									scores[i+ind1+1,ind1+1,1] >= scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&
									scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
									resolved[1,1] <- 1
								}else{
									resolved[1,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&
									scores[i+ind1+1,ind1+1,1] >= scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&
									scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
									resolved[1,2] <- 1
								}else{
									resolved[1,2] <- 0
								}
							}else if(seq.length - i < ind1){
								if(ind2 > 0){
									if(((seq.matrix[2,i] == seq.matrix[1,i - ind1] &&
										scores[i - ind1 + 1,ind1+1,1] >= scores[i - ind1 + 1,ind1+1,2])||
									(seq.matrix[2,i] == seq.matrix[2,i - ind1] &&
										scores[i - ind1 + 1,ind1+1,1] <= scores[i - ind1 + 1,ind1+1,2]))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[1,i] == seq.matrix[1,i - ind1] &&
										scores[i - ind1 + 1,ind1+1,1] >= scores[i - ind1 + 1,ind1+1,2]) ||
									(seq.matrix[1,i] == seq.matrix[2,i - ind1] &&
										scores[i - ind1 + 1,ind1+1,1] <= scores[i - ind1 + 1,ind1+1,2]))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}else{
									if(((seq.matrix[1,i] == seq.matrix[2,i-ind1] &&
										phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[1,i] == seq.matrix[1,i-ind1] &&
										phase.shift.matrix[1,i-ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i] == seq.matrix[2,i-ind1] &&
										phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[2,i] == seq.matrix[1,i-ind1] &&
										phase.shift.matrix[1,i-ind1+1]!=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}
							}else{
								if(ind2 > 0){
									if(((seq.matrix[2,i] == seq.matrix[1,i-ind1] && 
										scores[i-ind1+1,ind1+1,1] >= scores[i-ind1+1,ind1+1,2])||
									(seq.matrix[2,i] == seq.matrix[2,i-ind1]&&
										scores[i-ind1+1,ind1+1,1]<= scores[i-ind1+1,ind1+1,2])) &&
									((seq.matrix[1,i] == seq.matrix[2,i+ind1] &&
										scores[i+ind1+1,ind1+1,1]>=scores[i+ind1+1,ind1+1,2]) ||
									(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&
										scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[1,i] == seq.matrix[1,i-ind1] && 
										scores[i-ind1+1,ind1+1,1] >= scores[i-ind1+1,ind1+1,2])||
									(seq.matrix[1,i]==seq.matrix[2,i-ind1] && 
										scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2])) &&
									((seq.matrix[2,i]==seq.matrix[2,i+ind1] &&
										scores[i+ind1+1,ind1+1,1] >= scores[i+ind1+1,ind1+1,2])||
									(seq.matrix[2,i]==seq.matrix[1,i+ind1] &&
										scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}else{
									if(((seq.matrix[1,i] == seq.matrix[2,i-ind1] &&phase.shift.matrix[1,i-ind1+1] !=2) || 
									(seq.matrix[1,i] == seq.matrix[1,i-ind1] && phase.shift.matrix[1,i-ind1+1] !=1)) &&
									((seq.matrix[2,i] == seq.matrix[1,i+ind1] && phase.shift.matrix[1,i+ind1+1] != 2) || 
									(seq.matrix[2,i] == seq.matrix[2,i+ind1] && phase.shift.matrix[1,i+ind1+1] !=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i] == seq.matrix[2,i-ind1] && phase.shift.matrix[1,i-ind1+1] !=2) ||
										(seq.matrix[2,i] == seq.matrix[1,i-ind1] && phase.shift.matrix[1,i-ind1+1]!=1)) &&
									((seq.matrix[1,i] == seq.matrix[1,i+ind1] && phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[1,i] == seq.matrix[2,i+ind1] && phase.shift.matrix[1,i+ind1+1] !=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}
							}
							if(resolved[1,1] ==1 &&resolved[1,2] ==0){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(resolved[1,1] == 0 && resolved[1,2] ==1){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}
						}
					}
					if(phase.shift.matrix[1,i+1] != 0){
						j <- 1
					}
				}
			}


			if(j!=1){
				break
			}
		}

	}


	if(!is.right&&is.align){
		phase.shift.matrix<- alignLeft(phase.shift.matrix,seq.length,seq.matrix)$matr
		phase.shift.matrix <- alignRight(phase.shift.matrix,seq.length,seq.matrix)$matr
	}else if(is.right&&is.align){
		phase.shift.matrix <- alignRight(phase.shift.matrix,seq.length,seq.matrix)$matr
		phase.shift.matrix <- alignLeft(phase.shift.matrix,seq.length,seq.matrix)$matr
	}

	#resolve 3-fold degenerate bases

	if(bd){
		for(i in 1:seq.length){
			if(any(seq.matrix[1,i] == c("B","D","H","V","N"))){
				if(i>phase.shift.matrix[2,i+1] && i<=seq.length-phase.shift.matrix[2,i+1]){
					if(phase.shift.matrix[1,i-phase.shift.matrix[2,i+1]+1]>0 && phase.shift.matrix[1,i+phase.shift.matrix[2,i+1]+1]>0){
						if(any(seq.matrix[phase.shift.matrix[1,i-phase.shift.matrix[2,i+1]+1],i-phase.shift.matrix[2,i+1]] ==c("A","c","G","T"))){
							if(any(seq.matrix[3-phase.shift.matrix[1,i+phase.shift.matrix[2,i+1]+1],i+phase.shift.matrix[2,i+1]]==c("A","C","G","T"))){
								if(any(seq.matrix[1,i] ==c("B","D","H","V") & 
									seq.matrix[phase.shift.matrix[1,i-phase.shift.matrix[2,i+1]+1],i-phase.shift.matrix[2,i+1]]!=c("A","C","G","T")
								& seq.matrix[ 3 - phase.shift.matrix[1,i + phase.shift.matrix[2,i+1]+1],i+ phase.shift.matrix[2,i]]!=c("A","C","G","T") &
								seq.matrix[phase.shift.matrix[1,i - phase.shift.matrix[2,i+1]+1],i - phase.shift.matrix[2,i+1]] !=
								seq.matrix[ 3 - phase.shift.matrix[1,i + phase.shift.matrix[2,i+1]+1],i + phase.shift.matrix[2,i+1]])){
									seq.matrix[1,i] <- seq.matrix[3 - phase.shift.matrix[1,i + phase.shift.matrix[2,i+1]+1],i + phase.shift.matrix[2,i+1]]
									seq.matrix[2,i] <- seq.matrix[phase.shift.matrix[1,i - phase.shift.matrix[2,i+1]+1],i - phase.shift.matrix[2,i+1]]
								}else if(seq.matrix[1,i]=="N"&&seq.matrix[phase.shift.matrix[1,i-phase.shift.matrix[2,i+1]+1],i-phase.shift.matrix[2,i+1]+1]!=
									seq.matrix[3-phase.shift.matrix[1,i+phase.shift.matrix[2,i+1]+1],i+phase.shift.matrix[2,i+1]]){
									seq.matrix[1,i] <- seq.matrix[3-phase.shift.matrix[1,i+phase.shift.matrix[2,i+1]+1],i+phase.shift.matrix[2,i+1]]
									seq.matrix[2,i] <- seq.matrix[phase.shift.matrix[1,i-phase.shift.matrix[2,i+1]+1],i-phase.shift.matrix[2,i+1]]
								}
							}
						}
					}
				}
			}
		}
	}

	#based on phase.shift.matrix and seq.matrix reconstruct two allelic sequences
	temp1 <- c()
	temp2 <- c()
	temp3 <- c()
	ambig2 <- 0
	for(i in 1:seq.length){
		if(phase.shift.matrix[1,i+1]==1){
			temp1 <- c(temp1,seq.matrix[1,i])
			temp2 <- c(temp2,seq.matrix[2,i])
		}else if(phase.shift.matrix[1,i+1] == 2){
			temp1 <- c(temp1, seq.matrix[2,i])
			temp2 <- c(temp2, seq.matrix[1,i])
		}else{
			ambig2 <- ambig2 +1
			if(seq.matrix[1,i]=="A" && seq.matrix[2,i]=="G"){
				temp1 <- c(temp1,"R")
				temp2 <- c(temp2,"R")
			}else if(seq.matrix[1,i]=="C" && seq.matrix[2,i]=="T"){
				temp1 <- c(temp1,"Y")
				temp2 <- c(temp2,"Y")
			}else if(seq.matrix[1,i]=="G" && seq.matrix[2,i]=="C"){
				temp1 <- c(temp1,"S")
				temp2 <- c(temp2,"S")
			}else if(seq.matrix[1,i]=="A" && seq.matrix[2,i]=="T"){
				temp1 <- c(temp1,"W")
				temp2 <- c(temp2,"W")
			}else if(seq.matrix[1,i]=="A" && seq.matrix[2,i]=="C"){
				temp1 <- c(temp1,"M")
				temp2 <- c(temp2,"M")
			}else if(seq.matrix[1,i]=="G" && seq.matrix[2,i]=="T"){
				temp1 <- c(temp1,"K")
				temp2 <- c(temp2,"K")
			}
		}
		temp3 <- c(temp3,phase.shift.matrix[2,i+1])
	}


	#Align reconstructed allelic sequences (indicate gaps with dots)

	SC <- 0
	temp3 <- c()
	temp4 <- c()
	x <- c()
	for(i in 1:seq.length){
		if(phase.shift.matrix[2,i+1]>SC){
			dots <- character(phase.shift.matrix[2,i+1]-SC)
			dots[] <- "."
			temp3 <- c(temp3,dots)
			SC <- phase.shift.matrix[2,i+1]
		}else if(phase.shift.matrix[2,i+1]<SC){
			dots <- character(SC-phase.shift.matrix[2,i+1])
			dots[] <- "."
			temp4 <- c(temp4,dots)
			SC <- phase.shift.matrix[2,i+1]
		}
		temp3 <- c(temp3,temp1[i])
		temp4 <- c(temp4,temp2[i])
		x <- c(x,abs(phase.shift.matrix[2,i+1]))
	}
	if(SC>0){
		dots <- character(SC)
		dots[] <- "."
		temp4 <-c(temp4,dots)
	}else{
		dots <- character(abs(SC))
		dots[]<- "."
		temp3 <- c(temp3,dots)
	}

	#Outout reconstructed sequences with mismatches highlighted by red color

	j <- length(temp4)
	temp1 <- ""
	temp2 <- ""
	temp5 <- ""
	SCMax <- 0
	SCd <- 0
	mismatches <- 0
	ambiguities <- c("R","Y","K","M","S","W","B","D","H","V","N")
	ambiguities.resolved <- c("AG","CT","GT","AC","CG","AT","CGT","AGT","ACT","ACG","ACGT")
	for(i in 1:j){
		SCa <- temp3[i]
		SCb <- temp4[i]
		SCc <- ""
		if(!is.na(SCa)&&!is.na(SCb)){
			if(SCa=="." && SCb=="."){
			}else if(SCa==SCb && (any(SCa==c("A","C","G","T")))){
				temp1 <- temp1 %+% SCa
				temp2 <- temp2 %+% SCb
			}else{
				if(SCa=="."){
					temp1 <- temp1 %+% SCa
				}else{
					temp1 <- temp1 %+% red(SCa)
				}
				if(SCb=="."){
					temp2 <- temp2 %+% SCb
				}else{
					temp2 <- temp2 %+% red(SCb)
				}
			}

			if(SCa == "." &&SCb=="."){
			}else if(SCa=="."){
				temp5 <- temp5 %+% red(SCb)
				if(all(SCb!=c("A","C","G","T"))){
					SCd <- SCd +1
				}
			}else if(SCb=="."){
				temp5 <- temp5 %+% red(SCa)
				if(all(SCa!=c("A","C","G","T"))){
					SCd <- SCd + 1
				}
			}else if(SCa==SCb&&any(SCa==c("A","C","G","T"))){
				temp5 <- temp5 %+% SCa
				SCMax <- SCMax + 1
			}else{
				SCd <- SCd +1 
				if(any(SCa==ambiguities)){
					SCa <- ambiguities.resolved[SCa==ambiguities]
				}
				if(any(SCa==ambiguities)){
					SCa <- ambiguities.resolved[SCa==ambiguities]
				}
				if((countCharOccurrences("A",SCa) >0 && countCharOccurrences("A",SCb) >0) ||
					(countCharOccurrences("C",SCa) >0 && countCharOccurrences("C",SCb) >0) ||
					(countCharOccurrences("G",SCa) >0 && countCharOccurrences("G",SCb) >0) ||
					(countCharOccurrences("T",SCa) >0 && countCharOccurrences("T",SCb) >0)){
					SCMax <- SCMax +1 
				}else{
					mismatches <- mismatches +1
				}
				if(countCharOccurrences("A",SCa) >0||countCharOccurrences("A",SCb) >0){
					SCc <- SCc %+% "A"
				}
				if(countCharOccurrences("C",SCa) >0||countCharOccurrences("C",SCb) >0){
					SCc <- SCc %+% "C"
				}
				if(countCharOccurrences("G",SCa) >0||countCharOccurrences("G",SCb) >0){
					SCc <- SCc %+% "G"
				}
				if(countCharOccurrences("T",SCa) >0||countCharOccurrences("T",SCb) >0){
					SCc <- SCc %+% "T"
				}
				if(countCharOccurrences("G",SCc) >0 && countCharOccurrences("C",SCc) >0 && 
				countCharOccurrences("T",SCc) >0 && countCharOccurrences("A",SCc) >0){
					SCc <- "N"
				}else if(countCharOccurrences("G",SCc) >0&&countCharOccurrences("C",SCc) >0&&
					countCharOccurrences("T",SCc) >0){
					SCc <- "B"
				}else if(countCharOccurrences("G",SCc) >0 && countCharOccurrences("A",SCc) >0 &&
					countCharOccurrences("T",SCc) >0){
					SCc <- "D"
				}else if(countCharOccurrences("A",SCc) >0 && countCharOccurrences("C",SCc) >0 &&
					countCharOccurrences("T",SCc) >0){
					SCc <- "H"
				}else if(countCharOccurrences("G",SCc) >0 && countCharOccurrences("C",SCc) >0 &&
					countCharOccurrences("A",SCc) >0){
					SCc <- "V"
				}else if(countCharOccurrences("A",SCc) >0 && countCharOccurrences("G",SCc) >0){
					SCc <- "R"
				}else if(countCharOccurrences("C",SCc) >0 && countCharOccurrences("T",SCc) >0){
					SCc <- "Y"
				}else if(countCharOccurrences("G",SCc) >0 && countCharOccurrences("T",SCc) >0){
					SCc <- "K"
				}else if(countCharOccurrences("A",SCc) >0 && countCharOccurrences("C",SCc) >0){
					SCc <- "M"
				}else if(countCharOccurrences("C",SCc) >0 && countCharOccurrences("G",SCc) >0){
					SCc <- "S"
				}else if(countCharOccurrences("A",SCc) >0 && countCharOccurrences("T",SCc) >0){
					SCc <- "W"
				}
			}
			temp5 <- temp5 %+% red(SCc)
		}
	}

	cat("Reconstructed Sequences:\n",temp1,"\n",temp2,"\nCombined Sequence:\n",temp5,"\n")

}
		




		       
	       			       




