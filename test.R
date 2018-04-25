suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))

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
		[default \"%default\"]")
	)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#set variables
max.shift <- opt$max
penalty <- opt$penalty
if(!is.null(opt$fix)){
	fixed.shifts <- as.numeric(unlist(strsplit(opt$fix, split=",")))
	}else{fixed.shifts <- NULL}
is.right <- opt$right
is.align <- opt$align
is.longindel <- opt$long



ambig <- 0

seq.number <- opt$number



#check if input or string is set
if(!is.null(opt$input)){
	#parse input file
	sequence <- trimws(readLines(opt$input))

	#get combined sequence
	sequence <- sequence[grep("Combined",sequence)+1]

	if(length(sequence)<seq.number){
		print_help(opt_parser)
		stop("File does not contain ",seq.number," sequences, but only ",
		length(sequence),".",call.=FALSE)
	}else{
		sequence <- sequence[seq.number]
	}
}else if (!is.null(opt$string)){
	sequence <- opt$string
}else{
	print_help(opt_parser)
	stop("Please provide an input file or string.\n",call.=FALSE)
}	



#converte input string to matrix of sequences, where 2-fold ambiguities are split 
seq.toMatrix <- function(sequence,max.shift){

	#convert sequence to all uppercase
	sequence <- toupper(sequence)

	#split sequence letter by letter
	seq.asvector <- strsplit(sequence,"")[[1]]
	
	#set sequence length
	seq.length <- nchar(sequence)
	
	#maximum shift is limited to half of sequence length
	if(max.shift*2>seq.length){
		max.shift <- floor(nchar(sequence)/2)
	}

	#initialize matrix seq.matrix
	seq.matrix <- matrix(,nrow=2,ncol=seq.length)


	if (seq.length!=0){ #if sequence exists
	

	
	

	pos.of <- seq.asvector=="R"
	seq.matrix[1,pos.of] <- "A"
	seq.matrix[2,pos.of] <- "G"
	pos.of <- seq.asvector=="Y"
	seq.matrix[1,pos.of] <- "C"
	seq.matrix[2,pos.of] <- "T"
	pos.of <- seq.asvector=="S"
	seq.matrix[1,pos.of] <- "G"
	seq.matrix[2,pos.of] <- "C"
	pos.of <- seq.asvector=="W"
	seq.matrix[1,pos.of] <- "A"
	seq.matrix[2,pos.of] <- "T"
	pos.of <- seq.asvector=="M"
	seq.matrix[1,pos.of] <- "A"
	seq.matrix[2,pos.of] <- "C"
	pos.of <- seq.asvector=="K"
	seq.matrix[1,pos.of] <- "G"
	seq.matrix[2,pos.of] <- "T"
	#fill remaining characters
	pos.of <- is.na(seq.matrix[1,])
	seq.matrix[1,pos.of] <- seq.asvector[pos.of]
	seq.matrix[2,pos.of] <- seq.asvector[pos.of]
	#check for 3-fold degenerate bases

	ambig1 <- 0
	if (!(seq.matrix=="A" || seq.matrix=="T" || seq.matrix=="G" || seq.matrix=="C")){
		bd <- TRUE
		ambig1 <- sum(!(seq.matrix=="A" | seq.matrix=="T" | seq.matrix=="G" | seq.matrix=="C"))
	}else{
		bd <- FALSE
	}


	

	return(list(matrix=seq.matrix,is.three.fold=bd,ambig=ambig1,shift=max.shift))
	}
}


#calculates the forward score of the sequence
calculate.forward <- function(seq.matrix,max.shift,penalty){
	seq.length <- length(seq.matrix[1,])
	#set scores as in indelligent
	match.score <- 1
	mismatch.score <- 0

	#init forward matrix
	#forward saves score for forward[i,j,k]. where i is the position in the sequence, 
	#j is the shift and k is the interpretation of ambiguity

	forward <- array(0,dim=c(seq.length+2,max.shift+1,2))

	#set starting conditions for forward calculation
	for ( i in 2:(max.shift+1)){
		forward[i,i:(max.shift+1),] <- (i-1)
	}

	#calculate score in forward direction

	#iterate over seqence length
	for ( i in 1:seq.length){
		jMax <- max.shift

		#check if iteration is smaller than possible shift
		if ((i-1) < max.shift){	
			jMax <- i-1
		}

		#iterate over possible shifts
		for (j in 1:(jMax+1)){
			for (z in 1:2){
				SCMax <- -111111

				#calculate score for eveery shift
				for (x in 1:(jMax+1)){
					if (j==1 && x==1 && seq.matrix[1,i] != seq.matrix[2,i]){
						
						SC <- max(forward[i,x,1],forward[i,x,2]) - mismatch.score
					}else if (j==1 && seq.matrix[1,i]==seq.matrix[2,i]){
						SC <- max(forward[i,x,1],forward[i,x,2]) + match.score
					}else if (j != x){
						SC <- max(forward[i,x,1],forward[i,x,2]) + match.score
					}else if (((seq.matrix[3-z,i]==seq.matrix[1,i-(j-1)]) && 
					(forward[i-(j-2),j,1] >= forward[i-(j-2),j,2])) || 
					((seq.matrix[3-z,i]==seq.matrix[2,i-(j-1)]) &&
					( forward[i-(j-2),j,1] <= forward[i-(j-2),j,2]))){
						SC  <- max(forward[i,j,1],forward[i,j,2])+match.score
					} else {
						SC  <-  max(forward[i,j,1],forward[i,j,2])- mismatch.score
					}
					if (x!=j){
						SC <- SC - penalty - abs(x-j)
					}
					if (SC > SCMax){
						SCMax  <-  SC
					}
				}
				forward[i+1,j,z]  <- SCMax
			}
		}
	}
	return(forward)
}





calculate.forward.modified <- function(seq.matrix,max.shift,penalty){

	seq.length <- length(seq.matrix[1,])
	#set scores as in indelligent
	match.score <- 1
	mismatch.score <- 0

	#init forward matrix
	#forward saves score for forward[i,j,k]. where i is the position in the sequence, 
	#j is the shift and k is the interpretation of ambiguity

	forward <- array(0,dim=c(seq.length+2,max.shift+1,2))

	#set starting conditions for forward calculation
	for ( i in 2:(max.shift+1)){
		forward[i,i:(max.shift+1),] <- (i-1)
	}

	#calculate score in forward direction

	#iterate over seqence length
	for ( i in 1:seq.length){
		jMax <- max.shift

		#check if iteration is smaller than possible shift
		if ((i-1) < max.shift){	
			jMax <- i-1
		}

		#iterate over possible shifts
		for (j in 1:(jMax+1)){
			for (z in 1:2){
				SCMax <- -111111

				#calculate score for eveery shift

				score=vector("numeric",seq.length)

				if(j==1 ){#&& seq.matrix[1,i] == seq.matrix[2,i]){
					index <- seq.matrix[1,] == seq.matrix[2,]
					score[index] <- apply(forward[index,,],1,max) + match.score
					score[!index] <- apply(forward[!index,,],1,max) - mismatch.score #max(forward[i,1,1],forward[i,1,2]) - mismatch.score
				}else{

				}




				for (x in 1:(jMax+1)){
					if (j==1 && x==1 && seq.matrix[1,i] != seq.matrix[2,i]){
						
						SC <- max(forward[i,x,1],forward[i,x,2]) - mismatch.score
					}else if (j==1 && seq.matrix[1,i]==seq.matrix[2,i]){
						SC <- max(forward[i,x,1],forward[i,x,2]) + match.score
					}else if (j != x){
						SC <- max(forward[i,x,1],forward[i,x,2]) + match.score
					}else if (((seq.matrix[3-z,i]==seq.matrix[1,i-(j-1)]) && 
					(forward[i-(j-2),j,1] >= forward[i-(j-2),j,2])) || ((seq.matrix[3-z,i]==seq.matrix[2,i-(j-1)]) &&
					( forward[i-(j-2),j,1] <= forward[i-(j-2),j,2]))){
						SC  <- max(forward[i,j,1],forward[i,j,2])+match.score
					} else {
						SC  <-  max(forward[i,j,1],forward[i,j,2])- mismatch.score
					}
					if (x!=j){
						SC <- SC - penalty - abs(x-j)
					}
					if (SC > SCMax){
						SCMax  <-  SC
					}
				}
				forward[i+1,j,z]  <- SCMax
			}
		}
	}
}

#calculates the backward score of the sequence
calculate.backward <- function(seq.matrix,max.shift,penalty){

	seq.length <- length(seq.matrix[1,])
	#set scores as in indelligent
	match.score <- 1
	mismatch.score <- 0


	#init backward matrix
	backward <- array(0,dim=c(seq.length+2,max.shift+1,2))


	#set starting conditions for backward matrix
	for ( i in 2:(max.shift+1)){
		backward[(seq.length-i+3),i:(max.shift+1),] <- (i-1)
	}

	#calculate score in backward direction

	#iterate over seqence length
	for (i in seq.length:1){
		jMax <- max.shift

		#check if iteration is smaller than possible shift
		if((seq.length-i) < max.shift){
			jMax <- seq.length - i
		}

		#iterate over possible shifts
		for ( j in 1:(jMax+1)){
			for (z in 1:2){
				SCMax <- -111111

				#calculate score for every shift
				for (x in 1:(jMax+1)){
					if(j==1 && x == 1 && seq.matrix[1,i] != seq.matrix[2,i]){
						SC <- max(backward[i+2,x,1],backward[i+2,x,2]) -mismatch.score
					}else if(j==1 && seq.matrix[1,i]==seq.matrix[2,i]){
						SC <- max(backward[i+2,x,1],backward[i+2,x,2])+match.score
					}else if(j!=x){
						SC=max(backward[i+2,x,1],backward[i+2,x,2]) + match.score
					}else if(((seq.matrix[z,i]==seq.matrix[2,i+(j-1)]) && (backward[i+j,j,1] >= backward[i+j,j,2])) || 
						((seq.matrix[z,i]==seq.matrix[1,i+(j-1)]) && (backward[i+j,j,1] <=backward[i+j,j,2]))){
						SC <- max(backward[i+2,j,1],backward[i+2,j,2]) +match.score
					}else{
						SC <- max(backward[i+2,j,1],backward[i+2,j,2]) - mismatch.score
					}
					if (x!=j){
						SC <- SC-penalty-abs(x-j)
					}
					if(SC > SCMax){
						SCMax <- SC
					}
				}
				backward[i+1,j,z]=SCMax
			}
		}
	}
	return(backward)
}


#calculates the whole sequence score
calculate.scores <- function(seq.matrix,max.shift,penalty){
	forward <- calculate.forward(seq.matrix,max.shift,penalty)
	backward <- calculate.backward(seq.matrix,max.shift,penalty)
	return(forward+backward)
}

#calculates all possible indels
calculate.poss.indels <- function(combined.matrix,seq.length,max.shift,fixed.shifts){
	#init indel object
	all.indels <- c()

	#iterate over combined.matrix[i,,]
	for (i in 2:(seq.length+1)){
		score <- -111111
		indel <- 0

		#iterate over shifts
		for (j in 0:max.shift){
			#check for fixed shifts
			if(is.null(fixed.shifts) || any(fixed.shifts==j)){
				#check if indel is possible
				if (combined.matrix[i,j+1,1]> score){
					score <- combined.matrix[i,j+1,1]
					indel <- j
				}
				if (combined.matrix[i,j+1,2]>score){
					score <- combined.matrix[i,j+1,2]
					indel <- j
				}
			}
		}
		if(!any(all.indels==indel)){
			all.indels <- c(all.indels,indel)
		}
	}
	return(all.indels)
}


calculate.phase.shift <- function(scores,sequences,max.shift,fixed.shifts,possible.indels){
	#apply scores from scores and resolve positions in sequences. Results are stored in result.
	#result[1,i] the value indicating which base from sequences goes into the first string. 
	#If result[1,i]==0 then the position is ambiguous.
	#result[2,i] the value indicating the phase shift which resolves the position.

	removed.short.indel <- 0
	result <- matrix(0,nrow=4,ncol=(seq.length+1))
	repeat{
		for(i in 1:seq.length){
			score <- -111111
			pos.one <- 0	#possibility to solve in position 1
			pos.two <- 0	#possibility to solve in position 2
			shift <- 0	#shift size
			
			for(j in 0:max.shift){
				if((is.null(fixed.shifts) && removed.short.indel==0) || 
					(any(fixed.shifts==j) && removed.short.indel==0) || 
					(any(x==j) && removed.short.indel==1)){
					if(scores[i+1,j+1,1] > score){
						score <- scores[i+1,j+1,1]
						pos.one <- 0
						pos.two <- 0
						shift <- j
					}
					if(scores[i+1,j+1,2] > score){
						score <- scores[i+1,j+1,2]
						pos.one <- 0
						pos.two <- 0
						shift <- j
					}
					if(((scores[i+1,j+1,1]==score) || (scores[i+1,j+1,2]==score)) && 
						(0+j)==result[2,i]){
						shift <- j
					}
					if((scores[i+1,j+1,1]==score) && (any(x==j))){
						pos.one <- 1
					}
					if((scores[i+1,j+1,2]==score) && (any(x==j))){
						pos.two <- 1
					}
				}
			}
		
			if(sequences[1,i]==sequences[2,i]){
				result[1,i+1] <- 1
			}else if (pos.one==1 && pos.two==0){
				result[1,i+1] <- 1
			}else if (pos.two==1 && pos.one ==0){
				result[1,i+1] <- 2
			}else{
				result[1,i+1] <- 0
			}
			result[2,i+1] <- shift
		}
		#remove phase shifts recovered at the number of consecutive positions smaller than the phase shift magnitude

		
		u <- c()
		pos.one <- 0

		for(i in 1:seq.length){
			if(result[2,i+1] != result[2,i]){
				for(z in 1:result[2,i+1]){
					if(i+z > seq.length){
						break
					}
					if(result[2,i+z+1]!=result[2,i+1]){
						pos.one <- 1
						break
					}
					if (z==result[2,i+1]){
						z <- z+1
					}
				}
			}
			if(result[2,i+1] != result[2,2]){
				shift <- 1
			}
			if(z>result[2,i+1] && !any(u==result[2,i+1])){
				u <- c(u,result[2,i+1])
			}
		}
		if(!all(u == x) && removed.short.indel==0){
			   pos.one <- 1
		}
		x <- u
		if(pos.one==1 || removed.short.indel==1){
			removed.short.indel <- removed.short.indel+1
		}
		if(removed.short.indel!=1){		#recalculate result if short indel was removed
			break
		}
	}

	return(result)
}

mark.longindel <- function(is.longindel,sequences,scores,shifts){
	#mark and rotate positions for long indels 
	if (is.longindel){
		SCa = shifts[2,2]
		for(i in 1:seq.length){
			if(shifts[2,i+1]>SCa){
				for(j in 0:(SCa-1)){
					shifts[2,i+j+1] <- SCa
					shifts[1,i+j+1] <- 0
				}
				i <- i+SCa
				SCa <- shifts[2,i+1]
			}else if(shifts[2,i+1]<SCa){
				SCa <- shifts[2,i+1]
				for(j in 1:SCa){
					shifts[1,i-j+1] <- 0
					shifts[2,i-j+1] <- SCa
				}
			}
		}
		SCa <- shifts[2,2]
		SCd <- 0
		for(i in 1:seq.length){
			if(shifts[2,i+1] != SCa){
				for(j in max(1,i-10):min(seq.length,i+10)){
					if(scores[j,SCa,2]==scores[j,shifts[2,i+1],2] && scores[j,SCa,2] >= scores[j,SCa,1] && 
						scores[j,shifts[2,i+1],2]>= scores[j,shifts[2,i+1],1]){
						shifts[1,j+1] <- 0
					}else if(scores[j,SCa,1] == scores[j,shifts[2,i+1],1] && 
						scores[j,SCa,2] <= scores[j,SCa,1] && 
						scores[j,shifts[2,i+1],2] <=scores[j,shifts[2,i+1],1]){
						shifts[1,j] <- 0
					}
					if(sequences[1,j] == sequences[2,j]){
						shifts[1,j+1] <- 1
					}
				}
				if(SCa != 0){
					SCd <- abs(Scd-1)
				}
				SCa <- shifts[2,i+1]
			}
			if(SCd==1){
				shifts[2,i+1] <- -shifts[2,i+1]
				if(shifts[1,i+1] >0 && sequences[1,i] != sequences[2,i]){
					shifts[1,i+1]  <-  3-shifts[1,i+1]
				}
			}
		}
	}
}

resolve.ambiguities <- function(sequences,scores,shifts){
	seq.length <- length(sequences[1,])
	if(any(shifts[2,]!=shifts[2,2])){
		repeat{
			alignLeft(1)
			alignRight(1)

			#calculate and mark the ambiguities that could potentially be resolved

			for(i in 1:seq.length){
				if(shifts[1, i+1]==0){
					j <- abs(shifts[2,i+1])

					if(i<=j){
						shifts[3,i+1] <- 1	#1 - could be resolved, 0 - cannot be resolved
					}else if(j==0){
						shifts[3,i+1] <- 1
					}else if(shifts[3,i-j+1] && abs(shifts[2,i-j+1])==j && shifts[1,i-j+1]==0){
						shifts[3,i+1] <- 0
					}else if(shifts[3,i-j+1]==1 && abs(shifts[2,i-j+1])==j && shifts[1,i-j+1]==0){
						shifts[3,i+1] <- 1
					}else{
						SCc <- 1
						while(shifts[1,j*SCc+i+1]==0 && abs(shifts[2,j*SCc+1+i])==j && 
							j*SCc+i+1<=seq.length){
							SCc <- SCc+1
						}
						if(abs(shifts[2,i-j+1]) < j && shifts[2,i+1] >0){
							shifts[3,i+1] <- 1
						}else if(j*SCc +1 > seq.length){ # evtl +1
							shifts[3, i+1] <- 1
						}else if(abs(shifts[2,j*SCc+i+1]) != j){
							shifts[3,i+1] <- 1
						}else if(abs(shifts[2,j*SCc+i+1+shifts[3,i+1]]) != j){
							shifts[3,i+1] <- 1
						}else if(shifts[2,i+1] > 0 && (SCc / 2 )!=floor(SCc/2) && 
							((sequences[shifts[1,i-j+1],i-j] == sequences[2,i] && 
								sequences[1,i] == sequences[3-shifts[1,j*SCc*i+1],j*SCc+1]) || 
							(sequences[shifts[1,i-j+1],i-j]==sequences[1,i] && 
								sequences[2,i]==sequences[3-shifts[1,j*SCc+i+1],j*SCc+i]))){
							shifts[3,i+1] <- 1
						}else if(shifts[2,i+1] >0 && (SCc/2) == floor(SCc/2) && 
							((sequences[shifts[1,i-j+1],i-j]==sequences[1,i] && 
								sequences[2,i] == sequences[3-shifts[1,j*SCc+i+1],j*SCc+i])||
							(sequences[shifts[1,i-j+1],i-j] == sequences[1,i] && 
								sequences[1,i] == sequences[3- shifts[1,j*SCc+i+1],j*SCc+i]))){
							shifts[3,i+1] <- 1
						}else if(abs(shifts[2,i-j+1]) != j && shifts[2,i+1] < 0){
							shifts[3,i+1] <- 1
						}else if(shifts[1,i-j+1] ==0){
							shifts[3,i+1] <- 0
						}else if(shifts[2,i+1] < 0 && (SCc/2) != floor(SCc/2) && 
							((sequences[3-shifts[1,i-j+1],i-j] == sequences[1,i] && 
								sequences[2,i]== sequences[shifts[1,j*SCc+i+1],j*SCc+i]) || 
							(sequences[3-shifts[1,i-j+1],i-j] == sequences[2,i] && 
								sequences[1,i]==sequences[shifts[1,j*SCc+i+1],j*SCc+i]))){
							shifts[3,i+1] <- 1
						}else if(shifts[2,i+1] < 0 && (SCc/2) == floor(SCc/2) && 
							(( sequences[3-shifts[1,i-j+1],i-j]==sequences[1,i-j] && 
								sequences[1,i] == sequences[shifts[1,j*SCc+i+1],j*SCc+i]) || 
							(sequences[3-shifts[1,i-j+1],i-j] == sequences[2,i] &&
								sequences[2,i]==sequences[shifts[1,j*SCc+i+1],j*SCc+i]))){
							shifts[3,i+1] <- 1
						}else{
							shifts[3,i+1] <- 0
						}
					}
				}
			}
			j <- 0

			for(i in 1:seq.length){
				if(shifts[1,i+1]==0 && shifts[3,i+1]==1){
					for(z in 1:(max.shift+shifts[4,1])){
						if(z>=i || z> seq.length -i){
							break
						}
						if(shifts[2,i-z+1] != shifts[2,i+z+1]){
							break
						}
					}
					if(z<=i && z-shifts[3,i+1] <= max.shift+1 && z <= seq.length -(i+1)){
						ind1 <- abs(shifts[2,i-z+1])
						ind2 <- abs(shifts[2,i+z+1])
						if(!is.longindel){
							if((i>= ind1) && (i>= ind2) && (seq.length -i)>ind1 && 
								(seq.length-i)>ind2){
								if(((sequences[2,i]==sequences[1,i-ind1] && 
									scores[i-ind+1,ind1,1]>=scores[i-ind1+1,ind1,2]) || 
								(sequences[2,i] == sequences[2,i-ind1] && 
									scores[i-ind1+1,ind1,1] <= scores[i-ind1+1,ind1,2])) && 
								((sequences[1,i] == sequences[2,i+ind1] && 
									scores[i+ ind1+1,ind1,1] >= scores[i+ind1+1,ind1,2])||
								(sequences[1,i] == sequences[1,i+ind1] && 
									scores[i+ind1+1,ind1,1] <= scores[i+ind1+1,ind1,2]))){
									MX5[1,1] <- 1
								}else{
									MX5[1,1] <- 0
								}
								if(((sequences[1,i]==sequences[1,i-ind1] && 
									scores[i-ind1+1,ind1,1]>=scores[i-ind1+1,ind1,2])||
								(sequences[1,i] == sequences[2,i-ind] &&
									scores[i-ind1+1,ind1,1] <= scores[i-ind1+1,ind1,2]))&&
								((sequences[2,i]==sequences[2,i+ind1] && 
									scores[i+ind1+1,ind1,1]>=scores[i+ind1+1,ind1,2])||
								(sequences[2,i]==sequences[1,i+ind1] && 
									scores[i+ind1+1,ind1,1]<=scores[i+ind1+1,ind1,2]))){
									MX5[1,2] <- 1
								}else{
									MX5[1,2] <- 0
								}
								if(((sequences[2,i]==sequences[1,i-ind2]&&
									scores[i-ind2+1,ind2,1]>=scores[i-ind2+1,ind2,2])||
								(sequences[2,i]==sequences[2,i-ind2]&&
									scores[i-ind2+1,ind2,1]<= scores[i-ind2+1,ind2,2]))&&
								((sequences[1,i]==sequences[2,i+ind2]&&
									scores[i+ind2+1,ind2,1]>=scores[i+ind2+1,ind2,2])||
								(sequences[1,i]==sequences[1,i+ind2]&&
									scores[i+ind2+1,ind2,1]<=scores[i+ind2+1,ind2,2]))){
									MX5[2,1] <- 1
								}else{
									MX5[2,1] <- 0
								}
								if(((sequences[1,i]==sequences[1,i-ind2] && 
									scores[i-ind2+1,ind2,1]>=scores[i-ind2+1,ind2,2])||
								(sequences[1,i]==sequences[2,i-ind2]&&
									scores[i-ind2+1,ind2,1]<=scores[i-ind2+1,ind2,2]))&&
								((sequences[2,i]==sequences[2,i+ind2]&&
									scores[i+ind2+1,ind2,1]>=scores[i+ind2+1,ind2,2])||
								(sequences[2,i]==sequences[1,i+ind2]&&
									scores[i+ind2+1,ind2,1]<=scores[i+ind2+1,ind2,2]))){
									MX5[2,2] <- 1
								}else{
									MX5[2,2] <- 0
								}
								if((ind1 < ins2) &&
									((sequences[2,i]==sequences[1,i-ind2]&&
										scores[i-ind1+1,ind1,1]>=scores[i-ind1+1,ind1,2])||
									(sequences[2,i]==sequences[2,i-ind1]&&
										scores[i-ind1+1,ind1,1]<=scores[i-ind1+1,ind1,2]))&&
									((sequences[1,i]==sequences[2,i+ind2]&&
										scores[i+ind2+1,ind2,1]>=scores[i+ind2+1,ind2,2])||
									(sequences[1,i]==sequences[1,i+ind2]&&
										scores[i+ind2+1,ind2,1]<=scores[i+ind2+1,ind2,2]))){
									MX5[4,1] <- 1
								}else{
									MX5[4,1] <- 0
								}
								if((ind1<ind2)&&
									((sequences[1,i]==sequences[1,i-ind]&&
										scores[i-ind1+1,ind1,1]>=scores[i-ind1+1,ind1,2])||
									(sequences[1,i]==sequences[2,i-ind1]&&
										scores[i-ind1+1,ind1,1]<=scores[i-ind1+1,ind1,2]))&&
									((sequences[2,i]==sequences[2,i+ind2]&&
										scores[i+ind2+1,ind2,1]>=scores[i+ind2+1,ind2,2])||
									(sequences[2,i]==sequences[1,i+ind2]&&
										scores[i+ind2+1,ind2,1]<=scores[i+ind2+1,ind2,2]))){
									MX5[4,2] <- 1
								}else{
									MX5[4,1] <- 0
								}
							}else if(i<=ind1){
								if(((sequences[1,i]==sequences[2,i+ind1]&&
									scores[i+ind1+1,ind1,1]>=scores[i+ind1+1,ind1,2])||
								(sequences[1,i]==sequences[1,i+ind1]&&
									scores[i+ind1+1,ind1,1]<=scores[i+ind1+1,ind1,2]))){
									MX5[1,1] <- 1
								}else{
									MX5[1,1] <- 0
								}
								if(((sequences[2,i]==sequences[2,i+ind1]&&
									scores[i+ind1+1,ind1,1]>=scores[i+ind1+1,ind1,2])||
								(sequences[2,i]==sequences[1,i+ind1]&&
									scores[i+ind1+1,ind1,1]<=scores[i+ind1+1,ind1,2]))){
									MX5[1,2] <- 1
								}else{
									MX5[1,2] <- 0
								}
								MX5[2,] <- 0
							}else if(seq.length-i< ind2){
								MX5[1,] <- 0
								if(((sequences[2,i]==sequences[1,i-ind2]&&
									scores[i-ind2+1,ind2,1]>=scores[i-ind2+1,ind2,2])||
								(sequences[2,i]==sequences[2,i-ind2]&&
									scores[i-ind2+1,ind2,1]<=scores[i.ind2+1,ind2,2]))){
									MX5[2,1] <- 1
								}else{
									MX5[2,1] <- 0
								}
								if(((sequences[1,i]==sequences[1,i-ind2]&&
									scores[i-ind2+1,ind2,1]>=scores[i-ind2+1,ind2,2])||
								(sequences[1,i]==sequences[2,i-ind2]&&
									scores[i-ind2+1,ind2,1]<=scores[i-ind2+1,ind2,2]))){
									MX5[2,2] <- 1
								}else{
									MX5[2,2] <- 0
								}
							}
							if(shifts[2,i-z+1]<shifts[2,i+z+1]||i<=ind1){
								if(((sequences[1,i]==sequences[2,i+ind2]&&
									scores[i+ind2+1,ind2,1]>=scores[i+ind2+1,ind2,2])||
								(sequences[1,i]==sequences[1,i+ind2]&&
									scores[i+ind2+1,ind2,1]<=scores[i+ind2+1,ind2,2]))){
									MX5[3,1] <- 1
								}else{
									MX5[3,1] <- 0
								}
								if(((sequences[2,i]==sequences[1,i+ind2]&&
									scores[i+ind2+1,ind2,1]>=scores[i+ind2+1,ind2,2])||
								(sequences[2,i]==sequences[1,i+ind2]&&
									scores[i+ind2+1,ind2,1]<=scores[i+ind2+1,ind2,2]))){
									MX5[3,2] <- 1
								}else{
									MX5[3,2] <- 0
								}
							}else{
								if(((sequences[2,i]==sequences[1,i-ind1]&&
									scores[i-ind1+1,ind1,1]>=scores[i-ind1+1,ind1,2])||
								(sequences[2,i]==sequences[2,i-ind1]&&
									scores[i-ind1+1,ind1,1]<=scores[i-ind1+1,ind1,2]))){
									MX5[3,1] <- 1
								}else{
									MX5[3,1] <- 0
								}
								if(((sequences[1,i]==sequences[1,i-ind1]&&
									scores[i-ind1+1,ind1,1]>=scores[i-ind1+1,ind1,2])||
								(sequences[2,i]==sequences[2,i-ind1]&&
									scores[i-ind1+1,ind1,1]<=scores[i-ind1+1,ind1,2]))){
									MX5[3,2] <- 1
								}else{
									MX5[3,2] <- 0
								}
							}
							if(ind1==0||ind2==0){
								MX5[4,] <- 0
							}	
						}else{  #long indel
							if(shifts[2,i-z+1]>0){
								if(((sequences[2,i]==sequences[1,i-ind1]&&
									shifts[1,i-ind1+1] !=2)||
								(sequences[2,i]==sequences[2,i-ind2]&&shifts[1,i-ind1+1]!=1))){
									MX5[3,1] <- 1
								}else{
									MX5[3,1] <- 0
								}
								if(((sequences[1,i]==sequences[1,i-ind1]&&shifts[1,i-ind1+1]!=2)||
									(sequences[1,i]==sequences[2,i-ind1]&&shifts[1,i-ind1]!=1))){
									MX5[3,2] <- 1
								}else{
									MX5[3,2] <- 0
								}
							}else{
								if(((sequences[1,i]==sequences[2,i-ind1] && shifts[1,i-ind1+1]!=2)||
									(sequences[1,i]==sequences[1,i-ind1]&&shifts[1,i-ind1+1]!=1))){
									MX5[4,1] <- 1
								}else{
									MX5[4,1] <- 0
								}
								if(((sequences[2,i]==sequences[2,i-ind1]&&shifts[1,i-ind1+1]!=2)||
									(sequences[2,i]==sequences[1,i-ind1]&&shifts[1,i-ind1+1]!=1))){
									MX5[4,2] <- 1
								}else{
									MX5[4,2] <- 0
								}
							}
							if(shifts[2,i+z+1]>0){
								if(((sequences[1,i]==sequences[2,i+ind2]&&shifts[1,i+ind2+1]!=2)||
									(sequences[1,i]==sequences[1,i+ind2]&&shifts[1,i+ind2+2]!=1))){
									MX5[4,1] <- 1
								}else{
									MX5[4,1] <- 0
								}
								if(((sequences[2,i]==sequences[2,i+ind2]&&shifts[1,i+ind2+1]!=2)||
									(sequences[2,i]==sequences[1,i+ind2]&&shifts[1,i+ind2+1]!=1))){
									MX5[4,2] <- 1
								}else{
									MX5[4,2] <- 0
								}
							}else{
								if(((sequences[2,i]==sequences[1,i+ind2]&&shifts[1,i+ind2+1]!=2)||
									(sequences[2,i]==sequences[2,i+ind2]&&shifts[1,i+ind2+1]!=1))){
									MX5[4,1] <- 1
								}else{
									MX5[4,1] <- 0
								}
								if(((sequences[1,i]==sequences[1,i+ind2]&&shifts[1,i+ind2+1]!=2)||
									(sequences[1,i]==sequences[2,i+ind2]&&shifts[1,i+ind2+1]!=1))){
									MX5[4,2] <- 1
								}else{
									MX5[4,2] <- 0
								}
							}
							if(i >= ind1 && i>= ind2 && (seq.length-i)>ind1&&
								(seq.length-i)<ind2){
								if(shifts[2,i-z+1]>0){
									if(((sequences[2,i]==sequences[1,i-ind1]&&shifts[1,i-ind1+1]!=2)||
										(sequences[2,i]==sequences[2,i-ind1]&&shifts[1,i-ind1+1]!=1))&&
									((sequences[1,i]==sequences[2,i+ind1]&&shifts[1,i+ind1+1]!=2)||
										(sequences[1,i]==sequences[1,i+ind1]&&shifts[1,i+ind1+1]!=1))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((sequences[1,i]==sequences[1,i-ind1]&&shifts[1,i-ind1+1]!=2)||
										(sequences[1,i]==sequences[2,i-ind1]&&shifts[1,i-ind1]!=1))&&
									((sequences[2,i]==sequences[2,i+ind1]&&shifts[1,i+ind1+1]!=2)||
										(sequences[2,i]==sequences[1,i+ind1]&&shifts[1,i+ind1+1]!=1))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
								}else{
									if(((sequences[1,i]==sequences[2,i-ind1]&&shifts[1,i-ind1]!=2)||
										(sequences[1,i]==sequences[1,i-ind1]&&shifts[1,i-ind1+1]!=1))&&
									((sequences[2,i]==sequences[1,i+ind1]&&shifts[1,i+ind1+1]!=2)||
										(sequences[2,i]==sequences[1,i+ind1]&&shifts[1,i+ind1+1]!=1))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((sequences[2,i]==sequences[2,i-ind1]&&shifts[1,i-ind1+1]!=2)||
										(sequences[2,i]==sequences[1,i-ind1]&&shifts[1,i-ind1+1]!=1))&&
									((sequences[1,i]==sequences[1,i+ind1]&&shifts[1,i+ind1]!=2)||
										(sequences[1,i]==sequences[2,i+ind1]&&shifts[1,i+ind1+1]!=1))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
								}
								if(shifts[2,i+z+1]>0){
									if(((sequences[2,i]==sequences[1,i-ind2]&&shifts[1,i-ind2+1]!=2)||
										(sequences[2,i]==sequences[2,i-ind2]&&shifts[1,i-ind2+1]!=1))&&
									((sequences[1,i]==sequences[2,i+ind2]&&shifts[1,i+ind2+1]!=2)||
										(sequences[1,i]==sequences[1,i+ind2]&&shifts[1,i+ind2+1]!=1))){
										MX5[2,1] <- 1
									}else{
										MX5[2,1] <- 0
									}
									if(((sequences[1,i]==sequences[1,i-ind2]&&shifts[1,i-ind2+1]!=2)||
										(sequences[1,i]==sequences[2,i-ind2]&&shifts[1,i-ind2+1]!=1))&&
									((sequences[2,i]==sequences[2,i+ind2]&&shifts[1,i+ind2+1]!=2)||
										(sequences[2,i]==sequences[1,i+ind2]&&shifts[1,i+ind2+1]!=1))){
										MX5[2,2] <- 1
									}else{
										MX5[2,2] <- 0
									}
								}else{
									if(((sequences[1,i]==sequences[2,i-ind2]&&shifts[1,i-ind2+1]!=2)||
										(sequences[1,i]==sequences[1,i-ind2]&&shifts[1,i-ind2+1]!=1))&&
									((sequences[2,i]==sequences[1,i+ind2]&&shifts[1,i+ind2]!=2)||
										(sequences[2,i]==sequences[2,i+ind2]&&shifts[1,i+ind2+1]!=1))){
										MX5[2,1] <- 1
									}else{
										MX5[2,1] <- 0
									}
									if(((sequences[2,i]==sequences[2,i-ind2]&&shifts[1,i-ind2+1]!=2)||
										(sequences[1,i]==sequences[1,i-ind2]&&shifts[1,i-ind2]!=1))&&
									((sequences[1,i]==sequences[1,i+ind2]&&shifts[1,i+ind2+1]!=2)||
										(sequences[1,i]==sequences[2,i+ind2]&&shifts[1,i+ind2+1]!=1))){
										MX5[2,2] <- 1
									}else{
										MX5[2,2] <- 0
									}
								}
								######################################################
							}else if(i<=ind1){
									if(((sequences[1,i]==sequences[2,i+ind1]&&shifts[1,i+ind1+1]!=2)||(sequences[1,i]==sequences[1,i+ind1]&&shifts[1,i+ind1+1]!=1))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((sequences[2,i]==sequences[2,i+ind1]&&shifts[1,i+ind1]!=2)||(sequences[2,i]==sequences[1,i+ind1]&&shifts[1,i+ind1+1]))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
									MX5[2,] <- 0
							}else if((seq.length-i)<ind2){
									MX5[1,] <- 0
									if(((sequences[2,i]==sequences[1,i-ind2]&&shifts[1,i-ind1]!=2)||(sequences[2,i]==sequences[2,i-ind2]&&shifts[1,i-ind1]!=1))){
										MX5[2,1] <- 1
									}else{
										MX5[2,1] <- 0
									}
									if(((sequences[1,i]==sequences[1,i-ind2]&&shifts[1,i-ind1]!=2)||(sequences[1,i]==sequences[2,i-ind2]&&shifts[1,i-ind1]!=1))){
										MX5[2,2] <- 1
									}else{
										MX5[2,2] <- 0
									}
							}
						}

						if(ind1==0){
							MX5[1,] <- 0
						}
						if(ind2==0){
							MX5[2,] <- 0
						}
						if(!is.longindel){
							if(ind1>ind2){
								if(shifts[2,i]==ind1){
									if(z<=ind1){
										MX5[1,] <- 0
										if(z<=ind2){
											MX5[3,] <- 0
										}
									}

									if((shifts[4,i+1]-z)<0){
										MX5[2,] <- 0
									}
								}else{
									MX5[1,] <- 0
									MX5[3,] <- 0
								}
							}else{ #ind1<ind2
								if(shifts[2,i+1]==ind1){
									if(shifts[3,i+1]<=ind2){
										MX5[2,] <- 0
										if(ind1 != 0 && z > shifts[4,i+1]){
											MX5[3,] <- 0
										}
									}
									if(z>ind2){
										MX5[2,] <- 0
									}
								}else{
									if(z+shifts[4,i+1]<=ind2){
										MX5[2,] <- 0
									}
									if(z>ind1){
										MX5[4,] <- 0
									}
									MX5[1,] <- 0
								}
							}
						}else{
							if(ind1>ind2){
								if(abs(shifts[2,i+1])==ind1){
									if(z<=ind1){
										MX5[1,] <- 0
									}
									if((shifts[4,i+1]-z)<0){
										MX5[c(2,4),] <- 0
									}
								}else{
									MX5[c(1,3),] <- 0
								}
							}else{
								if(shifts[2,i+1] <= ind2){
									if(shifts[4,i+1]<=ind2){
										MX5[2,] <- 0
									}
									if(ind1!=0 &&z > shifts[4,i+1]){
										MX5[4,] <- 0
									}
									if(z>ind2){
										MX5[2,] <- 0
									}
								}else{
									if((z+shifts[4,i+1])<=ind2){
										MX5[2,] <- 0
									}
									MX5[c(1,3),] <- 0
								}
							}
						}
						if(!is.longindel){
							if(all(MX5[1,]==1)){
							}else if(all(MX5[2,]==1)){
							}else if(all(MX5[c(1,2),1]==1)&&all(MX5[c(1,2),2]==0)){
								shifts[1,i+1] <- 1
								scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1
								scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
							}else if(all(MX5[c(1,2),1]==0)&&all(MX5[c(1,2),2]==1)){
								shifts[1,i+1] <- 2
								scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1
								scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
							}else if((MX5[1,1]==1&&MX5[2,2]==1)||(MX5[1,2]==1&&MX5[2,1]==1)){
							}else if(all(MX5[c(1,3),1]==1)&&all(MX5[c(1,3),2]==0)){
								shifts[1,i+1] <- 1
								scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1
								scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
							}else if(all(MX5[c(1,3),1]==0)&&all(MX5[c(1,3),2]==1)){
								shifts[1,i+1] <- 2
								scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1
								scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
							}else if(((MX5[1,1]==1&&MX5[3,2]==1)||(MX5[1,2]==1&&MX5[3,1]==1))
								&&z< ind1 &&shifts[4,i+1]>ind1){
							}else if(((MX5[2,1]==1&&MX5[3,2]==1)||(MX5[2,2]==1&&MX5[3,1]==1))
								&&z<=ind1 &&shifts[4,i+1]>ind2 &&ind2 >ind1){
							}else if(((MX5[2,1]==1&&MX5[3,2]==1)||(MX5[2,2]==1&&MX5[3,1]==1))
								&&shifts[2,i+1]==ind2&&shifts[1,i+ind2+1]==0){
							}else if(((MX5[2,1]==1&&MX5[3,2]==1)||(MX5[2,2]==1&&MX5[3,1]==1))
								&&shifts[2,i+1]==ind1&&shifts[4,i+1]>=z){
							}else if(((MX5[1,1]==1&&MX5[3,2]==1)||(MX5[1,2]==1&&MX5[3,1]==1))
								&&shifts[2,i+1]==ind1&&ind2>ind1&&shifts[4,i+1]>=(z+ind1)){
							}else if(xor(MX5[1,1]==1,MX5[1,2]==1)){
								if(MX5[1,1]==1){
									shifts[1,i+1] <- 1
									scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1 
									scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
								}else{
									shifts[1,i+1] <- 2
									scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1 
									scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
								}
							}else if(xor(MX5[2,1]==1,MX5[2,2]==1)){
								if(MX5[2,1]==1){
									shifts[1,i+1] <- 1
									scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1 
									scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
								}else{
									shifts[1,i+1] <- 2
									scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1 
									scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
								}
							}else if(all(MX5[3,]==1)&&all(MX5[4,]==0)){
							}else if(all(MX5[3,]==1)&&all(MX5[4,]==c(1,0))){
								shifts[1,i+1] <- 1
								scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1 
								scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
							}else if(all(MX5[3,]==1)&&all(MX5[4,]==c(0,1))){
								shifts[1,i+1] <- 2
								scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1 
								scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
							}else if(MX5[3,1]==1){
								shifts[1,i+1] <- 1
								scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1 
								scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
							}else if(MX5[3,2]==1){
								shifts[1,i+1] <- 2
								scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1 
								scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
							}
						}else{

							if(all(MX5[1,]==1)){
							}else if(all(MX5[2,]==1)){
							}else if(((MX5[1,1]==1&&MX5[3,2]==1)||(MX5[1,2]==1&&MX5[3,1]==1))){
							}else if(((MX5[2,1]==1&&MX5[4,2]==1)||(MX5[2,2]==1&&MX5[4,1]==1))){
							}else if(((MX5[1,1]==1&&MX5[4,2]==1)||(MX5[1,2]==1&&MX5[4,1]==1))){
							}else if(((MX5[2,1]==1&&MX5[3,2]==1)||(MX5[2,2]==1&&MX5[3,1]==1))){
							}else if(((MX5[4,1]==1&&MX5[3,2]==1)||(MX5[4,2]==1&&MX5[3,1]==1))){
							}else if(xor(MX5[1,1]==1,MX5[1,2]==1)){
								if(MX5[1,1]==1){
									shifts[1,i+1] <- 1
									scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1 
									scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
								}else{
									shifts[1,i+1] <- 2
									scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1 
									scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
								}
							}else if(xor(MX5[2,1]==1,MX5[2,2]==1)){
								if(MX5[2,1]==1){
									shifts[1,i+1] <- 1
									scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1 
									scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
								}else{
									shifts[1,i+1] <- 2
									scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1 
									scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
								}
							}else if(all(MX5[3,]==1)){
							}else if(all(MX5[4,]==1)){
							}else if(all(MX5[4,]==c(1,0))){
								X4[1,i+1] <- 1
								scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1 
								scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
							}else if(all(MX5[4,]==c(0,1))){
								shifts[1,i+1] <- 2
								scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1 
								scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
							}else if(all(MX5[3,]==c(1,0))){
								X4[1,i+1] <- 1
								scores[i+1,ind2,1] <- max(scores[i+1,ind2,])+1 
								scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
							}else if(all(MX5[3,]==c(0,1))){
								shifts[1,i+1] <- 2
								scores[i+1,ind2,2] <- max(scores[i+1,ind2,])+1 
								scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
							}
						}
					}else{
						if(shifts[2,i]==shifts[2,i+2]||i==1||i==seq.length){
							if(i!=seq.length){
								ind1 <- abs(shifts[2,i+2])
								ind2 <- shifts[2,i+2]
							}else{
								ind1 <- abs(shifts[2,i])
								ind2 <- shifts[2,i]
							}
							if(i<=ind1){
								if(((sequences[1,i]==sequences[2,i+ind1]&&
									scores[i+ind1+1,ind1,1] >= scores[i+ind1+1,ind1,2])||
								(sequences[1,i]==sequences[1,i+ind1]&&
									scores[i+ind1+1,ind1,1] <= scores[i+ind1+1,ind1,2]))){
									MX5[1,1] <- 1
								}else{
									MX5[1,1] <- 0
								}
								if(((sequences[2,i]==sequences[2,i+ind1]&&
									scores[i+ind1+1,ind1,1] >= scores[i+ind1+1,ind1,2])||
								(sequences[2,i]==sequences[2,i-ind2]&&
									scores[i+ind1+1,ind1,1] <= scores[i+ind1+1,ind1,2]))){
									MX5[1,2] <- 1
								}else{
									MX5[1,2] <- 0
								}
							}else if(seq.length - i < ind1){
								if(ind2 > 0){
									if(((sequences[2,i] == sequences[1,i - ind1] &&
										scores[i - ind1 + 1,ind1,1] >= scores[i - ind1 + 1,ind1,2])||
									(sequences[2,i] == sequences[2,i - ind1] &&
										scores[i - ind1 + 1,ind1,1] <= scores[i - ind1 + 1,ind1,2]))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((sequences[1,i] == sequences[1,i - ind1] &&
										scores[i - ind1 + 1,ind1,1] >= scores[i - ind1 + 1,ind1,2]) ||
									(sequences[1,i] == sequences[2,i - ind1] &&
										scores[i - ind1 + 1,ind1,1] <= scores[i - ind1 + 1,ind1,2]))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
								}else{
									if(((sequences[1,i] == sequences[2,i-ind1] &&
										shifts[1,i-ind1+1]!=2)||
									(sequences[1,i] == sequences[1,i-ind1] &&
										shifts[1,i-ind1+1]!=1))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((sequences[2,i] == sequences[2,i-ind1] &&
										shifts[1,i-ind1+1]!=2)||
									(sequences[2,i] == sequences[1,i-ind1] &&
										shifts[1,i-ind1+1]!=1))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
								}
							}else{
								if(ind2 > 0){
									if(((sequences[2,i] == sequences[1,i-ind1] && 
										scores[i-ind1+1,ind1,1] >= scores[i-ind1+1,ind1,2])||
									(sequences[2,i] == sequences[2,i-ind1]&&
										scores[i-ind1+1,ind1,1]<= scores[i-ind1+1,ind1,2])) &&
									((sequences[1,i] == sequences[2,i+ind1] &&
										scores[i+ind1+1,ind1,1]>=scores[i+ind1+1,ind1,2]) ||
									(sequences[1,i]==sequences[1,i+ind1]&&
										scores[i+ind1+1,ind1,1] <= scores[i+ind1+1,ind1,2]))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((sequences[1,i] == sequences[1,i-ind1] && 
										scores[i-ind1+1,ind1,1] >= scores[i-ind1+1,ind1,2])||
									(sequences[1,i]==sequences[2,i-ind1] && 
										scores[i-ind1+1,ind1,1]<=scores[i-ind1+1,ind1,2])) &&
									((sequences[2,i]==sequences[2,i+ind1] &&
										scores[i+ind1+1,ind1,1] >= scores[i+ind1+1,ind1,2])||
									(sequences[2,i]==sequences[1,i+ind1] &&
										scores[i+ind1+1,ind1,1] <= scores[i+ind1+1,ind1,2]))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
								}else{
									if(((sequences[1,i] == sequences[2,i-ind1] &&shifts[1,i-ind1+1] !=2) || 
									(sequences[1,i] == sequences[1,i-ind1] && shifts[1,i-ind1+1] !=1)) &&
									((sequences[2,i] == sequences[1,i+ind1] && shifts[1,i+ind1+1] != 2) || 
									(sequences[2,i] == sequences[2,i+ind1] && shifts[1,i+ind1+1] !=1))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((sequences[2,i] == sequences[2,i-ind1] && shifts[1,i-ind1+1] !=2) ||
										(sequences[2,i] == sequences[1,i-ind1] && shifts[1,i-ind1+1]!=1)) &&
									((sequences[1,i] == sequences[1,i+ind1] && shifts[1,i+ind1+1]!=2)||
										(sequences[1,i] == sequences[2,i+ind1] && shifts[1,i+ind1+1] !=1))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
								}
							}
							if(MX5[1,1] ==1 &&MX5[1,2] ==0){
								shifts[1,i+1] <- 1
								scores[i+1,ind1,1] <- max(scores[i+1,ind1,])+1
							}else if(MX5[1,1] == 0 && MX5[1,2] ==1){
								shifts[1,i+1] <- 2
								scores[i+1,ind1,2] <- max(scores[i+1,ind1,])+1
							}
						}
					}
					if(shifts[1,i+1] != 0){
						j <- 1
					}
				}
			}
			if(j!=1){
				break
			}
		}

	}
}






#############################################################################################


cat("sequence: ",sequence,"\n")

first.calc <- seq.toMatrix(sequence,max.shift)
#save generated sequence matrix
sequences <- first.calc$matrix

cat("sequences matrix: ")
print(sequences)
#save if three fold degenerate bases exist in sequence
three.fold <- first.calc$is.three.fold

#save ambiguities (number of three fold degenerates)
ambig1 <- first.calc$ambig

#save real shift
shift <- first.calc$shift

#calculate forward matirx
second.calc <- calculate.forward(sequences,shift,penalty)

#calculate backward matrix
third.calc <- calculate.backward(sequences,shift,penalty)

scores <- calculate.scores(sequences,shift,penalty)

cat("forward: ")
print(second.calc)
cat("backward: ")
print(third.calc)
cat("combined: ")
print(scores)