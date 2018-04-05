suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))

#make command line options 
option_list <- list(
	make_option(c("-i","--input-file"), type="character",default=NULL,
		    help="file path to sequence data. If set -s is ignored."),
	make_option(c("-s","--string"), type="character",default=NULL,
		    help="sequence string"),
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
	make_option(c("-n","--not-align"),action="store_false",dest="align",
		help="Option for not displaying reconstructed allelic sequences."),
	make_option(c("-g","--long-indel"),action="store_true",default=FALSE,
		help="Display alternative, longer variants of reconstructed insertions.
		[default \"%default\"]")
	)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)





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


#check if input or string is set
if(!is.null(opt$input)){
	sequence <- trimws(read_file(opt$input))
}else if (!is.null(opt$string)){
	sequence <- opt$string
}else{
	print_help(opt_parser)
	stop("Please provide an input file or string.\n",call.=FALSE)
}	


#convert sequence to all uppercase
sequence <- toupper(sequence)
seq.asvector=strsplit(sequence,"")[[1]]
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

seq.length <- nchar(sequence)
#maximum shift is limited to half of sequence length
if(max.shift*2>seq.length){
	max.shift=floor(nchar(sequence)/2)
}

if (seq.length!=0){ #if sequence exists
	

	#set scores as in indelligent
	match.score=1
	mismatch.score=0
	

	#initialize matric MX1
	MX1 <- matrix(,nrow=2,ncol=seq.length)
	#fill MX1
	pos.of=seq.asvector=="R"
	MX1[1,pos.of]="A"
	MX1[2,pos.of]="G"
	pos.of=seq.asvector=="Y"
	MX1[1,pos.of]="C"
	MX1[2,pos.of]="T"
	pos.of=seq.asvector=="S"
	MX1[1,pos.of]="G"
	MX1[2,pos.of]="C"
	pos.of=seq.asvector=="W"
	MX1[1,pos.of]="A"
	MX1[2,pos.of]="T"
	pos.of=seq.asvector=="M"
	MX1[1,pos.of]="A"
	MX1[2,pos.of]="C"
	pos.of=seq.asvector=="K"
	MX1[1,pos.of]="G"
	MX1[2,pos.of]="T"
	#fill remaining characters
	pos.of=is.na(MX1[1,])
	MX1[1,pos.of]=seq.asvector[pos.of]
	MX1[2,pos.of]=seq.asvector[pos.of]
	#check for 3-fold degenerate bases
	if (!(MX1=="A" || MX1=="T" || MX1=="G" || MX1=="C")){
		bd <- TRUE
		ambig1 <- sum(!(MX1=="A" | MX1=="T" | MX1=="G" | MX1=="C"))
	}else{
		bd <- FALSE
	}
#	print(MX1)
	#init calculation matrices
	MX2 <- array(0,dim=c(seq.length+2,max.shift+1,2))
	MX3 <- array(0,dim=c(seq.length+2,max.shift+1,2))
	
	#init conditions for MX2 and MX3 (first 1 to max.shift+1 as the corresponding index 
	#in MX2, last max.shift to seq.length+1 as falling cooresponding index e.g. 10 9 8 ... 0)
	for ( i in 2:(max.shift+1)){
		MX2[i,i:(max.shift+1),] <- (i-1)
		MX3[(seq.length-i+3),i:(max.shift+1),] <- (i-1)
	}
		
#	print(MX2)
#	print(MX3)	
	#calculate MX2 in forward direction
	for ( i in 1:seq.length){
		jMax <- max.shift
		if ((i-1) < max.shift){	#hier evtl fehler?
			jMax <- i-1
		}
		for (j in 1:(jMax+1)){
			for (z in 1:2){
				SCMax <- -111111
				for (x in 1:(jMax+1)){
					if (j==1 && x==1 && MX1[1,i] != MX1[2,i]){
						
						SC <- max(MX2[i,x,1],MX2[i,x,2]) - mismatch.score
					}else if (j==1 && MX1[1,i]==MX1[2,i]){
						SC <- max(MX2[i,x,1],MX2[i,x,2]) + match.score
					}else if (j != x){
						SC <- max(MX2[i,x,1],MX2[i,x,2]) + match.score
					}else if (((MX1[3-z,i]==MX1[1,i-(j-1)]) && 
					(MX2[i-(j-2),j,1] >= MX2[i-(j-2),j,2])) || ((MX1[3-z,i]==MX1[2,i-(j-1)]) &&
					( MX2[i-(j-2),j,1] <= MX2[i-(j-2),j,2]))){
						SC  <- max(MX2[i,j,1],MX2[i,j,2])+match.score
					} else {
						SC  <-  max(MX2[i,j,1],MX2[i,j,2])- mismatch.score
					}
					if (x!=j){
						SC <- SC - penalty - abs(x-j)
					}
					if (SC > SCMax){
						SCMax  <-  SC
					}
				}
				MX2[i+1,j,z]  <- SCMax
			}
		}
	}

#	print(MX2)
	#calculate MX3 in reverse direction
	for (i in seq.length:1){
		jMax <- max.shift
		if((seq.length-i) < max.shift){
			jMax <- seq.length - i
		}
		for ( j in 1:(jMax+1)){
			for (z in 1:2){
				SCMax <- -111111
				for (x in 1:(jMax+1)){
					if(j==1 && x == 1 && MX1[1,i] != MX1[2,i]){
						SC <- max(MX3[i+2,x,1],MX3[i+2,x,2]) -mismatch.score
					}else if(j==1 && MX1[1,i]==MX1[2,i]){
						SC <- max(MX3[i+2,x,1],MX3[i+2,x,2])+match.score
					}else if(j!=x){
						SC=max(MX3[i+2,x,1],MX3[i+2,x,2]) + match.score
					}else if(((MX1[z,i]==MX1[2,i+(j-1)]) && (MX3[i+j,j,1] >= MX3[i+j,j,2])) || 
						((MX1[z,i]==MX1[1,i+(j-1)]) && (MX3[i+j,j,1] <=MX3[i+j,j,2]))){
						SC <- max(MX3[i+2,j,1],MX3[i+2,j,2]) +match.score
					}else{
						SC <- max(MX3[i+2,j,1],MX3[i+2,j,2]) - mismatch.score
					}
					if (x!=j){
						SC <- SC-penalty-abs(x-j)
					}
					if(SC > SCMax){
						SCMax <- SC
					}
				}
				MX3[i+1,j,z]=SCMax
			}
		}
	}
#	print(MX3)
	#add MX2 and MX3 store result in MX3
	MX3 <- MX2+MX3

	#make list of possible indels
	x <- c()
	for (i in 2:(seq.length+1)){
		SC <- -111111
		SCc <- 0
		for (j in 0:max.shift){
			if(is.null(fixed.shifts) || any(fixed.shifts==j)){
				if (MX3[i,j+1,1]> SC){
					SC <- MX3[i,j+1,1]
					SCc <- j
				}
				if (MX3[i,j+1,2]>SC){
					SC <- MX3[i,j+1,2]
					SCc <- j
				}
			}
		}
		if(!any(x==SCc)){
			x <- c(x,SCc)
		}
	}
#	print(x)
	#apply scores from MX3 and resolve positions in MX1. Results are stored in MX4.
	#MX4[1,i] the value indicating which base from MX1 goes into the first string. 
	#If MX4[1,i]==0 then the position is ambiguous.
	#MX4[2,i] the value indicating the phase shift which resolves the position.

	SCd <- 0
	MX4 <- matrix(0,nrow=4,ncol=(seq.length+1))
	repeat{
		for(i in 1:seq.length){
			SC <- -111111
			SCa <- 0	#possibility to solve in position 1
			SCb <- 0	#possibility to solve in position 2
			SCc <- 0	#shift size
			
			for(j in 0:max.shift){
				if((is.null(fixed.shifts) && SCd==0) || (any(fixed.shifts==j) && SCd==0) || 
					(any(x==j) && SCd==1)){
					if(MX3[i+1,j+1,1] > SC){
						SC <- MX3[i+1,j+1,1]
						SCa <- 0
						SCb <- 0
						SCc <- j
					}
					if(MX3[i+1,j+1,2] > SC){
						SC <- MX3[i+1,j+1,2]
						SCa <- 0
						SCb <- 0
						SCc <- j
					}
					if(((MX3[i+1,j+1,1]==SC) || (MX3[i+1,j+1,2]==SC)) && (0+j)==MX4[2,i]){
						SCc <- j
					}
					if((MX3[i+1,j+1,1]==SC) && (any(x==j))){
						SCa <- 1
					}
					if((MX3[i+1,j+1,2]==SC) && (any(x==j))){
						SCb <- 1
					}
				}
			}
		
			if(MX1[1,i]==MX1[2,i]){
				MX4[1,i+1] <- 1
			}else if (SCa==1 && SCb==0){
				MX4[1,i+1] <- 1
			}else if (SCb==1 && SCa ==0){
				MX4[1,i+1] <- 2
			}else{
				MX4[1,i+1] <- 0
			}
			MX4[2,i+1] <- SCc
		}
		#remove phase shifts recovered at the number of consecutive positions smaller than the phase shift magnitude

		
		u <- c()
		SCa <- 0

		for(i in 1:seq.length){
			if(MX4[2,i+1] != MX4[2,i]){
				for(z in 1:MX4[2,i+1]){
					if(i+z > seq.length){
						break
					}
					if(MX4[2,i+z+1]!=MX4[2,i+1]){
						SCa <- 1
						break
					}
					if (z==MX4[2,i+1]){
						z <- z+1
					}
				}
			}
			if(MX4[2,i+1] != MX4[2,2]){
				SCc <- 1
			}
			if(z>MX4[2,i+1] && !any(u==MX4[2,i+1])){
				u <- c(u,MX4[2,i+1])
			}
		}
		if(!all(u == x) && SCd==0){
			   SCa <- 1
		}
		x <- u
		if(SCa==1 || SCd==1){
			SCd <- SCd+1
		}
		if(SCd!=1){		#recalculate MX4 if short indel was removed
			break
		}
	}
	print(MX4)

	#resolve the remaining ambiguities.
	#MX5 is used to determine the possibility of resolving an ambiguous position with the 
	#same phase shift as the preceding positions or eith the same shift as the folowing positions
	MX5 <- matrix(,nrow=2,ncol=4)
	if (is.longindel){
		SCa = MX4[2,2]
		for(i in 1:seq.length){
			if(MX4[2,i+1]>SCa){
				for(j in 0:(SCa-1)){
					MX4[2,i+j+1] <- SCa
					MX4[1,i+j+1] <- 0
				}
				i <- i+SCa
				SCa <- MX4[2,i+1]
			}else if(MX4[2,i+1]<SCa){
				SCa <- MX4[2,i+1]
				for(j in 1:SCa){
					MX4[1,i-j+1] <- 0
					MX4[2,i-j+1] <- SCa
				}
			}
		}
		SCa <- MX4[2,2]
		SCd <- 0
		for(i in 1:seq.length){
			if(MX4[2,i+1] != SCa){
				for(j in max(1,i-10):min(seq.length,i+10)){
					if(MX3[j,SCa,2]==MX3[j,MX4[2,i+1],2] && MX3[j,SCa,2] >= MX3[j,SCa,1] && 
						MX3[j,MX4[2,i+1],2]>= MX3[j,MX4[2,i+1],1]){
						MX4[1,j+1] <- 0
					}else if(MX3[j,SCa,1] == MX3[j,MX4[2,i+1],1] && 
						MX3[j,SCa,2] <= MX3[j,SCa,1] && 
						MX3[j,MX4[2,i+1],2] <=MX3[j,MX4[2,i+1],1]){
						MX4[1,j] <- 0
					}
					if(MX1[1,j] == MX1[2,j]){
						MX4[1,j+1] <- 1
					}
				}
				if(SCa != 0){
					SCd <- abs(Scd-1)
				}
				SCa <- MX4[2,i+1]
			}
			if(SCd==1){
				MX4[2,i+1] <- -MX4[2,i+1]
				if(MX4[1,i+1] >0 && MX1[1,i] != MX1[2,i]){
					MX4[1,i+1]  <-  3-MX4[1,i+1]
				}
			}
		}
	}
	
	if(SCc==1){
		repeat{
			alignLeft()
			alignRight()

			#calculate and mark the ambiguities that could potentially be resolved

			for(i in 1:seq.length){
				if(MX4[1, i+1]==0){
					j <- abs(MX4[2,i+1])

					if(i<=j){
						MX4[3,i+1] <- 1	#1 could be resolved, 0- cannot be resolved
					}else if(j==0){
						MX4[3,i+1] <- 1
					}else if(MX4[3,i-j+1] && abs(MX4[2,i-j+1])==j && MX4[1,i-j+1]==0){
						MX4[3,i+1] <- 0
					}else if(MX4[3,i-j+1]==1 && abs(MX4[2,i-j+1])==j && MX4[1,i-j+1]==0){
						MX4[3,i+1] <- 1
					}else
					{
						SCc <- 1
						while(MX4[1,j*SCc+i+1]==0 && abs(MX4[2,j*SCc+1+i])==j && 
							j*SCc+i+1<=seq.length){
							SCc <- SCc+1
						}
						if(abs(MX4[2,i-j+1]) < j && MX4[2,i+1] >0){
							MX4[3,i+1] <- 1
						}else if(j*SCc +1 > seq.length){ # evtl +1
							MX4[3, i+1] <- 1
						}else if(abs(MX4[2,j*SCc+i+1]) != j){
							MX4[3,i+1] <- 1
						}else if(abs(MX4[2,j*SCc+i+1+MX4[3,i+1]]) != j){
							MX4[3,i+1] <- 1
						}else if(MX4[2,i+1] > 0 && (SCc / 2 )!=floor(SCc/2) && 
							((MX1[MX4[1,i-j+1],i-j] == MX1[2,i] && 
								MX1[1,i] == MX1[3-MX4[1,j*SCc*i+1],j*SCc+1]) || 
							(MX1[MX4[1,i-j+1],i-j]==MX1[1,i] && 
								MX1[2,i]==MX1[3-MX4[1,j*SCc+i+1],j*SCc+i]))){
							MX4[3,i+1] <- 1
						}else if(MX4[2,i+1] >0 && (SCc/2) == floor(SCc/2) && 
							((MX1[MX4[1,i-j+1],i-j]==MX1[1,i] && 
								MX1[2,i] == MX1[3-MX4[1,j*SCc+i+1],j*SCc+i])||
							(MX1[MX4[1,i-j+1],i-j] == MX1[1,i] && 
								MX1[1,i] == MX1[3- MX4[1,j*SCc+i+1],j*SCc+i]))){
							MX4[3,i+1] <- 1
						}else if(abs(MX4[2,i-j+1]) != j && MX4[2,i+1] < 0){
							MX4[3,i+1] <- 1
						}else if(MX4[1,i-j+1] ==0){
							MX4[3,i+1] <- 0
						}else if(MX4[2,i+1] < 0 && (SCc/2) != floor(SCc/2) && 
							((MX1[3-MX4[1,i-j+1],i-j] == MX1[1,i] && 
								MX1[2,i]== MX1[MX4[1,j*SCc+i+1],j*SCc+i]) || 
							(MX1[3-MX4[1,i-j+1],i-j] == MX1[2,i] && 
								MX1[1,i]==MX1[MX4[1,j*SCc+i+1],j*SCc+i]))){
							MX4[3,i+1] <- 1
						}else if(MX4[2,i+1] < 0 && (SCc/2) == floor(SCc/2) && 
							(( MX1[3-MX4[1,i-j+1],i-j]==MX1[1,i-j] && 
								MX1[1,i] == MX1[MX4[1,j*SCc+i+1],j*SCc+i]) || 
							(MX1[3-MX4[1,i-j+1],i-j] == MX1[2,i] &&
								MX1[2,i]==MX1[MX4[1,j*SCc+i+1],j*SCc+i]))){
							MX4[3,i+1] <- 1
						}else{
							MX4[3,i+1] <- 0
						}
					}
				}
			}
			j <- 0

			for(i in 1:seq.length){
				if(MX4[1,i+1]==0 && MX4[3,i+1]==1){
					for(z in 1:(max.shift+MX4[4,1])){
						if(z>=i || z> seq.length -i){
							break
						}
						if(MX4[2,i-z+1] != MX4[2,i+z+1]){
							break
						}
					}
					if(z<=i && z-MX4[3,i+1] <= max.shift+1 && z <= seq.length -(i+1)){
						ind1 <- abs(MX4[2,i-z+1])
						ind2 <- abs(MX4[2,i+z+1])
						if(!is.longindel){
							if((i>= ind1) && (i>= ind2) && (seq.length -i)>ind1 && 
								(seq.length-i)>ind2){
								if(((MX1[2,i]==MX1[1,i-ind1] && 
									MX3[i-ind+1,ind1,1]>=MX3[i-ind1+1,ind1,2]) || 
								(MX1[2,i] == MX1[2,i-ind1] && 
									MX3[i-ind1+1,ind1,1] <= MX3[i-ind1+1,ind1,2])) && 
								((MX1[1,i] == MX1[2,i+ind1] && 
									MX3[i+ ind1+1,ind1,1] >= MX3[i+ind1+1,ind1,2])||
								(MX1[1,i] == MX1[1,i+ind1] && 
									MX3[i+ind1+1,ind1,1] <= MX3[i+ind1+1,ind1,2]))){
									MX5[1,1] <- 1
								}else{
									MX5[1,1] <- 0
								}
								if(((MX1[1,i]==MX1[1,i-ind1] && 
									MX3[i-ind1+1,ind1,1]>=MX3[i-ind1+1,ind1,2])||
								(MX1[1,i] == MX1[2,i-ind] &&
									MX3[i-ind1+1,ind1,1] <= MX3[i-ind1+1,ind1,2]))&&
								((MX1[2,i]==MX1[2,i+ind1] && 
									MX3[i+ind1+1,ind1,1]>=MX3[i+ind1+1,ind1,2])||
								(MX1[2,i]==MX1[1,i+ind1] && 
									MX3[i+ind1+1,ind1,1]<=MX3[i+ind1+1,ind1,2]))){
									MX5[1,2] <- 1
								}else{
									MX5[1,2] <- 0
								}
								if(((MX1[2,i]==MX1[1,i-ind2]&&
									MX3[i-ind2+1,ind2,1]>=MX3[i-ind2+1,ind2,2])||
								(MX1[2,i]==MX1[2,i-ind2]&&
									MX3[i-ind2+1,ind2,1]<= MX3[i-ind2+1,ind2,2]))&&
								((MX1[1,i]==MX1[2,i+ind2]&&
									MX3[i+ind2+1,ind2,1]>=MX3[i+ind2+1,ind2,2])||
								(MX1[1,i]==MX1[1,i+ind2]&&
									MX3[i+ind2+1,ind2,1]<=MX3[i+ind2+1,ind2,2]))){
									MX5[2,1] <- 1
								}else{
									MX5[2,1] <- 0
								}
								if(((MX1[1,i]==MX1[1,i-ind2] && 
									MX3[i-ind2+1,ind2,1]>=MX3[i-ind2+1,ind2,2])||
								(MX1[1,i]==MX1[2,i-ind2]&&
									MX3[i-ind2+1,ind2,1]<=MX3[i-ind2+1,ind2,2]))&&
								((MX1[2,i]==MX1[2,i+ind2]&&
									MX3[i+ind2+1,ind2,1]>=MX3[i+ind2+1,ind2,2])||
								(MX1[2,i]==MX1[1,i+ind2]&&
									MX3[i+ind2+1,ind2,1]<=MX3[i+ind2+1,ind2,2]))){
									MX5[2,2] <- 1
								}else{
									MX5[2,2] <- 0
								}
								if((ind1 < ins2) &&((MX1[2,i]==MX1[1,i-ind2]&&MX3[i-ind1+1,ind1,1]>=MX3[i-ind1+1,ind1,2])||(MX1[2,i]==MX1[2,i-ind1]&&MX3[i-ind1+1,ind1,1]<=MX3[i-ind1+1,ind1,2]))&&((MX1[1,i]==MX1[2,i+ind2]&&MX3[i+ind2+1,ind2,1]>=MX3[i+ind2+1,ind2,2])||(MX1[1,i]==MX1[1,i+ind2]&&MX3[i+ind2+1,ind2,1]<=MX3[i+ind2+1,ind2,2]))){
									MX5[4,1] <- 1
								}else{
									MX5[4,1] <- 0
								}
								if((ind1<ind2)&&((MX1[1,i]==MX1[1,i-ind]&&MX3[i-ind1+1,ind1,1]>=MX3[i-ind1+1,ind1,2])||(MX1[1,i]==MX1[2,i-ind1]&&MX3[i-ind1+1,ind1,1]<=MX3[i-ind1+1,ind1,2]))&&((MX1[2,i]==MX1[2,i+ind2]&&MX3[i+ind2+1,ind2,1]>=MX3[i+ind2+1,ind2,2])||(MX1[2,i]==MX1[1,i+ind2]&&MX3[i+ind2+1,ind2,1]<=MX3[i+ind2+1,ind2,2]))){
									MX5[4,2] <- 1
								}else{
									MX5[4,1] <- 0
								}
							}else if(i<=ind1){
								if(((MX1[1,i]==MX1[2,i+ind1]&&MX3[i+ind1+1,ind1,1]>=MX3[i+ind1+1,ind1,2])||(MX1[1,i]==MX1[1,i+ind1]&&MX3[i+ind1+1,ind1,1]<=MX3[i+ind1+1,ind1,2]))){
									MX5[1,1] <- 1
								}else{
									MX5[1,1] <- 0
								}
								if(((MX1[2,i]==MX1[2,i+ind1]&&MX3[i+ind1+1,ind1,1]>=MX3[i+ind1+1,ind1,2])||(MX1[2,i]==MX1[1,i+ind1]&&MX3[i+ind1+1,ind1,1]<=MX3[i+ind1+1,ind1,2]))){
									MX5[1,2] <- 1
								}else{
									MX5[1,2] <- 0
								}
								MX5[2,1] <- 0
								MX5[2,2] <- 0
							}else if(seq.length-i< ind2){
								MX5[1,1] <- 0
								MX5[1,2] <- 0
								if(((MX1[2,i]==MX1[1,i-ind2]&&MX3[i-ind2+1,ind2,1]>=MX3[i-ind2+1,ind2,2])||(MX1[2,i]==MX1[2,i-ind2]&&MX3[i-ind2+1,ind2,1]<=MX3[i.ind2+1,ind2,2]))){
									MX5[2,1] <- 1
								}else{
									MX5[2,1] <- 0
								}
								if(((MX1[1,i]==MX1[1,i-ind2]&&MX3[i-ind2+1,ind2,1]>=MX3[i-ind2+1,ind2,2])||(MX1[1,i]==MX1[2,i-ind2]&&MX3[i-ind2+1,ind2,1]<=MX3[i-ind2+1,ind2,2]))){
									MX5[2,2] <- 1
								}else{
									MX5[2,2] <- 0
								}
							}
							if(MX4[2,i-z+1]<MX4[2,i+z+1]||i<=ind1){
								if(((MX1[1,i]==MX1[2,i+ind2]&&MX3[i+ind2+1,ind2,1]>=MX3[i+ind2+1,ind2,2])||(MX1[1,i]==MX1[1,i+ind2]&&MX3[i+ind2+1,ind2,1]<=MX3[i+ind2+1,ind2,2]))){
									MX5[3,1] <- 1
								}else{
									MX5[3,1] <- 0
								}
								if(((MX1[2,i]==MX1[1,i+ind2]&&MX3[i+ind2+1,ind2,1]>=MX3[i+ind2+1,ind2,2])||(MX1[2,i]==MX1[1,i+ind2]&&MX3[i+ind2+1,ind2,1]<=MX3[i+ind2+1,ind2,2]))){
									MX5[3,2] <- 1
								}else{
									MX5[3,2] <- 0
								}
							}else{
								if(((MX1[2,i]==MX1[1,i-ind1]&&MX3[i-ind1+1,ind1,1]>=MX3[i-ind1+1,ind1,2])||(MX1[2,i]==MX1[2,i-ind1]&&MX3[i-ind1+1,ind1,1]<=MX3[i-ind1+1,ind1,2]))){
									MX5[3,1] <- 1
								}else{
									MX5[3,1] <- 0
								}
								if(((MX1[1,i]==MX1[1,i-ind1]&&MX3[i-ind1+1,ind1,1]>=MX3[i-ind1+1,ind1,2])||(MX1[2,i]==MX1[2,i-ind1]&&MX3[i-ind1+1,ind1,1]<=MX3[i-ind1+1,ind1,2]))){
									MX5[3,2] <- 1
								}else{
									MX5[3,2] <- 0
								}
							}
							if(ind1==0||ind2==0){
								MX5[4,1] <- 0
								MX5[4,2] <- 0
							}	
						}else{
							if(MX4[2,i-z+1]>0){
								if(((MX1[2,i]==MX1[1,i-ind1]&&MX4[1,i-ind1+1] !=2)||(MX1[2,i]==MX1[2,i-ind2]&&MX4[1,i-ind1+1]!=1))){
									MX5[3,1] <- 1
								}else{
									MX5[3,1] <- 0
								}
								if(((MX1[i,1]==MX1[1,i-ind1]&&MX4[1,i-ind1+1]!=2)||(MX1[1,i]==MX1[2,i-ind1]&&MX4[1,i-ind1]!=1))){
									MX5[3,2] <- 1
								}else{
									MX5[3,2] <- 0
								}
							}else{
								if(((MX1[1,i]==MX1[2,i-ind1] && MX4[1,i-ind1+1]!=2)||(MX1[1,i]==MX1[1,i-ind1]&&MX4[1,i-ind1+1]!=1))){
									MX5[4,1] <- 1
								}else{
									MX5[4,1] <- 0
								}
								if(((MX1[2,i]==MX1[2,i-ind1]&&MX4[1,i-ind1+1]!=2)||(MX1[2,i]==MX1[1,i-ind1]&&MX4[1,i-ind1+1]!=1))){
									MX5[4,2] <- 1
								}else{
									MX5[4,2] <- 0
								}
							}
							if(MX4[2,i+z+1]>0){
								if(((MX1[1,i]==MX1[2,i+ind2]&&MX4[1,i+ind2+1]!=2)||(MX1[1,i]==MX1[1,i+ind2]&&MX4[1,i+ind2+2]!=1))){
									MX5[4,1] <- 1
								}else{
									MX5[4,1] <- 0
								}
								if(((MX1[2,i]==MX1[2,i+ind2]&&MX4[1,i+ind2+1]!=2)||(MX1[2,i]==MX1[1,i+ind2]&&MX4[1,i+ind2+1]!=1))){
									MX5[4,2] <- 1
								}else{
									MX5[4,2] <- 0
								}
							}else{
								if(((MX1[2,i]==MX1[1,i+ind2]&&MX4[1,i+ind2+1]!=2)||(MX1[2,i]==MX1[2,i+ind2]&&MX4[1,i+ind2+1]!=1))){
									MX5[4,1] <- 1
								}else{
									MX5[4,1] <- 0
								}
								if(((MX1[1,i]==MX1[1,i+ind2]&&MX4[1,i+ind2+1]!=2)||(MX1[1,i]==MX1[2,i+ind2]&&MX4[1,i+ind2+1]!=1))){
									MX5[4,2] <- 1
								}else{
									MX5[4,2] <- 0
								}
							}
							if(i >= ind1 && i>= ind2 && (seq.length-i)>ind1&&(seq.length-i)<ind2){
								if(MX4[2,i-z+1]>0){
									if(((MX1[2,i]==MX1[1,i-ind1]&&MX4[1,i-ind1+1]!=2)||(MX1[2,i]==MX1[2,i-ind1]&&MX4[1,i-ind1+1]!=1))&&((MX1[1,i]==MX1[2,i+ind1]&&MX4[1,i+ind1+1]!=2)||(MX1[1,i]==MX1[1,i+ind1]&&MX4[1,i+ind1+1]!=1))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((MX1[1,i]==MX1[1,i-ind1]&&MX4[1,i-ind1+1]!=2)||(MX1[1,i]==MX1[2,i-ind1]&&MX4[1,i-ind1]!1))&&((MX1[2,i]==MX1[2,i+ind1]&&MX4[1,i+ind1+1]!=2)||(MX1[2,i]==MX1[1,i+ind1]&&MX4[1,i+ind1+1]!=1))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
								}else{
									if(((MX1[1,i]==MX1[2,i-ind1]&&MX4[1,i-ind1]!=2)||(MX1[1,i]==MX1[1,i-ind1]&&MX4[1,i-ind1+1]!=1))&&((MX1[2,i]==MX1[1,i+ind1]&&MX4[1,i+ind1+1]!=2)||(MX1[2,i]==MX1[1,i+ind1]&&MX4[1,i+ind1+1]!=1))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((MX1[2,i]==MX1[2,i-ind1]&&MX4[1,i-ind1+1]!=2)||(MX1[2,i]==MX1[1,i-ind1]&&MX4[1,i-ind1+1]!=1))&&((MX1[1,i]==MX1[1,i+ind1]&&MX4[1,i+ind1]!=2)||(MX1[1,i]==MX1[2,i+ind1]&&MX4[1,i+ind1+1]!=1))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
								}else if(i<=ind1){
									if(((MX1[1,i]==MX1[2,i+ind1]&&MX4[1,i+ind1+1]!=2)||(MX1[1,i]==MX1[1,i+ind1]&&MX4[1,i+ind1+1]!=1))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((MX1[2,i]==MX1[2,i+ind1]&&MX4[1,i+ind1]!=2)||(MX1[2,i]==MX1[1,i+ind1]&&MX4[1,i+ind1+1]))){
										MX5[1,2] <- 1
									}else{
										MX5[1,2] <- 0
									}
									MX5[2,] <- 0
								}else if((seq.length-i)<ind2){
									MX5[1,] <- 0
									if(((MX1[2,i]==MX1[1,i-ind2]&&MX4[1,i-ind1]!=2)||(MX1[2,i]==MX1[2,i-ind2]&&MX4[1,i-ind1]!=1))){
										MX5[2,1] <- 1
									}else{
										MX5[2,1] <- 0
									}
									if(((MX1[1,i]==MX1[1,i-ind2]&&MX4[1,i-ind1]!=2)||(MX1[1,i]==MX1[2,i-ind2]&&MX4[1,i-ind1]!=1))){
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
									if(MX4[2,i]==ind1){
										if(z<=ind1){
											MX5[1,] <- 0
											if(z<=ind2){
												MX5[3,] <- 0
											}
										}

										if((MX4[4,i+1]-z)<0){
											MX5[2,] <- 0
										}
									}else{
										MX5[1,] <- 0
										MX5[3,] <- 0
									}
								}else{ #ind1<ind2
									if(MX4[2,i+1]==ind1){
										if(MX4[3,i+1]<=ind2){
											MX5[2,] <- 0
											if(ind1 != 0 && z > MX4[4,i+1]){
												MX5[3,] <- 0
											}
										}
										if(z>ind2){
											MX5[2,] <- 0
										}
									}else{
										if(z+MX4[4,i+1]<=ind2){
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
									if(abs(MX4[2,i+1])==ind1){
										if(z<=ind1){
											MX5[1,] <- 0
										}
										if((MX4[4,i+1]-z)<0){
											MX5[c(2,4),] <- 0
										}
									}else{
										MX5[c(1,3),] <- 0
									}
								}else{
									if(MX4[2,i+1] <0 ind2){
										if(MX4[4,i+1]<=ind2){
											MX5[2,] <- 0
										}
										if(ind1!=0 &&z > MX4[4,i+1]){
											MX5[4,] <- 0
										}
										if(z>ind2){
											MX5[2,] <- 0
										}
									}else{
										if((z+MX4[4,i+1])<=ind2){
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
									MX4[1,i+1] <- 1
									MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1
									MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
								}else if(all(MX5[c(1,2),1]==0)&&all(MX5[c(1,2),2]==1)){
									MX4[1,i+1] <- 2
									MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1
									MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
								}else if((MX5[1,1]==1&&MX5[2,2]==1)||(MX5[1,2]==1&&MX5[2,1]==1)){
								}else if(all(MX5[c(1,3),1]==1)&&all(MX5[c(1,3),2]==0)){
									MX4[1,i+1] <- 1
									MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1
									MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
								}else if(all(MX5[c(1,3),1]==0)&&all(MX5[c(1,3),2]==1))){
									MX4[1,i+1] <- 2
									MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1
									MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
								}else if(((MX5[1,1]==1&&MX5[3,2]==1)||(MX5[1,2]==1&&MX5[3,1]==1))
									&&z< ind1 &&MX4[4,i+1]>ind1{
								}else if(((MX5[2,1]==1&&MX5[3,2]==1)||(MX5[2,2]==1&&MX5[3,1]==1))
									&&z<=ind1 &&MX4[4,i+1]>ind2 &&ind2 >ind1){
								}else if(((MX5[2,1]==1&&MX5[3,2]==1)||(MX5[2,2]==1&&MX5[3,1]==1))
									&&MX4[2,i+1]==ind2&&MX4[1,i+ind2+1]==0){
								}else if(((MX5[2,1]==1&&MX5[3,2]==1)||(MX5[2,2]==1&&MX5[3,1]==1))
									&&MX4[2,i+1]==ind1&&MX4[4,i+1]>=z){
								}else if(((MX5[1,1]==1&&MX5[3,2]==1)||(MX5[1,2]==1&&MX5[3,1]==1))
									&&MX4[2,i+1]==ind1&&ind2>ind1&&MX4[4,i+1]>=(z+ind1)){
								}else if(xor(MX5[1,1]==1,MX5[1,2]==1){
									if(MX5[1,1]==1){
										MX4[1,i+1] <- 1
										MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1 
										MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
									}else{
										MX4[1,i+1] <- 2
										MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1 
										MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
									}
								}else if(xor(MX5[2,1]==1,MX5[2,2]==1)){
									if(MX5[2,1]==1){
										MX4[1,i+1] <- 1
										MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1 
										MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
									}else{
										MX4[1,i+1] <- 2
										MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1 
										MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
									}
								}else if(all(MX5[3,]==1)&&all(MX5[4,]==0)){
								}else if(all(MX5[3,]==1)&&all(MX5[4,]==c(1,0))){
									MX4[1,i+1] <- 1
									MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1 
									MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
								}else if(all(MX5[3,]==1)&&all(MX5[4,]==c(0,1))){
									MX4[1,i+1] <- 2
									MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1 
									MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
								}else if(MX5[3,1]==1){
									MX4[1,i+1] <- 1
									MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1 
									MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
								}else if(MX5[3,2]==1){
									MX4[1,i+1] <- 2
									MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1 
									MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
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
										MX4[1,i+1] <- 1
										MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1 
										MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
									}else{
										MX4[1,i+1] <- 2
										MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1 
										MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
									}
								}else if(xor(MX5[2,1]==1,MX5[2,2]==1)){
									if(MX5[2,1]==1){
										MX4[1,i+1] <- 1
										MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1 
										MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
									}else{
										MX4[1,i+1] <- 2
										MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1 
										MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
									}
								}else if(all(MX5[3,]==1)){
								}else if(all(MX5[4,]==1)){
								}else if(all(MX5[4,]==c(1,0))){
									X4[1,i+1] <- 1
									MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1 
									MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
								}else if(all(MX5[4,]==c(0,1))){
									MX4[1,i+1] <- 2
									MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1 
									MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
								}else if(all(MX5[3,]==c(1,0))){
									X4[1,i+1] <- 1
									MX3[i+1,ind2,1] <- max(MX3[i+1,ind2,])+1 
									MX3[i+1,ind1,1] <- max(MX3[i+1,ind1,])+1
								}else if(all(MX5[3,]==c(0,1))){
									MX4[1,i+1] <- 2
									MX3[i+1,ind2,2] <- max(MX3[i+1,ind2,])+1 
									MX3[i+1,ind1,2] <- max(MX3[i+1,ind1,])+1
								}
							}
						}else{
							if(MX4[2,i]==MX4[2,i+2]||i==1||i==seq.length){
								if(i!=seq.length){
									ind1 <- abs(MX4[2,i+2])
									ind2 <- MX4[2,i+2]
								}else{
									ind1 <- abs(MX4[2,i])
									ind2 <- MX4[2,i]
								}
								if(i<=ind1){
									if(((MX1[1,i]==MX1[2,i+ind1]&&
										MX3[i+ind1+1,ind1,1] >= MX3[i+ind1+1,ind1,2])||
									(MX1[1,i]==MX1[1,i+ind1]&&
										MX3[i+ind1+1,ind1,1] <= MX3[i+ind1+1,ind1,2]))){
										MX5[1,1] <- 1
									}else{
										MX5[1,1] <- 0
									}
									if(((MX1[2,i]==MX1[2,i+ind1]&& )))
								}
							}
						}
					}
				}
			}

		}
	}
 

}

		       
	       			       




