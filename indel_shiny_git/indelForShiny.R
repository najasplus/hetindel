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

	seq.length <- length(sequences[1,])
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
					(any(possible.indels==j) && removed.short.indel==1)){
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
					if((scores[i+1,j+1,1]==score) && (any(possible.indels==j))){
						pos.one <- 1
					}
					if((scores[i+1,j+1,2]==score) && (any(possible.indels==j))){
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
		z <- 0
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
		if(!all(u == possible.indels) && removed.short.indel==0){
			   pos.one <- 1
		}
		possible.indels <- u
		if(pos.one==1 || removed.short.indel==1){
			removed.short.indel <- removed.short.indel+1
		}
		if(removed.short.indel!=1){		#recalculate result if short indel was removed
			break
		}
	}

	return(list(result=result,shift.score=shift))
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

	return(shifts)
}

resolve.all <- function(phase.shift.matrix, seq.matrix, scores, shift.score, max.shift, is.longindel){

	resolved <- matrix(,nrow=4,ncol=2)

	seq.length <- length(seq.matrix[1,])
	if(shift.score==1){
		repeat{
			phase.shift.matrix <- alignLeft(phase.shift.matrix,seq.length,seq.matrix)$matr
			phase.shift.matrix <- alignRight(phase.shift.matrix,seq.length,seq.matrix)$matr

			#print(phase.shift.matrix)
			#calculate and mark the ambiguities that could potentially be resolved




			for(i in 1:seq.length){
				if(phase.shift.matrix[1, i+1]==0){
					j <- abs(phase.shift.matrix[2,i+1])


					#print(j*shift.score+i+1)
					#print(i)
					if(i<=j){
						phase.shift.matrix[3,i+1] <- 1	#1 could be resolved, 0- cannot be resolved
					}else if(j==0){
						phase.shift.matrix[3,i+1] <- 1
					}else if(phase.shift.matrix[3,i-j+1] == 0 && abs(phase.shift.matrix[2,i-j+1])==j && phase.shift.matrix[1,i-j+1]==0){
						phase.shift.matrix[3,i+1] <- 0
					}else if(phase.shift.matrix[3,i-j+1]==1 && abs(phase.shift.matrix[2,i-j+1])==j && phase.shift.matrix[1,i-j+1]==0){
						phase.shift.matrix[3,i+1] <- 1
					}else if(j*shift.score+i+1 <=seq.length){
						shift.score <- 1
						while(phase.shift.matrix[1,j*shift.score+i+1]==0 && abs(phase.shift.matrix[2,j*shift.score+1+i])==j && 
							j*shift.score+i+1<=seq.length){
							shift.score <- shift.score+1
						}
						#cat("erste Stelle",3-phase.shift.matrix[1,j*shift.score+i+1],"\n")
						#cat("zweite Stelle",j*shift.score+i,"\n")
						#cat("Wert",seq.matrix[3-phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i],"\n")


						if(abs(phase.shift.matrix[2,i-j+1]) < j && phase.shift.matrix[2,i+1] >0){
							phase.shift.matrix[3,i+1] <- 1
						}else if((j*shift.score +i) > seq.length){ # evtl +1
							phase.shift.matrix[3, i+1] <- 1
						}else if(abs(phase.shift.matrix[2,j*shift.score+i+1]) != j){
							phase.shift.matrix[3,i+1] <- 1
						}else if(abs(phase.shift.matrix[2,j*shift.score+i+1+phase.shift.matrix[4,i+1]]) != j){
							phase.shift.matrix[3,i+1] <- 1


						}else if(phase.shift.matrix[1,i-j+1]>0 && (3-phase.shift.matrix[1,j*shift.score+i+1]<3 )&& 
							(3-phase.shift.matrix[1,j*shift.score+i+1])>0){

							if((phase.shift.matrix[2,i+1] > 0) && ((shift.score / 2 ) != floor(shift.score/2)) && 
							((seq.matrix[phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] && 
									seq.matrix[1,i] == seq.matrix[3-phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]) || 
								(seq.matrix[phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i] && 
									seq.matrix[2,i]==seq.matrix[3-phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]))){
								

								phase.shift.matrix[3,i+1] <- 1
							}
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*shift.score+i+1])<3 && 
							3-phase.shift.matrix[1,j*shift.score+i+1]>0){

							if(phase.shift.matrix[2,i+1] >0 && (shift.score/2) == floor(shift.score/2) && 
								((seq.matrix[phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i] && 
									seq.matrix[2,i] == seq.matrix[3-phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i])||
								(seq.matrix[phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[1,i] && 
									seq.matrix[1,i] == seq.matrix[3- phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]))){
								phase.shift.matrix[3,i+1] <- 1
							}
						}else if(abs(phase.shift.matrix[2,i-j+1]) != j && phase.shift.matrix[2,i+1] < 0){
							phase.shift.matrix[3,i+1] <- 1
						}else if(phase.shift.matrix[1,i-j+1] ==0){
							phase.shift.matrix[3,i+1] <- 0
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*shift.score+i+1])<3 && 
							3-phase.shift.matrix[1,j*shift.score+i+1]>0){

							if(phase.shift.matrix[2,i+1] < 0 && (shift.score/2) != floor(shift.score/2) && 
								((seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[1,i] && 
									seq.matrix[2,i]== seq.matrix[phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]) || 
								(seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] && 
									seq.matrix[1,i]==seq.matrix[phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]))){
								phase.shift.matrix[3,i+1] <- 1

						}
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*shift.score+i+1])<3 && 
							3-phase.shift.matrix[1,j*shift.score+i+1]>0){

							if(phase.shift.matrix[2,i+1] < 0 && (shift.score/2) == floor(shift.score/2) && 
								(( seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i-j] && 
									seq.matrix[1,i] == seq.matrix[phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]) || 
								(seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] &&
									seq.matrix[2,i]==seq.matrix[phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]))){
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
										scores[i-ind1+1,ind1+11,1]>=scores[i-ind1+1,ind1+1,2])||
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
								######################################################
							}else if(i<=ind1){
									if(((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]))){
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
#############################################################################################
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
								X4[1,i+1] <- 1
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
	return(list(seq.matrix=seq.matrix, phase.shift.matrix=phase.shift.matrix, resolved=resolved, scores=scores))

}

align <- function(is.right,is.align,phase.shift.matrix,seq.matrix){
	seq.length <- length(seq.matrix[1,])
	if(!is.right&&is.align){
		phase.shift.matrix<- alignLeft(phase.shift.matrix,seq.length,seq.matrix)$matr
		phase.shift.matrix <- alignRight(phase.shift.matrix,seq.length,seq.matrix)$matr
	}else if(is.right&&is.align){
		phase.shift.matrix <- alignRight(phase.shift.matrix,seq.length,seq.matrix)$matr
		phase.shift.matrix <- alignLeft(phase.shift.matrix,seq.length,seq.matrix)$matr
	}
	return(phase.shift.matrix)
}

resolve.three.fold <- function(bd, seq.matrix, phase.shift.matrix){

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

	return(seq.matrix)
}

reconstruct <- function(seq.matrix, phase.shift.matrix){

	#based on phase.shift.matrix and seq.matrix reconstruct two allelic sequences
	seq.length <- length(seq.matrix[1,])

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
			temp5 <- temp5 %+% red(SCc)
		}
	}

	return(list(seq1=temp1, seq2=temp2, combined=temp5))


}

calculate.indel <- function(sequence, max.shift=15){

	temp <- seq.toMatrix(sequence, max.shift)
	seq.matrix <- temp$matrix
	is.three.fold <- temp$is.three.fold
	ambig <- temp$ambig
	shift <- temp$shift

	penalty <- 2
	scores <- calculate.scores(seq.matrix, shift, penalty)

	fixed.shifts <- NULL

	poss.indel <- calculate.poss.indels(scores,length(sequence), shift, fixed.shifts)

	phase.shift <- calculate.phase.shift(scores, seq.matrix, shift, fixed.shifts, poss.indel)
	phase.shift.matrix <- phase.shift$result
	shift.score <- phase.shift$shift.score

	resolve <- resolve.all(phase.shift.matrix, seq.matrix, scores, shift.score, shift, FALSE)

	seq.matrix <- resolve$seq.matrix
	phase.shift.matrix <- resolve$phase.shift.matrix
	resolved <- resolve$resolved
	scores <- resolve$scores
	
	is.right <- FALSE
	is.align <- TRUE

	phase.shift.matrix <- align(is.right,is.align,phase.shift.matrix, seq.matrix)
	
	seq.matrix <- resolve.three.fold(is.three.fold, seq.matrix, phase.shift.matrix)


	reconstructed <- reconstruct(seq.matrix, phase.shift.matrix)

	seq1 <- reconstructed$seq1
	seq2 <- reconstructed$seq2
	combined <- reconstructed$combined

	return(list(seq1=seq1,seq2=seq2))
	#cat("Reconstructed Sequences:\n",seq1,"\n",seq2,"\nCombined Sequence:\n",combined,"\n")
}