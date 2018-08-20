seq.toMatrix <-
function(sequence,max.shift){

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
