calculate.forward <-
function(seq.matrix,max.shift,penalty){
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
