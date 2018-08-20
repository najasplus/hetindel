calculate.backward <-
function(seq.matrix,max.shift,penalty){

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
