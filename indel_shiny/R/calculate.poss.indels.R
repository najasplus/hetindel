calculate.poss.indels <-
function(combined.matrix,seq.length,max.shift,fixed.shifts){
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
