calculate.phase.shift <-
function(scores,sequences,max.shift,fixed.shifts,possible.indels){
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
