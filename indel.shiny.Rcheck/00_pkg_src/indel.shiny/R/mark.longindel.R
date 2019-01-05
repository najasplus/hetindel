mark.longindel <-
function(is.longindel,sequences,scores,shifts){
	#mark and rotate positions for long indels 
	seq.length <- length(sequences[1,])
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
					SCd <- abs(SCd-1)
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
