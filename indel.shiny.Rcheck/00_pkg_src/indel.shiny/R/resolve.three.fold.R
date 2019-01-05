resolve.three.fold <-
function(bd, seq.matrix, phase.shift.matrix){

#resolve 3-fold degenerate bases
	seq.length <- length(seq.matrix[1,])
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
