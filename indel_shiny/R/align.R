align <-
function(is.right,is.align,phase.shift.matrix,seq.matrix){
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
