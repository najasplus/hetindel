calculate.scores <-
function(seq.matrix,max.shift,penalty){
	forward <- calculate.forward(seq.matrix,max.shift,penalty)
	backward <- calculate.backward(seq.matrix,max.shift,penalty)
	return(forward+backward)
}
