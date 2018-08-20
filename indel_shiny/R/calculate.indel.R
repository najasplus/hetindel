calculate.indel <-
function(sequence, max.shift=15){

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
