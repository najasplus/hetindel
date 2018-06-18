library(crayon)

set.seed(0)

randomBase <- function(dist=runif(1,0,4)){
	randVar <- dist
	if(randVar<1){
		return("A")
	}else if(randVar<2){
		return("C")
	}else if(randVar<3){
		return("G")
	}else{
		return("T")
	}
}

generateSeq <- function(length,dist=runif(1,0,4)){
	seqence <- c()
	for (i in 1:length) {
		seqence <- c(seqence,randomBase())
	}
	return(seqence)
}

iupacNotation <- function(chars){
	if(any(chars[]==" ")&&any(chars[]=="A")){
		return("A")
	}else if(any(chars[]==" ")&&any(chars[]=="C")){
		return("C")
	}else if(any(chars[]==" ")&&any(chars[]=="G")){
		return("G")
	}else if(any(chars[]==" ")&&any(chars[]=="T")){
		return("T")
	}else if(any(chars[]=="A")&&any(chars[]=="G")){
		return("R")
	}else if(any(chars[]=="C")&&any(chars[]=="T")){
		return("Y")
	}else if(any(chars[]=="G")&&any(chars[]=="T")){
		return("K")
	}else if(any(chars[]=="A")&&any(chars[]=="C")){
		return("M")
	}else if(any(chars[]=="C")&&any(chars[]=="G")){
		return("S")
	}else if(any(chars[]=="A")&&any(chars[]=="T")){
		return("W")
	}else{
		if(chars[1]=="C"){
			return("C")
		}else if(chars[1]=="A"){
			return("A")
		}else if(chars[1]=="T"){
			return("T")
		}else {
			return("G")
		}
	}
}

combineSeq <- function(seq1,seq2){
	result <- rbind(seq1,seq2)
	result <- apply(result,2,iupacNotation)
	return(result)
}


introduceShift <- function(seqence,exact=TRUE,occurrence="m"){
	
	length <- length(seqence)
	if(length==150){
		lower <- 0
		upper <- 10
	}else if(length==300){
		lower <- 10
		upper <- 20
	}else if(length==700){
		lower <- 35
		upper <- 45
	}

	shiftLength <- floor(runif(1,lower,upper+1))
	if(exact){
		shiftStart <- 0
	}else{
	shiftStart <- floor(runif(1,-shiftLength,shiftLength))
	}

	print(shiftLength)
	shiftedSeq <- character(length+shiftLength)
	shiftedSeq[] <- " "
	if(occurrence=="b"){
		shiftedSeq[(shiftLength+1):(length+shiftLength)] <- seqence[1:length]	
	}else if(occurrence=="m"){
		shiftStart <- floor(length/2) + shiftStart
		shiftedSeq[1:shiftStart] <- seqence[1:shiftStart]
		shiftedSeq[(shiftStart+1+shiftLength):(length+shiftLength)] <- seqence[(shiftStart+1):length]
	}else if(occurrence=="q1"){
		shiftStart <- floor(length/4) + shiftStart
		shiftedSeq[1:shiftStart] <- seqence[1:shiftStart]
		shiftedSeq[(shiftStart+shiftLength+1):(length+shiftLength)] <- seqence[(shiftStart+1):length]
	}else if(occurrence=="q3"){
		shiftStart <- floor(3*length/4) + shiftStart
		shiftedSeq[1:shiftStart] <- seqence[1:shiftStart]
		shiftedSeq[(shiftStart+shiftLength+1):(length+shiftLength)] <- seqence[(shiftStart+1):length]
	}
		print(shiftedSeq)
		insert <- character(shiftLength)
		insert[] <- " "
		seqence <- c(seqence, insert)
		print(seqence)
		combinedSeq <- combineSeq(seqence,shiftedSeq)
	return(combinedSeq)


}

generateRandBioString <- function(length,dist=runif(1,0,4),exact=TRUE,occurrence="m"){
	seqence <- generateSeq(length,dist)
	seqence <- introduceShift(seqence,exact,occurrence)
	return(seqence)
}