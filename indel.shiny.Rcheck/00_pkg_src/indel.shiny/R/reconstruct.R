reconstruct <-
function(seq.matrix, phase.shift.matrix){

	#based on phase.shift.matrix and seq.matrix reconstruct two allelic sequences
	seq.length <- length(seq.matrix[1,])

	temp1 <- c()
	temp2 <- c()
	temp3 <- c()
	ambig2 <- 0
	for(i in 1:seq.length){
		if(phase.shift.matrix[1,i+1]==1){
			temp1 <- c(temp1,seq.matrix[1,i])
			temp2 <- c(temp2,seq.matrix[2,i])
		}else if(phase.shift.matrix[1,i+1] == 2){
			temp1 <- c(temp1, seq.matrix[2,i])
			temp2 <- c(temp2, seq.matrix[1,i])
		}else{
			ambig2 <- ambig2 +1
			if(seq.matrix[1,i]=="A" && seq.matrix[2,i]=="G"){
				temp1 <- c(temp1,"R")
				temp2 <- c(temp2,"R")
			}else if(seq.matrix[1,i]=="C" && seq.matrix[2,i]=="T"){
				temp1 <- c(temp1,"Y")
				temp2 <- c(temp2,"Y")
			}else if(seq.matrix[1,i]=="G" && seq.matrix[2,i]=="C"){
				temp1 <- c(temp1,"S")
				temp2 <- c(temp2,"S")
			}else if(seq.matrix[1,i]=="A" && seq.matrix[2,i]=="T"){
				temp1 <- c(temp1,"W")
				temp2 <- c(temp2,"W")
			}else if(seq.matrix[1,i]=="A" && seq.matrix[2,i]=="C"){
				temp1 <- c(temp1,"M")
				temp2 <- c(temp2,"M")
			}else if(seq.matrix[1,i]=="G" && seq.matrix[2,i]=="T"){
				temp1 <- c(temp1,"K")
				temp2 <- c(temp2,"K")
			}
		}
		temp3 <- c(temp3,phase.shift.matrix[2,i+1])
	}


	#Align reconstructed allelic sequences (indicate gaps with dots)

	SC <- 0
	temp3 <- c()
	temp4 <- c()
	x <- c()
	for(i in 1:seq.length){
		if(phase.shift.matrix[2,i+1]>SC){
			dots <- character(phase.shift.matrix[2,i+1]-SC)
			dots[] <- "."
			temp3 <- c(temp3,dots)
			SC <- phase.shift.matrix[2,i+1]
		}else if(phase.shift.matrix[2,i+1]<SC){
			dots <- character(SC-phase.shift.matrix[2,i+1])
			dots[] <- "."
			temp4 <- c(temp4,dots)
			SC <- phase.shift.matrix[2,i+1]
		}
		temp3 <- c(temp3,temp1[i])
		temp4 <- c(temp4,temp2[i])
		x <- c(x,abs(phase.shift.matrix[2,i+1]))
	}
	if(SC>0){
		dots <- character(SC)
		dots[] <- "."
		temp4 <-c(temp4,dots)
	}else{
		dots <- character(abs(SC))
		dots[]<- "."
		temp3 <- c(temp3,dots)
	}

	#Outout reconstructed sequences with mismatches highlighted by red color

	j <- length(temp4)
	temp1 <- ""
	temp2 <- ""
	temp5 <- ""
	SCMax <- 0
	SCd <- 0
	mismatches <- 0
	ambiguities <- c("R","Y","K","M","S","W","B","D","H","V","N")
	ambiguities.resolved <- c("AG","CT","GT","AC","CG","AT","CGT","AGT","ACT","ACG","ACGT")
	for(i in 1:j){
		SCa <- temp3[i]
		SCb <- temp4[i]
		SCc <- ""
		if(SCa=="." && SCb=="."){
		}else if(SCa==SCb && (any(SCa==c("A","C","G","T")))){
			temp1 <- temp1 %+% SCa
			temp2 <- temp2 %+% SCb
		}else{
			if(SCa=="."){
				temp1 <- temp1 %+% SCa
			}else{
				temp1 <- temp1 %+% red(SCa)
			}
			if(SCb=="."){
				temp2 <- temp2 %+% SCb
			}else{
				temp2 <- temp2 %+% red(SCb)
			}
		}

		if(SCa == "." &&SCb=="."){
		}else if(SCa=="."){
			temp5 <- temp5 %+% red(SCb)
			if(all(SCb!=c("A","C","G","T"))){
				SCd <- SCd +1
			}
		}else if(SCb=="."){
			temp5 <- temp5 %+% red(SCa)
			if(all(SCa!=c("A","C","G","T"))){
				SCd <- SCd + 1
			}
		}else if(SCa==SCb&&any(SCa==c("A","C","G","T"))){
			temp5 <- temp5 %+% SCa
			SCMax <- SCMax + 1
		}else{
			SCd <- SCd +1 
			if(any(SCa==ambiguities)){
				SCa <- ambiguities.resolved[SCa==ambiguities]
			}
			if(any(SCa==ambiguities)){
				SCa <- ambiguities.resolved[SCa==ambiguities]
			}
			if((countCharOccurrences("A",SCa) >0 && countCharOccurrences("A",SCb) >0) ||
				(countCharOccurrences("C",SCa) >0 && countCharOccurrences("C",SCb) >0) ||
				(countCharOccurrences("G",SCa) >0 && countCharOccurrences("G",SCb) >0) ||
				(countCharOccurrences("T",SCa) >0 && countCharOccurrences("T",SCb) >0)){
				SCMax <- SCMax +1 
			}else{
				mismatches <- mismatches +1
			}
			if(countCharOccurrences("A",SCa) >0||countCharOccurrences("A",SCb) >0){
				SCc <- SCc %+% "A"
			}
			if(countCharOccurrences("C",SCa) >0||countCharOccurrences("C",SCb) >0){
				SCc <- SCc %+% "C"
			}
			if(countCharOccurrences("G",SCa) >0||countCharOccurrences("G",SCb) >0){
				SCc <- SCc %+% "G"
			}
			if(countCharOccurrences("T",SCa) >0||countCharOccurrences("T",SCb) >0){
				SCc <- SCc %+% "T"
			}
			if(countCharOccurrences("G",SCc) >0 && countCharOccurrences("C",SCc) >0 && 
			countCharOccurrences("T",SCc) >0 && countCharOccurrences("A",SCc) >0){
				SCc <- "N"
			}else if(countCharOccurrences("G",SCc) >0&&countCharOccurrences("C",SCc) >0&&
				countCharOccurrences("T",SCc) >0){
				SCc <- "B"
			}else if(countCharOccurrences("G",SCc) >0 && countCharOccurrences("A",SCc) >0 &&
				countCharOccurrences("T",SCc) >0){
				SCc <- "D"
			}else if(countCharOccurrences("A",SCc) >0 && countCharOccurrences("C",SCc) >0 &&
				countCharOccurrences("T",SCc) >0){
				SCc <- "H"
			}else if(countCharOccurrences("G",SCc) >0 && countCharOccurrences("C",SCc) >0 &&
				countCharOccurrences("A",SCc) >0){
				SCc <- "V"
			}else if(countCharOccurrences("A",SCc) >0 && countCharOccurrences("G",SCc) >0){
				SCc <- "R"
			}else if(countCharOccurrences("C",SCc) >0 && countCharOccurrences("T",SCc) >0){
				SCc <- "Y"
			}else if(countCharOccurrences("G",SCc) >0 && countCharOccurrences("T",SCc) >0){
				SCc <- "K"
			}else if(countCharOccurrences("A",SCc) >0 && countCharOccurrences("C",SCc) >0){
				SCc <- "M"
			}else if(countCharOccurrences("C",SCc) >0 && countCharOccurrences("G",SCc) >0){
				SCc <- "S"
			}else if(countCharOccurrences("A",SCc) >0 && countCharOccurrences("T",SCc) >0){
				SCc <- "W"
			}
			temp5 <- temp5 %+% red(SCc)
		}
	}

	return(list(seq1=temp1, seq2=temp2, combined=temp5))


}
