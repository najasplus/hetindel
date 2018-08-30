
countCharOccurrences <- function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
}

compare <- function(target, lead){
	threshold <- 5
	errors <- 0
	TP <- 0
	FP <- 0
	TN <- 0
	FN <- 0
	indel.length <- countCharOccurrences("\\.",target)
	found.indel.length <- countCharOccurrences("\\.",lead)
	diff <- abs(nchar(target) - nchar(lead))
	target <- trimws(target)
	target <- paste(target, substr(target,nchar(target)-diff+1,nchar(target)),sep="")
	for(i in 1:nchar(target)){
		if(substr(target,i,i)!=substr(lead,i,i)){
			errors <- errors + 1
			if(substr(target,i,i)=="."){
				FN <- FN +1
			}else if(substr(lead,i,i)=="."){
				FP <- FP +1
			}else{
				TN <- TN +1 
			}
		}else{
			if(substr(target,i,i)=="."){
				TP <- TP +1
			}else{
				TN <- TN +1
			}
		}
	}
	indel.found <- ifelse(indel.length == found.indel.length && errors >= threshold, 1, 0)
	errorless <- ifelse(errors < threshold, 1,0)
	return(list(errors=errors,basenumber=nchar(target),indel.found=indel.found,errorless=errorless,TP=TP,TN=TN,FP=FP,FN=FN))
}

len <- c("150","300","700")
pos <- c("middle", "first_quater","last_quater")
tests <- list(short.indel.error=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			long.indel.error=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			short.indel.perc=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			long.indel.perc=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			short.indel.found=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			long.indel.found=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			short.indel.complete=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			long.indel.complete=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			short.indel.false=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos)),
			long.indel.false=
			matrix(0, nrow=3,ncol=3,dimnames=list( len,pos))
			)
runs <- 1000
lengths <- c("150", "300", "700")
positions <- c("m", "q1", "q3")
indel.lengths <- c("9","50")
lengths.text <- c("short indels","long indels")
types <- c("per base error","percentage per base error", "found indels", "reconstructed sequence", "incorrect sequence")

reconstructed.seq <- 0
correct.indels <- 0
incorrect.indels <- 0

TP <- 0 
TN <- 0
FP <- 0
FN <- 0
for(indel.length in 1:length(indel.lengths)){
	for(position in 1:length(positions)){
		for(length in 1:length(lengths) ){
			error <- 0
			found <- 0
			no.err <- 0
			total.bases <- 0
			for(i in 1:runs){
				cat(paste0("Run ",i," of ",runs," for ",lengths.text[indel.length],", ",positions[position]," position, and ",lengths[length]," length\n"))
				output <- system(paste("Rscript indelScript.R -m ",indel.lengths[indel.length]," -G ",lengths[length],",TRUE,",positions[position],sep=""), intern=TRUE)
				generated <- output[1]
				lead <- output[2]
				result <- compare(generated, lead)

				error <- error + result$errors
				total.bases <- total.bases + result$basenumber
				found <- found + result$indel.found
				no.err <- no.err + result$errorless
				TP <- TP + result$TP
				FP <- FP + result$FP
				TN <- TN + result$TN
				FN <- FN + result$FN

			}

			tests[[indel.length]][length,position] <- error
			tests[[indel.length+2]][length,position] <- round(error / total.bases *100 ,4)
			tests[[indel.length+4]][length,position] <- found 
			tests[[indel.length+6]][length,position] <- no.err
			tests[[indel.length+8]][length,position] <- runs - (found + no.err)

			if(result$errorless){
				reconstructed.seq <- reconstructed.seq + 1
			}else{
				if(result$indel.found){
					correct.indels <- correct.indels + 1
				}else{
					incorrect.indels <- incorrect.indels + 1
				}
			}
		}
	}
}

sensitivity <- TP / (TP+FN)
specificity <- TN / (TN+FP)
accuracy <- (TP + TN) / (TP + TN + FP + FN)
cat("Results for ")
cat(runs)
cat(" runs of indel per option:\n")

#cat(tests)
for(j in 1:length(types)){
	for(i in 1:length(indel.lengths)){
		cat(paste0("For ",lengths.text[i]," and ",types[j],"\n"))
		cat(paste0("\t\tPositions\n"))
		cat(paste0("Lengths\t"))
		for(p in 1:length(positions)){
			cat(paste0(pos[p],"\t"))
		}
		cat("\n")
		for(l in 1:length(lengths)){
			cat(paste0(lengths[l],"\t"))
			for(p in 1:length(positions) ){
				cat(paste0(tests[[i+(2*(j-1))]][l,p],"\t\t"))
			}
			cat("\n")
		}
	}
	cat("\n")
}
cat("\n")
cat("Number of perfectly reconstructed sequences: ")
cat(reconstructed.seq)
cat("\n")
cat("number of correctly found indel lengths: ")
cat(correct.indels)
cat("\n")
cat("number of failed reconstructions: ")
cat(incorrect.indels)
cat("\n")
cat(paste0("Sensitivity for per Base analysis was: ", as.character(sensitivity),"\n"))
cat(paste0("Specificity for per Base analysis was: ", as.character(specificity),"\n"))
cat(paste0("Accuracy for per Base analysis was: ", as.character(accuracy),"\n"))
