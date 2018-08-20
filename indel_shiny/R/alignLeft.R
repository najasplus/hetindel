alignLeft <-
function(matrix,seq.length,seq.matrix){
	z1 <- 0
	for(i in seq.length:2){
		matrix[4,i] <- 0
		if((seq.length - i-1 - matrix[2,i]) >0){
			j <- matrix[2,i+1]
			x <- matrix[2,i]
			if(x<j){
				matrix[4,i] <- 1
				if(any(matrix[1,i]==c(0,1))&&any(matrix[1,i+j]==c(0,1))&&
					seq.matrix[1,i-1]==seq.matrix[2,i+j-1]){
					##########################################
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(2,0))&&
					seq.matrix[2,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(1,0))&&
					seq.matrix[2,i-1]==seq.matrix[2,i+j-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(2,0))&&
					seq.matrix[1,i-1]==seq.matrix[1,i+j-1]){
					matrix[2,i] <- j
				}else{
					matrix[4,i] <- 0
				}
			}else if(x>0&&j<0){
				matrix[4,i] <- 1
				if(i<j){
					matrix[4,i] <- 0
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i-j]==c(1,0))&&
					seq.matrix[1,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i-j]==c(2,0))&&
					seq.matrix[2,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i-j]==c(1,0))&&
					seq.matrix[1,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i-j]==c(2,0))&&
					seq.matrix[2,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <-j
				}else{
					matrix[4,i] <- 0
				}
			}else if( x>j){
				matrix[4,i] <- 1
				if(i<j){
					matrix[4,i] <- 0
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(1,0))&&
					any(matrix[1,i-j]==c(1,0))&&seq.matrix[1,i-1]==seq.matrix[2,i+j-1]&&
					seq.matrix[1,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(1,0))&&
					any(matrix[1,i-j]==c(2,0))&&seq.matrix[1,i-1]==seq.matrix[2,i+j-1]&&
					seq.matrix[2,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(2,0))&&
					any(matrix[1,i-j]==c(1,0))&&seq.matrix[2,i-1]==seq.matrix[1,i+j-1]&&
					seq.matrix[1,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(2,0))&&
					any(matrix[1,i-j]==c(2,0))&&seq.matrix[2,i-1]==seq.matrix[1,i+j-1]&&
					seq.matrix[2,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(1,0))&&
					any(matrix[1,i-j]==c(1,0))&&seq.matrix[2,i-1]==seq.matrix[2,i+j-1]&&
					seq.matrix[1,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(2,0))&&any(matrix[1,i+j]==c(1,0))&&
					any(matrix[1,i-j]==c(2,0))&&seq.matrix[2,i-1]==seq.matrix[2,i+j-1]&&
					seq.matrix[2,i-j-1]==seq.matrix[1,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(2,0))&&
					any(matrix[1,i-j]==c(1,0))&&seq.matrix[1,i-1]==seq.matrix[1,i+j-1]&&
					seq.matrix[1,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else if(any(matrix[1,i]==c(1,0))&&any(matrix[1,i+j]==c(2,0))&&
					any(matrix[1,i-j]==c(2,0))&&seq.matrix[1,i-1]==seq.matrix[1,i+j-1]&&
					seq.matrix[2,i-j-1]==seq.matrix[2,i-1]){
					matrix[2,i] <- j
				}else{
					matrix[4,i] <- 0
				}
			}
		}
		if(matrix[4,i] > 0){
			z1 <- z1+1
		}
	}
	return(list(matr=matrix,val=z1))
}
