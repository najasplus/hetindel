resolve.all <-
function(phase.shift.matrix, seq.matrix, scores, shift.score, max.shift, is.longindel){

	resolved <- matrix(,nrow=4,ncol=2)

	seq.length <- length(seq.matrix[1,])
		if(shift.score==1){
		repeat{
			phase.shift.matrix <- alignLeft(phase.shift.matrix,seq.length,seq.matrix)$matr
			phase.shift.matrix <- alignRight(phase.shift.matrix,seq.length,seq.matrix)$matr

			#print(phase.shift.matrix)
			#calculate and mark the ambiguities that could potentially be resolved




			for(i in 1:seq.length){
				if(phase.shift.matrix[1, i+1]==0){
					j <- abs(phase.shift.matrix[2,i+1])


					#print(j*shift.score+i+1)
					#print(i)
					if(i<=j){
						phase.shift.matrix[3,i+1] <- 1	#1 could be resolved, 0- cannot be resolved
					}else if(j==0){
						phase.shift.matrix[3,i+1] <- 1
					}else if(phase.shift.matrix[3,i-j+1] == 0 && abs(phase.shift.matrix[2,i-j+1])==j && phase.shift.matrix[1,i-j+1]==0){
						phase.shift.matrix[3,i+1] <- 0
					}else if(phase.shift.matrix[3,i-j+1]==1 && abs(phase.shift.matrix[2,i-j+1])==j && phase.shift.matrix[1,i-j+1]==0){
						phase.shift.matrix[3,i+1] <- 1
					}else if(j*shift.score+i+1 <=seq.length){
						shift.score <- 1
						while(phase.shift.matrix[1,j*shift.score+i+1]==0 && abs(phase.shift.matrix[2,j*shift.score+1+i])==j && 
							j*shift.score+i+1<=seq.length){
							shift.score <- shift.score+1
						}
						#cat("erste Stelle",3-phase.shift.matrix[1,j*shift.score+i+1],"\n")
						#cat("zweite Stelle",j*shift.score+i,"\n")
						#cat("Wert",seq.matrix[3-phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i],"\n")


						if(abs(phase.shift.matrix[2,i-j+1]) < j && phase.shift.matrix[2,i+1] >0){
							phase.shift.matrix[3,i+1] <- 1
						}else if((j*shift.score +i) > seq.length){ # evtl +1
							phase.shift.matrix[3, i+1] <- 1
						}else if(abs(phase.shift.matrix[2,j*shift.score+i+1]) != j){
							phase.shift.matrix[3,i+1] <- 1
						}else if(abs(phase.shift.matrix[2,j*shift.score+i+1+phase.shift.matrix[4,i+1]]) != j){
							phase.shift.matrix[3,i+1] <- 1


						}else if(phase.shift.matrix[1,i-j+1]>0 && (3-phase.shift.matrix[1,j*shift.score+i+1]<3 )&& 
							(3-phase.shift.matrix[1,j*shift.score+i+1])>0){

							if((phase.shift.matrix[2,i+1] > 0) && ((shift.score / 2 ) != floor(shift.score/2)) && 
							((seq.matrix[phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] && 
									seq.matrix[1,i] == seq.matrix[3-phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]) || 
								(seq.matrix[phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i] && 
									seq.matrix[2,i]==seq.matrix[3-phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]))){
								

								phase.shift.matrix[3,i+1] <- 1
							}
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*shift.score+i+1])<3 && 
							3-phase.shift.matrix[1,j*shift.score+i+1]>0){

							if(phase.shift.matrix[2,i+1] >0 && (shift.score/2) == floor(shift.score/2) && 
								((seq.matrix[phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i] && 
									seq.matrix[2,i] == seq.matrix[3-phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i])||
								(seq.matrix[phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[1,i] && 
									seq.matrix[1,i] == seq.matrix[3- phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]))){
								phase.shift.matrix[3,i+1] <- 1
							}
						}else if(abs(phase.shift.matrix[2,i-j+1]) != j && phase.shift.matrix[2,i+1] < 0){
							phase.shift.matrix[3,i+1] <- 1
						}else if(phase.shift.matrix[1,i-j+1] ==0){
							phase.shift.matrix[3,i+1] <- 0
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*shift.score+i+1])<3 && 
							3-phase.shift.matrix[1,j*shift.score+i+1]>0){

							if(phase.shift.matrix[2,i+1] < 0 && (shift.score/2) != floor(shift.score/2) && 
								((seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[1,i] && 
									seq.matrix[2,i]== seq.matrix[phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]) || 
								(seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] && 
									seq.matrix[1,i]==seq.matrix[phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]))){
								phase.shift.matrix[3,i+1] <- 1

						}
						}else if(phase.shift.matrix[1,i-j+1]!=0 && (3-phase.shift.matrix[1,j*shift.score+i+1])<3 && 
							3-phase.shift.matrix[1,j*shift.score+i+1]>0){

							if(phase.shift.matrix[2,i+1] < 0 && (shift.score/2) == floor(shift.score/2) && 
								(( seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j]==seq.matrix[1,i-j] && 
									seq.matrix[1,i] == seq.matrix[phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]) || 
								(seq.matrix[3-phase.shift.matrix[1,i-j+1],i-j] == seq.matrix[2,i] &&
									seq.matrix[2,i]==seq.matrix[phase.shift.matrix[1,j*shift.score+i+1],j*shift.score+i]))){
								phase.shift.matrix[3,i+1] <- 1
							}
						}else{

							phase.shift.matrix[3,i+1] <- 0
						}
					}
				}
			}

			j <- 0

			for(i in 1:seq.length){
				

				if(phase.shift.matrix[1,i+1]==0 && phase.shift.matrix[3,i+1]==1){

					for(z in 1:(max.shift+phase.shift.matrix[4,1])){
						if(z>=i || z> seq.length -i){
							break
						}
						if(phase.shift.matrix[2,i-z+1] != phase.shift.matrix[2,i+z+1]){
							break
						}
					}

					if(z<=i && z-phase.shift.matrix[3,i+1] <= max.shift+1 && z <= seq.length -(i+1)){

						ind1 <- abs(phase.shift.matrix[2,i-z+1])
						ind2 <- abs(phase.shift.matrix[2,i+z+1])
						if(!is.longindel){
							if((i> ind1) && (i>ind2) && (seq.length -i)>ind1 && 
								(seq.length-i)>ind2){
								if(((seq.matrix[2,i]==seq.matrix[1,i-ind1] && 
									scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2]) || 
								(seq.matrix[2,i] == seq.matrix[2,i-ind1] && 
									scores[i-ind1+1,ind1+1,1] <= scores[i-ind1+1,ind1+1,2])) && 
								((seq.matrix[1,i] == seq.matrix[2,i+ind1] && 
									scores[i+ ind1+1,ind1+1,1] >= scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[1,i] == seq.matrix[1,i+ind1] && 
									scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
									resolved[1,1] <- 1
								}else{
									resolved[1,1] <- 0
								}
								if(((seq.matrix[1,i]==seq.matrix[1,i-ind1] && 
									scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
								(seq.matrix[1,i] == seq.matrix[2,i-1] &&
									scores[i-ind1+1,ind1+1,1] <= scores[i-ind1+1,ind1+1,2]))&&
								((seq.matrix[2,i]==seq.matrix[2,i+ind1] && 
									scores[i+ind1+1,ind1+1,1]>=scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[1,i+ind1] && 
									scores[i+ind1+1,ind1+1,1]<=scores[i+ind1+1,ind1+1,2]))){
									resolved[1,2] <- 1
								}else{
									resolved[1,2] <- 0
								}


								if(((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]>=scores[i-ind2+1,ind2+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]<= scores[i-ind2+1,ind2+1,2]))&&
								((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
								(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[2,1] <- 1
								}else{
									resolved[2,1] <- 0
								}


								if(((seq.matrix[1,i]==seq.matrix[1,i-ind2] && 
									scores[i-ind2+1,ind2+1,1]>=scores[i-ind2+1,ind2+1,2])||
								(seq.matrix[1,i]==seq.matrix[2,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]<=scores[i-ind2+1,ind2+1,2]))&&
								((seq.matrix[2,i]==seq.matrix[2,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
								(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[2,2] <- 1
								}else{
									resolved[2,2] <- 0
								}
								if((ind1 < ind2) &&
									((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&
										scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
									(seq.matrix[2,i]==seq.matrix[2,i-ind1]&&
										scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2]))&&
									((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&
										scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
									(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&
										scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[4,1] <- 1
								}else{
									resolved[4,1] <- 0
								}
								if((ind1<ind2)&&
									((seq.matrix[1,i]==seq.matrix[1,i-ind1]&&
										scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
									(seq.matrix[1,i]==seq.matrix[2,i-ind1]&&
										scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2]))&&
									((seq.matrix[2,i]==seq.matrix[2,i+ind2]&&
										scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
									(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&
										scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[4,2] <- 1
								}else{
									resolved[4,1] <- 0
								}
							}else if(i<=ind1){
								if(((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&
									scores[i+ind1+1,ind1+1,1]>=scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&
									scores[i+ind1+1,ind1+1,1]<=scores[i+ind1+1,ind1+1,2]))){
									resolved[1,1] <- 1
								}else{
									resolved[1,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&
									scores[i+ind1+1,ind1+1,1]>=scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&
									scores[i+ind1+1,ind1+1,1]<=scores[i+ind1+1,ind1+1,2]))){
									resolved[1,2] <- 1
								}else{
									resolved[1,2] <- 0
								}
								resolved[2,] <- 0
							}else if(seq.length-i< ind2){
								resolved[1,] <- 0
								if(((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]>=scores[i-ind2+1,ind2+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]<=scores[i-ind2+1,ind2+1,2]))){
									resolved[2,1] <- 1
								}else{
									resolved[2,1] <- 0
								}
								if(((seq.matrix[1,i]==seq.matrix[1,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]>=scores[i-ind2+1,ind2+1,2])||
								(seq.matrix[1,i]==seq.matrix[2,i-ind2]&&
									scores[i-ind2+1,ind2+1,1]<=scores[i-ind2+1,ind2+1,2]))){
									resolved[2,2] <- 1
								}else{
									resolved[2,2] <- 0
								}
							}							
							if(phase.shift.matrix[2,i-z+1]<phase.shift.matrix[2,i+z+1]||i<=ind1){
								
								if(((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
								(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[3,1] <- 1
								}else{
									resolved[3,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]>=scores[i+ind2+1,ind2+1,2])||
								(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&
									scores[i+ind2+1,ind2+1,1]<=scores[i+ind2+1,ind2+1,2]))){
									resolved[3,2] <- 1
								}else{
									resolved[3,2] <- 0
								}
							}else{

								if(ind1==0||((seq.matrix[2,i]==seq.matrix[1,i-ind1]&&
									scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind1]&&
									scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2]))){
									resolved[3,1] <- 1
								}else{
									resolved[3,1] <- 0
								}
								if(ind1==0||((seq.matrix[1,i]==seq.matrix[1,i-ind1]&&
									scores[i-ind1+1,ind1+1,1]>=scores[i-ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind1]&&
									scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2]))){
									resolved[3,2] <- 1
								}else{
									resolved[3,2] <- 0
								}
							}
							if(ind1==0||ind2==0){
								resolved[4,] <- 0
							}	

						}else{  #long indel
							if(phase.shift.matrix[2,i-z+1]>0){
								if(((seq.matrix[2,i]==seq.matrix[1,i-ind1]&&
									phase.shift.matrix[1,i-ind1+1] !=2)||
								(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind1+1]!=1))){
									resolved[3,1] <- 1
								}else{
									resolved[3,1] <- 0
								}
								if(((seq.matrix[1,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[1,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1]!=1))){
									resolved[3,2] <- 1
								}else{
									resolved[3,2] <- 0
								}
							}else{
								if(((seq.matrix[1,i]==seq.matrix[2,i-ind1] && phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[1,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))){
									resolved[4,1] <- 1
								}else{
									resolved[4,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[2,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))){
									resolved[4,2] <- 1
								}else{
									resolved[4,2] <- 0
								}
							}
							if(phase.shift.matrix[2,i+z+1]>0){
								if(((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
									(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+2]!=1))){
									resolved[4,1] <- 1
								}else{
									resolved[4,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
									(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
									resolved[4,2] <- 1
								}else{
									resolved[4,2] <- 0
								}
							}else{
								if(((seq.matrix[2,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
									(seq.matrix[2,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
									resolved[4,1] <- 1
								}else{
									resolved[4,1] <- 0
								}
								if(((seq.matrix[1,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
									(seq.matrix[1,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
									resolved[4,2] <- 1
								}else{
									resolved[4,2] <- 0
								}
							}
							if(i >= ind1 && i>= ind2 && (seq.length-i)>ind1&&
								(seq.length-i)<ind2){
								if(phase.shift.matrix[2,i-z+1]>0){
									if(((seq.matrix[2,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))&&
									((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[1,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1]!=1))&&
									((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}else{
									if(((seq.matrix[1,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))&&
									((seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i]==seq.matrix[2,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i-ind1]&&phase.shift.matrix[1,i-ind1+1]!=1))&&
									((seq.matrix[1,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}
								if(phase.shift.matrix[2,i+z+1]>0){
									if(((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=1))&&
									((seq.matrix[1,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
										resolved[2,1] <- 1
									}else{
										resolved[2,1] <- 0
									}
									if(((seq.matrix[1,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=1))&&
									((seq.matrix[2,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
										resolved[2,2] <- 1
									}else{
										resolved[2,2] <- 0
									}
								}else{
									if(((seq.matrix[1,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=1))&&
									((seq.matrix[2,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2]!=2)||
										(seq.matrix[2,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
										resolved[2,1] <- 1
									}else{
										resolved[2,1] <- 0
									}
									if(((seq.matrix[2,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind2]!=1))&&
									((seq.matrix[1,i]==seq.matrix[1,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i+ind2]&&phase.shift.matrix[1,i+ind2+1]!=1))){
										resolved[2,2] <- 1
									}else{
										resolved[2,2] <- 0
									}
								}
							}else if(i<=ind1){
									if(((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&phase.shift.matrix[1,i+ind1]!=2)||
										(seq.matrix[2,i]==seq.matrix[1,i+ind1]&&phase.shift.matrix[1,i+ind1+1]!=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
									resolved[2,] <- 0
							}else if((seq.length-i)<ind2){
									resolved[1,] <- 0
									if(((seq.matrix[2,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind1]!=2)||
										(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind1]!=1))){
										resolved[2,1] <- 1
									}else{
										resolved[2,1] <- 0
									}
									if(((seq.matrix[1,i]==seq.matrix[1,i-ind2]&&phase.shift.matrix[1,i-ind1]!=2)||
										(seq.matrix[1,i]==seq.matrix[2,i-ind2]&&phase.shift.matrix[1,i-ind1]!=1))){
										resolved[2,2] <- 1
									}else{
										resolved[2,2] <- 0
									}
							}
						}

						if(ind1==0){
							resolved[1,] <- 0
						}
						if(ind2==0){
							resolved[2,] <- 0
						}
						if(!is.longindel){
							if(ind1>ind2){
								if(phase.shift.matrix[2,i]==ind1){
									if(z<=ind1){
										resolved[1,] <- 0
										if(z<=ind2){
											resolved[3,] <- 0
										}
									}

									if((phase.shift.matrix[4,i+1]-z)<0){
										resolved[2,] <- 0
									}
								}else{
									resolved[1,] <- 0
									resolved[3,] <- 0
								}
							}else{ #ind1<ind2
								if(phase.shift.matrix[2,i+1]==ind1){
									if(phase.shift.matrix[3,i+1]<=ind2){
										resolved[2,] <- 0
										if(ind1 != 0 && z > phase.shift.matrix[4,i+1]){
											resolved[3,] <- 0
										}
									}
									if(z>ind2){
										resolved[2,] <- 0
									}
								}else{
									if(z+phase.shift.matrix[4,i+1]<=ind2){
										resolved[2,] <- 0
									}
									if(z>ind1){
										resolved[4,] <- 0
									}
									resolved[1,] <- 0
								}
							}
						}else{
							if(ind1>ind2){
								if(abs(phase.shift.matrix[2,i+1])==ind1){
									if(z<=ind1){
										resolved[1,] <- 0
									}
									if((phase.shift.matrix[4,i+1]-z)<0){
										resolved[c(2,4),] <- 0
									}
								}else{
									resolved[c(1,3),] <- 0
								}
							}else{
								if(phase.shift.matrix[2,i+1] <= ind2){
									if(phase.shift.matrix[4,i+1]<=ind2){
										resolved[2,] <- 0
									}
									if(ind1!=0 &&z > phase.shift.matrix[4,i+1]){
										resolved[4,] <- 0
									}
									if(z>ind2){
										resolved[2,] <- 0
									}
								}else{
									if((z+phase.shift.matrix[4,i+1])<=ind2){
										resolved[2,] <- 0
									}
									resolved[c(1,3),] <- 0
								}
							}
						}
						if(!is.longindel){
							if(all(resolved[1,]==1)){
							}else if(all(resolved[2,]==1)){
							}else if(all(resolved[c(1,2),1]==1)&&all(resolved[c(1,2),2]==0)){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[c(1,2),1]==0)&&all(resolved[c(1,2),2]==1)){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}else if((resolved[1,1]==1&&resolved[2,2]==1)||(resolved[1,2]==1&&resolved[2,1]==1)){
							}else if(all(resolved[c(1,3),1]==1)&&all(resolved[c(1,3),2]==0)){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[c(1,3),1]==0)&&all(resolved[c(1,3),2]==1)){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}else if(((resolved[1,1]==1&&resolved[3,2]==1)||(resolved[1,2]==1&&resolved[3,1]==1))
								&&z< ind1 &&phase.shift.matrix[4,i+1]>ind1){
							}else if(((resolved[2,1]==1&&resolved[3,2]==1)||(resolved[2,2]==1&&resolved[3,1]==1))
								&&z<=ind1 &&phase.shift.matrix[4,i+1]>ind2 &&ind2 >ind1){
							}else if(((resolved[2,1]==1&&resolved[3,2]==1)||(resolved[2,2]==1&&resolved[3,1]==1))
								&&phase.shift.matrix[2,i+1]==ind2&&phase.shift.matrix[1,i+ind2+1]==0){
							}else if(((resolved[2,1]==1&&resolved[3,2]==1)||(resolved[2,2]==1&&resolved[3,1]==1))
								&&phase.shift.matrix[2,i+1]==ind1&&phase.shift.matrix[4,i+1]>=z){
							}else if(((resolved[1,1]==1&&resolved[3,2]==1)||(resolved[1,2]==1&&resolved[3,1]==1))
								&&phase.shift.matrix[2,i+1]==ind1&&ind2>ind1&&phase.shift.matrix[4,i+1]>=(z+ind1)){
							}else if(xor(resolved[1,1]==1,resolved[1,2]==1)){
								if(resolved[1,1]==1){
									phase.shift.matrix[1,i+1] <- 1
									scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
								}else{
									phase.shift.matrix[1,i+1] <- 2
									scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
								}
							}else if(xor(resolved[2,1]==1,resolved[2,2]==1)){
								if(resolved[2,1]==1){
									phase.shift.matrix[1,i+1] <- 1
									scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
								}else{
									phase.shift.matrix[1,i+1] <- 2
									scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
								}
							}else if(all(resolved[3,]==1)&&all(resolved[4,]==0)){
							}else if(all(resolved[3,]==1)&&all(resolved[4,]==c(1,0))){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[3,]==1)&&all(resolved[4,]==c(0,1))){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}else if(resolved[3,1]==1){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(resolved[3,2]==1){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}
						}else{

							if(all(resolved[1,]==1)){
							}else if(all(resolved[2,]==1)){
							}else if(((resolved[1,1]==1&&resolved[3,2]==1)||(resolved[1,2]==1&&resolved[3,1]==1))){
							}else if(((resolved[2,1]==1&&resolved[4,2]==1)||(resolved[2,2]==1&&resolved[4,1]==1))){
							}else if(((resolved[1,1]==1&&resolved[4,2]==1)||(resolved[1,2]==1&&resolved[4,1]==1))){
							}else if(((resolved[2,1]==1&&resolved[3,2]==1)||(resolved[2,2]==1&&resolved[3,1]==1))){
							}else if(((resolved[4,1]==1&&resolved[3,2]==1)||(resolved[4,2]==1&&resolved[3,1]==1))){
							}else if(xor(resolved[1,1]==1,resolved[1,2]==1)){
								if(resolved[1,1]==1){
									phase.shift.matrix[1,i+1] <- 1
									scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
								}else{
									phase.shift.matrix[1,i+1] <- 2
									scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
								}
							}else if(xor(resolved[2,1]==1,resolved[2,2]==1)){
								if(resolved[2,1]==1){
									phase.shift.matrix[1,i+1] <- 1
									scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
								}else{
									phase.shift.matrix[1,i+1] <- 2
									scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
									scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
								}
							}else if(all(resolved[3,]==1)){
							}else if(all(resolved[4,]==1)){
							}else if(all(resolved[4,]==c(1,0))){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[4,]==c(0,1))){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[3,]==c(1,0))){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind2+1,1] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(all(resolved[3,]==c(0,1))){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind2+1,2] <- max(scores[i+1,ind2+1,])+1 
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}
						}
					}else{
						if(i==1||i==seq.length||phase.shift.matrix[2,i]==phase.shift.matrix[2,i+2]){
							if(i!=seq.length){
								ind1 <- abs(phase.shift.matrix[2,i+2])
								ind2 <- phase.shift.matrix[2,i+2]
							}else{
								ind1 <- abs(phase.shift.matrix[2,i])
								ind2 <- phase.shift.matrix[2,i]
							}
							if(i<=ind1){
								if(((seq.matrix[1,i]==seq.matrix[2,i+ind1]&&
									scores[i+ind1+1,ind1+1,1] >= scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&
									scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
									resolved[1,1] <- 1
								}else{
									resolved[1,1] <- 0
								}
								if(((seq.matrix[2,i]==seq.matrix[2,i+ind1]&&
									scores[i+ind1+1,ind1+1,1] >= scores[i+ind1+1,ind1+1,2])||
								(seq.matrix[2,i]==seq.matrix[2,i-ind2]&&
									scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
									resolved[1,2] <- 1
								}else{
									resolved[1,2] <- 0
								}
							}else if(seq.length - i < ind1){
								if(ind2 > 0){
									if(((seq.matrix[2,i] == seq.matrix[1,i - ind1] &&
										scores[i - ind1 + 1,ind1+1,1] >= scores[i - ind1 + 1,ind1+1,2])||
									(seq.matrix[2,i] == seq.matrix[2,i - ind1] &&
										scores[i - ind1 + 1,ind1+1,1] <= scores[i - ind1 + 1,ind1+1,2]))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[1,i] == seq.matrix[1,i - ind1] &&
										scores[i - ind1 + 1,ind1+1,1] >= scores[i - ind1 + 1,ind1+1,2]) ||
									(seq.matrix[1,i] == seq.matrix[2,i - ind1] &&
										scores[i - ind1 + 1,ind1+1,1] <= scores[i - ind1 + 1,ind1+1,2]))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}else{
									if(((seq.matrix[1,i] == seq.matrix[2,i-ind1] &&
										phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[1,i] == seq.matrix[1,i-ind1] &&
										phase.shift.matrix[1,i-ind1+1]!=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i] == seq.matrix[2,i-ind1] &&
										phase.shift.matrix[1,i-ind1+1]!=2)||
									(seq.matrix[2,i] == seq.matrix[1,i-ind1] &&
										phase.shift.matrix[1,i-ind1+1]!=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}
							}else{
								if(ind2 > 0){
									if(((seq.matrix[2,i] == seq.matrix[1,i-ind1] && 
										scores[i-ind1+1,ind1+1,1] >= scores[i-ind1+1,ind1+1,2])||
									(seq.matrix[2,i] == seq.matrix[2,i-ind1]&&
										scores[i-ind1+1,ind1+1,1]<= scores[i-ind1+1,ind1+1,2])) &&
									((seq.matrix[1,i] == seq.matrix[2,i+ind1] &&
										scores[i+ind1+1,ind1+1,1]>=scores[i+ind1+1,ind1+1,2]) ||
									(seq.matrix[1,i]==seq.matrix[1,i+ind1]&&
										scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[1,i] == seq.matrix[1,i-ind1] && 
										scores[i-ind1+1,ind1+1,1] >= scores[i-ind1+1,ind1+1,2])||
									(seq.matrix[1,i]==seq.matrix[2,i-ind1] && 
										scores[i-ind1+1,ind1+1,1]<=scores[i-ind1+1,ind1+1,2])) &&
									((seq.matrix[2,i]==seq.matrix[2,i+ind1] &&
										scores[i+ind1+1,ind1+1,1] >= scores[i+ind1+1,ind1+1,2])||
									(seq.matrix[2,i]==seq.matrix[1,i+ind1] &&
										scores[i+ind1+1,ind1+1,1] <= scores[i+ind1+1,ind1+1,2]))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}else{
									if(((seq.matrix[1,i] == seq.matrix[2,i-ind1] &&phase.shift.matrix[1,i-ind1+1] !=2) || 
									(seq.matrix[1,i] == seq.matrix[1,i-ind1] && phase.shift.matrix[1,i-ind1+1] !=1)) &&
									((seq.matrix[2,i] == seq.matrix[1,i+ind1] && phase.shift.matrix[1,i+ind1+1] != 2) || 
									(seq.matrix[2,i] == seq.matrix[2,i+ind1] && phase.shift.matrix[1,i+ind1+1] !=1))){
										resolved[1,1] <- 1
									}else{
										resolved[1,1] <- 0
									}
									if(((seq.matrix[2,i] == seq.matrix[2,i-ind1] && phase.shift.matrix[1,i-ind1+1] !=2) ||
										(seq.matrix[2,i] == seq.matrix[1,i-ind1] && phase.shift.matrix[1,i-ind1+1]!=1)) &&
									((seq.matrix[1,i] == seq.matrix[1,i+ind1] && phase.shift.matrix[1,i+ind1+1]!=2)||
										(seq.matrix[1,i] == seq.matrix[2,i+ind1] && phase.shift.matrix[1,i+ind1+1] !=1))){
										resolved[1,2] <- 1
									}else{
										resolved[1,2] <- 0
									}
								}
							}
							if(resolved[1,1] ==1 &&resolved[1,2] ==0){
								phase.shift.matrix[1,i+1] <- 1
								scores[i+1,ind1+1,1] <- max(scores[i+1,ind1+1,])+1
							}else if(resolved[1,1] == 0 && resolved[1,2] ==1){
								phase.shift.matrix[1,i+1] <- 2
								scores[i+1,ind1+1,2] <- max(scores[i+1,ind1+1,])+1
							}
						}
					}
					if(phase.shift.matrix[1,i+1] != 0){
						j <- 1
					}
				}
			}


			if(j!=1){
				break
			}
		}

	}
	return(list(seq.matrix=seq.matrix, phase.shift.matrix=phase.shift.matrix, resolved=resolved, scores=scores))

}
