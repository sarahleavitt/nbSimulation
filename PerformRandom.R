#Sarah V. Leavitt
#Boston University Dissertation
#Simulation

################################################################################
# This program calculates probability of transmission only using
# generation interval
################################################################################


#### Inputs ####

#df = pair-level dataset ordered by observationDate
#indIDVar = name (in quotes) of the column with the individual ID.
#(data frame orderedPair must have columns called <indIDVar>.1 and <indIDVar>.2).
#pairIDVar = 	name (in quotes) of the column with the unique pair ID variable. 


#### Outputs ####

#One data frame with the predicted probabilities for each pair


########### Interval Method ##########

performRandom <- function(df, indIDVar = "individualID", pairIDVar = "edgeID"){
  
  indIDVar2 <- paste0(indIDVar, ".2")

  probs <- as.data.frame(df)
  probs$pAvg = runif(nrow(probs), 0, 1)
  
  #Calculating scaled probabilities
  #Calculating the total of all probabilities per infectee
  totalP <- stats::aggregate(probs$pAvg, by = list(probs[, indIDVar2]), sum, na.rm = TRUE)
  names(totalP) <- c(indIDVar2, "pTotal")
  
  #Calculating the scaled probabilities
  probs2 <- merge(probs, totalP, by = indIDVar2)
  probs2$pScaled <- ifelse(probs2$pTotal != 0, probs2$pAvg / probs2$pTotal, 0)
  #Ranking the probabilities for each possible infector
  #Ties are set to the minimum rank of that group
  probs2 <- probs2[order(probs2[, indIDVar2], -probs2$pScaled), ]
  probs2$pRank <- stats::ave(-probs2$pScaled, probs2[, indIDVar2], 
                             FUN = function(x){
                               rank(x, ties.method = "min") 
                             })
  
  #Only keeping columns of interest
  probs2$label <- "Random"
  probs2 <- probs2[, c("label", pairIDVar, "pAvg", "pScaled", "pRank")]
  
  return(probs2)
}

