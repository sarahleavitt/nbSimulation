#Sarah V. Leavitt
#Boston University Dissertation
#Simulation


################################################################################
# This program calculates probability of transmission only using
# generation interval
################################################################################


#### Inputs ####

#df = pair-level dataset ordered by observationDate
#genShape = shape parameter for the gamma distribution describing the generation time in years
#genScale = scale parameter for the gamma distribution describing the generation time in years
#shift = number of years to shift the gamma distribution describing the generation time
#observationDiff = variable name that denotes the observation time difference between cases in days
#indIDVar = name (in quotes) of the column with the individual ID.
#(data frame orderedPair must have columns called <indIDVar>.1 and <indIDVar>.2).
#pairIDVar = 	name (in quotes) of the column with the unique pair ID variable. 


#### Outputs ####

#One data frame with the predicted probabilities for each pair with:
  #label = denotes that this is "Random" probabilities 
  #pAvg = the random probability (called pAvg to match with the nbProbabilities function)
  #pScaled = the probability scaled to add to one for all possible infectors for
    #each infectee
  #pRank = the rank of the probability of that infector for the infectee.


performInterval <- function(df, genShape, genScale, shift, observationDiff,
                            indIDVar = "individualID", pairIDVar = "edgeID"){
  
  indIDVar2 <- paste0(indIDVar, ".2")
  
  probs <- as.data.frame(df)
  probs$timeDiff <- round(probs[, observationDiff] / 365, 3)
  probs$pAvg <- ifelse(probs$timeDiff < shift, 0,
                       pgamma(probs$timeDiff - shift,
                              shape = genShape, scale = genScale) -
                       pgamma(probs$timeDiff - (shift + 0.001),
                              shape = genShape, scale = genScale))
  
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
  probs2$label <- paste0("GenInt", round(genShape, 1), "/", round(genScale, 1))
  probs2 <- probs2[, c("label", pairIDVar, "pAvg", "pScaled", "pRank")]
  
  return(probs2)
}

