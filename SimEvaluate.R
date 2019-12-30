#Sarah V. Leavitt
#Boston University Dissertation
#Simulation

################################################################################
# This program contains a function that calculates various metrics to evaluate
# the performance of the nb transmission method
################################################################################


simEvaluate <- function(probs, truthVar = "transmission", pVar = "pScaled"){

  probs <- as.data.frame(probs)
  probs$p <- probs[, pVar]
  probs$truth <- probs[, truthVar]
  
  #If outbreakID is missing, setting it to 1 (only need this variable for simulations)
  if(!"outbreakID.2" %in% names(probs)){
    probs$outbreakID.1 = 1
    probs$outbreakID.2 = 1
  }
  
  #Extracting the probabilities for transmission pairs and non-transmission pairs
  pT <- probs[!is.na(probs$p) & probs$truth == TRUE, "p"]
  pF <- probs[!is.na(probs$p) & probs$truth == FALSE, "p"]
  
  #Only executing function if there are non-missing values for the probabilities
  if(sum(is.na(probs[, "p"])) != nrow(probs)){
    
    #Finding the number of individual outbreaks
    nOutbreaks = sum(!is.na(unique(c(probs$outbreakID.1, probs$outbreakID.2))))
    #Finding the number of cases
    nCases = length(unique(probs$individualID.1)) + 1
    
    #Finding the AUC from the ROC curve
    aucVal <- as.numeric(pROC::auc(response = probs$truth, predictor = probs$p))
    
    #Ranking the probabilities for each possible infector
    #Ties are set to the minimum rank of that group
    probs <- (probs
              %>% group_by(individualID.2)
              %>% arrange(desc(p))
              %>% mutate(pRank = rank(desc(p), ties.method = "min"))
              %>% ungroup()
    )
    
    #Finding the median probabilities for transmission and non-transmission pairs
    medT <- median(pT)
    medF <- median(pF)
    
    #Finding rank and probability for the true infector for each infectee
    trueInfector <- probs %>% filter(truth == TRUE)
    trueInfector <- (trueInfector
                     %>% mutate(pInfector = p,
                                rankInfector = pRank)
                     %>% select(individualID.2, pInfector, rankInfector)
    )
    #Finding the number of possible infectors per infectee
    nInfectors <- (probs
                   %>% group_by(individualID.2)
                   %>% dplyr::summarize(nInfectors = n())
    )
    #Finding the maximum probability for each infectee
    allInd <- (probs
               %>% full_join(trueInfector, by = "individualID.2")
               %>% full_join(nInfectors, by = "individualID.2")
               #Calculating the top 5 and 10 percent of possible infectors
               %>% mutate(top5p = ceiling(nInfectors * 0.05),
                          top10p = ceiling(nInfectors * 0.1),
                          top15p = ceiling(nInfectors * 0.15),
                          top20p = ceiling(nInfectors * 0.2),
                          top25p = ceiling(nInfectors * 0.25),
                          top50p = ceiling(nInfectors * 0.5),
                          pPercent = pRank / nInfectors)
               %>% filter(!duplicated(individualID.2))
    )
    
    #Proportion of times that true infector has highest non-zero probability
    nTot <- sum(!is.na(allInd$rankInfector))
    pCorrect <- sum(allInd$rankInfector == 1 & allInd$pInfector != 0, na.rm = TRUE) / nTot
    #Proportion of times that the true infector is in the top n possible infectors
    pTop5 <- sum(allInd$rankInfector <= allInd$top5p & allInd$pInfector != 0, na.rm = TRUE) / nTot
    pTop10 <- sum(allInd$rankInfector <= allInd$top10p & allInd$pInfector != 0, na.rm = TRUE) / nTot
    pTop15 <- sum(allInd$rankInfector <= allInd$top15p & allInd$pInfector != 0, na.rm = TRUE) / nTot
    pTop20 <- sum(allInd$rankInfector <= allInd$top20p & allInd$pInfector != 0, na.rm = TRUE) / nTot
    pTop25 <- sum(allInd$rankInfector <= allInd$top25p & allInd$pInfector != 0, na.rm = TRUE) / nTot
    pTop50 <- sum(allInd$rankInfector <= allInd$top50p & allInd$pInfector != 0, na.rm = TRUE) / nTot
    #Average rank of true infector
    avgRank <- mean(allInd$rankInfector, na.rm = TRUE)
    
    
    #Calculating sensitivity and specificity
    #Threshold = 1/average number of infectors
    threshold <- 1 / mean(nInfectors$nInfectors)
    classData <- (probs
                  %>% select(p, truth)
                  %>% mutate(est = as.factor(p > threshold),
                             truth = as.factor(truth))
    )
    
    table(classData$truth, classData$est)
    
    sens <- sensitivity(classData$est, classData$truth, positive = "TRUE", negative = "FALSE")
    spec <- specificity(classData$est, classData$truth, positive = "TRUE", negative = "FALSE")
    ppv <- posPredValue(classData$est, classData$truth, positive = "TRUE", negative = "FALSE")
    npv <- negPredValue(classData$est, classData$truth, positive = "TRUE", negative = "FALSE")

    
  }else{
    aucVal <- NA
    sens <- NA
    spec <- NA
    ppv <- NA
    npv <- NA
    medT <- NA
    medF <- NA
    pCorrect <- NA
    avgRank <- NA
    pTop5 <- NA
    pTop10 <- NA
    pTop15 <- NA
    pTop20 <- NA
    pTop25 <- NA
    pTop50 <- NA
  }
  
  return(cbind.data.frame(nOutbreaks, nCases, aucVal, sens, spec, ppv, npv,
                          pCorrect, pTop5, pTop10, pTop15, pTop20, pTop25, pTop50,
                          avgRank, medT, medF, stringsAsFactors = FALSE))
}
