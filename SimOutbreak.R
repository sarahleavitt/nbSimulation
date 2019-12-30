#Sarah V. Leavitt
#Boston University Dissertation
#Simulation

##################################################################################
# This program simulates an outbreak with one or multiple transmission chains
##################################################################################


#### Inputs  ####

#neg = effective population size times pathogen generation time in years
#pi = proportion of cases sampled
#off.r = first parameter of the negative binomial distribution for offspring number
  #(if off.p = 0.5 then this value is the reproductive number)
#off.p = second parameter of the binomial distribution for offspring number (default is 0.5)
#w.shape = shape parameter of the gamma probability density function representing the generation time in years
  #(time from infection of primary case to infection of secondary case)
#w.scale = scale parameter of the gamma probability density function representing the generation time in years
#w.shift = number of years to shift the gamma distribution of the generation time
#ws.shape = shape parameter of the gamma probability density function representing the sampling time in years
#(time from DNA sampling of primary case to DNA sampling of secondary case)
#ws.scale = scale parameter of the gamma probability density function representing the sampling time in years
#ws.shift = number of years to shift the gamma distribution of the sampling time
#sampleSize = total number of cases to be involved in the outbreak/outbreaks
#startDate = year to start the outbreak
#time = number of years the outbreak should run
#multOutbreaks = do you want multiple simultaneous outbreaks (default is FALSE)
#length = length of sequence to simulate
#rate = mutation rate in mutations/site/year
#rootseq = a vector of length 1 containing the root sequence
  #if NULL, root sequence is randomly generated


#### Outputs ####

#List containing two dataframes
  #1: indData with the information about the individuals in the outbreak
  #2: pairData with information on all pairs of individuals (including SNP distance)


############# Function to simulate the outbreak ###############

simOutbreak <- function(neg, pi, off.r, off.p, w.scale, w.shape, w.shift,
                        ws.scale, ws.shape, ws.shift, sampleSize,
                        startDate = 1980, time, multOutbreaks = FALSE,
                        length, rate, rootseq = NULL){
  
  
  #### Simulate single outbreak #####
  
  #If multOutbreaks = FALSE, simulating single outbreak with the desired sample size
  if(multOutbreaks == FALSE){
    
    #Force outbreaks to be more than 2 people so genetic sequence simulation works
    nSampled <- 0
    while(nSampled < 2){
      simu <- simulateOutbreakS(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                                w.shape = w.shape, w.scale = w.scale, w.shift = w.shift,
                                ws.shape = ws.shape, ws.scale = ws.scale, ws.shift = ws.shift,
                                nSampled = sampleSize,
                                dateStartOutbreak = startDate, dateT = startDate + time)
      #Extracting the transmission information
      indData <- as.data.frame(extractTTree(simu)$ttree)
      nSampled <- nrow(indData)
    }
    
    indData <- (indData
                #Creating an ID variable for the outbreak and the individual
                %>% mutate(outbreakID = 1,
                           individualID = 10000 * outbreakID + 1:nrow(.),
                           infector = 10000 * outbreakID + V3)
                #Formating the dates
                %>% mutate(infectionDate = date_decimal(V1),
                           sampleDate = date_decimal(V2))
                %>% select(-V1, -V2, -V3)
                %>% arrange(individualID)
    )
    
    #Generating sequences
    #Extract the phylogenetic tree and convert it to a phylo object
    p <- phyloFromPTree(extractPTree(simu))
    #Simulating the sequences
    seq <- simSeq(p, l = length, rate = rate, rootseq = rootseq)
    names(seq) <- 10000 + as.numeric(names(seq))
    #Making the sequences into a DNAbin object
    seqBin <- as.DNAbin(seq)
  }
  
  
  
  #### Simulate multiple outbreaks ####
  
  startDate = startDate
  endDate = startDate + time
  
  #If more than one outbreak is desired, simulate multiple outbreaks and combine together
  if (multOutbreaks == TRUE){
    
    #### First outbreak ####
    
    #Simulate initial outbreak (to initialize objects incorrect forms)
    #Force outbreaks to be more than 2 people so sequence simulation works
    nSampled <- 0
    while(nSampled < 2){
      simu <- simulateOutbreakS(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                               w.shape = w.shape, w.scale = w.scale, w.shift = w.shift,
                               ws.shape = ws.shape, ws.scale = ws.scale, ws.shift = ws.shift,
                               dateStartOutbreak = startDate, dateT = endDate)
      #Extracting the transmission information
      indData <- as.data.frame(extractTTree(simu)$ttree)
      nSampled <- nrow(indData)
    }
    
    indData <- (indData
                #Creating an ID variable for the outbreak and the individual
                %>% mutate(outbreakID = 1,
                           individualID = 10000 * outbreakID + 1:nrow(.),
                           infector = 10000 * outbreakID + V3)
                #Formating the dates
                %>% mutate(infectionDate = date_decimal(V1),
                           sampleDate = date_decimal(V2))
                %>% select(-V1, -V2, -V3)
                %>% arrange(individualID)
    )
    
    #Generating sequences
    #Extract the phylogenetic tree and convert it to a phylo object
    p <- phyloFromPTree(extractPTree(simu))
    #Simulating the sequences
    seq <- simSeq(p, l = length, rate = rate, rootseq = rootseq)
    names(seq) <- 10000 + as.numeric(names(seq))
    #Making the sequences into a DNAbin object
    seqBin <- as.DNAbin(seq)
    
    
    #### Additional Outbreaks ####
    
    #Initilize counters
    size <- nrow(indData)
    i <- 2
    
    while(size < sampleSize){
      
      #Have the outbreaks vary their start and end times by adding or subtracting
      #half of the time for the outbreak from fixed start and end date.
      nSampled <- 0
      while(nSampled < 2){
        dateShift <- round(runif(1, -time/2, time/2))
        simu <- simulateOutbreakS(neg = neg, pi = pi, off.r = off.r, off.p = off.p,
                                 w.shape = w.shape, w.scale = w.scale, w.shift = w.shift,
                                 ws.shape = ws.shape, ws.scale = ws.scale, ws.shift = ws.shift,
                                 dateStartOutbreak = startDate + dateShift,
                                 dateT = endDate + dateShift)
        #Extracting the transmission information
        indDataTemp <- as.data.frame(extractTTree(simu)$ttree)
        nSampled <- nrow(indDataTemp)
      }

      indDataTemp <- (indDataTemp
                  #Creating an ID variable for the outbreak and the individual
                  %>% mutate(outbreakID = i,
                             individualID = 10000 * outbreakID + 1:nrow(.),
                             infector = 10000 * outbreakID + V3)
                  #Formating the dates
                  %>% mutate(infectionDate = date_decimal(V1),
                             sampleDate = date_decimal(V2))
                  %>% select(-V1, -V2, -V3)
                  %>% arrange(individualID)
      )

      #Generating sequences
      #Extract the phylogenetic tree and convert it to a phylo object
      p <- phyloFromPTree(extractPTree(simu))
      #Simulating the sequences
      seq <- simSeq(p, l = length, rate = rate, rootseq = rootseq)
      names(seq) <- 10000 * i + as.numeric(names(seq))
      #Making the sequences into a DNAbin object
      seqBinTemp <- as.DNAbin(seq)
      
      #Combining with other info
      indData <- bind_rows(indData, indDataTemp)
      seqBin <- rbind(seqBin, seqBinTemp)
      
      size = nrow(indData)
      i <- i + 1
    }
    
  }
  
  
  #### Combine all datasets and find pair info ####
  
  #This function finds the number of SNPs different for each pair
  snpDist <- dist.gene(seqBin)
  #Using melt to create a dataset with one row per pair and the SNP distance
  snpDistDf <- reshape2::melt(as.matrix(snpDist), varnames = c("individualID.1", "individualID.2"))
  snpDistDf <- snpDistDf %>% rename(snpDist = value)
  
  #Finding all pairs of IDs (order matters)
  pairs <- expand.grid(indData$individualID, indData$individualID)
  pairs2 <- pairs %>% rename(individualID.1 = Var1,  individualID.2 = Var2)
  
  #Creating a dataset of the transmission pairs
  pairData <- (indData
               %>% select(infector, individualID)
               %>% rename(individualID.1 = infector,
                          individualID.2 = individualID)
               #Creating a variable that defines these pairs as transmission pairs
               #Creating a rowID to merge with pairs dataset
               %>% mutate(transmission = TRUE)
               #Combining with full permutation to get all pairs
               %>% full_join(pairs2, by = c("individualID.1", "individualID.2"))
               #Combining with the SNP distance matrix
               %>% full_join(snpDistDf, by = c("individualID.1", "individualID.2"))
               #Making the transmission variable FALSE for all transmission pairs
               %>% replace_na(list(transmission = FALSE))
               #Joining with the individual data to get the information about the infector
               #Making ID.1 the infector information and ID.2 the infectee information
               %>% full_join(indData, by = c("individualID.1" = "individualID"))
               %>% full_join(indData, by = c("individualID.2" = "individualID"),
                             suffix = c(".1", ".2"))
               #Creating time difference variables
               #And a variable that indicates if the cases are in the same outbreak
               %>% mutate(sampleDiff = round(as.numeric(difftime(sampleDate.2, sampleDate.1, units = "days")), 0),
                          infectionDiff = round(as.numeric(difftime(infectionDate.2, infectionDate.1, units = "days")), 0),
                          sameOutbreak = outbreakID.1 == outbreakID.2)
               #Creating an edgeID
               %>% unite(edgeID, individualID.1, individualID.2, remove = FALSE)
               %>% select(edgeID, individualID.1, individualID.2,  outbreakID.1, outbreakID.2,
                          transmission, snpDist, sampleDiff, infectionDiff,
                          infectionDate.1, infectionDate.2, sampleDate.1, sampleDate.2, sameOutbreak)
               %>% filter(individualID.1 != individualID.2)
  )
  
  return(list(indData, pairData))
  
}
