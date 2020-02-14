# nbSimulation

This directory contains a set of programs used to simulate outbreak data.
These programs are used by other programs in the npPaper1, nbPaper2, and nbPaper3
directories which create the outputs for the papers developing the naive Bayes transmission
method.


### SimOutbreak.R

This program contains the function "simOutbreak" which simulates an infectious
disease outbreak from various input parameters including reproductive number,
generation interval, and sample size. The outbreak is simulated using an adapted 
version of the "simulateOutbreak" function from the TransPhylo package (see below). 
The outbreak could be composed of one transmission chains starting with one case 
or multiple parallel transmission chains. The function also simulates pathogen 
genetic sequences using the phangorn package. The output is a list with 
individidual-level and pair-level datasets where the pathogen genetic sequence data 
are captured in SNP distances between pairs.



### SimulateOutbreakS.R

This program contains a set of functions taken directly from the TransPhylo package 
(https://github.com/xavierdidelot/TransPhylo/tree/master/R) but adapting them to
use a shifted gamma distribution for the generation and sampling times. The main
function of interest from this program is "simulateOutbreakS" which is the same as the
TransPhylo function "simulateOutbreak" but adding the shift parameter. This function
is used in SimOutbreak.R.



### SimCovariates.R

This program contains the function "simCovariates" which simulates four covariates 
at the individual-level and transforms them into pair-level covariates. The details 
of the covariates is described at the beginning of the program. The function also 
creates an "observationDate" that is a copy of the variable which contains the 
date to be used in all downstream analyses (creating ordered pairs, estimating
the reproductive number, estimating the generation interval).



### SimEvaluate.R

This program contains the function "simEvaluate" which evaluates the performance 
of the naive Bayes transmission probabilities by calculating various metrics 
including the AUC, the percentage of time the true infector has the highest 
probability, the precentage of time the true infector is in the top n% of 
the possible infectors, and the sensitivity, specificity, ppv, and npv at one 
possible threshold.



### PerformRandom.R

This program contains the function "performRandom" which simulates random 
probabilities for a dataset of ordered pairs of cases. The output is a data frame 
with a label, the pairID variable, the random probability, the scaled probability,
and the rank of the probability. 



### PerformInterval.R

This program contains the function "performInterval" which estimates the probability 
pair is linked based on the time between the obervation dates and a prespecified
generation/serial interval. The output is a data frame with a label, the pairID 
variable, the random probability, the scaled probability, and the rank of the 
probability. 







