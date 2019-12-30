#Sarah V. Leavitt
#Boston University Dissertation
#Simulation


## Adapted from TransPhylo ##

########################################################################################
# This programs adapts the SimulateOutbreak function from the TransPhylo package
# so that it draws the generation and sampling times from a shifted gamma distribution
########################################################################################

#' Simulate an outbreak
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param neg the within-host effective population size (Ne) times generation duration in years (g)
#' @param nSampled number of sampled infected individuals, or NA for any
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time in years
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time in years
#' ## Added by SVL ##
#' @param w.shift amount to shift the Gamma probability density function for both generation time
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time in years
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time in years
#' ## Added by SVL ##
#' @param ws.shift amount to shift the Gamma probability density function for both sampling time
#' @param dateStartOutbreak Date when index case becomes infected
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return Combined phylogenetic and transmission tree
#' @examples
#' plotCTree(simulateOutbreak())
#' @export

simulateOutbreakS = function(off.r=1,off.p=0.5,neg=0.25,nSampled=NA,pi=0.5,w.shape=2,w.scale=1,w.shift=0,
                             ws.shape=w.shape,ws.scale=w.scale,ws.shift=w.shift,dateStartOutbreak=2000,dateT=Inf) {
  #Create a transmission tree with nSampled infected sampled individuals
  nsam<-0
  nh<-0
  rejected=-1
  while (is.na(nSampled)||nsam!=nSampled) {
    ttree=NULL
    while (is.null(ttree)) {
      mtt<-makeTTreeS(off.r,off.p,pi,w.shape,w.scale,w.shift,ws.shape,ws.scale,ws.shift,dateT-dateStartOutbreak,nSampled)
      rejected=rejected+1
      ttree<-mtt$ttree
      probttree<-mtt$prob
      if (mtt$pruned>0) {
        dateStartOutbreak=dateStartOutbreak+mtt$pruned
        cat(sprintf('Note that simulated outbreak was pruned: in order to have %d sampled by present date %f, the start date was set to %f\n',nSampled,dateT,dateStartOutbreak))
      }
    }
    nsam<-length(which(!is.na(ttree[,2])))
    nh=nrow(ttree)-nsam
    if (is.na(nSampled)) nSampled=nsam
  }
  if (rejected>0) cat(sprintf('Note that rejection sampling was used %d times to simulate outbreak with %d sampled individuals\n',rejected,nSampled))
  n<-nsam+nh
  
  #Create a within-host phylogenetic tree for each infected host
  wtree<-vector('list',n)
  probwithin=0
  for (i in (1:n)) {
    if (is.na(ttree[i,2])) {times<-c(           ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    else {times<-c(ttree[i,2],ttree[which(ttree[,3]==i),1])-ttree[i,1]}
    a<-withinhost(times,neg)
    wtree[[i]]=a$nodes
    probwithin=probwithin+a$prob
  }
  
  #Glue these trees together
  truth<-.glueTrees(ttree,wtree)
  truth[,1]<-truth[,1]+dateStartOutbreak
  return(list(ctree=truth,nam=mtt$nam,probttree=probttree,probwithin=probwithin))
}  



#' Simulate a transmission tree
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' ## Added by SVL ##
#' @param w.shift amount to shift the Gamma probability density function for both generation time
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time in years
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time in years
#' ## Added by SVL ##
#' @param ws.shift amount to shift the Gamma probability density function for both sampling time
#' @param maxTime Duration of simulation (can be Inf)
#' @param nSampled Number of sampled individuals (can be NA for any)
#' @return A N*3 matrix in the following format with one row per infected host,
#'  first column is time of infection, second column is time of sampling, third column is infector
#' @export

makeTTreeS <-function(off.r,off.p,pi,w.shape,w.scale,w.shift,ws.shape,ws.scale,ws.shift,maxTime,nSampled) { 
  ttree<-matrix(0,1,3)
  prob<-0
  todo<-1
  while (length(todo)>0) {
    rand<-runif(1)
    if (rand<pi) {
      #This individual is sampled
      prob<-prob+log(pi)
      draw<-rgamma(1,shape=ws.shape,scale=ws.scale) + ws.shift
      if (ttree[todo[1],1]+draw<maxTime)
        ttree[todo[1],2]<-ttree[todo[1],1]+draw
      else ttree[todo[1],2]<-NA
      prob<-prob+log(dgamma(draw - ws.shift,shape=ws.shape,scale=ws.scale))}
    else {
      #This individual is not sampled
      prob<-prob+log(1-pi)
      ttree[todo[1],2]<-NA}
    
    offspring<-rnbinom(1,off.r,1-off.p)
    prob<-prob+log(dnbinom(offspring,off.r,1-off.p))
    if (offspring>0) {
      for (i in 1:offspring) {
        draw<-rgamma(1,shape=w.shape,scale=w.scale) + w.shift
        prob<-prob+log(dgamma(draw - w.shift,shape=w.shape,scale=w.scale))
        if (ttree[todo[1],1]+draw>maxTime) next
        ttree<-rbind(ttree,c(ttree[todo[1],1]+draw,0,todo[1]))
        todo<-c(todo,nrow(ttree))
      }
    }
    todo<-todo[-1] 
  }
  
  #Prune down number of sampled individuals if needed (only for ongoing outbreaks)
  pruned=0
  samTimes=ttree[which(!is.na(ttree[,2])),2]
  if (!is.na(nSampled)&&maxTime<Inf&&nSampled<length(samTimes)) {
    samTimes=sort(samTimes)
    cutoff=(samTimes[nSampled]+samTimes[nSampled+1])/2
    pruned=maxTime-cutoff
    ttree[ttree[,2]>cutoff,2]=NA
  }
  
  #Remove infected individuals who are not sampled and have no children
  while (TRUE) {
    if (nrow(ttree)==1 && is.na(ttree[1,2])) {return(list(ttree=NULL,prob=NULL,pruned=0))} #Nothing left
    torem=c()
    for (i in 1:nrow(ttree)) {if (is.na(ttree[i,2])&&length(which(ttree[,3]==i))==0) {torem=c(torem,i)}}
    if (length(torem)==0) {break}
    ttree<-ttree[setdiff(1:nrow(ttree),torem),,drop=FALSE]
    for (i in 1:nrow(ttree)) {ttree[i,3]=ttree[i,3]-length(which(torem<ttree[i,3]))}
  }  
  
  #Remove root if not sampled and only has one child
  while (TRUE) {
    if (nrow(ttree)==1 && is.na(ttree[1,2])) {return(list(ttree=NULL,prob=NULL,pruned=0))} #Nothing left
    if (is.na(ttree[1,2])&&length(which(ttree[,3]==1))==1) {
      ttree[,2]<-ttree[,2]-ttree[2,1]
      ttree[,1]<-ttree[,1]-ttree[2,1]
      ttree<-ttree[-1,,drop=FALSE]
      ttree[,3]=ttree[,3]-1
    } else {break}
  }
  
  #Reorder so that sampled hosts are first
  order<-c(which(!is.na(ttree[,2])),which(is.na(ttree[,2])))
  invorder=1:length(order);invorder[order]=1:length(order)
  ttree<-ttree[order,,drop=FALSE]
  ttree[ttree[,3]>0,3]=invorder[ttree[ttree[,3]>0,3]]
  return(list(ttree=ttree,nam=sprintf('%d',seq(1:length(which(!is.na(ttree[,2]))))),prob=prob,pruned=pruned))
}



#Glue together some within-host trees using a transmission tree 
.glueTrees = function(ttree,wtree)  {
  nsam <- length(which(!is.na(ttree[,2])))
  n <- nrow(ttree)
  nh <- n-nsam
  no <- nsam + 1 
  ma<-0
  for (i in (1:n)) {ma<-max(ma,nrow(wtree[[i]])/2)} 
  lab<-matrix(0,n,ma)
  labh<-rep(0,n)
  for (i in (n:1)) {#Consider all hosts 
    ni <- nrow(wtree[[i]])/2 
    for (j in (1:ni)) { 
      lab[i,j] <- no 
      no <- no + 1 
    } 
    labh[i] <- no-1 
  } 
  leaves <- cbind() 
  intnodes <- cbind() 
  for (i in (1:n)) { 
    tree <- wtree[[i]]
    ni <- nrow(wtree[[i]])/2 
    tree[ ,1] <- tree[ ,1] + ttree[i,1];#Add infection time to all nodes 
    if (!is.na(ttree[i,2])) leaves <- rbind(leaves,tree[1, ]) 
    tree <- tree[(ni+1):nrow(tree), ,drop=FALSE] #keep only internal nodes of current within-host tree
    f <- which( ttree[ ,3] == i )#Infected by current nodes 
    for (j in (1:nrow(tree))) for (k in (2:3)) { 
      if (tree[j,k] == 0)  {next} 
      if (!is.na(ttree[i,2])&&tree[j,k] == 1)  {tree[j,k] <- i}
      else if (tree[j,k] <= ni)  { 
        tree[j,k] <- labh[f[tree[j,k]-!is.na(ttree[i,2])]] 
      } else { tree[j,k] <- lab[i,tree[j,k]-ni] 
      } 
    } 
    intnodes <- rbind(tree,intnodes) 
  } 
  tree <- rbind(leaves,intnodes) 
  ind = order(tree[(nsam+1):nrow(tree),1],decreasing=TRUE)
  for (i in (nsam+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > nsam) tree[i,j] <- nsam + which( ind == tree[i,j]-nsam ) 
  tree <- tree[c(1:nsam,nsam + ind), ] 
  tree <- cbind(tree,.computeHost(tree)) 
  return(tree)
} 



#Build vector 'host' indicating in which host each node of the ctree is found 
.computeHost = function(ctree)  {
  fathers <- rep(NA, nrow(ctree)) 
  fathers[ctree[ ,2] + 1] <- 1:nrow(ctree) 
  fathers[ctree[ ,3] + 1] <- 1:nrow(ctree) 
  fathers <- fathers[-1] 
  host <- rep(0, nrow(ctree)) 
  nsam <- sum(ctree[ ,2] == 0&ctree[ ,3] == 0) 
  for (i in 1:nsam) { 
    j <- i 
    while (1)  { 
      if (host[j]>0) print('Error: two leaves in same host')
      host[j] <- i 
      j <- fathers[j] 
      if (ctree[j,3] == 0) break 
    } 
  } 
  if (nsam==1) return(host)
  
  dispo=nsam+1
  f <- which( ctree[,3] == 0 & ctree[,2]>0 & (1:nrow(ctree))<nrow(ctree) ) #transmission events other than root
  for (i in 1:length(f)) { 
    j <- f[i] 
    tocol <- c() 
    while (ctree[fathers[j],3]>0&&fathers[j]<nrow(ctree))  { 
      tocol <- c(tocol,j) 
      j <- fathers[j] 
    } 
    if (host[j]==0) {host[j]=dispo;dispo=dispo+1}
    host[tocol] <- host[j] 
  } 
  return(host)
} 


# Simulates the within-host coalescent model
# @param times times at which N samples are taken(counted forward in time from infection time) 
# @param neg is the product of the within-host effective population size and the generation duration in days 
# @return array of size(2N)*3 where each row is a node,the first column indicate the date of the node and the
#last two columns indicate the two children. This array has internal nodes sorted in order of most recent to
#most ancient node(and remains so during the algorithm). The last node corresponds to infection time and only
#has one child 
withinhost = function(times,neg)  {
  prob <- 0 
  #log1overneg=log(1/neg)
  ind=order(times,decreasing=T);tim=times[ind]
  n <- length(tim) 
  nodes <- cbind(0,ind[1],0);#Start with one node at time 0 and with the first isolate connected to it 
  i <- 2 
  while (i <= n) {#Graft branches one by one 
    curt <- tim[i]
    fi <- which( nodes[ ,1] < curt );fi<-fi[1]
    trunc=0
    for (j in (seqML(fi,nrow(nodes))))  {
      trunc=trunc+(curt-nodes[j,1]) * (i-j) 
      curt <- nodes[j,1] 
    }
    r=rexpT(1/neg,trunc)
    prob=prob+r$prob
    r=r$x
    curt <- tim[i]#Current time:start with date of isolate and go back in time until coalescence happens 
    fi <- which( nodes[ ,1] < curt );fi<-fi[1]
    for (j in (seqML(fi,nrow(nodes))))  {
      if (r > (curt-nodes[j,1]) * (i-j))  { 
        r <- r-(curt-nodes[j,1]) * (i-j) 
        curt <- nodes[j,1] 
      } else { 
        curt <- curt-r/(i-j)#Found the time for grafting
        r <- 0 
        break 
      } 
    } 
    if (r>0) print('error: this shound not happen')#{next} 
    #Create new node 
    a <- nodes[ ,2:3];a[a >= j + n] <- a[a >= j + n] + 1;nodes[ ,2:3] <- a;#Renumbering according to table insertion in next line 
    nodes <- rbind(nodes[seqML(1,j-1), ],c(curt,ind[i],0),nodes[seqML(j,nrow(nodes)),]) 
    #Now choose on which branch to regraft amongst the branches alive at time curt 
    no <- j 
    side <- 2 
    #prob <- prob + log(1/(nrow(nodes)-j))
    w <- 1 + floor(runif(1) * (nrow(nodes)-j)) 
    while (w > 0)  { 
      no <- no + side-1 
      side <- 3-side 
      if (nodes[no,side + 1] <= n ||(nodes[no,side + 1] > n && nodes[nodes[no,side + 1]-n,1] > curt))  { 
        w <- w-1 
      } 
    } 
    nodes[j,3] <- nodes[no,side + 1] 
    nodes[no,side + 1] <- n + j 
    i <- i + 1 
  } 
  nodes <- rbind(matrix(0, nrow = n, ncol = 3),nodes) 
  nodes[1:n,1] <- times 
  return(list(nodes = nodes,prob = prob))
} 

seqML <- function(from, to, by=1) {if (from > to) integer(0) else seq.int(from, to, by)}

#truncated exponential distribution
rexpT = function(rate,trunc) {
  x=-log(1-runif(1)*(1-exp(-rate*trunc)))/rate
  prob=log(rate)-rate*x-log(1-exp(-rate*trunc))
  return(list(x=x,prob=prob))
}

withinhost = function(times,neg)  {
  prob <- 0 
  #log1overneg=log(1/neg)
  ind=order(times,decreasing=T);tim=times[ind]
  n <- length(tim) 
  nodes <- cbind(0,ind[1],0);#Start with one node at time 0 and with the first isolate connected to it 
  i <- 2 
  while (i <= n) {#Graft branches one by one 
    curt <- tim[i]
    fi <- which( nodes[ ,1] < curt );fi<-fi[1]
    trunc=0
    for (j in (seqML(fi,nrow(nodes))))  {
      trunc=trunc+(curt-nodes[j,1]) * (i-j) 
      curt <- nodes[j,1] 
    }
    r=rexpT(1/neg,trunc)
    prob=prob+r$prob
    r=r$x
    curt <- tim[i]#Current time:start with date of isolate and go back in time until coalescence happens 
    fi <- which( nodes[ ,1] < curt );fi<-fi[1]
    for (j in (seqML(fi,nrow(nodes))))  {
      if (r > (curt-nodes[j,1]) * (i-j))  { 
        r <- r-(curt-nodes[j,1]) * (i-j) 
        curt <- nodes[j,1] 
      } else { 
        curt <- curt-r/(i-j)#Found the time for grafting
        r <- 0 
        break 
      } 
    } 
    if (r>0) print('error: this shound not happen')#{next} 
    #Create new node 
    a <- nodes[ ,2:3];a[a >= j + n] <- a[a >= j + n] + 1;nodes[ ,2:3] <- a;#Renumbering according to table insertion in next line 
    nodes <- rbind(nodes[seqML(1,j-1), ],c(curt,ind[i],0),nodes[seqML(j,nrow(nodes)),]) 
    #Now choose on which branch to regraft amongst the branches alive at time curt 
    no <- j 
    side <- 2 
    #prob <- prob + log(1/(nrow(nodes)-j))
    w <- 1 + floor(runif(1) * (nrow(nodes)-j)) 
    while (w > 0)  { 
      no <- no + side-1 
      side <- 3-side 
      if (nodes[no,side + 1] <= n ||(nodes[no,side + 1] > n && nodes[nodes[no,side + 1]-n,1] > curt))  { 
        w <- w-1 
      } 
    } 
    nodes[j,3] <- nodes[no,side + 1] 
    nodes[no,side + 1] <- n + j 
    i <- i + 1 
  } 
  nodes <- rbind(matrix(0, nrow = n, ncol = 3),nodes) 
  nodes[1:n,1] <- times 
  return(list(nodes = nodes,prob = prob))
} 

seqML <- function(from, to, by=1) {if (from > to) integer(0) else seq.int(from, to, by)}

#truncated exponential distribution
rexpT = function(rate,trunc) {
  x=-log(1-runif(1)*(1-exp(-rate*trunc)))/rate
  prob=log(rate)-rate*x-log(1-exp(-rate*trunc))
  return(list(x=x,prob=prob))
}