
# R code to generate simulated data for the chapter 
# What makes for reliable kinship assignments; how sample design, genetic markers and programs used shape assignment accuracies.

#create families 

library("dplyr", lib.loc="/usr/local/R/3.5.1-gcc7/lib64/R/library")
binomparent <- vector()
parentlist = list()
offlist = list()
dataframeoff<- (data.frame(matrix(ncol = 40000, nrow = 20)))
dataframeparent <- (data.frame(matrix(ncol = 40000, nrow = 2)))

# Create a truncated normal distribution which max of 1 and min of 0
# Note it does actually calculate normal distribution values below 0 and above 1 it just calculates extra values and filters out those 
# too low or high. Same effect

normallimit <- function(n, m, s, lwr, upr, nnorm) {
  samp <- rnorm(nnorm, m, s)
  samp <- samp[samp >= lwr & samp <= upr]
  if (length(samp) >= n) {
    return(sample(samp, n))
  }  
  stop(simpleError("Not enough values to sample from"))
}

# full dataframe has 40000 columns (40,000 alleles, 20,000 loci) and 5500 rows (500 adult parents + 5000 offspring)
fulldataframe <-  data.frame(matrix(ncol = 40000, nrow = 5500))
mysamp <- normallimit(n=40000, m=0.5, s=0.20, lwr=0, upr=1, nnorm=60000)



for (k in c(1:250)) {
  
  for (p in c(1:2)) {
    
    for (i in c(1:40000)) {
      
      binomparent[i] <- (rbinom(1,1,(mysamp[i]))+1)
      
      
    }
    
    dataframeparent[p,] <- binomparent
    
  }
  
  
  for (j in c(1:20)) {
    # This part is used to randomise which alleles from which parent is passed onto each child. This is achieved through a simple binom probability 
  
    #dataframeparent <-  as.data.frame(matrix(unlist(parentlist),nrow =2))
    
    alleleintpar1 <- as.numeric(dataframeparent[1,])
    allleleintpar2 <- as.numeric(dataframeparent[2,])
    Allele1par1 <- alleleintpar1[seq(1, length(alleleintpar1), 2)]
    Allele2par1 <- alleleintpar1[seq(2, length(alleleintpar1), 2)]
    Allele1par2 <- allleleintpar2[seq(1, length(allleleintpar2), 2)]
    Allele2par2 <- allleleintpar2[seq(2, length(allleleintpar2), 2)]
    binomprob1 <- rbinom(length(Allele1par1),1,0.5)
    binomprob2 <- rbinom(length(Allele1par2),1,0.5) 
    gampar1 <- ifelse(binomprob1[1:length(binomprob1)]==0, Allele1par1[seq(1, length(Allele1par1),1)], Allele2par1[seq(1, length(Allele2par1),1)])   
    gampar2 <- ifelse(binomprob2[1:length(binomprob2)]==0, Allele1par2[seq(1, length(Allele1par2),1)], Allele2par2[seq(1, length(Allele2par2),1)])
    
    
    #produces offspring but alleles don't segregate on their own (allele from gamete 1 is always allele 1 for offspring)
    off <-  c(rbind(gampar1, gampar2))
    
    #rerandomises alleles across allele 1 and 2
    binomprob3 <-rbinom(length(Allele1par1),1,0.5)
    offrandom <- ifelse(binomprob3[1:length(binomprob3)]==0, off[seq(1, length(off), 2)], off[seq(2, length(off), 2)])
    
    offrandom2 <- ifelse(binomprob3[1:length(binomprob3)]==1, off[seq(1, length(off), 2)], off[seq(2, length(off), 2)])
    
    offfullrandom <-  c(rbind(offrandom, offrandom2))
    
    dataframeoff[j,] <- offfullrandom
    
  }
  
  
  colnames(dataframeparent) <- colnames(dataframeoff)
  colnames(fulldataframe) <- colnames(dataframeoff)
  dataframeparentpair <- rbind(dataframeparent,dataframeoff)
  
  
  fulldataframe[(1+(k*(nrow(dataframeparentpair)))-nrow(dataframeparentpair)):((nrow(dataframeparentpair))*k), ] <- dataframeparentpair[1:22, ]
  
  
}


dataframe_names <-  data.frame(matrix(ncol = 1, nrow = 5500))
#rename stuff
# Downstream for Colony analysis, a consistent naming scheme allows for mass analysing of results in ways that manually relooking up families would
# be too time consuming. 
# Generate a key of values for names e.g., P for parent, O for offspring, then the number the ind is, if a male then MF then the number of the partner
# if female FM then the number of the partner, then G for which family group it belongs to (e.g., 1-250 for families created)
for (i in c(1:250)) {
  
  dataframe_names[(((i*22)-22)+1),1] <- paste("P", (((i*22)-22)+1) , "MF",(((i*22)-22)+2),"G", i , sep="" )
  dataframe_names[(((i*22)-22)+2),1] <- paste("P", (((i*22)-22)+2) , "FM",(((i*22)-22)+1),"G", i , sep="" )
  
  for(j in c(1:20)) {
    
    dataframe_names[((((i*22)-22)+2)+j),1] <- paste("O", ((((i*22)-22)+2)+j),"OM",(((i*22)-22)+1),"F",(((i*22)-22)+2),"G",i,sep="")
    
    
    
  }
}
# repeat this for a random subsample of adults to rebreed for half siblings
rownames(fulldataframe)<- dataframe_names$matrix.ncol...1..nrow...5500.
parents <- subset(fulldataframe, grepl("P",rownames(fulldataframe)))
maleparents <- subset(parents, grepl("MF",rownames(parents)))
femaleparents <- subset(parents, grepl("FM",rownames(parents)))
dataframefull <- data.frame(matrix(ncol = 40000, nrow = 1000))
dataframenames2 <-  data.frame(matrix(ncol = 1, nrow = 1000))
dataframeoff2<- (data.frame(matrix(ncol = 40000, nrow = 20)))

for (p in c(1:50)) {
  
  males <- maleparents[sample(nrow(maleparents),50, replace = F),]
  
  females <- femaleparents[sample(nrow(femaleparents),50, replace = F),]     
 
  alleleintpar1 <- as.numeric(males[p,])
  allleleintpar2 <- as.numeric(females[p,])

    #convert parent vectors to separate alleles
    Allele1par1 <- alleleintpar1[seq(1, length(alleleintpar1), 2)]
    Allele2par1 <- alleleintpar1[seq(2, length(alleleintpar1), 2)]
    
    Allele1par2 <- allleleintpar2[seq(1, length(allleleintpar2), 2)]
    Allele2par2 <- allleleintpar2[seq(2, length(allleleintpar2), 2)]

 
  for (i in c(1:20)) {
    
        

    binomprob1 <- rbinom(20000,1,0.5)
    binomprob2 <- rbinom(20000,1,0.5) 
    
    
    gampar1 <- ifelse(binomprob1[1:length(binomprob1)]==0, Allele1par1[seq(1, length(Allele1par1),1)], Allele2par1[seq(1, length(Allele2par1),1)])   
    gampar2 <- ifelse(binomprob2[1:length(binomprob2)]==0, Allele1par2[seq(1, length(Allele1par2),1)], Allele2par2[seq(1, length(Allele2par2),1)])
    
  
    off <- c(rbind(gampar1, gampar2))
    binomprob3 <-rbinom(20000,1,0.5)
    offrandom <- ifelse(binomprob3[1:length(binomprob3)]==0, off[seq(1, length(off), 2)], off[seq(2, length(off), 2)])
    offrandom2 <- ifelse(binomprob3[1:length(binomprob3)]==1, off[seq(1, length(off), 2)], off[seq(2, length(off), 2)])
    offfullrandom2 <-  c(rbind(offrandom, offrandom2))
    
    dataframeoff2[i,] <- offfullrandom2
    
  }

  #make names
  
  sub1_1 <- sub("M.*", "", rownames(males)[p])
  sub1_2 <- sub("P", "M", sub1_1)
  
  sub2_1 <- sub("F.*", "", rownames(females)[p])
  sub2_2 <- sub("P", "F", sub2_1)
  for (l in c(1:20)) {
    
    dataframenames2[((p*20)-20)+l,] <- paste("Off",l,sub1_2,sub2_2, sep="")
    
    
  }
  dataframefull[(((p*20)-20)+1):(p*20),] <-  dataframeoff2[1:20,]
}
rownames(dataframefull) <- dataframenames2$matrix.ncol...1..nrow...1000.
complete_simulation_datasethalf <- rbind(fulldataframe,dataframefull)


dataframeotherM <- (data.frame(matrix(ncol = 40000, nrow = 250)))
dataframeotherF <- (data.frame(matrix(ncol = 40000, nrow = 250)))
dataframenames3 <- data.frame(matrix(ncol = 1, nrow = 250))
dataframenames4 <- data.frame(matrix(ncol = 1, nrow = 250))
dataframeother <- (data.frame(matrix(ncol = 40000, nrow = 500)))
 binomotherM <- vector()
 binomotherF <- vector()

for (k in c(1:250)) {
      
    for (i in c(1:40000)) {
      
      binomotherM[i] <- (rbinom(1,1,(mysamp[i]))+1)
      
      
    }
    
  dataframeotherM[k,] <- binomotherM
  dataframenames3[k,] <- paste("Rand","MM",k,sep="")  
  }

for (k in c(1:250)) {
      
    for (i in c(1:40000)) {
      
      binomotherF[i] <- (rbinom(1,1,(mysamp[i]))+1)
      
      
    }
    
  dataframeotherF[k,] <- binomotherF
  dataframenames4[k,] <- paste("Rand","FF",k,sep="")  
  }

rownames(dataframeotherM) <- dataframenames3$matrix.ncol...1..nrow...250.
rownames(dataframeotherF) <- dataframenames4$matrix.ncol...1..nrow...250.

  dataframeother <- rbind(dataframeotherM,dataframeotherF)

complete_simulation_dataset<- rbind(complete_simulation_datasethalf,dataframeother)





rm(dataframe_names,dataframefull,dataframenames,dataframeoff,dataframeoff2,dataframeparent,dataframeparentpair,maleparents,femaleparents,offlist,parentlist,parents,Allele1par1,Allele1par2,Allele2par1,Allele2par2, binomparent,binomprob1,binomprob2,binomprob3,gampar1,gampar2,i,j,k,l,off,offfullrandom,offrandom,offrandom2,p,sub1_1,sub1_2,sub2_1,sub2_2,fulldataframe,alleleintpar1,allleleintpar2,males,females,dataframeother,dataframeotherM, dataframeotherF,dataframenames3, dataframenames4,binomotherM,binomotherF,complete_simulation_datasethalf)

# save environment as rdata to desired location
save.image("/home/group/murphylab/James/simulated_data/correctedfullsimulationdatanamed.RData")

#introduce NAs
# instead of rolling through every one of the 280 million cells and running a probability of success fail for introducing Nas (which would take
# forever) I created a vector of length 280 million for index numbers referencing the cells, randomly sampled the required NA and typing error percentages
# and changed the indext cell values which were sampled 
# note the rate is set at 2% here because there is a 50% chance that it is changed to its already existing value which gives an overal
#1%  error rate.
complete_simulation_matrix <- as.matrix(complete_simulation_dataset)
sampleforNA <- sample(1:280000000, (0.1*280000000), replace=F)
samplefortyping <- sample(1:280000000, (0.02*280000000), replace=F)


for (j in c(1:(0.02*280000000))) {
  
  
i <- samplefortyping[j]
    
    complete_simulation_matrix[i] <- ifelse(rbinom(1,1,0.5)==1,1,2)
    
}


save.image("/home/group/murphylab/James/simulated_data/typingerrorscomplete.RData")


for (k in c(1:(0.05*280000000)))     {
  
  i <- sampleforNA[k]
  complete_simulation_matrix[i] <- '0'
  

}
  

dim(complete_simulation_matrix) <- c(7000,40000)  
rownames(complete_simulation_matrix) <- rownames(complete_simulation_dataset)
colnames(complete_simulation_matrix) <- colnames(complete_simulation_dataset)
rm(complete_simulation_dataset,mysamp, samplefortyping,sampleforNA,i,j,k)

# convert matrix of values to a genid object for genetic analysis
# geno to genind
library("adegenet", lib.loc="/usr/local/R/3.5.1-gcc7/lib64/R/library")
library("pegas", lib.loc="/usr/local/R/3.5.1-gcc7/lib64/R/library")

complete_simulation_matrix <- as.data.frame(complete_simulation_matrix)

allnames <- rownames(complete_simulation_matrix)
complete_simulation_matrix <- alleles2loci(complete_simulation_matrix, ploidy = 2, population = NULL,
             phased = FALSE)
save.image("/home/group/murphylab/James/simulated_data/dfallindivsfullsimerrors.RData")

complete_simulation_matrix[] <- lapply(complete_simulation_matrix, function(x) replace(x, grep("[0]", x), NA))

genidallindivs <- df2genind(complete_simulation_matrix, sep = "/", ncode = NULL, ind.names = NULL,
  loc.names = NULL, pop = NULL, NA.char = "NA", ploidy = 2,
  type = c("codom"), strata = NULL, hierarchy = NULL)

rm(complete_simulation_matrix)
save.image("/home/group/murphylab/James/simulated_data/genidallindivsfullsimerrors.RData")
