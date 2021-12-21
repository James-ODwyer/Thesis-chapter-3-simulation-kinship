# R code to filter Macquarie perch data for the chapter 
# What makes for reliable kinship assignments; how sample design, genetic markers and programs used shape assignment accuracies.

library("ggplot2", lib.loc="~/R/win-library/3.5")
library("lattice", lib.loc="~/R/win-library/3.5")
library("MASS", lib.loc="~/R/win-library/3.5")
library("PopGenReport", lib.loc="~/R/win-library/3.5")
library("poppr", lib.loc="~/R/win-library/3.5")
library("dartR", lib.loc="~/R/win-library/3.5")
library("adegenet", lib.loc="~/R/win-library/3.5")
library("devtools", lib.loc="~/R/win-library/3.5")
library("diveRsity", lib.loc="~/R/win-library/3.5")
library("StAMPP", lib.loc="~/R/win-library/3.5")
library("radiator", lib.loc="~/R/win-library/3.5")
devtools::install_github('green-striped-gecko/dartR', dependencies = TRUE)



if (!require("devtools")) install.packages("devtools")
devtools::install_github("thierrygosselin/radiator")



setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset")
setwd("D:/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset")
glbase <- gl.read.dart(filename="Report_DMacq17-3037_3_moreOrders_SNP_singlerow_2_withoutCR_allchrom.csv", covfilename = "Ovens_Stocking_and Monitoring_COV_File.csv")


#Order 
##Rep av
##monomorphic
##callrate loci
##secondaries
##callrate ind
##HWE 
##hamming distance

#Recoding ind names based on colony run and known relationships
?gl.make.recode.ind

#Don't rewrite over the correct one a second time
#gl.make.recode.ind(glbase, outfile="new_ind_assignment2s.csv",outpath=getwd())
# after written had to manually change based on sorting and prev colony runs
indNames(glbase)
nInd(glbase)
?gl.recode.ind
glnew3 <- gl.recode.ind(glbase, ind.recode="new_ind_assignments2.csv")
getwd()

# Inds which were not confidently assignmed in large run
glreduced <- gl.drop.ind(glnew3, ind.list=c("REMOVE","MP-OV012U","MP-OV018U","MP-OV040U","MP-OV048U","MP-OV050U","MP-OV056U","MP-OV062U","MP-OV066U","MP-OV067U","MP-OV076U","MP-OV089U","MP-OV096U","MP-OV108U","MP-OV116U","MP-OV123U","O1OM101XF3X"
                                              , "O2OM101XF3X", "O1OM1XF2X", "O1OM2XF4X", "O2OM2XF4X","O10M20XF7X", "O20M20XF7X", "O1OM15XF15X", "O1OM11XF13X", "O1OM18XF18","O1OM3F5X","O2OM15XF15X", "O1OM101XF2X", "O1OM5XF2X" ))

# (Note, test for sex linkage before filtering)
indNames(glreduced)
nLoc(glbase)
nInd(glbase)
gl_2_repavg <- gl.filter.reproducibility(glreduced, t=0.99)
nLoc(gl_2_repavg)
gl_3_monomorphs <- gl.filter.monomorphs(gl_2_repavg, v=1)
nLoc(gl_3_monomorphs)
par(mar = c(2, 2, 2, 2))
gl_4_callrate <- gl.filter.callrate(gl_3_monomorphs, method = "loc", threshold = 0.90)
nLoc(gl_4_callrate)
gl_5_secondaries <- gl.filter.secondaries(gl_4_callrate)
nLoc(gl_5_secondaries)
nInd(gl_5_secondaries)
gl_6_callrate_ind <- gl.filter.callrate(gl_5_secondaries, method = "ind", threshold = 0.90)
nLoc(gl_6_callrate_ind)
nInd(gl_6_callrate_ind)
gl_final <- gl_6_callrate_ind

# change to  0.15 or 0.02 for the two MAFs
glfilterdmaf <- gl.filter.maf(gl_final, threshold = 0.15)
#glfilterdmaf <- gl.filter.maf(gl_final, threshold = 0.02)

filtered_freqs <- glMean(glfilterdmaf)
hist(filtered_freqs)
mean(filtered_freqs)
max(filtered_freqs)
sd(filtered_freqs)
gl.report.maf(glfilterdmaf)
nLoc(glfilterdmaf)
##### gl_7_HWE <- gl.filter.hwe(p=0.05, gl_6_callrate_ind) (Don't do HWE first)



tiff("MP MAF graphs off_adults_015MAF.tif",units='cm', width =30,height=30, res=400)
gl.report.maf(glfilterdmaf)
dev.off()


gl_6_callrate_ind <- gl.filter.callrate(glfilterdmaf, method = "ind", threshold = 0.85)
nLoc(gl_6_callrate_ind)
gl_7_monomorphs <- gl.filter.monomorphs(gl_6_callrate_ind, v=1)
gl_final <- gl_7_monomorphs

nLoc(gl_final)
nInd(gl_7_monomorphs)


indNames(gl_final)
gl.report.maf(glfinal)
glfinal <- gi2gl(genidallindivs)
mean(glfinal@other$loc.metrics$maf)
sd(gl_final@other$loc.metrics$maf)


# remove individuals with uncertain parentage

# This is moved to the top now it's been tested. THis has been tested and so the ind dropping should take place at the start of the filtering steps

#glreduced <- gl.drop.ind(gl_final, ind.list=c("REMOVE","MP-OV012U","MP-OV018U","MP-OV040U","MP-OV048U","MP-OV050U","MP-OV056U","MP-OV062U","MP-OV066U","MP-OV067U","MP-OV076U","MP-OV089U","MP-OV096U","MP-OV108U","MP-OV116U","MP-OV123U","O1OM101XF3X"
#                                              , "O2OM101XF3X", "O1OM1XF2X", "O1OM2XF4X", "O2OM2XF4X","O10M20XF7X", "O20M20XF7X", "O1OM15XF15X", "O1OM11XF13X", "O1OM18XF18","O1OM3F5X","O2OM15XF15X", "O1OM101XF2X", "O1OM5XF2X" ))

#Test, remove all individuals with no full sibling relationships
#glreduced <- gl.drop.ind(glreduced, ind.list=c("O1OM10XF4X", "O1OM13XF14X", "O1OM14XF15X", "O1OM15F16", "O1OM16XF9", "O1OM16XF9", "O1OM17XF15X", "O1OM17XF16X", "O1OM17XF5X", "O1OM17XF7X", "O1OM23XF18X", "O1OM23XF19X", "O1OM23XF9X", "O1OM24XF4X", "O1OM25XF21X", "O1OM2F2","O1OM3F100X", "O1OM5XF8X", "O1OM7XF11X"))





# confirm secondaries removed correctly (function already loaded on)

check_dartR_locnames(gl_final)



#number of loci after all restrictions
#number of individuals after all restrictions
nLoc(gl_final)
nInd(gl_final)


#sex linkage analysis
gl.sexlinkage(gl_final)


# create gi file for colony on HPC


rm(gl_2_repavg, gl_3_monomorphs, gl_4_callrate, gl_5_secondaries,gl_6_callrate_ind, gl_7_monomorphs,glbase,glfilterdmaf,glreduced,glnew3)

gi_final <- gl2gi(gl_final)
save.image("filtered_recoded_MP_genlight_genind_15MAF_99repro.RData")



# randomise seed

set.seed(12654421)

gl.drop.pop2 <- function (x, pop.list, recalc = FALSE, mono.rm = FALSE, v = 2) {
  if (class(x) != "genlight") {
    cat("Fatal Error: genlight object required for gl.drop.pop.r!\n")
    stop("Execution terminated\n")
  }
  if (v > 0) {
    cat("Starting gl.drop.pop: Deleting selected populations\n")
  }
  if (v > 1) {
    cat("  Deleting populations", pop.list, "\n")
  }
  x <- x[!x$pop %in% pop.list]
  if (mono.rm) {
    x <- gl.filter.monomorphs(x, v = 0)
  }
  if (recalc) {
    gl.recalc.metrics(x, v = v)
  }
  if (v > 2) {
    cat("Summary of recoded dataset\n")
    cat(paste("  No. of loci:", nLoc(x), "\n"))
    cat(paste("  No. of individuals:", nInd(x), "\n"))
    cat(paste("  No. of populations: ", length(levels(factor(pop(x)))), 
              "\n"))
  }
  if (v > 1) {
    if (!recalc) {
      cat("Note: Locus metrics not recalculated\n")
    }
    if (!mono.rm) {
      cat("Note: Resultant monomorphic loci not deleted\n")
    }
  }
  if (v > 0) {
    cat("Completed gl.drop.pop\n\n")
  }
  return <- x
}

gi2gl2<- function(gi, parallel) {
  
  locna <- gi@loc.n.all
  ccc <- 1
  for (i in 2:length(locna)) {
    if (locna[i - 1] == 1) 
      ccc[i] <- ccc[i - 1] + 1
    else ccc[i] <- ccc[i - 1] + 2
  }
  gl <- new("genlight", gi@tab[, ccc], pop = pop(gi), 
            other = gi@other, ploidy = 2, loc.names = locNames(gi), 
            ind.names = indNames(gi), parallel = parallel)
  return(gl)
}

#a = object, b = loc proportion, c = ind proportion 
subset_data <- function (a, b, c) {
  
  
  locindex <- sample(nLoc(a), (nLoc(a)*b), replace = F)
  indindex <- sample(nInd(a), (nInd(a)*c), replace = F)
  gl <- a[indindex, locindex]
  gl@other$ind.metrics <- a@other$ind.metrics[indindex,]
  gl@other$loc.metrics <- a@other$loc.metrics[locindex,]
  x <- as.matrix(gl[, ])
  for (i in 1:nrow(x)) {
    for (ii in 1:ncol(x)) {
      inp <- x[i, ii]
      if (!is.na(inp)) {
        if (inp == 0) 
          x[i, ii] <- "A/A"
        else if (inp == 1) 
          x[i, ii] <- "A/B"
        else if (inp == 2) 
          x[i, ii] <- "B/B"
      }
    }
    
  }
  
  gen <- df2genind(x[, ], sep = "/", ncode = 1, ind.names = gl@ind.names, 
                   pop = gl@pop, ploidy = 2)
  gen@other <- gl@other
  
  gen
}

genidallindivs <- gi_final  

vector1 <- rep('A', times = nInd(genidallindivs)-1, each = 1)
vector2 <- rep('B', times = nInd(genidallindivs)+1, each = 1)         
vector3 <- c(vector1,vector2)
vector4 <- vector3[1:nInd(genidallindivs)]
vector4 <- as.factor(vector4)
genidallindivs@pop 
genidallindivs@pop <- as.factor(vector4)
genidallindivs@pop<-  as.factor(genidallindivs@pop)
rm(vector1,vector2,vector3,vector4)

genlightallindivs <- gi2gl2(genidallindivs, FALSE) 

# Need to have a pop column which is adults vs offspring.
# Filter to offspring section
for (i in c(1:length(genlightallindivs@ind.names))) {
  
  if  (grepl('O', genlightallindivs@ind.names[i])) {
    genlightallindivs@pop[i] <- 'A'
  }
  else {
    genlightallindivs@pop[i] <- 'B'
  }
}



genlightoffspring <- gl.drop.pop2(genlightallindivs, pop.list = c('B'))


genidoffspring <- subset_data(genlightoffspring, 1, 1)
#genidadults <- subset_data(genlightadults, 1, 0.5)

# end filter to offspring section

genlightadults 

genlightoffspring



nLoc(gi_final)

# Loci numbers to filter down to for each MAF

# for MAF 2%
# 2396, 2000, 1600 1200, 1000, 800, 600, 400, 300, 200, 100

2396*0.8348 #(2000 markers)

2396*0.6678 #(1600 markers)

2396*0.5009 #(1200 markers)

2396*0.41739 #(1000 markers)

2396*0.3339 #(800 markers)

2396*0.2505 #(600 markers)

2396*0.167 #(400 markers)

2396*0.1253 #(300 markers)

2396*0.0836 #(200 markers)

2396*0.0418 #(100 markers)


# for MAF 15%
# 693 loci = 100%
 693*0.8659 #= 600 markers
 693*0.722 #= 500 markers
 693*0.5776 #= 400 markers
 693*0.433 #= 300 markers
 693*0.289 #= 200 markers
 693*0.1445 #= 100 markers
 693*0.0724 #= 50 markers

693*0.0724

472 =400*x 
x =472/400 
472/400
400*1.18
(1.18/472)*100
0.106*472





# all inds 
#genidreduced <- subset_data(genlightallindivs,0.167, 1)
# offspring only one
genidreduced <- subset_data(genlightoffspring,0.0724, 1)

nLoc(genidreduced)
nInd(genidreduced)
#rm(genidallindivs,genlightadults,genlightoffspring,genidadults,genidoffspring)
getwd()



tidydata <- tidy_genind(
  genidreduced,
  keep.allele.names = TRUE,
  tidy = TRUE,
  gds = FALSE,
  write = FALSE,
  verbose = FALSE
)

# save reduced genid data for later incorporation into related
save.image("D:/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/genidreduced100%ind7%loc15MAF_99repoff.RData")
#save.image("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/genidreduced100%ind100%loc02MAF_99rep.RData")


# create colony file 
write_colony(tidydata)





