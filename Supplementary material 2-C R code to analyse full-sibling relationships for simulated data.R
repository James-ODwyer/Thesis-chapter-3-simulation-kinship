# Supplementary material 2-C R code to analyse full-sibling relationships for simulated data
# What makes for reliable kinship assignments; how sample design, genetic markers and programs used shape assignment accuracies.


#Full sib completed code 
#change based on which run is being analysed
# I have run this analysis on two different computers through a dropbox so the file pathways are slightly different
setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults")
setwd("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults")
#change name of file based on which run is being analysed


full_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/25%ind/4%loc/simulation25ind4locfpls.FullSibDyad" ,sep=",")
#full_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/1%ind/10%loc/simulation1ind10locpwlong.FullSibDyad" ,sep=",")
full_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/1%ind/10%loc/simulation1ind10locoffpwlong.FullSibDyad",sep=",")
full_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/25%ind/4%loc/simulation25ind4locofffpls.FullSibDyad",sep=",")

#full_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/variedparentsproportion/25%parents/1%off/2%loc/simulation1ind5loc25parfpls.FullsibDyad",sep=",")
#full_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/1%ind/1%loc/simulation1ind1locoffpwlong.FullSibDyad",sep=",")
#full_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/25%ind/005%loc/simulation25ind005locpwlong_2.FullSibDyad" ,sep=",")
#full_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/50%ind/1%loc/redo2/simulation50ind1locofffpls.FullSibDyad",sep=",")


#For the partial parents data
#full_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/variedparentsproportion/redo/25%parents/1%off/005%loc/simulation1ind005loc25parfpls.FullsibDyad",sep=",")
#full_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/variedparentsproportion/25%parents/1%off/005%loc/simulation1ind005loc25parfpls.FullsibDyad" ,sep=",")
#full_sibtable <- read.csv("simulation5ind5loc50parfpls.FullSibDyad" ,sep=",")



#create starting new columns to 
full_sibtable$M1[1] <- NA
full_sibtable$F1[1] <- NA
full_sibtable$M2[1] <- NA
full_sibtable$F2[1] <- NA
full_sibtable$off1fam[1] <- NA
full_sibtable$off2fam[1] <- NA




for (j in c(1:nrow(full_sibtable))) { 
  sub1_1 <- sub(".*G", "G",full_sibtable[j,1])
  
  sub2_1 <- sub(".*G", "G",full_sibtable[j,2])

  full_sibtable$off1fam[j] <- sub1_1
  full_sibtable$off2fam[j] <- sub2_1
  
}
#half sib dyads, change the starting row based on how many full sib dyads exist
for (i in c(1:nrow(full_sibtable))) {
  sub1_1 <- sub(".*M", "M",full_sibtable[i,1])
  sub1_1_2 <- sub("F.*", "",sub1_1)
  sub1_2 <- sub (".*F", "F", sub1_1)
  sub1_2_1 <- sub ("G.*", "", sub1_2)
  #sub1_2_2 <- substr (sub1_2_1, 1,4)
  
  
  sub2_1 <- sub(".*M", "M",full_sibtable[i,2])
  sub2_1_2 <- sub("F.*", "",sub2_1)
  sub2_2 <- sub (".*F", "F", sub2_1)
  sub2_2_1 <- sub ("G.*", "", sub2_2)
  #sub2_2_2 <- substr (sub2_2_1, 1,4)
  
  
  full_sibtable$M1[i] <- sub1_1_2
  full_sibtable$F1[i] <- sub1_2_1
  full_sibtable$M2[i] <- sub2_1_2
  full_sibtable$F2[i] <- sub2_2_1
  
}


#create table to store data 
summary_table_fullsibs <- matrix(data = NA, nrow = 8, ncol = 10) 
colnames(summary_table_fullsibs) = c("fullsibs correct","mother correct", "father correct", "fullsib correct probability over 95%", "fullsibs neither correct", "fullsib incorrect over 95%","total dyad pairs", "halfsib relationships listed as full","percentage correctly identified over incorrectly for over 95%",
                                     "percent correctly identified over 95% over total percent correctly identified")
#rownames(summary_table_fullsibs) = c("pairwise full data", "pairwise all ind 75% loc", "pairwise all ind 50% loc", "pairwise all ind 25% loc")
#results
#input full sib and half sib starting points, p = where fullsib stops q= where half sib starts (q=p+1) control find 'off' in fullsibtable
# r = area where intermixxing starts 
#p = 47900
#q = 47901
#r =43608

summary_table_fullsibs[1,1] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] == full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] == full_sibtable$M2[1:nrow(full_sibtable)] )
summary_table_fullsibs[1,2] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] == full_sibtable$F2[1:nrow(full_sibtable)])
summary_table_fullsibs[1,3] <-sum(full_sibtable$M1[1:nrow(full_sibtable)] == full_sibtable$M2[1:nrow(full_sibtable)])
summary_table_fullsibs[1,4] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] == full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] == full_sibtable$M2[1:nrow(full_sibtable)] & full_sibtable$Probability[1:nrow(full_sibtable)] >=0.95 )
summary_table_fullsibs[1,5] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] != full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] != full_sibtable$M2[1:nrow(full_sibtable)] )
summary_table_fullsibs[1,6] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] != full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] != full_sibtable$M2[1:nrow(full_sibtable)] & full_sibtable$Probability[1:nrow(full_sibtable)] >=0.95  )
#summary_table_fullsibs[1,7] <-sum(full_sibtable$off1fam[1:p] ==full_sibtable$off2fam[1:p])
#summary_table_fullsibs[1,8] <-sum(full_sibtable$off1fam[1:p] ==full_sibtable$off2fam[1:p] & full_sibtable$Probability[1:p] >=0.95)
#summary_table_fullsibs[1,9] <-sum(full_sibtable$off1fam[1:p] !=full_sibtable$off2fam[1:p])
#summary_table_fullsibs[1,10] <-sum(full_sibtable$off1fam[1:p] !=full_sibtable$off2fam[1:p] & full_sibtable$Probability[1:p] >=0.95)
summary_table_fullsibs[1,7] <- nrow(full_sibtable)
summary_table_fullsibs[1,8] <- sum((full_sibtable$M1[1:nrow(full_sibtable)] == full_sibtable$M2[1:nrow(full_sibtable)] & full_sibtable$F1[1:nrow(full_sibtable)] != full_sibtable$F2[1:nrow(full_sibtable)]) | (full_sibtable$F1[1:nrow(full_sibtable)] == full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] != full_sibtable$M2[1:nrow(full_sibtable)]))
summary_table_fullsibs[1,9] <- ((summary_table_fullsibs[1,4])/(summary_table_fullsibs[1,4]+summary_table_fullsibs[1,6])*100)
summary_table_fullsibs[1,10] <- ((summary_table_fullsibs[1,4]/summary_table_fullsibs[1,1])*100)
# 12 = M1 = M2 & F1 != F2 | F1 = F2 & M1 != M2

#Placeholder for R (halfsibs listed as full)

#halfsiblings misidentified as full siblings |

# Un comment this part and run once initially when starting
#multiple_analysis_summary_table_fullsibs <- matrix(data = NA, nrow = 11, ncol = 10) 
#colnames(multiple_analysis_summary_table_fullsibs) = c("fullsibs correct","mother correct", "father correct", "fullsib correct probability over 95%", "fullsibs neither correct", "fullsib incorrect over 95%","total dyad pairs", "halfsib relationships listed as full","percentage correctly identified over incorrectly for over 95%",
#                                    "percent correctly identified over 95% over total percent correctly identified")

# Current code is manual per read in (although this can be automated fairly easily, I just didn't take the time to do it)
# So each loading of colony file then run the main body, then put it into the correct row for the multiple analysis summary table, increase row
# num in the multiple analysis summary table, and load in next colony results
multiple_analysis_summary_table_fullsibs[11,] <- summary_table_fullsibs[1,]


# Save final results
#setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/variedparentsproportion/25%parents")
#setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/results")
#setwd("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults")
#setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/variedparentsproportion/redo")
#write.csv(multiple_analysis_summary_table_fullsibs, 'multiple_analysis_summary_table_fullsibs_1indivs_25par_fpls.csv')
write.csv(multiple_analysis_summary_table_fullsibs, 'multiple_analysis_summary_table_fullsibs_1indivsoff_fpls.csv')

