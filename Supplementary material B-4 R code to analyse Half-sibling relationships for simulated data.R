# Supplementary material 2-D R code to analyse Half-sibling relationships for simulated data
# What makes for reliable kinship assignments; how sample design, genetic markers and programs used shape assignment accuracies.


# I have run this analysis on two different computers through a dropbox so the file pathways are slightly different

#setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults")

setwd("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/variedparentsproportion/100%parents")
read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/variedparentsproportion/100%parents/10%off/5%loc/simulation10ind5locallparfpls.HalfsibDyad" ,sep=",")
#half_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/10%ind/1%loc/simulation10ind1locpwlong.HalfSibDyad",sep=",")
#half_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/2%ind/1%loc/simulation2ind1locoffpwlong.HalfSibDyad",sep=",")

#half_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/25%ind/4%loc/simulation25ind4locfpls.HalfSibDyad" ,sep=",")
#half_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/1%ind/10%loc/simulation1ind10locpwlong.HalfSibDyad" ,sep=",")
half_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/1%ind/10%loc/simulation100ind10locoffpwlong.HalfSibDyad",sep=",")
half_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/25%ind/4%loc/simulation25ind4locofffpls.HalfSibDyad",sep=",")



#setwd("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/100%ind")
#half_sibtable <- read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/1%ind/005%loc/simulation1ind005locfpls.HalfSibDyad" ,sep=",")


for (i in c(1:nrow(half_sibtable))) {
  
  sub1_1_1 <- sub(".*M", "M", half_sibtable[i,1])
  sub1_1_2 <- sub("F.*", "", sub1_1_1)

  sub1_2_1 <- sub(".*F", "F",half_sibtable[i,1])
  sub1_2_2 <- sub("G.*","", sub1_2_1)

  sub2_1_1 <- sub(".*M", "M", half_sibtable[i,2])
  sub2_1_2 <- sub("F.*", "", sub2_1_1)

  sub2_2_1 <- sub(".*F", "F",half_sibtable[i,2])
  sub2_2_2 <- sub("G.*","", sub2_2_1)
  
half_sibtable$M1[i] <- sub1_1_2
half_sibtable$F1[i] <- sub1_2_2
half_sibtable$M2[i] <- sub2_1_2
half_sibtable$F2[i] <- sub2_2_2

}


summary_table_halfsibs <- matrix(data = NA, nrow = 10, ncol = 8)
colnames(summary_table_halfsibs) = c("correctly identified half siblings",
                                     "correctly identified half siblings over 95%",
                                     "incorrectly indentified half siblings", 
                                     "incorrectly indentified half siblings over 95%",
                                     "fullsiblings misidentified as halfsiblings", 
                                     "fullsibings misidentified as halfsiblings over 95%",
                                     "percentage correctly identified over incorrectly for over 95%",
                                     "percent correctly identified over 95% over total percent correctly identified")

summary_table_halfsibs[1,1] <-(sum(half_sibtable$M1[1:nrow(half_sibtable)] == half_sibtable$M2[1:nrow(half_sibtable)] |  half_sibtable$F1[1:nrow(half_sibtable)] == half_sibtable$F2[1:nrow(half_sibtable)]) - sum(half_sibtable$M1[1:nrow(half_sibtable)] == half_sibtable$M2[1:nrow(half_sibtable)] &  half_sibtable$F1[1:nrow(half_sibtable)] == half_sibtable$F2[1:nrow(half_sibtable)]))
summary_table_halfsibs[1,2] <- sum((half_sibtable$M1[1:nrow(half_sibtable)] == half_sibtable$M2[1:nrow(half_sibtable)] & half_sibtable$Probability[1:nrow(half_sibtable)] >=0.95) | (half_sibtable$F1[1:nrow(half_sibtable)] == half_sibtable$F2[1:nrow(half_sibtable)] & half_sibtable$Probability[1:nrow(half_sibtable)] >=0.95)) -sum(half_sibtable$M1[1:nrow(half_sibtable)] == half_sibtable$M2[1:nrow(half_sibtable)] &  half_sibtable$F1[1:nrow(half_sibtable)] == half_sibtable$F2[1:nrow(half_sibtable)] & half_sibtable$Probability[1:nrow(half_sibtable)] >=0.95)
summary_table_halfsibs[1,3] <-(sum(half_sibtable$M1[1:nrow(half_sibtable)] != half_sibtable$M2[1:nrow(half_sibtable)] &  half_sibtable$F1[1:nrow(half_sibtable)] != half_sibtable$F2[1:nrow(half_sibtable)]) + sum(half_sibtable$M1[1:nrow(half_sibtable)] == half_sibtable$M2[1:nrow(half_sibtable)] &  half_sibtable$F1[1:nrow(half_sibtable)] == half_sibtable$F2[1:nrow(half_sibtable)]))
summary_table_halfsibs[1,4] <-(sum(half_sibtable$M1[1:nrow(half_sibtable)] != half_sibtable$M2[1:nrow(half_sibtable)] &  half_sibtable$F1[1:nrow(half_sibtable)] != half_sibtable$F2[1:nrow(half_sibtable)] & half_sibtable$Probability[1:nrow(half_sibtable)] >=0.95) + sum(half_sibtable$M1[1:nrow(half_sibtable)] == half_sibtable$M2[1:nrow(half_sibtable)] &  half_sibtable$F1[1:nrow(half_sibtable)] == half_sibtable$F2[1:nrow(half_sibtable)] & half_sibtable$Probability[1:nrow(half_sibtable)] >=0.95))
summary_table_halfsibs[1,5] <-sum(half_sibtable$M1[1:nrow(half_sibtable)] == half_sibtable$M2[1:nrow(half_sibtable)] &  half_sibtable$F1[1:nrow(half_sibtable)] == half_sibtable$F2[1:nrow(half_sibtable)])
summary_table_halfsibs[1,6] <-sum(half_sibtable$M1[1:nrow(half_sibtable)] == half_sibtable$M2[1:nrow(half_sibtable)] &  half_sibtable$F1[1:nrow(half_sibtable)] == half_sibtable$F2[1:nrow(half_sibtable)] & half_sibtable$Probability[1:nrow(half_sibtable)] >=0.95)
summary_table_halfsibs[1,7] <-(summary_table_halfsibs[1,2] / sum(summary_table_halfsibs[1,2]+summary_table_halfsibs[1,4])*100)
summary_table_halfsibs[1,8] <- ((summary_table_halfsibs[1,2]/summary_table_halfsibs[1,1])*100)
sum(half_sibtable$Probability[1:nrow(half_sibtable)] >=0.95)

#Check the output of summary table 1,7 when redoing it, it seems like it should be double counting the full siblings
# But some of the higher ups are 100% accurate so therefore it cant't be? Just double check...
# Will need to reanalyse anyway with new data so all good. 
# read above (30/11/2020)
# Full sibs checked and no errors in it

# double checked, it is double counting. can lead to a few % points off unfortunately. 
# Fix with rerun, tested 1/12/2020
# old code summary_table_halfsibs[1,7] <-(summary_table_halfsibs[1,2] / sum(summary_table_halfsibs[1,2]+summary_table_halfsibs[1,4]+summary_table_halfsibs[1,6])*100)
# new code summary_table_halfsibs[1,7] <-(summary_table_halfsibs[1,2] / sum(summary_table_halfsibs[1,2]+summary_table_halfsibs[1,4])*100)


# Un comment this part and run once initially when starting
#multiple_analysis_summary_table_halfsibs <- matrix(data = NA, nrow = 11, ncol = 8) 
colnames(multiple_analysis_summary_table_halfsibs) =  c("correctly identified half siblings",
                                                        "correctly identified half siblings over 95%",
                                                        "incorrectly indentified half siblings", 
                                                        "incorrectly indentified half siblings over 95%",
                                                        "fullsiblings misidentified as halfsiblings", 
                                                        "fullsibings misidentified as halfsiblings over 95%",
                                                        "percentage correctly identified over incorrectly for over 95%",
                                                        "percent correctly identified over 95% over total percent correctly identified")


# Current code is manual per read in (although this can be automated fairly easily, I just didn't take the time to do it)
# So each loading of colony file then run the main body, then put it into the correct row for the multiple analysis summary table, increase row
# num in the multiple analysis summary table, and load in next colony results

multiple_analysis_summary_table_halfsibs[11,] <- summary_table_halfsibs[1,]
#setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/results")
setwd("D:/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults")
write.csv(multiple_analysis_summary_table_halfsibs, 'multiple_analysis_summary_table_halfsibs_1indivs_offs_fpls.csv')
