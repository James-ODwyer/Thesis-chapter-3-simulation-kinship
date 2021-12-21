


#Supplementary material 2-F R code to analyse half-sibling relationships for Macquarie perch data
# What makes for reliable kinship assignments; how sample design, genetic markers and programs used shape assignment accuracies.




#Half sib completed code 
#change based on which run is being analysed

setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/offspring only")
#half_sibtable <-read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/allindivs/36.5%loc/MP100ind36loc.HalfSibDyad",sep=",")
#half_sibtable <-read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/offspring only/6.09%loci/MP100ind6locoff.HalfSibDyad",sep=",")
#half_sibtable <-read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/halfpar/rerun/6.09%loc/MP100ind6loc50par.HalfSibDyad",sep=",")
#half_sibtable <-read.csv("MP100ind18locofffsonly.HalfSibDyad")
#half_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/10%ind/1%loc/simulation10ind1locpwlong.HalfSibDyad",sep=",")
#half_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/2%ind/1%loc/simulation2ind1locoffpwlong.HalfSibDyad",sep=",")
#half_sibtable <-read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/15MAF/11%loc/MP100ind11loc15MAF.HalfSibDyad",sep=",")
#half_sibtable <-read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/offspring15MAF/11%loc/MP100ind11loc15MAFoff.HalfSibDyad",sep=",")
half_sibtable <-read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/MP99/15MAF/offspring only/7loc/MP100ind7loc15MAF99repoff.HalfSibDyad",sep=",")


for (i in c(1:nrow(half_sibtable))) {
  
  sub1_1_1 <- sub(".*M", "M", half_sibtable[i,1])
  sub1_1_2 <- sub("F.*", "", sub1_1_1)

  sub1_2_1 <- sub(".*F", "F",half_sibtable[i,1])
  sub1_2_2 <- sub("T.*","", sub1_2_1)

  sub2_1_1 <- sub(".*M", "M", half_sibtable[i,2])
  sub2_1_2 <- sub("F.*", "", sub2_1_1)

  sub2_2_1 <- sub(".*F", "F",half_sibtable[i,2])
  sub2_2_2 <- sub("T.*","", sub2_2_1)
  
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
summary_table_halfsibs[1,7] <-(summary_table_halfsibs[1,2] / sum(summary_table_halfsibs[1,2]+summary_table_halfsibs[1,4]+summary_table_halfsibs[1,6])*100)
summary_table_halfsibs[1,8] <- ((summary_table_halfsibs[1,2]/summary_table_halfsibs[1,1])*100)
sum(half_sibtable$Probability[1:nrow(half_sibtable)] >=0.95)


# Un comment this part and run once initially when starting
#multiple_analysis_summary_table_halfsibs <- matrix(data = NA, nrow = 14, ncol = 8) 
#colnames(multiple_analysis_summary_table_halfsibs) =  c("correctly identified half siblings",
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
multiple_analysis_summary_table_halfsibs[8,] <- summary_table_halfsibs[1,]
#setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/results")
write.csv(multiple_analysis_summary_table_halfsibs, 'multiple_analysis_summary_table_halfsibs_MP_100indivs15MAF99repoff.csv')


