
#Supplementary material 2-E R code to analyse full-sibling relationships for Macquarie perch data
# What makes for reliable kinship assignments; how sample design, genetic markers and programs used shape assignment accuracies.


#Full sib completed code 
#change based on which run is being analysed


#change name of file based on which run is being analysed
#setwd("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/offspring only")
setwd("D:/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/halfpar/rerun")
#full_sibtable <-read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/offspring only/6.09%loci/MP100ind6locoff.FullSibDyad",sep=",")
#full_sibtable <-read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/halfpar/6.09%loc/MP100ind6loc50par.FullSibDyad",sep=",")
#full_sibtable <-read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/15MAF/85%loc/MP100ind85loc15MAF.FullSibDyad",sep=",")

full_sibtable <-read.csv("D:/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/MP99/15MAF/offspring only/7loc/MP100ind7loc15MAF99repoff.FullSibDyad",sep=",")

#full_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/1%ind/1%loc/simulation1ind1locpwlong.FullSibDyad",sep=",")
#full_sibtable <- read.csv("C:/Users/18088076/Dropbox/PhD/Chapters/Power analysis/Simulation/Colony runs/Simulated dataset/addedotheradults/offspring only/1%ind/1%loc/simulation1ind1locoffpwlong.FullSibDyad",sep=",")



#create starting new columns to 
full_sibtable$M1[1] <- NA
full_sibtable$F1[1] <- NA
full_sibtable$M2[1] <- NA
full_sibtable$F2[1] <- NA
full_sibtable$off1fam[1] <- NA
full_sibtable$off2fam[1] <- NA




for (j in c(1:nrow(full_sibtable))) { 
  sub1_1 <- sub(".*T", "T",full_sibtable[j,1])
  
  sub2_1 <- sub(".*T", "T",full_sibtable[j,2])

  full_sibtable$off1fam[j] <- sub1_1
  full_sibtable$off2fam[j] <- sub2_1
  
}
#half sib dyads, change the starting row based on how many full sib dyads exist
for (i in c(1:nrow(full_sibtable))) {
  sub1_1 <- sub(".*M", "M",full_sibtable[i,1])
  sub1_1_2 <- sub("F.*", "",sub1_1)
  sub1_2 <- sub (".*F", "F", sub1_1)
  sub1_2_1 <- sub ("T.*", "", sub1_2)
  #sub1_2_2 <- substr (sub1_2_1, 1,4)
  
  
  sub2_1 <- sub(".*M", "M",full_sibtable[i,2])
  sub2_1_2 <- sub("F.*", "",sub2_1)
  sub2_2 <- sub (".*F", "F", sub2_1)
  sub2_2_1 <- sub ("T.*", "", sub2_2)
  #sub2_2_2 <- substr (sub2_2_1, 1,4)
  
  
  full_sibtable$M1[i] <- sub1_1_2
  full_sibtable$F1[i] <- sub1_2_1
  full_sibtable$M2[i] <- sub2_1_2
  full_sibtable$F2[i] <- sub2_2_1
  
}


#create table to store data 
summary_table_fullsibs <- matrix(data = NA, nrow = 8, ncol = 14) 
colnames(summary_table_fullsibs) = c("fullsibs correct","mother correct", "father correct", "fullsib correct probability over 95%", "fullsibs neither correct", "fullsib incorrect over 95%","total dyad pairs", "halfsib relationships listed as full","percentage correctly identified over incorrectly for over 95%",
                                     "percent correctly identified over 95% over total percent correctly identified", "number of known relationships correctly identified", "number of known relationships incorrectly identified","number of known relationships correctly identified over 95% confidence", "number of known relatinoships incorrectly identified over 95% confidence")
#rownames(summary_table_fullsibs) = c("pairwise full data", "pairwise all ind 75% loc", "pairwise all ind 50% loc", "pairwise all ind 25% loc")


summary_table_fullsibs[1,1] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] == full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] == full_sibtable$M2[1:nrow(full_sibtable)] )
summary_table_fullsibs[1,2] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] == full_sibtable$F2[1:nrow(full_sibtable)])
summary_table_fullsibs[1,3] <-sum(full_sibtable$M1[1:nrow(full_sibtable)] == full_sibtable$M2[1:nrow(full_sibtable)])
summary_table_fullsibs[1,4] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] == full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] == full_sibtable$M2[1:nrow(full_sibtable)] & full_sibtable$Probability[1:nrow(full_sibtable)] >=0.95 )
summary_table_fullsibs[1,5] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] != full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] != full_sibtable$M2[1:nrow(full_sibtable)] )
summary_table_fullsibs[1,6] <-sum(full_sibtable$F1[1:nrow(full_sibtable)] != full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] != full_sibtable$M2[1:nrow(full_sibtable)] & full_sibtable$Probability[1:nrow(full_sibtable)] >=0.95  )
summary_table_fullsibs[1,7] <- nrow(full_sibtable)
summary_table_fullsibs[1,8] <- sum((full_sibtable$M1[1:nrow(full_sibtable)] == full_sibtable$M2[1:nrow(full_sibtable)] & full_sibtable$F1[1:nrow(full_sibtable)] != full_sibtable$F2[1:nrow(full_sibtable)]) | (full_sibtable$F1[1:nrow(full_sibtable)] == full_sibtable$F2[1:nrow(full_sibtable)] & full_sibtable$M1[1:nrow(full_sibtable)] != full_sibtable$M2[1:nrow(full_sibtable)]))
summary_table_fullsibs[1,9] <- ((summary_table_fullsibs[1,4])/(summary_table_fullsibs[1,4]+summary_table_fullsibs[1,6])*100)
summary_table_fullsibs[1,10] <- ((summary_table_fullsibs[1,4]/summary_table_fullsibs[1,1])*100)
summary_table_fullsibs[1,11] <- sum(full_sibtable$off1fam[1:nrow(full_sibtable)] == full_sibtable$off2fam[1:nrow(full_sibtable)] & grepl('T',full_sibtable$off1fam[1:nrow(full_sibtable)]) & grepl('T',full_sibtable$off2fam[1:nrow(full_sibtable)])  )            
summary_table_fullsibs[1,12] <- sum(full_sibtable$off1fam[1:nrow(full_sibtable)] != full_sibtable$off2fam[1:nrow(full_sibtable)] & grepl('T',full_sibtable$off1fam[1:nrow(full_sibtable)]) & grepl('T',full_sibtable$off2fam[1:nrow(full_sibtable)])  )            
summary_table_fullsibs[1,13] <-sum(full_sibtable$off1fam[1:nrow(full_sibtable)] == full_sibtable$off2fam[1:nrow(full_sibtable)] & grepl('T',full_sibtable$off1fam[1:nrow(full_sibtable)]) & grepl('T',full_sibtable$off2fam[1:nrow(full_sibtable)]) &  full_sibtable$Probability[1:nrow(full_sibtable)] >=0.95)            
summary_table_fullsibs[1,14] <- sum(full_sibtable$off1fam[1:nrow(full_sibtable)] != full_sibtable$off2fam[1:nrow(full_sibtable)] & grepl('T',full_sibtable$off1fam[1:nrow(full_sibtable)]) & grepl('T',full_sibtable$off2fam[1:nrow(full_sibtable)]) &  full_sibtable$Probability[1:nrow(full_sibtable)] >=0.95  )            


# Un comment this part and run once initially when starting
#multiple_analysis_summary_table_fullsibs <- matrix(data = NA, nrow = 13, ncol = 14) 
#colnames(multiple_analysis_summary_table_fullsibs) = c("fullsibs correct","mother correct", "father correct", "fullsib correct probability over 95%", "fullsibs neither correct", "fullsib incorrect over 95%","total dyad pairs", "halfsib relationships listed as full","percentage correctly identified over incorrectly for over 95%",
#                                    "percent correctly identified over 95% over total percent correctly identified", "number of known relationships correctly identified", "number of known relationships incorrectly identified","number of known relationships correctly identified over 95% confidence", "number of known relatinoships incorrectly identified over 95% confidence")



# Current code is manual per read in (although this can be automated fairly easily, I just didn't take the time to do it)
# So each loading of colony file then run the main body, then put it into the correct row for the multiple analysis summary table, increase row
# num in the multiple analysis summary table, and load in next colony results

multiple_analysis_summary_table_fullsibs[8,] <- summary_table_fullsibs[1,]
setwd("D:/Dropbox/PhD/Chapters/Power analysis/Power analysis MP/MP dataset/colonyruns/MP99")
write.csv(multiple_analysis_summary_table_fullsibs, 'multiple_analysis_summary_table_fullsibs_MP_100indivs15MAF99repoff.csv')


