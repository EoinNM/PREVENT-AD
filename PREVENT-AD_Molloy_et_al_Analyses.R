##########################################################################################################################

#Analysis of PREVENT-AD data - baseline fMRI and follow-up PET with genotype markers & cognition
#Code by Eóin N. Molloy, PhD - MultiModal Neuroimaging Lab, DZNE Magdeburg

#########################################################################################################################

#Import packages:
library(dplyr)
library(lme4)
library(car)
library(lmerTest)
library(MuMIn)
library(sjstats)
library(ggpubr)
library(devtools)
library(Lahman)
library(ppcor)

#########################################################################################################################

#0A Set the working directory:

setwd('/Users/eoin/Documents/PREVENT-AD/Analyses/PREVENT-AD_Data_Analayses')

#0B Read in the data first:

#First the longitudinal precuneus ROI data:
prec <-read.csv("Precunues.csv", header = T)
#Second the RBANS cognitive data:
rbans <-read.csv("RBANS_Delayed_Memory_All_Times.csv", header = T)
#However, just keep it to the times matched by the rest of the data:
rbans <- subset(rbans, Time != "FU36" & Time != "FU60" & Time != "FU72" & Time != "FU84")
#Third the corrected hit rate behavioural data from the scan sessions:
corrhr <-read.csv("CorrHR.csv", header = T)
#Fourth, our main dataframe with all fMRI ROI and PET data:
base <-read.csv("PAD_Masterfile_BasefMRI_PET.csv", header = T)

#########################################################################################################################

#SECTION: Assessment of baseline associations between precuneus activation and cognition

#Note --> NaN have been subset out for subjects with missing data points

#1A RBANS vs Precuneus BOLD Retrieval
nonan <- subset(base, Subject != "sub-MTL0384" & Subject != "sub-MTL0390")
prec_v_RBANS <- nonan[,c("RBANS_delayed_memory_index_score_Baseline", "fMRI_Precuneus_Baseline", "Age_Baseline_Months", "Sex", "Education_Baseline_Years")]
pcor(prec_v_RBANS, method = "spearman") #no correlation (rho=-0.01998261, p = 0.80) 

#1B CorrHR vs Precuneus BOLD Retrieval
prec_v_CorrHR <- nonan[,c("CorrHR_Baseline", "fMRI_Precuneus_Baseline", "Age_Baseline_Months", "Sex", "Education_Baseline_Years")]
pcor(prec_v_CorrHR, method = "spearman") #no correlation (rho= 0.13256168, p = 0.09)

#1C MoCA vs Precuneus BOLD Retrieval
nonan2 <- subset(base, Subject != "sub-MTL0550")
prec_v_moca <- nonan2[,c("MOCA_Total_Score", "fMRI_Precuneus_Baseline", "Age_Baseline_Months", "Sex", "Education_Baseline_Years")]
pcor(prec_v_moca, method = "spearman") #no correlation (rho= 0.07525192  , p = 0.34)

#########################################################################################################################

#SECTION: Assessment of longitudinal precuneus activation and longitudinal episodic memory retrieval performance

#2A Change in precuneus activation over time:
PREC <- lmer(PREC_Bilat ~ Time + Age + Sex + Education + (1 |Subject), data=prec, REML = F)
anova(PREC, type = 3) #Significant effect of time

#2B Change in "Corrected Hit Rate" retrieval performance from the scanner:
CorrHR <- lmer(HR_corr ~ Time + Age + Sex + Education + (1 |Subject), data=corrhr, REML = F)
anova(CorrHR, type = 3) #Significant effect of time and sex

#2C Change in delayed memory score on the RBANS:
RBANS <- lmer(delayed_memory_index_score ~ Time + Age + Sex + Education + (1 |Subject), data=rbans, REML = F)
anova(RBANS, type = 3) #Significant effect of time, sex, and education

#########################################################################################################################

#SECTION: Assessment of baseline and longitudinal precuneus activation effects on AD pathology

#First we need to account for the non-normal amyloid data distributions:
#Box-Cox correction on the whole brain amyloid data:

#compute a linear model and pass it to the boxcox function:
mod <- boxcox(lm(base$Amyloid_WB ~ 1))
#calculate lambda
lambda <- mod$x[which.max(mod$y)]
lambda
#apply lambda to data
Amyloid_WB_BoxCox <- (base$Amyloid_WB ^ lambda - 1) / lambda
#how does it look?
hist(Amyloid_WB_BoxCox)

#3A Baseline prec on amyloid:
PRECvsAMY <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Baseline + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsAMY) # effect precuneus

#3B Baseline prec on tau:
PRECvsTAU <- lm(Entorhinal_Tau_PET_LR_MEAN ~ fMRI_Precuneus_Baseline + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU) #significant sex difference

#3C Slope (change in) prec on amyloid: 
PRECvsAMY_SLOPE <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Slope + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsAMY_SLOPE) #  effect slope

#3D Slope (change in) prec on tau: 
PRECvsTAU_SLOPE <- lm(Entorhinal_Tau_PET_LR_MEAN ~ fMRI_Precuneus_Slope + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU_SLOPE) #significant sex difference

#########################################################################################################################

#SECTION: Assessment of baseline and longitudinal precuneus activity, AD pathology and APOE genotype:

#Repeat the above, but now with carrier status as a grouping factor:
#Set carrier status as factor:
base$APOE4_Group <- as.factor(base$APOE4_Group) #set as factor

#4A Baseline prec on amyloid:
PRECvsAMY_GEN <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Baseline*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsAMY_GEN) # effect precuneus, effect genotype, interaction

#4B Baseline prec on tau:
PRECvsTAU_GEN <- lm(Entorhinal_Tau_PET_LR_MEAN ~ fMRI_Precuneus_Baseline*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU_GEN) #significant sex difference

#4C Slope (change in) prec on amyloid: 
PRECvsAMY_SLOPE_GEN <- lm(Amyloid_WB_BoxCox ~ fMRI_Precuneus_Slope*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsAMY_SLOPE_GEN) # effect precuneus, effect genotype, interaction

#4D Slope (change in) prec on tau: 
PRECvsTAU_SLOPE_GEN <- lm(Entorhinal_Tau_PET_LR_MEAN ~ fMRI_Precuneus_Slope*APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex + BaseMRI_to_Amy_PET_Months + Prec_GM_Volume, data=base)
anova(PRECvsTAU_SLOPE_GEN) #significant sex difference

#Now we test for genotype group differences in precuneus activity and whole brain amyloid:

#4E Group differences in baseline activation?
baseprec <- lm(fMRI_Precuneus_Baseline ~ APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex + Prec_GM_Volume, data=base)
summary(baseprec) #significant sex difference

#4F Group differences in activation over time?
slopeprec <- lm(fMRI_Precuneus_Slope ~ APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex + Prec_GM_Volume, data=base)
summary(slopeprec) #no effects

#4G Group differences in PET Amyloid?
samy <- lm(Amyloid_WB_BoxCox ~ APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex, data=base)
summary(samy) #Significant group difference

#4H Group differences in PET Entorhinal Tau?
tau <- lm(Entorhinal_Tau_PET_LR_MEAN ~ APOE4_Group + Age_Baseline_Months + Education_Baseline_Years + Sex, data=base)
summary(tau) #Significant group and sex difference

#########################################################################################################################

#SECTION: Relationship between baseline precuneus activity, APOE genotype, and cognition:

#5A Baseline activity, genotype, and time on RBANS delayed memory:
PREC_GEN_TIME_RBANS <- lmer(delayed_memory_index_score ~ Time*APO4_Group*PREC_Bilat_Baseline + Age + Sex + Education + (1 |Subject), data=rbans, REML = F)
Anova(PREC_GEN_TIME_RBANS) #Effect of time, sex, education, and sig 3-way interaction

#5B Baseline activity, genotype, and time on Corrected Hit Rate performance:
PREC_GEN_TIME_CorrHR <- lmer(HR_corr ~ Time*APO4_Group*PREC_Bilat_Baseline + Age + Sex + Education + (1 |Subject), data=corrhr, REML = F)
Anova(PREC_GEN_TIME_CorrHR) #Effect of time, sex

#########################################################################################################################