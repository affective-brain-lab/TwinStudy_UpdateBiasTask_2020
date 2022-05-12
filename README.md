# TwinStudy_UpdateBiasTask_2020

Data

Data contains the scores of the variables (priors, update, subjective ratings, learning, estimation error and memory score), zygosity, gender and age for each twin.
Twins pairs are identified by having the first digits equal but the last digit different. TwinA's last digit is 1, twinB's last digit is 2. Fr example, 201 refers to twinA, 202 refers to twinB.

Corr_FE_ratings contains the correlation between priors and each ratings and memory score for each subject.

DataMixed contains scores at trial level of priors and subjective ratings.

Code

ACE assesses the genetic contribution of the cognitive and affective factors using the ACE model which is an epidemiological model stating that genetic factors (A), common environmental factors (C) and specific environmental factors (E) explain individual differences in a given phenotype. The model enables us to quantify the contribution of each factor by comparing the covariance of scores in the phenotype of interest across monozygotic and dizygotic twins. The script run on R-3.5.0 and requires data in .csv. 
