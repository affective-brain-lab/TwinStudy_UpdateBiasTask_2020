# TwinStudy_UpdateBiasTask_2020

Data

Data contains the scores of the variables (priors, update, subjective ratings, learning, estimation error and memory score), zygosity and age for each twin.
Twins pairs are identified by having the first digits equal but the last digit different. TwinA's last digit is 1, twinB's last digit is 2. Fr example, 201 refers to twinA, 202 refers to twinB.

Corr_FE_ratings contains the correlation between priors and each ratings and memory score for each subject.

DataMixed contained scores at trial level of priors and subjective ratings.

Code

ComputeFalconers takes correlation coefficients from MZ and DZ pairs and compute the Falconer's formula

Compute Part Corr compute partial correlation between TwinA and TwinB scores. Since assignment of twins in each pair to the abscissa (i.e., twin A) or ordinate (i.e., twin B) is arbitrary and this assignment influences the correlation value, we performed a permutation analysis. Specifically, we randomly reassigned each twin in the pair as twin A or twin B 10,000 times, each time obtaining a correlation coefficient for the monozygotic samples and for the dizygotic sample.

ACE assesses the genetic contribution of the cognitive and affective factors using the ACE model which is an epidemiological model stating that genetic factors (A), common environmental factors (C) and specific environmental factors (E) explain individual differences in a given phenotype. The model enables us to quantify the contribution of each factor by comparing the covariance of scores in the phenotype of interest across monozygotic and dizygotic twins. The script run on R-3.5.0 and requires data in .csv. 
