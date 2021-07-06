# Comparison of MR-MOE vs two-sample MR in Identifying Causal Risk Factors for Alzheimer’s Disease

**Background and Objective:**
The global prevalence of Alzheimer’s Disease is currently estimated to affect forty-four million people, with numbers doubling each year, primarily for those in the 65+ age group. Previous studies have demonstrated a strong relationship exists between Alzheimer’s Disease (AD) and genetics. Some important findings indicate that variance in the APOE, CLU, and SORL1 genes (all responsible for lipid metabolism) are associated with AD. Further, AD has been found to be common among people with cardiovascular disease, smoking, and obesity, among other common factors. Identifying causal risk factors for AD is the next step in creating interventions that can prevent the emergence of the disease. A methods gain in populaty for identifying such risk factors is Mendelian Randomization (MR), which involves using data from previously conducted GWAS to determine if a relationship exists between disease and the risk factor. The genotype chosen must affect the disease indirectly through its exposure. Unlike epidemiological studies (which involves studying subjects over a long period of time to see if a person with a certain condition developed a disease), since genotypes are randomly assigned, the distribution itself is free of confounding factors. This means is easier to directly identify variants that lead to AD.

*Mendelian Randomization – Mixture of Experts (MR-MOE)*, recently developed by *Hemani et. al in 2017*, is a method used to determine the most suitable MR test through a machine learning algorithm. Data is provided to several methods and the one ‘expert’ that is most reliable is the one selected using a performance estimate. This method can improve power and lower rates of false discovery associated with single method studies. 

Using conclusions and data from a previously conducted two-sample MR (Andrews et. al, 2021), this project will entail using MR-MOE to determine the effectiveness of the program for potential use in identifying risk factors. In particular, education and cholesterol levels, using summary statistics from Social Science Genetic Association Consortium (SSGAS) and Global Lipids Genetics Consortium (GLGC) respectively as well as data from a GWAS the million veterans program, will be the measurements used to begin conducting the MR-MOE analysis.

MR-MOE will then be used to evaluate the causal associations between other key AD risk factors, such as Blood Pressure, Cigarettes per Day, and Alcohol Consumption. If time permits, some additional exposures of interest are the genetics of mental illnesses, such as loneliness and depression, as potential risk factors for Alzheimer’s Disease, which will involve using the MR-MOE approach. Another possible exposure is drug addiction.  

**Planned analysis by Aims:**
GWAS Summary statistics used for analysis will be the same as those previously used to conduct two-sample MR. For example, total cholesterol levels data will come from the GLGC and for education, from the SSGAS.

## Aim 1: Use MR-MoE to evaluate causal associations between education and total cholesterol and AD

**Hypothesis 1a:**
High levels of education is associated with a reduced risk of Alzheimer’s Disease. 

**Hypothesis 1b:**
High levels of total cholesterol is associated with an increased risk of Alzheimer’s Disease. 

**Rationale:**
In previously conducted two-sample MR, education was associated with a reduced risk of AD using the IVW, MR-Egger, Weight Median, and Weight Mode two-sample MR methods. High levels of total cholesterol are associated with an increased risk of neuritic plaque buildup and VBI used the same two-sample MR methods mentioned above. MR-MOE will identify the best approach to use for future studies that best corresponds with the data of interest and gives the highest statistical power. 

**Primary Proposed Analysis:**
MR-MOE will report the most appropriate MR-approach. Whichever method is chosen by the algorithm will then be used to conduct the MR analysis. The significance levels are then compared to the other techniques used previously (IVW, MR-Egger, Weight Median, Weight Mode). 

**Anticipated Problems and Alternative Approaches:**
The following applies to all MR-MoE analyses: they may be slower to run and the program does not give a reason as to why a particular method was selected over another. Therefore, the alternative is to use a combination  of 2-sample MR methods (IVW, MR-Egger, Weight Median, Weight Mode) thought to best apply to a set of risk factors by the researcher. 

## Aim 2:Use MR-MoE to evaluate the causal association between key AD risk factors and AD

**Hypothesis:**
Other key risk factors, such as smoking initiation, blood pressure, and sleep duration, are also associated with AD. In previous studies, Smoking initiation has been associated with decreased cortical thickness, blood pressure with increased risk for VBI, and sleep duration with increased cortical thickness.

**Rationale:**
MR-MOE will identify the best approach that best corresponds with the data of interest and gives the highest statistical power. The method can then be used to determine if there is an association between the risk factor and AD.

**Primary Proposed Analysis:**
MR-MOE will report the most appropriate MR-approach. Whichever method is chosen by the algorithm will then be used to conduct the MR analysis. The significance levels are then compared to the other techniques used previously (IVW, MR-Egger, Weight Median, Weight Mode). 

**Anticipated Problems and Alternative Approaches:**
The following applies to all MR-MoE analyses: they may be slower to run and the program does not give a reason as to why a particular method was selected over another. Therefore, the alternative is to use a combination  of 2-sample MR methods (IVW, MR-Egger, Weight Median, Weight Mode) thought to best apply to a set of risk factors by the researcher. 

## Aim 3: Expand MR-MOE analysis to Endophenotypes of Interest: Social Isolation

**Hypothesis:**
Social isolation leads to an increased risk of AD. 

**Rationale:**
Various levels of socialization are distributed continuously throughout the population. Particularly, decreased sociability has an demonstrated a relationship exists with several neurodegenerative disorders, including AD. 

**Primary Proposed Analysis:** 
MR-MOE will report the most appropriate MR-approach. Whichever method is chosen by the algorithm will then be used to conduct the MR analysis. The significance levels are then analyzed to determine what SNPs, if any, are associated with AD. 

**Anticipated Problems and Alternative Approaches:**
The following applies to all MR-MoE analyses: they may be slower to run and the program does not give a reason as to why a particular method was selected over another. Therefore, the alternative is to use a combination  of 2-sample MR methods (IVW, MR-Egger, Weight Median, Weight Mode) thought to best apply to a set of risk factors by the researcher. 

## Gantt Chart 
![image](https://user-images.githubusercontent.com/85951551/122097091-49874380-cddd-11eb-846b-c52b208fc6dd.png)

## References

Andrews, S. J. et al. Causal Associations Between Modifiable Risk Factors and the Alzheimer’s Phenome. Ann Neurol 89, 54–65 (2021).

Bralten, J., Mota, N.R., Klemann, C.J.H.M. et al. Genetic underpinnings of sociability in the general population. Neuropsychopharmacol. (2021). https://doi.org/10.1038/s41386-021-01044-z

Hemani, G. et al. Automating Mendelian randomization through machine learning to construct a putative causal map of the human phenome. (2017).

Kivipelto, M., Mangialasche, F. & Ngandu, T. Lifestyle interventions to prevent cognitive impairment, dementia and Alzheimer disease. Nat Rev Neurol 14, 653–666 (2018).











