# Comparison of MR-MOE vs Two-Sample MR in Identifying Causal Risk Factors for Alzheimer’s Disease

Esimate the causal relationships between potentially modifiable risk factors and AD using Mendelian randomization and comparing approaches.

**MR Analysis**

Mendelian randomization analysis was conducted using the TwoSampleMR package. The MR-MoE approach was used to determine the best method for Two-Sample MR.

**Risk Factors**

Risk factors of choice include education, high total cholesterol levels, and social isolation.

**Results**

Education - Fixed Estimates IVW produced the highest AUROC score, with an odds ratio of 0.69. This means that a person with higher education is likely to not get AD.

**Background and Objective:**
The global prevalence of Alzheimer’s Disease is currently estimated to affect forty-four million people, with numbers doubling each year, primarily for those in the 65+ age group. Previous studies have demonstrated a strong relationship exists between Alzheimer’s Disease (AD) and genetics. Some important findings indicate that variance in the APOE, CLU, and SORL1 genes (all responsible for lipid metabolism) are associated with AD. Further, AD has been found to be common among people with cardiovascular disease, smoking, and obesity, among other common factors. Identifying causal risk factors for AD is the next step in creating interventions that can prevent the emergence of the disease. A methods gain in populaty for identifying such risk factors is Mendelian Randomization (MR), which involves using data from previously conducted GWAS to determine if a relationship exists between disease and the risk factor. The genotype chosen must affect the disease indirectly through its exposure. Unlike epidemiological studies (which involves studying subjects over a long period of time to see if a person with a certain condition developed a disease), since genotypes are randomly assigned, the distribution itself is free of confounding factors. This means is easier to directly identify variants that lead to AD.

*Mendelian Randomization – Mixture of Experts (MR-MOE)*, recently developed by *Hemani et. al in 2017*, is a method used to determine the most suitable MR test through a machine learning algorithm. Data is provided to several methods and the one ‘expert’ that is most reliable is the one selected using a performance estimate. This method can improve power and lower rates of false discovery associated with single method studies. 

Using conclusions and data from a previously conducted two-sample MR (Andrews et. al, 2021), this project will entail using MR-MOE to determine the effectiveness of the program for potential use in identifying risk factors. In particular, education and cholesterol levels, using summary statistics from Social Science Genetic Association Consortium (SSGAS) and Global Lipids Genetics Consortium (GLGC) respectively as well as data from a GWAS the million veterans program, will be the measurements used to begin conducting the MR-MOE analysis.

MR-MOE will then be used to evaluate the causal associations between other key AD risk factors, such as Blood Pressure, Cigarettes per Day, and Alcohol Consumption. If time permits, some additional exposures of interest are the genetics of mental illnesses, such as loneliness and depression, as potential risk factors for Alzheimer’s Disease, which will involve using the MR-MOE approach. Another possible exposure is drug addiction. 
