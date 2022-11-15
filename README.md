# Zero-Covid
Code for: Criteria for the maintenance of zero COVID in China: a large-scale modelling study



## Citation
Criteria for the maintenance of zero COVID in China: a large-scale modelling study



## Abstract
During the past three years of the COVID-19 pandemic, many countries have experienced disruptions and collapses in their health-care systems. Rethinking population health, broader than health-care system alone, is a prerequisite for health preparedness. As China continues to be challenged by successive SARS-CoV-2 variants, the feasibility of maintaining the Zero-COVID strategy and the disease burden have not been quantitatively assessed. We developed an age-stratified metapopulation model with testing, contact tracing and isolation based on 419 million travel movements among 366 Chinese cities to explore the feasibility of maintaining zero COVID and the corresponding disease burden with consideration of the prevalence of underlying health conditions in the Chinese population, including hypertension, diabetes mellitus, chronic kidney disease, obesity, chronic obstructive pulmonary disease, asthma, chronic liver disease and cancer. We find that, if the Zero-COVID strategy aims at suppressing SARS-CoV-2 and COVID-19, including Omicron variants, population-wide viral testing at city level should be carried out at an interval of ≤ 3 days, given current levels of vaccination, even without lockdown for most cities. Satisfying that criterion would not only eliminate infection but also keep hospital bed occupancy at less than 6% of capacity for respiratory illness and ICU bed occupancy at less than 44% of capacity during an outbreak. Relaxing the testing interval to 4 days would keep the COVID-19 burden within hospital capacity, but with an endemic. Without testing, hospital bed capacity could be exceeded 3.5-fold and ICU bed capacity 25-fold. In addition, under a 3-day testing interval, the underlying health conditions among the SARS-CoV-2 infections would lead to an extra 3.81 million hospital admissions and 0.51 million ICU admissions, corresponding to a 2.3-fold and 2.8-fold of that under the scenario without underlying conditions, respectively. The major proportion (≥ 50%) of extra COVID-19 burden would come from the ≥ 60 age group. Among the underlying health conditions, obesity is expected to contribute to most extra hospital admissions, while hypertension could contribute to most extra ICU admissions. 

## System requirements
### Versions of the software
The COVID-19 burden analysis was performed using R4.0.5, and the fitting/simulation was performed using MATLAB R2021b on a Microsoft Windows Server version 1607(14393.447).  
It takes about 3 hours to fit the model and 20 minitures to run a simulation on a machine with the recommended specs.

### Matlab dependencies
To run the fitting, you need a Matlab toolbox called "mcmcstat" available from https://mjlaine.github.io/mcmcstat/. It implements DRAM algorithm, which is a combination of two ideas for improving the efficiency of Metropolis-Hastings type Markov chain Monte Carlo (MCMC) algorithms, Delayed Rejection and Adaptive Metropolis. The mcmcstat toolbox documentation can be found at https://mjlaine.github.io/mcmcstat/.

### R dependencies
Users should install the following packages prior to performing COVID-19 burden analysis, from an R terminal:
-install.packages(c('ggplot2', 'abind', 'irlba', 'knitr', 'rmarkdown', 'latex2exp', 'MASS', 'randomForest'))
  
## Notes on the code


About code folder: We used ”Fitting2020” folder to estimate parameters for the first wave of COVID-19 in 2020 in mainland China, “Simulation2020” folder to perform simulations evaluating the effect of zero-COVID strategy implemented in 2020 (social distancing, travel ban) on Omicron wave, “SimulationTesting” folder to perform simulations investigating the effect of population-level testing on Omicron wave, and "COVID-19 burden" folder to quantitatively assesse disease burden burden of COVID-19 on Omicron wave taking unerlying health condition into account. The setting of the fixed parameters can be found in the main text and supplementary material.

## Data
### Epidemiological data
We collected the daily official case reports from the health commission of 34 provincial-level administrative units and 366 city-level units, the website’s links are provided in Tian, H. et al, 2020. Demographic data for each city were collected from the China City Statistical Yearbook 2019 (http://olap.epsnet.com.cn/).


### Human mobility data
Human movement in China can be observed directly using mobile phone data with Baidu location-based services. We obtained both the recorded movements and relative volume of inflows, outflows, and internal movement among cities (n=366) from the migration flows database (http://qianxi.baidu.com/) from 1 January 2019 to 1 April 2020. We considered the averaged travel flow in 2019 as the flow at baseline and constructed a flow network across 366 cities in China to simulate SARS-CoV-2 transmission across Chinese cities.
