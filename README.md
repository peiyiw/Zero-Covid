# Zero-Covid
Code for: Marginal effects of public health measures and COVID-19 disease burden in China: a large-scale modelling study



## Citation
Marginal effects of public health measures and COVID-19 disease burden in China: a large-scale modelling study



## Abstract
China had conducted some of the most stringent public health measures to control the spread of successive SARS-CoV-2 variants. However, the effectiveness of these measures and their impacts on the associated disease burden have rarely been quantitatively assessed at the national level. To address this gap, we developed a stochastic age-stratified metapopulation model that incorporates testing, contact tracing and isolation, based on 419 million travel movements among 366 Chinese cities. The study period for this model began from September 2022. The COVID-19 disease burden was evaluated, considering 8 types of underlying health conditions in the Chinese population. We identified the marginal effects between the testing speed and reduction in the epidemic duration. The findings suggest that assuming a vaccine coverage of 89%, the Omicron-like wave could be suppressed by 3-day interval population-level testing (PLT), while it would become endemic with 4-day interval PLT, and without testing, it would result in an epidemic. PLT conducted every 3 days would not only eliminate infections but also keep hospital bed occupancy at less than 29.46% (95% CI, 22.73-38.68%) of capacity for respiratory illness and ICU bed occupancy at less than 58.94% (95% CI, 45.70-76.90%) during an outbreak. Furthermore, the underlying health conditions would lead to an extra 2.35 (95% CI, 1.89-2.92) million hospital admissions and 0.16 (95% CI, 0.13-0.2) million ICU admissions. Our study provides insights into health preparedness to balance the disease burden and sustainability for a country with a population of billions.

## System requirements
### Versions of the software
The COVID-19 burden analysis was performed using R 4.0.5, and the fitting/simulation was performed using MATLAB R2021b on a Microsoft Windows Server version 1607(14393.447).  
It takes about 3 hours to fit the model and 8 hours to run the simulation on a machine with the recommended specs.

### Matlab dependencies
To run the fitting, you need a Matlab toolbox called "mcmcstat" available from https://mjlaine.github.io/mcmcstat/. It implements DRAM algorithm, which is a combination of two ideas for improving the efficiency of Metropolis-Hastings type Markov chain Monte Carlo (MCMC) algorithms, Delayed Rejection and Adaptive Metropolis. The mcmcstat toolbox documentation can be found at https://mjlaine.github.io/mcmcstat/.

### R dependencies
Users should install the following packages prior to performing COVID-19 burden analysis, from an R terminal:
```
install.packages(c('ggplot2', 'openxlsx', 'lubridate', 'reshape2', 'scales', 'Rmisc', 'ggpubr', 'dplyr', 'ggforce', 'RColorBrewer'))
```

## Notes on the code


About code folder: We used ”Fitting2020” folder to estimate parameters for the first wave of COVID-19 in 2020 in mainland China, “Simulation2020” folder to perform simulations evaluating the effect of public health measures implemented in 2020 (social distancing, travel ban) on Omicron wave, “SimulationTesting” folder to perform simulations investigating the effect of population-level testing on Omicron wave, and "COVID-19 burden" folder to quantitatively assesse disease burden burden of COVID-19 on Omicron wave taking unerlying health condition into account. The detailed information on the epidemiological data of the first wave of COVID-19 in 2020 in mainland China were stroed in the "Data" folder.The setting of the fixed parameters can be found in the main text and supplementary material. The chain of MCMC (200,000 were reserved) of fitting and the results of simulation of testing were stroed in the "Output" folder.


## Data
### Epidemiological data
We collected the daily official case reports from the health commission of 34 provincial-level administrative units and 366 city-level units, the website’s links are provided in Tian, H. et al, 2020. Demographic data for each city were collected from the China City Statistical Yearbook 2019 (http://olap.epsnet.com.cn/).


### Human mobility data
Human movement in China can be observed directly using mobile phone data with Baidu location-based services. We obtained both the recorded movements and relative volume of inflows, outflows, and internal movement among cities (n=366) from the migration flows database (http://qianxi.baidu.com/) from 1 January 2019 to 1 April 2020. We considered the averaged travel flow in 2019 as the flow at baseline and constructed a flow network across 366 cities in China to simulate SARS-CoV-2 transmission across Chinese cities.
