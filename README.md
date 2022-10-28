# Zero-Covid
Code for: Quantitative assessment of zero-COVID strategy in mainland China

## Citation
Criteria for maintenance of the Zero-COVID strategy in mainland China: a large-scale modelling study


## Abstract
Background: During the past three years of the COVID-19 pandemic, many countries have experienced disruptions and collapses in their health-care systems. Rethinking population heath, broader than health-care system alone, is a prerequisite for health preparedness. As China continues to be challenged by successive SARS-CoV-2 variants, the feasibility of successfully maintaining the Zero-COVID strategy and the disease burden have not been quantitatively assessed.

Methods: We developed an age-stratified metapopulation model with testing, contact tracing and isolation based on 419 million travel movements among 366 Chinese cities to explore the feasibility of maintaining the Zero-COVID strategy and the corresponding disease burden with consideration of the prevalence of underlying health conditions in the Chinese population, including hypertension, diabetes mellitus, chronic kidney disease, obesity, chronic obstructive pulmonary disease, asthma, chronic liver disease and cancer. 

Findings: We find that, if the Zero-COVID strategy is to be successful in suppressing SARS-CoV-2 and COVID-19, including Omicron variants, population-wide viral testing should be carried out at an interval of ≤ 3 days, given current levels of vaccination, even without lockdown for most cities. Satisfying that criterion would keep hospital bed occupancy at less than 6% of the capacity for respiratory illness, and ICU bed occupancy at less than 44% of the capacity for respiratory illness, during any COVID-19 outbreak. Relaxing the testing interval to 4 days would still keep the COVID-19 burden under the capacity, but lead to the endemic of COVID-19.  Without testing, bed capacity could be exceeded 3.5-fold and ICU bed capacity 25-fold. In addition, the underlying health conditions among the SARS-CoV-2 infections would lead to an extra 3.81 million hospital admissions (2276% increase over the scenario without underlying conditions) and 0.51 million ICU admissions (2830% increase over the scenario without underlying conditions  )   under the scenario of the least stringent Zero-COVID strategy, with the major contribution from the ≥ 60 age group compared with other age groups. Among the underlying health conditions, obesity is expected to contribute to most extra hospital admissions, while hypertension could contribute to most extra ICU admissions. 

Interpretation: Even with China’s high vaccination coverage, SARS-CoV-2 and COVID-19 are unlikely to be contained without frequent population-wide testing and isolation linked to every outbreak. Our findings could help guide future COVID-19 control policies in China.


## Notes on the code
To run the fitting, you need a Matlab toolbox called "mcmcstat" available from https://mjlaine.github.io/mcmcstat/. It implements DRAM algorithm, which is a combination of two ideas for improving the efficiency of Metropolis-Hastings type Markov chain Monte Carlo (MCMC) algorithms, Delayed Rejection and Adaptive Metropolis. The mcmcstat toolbox documentation can be found at https://mjlaine.github.io/mcmcstat/.

About code folder: We used ”Fitting2020” folder to estimate parameters for the first wave of COVID-19 in 2020 in mainland China, “Simulation2020” folder to perform simulations evaluating the effect of zero-COVID strategy implemented in 2020 (social distancing, travel ban) on Omicron wave, “SimulationTesting” folder to perform simulations investigating the effect of population-level testing on Omicron wave, and "COVID-19 burden" folder to quantitatively assesse disease COVID-19 burden taking unerlying health condition into account. The setting of the fixed parameters can be found in the main text and supplementary material.

## Data
### Epidemiological data
We collected the daily official case reports from the health commission of 34 provincial-level administrative units and 366 city-level units, the website’s links are provided in Tian, H. et al, 2020. Demographic data for each city were collected from the China City Statistical Yearbook 2019 (http://olap.epsnet.com.cn/).


### Human mobility data
Human movement in China can be observed directly using mobile phone data with Baidu location-based services. We obtained both the recorded movements and relative volume of inflows, outflows, and internal movement among cities (n=366) from the migration flows database (http://qianxi.baidu.com/) from 1 January 2019 to 1 April 2020. We considered the averaged travel flow in 2019 as the flow at baseline and constructed a flow network across 366 cities in China to simulate SARS-CoV-2 transmission across Chinese cities.
