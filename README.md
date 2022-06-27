# Zero-Covid
Code for: Quantitative assessment of zero-COVID strategy in mainland China

## Citation
Quantitative assessment of zero-COVID strategy in mainland China

## Abstract
Although the zero-COVID strategy in China has been successful in combating series of SARS-CoV-2 variants, arguments for its benefit and negative consequence are increasing but seldom quantitatively assessed. Here, we develop a metapopulation model, based on 419 million travel movements among 366 Chinese cities with population-level testing, contact contracting and isolating, to explore sustainability of current strategy. We find that original control measures hardly contain Omicron outbreak even with full-dose vaccination and lockdown  . As a pillar of zero-COVID strategy, instituting population-level testing with a less than 3 days interval is critical to achieve zero-COVID under Omicron wave even lifting lockdown  . The expected value of perfect information analysis indicates that the 2-day testing interval is optimal after considering the uncertainty in transmissibility of the future wave. Our findings provide important information for future self-monitoring based pandemic responsiveness.

## Notes on the code
To run the fitting, you need a Matlab toolbox called "DRAM": DRAM is a combination of two ideas for improving the efficiency of Metropolis-Hastings type Markov chain Monte Carlo (MCMC) algorithms, Delayed Rejection and Adaptive Metropolis. This page explains the basic ideas behind DRAM and provides examples and Matlab code for the computations.(see http://helios.fmi.fi/~lainema/dram/)

About code folder: We used ”Fitting” folder to estimate parameters for the first wave of COVID-19 in 2020 in mainland China, “Simulation2020” folder to perform simulations evaluating the effect of zero-COVID strategy implemented in 2020 (social distancing, travel ban) on Omicron wave, and “SimulationTesting” folder to perform simulations investigating the effect of population-level testing on Omicron and future SARS-CoV-2 wave. The setting of the fixed parameters can be found in the main text and supplementary material.

## Data
### Epidemiological data
We collected the daily official case reports from the health commission of 34 provincial-level administrative units and 36605 city-level units, the website’s links are provided 16. Demographic data for each city were collected from the China City Statistical Yearbook 2019 (http://olap.epsnet.com.cn/).


### Human mobility data
Human movement in China can be observed directly from mobile phone data, through the Baidu location-based services (LBS). Both the recorded movements and relative volume of inflows, outflows, inner movement among cities (n=366), were obtained from the migration flows database (http://qianxi.baidu.com/) from 1 January 2019 to 1 April 2020. We consider the averaged travel flow in 2019 as the flow of baseline, and construct a flow network across 366 cities in China to simulate the SARS-CoV-2 transmission across cities in China.
