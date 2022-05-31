# Zero-Covid
Code for: Effect of changing zero-COVID strategy on Omicron and future SARS-CoV-2 wave in mainland China

## Citation
Effect of changing zero-COVID strategy on Omicron and future SARS-CoV-2 wave in mainland China

## Abstract
To adapt to the evolution of SARS-CoV-2, the zero-COVID strategy in China is changing from strict lockdown to population-level testing and has been successful in combating the SARS-CoV-2 spread. However, it is unknown whether and how long this strategy still works for Omicron and future wave. To address this question, we develop a meta-population transmission dynamic model including population-level testing based on 398 million travel movements among 305 Chinese cities. We show that Omicron outbreak is difficult to contain by the control measures implemented during the first COVID-19 wave and vaccines. Travel ban has limited effects on preventing the dissemination of SARS-CoV-2, but it can significantly reduce the epidemic duration if carried out in an appropriate time window. As an active surveillance measure, population-level testing is critical to curtail the nationwide Omicron wave and future SARS-CoV-2 wave. Before effective vaccines and anti-virus treatments coming into the market, zero-COVID strategy is still necessary to reduce the COVID-19 burden in China. 

## Notes on the code
To run the fitting, you need a Matlab toolbox called "DRAM": DRAM is a combination of two ideas for improving the efficiency of Metropolis-Hastings type Markov chain Monte Carlo (MCMC) algorithms, Delayed Rejection and Adaptive Metropolis. This page explains the basic ideas behind DRAM and provides examples and Matlab code for the computations.(see http://helios.fmi.fi/~lainema/dram/)

About code folder: We used ”Fitting” folder to estimate parameters for the first wave of COVID-19 in 2020 in mainland China, “SimulationPrevious” folder to perform simulations evaluating the effect of previous zero-COVID strategy (social distancing, travel ban and vaccination) during the Omicron wave in China, and “SimulationTesting” folder to perform simulations investigating the effect of population-level testing with or without travel ban on Omicron and future SARS-CoV-2 wave. The setting of the fixed parameters can be found in the main text and supplemented material.

## Data
### The COVID-19 case data
We collected the daily official case reports from the health commission of 34 provincial-level administrative units and 305 city-level units, the website’s links are provided 9. Demographic data for each city were collected from the China City Statistical Yearbook 2019 (http://olap.epsnet.com.cn/).

### Human mobility data
Human movement in China can be observed directly from mobile phone data, through the Baidu location-based services (LBS). Both the recorded movements and relative volume of inflows, outflows, inner movement among cities (n=30558), were obtained from the migration flows database (http://qianxi.baidu.com/) from 1 January 2019 to 1 April 2020. 
