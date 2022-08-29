# Zero-Covid
Code for: Quantitative assessment of zero-COVID strategy in mainland China

## Citation
Quantitative assessment of zero-COVID strategy in mainland China

## Abstract
Although the zero-COVID strategy in China has temporarily confronted several SARS-CoV-2 variants, arguments for its benefits and disadvantages are increasing but have rarely been quantitatively assessed. Here, we developed a metapopulation model with testing-contact tracing-isolation, including 419 million travel movements among 366 Chinese cities, to explore sustainability of the zero-COVID strategy. The analysis shows that instituting faster testing with diminishing marginal utility, particularly in less than 3-day interval, is critical to achieving zero-COVID under the benchmark of epidemic duration, even with lockdowns lifted. Strikingly, developing regions may be more vulnerable to current strategy, given higher proportion of people isolated, than developed ones. Despite high vaccination coverage, control measures without testing could scarcely contain the Omicron and more transmissible variant. Our findings provide evidence for future pandemic responsiveness based on self-monitoring.

## Notes on the code
To run the fitting, you need a Matlab toolbox called "mcmcstat" available from https://mjlaine.github.io/mcmcstat/. It implements DRAM algorithm, which is a combination of two ideas for improving the efficiency of Metropolis-Hastings type Markov chain Monte Carlo (MCMC) algorithms, Delayed Rejection and Adaptive Metropolis. The mcmcstat toolbox documentation can be found at https://mjlaine.github.io/mcmcstat/.

About code folder: We used ”Fitting2020” folder to estimate parameters for the first wave of COVID-19 in 2020 in mainland China, ”FittingTesting” folder to demostrate the ability of metapopulation model with testing-contact tracing-isolation by reconstructing transmission dynamic in the first city implementing population-level testing, “Simulation2020” folder to perform simulations evaluating the effect of zero-COVID strategy implemented in 2020 (social distancing, travel ban) on Omicron wave, and “SimulationTesting” folder to perform simulations investigating the effect of population-level testing on Omicron and future SARS-CoV-2 wave. The setting of the fixed parameters can be found in the main text and supplementary material.

## Data
### Epidemiological data
We collected the daily official case reports from the health commission of 34 provincial-level administrative units and 366 city-level units, the website’s links are provided in Tian, H. et al, 2020. Demographic data for each city were collected from the China City Statistical Yearbook 2019 (http://olap.epsnet.com.cn/).


### Human mobility data
Human movement in China can be observed directly using mobile phone data with Baidu location-based services. We obtained both the recorded movements and relative volume of inflows, outflows, and internal movement among cities (n=366) from the migration flows database (http://qianxi.baidu.com/) from 1 January 2019 to 1 April 2020. We considered the averaged travel flow in 2019 as the flow at baseline and constructed a flow network across 366 cities in China to simulate SARS-CoV-2 transmission across Chinese cities.
