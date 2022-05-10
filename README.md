# SIMEX to correct for bias due to interval-censored events in survival analysis with time-varying exposure 

R code to implement for a given dataset the proposed methods in the article:  
"Abrahamowicz M, Beauchamp ME, Moura CS, Bernatsky S, Ferreira Guerra S, Danieli C. Adapting SIMEX to correct for bias due to interval-censored outcomes in survival analysis with time-varying exposure. *Biom J* 2022".

For questions or comments about the code contact Marie-Eve Beauchamp (marie-eve.beauchamp at rimuhc.ca).

R code to reproduce the simulations presented in the article can be found as Supporting Information at doi/10.1002/bimj.2021xxxxx. 

## Content

#### `Functions_SIMEXforIntervalCensoredOutcomes.R`
Code of the functions to implement the proposed SIMEX-like method for a given dataset, including the diagnostic plot for choosing the extrapolating function and estimating the bootstrap CIs.

#### `Example_Use.R`
Program illustrating the implementation of the proposed methods for a given dataset available in `Data_for_Example.RData`. This program calls `Functions_SIMEXforIntervalCensoredOutcomes.R`.

#### `Data_for_Example.RData`
File including the dataset (dat.end) used in the program `Example_Use.R`. It also includes an object (visits) listing the visit times of each subject. Both the dataset and visits were simulated.
