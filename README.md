# SIMEXforIntervalCensoredOutcomes-BiomJ2022
R code for article "Abrahamowicz M, Beauchamp ME, Moura CS, Bernatsky S, Ferreira Guerra S, Danieli C. Adapting SIMEX to correct for bias due to interval-censored outcomes in survival analysis with time-varying exposure. *Biom J* 2022".
## Description
For questions or comments about the code contact Marie-Eve Beauchamp (marie-eve.beauchamp at rimuhc.ca).
 
The code has been written using R with the following version information:<br/>
- R version 4.1.1 (2021-08-10)<br/> 
- Platform x86_64-w64-mingw32/x64 (64-bit)<br/> 
- Using R packages:<br/> 
  - PermAlgo version 1.1<br/>
  - survival version 3.2-11
## Content
### Directory: Simulations 
#### `Simulations_and_Results.R`
Code to reproduce the simulations, as well as to extract the results for the tables and create the figures, presented in the manuscript and Supporting Information.

#### `Functions_for_Simulations.R`
Includes the functions needed for the *simulations*. This program is called by 'Simulations_and_Results.R'. For functions to implement the proposed methods for a given dataset, see 'Functions_SIMEXforIntervalCensoredOutcomes.R' in the directory Use.

### Directory: Use 
#### 'Functions_SIMEXforIntervalCensoredOutcomes.R'
Includes the functions to implement the proposed methods for a given dataset.

#### 'Example_Use.R'
Program illustrating the implementation of the proposed methods for a given dataset. This program calls 'Functions_SIMEXforIntervalCensoredOutcomes.R'.

#### 'Data_for_Example.RData'
This file includes a dataset (dat.end) used in the example of use of the proposed methods presented in the program 'Example_Use.R'. It also includes an object (visits) listing the visit times of each subject. Both the dataset and visits were simulated.
