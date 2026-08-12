# mixture of varying coefficient compositional linear regressions
Code and Data for TT H and HJ C "Mixture of Varying-coefficient Compositional Linear Regression and Its Application to Industrial Structure - PM2.5 Data", in review.

## Directory Structure

industry-pm/                       # R code for the processing of industry structure-pm2.5 data
│   
├── auxiliary functions.R          # auxiliary functions used in the real data analysis
├── industry and others.txt        # the dataset contains economic and weather variables, including industrial structure, GDP, temperature, relative humidity, and wind speed.
├── pm25.txt                       # the pm2.5 data obtained from the Huiju Atmosphere (https://airwise.hjhj-e.com/)
├── pre analysis.R                 # preliminary analysis of the dataset implementing the mixture of scalar-on-composition regressions and plotting the triangaulr prim
└── real analysis.R                # estimate the coefficient functions and visualize the cities on map

simulation/                        # R code for the simulation studies presented in the manuscript
│   
├── case1/                         # the case 1 in simulation study 
│     ├── auxiliary functions1.R   # auxiliary functions for fitting the mixture of function-on-composition without within-subject correlation
│     ├── auxiliary functions2.R   # auxiliary functions for fitting the mixture of function-on-composition with within-subject correlation
│     ├── comparing models.R       # implementation of the non-mixture and MFoC-i models   
│     ├── data generation.R        # generate the simulated data in case 1
│     ├── existing models.R        # the code for implement the existing models including the mixture of scalar-on composition regressions and mixture of function-on-scalar regressions
│     └── mixfoc-c.R               # implementation of the proposed model
│
├── case2/                         # the case 2 in simulation study 
│     ├── auxiliary functions1.R   # auxiliary functions for the MFoC-i model
│     ├── auxiliary functions2.R   # auxiliary functions for the MFoC-c mdoel
│     ├── comparing models.R       # implementation of the non-mixture and MFoC-i models   
│     ├── data generation.R        # generate the simulated data in case 2
│     ├── existing models.R        # the code for implementing the existing models
│     └── mixfoc-c.R               # implementation of the proposed model
│
└── case3/                         # the case 3 in simulation study 
      ├── auxiliary functions1.R   # auxiliary functions for the MFoC-i model
      ├── auxiliary functions2.R   # auxiliary functions for the MFoC-c mdoel
      ├── comparing models.R       # implementation of the non-mixture and MFoC-i models   
      ├── data generation.R        # generate the simulated data in case 3
      ├── existing models.R        # the code for implementing the MSoC and MFoS models
      └── mixfoc-c.R               # implementation of the MFoC-c model


