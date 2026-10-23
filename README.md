# GEpiSim - *G*enetic *Epi*demiology *Sim*ulator
GEpiSim is an R package for simulating genetic epidemiology dynamics using compartmental models and also simulating breeding programmes accounting for epidemiology. It extends the [AlphaSimR](https://github.com/gaynorr/AlphaSimR) package to simulate infectious disease spread in structured populations with genetic variation.

## Features

 - **Compartmental Models**: Simulate epidemics using SIR, SIS, SEIR, SIDR, and SEIDR models.
 - **Population Dynamics**: Simulate disease progression, transmission, and recovery with genetic effects.
 - **Visualization**: Plot epidemic time series and Kaplan-Meier curves for risk estimation.
 - **Genetic Parameters**: Incorporate genetic variation in susceptibility, infectivity, tolerance, and detection.
 - **Basic Reproduction Number (R0)**: Estimate R0 for infectious diseases.
 - **Customisable Parameters**: Set transmission rates, latent/removal/detection periods, and simulation time.

## Usage
### 1. Initialise simulation parameter
```{r}
library(GEPISim)
# Create founder haplotypes
founderPop <- AlphaSimR::quickHaplo(nInd = 100, nChr = 5, segSites = 1000)

# Set simulation parameters for an SIR model
SP <- SimParamEpidemic$new(
   founderPop,
   model = "SIR",
   r_beta = 0.5,
   removal_period = 10
)
```
### 2. Create a Population with an Epidemic
```{r}
# Add genetic traits (e.g., susceptibility, infectivity, tolerance)
means <- rep(0, 3)
vars  <- rep(1, 3)
SP$addTraitA(10, mean = means, var = vars, name = c("sus", "inf", "tol"))

# Create a population with one initial infected individual per group
pop <- newPopEpidemic(founderPop, group_list = rep(1, 100), indCases = NULL)

# Set phenotypes (heritabilities)
pop <- setPheno(pop, h2 = c(0.3, 0.4, 0.2))
```

### 3. Run the Epidemic Simulation
```{r}
# Simulate the epidemic
pop <- runEpidemic(pop, SP)
```

### 4. Visualize Results
```{r}
# Plot epidemic time series
plotModel(pop, SP)

# Estimate R0
R0 <- getR0(pop, SP)
print(paste("R0:", R0))

# Plot survival curves (Kaplan-Meier)
plotKM(pop, SP)
```

### 5. Implement selection
Use a solver to obtain breeding values for underlying traits (such as [SIRE 2.0](https://github.com/theITEAM/SIRE2.0)) or other traits and select via selection functions in `AlphaSimR`.

## Models Supported

 - SIR: Susceptible → Infectious → Removed
 - SIS: Susceptible → Infectious → Susceptible
 - SEIR: Susceptible → Exposed → Infectious → Removed
 - SIDR: Susceptible → Infectious → Detected → Removed
 - SEIDR: Susceptible → Exposed → Infectious → Detected → Removed

## Future work

We are actively working on the following improvements:

 - **Enhanced Documentation**: Adding more examples, vignettes, and detailed function documentation.
 - **CRAN Compatibility**: Preparing the package for submission to CRAN to ensure long-term maintenance and accessibility.
 - **Performance Optimisation**: Rewriting backend simulation models in C++ using `Rcpp` for faster execution, especially for very large populations.
