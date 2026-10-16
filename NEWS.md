# GEPISim  1.0.0

* indexing added for PopEpidemic

* SimParamEpidemic to show a message on creation

* Removed epi-helper, group_inf from dynamics

* get_R0 renamed to getR0 and added option parameter aggfunc with default `median`

* To-do task
- [ ] Write models in C++
- [ ] Add fixed effect
- [ ] Include a vignettes with cattle data and SIS model
- [ ] Include a vignette with fish genome and CMS (SIR model)
- [ ] custom trait names for all model and rectify devtools::check()

# GEPISim 0.15.1.9002

* appropriate names for timings (now tdeath for SIS)

* model info is removed from PopEpidemic object (wouldn't work for SIS)

* Add model plots for all models

* Include viridis package for colors

* Informative message for when `runEpidemic` is called without approrpriate types

* Rates (removal, latency, detection) are now constant

* Version tracking

* NEWS tracking

* Models added SIR SIS SEIR SIDR SEIDR

* R0 calculator

* plotModel

* KM plots for all but SIS models

* Add unit tests
