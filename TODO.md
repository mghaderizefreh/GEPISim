## High Priority
- [x] Add model plots for all models
    - [x] SIR
    - [x] SIDR
    - [x] SEIR
    - [x] SEIDR
- [x] Include viridis package for colors
- [x] `addTrait` currently assumes all traits are of the same type. This forced some behaviour that is needed (naming) but pleiotropic with same QTL and same trait type is also a strong assumption. For naming, `runEpidemic` now checks for different models. Check with Gregor if it is a good idea to check in `runEpidemic` or in `PopEpidemic` constructor. Also, the non-zero mean may be checked with a t.test()$p.value but that can only be done in the constructor and forces the user to define traits
    - [x] SIR
    - [x] SIDR
    - [x] SEIR
    - [x] SEIDR
    - [ ] include t.test for nonzero checking -> consult with GG
- [ ] Include two vignettes
- [x] Include implicit conversion in `run` -> a warning message is produced instead
    - [x] Rename `run` to `runEpidemic`
- [x] Convert rates (removal, latency, detection) to constant values and remove all distribution-implementation (following discussion with Andrea on Sep 4th)
- [x] add versions
- [ ] other plots
    - [x] km plot

## Medium Priority
- [ ] SIS
    - [ ] model
    - [ ] R0
    - [ ] plot
### Remove all warnings from devtool::build()
    - [ ] pheno traits (as strings) in `runEpidemic` to be replaced with proper placeholders
        - [ ] SIR
        - [ ] SIDR
        - [ ] SEIR
        - [ ] SEIDR
    - [ ] timing names in `runEpidemic` to be replaced with placeholders
        - [ ] SIR
        - [ ] SIDR
        - [ ] SEIR
        - [ ] SEIDR
    - [ ] other
- [x] Include `addEpiTrait2` for "other" functions in AlphaSimR that don't start with `add` - NO LONGER NEEDED
- [x] Add unit tests
- [ ] Add fixed effect

## Low Priority
- [ ] Write models in C++
