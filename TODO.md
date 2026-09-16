## High Priority
- [ ] Add model plots for all models
    - [x] SIR
    - [ ] SIDR
    - [ ] SEIR
    - [ ] SEIDR
- [x] `addTrait` currently assumes all traits are of the same type. This forces some behaviour that is needed (naming) but pleitropic with same QTL and same trait type is also a strong assumption. Think of a way to keep the naming constraint but remove the other
    - [x] SIR
    - [ ] SIDR
    - [ ] SEIR
    - [ ] SEIDR
- [ ] Include two vignettes
- [ ] Include implicit conversion in `run`
    - [ ] Rename `run` to `runEpidemic`
- [ ] add versions

## Medium Priority
- [ ] Add SIS and R0 for SIS
- [ ] Remove all warnings from devtool::build()
    - [ ] pheno traits (as strings) in `runEpidemic` to be be replaced with proper placeholders
        - [ ] SIR
        - [ ] SIDR
        - [ ] SEIR
        - [ ] SEIDR
    - [ ] timing names in `runEpidemic` to be replace with placeholders
        - [ ] SIR
        - [ ] SIDR
        - [ ] SEIR
        - [ ] SEIDR
- [x] Include `addEpiTrait2` for "other" functions in AlphaSimR that don't start with `add` - NO LONGER NEEDED

## Low Priority
- [ ] Write model in C++
