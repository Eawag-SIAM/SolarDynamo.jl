# SolarDynamo.jl

[![Build Status](https://github.com/Eawag-SIAM/SolarDynamo.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Eawag-SIAM/SolarDynamo.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Eawag-SIAM/SolarDynamo.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Eawag-SIAM/SolarDynamo.jl)


Stochastic simulation of the number of sunspots.


## Installation


```Julia
] add https://github.com/Eawag-SIAM/SolarDynamo.jl
```

## Usage


```julia
using SolarDynamo

θ = [2.0,  # τ
     3,     # T
     8.5,   # Nd
     0.1,   # sigma
     10]    # Bmax

sn(θ, Tobs = 929, Twarmup = 200)


# providing a random seed to the solver
sn(θ, Tobs = 929, Twarmup = 200, seed=314)
```
