# SolarDynamo.jl 🌞

[![Build Status](https://github.com/Eawag-SIAM/SolarDynamo.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Eawag-SIAM/SolarDynamo.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Eawag-SIAM/SolarDynamo.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Eawag-SIAM/SolarDynamo.jl)


Stochastic simulation of the number of sunspots.


## Installation


```Julia
] add git@github.com:Eawag-SIAM/SolarDynamo.jl.git
```

## Usage


```julia
using SolarDynamo

θ = [2.0,   # τ
     3,     # T
     8.5,   # Nd
     0.1,   # sigma
     10]    # Bmax

# -- run model

sn(θ, Tobs = 929, Twarmup = 200)

# You can provide a random seed to the solver
sn(θ, Tobs = 929, Twarmup = 200, seed=314)


# -- Compute summary statics based on Fourier components

res = sn(θ, Tobs = 929)
summary_statistics(res)


# The window weights can also be precomputed:
Tobs = 929
window = hann_window(Tobs)
summary_statistics(res, window)
```

## Reference

Ulzega, S., Beer, J., Ferriz-Mas, A., Dirmeier, S., & Albert, C. (2025). Shedding light on the solar dynamo using data-driven Bayesian parameter inference. The Astrophysical Journal, 992(1), 61.
