module SolarDynamo


"""
Include the following lines in the main script:
---
include("./SDDESolarDynamo.jl")
using .SDDESolarDynamo
---
This gives access to functions 'bfield' and 'sn'

See here for the SDDE interface: https://docs.sciml.ai/DiffEqDocs/stable/tutorials/dde_example/
"""

export bfield, sn

using StochasticDelayDiffEq
using SpecialFunctions: erf
using StaticArrays

# --- Nonlinear function
ftilde(x, Bmin, Bmax) = x/4 * (1 + erf(x^2-Bmin^2)) * (1 - erf(x^2-Bmax^2))


# ---------------------------------
# B-field WITHOUT Jupiter
# ---------------------------------

# --- Model : B field

function f(u,h,p,t)     # Drift function
    #  u = [B, dB/dt]
    # du = [dB/dt, d^2B/dt^2]
    τ, T, Nd, sigma, Bmax = p
    # --- with Jupiter
    hist = h(p, t - T, idxs = 1)    # B[1](t-T)
    du1 = u[2]
    du2 = -u[1]/τ^2 - 2*u[2]/τ - Nd/τ^2*ftilde(hist, 1, Bmax)
    SA[du1, du2]
end

function g(u,h,p,t)     # Diffusion function
    τ, T, Nd, sigma, Bmax = p
    du1 = 0
    du2 = Bmax*sigma / (τ^(3/2))
    SA[du1, du2]
end


function bfield(θ, Tsim; kwargs...)

    τ, T, Nd, sigma, Bmax = θ

    # --- define initial values
    u0 = SA[Bmax, 0.0]
    # define inital values [B(t), dB(t)/dt] if t < t0
    # h(p, t) = (Bmax, 0.0)
    # we can speed things up by providing a call by index, see example linked above:
    h(p, t; idxs = nothing) = idxs == 1 ? Bmax : (Bmax, 0.0)

    # use constant lags
    lags = (T, )
    tspan = (0.0, Tsim)

    prob = SDDEProblem(f, g, u0, h, tspan, θ; constant_lags = lags)
    solve(prob, EM(); dt=0.1, saveat=1.0, kwargs...)
end



# --- SN from B field

# Data-generating model
function sn(θ; Twarmup = 200, Tobs = 929, kwargs...)

    Tsim = Twarmup + Tobs  # Total simulation steps

    sol = bfield(θ, Tsim; kwargs...)

    # square result and get rid of warmup points
    y = map(abs2, sol[1, (Twarmup + 2):end])

    return y
end


end
