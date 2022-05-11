

A (fast) implementation of parallel-tempered, ensemble MCMC.
Currently unregistered.

At the moment, the following are implemented:

    Sampling values of additional observables:
	- These should be packaged into a second output of 
	  the user's 'log_pdf' function, ie. 
	  'log_pdf(x)::(lp::Float64, obs::Vector{Float64})'

    Proposal Moves:
        - Ensemble "Stretch" Move (GW10)

    Parallel Tempering Schemes:
        - Fixed ladder of exponentially scaling temps, 
	  T = [1, dT^1, dT^2, ...]
          for user-input ladder size (`ntemps`) and dT (`dtemps`).
        - Control of the overall likelihood of proposing
          temperature swaps at any step via
          `prob_propose_swap`; default = 0.1. 

    Automatic Multi-Threading using Threads.@threads

While these other options are not yet:

	- Reproducible randomness
	  (reproducible Julia RNGs are not thread-safe)
	- Adaptive temperature scaling 

For example, to sample from Rosenbrock's banana function (2d), 

```julia
julia> using TemperedEnsembleMCMC
julia> log_pdf(x) = -1/20 * (100(x[2]-x[1]^2)^2 + (1-x[1])^2)

julia> nsteps = 50
julia> nwalkers = 300
# nsamples = nsteps * nwalkers
julia> init_samples = [ randn(2) for i=1:nwalkers ]

# split walkers into T=0 and T=1.41 "rungs"
julia> ntemps = 2
julia> dtemp = 1.41

# optionally, add burn-in steps which are not recorded
julia> nburnin_steps = 20

# also, note that our log_pdf function does not 
# produce any additional observables,
julia> save_observables = false

julia> out = TemperedEnsembleMCMC.run(
	log_pdf, samples_0;
	nsteps=nsteps, nwalkers=nwalkers, 
	nburnin_steps=20
	ntemps=2, dtemp=dtemp,
	save_observables=false)

julia> x, accepted, nswaps_accepted, nswaps_proposed, 
	log_probs, obs = out

# reshape the output "x" to (nsamples, ndim)
julia> samples = reduce(vcat, reduce(vcat, x)')

# plot!
julia> using Plots
julia> histogram2d(
	samples[:,1], samples[:,2],
	show_empty=true, bins=200, c=cgrad(:blues)
)
```

We can benchmark serial performance:

```julia
julia> using BenchmarkTools
julia> @benchmark TemperedEnsembleMCMC.run(
	log_pdf, init_samples;
	nsteps=nsteps, nwalkers=nwalkers, ntemps=2
)

BenchmarkTools.Trial: 7338 samples with 1 evaluation.
 Range (min … max):  632.351 μs …   3.058 ms  ┊ GC (min … max): 0.00% … 73.15%
 Time  (median):     647.424 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   677.091 μs ± 217.372 μs  ┊ GC (mean ± σ):  3.99% ±  9.00%

  █▃                                                            ▁
  ██▇▆▆▁▁▃▁▄▁█▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▃▃▄▄▆▇ █
  632 μs        Histogram: log(frequency) by time       2.33 ms <

```

and parallel:



References

Goodman & Weare (2010) https://msp.org/camcos/2010/5-1/camcos-v5-n1-p04-s.pdf

Foreman-Mackey et. al (2011) https://arxiv.org/abs/1202.3665
emcee: https://github.com/dfm/emcee

