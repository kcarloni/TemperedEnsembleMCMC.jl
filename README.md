An (fast) implementation of parallel-tempered, ensemble MCMC.

At the moment, the following are implemented:

    Proposal Moves:
        - Ensemble "Stretch" Move (GW10)

    Parallel Tempering Schemes:
        - Fixed ladder of exponentially scaling temps, 
	  T = [1, dT^1, dT^2, ...]
          for user-input ladder size (ntemps) and dT (dtemps).
        - Control of the overall likelihood of proposing 
	  temperature swaps at any step via
          "prob_propose_swap"; default = 0.1. 

    Automatic Multi-Threading using Threads.@threads

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

julia> x, accepted, nswaps_accepted, nswaps_proposed, log_probs, obs = out

# reshape the output "x" to (nsamples, ndim)
julia> samples = reduce(vcat, reduce(vcat, x)')

# plot!
julia> using Plots
julia> histogram2d(
	samples[:,1], samples[:,2],
	show_empty=true, bins=200, c=cgrad(:blues)
)
```


References

Goodman & Weare (2010) 


