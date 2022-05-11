An (fast) implementation of parallel-tempered, ensemble MCMC.

At the moment, the following are implemented:

    Proposal Moves:
        - Ensemble "Stretch" Move (GW10)

    Parallel Tempering Schemes:
        - Ladder of exponentially scaling temps, T = [1, dT^1, dT^2, ...]
            for user-input ladder size (ntemps) and dT (dtemps).
        - Control the overall likelihood of proposing temperature swaps at any step
            via "prob_propose_swap"; default = 0.1. 

    Automatic Multi-Threading using Threads.@threads

For example, to sample from Rosenbrock's banana function (2d), 

```julia
using TemperedEnsembleMCMC
log_pdf(x) = -1/20 * (100(x[2]-x[1]^2)^2 + (1-x[1])^2)

nsteps = 50
nwalkers = 300
# nsamples = nsteps * nwalkers
init_samples = [ randn(2) for i=1:nwalkers ]

# split walkers into T=0 and T=1.41 "rungs"
ntemps = 2
dtemp = 1.41

# optionally, add burn-in steps which are not recorded
nburnin_steps = 20

# also, note that our log_pdf function does not produce any additional observables,
save_observables = false

x, accepted, n_swaps_accepted, n_swaps_proposed, log_probs, obs = TemperedEnsembleMCMC.run(
	log_pdf, samples_0;
	nsteps=nsteps, nwalkers=nwalkers, 
	nburnin_steps=20
	ntemps=2, dtemp=dtemp,
	save_observables=false)

# reshape the output "x" to (nsamples, ndim)
samples = reduce(vcat, reduce(vcat, x)')

# plot!
using Plots
histogram2d(
	samples[:,1], samples[:,2],
	show_empty=true, bins=200, c=cgrad(:blues)
)
```


References
GW10

