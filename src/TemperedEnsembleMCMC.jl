module TemperedEnsembleMCMC

# main loop 
using Random
using StaticArrays 

include("util.jl")
include("moves.jl")

#observable_names=["",]
function run(
    log_pdf, init_samples;
    nsteps=50, nwalkers=100, seed=123,
    nburnin_steps=0, # not saved!
    ntemps=2, dtemp=1.41, prob_propose_swap=0.1,
    save_observables=false,
)
    @assert (length(init_samples) == nwalkers) "initial samples should have length = nwalkers"
    @assert (nwalkers % 2 == 0) "nwalkers should be even"
    @assert (nwalkers % ntemps == 0) "number of temperatures should evenly divide nwalkers"
    @assert ((nwalkers ÷ ntemps)% 2 == 0) "nwalkers/ntemps should be even"

    whalf = nwalkers ÷ 2
    n_per_temp = nwalkers ÷ ntemps
    prob_propose_swap = (ntemps==1) ? 0 : prob_propose_swap
    # temperature ladder...
    β(r) = dtemp^(-r) # r = 0:ntemps-1

    log_pdf_ = save_observables ? log_pdf : x -> (log_pdf(x), nothing)

    rng = MersenneTwister(seed)

    x0 = copy(init_samples)
    tmp = log_pdf_.(init_samples)
    lp0, o0 = getindex.(tmp, 1), getindex.(tmp, 2) 

    x = init_output_vector_of_vector(x0, nsteps, nwalkers) # save timesteps together
    lp = init_output_vector_of_vector(lp0, nsteps, nwalkers)
    o = init_output_vector_of_vector(o0, nsteps, nwalkers)

    naccept = zeros(Int, nwalkers) 
    nswaps_p = 0
    nswaps_a = 0
    perm = randperm(rng, n_per_temp) # allocate once!
    wsets = SVector(1:whalf, whalf+1:nwalkers) # divide walkers

    for t = 1-nburnin_steps:nsteps # t>0 = save_data
        # propose move...
        propose_move!(
            x0, o0, lp0, naccept, rng, 
            ntemps, n_per_temp, β,
            perm, wsets, whalf, log_pdf_,
            StretchMove
        )

        # propose swap...
        if rand(rng) < prob_propose_swap
            nswaps_p += nwalkers-ntemps # skip r=0
            nswaps_a += propose_swap!(
                (x0, o0, lp0, naccept), rng, 
                ntemps, n_per_temp, β, perm
            )
        end

        if t == 0         # next step save
            naccept .*= 0
            nswaps_p = 0
            nswaps_a = 0
        elseif t > 0        # save to output...
            save_copy!(x,x0,t)
            save_observables && save_copy!(o,o0,t)
            save_copy!(lp,lp0,t)
        end

    end
    accept_ratio = [na/nsteps for na in naccept]
    return x, accept_ratio, nswaps_a, nswaps_p, lp, o   
end

# PARALLEL
function propose_move!(
    x0, o0, lp0, naccept, rng, 
    ntemps, n_per_temp, β, 
    perm, wsets, whalf, log_pdf_,
    proposal_type::Type{pT}
) where {pT <: OneRefMove}
    # shuffle set of walkers at each temp
    shuffle!(rng, perm)
    for i=1:2
        Threads.@threads for w0 in wsets[i]
            # index manips to quickly
            # - shuffle each temperature rung
            # - get a random ref walker in the other half of the rung.

                # w ∈ 1:nwalkers, r ∈ 0:ntemps-1, k ∈ 1:n_per_temp.
                
                # set of all w in rung r (all at same fixed temp) is
                # {w | r(w) = r } = (1+r) + ntemps*(1:n_per_temp - 1)

                # for fixed r,
                # k(w) = 1 + (w - 1 - r) ÷ ntemps
                # w(k) = (1+r) + (k-1)*ntemps

                # - shuffle rung by sending k -> perm[k]

            r = (w0 - 1) % ntemps
            # - shuffle rung by sending k -> perm[k]
            w = (1 + r) + ntemps * (perm[1 + (w0 - r - 1)÷ntemps]-1)

            # - get random k from other half of rung 
            #   by drawing from (1:n_per_temp÷2)-1 and applying transf.
            kref = 1 + (whalf *(i%2) + ntemps*(rand(rng, 1:n_per_temp÷2)-1)) ÷ ntemps
            wref = (1 + r) + ntemps * (perm[kref]-1)

            x1, f = get_proposal(
                proposal_type, x0[w], x0[wref], rng)
            lp1, o1 = log_pdf_(x1)

            if log(rand(rng)) < f + β((w-1)%ntemps) * (lp1-lp0[w])
                @inbounds begin
                    naccept[w] += 1
                    x0[w] = x1
                    lp0[w] = lp1
                    o0[w] = o1
                end
            end         
        end
    end
end



function propose_swap!(state0, rng, ntemps, n_per_temp, β, perm)
    lp0 = state0[3]
    swaps_accepted = 0
    for r in (ntemps-1):-1:1 # skip 0, and go backwards
        shuffle!(rng, perm)
        dβ = β(r-1) - β(r)
        for k in 1:n_per_temp
            # w(r, k) = (1+r) + (k-1)*ntemps
            w_high = (1+r) + (perm[k]-1)*ntemps # w(r, perm[k])
            w_low = r + (k-1)*ntemps # w(r-1), k

            if log(rand(rng)) < dβ * (lp0[w_high] - lp0[w_low])
                swaps_accepted += 1
                for v in state0
                    swap_at!(v, w_high, w_low)
                end
            end
        end
    end
    return swaps_accepted
end

# module ends
end