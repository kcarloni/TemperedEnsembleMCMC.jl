# Moves.jl

abstract type Move end

abstract type EnsembleMove <: Move end
abstract type OneRefMove <: EnsembleMove end

abstract type IndependentMove <: Move end

struct GaussianMove <: IndependentMove end

Base.@kwdef struct StretchMove <: OneRefMove end

function get_proposal( # StretchMove 
    proposal_type::Type{T},
    old_sample::Vector{Float64},
    ref_sample::Vector{Float64}, rng, a=2,
    # c_set::Matrix{Float64},
    ) where (T<:StretchMove)

    # get a random ref sample from complementary set of samples
    # ref_sample = rand(c_set)

    # want to sample z ~ g(z) = N * 1/sqrt(z) for z in (1/a, a); 0 else.
    # sample t ~ [0, 1], and transform by
    # z = f(t) = ((a-1)t + 1)^2 / a
    # then g(z) = [0,1](f⁻¹(z)) * | d/dz (f⁻¹(z)) | 
    # normalization N = 1/(2*(sqrt(a)-1/sqrt(a)))
    z = ((a - 1)*rand(rng) + 1)^2 / a

    # stretch move:
    new_sample = @. ref_sample + z*(old_sample - ref_sample)

    # accept_prob = z^(ndims-1) p(new)/p(old)
    # log_accept_prob = (ndims-1) log(z) + logp(new) - logp(old)
    # "factor" = (ndims-1) log z 
    factor = (length(old_sample) - 1) * log(z)
    return new_sample, factor
end