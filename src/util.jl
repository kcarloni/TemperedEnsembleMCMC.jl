
# functions to initialize output vectors.
# useful for handling the "don't save" = nothing case.

function init_output_vector(v0s, nsamples)
    vs = Vector{typeof(v0s)}(undef, nsamples)
end
init_output_vector(v0s::Nothing, nsamples) = nothing

# nouter = nsteps, ninner = nwalkers
init_output_vector_of_vector(
    v0s, nouter, ninner, 
    init_fn=init_output_vector) = [init_fn(v0s[1], ninner) for i=1:nouter]

init_output_vector_of_vector(
    v0s::Vector{Nothing}, nouter, ninner,
    init_fn=init_output_vector) = nothing

# convience swap function 
function swap_at!(vec, i1::Int, i2::Int)
    vec[i1], vec[i2] = vec[i2], vec[i1]
end

# # convenience save function
function save_copy!(vv, v, t)
    @inbounds vv[t] = copy(v)
end
# save_copy!(vv, v::Vector{nothing}, t) = nothing

# # convenience save function
# function save_at!(vv, v, w, t)
#     @inbounds vv[w][t] .= v[w]
# end