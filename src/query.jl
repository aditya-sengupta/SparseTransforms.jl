"""
Methods for the query generator: specifically, to

1. generate sparsity coefficients b and subsampling matrices M
2. get the indices of a signal subsample
3. compute a subsampled and delayed Walsh-Hadamard transform.
"""

include("utils.jl")
include("input_signal.jl")
using StatsBase

"""
A semi-arbitrary fixed choice of the sparsity coefficient. 
See get_b for full signature.
"""
function get_b_simple(signal::InputSignal)
    return signal.n / 2
end

get_b_lookup = Dict(
    :simple => get_b_simple
)

"""
Get the sparsity coefficient for the signal.

Arguments
---------
signal : InputSignal
The signal whose WHT we want.

method : string
The method to use. All methods referenced must use this signature (minus "method".)

Returns
-------
b : Int64
The sparsity coefficient.
"""
function get_b(signal::InputSignal, method=:simple)
    return get_b_lookup.get(method)(signal)
end

"""
A semi-arbitrary fixed choice of the subsampling matrices. See get_Ms for full signature.
"""
function get_Ms_simple(n::Int64, b::Int64, num_to_get::Int64 = nothing)
    @assert n % b == 0, "b must be exactly divisible by n"
    if isnothing(num_to_get)
        num_to_get = n / b
    end

    Ms = Any[]
    for i in num_to_get-1:-1:-1
        M = zeros(n, b)
        M[(b * i):(b * (i + 1)), :] = I
        push!(Ms, M)
    end
    return Ms
end

get_Ms_lookup = Dict(
    :simple => get_Ms_simple
)

function get_Ms(signal::InputSignal, method=:simple)
    return get_Ms_lookup.get(method)(signal)
end

"""
Gets the delays matrix [0; I], of dimension (n+1, n). See get_D for full signature.
"""
function get_D_standard(n::Int64, kwargs...)
    return Bool.(vcat(zeros(1, n), Matrix{Bool}(I, n, n)))
end

"""
Gets a random delays matrix of dimension (num_delays, n). See get_D for full signature.
"""
function get_D_random(n::Int64; kwargs...)
    num_delays = kwargs[:num_delays]
    choices = sample(1:2^n, num_delays, replace=false)
    return map(c -> dec_to_bin(c, n), choices)
end

"""
Get a repetition code based (NSO-SPRIGHT) delays matrix. See get_D for full signature.
"""
function get_D_nso(n::Int64; kwargs...)
    num_delays = kwargs[:num_delays]
    p1 = num_delays / n
    random_offsets = get_D_random(n; num_delays=p1)
    D = Array{Bool}(undef, 0, n)
    identity_like = get_D_identity_like(n)
    for row in random_offsets
        modulated_offsets = mod(row + identity_like, 2)
        D = vstack(D, modulated_offsets)
    end
end

get_D_lookup = Dict(
    :standard => get_D_standard,
    :random => get_D_random
)

"""
Delay generator: gets a delays matrix.

Arguments
---------
n : Int64
number of bits: log2 of the signal length.

Returns
-------
D : num_delays×n Array{Bool,2}
The delays matrix; if num_delays is not specified in kwargs, see the relevant sub-function for a default.
"""
function get_D(n::Int64, method=:standard; kwargs...)
    return get_D_lookup.get(method)(n)
end

"""
Query generator: creates indices for signal subsamples.
Arguments
---------
M : n×b Array{Bool,2}
The subsampling matrix.

d : n-element Array{Bool,1}
The subsampling offset.

Returns
-------
indices : B-element Array{Int64,1}
The (decimal) subsample indices.
"""
function subsample_indices(M, d)
    L = binary_ints(size(M, 2))
    inds_binary = mod(M⋅L + d, 2)
    return bin_to_dec(inds_binary)
end

"""
Creates random delays, subsamples according to M and the random delays,
and returns the subsample WHT along with the delays.

Arguments
---------
signal : InputSignal
The signal to subsample, delay, and compute the WHT of.

M : n×b Array{Bool,2}
The subsampling matrix.

num_delays : Int64
The number of delays to apply; or, the number of rows in the delays matrix.

force_identity_like : Bool
Whether to make D = [0; I] like in the noiseless case; for debugging.
"""
function compute_delayed_wht(signal::InputSignal, M, D)
    inds = map(d -> subsample_indices(M, d), D)
    used_inds = Set(inds)
    samples_to_transform = signal.signal_t[inds]
    return fwht.(samples_to_transform), used_inds
end