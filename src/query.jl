include("utils.jl")
include("input_signal.jl")

get_b_lookup = Dict(
    :simple => get_b_simple
)

get_Ms_lookup = Dict(
    :simple => get_Ms_simple
)

function get_b_simple(signal::InputSignal)
    return signal.n / 2
end

function get_b(signal::InputSignal, method=:simple)
    return get_b_lookup.get(method)(signal)
end

function get_Ms_simple(n::Int64, b::Int64, num_to_get::Int64 = nothing)
    """
    A semi-arbitrary fixed choice of the subsampling matrices.
    """
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

function get_b(signal::InputSignal, method=:simple)
    return get_Ms_lookup.get(method)(signal)
end

function get_D_standard(n::Int64, kwargs...)
    """
    Gets the delays matrix [0; I], of dimension (n+1, n). See get_D for full signature.
    """
    return map(x -> Bool(x), vcat(zeros(1, n), diagm(repeat([Bool(1)], n))))
end

function get_D_random(n::Int64, kwargs...)
end