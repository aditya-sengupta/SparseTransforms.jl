"""
Methods for the reconstruction engine; specifically, to

1. carry out singleton detection
2. get the cardinalities of all bins in a subsampling group (debugging only).
"""

include("query.jl")

"""
Noiseless-case singleton detection. Assumes P = n + 1 and D = [0; I].
See singleton_detection for full signature.
"""
function singleton_detection_noiseless(U_slice; kwargs...)
    return -sign.(U_slice * U_slice[1])[2:length(U_slice)], true
end

"""
Singleton detection that uses MLE: looks at the residuals created by peeling off each possible singleton.
See singleton_detection for full signature.
"""
function singleton_detection_mle(U_slice; kwargs...)
    selection, S_slice, n = kwargs[:selection], kwargs[:S_slice], kwargs[:n]
    P = size(S_slice, 1)
    alphas = (1 / P) * S_slice' ⋅ U_slice
    residuals = norm.(U_slice - (alphas * S_slice)') # this is sketch
    k_sel = argmin(residuals)
    return dec_to_bin(selection[k_sel], n), sign(alphas[k_sel])
end

"""
Non-sample-optimal singleton detection.
"""
function singleton_detection_nso(U_slice; kwargs...)
    n = kwargs[:n]
    chunks = reshape(U_slice, (length(U_slice) / (n + 1), n + 1)) |> sign_spright
    chunks = (@pipe (chunks' + chunks[:,1])' |> mod.(_, 2) |> Bool)[:, 2:n]
    choices = vstack((sum(chunks, dims=1), sum(flip.(chunks), dims=1)))
    nso_k = argmin(choices, dims=1)
    return nso_k, 1
end

singleton_detection_lookup = Dict(
    :mle => singleton_detection_mle,
    :noiseless => singleton_detection_noiseless,
    :nso => singleton_detection_nso
)

"""
Finds the true index of a singleton, or the best-approximation singleton of a multiton.
Arguments
---------
U_slice: P-element Array{Float64,1}
The WHT component of a subsampled bin, with element i corresponding to delay i.

Returns
-------
k : n-element Array{Bool,1}
Index of the corresponding right node.
"""
function singleton_detection(U_slice, method=:mle; kwargs...)
    return singleton_detection_lookup.get(method)(U_slice; kwargs)
end