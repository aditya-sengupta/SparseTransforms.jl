using SparseTransforms
using NPZ

noise_sd = 0.1
synth_data = npzread("data/1kuh_10_exhaustive.npz")
signal = InputSignal(synth_data["y"], noise_sd)
transformed = spright(signal, [:simple, :identity_like, :mle])
