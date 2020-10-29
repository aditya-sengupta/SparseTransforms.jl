# SPRIGHT.jl

An implementation of and extension to the SPRIGHT (SParse Robust Iterative Graph-based Hadamard Transform) algorithm, in Julia.

**Setup**

In order to start using SPRIGHT, the Julia programming language is required, ideally with the Revise package. The specific steps to set this up are:

1. `git clone` the repo to wherever you want
2. Download Julia from julialang.org
3. Run `julia` from terminal, type `]`, and run `add Revise`
4. Backspace and run `using Revise` (the purpose of this is just to have it precompile once)
5. Go to `~/.julia/config/startup.jl` (or the equivalent for a different OS) and add `using Revise` there
6. `cd` to wherever you cloned the repo and run `julia`
7. Within the REPL type `]` and run `activate .`
8. Hit backspace and type using SPRIGHT  - this precompile will take a while but you'll only have to do it once (you may also get warnings about "replacing docs" - I'm working on trying to fix those)

You can now access the `transform` function, which carries out the sparse W-H transform of a signal; this should be the only interface needed for users of SPRIGHT.

**References**

[1] Li, X., Ramchandran, K. (2015). An Active Learning Framework using Sparse-Graph Codes for Sparse Polynomials and Graph Sketching. NeurIPS, http://papers.neurips.cc/paper/5697-an-active-learning-framework-using-sparse-graph-codes-for-sparse-polynomials-and-graph-sketching.pdf

[2] Li, X., Bradley, J., Pawar, S., Ramchandran, K. (2015). SPRIGHT: A Fast and Robust Framework for Sparse Walsh-Hadamard Transform. arXiv.org. cs.IT. https://arxiv.org/abs/1508.06336

Comments in the code referencing equation and page numbers refer to [2] unless otherwise specified.