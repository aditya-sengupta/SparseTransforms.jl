# SparseTransforms.jl

Implementations of and extensions to the SPRIGHT (SParse Robust Iterative Graph-based Hadamard Transform) algorithm in Julia, and potentially the FFAST (Fast Fourier Aliasing-based Sparse Transform) algorithm too.

**Authors and Acknowledgements**

Version 0.1.0 of this package was written by Aditya Sengupta (@aditya-sengupta), Tynan Sigg (@trsigg), Catherine Huang (@thecatherinehuang) and Ben Hoberman (@bhoberman) as the final project for EECS 229A at UC Berkeley in Fall 2020, with advice from Orhan Ocal, Amirali Aghazadeh, and Kannan Ramchandran.

**References**

[1] Li, X., Ramchandran, K. (2015). An Active Learning Framework using Sparse-Graph Codes for Sparse Polynomials and Graph Sketching. NeurIPS, http://papers.neurips.cc/paper/5697-an-active-learning-framework-using-sparse-graph-codes-for-sparse-polynomials-and-graph-sketching.pdf

[2] Li, X., Bradley, J., Pawar, S., Ramchandran, K. (2015). SPRIGHT: A Fast and Robust Framework for Sparse Walsh-Hadamard Transform. arXiv.org. cs.IT. https://arxiv.org/abs/1508.06336

[3] Pawar, S., Ramchandran, K. (2017). FFAST: An Algorithm for Computing an Exactly $k$-Sparse DFT in $O(k \log k)$ Time. IEEE, https://ieeexplore.ieee.org/abstract/document/8022929

Comments in the code referencing equation and page numbers refer to [2] unless otherwise specified.
