# tUPL-to-MPI
Repository for the master's thesis of Dennis Buurman, s2027100, at Leiden University. This project is supervised by Dr. K.F.D. Rietveld.

tUPL [1] is a framework for parallel programming. tUPL specifications are inherently parallel and do not allow explicit data structures to be specified.
Hommelberg et al. [2] defined transformations on tUPL specifications that allow the generation of optimized parallel C/C++ MPI implementations. 
However, this work did not include transformations for communication code resulting in an incomplete derivation process.
Within this project, the goal is to complete this derivation process by adding rigid transformations for the communication code.
Furthermore, we will investigate whether different derivations from the same tUPL specification lead to best performance on different hardware and MPI rank configurations.

[1] K. F. D. Rietveld and H. A. G. Wijshoff, “Forelem: A versatile optimization framework for tuple-based computations,” in CPC 2013: 17th Workshop on Compilers for Parallel Computing, Citeseer, 2013.

[2] A. Hommelberg, B. van Strien, K. F. D. Rietveld, and H. A. G. Wijshoff, “A new framework for expressing, parallelizing and optimizing big data applications,” arXiv preprint arXiv:2203.01081, 2022.
