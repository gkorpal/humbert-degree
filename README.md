# Refined Humbert Invariants in Supersingular Isogeny Degree Analysis

This repository contains the code and data accompanying the paper.

**ePrint:** [2025/1605](https://eprint.iacr.org/2025/1605)

**Authors:** Eda Kırımlı and Gaurish Korpal

## Scripts

All scripts are in [`/scripts/`](/scripts/) and are implemented in Julia using the [Oscar](https://www.oscar-system.org/) package and its dependency [Nemo/Hecke](https://thofma.github.io/Hecke.jl/stable/). Each script is self-documented with usage instructions in its header comments. The scripts form two pipelines.

**Pipeline 1: Polarizations and reducible RHIs**

- [`polz.jl`](/scripts/polz.jl) implements $\mathsf{dynamicPolz}(p)$ (Algorithm 2), enumerating all principal polarizations on $E \times E$ up to equivalence for a given prime $p \equiv 11 \pmod{12}$. Each polarization is a unimodular positive definite binary Hermitian form $\theta$ over the fixed maximal order $\mathcal{O} \subset B_p$. Following the approach of [Greenberg–Voight](http://arxiv.org/abs/1209.2460v1) and [Chisholm](http://hdl.handle.net/11023/1920), each polarization is converted to a rank-8 integral lattice equipped with four auxiliary trace forms $F_1, F_2, F_3, F_4$. The auxiliary forms are LLL-reduced and two candidates are identified as equivalent via exact simultaneous Plesken–Souvignier isometry. The enumeration bound grows dynamically until the number of non-isometric polarizations reaches the class number $\mathbf{H}(p)$. This is the most computationally intensive step.

- [`forms.jl`](/scripts/forms.jl) implements Algorithms 1 and 3. It reads the polarizations produced by `polz.jl` and computes, for each polarization $\theta$, the coefficient matrix of the reducible refined Humbert invariant $q_{(E \times E, \theta)}$ and the corresponding degree form $\deg$, where
$$q_{(E \times E,\, \theta)}(D) = X^2 + 4 \cdot \deg.$$

- [`deg.jl`](/scripts/deg.jl) reads the output of `forms.jl` and computes the distribution of $d = \max_{E_1, E_2} \min\{N : q_{E_1,E_2} = N\}$ over all computed isometry classes of reducible refined Humbert invariants.

**Pipeline 2: Genus of $q_4$ and degree maps**
- [`genus.jl`](/scripts/genus.jl) computes the genus of the lattice
$$q_4 = t_1^2 + t_1t_3 + t_2^2 + t_2t_4 + \tfrac{p+1}{4}t_3^2 + \tfrac{p+1}{4}t_4^2,$$
which is the degree form extracted from the canonical refined Humbert invariant
$$q_5 = t_0^2 + 4t_1^2 + 4t_1t_3 + 4t_2^2 + 4t_2t_4 + (p+1)t_3^2 + (p+1)t_4^2$$
corresponding to the identity polarization $\theta_0$. The genus of $q_4$ counts all isometry classes of reducible refined Humbert invariants.

- [`deg_gen.jl`](/scripts/deg_gen.jl) computes the $d$ value for each genus class of $q_4$, giving the complete distribution of minimum isogeny degrees.

## Data

All data is in [`/data/`](/data/) as human-readable `.txt` files, organized by computation type.

- [`/data/polarizations/`](/data/polarizations/) contains the principal polarizations produced by `polz.jl` for each prime $p$. Each file lists the parameters $(u, v, w, x, y, z)$ of each polarization $\theta$ with diagonal entries $u, v \in \mathbb{Z}_{>0}$ and off-diagonal entry $\alpha \in \mathcal{O}$ such that $\alpha = w + xi + y\frac{i+j}{2} + z\frac{1+k}{2}$.

- [`/data/RHI1/`](/data/RHI1/) contains the quadratic forms produced by `forms.jl`.

- [`/data/deg/`](/data/deg/) contains the $d$ values produced by `deg.jl` for the experimentally computed isometry classes.

- [`/data/Gen5/`](/data/Gen5/) contains the genus data for $q_4$ produced by `genus.jl`.

- [`/data/deg_gen/`](/data/deg_gen/) contains the $d$ values produced by `deg_gen.jl` for the genus classes of $q_4$.

The data covers 31 primes $p \equiv 11 \pmod{12}$ with $10 < p < 660$. For primes where `polz.jl` did not reach the full class number $\mathbf{H}(p)$ within the time limit, the RHI and degree data reflects only the polarizations found so far.

## System Requirements

The computations in `polz.jl` were run on the HPC facilities of the *Advanced Computing Research Centre, University of Bristol*, with a time limit of two weeks per prime. All other scripts were run on a local machine with the following specifications.

- **CPU:** 12th Gen Intel i7-12800HX (24 cores)
- **Memory:** 15843 MiB
- **OS:** Ubuntu 22.04.5 LTS on Windows 11 x86\_64 (WSL)
- **Julia:** Version 1.12.6
- **Oscar:** Version 1.8.0
