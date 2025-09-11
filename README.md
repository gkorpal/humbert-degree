# Refined Humbert Invariants in Supersingular Isogeny Degree Analysis

This directory contains the implementation of the algorithms and the results of the computations described in the paper.

**ePrint:** [2025/1605](https://eprint.iacr.org/2025/1605)

**Authors:** Eda Kırımlı and Gaurish Korpal

## Directory Organization

### A. Package

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)

This is a local package, not intended for the [General Registry](https://github.com/JuliaRegistries/General). Therefore, it lacks `git` configuration. It was created using [pkg](https://pkgdocs.julialang.org/v1/creating-packages/) (`pkg> generate Humbert`, `pkg> activate ./Humbert`, and `pkg> add Oscar@1.2.2`).

```
Humbert/
├── Project.toml
└── src/
    ├── Humbert.jl
    ├── main.jl
    └── utils.jl
```

#### A1. `Project.toml`

This file contains the name of the package, its unique UUID, its version, the authors and potential dependencies.

#### A2. `src/`

This directory contains all the functions that can be called to run various experiments. 

1. `Humbert.jl`: This file defines the [module](https://docs.julialang.org/en/v1/base/base/#Reflection) called `Humbert`. It is the heart of our package, which gives access to the following functions:

    ```julia
    julia> println(names(Humbert))
    [:Humbert, :RHI, :allRHI, :classNumbers, :degForm, :fileParser, :getEll, :minDeg, :polyForm, :polz, :polzSQI, :primeList]
    ```

    The documentation of each function is accessible using the [REPL help mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Help-mode).

2. `main.jl`: This file contains an implementation of algorithms `polz`, ,`polzSQI`, `RHI`, `degForm`, `allRHI`, and `minDeg` introduced in our paper. Here `allRHI(p, a, b)` generates a file `RHI_p_a_b.txt` containing all unique refined Humbert invariants and degree forms, along with a polarization $\theta = \begin{pmatrix}m & \alpha \\ \overline{\alpha} & m\end{pmatrix}$ where $a\leq m \leq b$. On the other hand, `minDeg(p, a, b)` parses the file `RHI_p_a_b.txt` and constructs a minimum degree frequency distribution dictionary.

3. `utils.jl`: This file contains helper functions `primeList`, `classNumbers`, `getEll`, `polyForm`, and `fileParser` which are useful for running the experiments.

### B. Scripts

This directory contains some example scripts that can be used for experimentation. 

```
expt/
├── getDeg.jl
├── getRHIs.jl
└── getHh.jl
```

#### B1. `getDeg.jl`

It runs `minDeg(p, ell, 2*ell*p, 2*ell*p+p)` for $n$ primes $p\equiv m\pmod 4$ greater than or equal to $N$. The outputs are saved in a file called `deg_N_lastp.txt`.

#### B2. `getRHIs.jl`

It contains a [parallelizable](https://docs.julialang.org/en/v1/manual/multi-threading/#The-@threads-Macro) `for` loop that runs `allRHI(p, ell, 2*ell*p, 2*ell*p+p)` for $n$ primes $p$ such that $p\geq N$ and $p\equiv m\pmod 4$.

#### B3. `getHh.jl`

It computes the class numbers $h = H_1(p,1)$ and $H=H_2(p,1)$ for any prime $p$.

### C. Data

The following data was generated using (different versions of) above scripts.

```
data/
├── RHI_1_4p_3mod4/
│   └── RHI_p_1_4*p.txt
├── RHI_1_p2_3mod4/
│   └── RHI_p_1_p^2.txt
├── RHI_2lp_2lp+p_1mod4/
│   └── RHI_p_2*l*p_(2*l+1)*p.txt
├── RHI_2lp_3lp_1mod4/
│   └── RHI_p_2*l*p_3*l*p.txt
├── RHI_2p_3p_1mod4/
│   └── RHI_p_2*p_3*p.txt
├── RHI_2p_3p_3mod4/
│   └── RHI_p_2*p_3*p.txt
├── RHI_5ary_deg_4ary/
│   └── RHI_p_ell_2*ell*p_2*ell*p+p.txt
├── deg/
|   ├── deg_13_197.txt
|   ├── deg_19_599.txt
│   └── deg_607_727.txt
└── class_numbers/
    └── class_numbers_19_64_3.txt
```

#### C1. `RHI_1_4p_3mod4/` and `RHI_1_p2_3mod4/`

These directories contain `allRHI(p, 1, 1, 4*p)` and `allRHI(p, 1, 1, p^2)` output to support our choice of searching between $2p$ and $3p$. Note that, for $p = 251$, we get the last new RHI exactly for $m=3\cdot 251=753$. Also, for $p=47,71$ the last $m$ is quite close to $3p$.

*Ignore the first line that says "Theory suggests that we should have `polzCount` quadratic forms."*

#### C2. `RHI_2p_3p_3mod4/`

This directory contains outputs of `allRHI(p, 1, 2p, 3p)` for $p\equiv 3 \pmod 4$ and $19\leq p \leq 727$ from `getRHIs.jl` script. 

*Ignore the first line that says "Theory suggests that we should have `polzCount` quadratic forms."*

#### C3. `RHI_2lp_3lp_1mod4/` and `RHI_2p_3p_1mod4/`

These directories contain `allRHI(p, ell, 2*ell*p, 3*ell*p)` and `allRHI(p, ell, 2*p, 3*p)` output to support our choice of searching between $2 \ell p$ and $2 \ell p + p$. Note that, for $p = 17$, we get the last new RHI exactly for $m=2\cdot 3 \cdot 17 + 17 = 119$. 

#### C4.  `RHI_2lp_2lp+p_1mod4/`

This directory contains outputs of `allRHI(p, ell, 2*ell*p, 2*ell*p+1)` for $p\equiv 1 \pmod 4$ and $13\leq p \leq 200$ from `getRHIs.jl` script. 

#### C5. `RHI_5ary_get_4ary`

In addition to `allRHI(p, ell, 2*ell*p, 2*ell*p+1)` for $13 \leq p \leq 59$ such that for each 5-ary form the automorphism group size and number of ways represents 4 are also recorded. 

Additionally, it contains `Representatives of 5-ary genus rep 1` that is the number of isometry class representatives of the genus of $x_0^2+4(x_1^2 + \ell x_2^2 + px_3^2 + \ell px_4^2)$ which represent 1.

```julia
# Representatives of 5-ary genus rep 1
function Genus5(p::Int, ell::Int)
    A = ZZ[1 0 0 0 0; 0 4 0 0 0; 0 0 4*ell 0 0; 0 0 0 4*p 0; 0 0 0 0 4*ell*p]
    L = integer_lattice(; gram = A)
    Llist = genus_representatives(L)
    count = 0
    for L in Llist
        if minimum(L) == 1
            count += 1
        end        
    end
    return count
end
```

It also contains `Representatives of 4-ary norm genus` that is the class number of genus of $x_1^2 + \ell x_2^2 + px_3^2 + \ell px_4^2$.

```julia
# Representatives of 4-ary norm genus
function Genus4(p::Int, ell::Int)
    A = ZZ[1 0 0 0; 0 ell 0 0; 0 0 p 0; 0 0 0 ell*p]
    L = integer_lattice(; gram = A)
    Llist = genus_representatives(L)
    return length(Llist)
end
```

#### C6. `deg/`

This directory contains the output of `getDeg(19, 54, 3)`, `getDeg(607, 10, 3)`, and `getDeg(13, 20, 1)` from the `getDeg.jl` script.

#### C7. `class_numbers/`

This directory contains output of function `getHh(19,64,3)` from the `getHh.jl` script.

## Recommended Workflow

To experiment with this code, we can follow [the recommended development workflow](https://modernjuliaworkflows.org/writing/#development_workflow). 

0. Clone the repository, and enter inside the directory. We will assume that our [Julia environment](https://pkgdocs.julialang.org/v1/environments/) has the [Oscar package](https://github.com/oscar-system/Oscar.jl/?tab=readme-ov-file#installation) installed.

1. Use `Pkg.develop` to add the functions to our Julia environment.
    ```julia 
    julia> using Pkg

    julia> Pkg.develop(path="./Humbert")
    ```
    This will generate a file called `Manifest.toml` which contains an absolute record of the state of the packages in the `./Humbert` environment. It includes exact information about (direct and indirect) dependencies of the project. It is not tracked by Git.

2.  Now we can use any of our function by using the local `Humbert` package.

    ```julia

    julia> using Humbert

    julia> @time allRHI(19, 1, 2*19, 3*19)
    working with prime 19
    saved data for prime 19
      1.550728 seconds (6.31 M allocations: 350.799 MiB, 14.78% gc time, 50.12% compilation time)

    julia> @time minDeg(19, 1, 2*19, 3*19)
      0.614780 seconds (1.52 M allocations: 77.198 MiB, 98.32% compilation time: 16% of which was recompilation)
    4

    julia> include("./expt/getRHIs.jl")
    working with prime 19
    saved data for prime 19
    working with prime 23
    saved data for prime 23

    julia> include("./expt/getDeg.jl")
    ```

All code have been formatted in [BlueStyle](https://domluna.github.io/JuliaFormatter.jl/dev/blue_style/).

## System Requirements

The code was written and tested on a system with the following specifications:

CPU: 12th Gen Intel i7-12800HX (24)

Memory: 15843MiB

OS: Ubuntu 22.04.5 LTS on Windows 11 x86_64 (WSL)

Software: Julia Version 1.11.3 with Oscar Version 1.2.2

*Most of the data was generated using the HPC facilities of the Advanced Computing Research Centre, University of Bristol.*