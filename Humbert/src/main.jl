"""
    polz(p, ell, m) and polzSQI(p, m)

Given prime p and m > 0 such that (m^2 - 1) is not divisible by p, this gives 
polarizations of the form 
[u       h]
[bar{h} v]
where u = v = m and 
h = w + x*i + y*j + z*k with 
nrd(h) = m^2 - 1 and  0 <= w <= x <= m-1 and 0 <= y <= z <= m-1.

# Examples
```jldoctest
julia> @time polzSQI(19, 7)
  0.000007 seconds (6 allocations: 352 bytes)
2-element Vector{Vector{Int64}}:
 [7, 7, 1, 3, 1, 1]
 [7, 7, 2, 5, 0, 1]

# References
[1] K. Hashimoto and T. Ibukiyama, On class numbers of positive definite binary quaternion Hermitian forms, J. Fac. Sci. Univ. Tokyo Sect. IA Math. **27** (1980), no. 3, 549--601; MR0603952
[2] N. Grieve, Reduced norms and the Riemann-Roch theorem for Abelian varieties, New York J. Math. **23** (2017), 1087--1110; MR3711271
[3] T. Ibukiyama, Supersingular abelian varieties and quaternion hermitian lattices, in *Theory and Applications of Supersingular Curves and Supersingular Abelian Varieties*, 17--37, RIMS Kokyuroku Bessatsu, B90, Res. Inst. Math. Sci. (RIMS), Kyoto, ; MR4521511
```
"""

function polz(p::Integer, ell::Integer, m::Integer)
    sols = Vector{Vector{Int}}()

    needed = m^2 - 1
    B = m - 1
    # Outermost loop: iterate over z from 0 to B.
    for z in 0:B
        # Early exit if minimal contribution p*z^2 already exceeds needed.
        if ell * p * z^2 > needed
            break
        end
        # Loop over y from 0 to B // with y ≤ z.
        for y in 0:B
            part_yz = p * (y^2 + ell * z^2)
            if part_yz > needed
                break
            end
            rem = needed - part_yz
            # Loop over x from 0 to B.
            for x in 0:B
                if ell * x^2 > rem
                    break
                end
                candidate = rem - ell * x^2
                # Directly solve for w: there is a valid w≤x when candidate is a perfect square.
                w = isqrt(candidate)
                if w^2 == candidate && w <= B
                    push!(sols, [m, m, w, x, y, z])
                end
            end
        end
    end
    return sols
end

function polzSQI(p::Integer, m::Integer)
    sols = Vector{Vector{Int}}()
    needed = m^2 - 1
    B = m - 1
    # Outermost loop: iterate over z from 0 to B.
    for z in 0:B
        # Early exit if minimal contribution p*z^2 already exceeds needed.
        if p * z^2 > needed
            break
        end
        # Loop over y with y ≤ z.
        for y in 0:z
            part_yz = p * (y^2 + z^2)
            if part_yz > needed
                break
            end
            rem = needed - part_yz
            # Loop over x from 0 to B.
            for x in 0:B
                if x^2 > rem
                    break
                end
                candidate = rem - x^2
                # Directly solve for w: there is a valid w when candidate is a perfect square.
                w = isqrt(candidate)
                if w^2 == candidate && w <= x
                    push!(sols, [m, m, w, x, y, z])
                end
            end
        end
    end
    return sols
end

"""
    RHI(p, ell, param)

Given Bp = (-ell, -p | Q) and polarization param = [u0, v0, w0, x0, y0, z0] computes the coefficient matrix of the 5-ary refined Humbert invariant which represent 1. 

Note that polarization is a divisor with zero index (u0 > 0 and v0 > 0) and two self-intersection numbers (u0*v0 - w0^2 - ell*x0^2 - p*y0^2 - ell*p*z0^2 == 1).  

# Examples
```jldoctest
julia> @time RHI(19, 1, [14,31,7,2,4,2])
  0.000723 seconds (2.87 k allocations: 203.000 KiB)
[2    0    0     0     0]
[0   24    8     8     8]
[0    8   32    -8    16]
[0    8   -8    48   -16]
[0    8   16   -16    64]

julia> @time RHI(19, 1, [2,217,7,2,4,2])
  0.000512 seconds (2.54 k allocations: 178.336 KiB)
[2   0   0     0     0]
[0   8   0     0     0]
[0   0   8     0     0]
[0   0   0   152     0]
[0   0   0     0   152]

julia> @time RHI(19, 1, [7,62,7,2,4,2])
  0.000537 seconds (2.82 k allocations: 201.219 KiB)
0
```

# References
[1] U. Fincke and M. E. Pohst, Improved methods for calculating vectors of short length in a lattice, including a complexity analysis, Math. Comp. **44** (1985), no. 170, 463--471; MR0777278
[2] D. Simon, Solving quadratic equations using reduced unimodular quadratic forms, Math. Comp. **74** (2005), no. 251, 1531--1543; MR2137016
[3] E. J. Kani, Elliptic subcovers of a curve of genus 2. I. The isogeny defect, Ann. Math. Qu\'e. **43** (2019), no.~2, 281--303; MR3996071
[4] E. J. Kani, Elliptic subcovers of a curve of genus 2 II. The refined Humbert invariant, J. Number Theory **193** (2018), 302--335; MR3846811
"""
function RHI(p::Integer, ell::Integer, param::Vector{<:Integer})
    u0, v0, w0, x0, y0, z0 = param

    A = zero_matrix(ZZ, 6, 6)

    # Fill diagonal
    A[1, 1] = 2 * v0^2
    A[2, 2] = 2 * u0^2
    A[3, 3] = 8 * (1 + w0^2)
    A[4, 4] = 8 * ell * (1 + ell * x0^2)
    A[5, 5] = 8 * p * (1 + p * y0^2)
    A[6, 6] = 8 * ell * p * (1 + ell * p * z0^2)

    # Fill upper triangular off-diagonals
    A[1, 2] = 2 * (u0 * v0 - 2)
    A[1, 3] = -4 * v0 * w0
    A[1, 4] = -4 * ell * v0 * x0
    A[1, 5] = -4 * p * v0 * y0
    A[1, 6] = 4 * ell * p * v0 * z0

    A[2, 3] = -4 * u0 * w0
    A[2, 4] = -4 * ell * u0 * x0
    A[2, 5] = -4 * p * u0 * y0
    A[2, 6] = 4 * ell * p * u0 * z0

    A[3, 4] = 8 * ell * w0 * x0
    A[3, 5] = 8 * p * w0 * y0
    A[3, 6] = -8 * ell * p * w0 * z0

    A[4, 5] = 8 * ell * p * x0 * y0
    A[4, 6] = -8 * ell^2 * p * x0 * z0

    A[5, 6] = -8 * ell * (p^2) * y0 * z0

    # Mirror the upper triangular part to the lower triangular part
    for i in 2:6
        for j in 1:(i - 1)
            A[i, j] = A[j, i]
        end
    end

    # Check A is positive semidefinite
    V = quadratic_space(QQ, A)
    D = diagonal(V)
    if !all(>=(0), D)
        return 0
    end

    # Transform A to reduce radical
    T = Hecke._complete_to_basis(kernel(A; side=:left))
    AA = T * A * transpose(T)

    # Check rank, symmetry, and last row
    if rank(AA) != 5 || !is_symmetric(AA) || any(!iszero, AA[6, :])
        return 0
    end

    # Work with top-left 5×5 submatrix
    B = @view AA[1:5, 1:5]
    if det(B) == (2^13) * ell^2 * p^2
        L = integer_lattice(; gram=B .÷ 2)
        if is_positive_definite(L)
            LL = lll(L)
            C = gram_matrix(LL)
            if minimum(LL) == 1
                # x1^2 + 4*deg
                return change_base_ring(ZZ, 2 * C)
            end
        end
    end
    return 0
end

"""
    degForm(M)

Given the coefficient matrix M of RHI, it computes the coefficient matrix of the corresponding 4-ary quadratic form called degree form. 

Note that, RHI = X^2 + 4*deg.

# Examples
```jldoctest
julia> M = RHI(19, 1, [14,31,7,2,4,2]);

julia> @time degForm(M)
  0.000013 seconds (64 allocations: 1.438 KiB)
[6    2    2    2]
[2    8   -2    4]
[2   -2   12   -4]
[2    4   -4   16]
```
"""
function degForm(M::ZZMatrix)
    # Check directly without allocating an intermediate matrix.
    if M[1, 1] == 2 && M[1, 2] == 0 && M[1, 3] == 0 && M[1, 4] == 0 && M[1, 5] == 0
        N = @view M[2:5, 2:5]  # Use a view to avoid copying the submatrix.
        return N .÷ 4
    end
    return nothing
end

"""
    allRHI(p,ell,a,b)

Given Bp = (-ell, -p | Q) and integers 0 < a <= b, it saves a list of all unique RHIs and their corresponding degree forms corresponding to polarizations with a <= u = v <= b.

# Examples
```jldoctest
julia> @time allRHI(19, 1, 2*19, 3*19)
working with prime 19
saved data for prime 19
  0.668610 seconds (5.93 M allocations: 329.830 MiB, 20.30% gc time)
```

# References
[1] W. Plesken and B. Souvignier, Computing isometries of lattices, J. Symbolic Comput. **24** (1997), no.~3-4, 327--334; MR1484483
"""
function allRHI(p::Integer, ell::Integer, a::Integer, b::Integer)
    start_time = time()
    println("working with prime ", p)

    # Open file once for writing.
    filename = "./RHI_$(p)_$(ell)_$(a)_$(b).txt"

    file = open(filename, "w")
    try
        # Write global header immediately.
        println(file, "p = ", p)
        println(file, "ℓ = ", ell, "\n")

        idx = 0   # Counting unique forms.
        count = 0 # total polarizations cycled through.
        total = 0 # Total RHIs computed.

        unique_forms = Vector{ZZMatrix}()
        pol_count = Dict{Int,Int}()

        # Process each m in a block.
        for m in a:b
            if rem(m^2, p) != 1
                if ell == 1
                    params = polzSQI(p, m)
                else
                    params = polz(p, ell, m)
                end
                count += length(params)

                # Create a temporary IOBuffer for this m
                m_buffer = IOBuffer()
                #println("m = ", m)
                #println(m_buffer, "m = ", m, "\n")

                for param in params
                    cmA = RHI(p, ell, param)  # ZZMatrix coefficient matrix.
                    if cmA != 0
                        total += 1
                        is_unique = true

                        for (k, cmB) in enumerate(unique_forms)
                            if cmA == cmB
                                is_unique = false
                                pol_count[k] = get(pol_count, k, 1) + 1
                                break
                            end
                        end

                        if is_unique
                            LA = integer_lattice(; gram=cmA .÷ 2)
                            for (k, cmB) in enumerate(unique_forms)
                                LB = integer_lattice(; gram=cmB .÷ 2)
                                if is_isometric(LA, LB)
                                    is_unique = false
                                    pol_count[k] = get(pol_count, k, 1) + 1
                                    break
                                end
                            end
                        end

                        if is_unique
                            idx += 1
                            println(m_buffer, "Type ", idx)
                            #println("Found no. ", idx)
                            println(
                                m_buffer,
                                "θ = [",
                                param[1],
                                "  ",
                                param[3],
                                "+",
                                param[4],
                                "i+",
                                param[5],
                                "j+",
                                param[6],
                                "k]",
                            )
                            println(
                                m_buffer,
                                "    [",
                                param[3],
                                "-",
                                param[4],
                                "i-",
                                param[5],
                                "j-",
                                param[6],
                                "k  ",
                                param[2],
                                "]",
                            )

                            push!(unique_forms, cmA)

                            q = polyForm(cmA)
                            println(m_buffer, "q(E1xE2,θ) = ", q)

                            L = integer_lattice(; gram=cmA .÷ 2)
                            ways = 2*length(short_vectors(L, 4, 4))
                            println(m_buffer, "ways to represent 4 = ", ways)

                            s = automorphism_group_order(L)
                            println(m_buffer, "|Aut(q)| = ", s)

                            cmC = degForm(cmA)
                            deg = polyForm(cmC)
                            println(m_buffer, "deg(E1xE2,θ) = ", deg, "\n")
                        end
                    end
                end
                # Write the chunk for m to the file and flush.
                write(file, String(take!(m_buffer)))
                flush(file)
            end
        end

        println(file, "total polarizations checked: ", count)
        println(file, "total RHI's computed: ", total, "\n")
        println(file, "polarization leading to same type: ", pol_count, "\n")

        end_time = time()
        elapsed_time = end_time - start_time
        hours = floor(elapsed_time / 3600)
        minutes = floor((elapsed_time % 3600) / 60)
        seconds = round(elapsed_time % 60)
        println(file, "Total run time: ", hours, " hrs ", minutes, " min ", seconds, " sec")
    finally
        close(file)
    end
    println("saved data for prime ", p)
    return nothing
end

"""
    minDeg(p,ell,a,b)

Retrieve all the degree forms from the file "RHI_p_a_b.txt" and compute minimum vector for the lattice corresponding to the gram matrix. Finally, returns the frequency distribution of all minimum vectors.

# Examples
```jldoctest 
julia> @time minDeg(19, 1, 2*19, 3*19)
  0.006532 seconds (14.66 k allocations: 759.359 KiB)
Dict{Int64, Int64} with 4 entries:
  4 => 1
  2 => 2
  3 => 2
  1 => 3
```
"""
function minDeg(p::Integer, ell::Integer, a::Integer, b::Integer)
    monos = [
        [2, 0, 0, 0, 0],
        [1, 1, 0, 0, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 0, 1, 0],
        [1, 0, 0, 0, 1],
        [0, 2, 0, 0, 0],
        [0, 1, 1, 0, 0],
        [0, 1, 0, 1, 0],
        [0, 1, 0, 0, 1],
        [0, 0, 2, 0, 0],
        [0, 0, 1, 1, 0],
        [0, 0, 1, 0, 1],
        [0, 0, 0, 2, 0],
        [0, 0, 0, 1, 1],
        [0, 0, 0, 0, 2],
    ]

    filename = "RHI_$(p)_$(ell)_$(a)_$(b).txt"
    prime, forms = fileParser(filename)
    #@assert prime == p "File mismatch"

    NDict = Dict{Int,Int}()

    for tp in keys(forms)
        _, q = forms[tp]
        coeffs = [coeff(q, mono) for mono in monos] #15 elements

        A = ZZ[
            2*coeffs[1] coeffs[2] coeffs[3] coeffs[4];
            coeffs[2] 2*coeffs[6] coeffs[7] coeffs[8];
            coeffs[3] coeffs[7] 2*coeffs[10] coeffs[11];
            coeffs[4] coeffs[8] coeffs[11] 2*coeffs[13]
        ]

        L = integer_lattice(; gram=A .÷ 2)
        N = Int(minimum(L))
        NDict[N] = get(NDict, N, 0) + 1
    end
    return NDict
end
