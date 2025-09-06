"""
    primeList(N, n, m)

Computes a list of n pseudoprimes bigger than N and congruent to m mod 4.
    
We take N > 11 so the moduli space of the supersingular surface has more than one element.

# Examples
```jldoctest
julia> @time primeList(19,2,3)
  0.000011 seconds (6 allocations: 240 bytes)
2-element Vector{Int64}:
 19
 23
```
"""
function primeList(N::Integer, n::Integer, m::Integer)
    primes = Vector{Int}()
    p = N
    # Adjust p so that it is congruent to m mod 4.
    rem = p % 4
    if rem != m
        p += (rem - m) % 4
    end
    while length(primes) < n
        if is_probable_prime(p)
            push!(primes, p)
        end
        p += 4  # Only check numbers congruent to m mod 4.
    end
    return primes
end

"""
    classNumbers(p)

Returns the class numbers h, (h(h+1)/2) and H. 

# Examples
```jldoctest
julia> @time classNumbers(19)
  0.000006 seconds (2 allocations: 80 bytes)
3-element Vector{Int64}:
  2
  3
 10
```

# References
[1] K. Hashimoto and T. Ibukiyama, On class numbers of positive definite binary quaternion Hermitian forms, J. Fac. Sci. Univ. Tokyo Sect. IA Math. **27** (1980), no. 3, 549--601; MR0603952
[2] T. Ibukiyama, T. Katsura and F. Oort, Supersingular curves of genus two and class numbers, Compositio Math. **57** (1986), no. 2, 127--152; MR0827350
[3] T. Katsura and F. Oort, Supersingular abelian varieties of dimension two or three and class numbers, *in Algebraic geometry, Sendai, 1985*, 253--281, Adv. Stud. Pure Math., 10, North-Holland, Amsterdam, ; MR0946242
[4] T. Ibukiyama and T. Katsura, On the field of definition of superspecial polarized abelian varieties and type numbers, Compositio Math. **91** (1994), no. 1, 37--46; MR1273924
[5] T. Ibukiyama, Principal polarizations of supersingular abelian surfaces, J. Math. Soc. Japan **72** (2020), no 4, 1161--1180; MR4165927
[6] T. Ibukiyama, Supersingular abelian varieties and quaternion hermitian lattices, in *Theory and Applications of Supersingular Curves and Supersingular Abelian Varieties*, 17--37, RIMS Kokyuroku Bessatsu, B90, Res. Inst. Math. Sci. (RIMS), Kyoto, 2022; MR4521511
"""
function classNumbers(p::Integer)
    # Handle small p quickly H - (h(h+1))/2
    if p == 2 || p == 3
        return [1, 1, 0]
    elseif p == 5
        return [2, 1, 1]
    end

    # Compute H
    # Precompute Jacobi values and remainders
    a = (1 - jacobi_symbol(-1, p))
    b = (1 - jacobi_symbol(-2, p))
    c = (1 - jacobi_symbol(-3, p))
    r5 = p % 5
    d = (r5 == 4) ? (4//5) : 0
    H =
        ((p - 1) * (p + 12) * (p + 23))//2880 +
        (a * (2p + 13))//96 +
        (c * (p + 11))//36 +
        (b//8) +
        ((a * c)//12) +
        d

    # Compute h
    # Use a small lookup table for offset in h
    # Only 1, 5, 7, 11 matter; default 0 otherwise
    r12 = p % 12
    offset = if r12 == 1
        0
    elseif r12 == 5 || r12 == 7
        1
    elseif r12 == 11
        2
    else
        0
    end
    h = ((p - r12)//12 + offset)

    # Return result
    return [Int(h), Int((h * (h + 1)) ÷ 2), Int(H)]
end

"""
    getEll(p,m)

Let Bp = (-ell, -p | Q). Then this function gives us a value of ell >= m.

# Examples
```jldoctest
julia> @time getEll(19,2)
  0.000012 seconds (4 allocations: 64 bytes)
5
```

# References
[1] S. Li, Y. Ouyang and Z. Xu, Endomorphism rings of supersingular elliptic curves over Fp, Finite Fields Appl. **62** (2020); MR4038249
"""
function getEll(p::Integer, m::Integer)
    if is_probable_prime(ZZ(p))
        if mod(p, 4) == 3
            if m > 1
                B = ceil(p*(log(p))^2)
                for q in m:B
                    if mod(q, 4) == 1 &&
                        is_probable_prime(ZZ(q)) &&
                        kronecker_symbol(ZZ(q), ZZ(p)) == 1
                        return Integer(q)
                    end
                end
            else
                return 1
            end
        elseif mod(p, 8) == 5
            if m > 2
                B = ceil(p*(log(p))^2)
                for q in m:B
                    if mod(q, 4) == 3 &&
                        is_probable_prime(ZZ(q)) &&
                        kronecker_symbol(ZZ(q), ZZ(p)) == -1
                        return Integer(q)
                    end
                end
            else
                return 2
            end
        elseif mod(p, 8) == 1
            B = ceil(p*(log(p))^2)
            for q in 3:B
                if mod(q, 4) == 3 &&
                    is_probable_prime(ZZ(q)) &&
                    kronecker_symbol(ZZ(q), ZZ(p)) == -1
                    return Integer(q)
                end
            end
        end
    else
        return 0
    end
end

"""
    polyForm(M)

Given the coefficient matrix M of a quadratic form, this function computes the polynomial form using f = 1/2 * (X^t * M * X).

# Examples
```jldoctest
julia> M = RHI(19, 1, [14,31,7,2,4,2]);

julia> @time polyForm(M)
  0.000080 seconds (194 allocations: 7.375 KiB)
x1^2 + 12*x2^2 + 8*x2*x3 + 8*x2*x4 + 8*x2*x5 + 16*x3^2 - 8*x3*x4 + 16*x3*x5 + 24*x4^2 - 16*x4*x5 + 32*x5^2

julia> N = degForm(M)
[6    2    2    2]
[2    8   -2    4]
[2   -2   12   -4]
[2    4   -4   16]

julia> @time polyForm(N)
  0.000066 seconds (148 allocations: 6.164 KiB)
3*x1^2 + 2*x1*x2 + 2*x1*x3 + 2*x1*x4 + 4*x2^2 - 2*x2*x3 + 4*x2*x4 + 6*x3^2 - 4*x3*x4 + 8*x4^2
```
"""
function polyForm(M::ZZMatrix)
    n = number_of_rows(M)
    R, x = polynomial_ring(ZZ, n)
    f = R(0)
    # Use the formula: f = 1/2 * Σᵢ M[i,i]*x[i]^2 + Σ₍ᵢ<ⱼ₎ M[i,j]*x[i]*x[j]
    for i in 1:n
        f += (M[i, i] ÷ 2) * x[i]^2
        for j in (i + 1):n
            f += M[i, j] * x[i] * x[j]
        end
    end
    return f
end

"""
    fileParser(filename)

Reads in txt file, identifies prime p, and identifies all the degree forms and stores their coefficient matrices.

# Examples
```jldoctest 
julia> prime, forms = fileParser("RHI_19_1_38_57.txt");

julia> types = sort(collect(keys(forms)));

julia> length(types)
8

julia> q_form, deg_form = forms[6];

julia> println("q(E1xE2,θ) = ", q_form)
q(E1xE2,θ) = x1^2 + 8*x2^2 + 8*x2*x4 + 12*x3^2 - 8*x3*x4 + 16*x4^2 + 76*x5^2

julia> println("deg(E1xE2,θ) = ", deg_form)
deg(E1xE2,θ) = 2*x1^2 + 2*x1*x3 + 3*x2^2 - 2*x2*x3 + 4*x3^2 + 19*x4^2
```
"""
function fileParser(filename::String)
    # Define the polynomial ring in 5 variables over ZZ.
    global R, (x1, x2, x3, x4, x5) = polynomial_ring(ZZ, 5)
    PolyElem = typeof(x1^2)

    prime = nothing
    types_data = Dict{Int,NTuple{2,PolyElem}}()

    current_type = 0
    current_q = nothing
    current_deg = nothing

    # regex patterns.
    prime_re = r"^\s*p\s*=\s*(\d+)"
    type_re = r"^\s*Type\s+(\d+)"
    q_re = r"^\s*q\(E1xE2,θ\)\s*=\s*(.+)"
    deg_re = r"^\s*deg\(E1xE2,θ\)\s*=\s*(.+)"

    for line in eachline(filename)
        # Check for the prime line.
        m = match(prime_re, line)
        if m !== nothing
            prime = parse(Int, m.captures[1])
            continue
        end

        m = match(type_re, line)
        if m !== nothing
            if current_type != 0
                types_data[current_type] = (current_q, current_deg)
            end
            current_type = parse(Int, m.captures[1])
            current_q = nothing
            current_deg = nothing
            continue
        end

        if current_type != 0
            m = match(q_re, line)
            if m !== nothing
                q_str = strip(m.captures[1])
                current_q = eval(Meta.parse(q_str))
                continue
            end

            m = match(deg_re, line)
            if m !== nothing
                deg_str = strip(m.captures[1])
                current_deg = eval(Meta.parse(deg_str))
                continue
            end
        end
    end

    if current_type != 0
        types_data[current_type] = (current_q, current_deg)
    end

    return prime, types_data
end
