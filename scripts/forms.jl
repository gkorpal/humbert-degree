using Oscar

"""
    read_polarizations(filename::AbstractString) -> (prime, polarizations)

Parse a `polarizations_p<p>.txt`-style file: `prime` is the integer
from a `p = ...`/`# p = ...` header line (or `nothing` if none was
found), `polarizations` is every `u v w x y z` data row as a
`Vector{Int}`. Malformed data rows (wrong column count, non-integer
entries) are skipped rather than raising.
"""
function read_polarizations(filename::AbstractString)
    prime = nothing
    polarizations = Vector{Vector{Int}}()

    prime_re = r"^\s*#?\s*p\s*=\s*(\d+)"
    data_re = r"^\s*\[\d+\]\s+(.+)$"

    for line in eachline(filename)
        raw = strip(line)
        isempty(raw) && continue

        m = match(prime_re, raw)
        if m !== nothing
            prime = parse(Int, m.captures[1])
            continue
        end
        startswith(raw, "#") && continue

        data = match(data_re, raw)
        raw = data === nothing ? raw : data.captures[1]

        tokens = split(raw)
        length(tokens) < 6 && continue
        values = tryparse.(Int, tokens[1:6])
        any(isnothing, values) && continue
        push!(polarizations, values)
    end

    return prime, polarizations
end

"""
    class_numbers(p::Integer) -> (h, h*(h+1)÷2, H)

Class numbers from Ibukiyama-Katsura-Oort, *"Supersingular curves of
genus two and class numbers"*: `h`, the number of reducible
polarizations' generating count, `h*(h+1)/2`, the number of reducible
polarizations, and `H`, the total number of principal polarizations.
Always returns a plain `NTuple{3,Int}`, including the `p in (2,3,5)`
fast paths.
"""
function class_numbers(p::Integer)
    if p == 2 || p == 3
        return (1, 1, 0)
    elseif p == 5
        return (2, 1, 1)
    end

    a = 1 - jacobi_symbol(-1, p)
    b = 1 - jacobi_symbol(-2, p)
    c = 1 - jacobi_symbol(-3, p)
    d = (p % 5 == 4) ? 4 // 5 : 0 // 1

    H = ((p - 1) * (p + 12) * (p + 23)) // 2880 +
        (a * (2p + 13)) // 96 +
        (c * (p + 11)) // 36 +
        b // 8 +
        (a * c) // 12 +
        d

    r12 = p % 12
    offset = r12 == 1 ? 0 : (r12 == 5 || r12 == 7) ? 1 : r12 == 11 ? 2 : 0
    h = (p - r12) // 12 + offset

    return Int(h), Int(h * (h + 1) ÷ 2), Int(H)
end

"""
    rhi1(p::Integer, param::AbstractVector{<:Integer}) -> Union{Nothing,ZZMatrix}

Build the auxiliary Gram form for polarization `param = [u0,v0,w0,x0,y0,z0]`
at prime `p`, and return its doubled top-left 5x5 block if it passes
the positive-semidefinite / rank / determinant / positive-definite-
minimum-1 gates, or `nothing` if any gate fails.
"""
function rhi1(p::Integer, param::AbstractVector{<:Integer})
    u0, v0, w0, x0, y0, z0 = param

    A = zero_matrix(QQ, 6, 6)

    A[1, 1] = 2 * v0^2
    A[2, 2] = 2 * u0^2
    A[3, 3] = 8 * w0^2 + 8 * w0 * z0 + 2 * z0^2 + 8
    A[4, 4] = 8 * x0^2 + 8 * x0 * y0 + 2 * y0^2 + 8
    A[5, 5] = 2 * x0^2 + 2 * x0 * y0 * p + 2 * x0 * y0 + QQ(1, 2) * y0^2 * p^2 + y0^2 * p + QQ(1, 2) * y0^2 + 2 * p + 2
    A[6, 6] = 2 * w0^2 - 2 * w0 * z0 * p + 2 * w0 * z0 + QQ(1, 2) * z0^2 * p^2 - z0^2 * p + QQ(1, 2) * z0^2 + 2 * p + 2

    A[1, 2] = 2 * u0 * v0 - 4
    A[1, 3] = -4 * v0 * w0 - 2 * v0 * z0
    A[1, 4] = -4 * v0 * x0 - 2 * v0 * y0
    A[1, 5] = -2 * v0 * x0 - v0 * y0 * p - v0 * y0
    A[1, 6] = -2 * v0 * w0 + v0 * z0 * p - v0 * z0

    A[2, 3] = -4 * u0 * w0 - 2 * u0 * z0
    A[2, 4] = -4 * u0 * x0 - 2 * u0 * y0
    A[2, 5] = -2 * u0 * x0 - u0 * y0 * p - u0 * y0
    A[2, 6] = -2 * u0 * w0 + u0 * z0 * p - u0 * z0

    A[3, 4] = 8 * w0 * x0 + 4 * w0 * y0 + 4 * x0 * z0 + 2 * y0 * z0
    A[3, 5] = 4 * w0 * x0 + 2 * w0 * y0 * p + 2 * w0 * y0 + 2 * x0 * z0 + y0 * z0 * p + y0 * z0
    A[3, 6] = 4 * w0^2 - 2 * w0 * z0 * p + 4 * w0 * z0 - z0^2 * p + z0^2 + 4

    A[4, 5] = 4 * x0^2 + 2 * x0 * y0 * p + 4 * x0 * y0 + y0^2 * p + y0^2 + 4
    A[4, 6] = 4 * w0 * x0 + 2 * w0 * y0 - 2 * x0 * z0 * p + 2 * x0 * z0 - y0 * z0 * p + y0 * z0

    A[5, 6] = 2 * w0 * x0 + w0 * y0 * p + w0 * y0 - x0 * z0 * p + x0 * z0 - QQ(1, 2) * y0 * z0 * p^2 + QQ(1, 2) * y0 * z0

    for i in 2:6, j in 1:(i-1)
        A[i, j] = A[j, i]
    end

    V = quadratic_space(QQ, A)
    all(>=(0), diagonal(V)) || return nothing

    Azz = map_entries(x -> ZZ(x // 2), A)
    AA = lll_gram(Azz)

    (rank(AA) == 5 && is_symmetric(AA) && all(iszero, AA[6, :])) || return nothing

    B = @view AA[1:5, 1:5]
    det(B) == (2^4) * p^2 || return nothing

    L = integer_lattice(; gram = B)
    is_positive_definite(L) && minimum(L) == 1 || return nothing

    return B .* 2  # det = 2^9 * p^2
end

"""
    poly_form(M::ZZMatrix) -> polynomial

`f = (1/2) * Σᵢ M[i,i]*xᵢ² + Σ_{i<j} M[i,j]*xᵢ*xⱼ`, the quadratic form
with Gram matrix `M` (an `n`-variable integer polynomial ring is
created for this call).
"""
function poly_form(M::ZZMatrix)
    n = number_of_rows(M)
    R, x = polynomial_ring(ZZ, n)
    f = R(0)
    for i in 1:n
        f += (M[i, i] ÷ 2) * x[i]^2
        for j in (i+1):n
            f += M[i, j] * x[i] * x[j]
        end
    end
    return f
end

"""
    deg_form(M::ZZMatrix) -> Union{Nothing,AbstractMatrix}

`M`'s degree-4 sublattice: the bottom-right 4x4 block divided by 4, if
row 1 of `M` is `[2 0 0 0 0]` (view, no copy); `nothing` otherwise.
"""
function deg_form(M::ZZMatrix)
    all(iszero, @view M[1, 2:5]) && M[1, 1] == 2 || return nothing
    return @view(M[2:5, 2:5]) .÷ 4
end

"""
    sign_str(x::Integer) -> String

`"+"` for `x >= 0`, `"-"` otherwise.
"""
sign_str(x::Integer) = x >= 0 ? "+" : "-"

"""
    all_rhi1(p::Integer)

Read `polarizations_p<p>.txt`, compute `rhi1` for every polarization,
group the results into isometry classes (exact `ZZMatrix` equality
first, then `is_isometric` as a fallback), and write a per-class report
plus summary stats to `./RHI1_<p>.txt`.
"""
function all_rhi1(p::Integer)
    start_time = time()
    println("working with prime ", p)

    filename = "./RHI1_$(p).txt"

    open(filename, "w") do file
        println(file, "p = ", p, "\n")

        h, h2, H = class_numbers(p)
        println(file, "h(h+1)/2 = ", h2)
        println(file, "H = ", H, "\n")

        idx = 0
        total = 0

        unique_forms = Vector{ZZMatrix}()
        pol_count = Dict{Int,Int}()

        prime, params = read_polarizations("polarizations_p$(p).txt")
        count = length(params)

        if prime == p
            for param in params
                cmA = rhi1(p, param)
                cmA === nothing && continue
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
                    LA = integer_lattice(gram = cmA .÷ 2)
                    for (k, cmB) in enumerate(unique_forms)
                        LB = integer_lattice(gram = cmB .÷ 2)
                        if is_isometric(LA, LB)
                            is_unique = false
                            pol_count[k] = get(pol_count, k, 1) + 1
                            break
                        end
                    end
                end
                is_unique || continue

                idx += 1
                push!(unique_forms, cmA)

                u, v, w, x, y, z = param
                println(file, "Type ", idx)
                println(
                    file,
                    "θ = [", u, "  ", w, sign_str(x), abs(x), "β₁",
                    sign_str(y), abs(y), "β₂", sign_str(z), abs(z), "β₃]",
                )
                println(
                    file,
                    "    [", w, sign_str(-x), abs(x), "β₁",
                    sign_str(-y), abs(y), "β₂", sign_str(-z), abs(z), "β₃  ", v, "]",
                )

                q = poly_form(cmA)
                println(file, "q(ExE,θ) = ", q)

                cmC = deg_form(cmA)
                deg = poly_form(cmC)
                println(file, "deg(ExE,θ) = ", deg, "\n")
            end
        end

        println(file, "total polarizations checked: ", count)
        println(file, "total RHI's computed: ", total, "\n")
        println(file, "polarization leading to same type: ", pol_count, "\n")

        elapsed = time() - start_time
        hours, rem = divrem(elapsed, 3600)
        minutes, seconds = divrem(rem, 60)
        println(file, "Total run time: ", floor(hours), " hrs ", floor(minutes), " min ", round(seconds), " sec")
    end
    println("saved data for prime ", p)
    return nothing
end

# for p in 11:12:660
#     if is_probable_prime(p)
#         all_rhi1(p)
#     end
# end