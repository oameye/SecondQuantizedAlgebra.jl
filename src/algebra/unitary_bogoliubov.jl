# === General exact bosonic Bogoliubov transformations ===

function bogoliubov_modes(modes::AbstractVector{Op})
    isempty(modes) && unitary_error("`Bogoliubov` needs at least one Fock mode")
    lowering_modes = Op[]
    sizehint!(lowering_modes, length(modes))
    seen = Set{SiteKey}()
    for mode in modes
        lowering = fock_or_throw(mode, "`Bogoliubov`")
        key = site_key(lowering)
        key in seen && unitary_error(
            "`Bogoliubov` received the same Fock mode more than once: `$lowering`",
        )
        push!(seen, key)
        push!(lowering_modes, lowering)
    end
    return lowering_modes
end

bogoliubov_modes(mode::Op) = bogoliubov_modes(Op[mode])
bogoliubov_modes(modes::Tuple{Vararg{Op}}) = bogoliubov_modes(Op[modes...])

function bogoliubov_basis(modes::Vector{Op})
    basis = Vector{Op}(undef, 2length(modes))
    n = length(modes)
    for i in 1:n
        basis[i] = modes[i]
        basis[n + i] = adjoint(modes[i])
    end
    return basis
end

function bogoliubov_matrix(S::AbstractMatrix, n::Int)
    size(S) == (2n, 2n) || unitary_error(
        "`Bogoliubov` needs a $(2n)×$(2n) Nambu matrix for $n modes; got $(size(S))",
    )
    matrix = Matrix{CNum}(undef, 2n, 2n)
    for j in 1:(2n), i in 1:(2n)
        matrix[i, j] = to_cnum(S[i, j])
    end
    return matrix
end

function bogoliubov_matrix(U::AbstractMatrix, V::AbstractMatrix, n::Int)
    size(U) == (n, n) || unitary_error(
        "`Bogoliubov` needs an $n×$n `U` block; got $(size(U))",
    )
    size(V) == (n, n) || unitary_error(
        "`Bogoliubov` needs an $n×$n `V` block; got $(size(V))",
    )
    matrix = Matrix{CNum}(undef, 2n, 2n)
    for j in 1:n, i in 1:n
        u = to_cnum(U[i, j])
        v = to_cnum(V[i, j])
        matrix[i, j] = u
        matrix[i, n + j] = v
        matrix[n + i, j] = conj_cnum(v)
        matrix[n + i, n + j] = conj_cnum(u)
    end
    return matrix
end

function exact_bogoliubov(modes::Vector{Op}, matrix::Matrix{CNum})
    basis = bogoliubov_basis(modes)
    action = AffineAction(
        BosonicNambu(), basis, matrix, fill(CNUM_ZERO, length(basis)),
    )
    return static_transform(action)
end

"""
    Bogoliubov(modes, S)

Construct an exact bosonic Bogoliubov transformation in Nambu ordering
`(a₁, …, aₙ, a₁', …, aₙ')`. `S` must preserve both adjoints and the bosonic
commutator form exactly. Symbolic matrices whose canonicality cannot be proven are rejected.
"""
function Bogoliubov(modes::AbstractVector{Op}, S::AbstractMatrix)
    lowering_modes = bogoliubov_modes(modes)
    return exact_bogoliubov(
        lowering_modes, bogoliubov_matrix(S, length(lowering_modes)),
    )
end

Bogoliubov(mode::Op, S::AbstractMatrix) = Bogoliubov(Op[mode], S)
Bogoliubov(modes::Tuple{Vararg{Op}}, S::AbstractMatrix) = Bogoliubov(Op[modes...], S)

"""
    Bogoliubov(modes, U, V)

Construct the exact bosonic map `a ↦ U*a + V*a'`. The implied Nambu matrix is
`[U V; conj(V) conj(U)]` and must satisfy the bosonic canonical relations exactly.
"""
function Bogoliubov(
        modes::AbstractVector{Op}, U::AbstractMatrix, V::AbstractMatrix,
    )
    lowering_modes = bogoliubov_modes(modes)
    return exact_bogoliubov(
        lowering_modes,
        bogoliubov_matrix(U, V, length(lowering_modes)),
    )
end

Bogoliubov(mode::Op, U::AbstractMatrix, V::AbstractMatrix) = Bogoliubov(Op[mode], U, V)
Bogoliubov(modes::Tuple{Vararg{Op}}, U::AbstractMatrix, V::AbstractMatrix) =
    Bogoliubov(Op[modes...], U, V)
