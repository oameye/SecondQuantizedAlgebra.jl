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
    n = length(modes)
    basis = Vector{Op}(undef, 2 * n)
    for i in 1:n
        basis[i] = modes[i]
        basis[n + i] = adjoint(modes[i])
    end
    return basis
end

function bogoliubov_matrix(S::AbstractMatrix, n::Int)
    dimension = 2 * n
    size(S) == (dimension, dimension) || unitary_error(
        "`Bogoliubov` needs a $dimension×$dimension Nambu matrix for $n modes; got $(size(S))",
    )
    matrix = Matrix{CNum}(undef, dimension, dimension)
    for j in 1:dimension, i in 1:dimension
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
    dimension = 2 * n
    matrix = Matrix{CNum}(undef, dimension, dimension)
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

function push_bogoliubov_residual!(residuals::Vector{CNum}, residual::CNum)
    reduced = reduce_constraint(residual)
    iszero_cnum(reduced) || append_constraint!(residuals, reduced)
    return residuals
end

# Return exact zero residuals required for the Nambu map to preserve adjoints and the bosonic
# commutator form. A symbolic nonzero residual is a condition, not automatically a failure;
# `unresolved_constraints` distinguishes it from a provably nonzero constant residual.
function bogoliubov_residuals(action::AffineAction)
    action.structure === AFFINE_BOSONIC_NAMBU || unitary_error(
        "internal Bogoliubov validation requires a bosonic Nambu affine action",
    )
    n = length(action.basis)
    half = n ÷ 2
    residuals = CNum[]

    for i in 1:half
        for j in 1:half
            residual = add_cnum(
                action.linear[half + i, half + j],
                neg_cnum(conj_cnum(action.linear[i, j])),
            )
            push_bogoliubov_residual!(residuals, residual)
            residual = add_cnum(
                action.linear[half + i, j],
                neg_cnum(conj_cnum(action.linear[i, half + j])),
            )
            push_bogoliubov_residual!(residuals, residual)
        end
    end

    inverse = inverse_linear(action.linear, action.structure, action.relations)
    for (left, right) in ((action.linear, inverse), (inverse, action.linear))
        for j in 1:n, i in 1:n
            residual = i == j ? CNUM_NEG1 : CNUM_ZERO
            for k in 1:n
                residual = add_cnum(residual, mul_cnum(left[i, k], right[k, j]))
            end
            push_bogoliubov_residual!(residuals, residual)
        end
    end
    return residuals
end

function bogoliubov_action(modes::Vector{Op}, matrix::Matrix{CNum})
    basis = bogoliubov_basis(modes)
    return AffineAction(
        BosonicNambu(), basis, matrix, fill(CNUM_ZERO, length(basis)),
    )
end

function exact_bogoliubov(modes::Vector{Op}, matrix::Matrix{CNum})
    action = bogoliubov_action(modes, matrix)
    residuals = unresolved_constraints(
        bogoliubov_residuals(action), action.relations; context = "`Bogoliubov`",
    )
    transform = canonical_transform(action)
    return isempty(residuals) ? transform : ConditionalTransform(transform, residuals)
end

# Native real matrices cannot produce an unresolved symbolic condition. Keeping this path
# separate gives the value-independent exact API a concrete inferred return type while the
# symbolic path is free to return `ConditionalTransform` when necessary.
const NativeBogoliubovReal = Union{Integer, Rational, AbstractFloat}

function strict_bogoliubov(modes::Vector{Op}, matrix::Matrix{CNum})
    action = bogoliubov_action(modes, matrix)
    residuals = unresolved_constraints(
        bogoliubov_residuals(action), action.relations; context = "`Bogoliubov`",
    )
    isempty(residuals) || unitary_error(
        "native Bogoliubov coefficients unexpectedly produced unresolved canonicality",
    )
    return canonical_transform(action)
end

"""
    Bogoliubov(modes, S)

Construct an exact bosonic Bogoliubov transformation in Nambu ordering
`(a₁, …, aₙ, a₁', …, aₙ')`. `S` must preserve both adjoints and the bosonic
commutator form. Proven exact matrices return an ordinary [`UnitaryTransform`](@ref),
provably noncanonical matrices are rejected, and unresolved symbolic canonicality returns a
[`ConditionalTransform`](@ref) whose [`constraints`](@ref) must vanish.
"""
function Bogoliubov(modes::AbstractVector{Op}, S::AbstractMatrix)
    lowering_modes = bogoliubov_modes(modes)
    return exact_bogoliubov(
        lowering_modes, bogoliubov_matrix(S, length(lowering_modes)),
    )
end

function Bogoliubov(
        modes::AbstractVector{Op}, S::AbstractMatrix{T},
    ) where {T <: NativeBogoliubovReal}
    lowering_modes = bogoliubov_modes(modes)
    return strict_bogoliubov(
        lowering_modes, bogoliubov_matrix(S, length(lowering_modes)),
    )
end

Bogoliubov(mode::Op, S::AbstractMatrix) = Bogoliubov(Op[mode], S)
Bogoliubov(modes::Tuple{Vararg{Op}}, S::AbstractMatrix) = Bogoliubov(Op[modes...], S)

"""
    Bogoliubov(modes, U, V)

Construct the exact bosonic map `a ↦ U*a + V*a'`. The implied Nambu matrix is
`[U V; conj(V) conj(U)]`. Proven exact blocks return an ordinary
[`UnitaryTransform`](@ref), while unresolved symbolic canonicality is returned explicitly as
a [`ConditionalTransform`](@ref).
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

function Bogoliubov(
        modes::AbstractVector{Op}, U::AbstractMatrix{TU}, V::AbstractMatrix{TV},
    ) where {TU <: NativeBogoliubovReal, TV <: NativeBogoliubovReal}
    lowering_modes = bogoliubov_modes(modes)
    return strict_bogoliubov(
        lowering_modes,
        bogoliubov_matrix(U, V, length(lowering_modes)),
    )
end

Bogoliubov(mode::Op, U::AbstractMatrix, V::AbstractMatrix) = Bogoliubov(Op[mode], U, V)
Bogoliubov(modes::Tuple{Vararg{Op}}, U::AbstractMatrix, V::AbstractMatrix) =
    Bogoliubov(Op[modes...], U, V)
