# Internal affine representation for exact total canonical transformations.
#
# `AffineAction` is the construction IR. `UnitaryTransform` remains the compiled execution
# representation used by `conjugate` and `transform`.
abstract type AffineStructure end

"""Fallback structure for an exact affine map with no stronger algebraic inverse formula."""
struct GenericAffine <: AffineStructure end

"""Bosonic Nambu ordering `(a₁,…,aₙ,a₁†,…,aₙ†)`."""
struct BosonicNambu <: AffineStructure end

"""Canonical phase-space ordering `(x₁,…,xₙ,p₁,…,pₙ)`."""
struct SymplecticPhaseSpace <: AffineStructure end

"""Real orthogonal action, used for spin/Pauli rotations."""
struct OrthogonalAction <: AffineStructure end

"""Unitary linear action in the chosen generator basis, e.g. matrix-unit conjugation."""
struct UnitaryLinearAction <: AffineStructure end

struct AffineAction{S <: AffineStructure}
    structure::S
    basis::Vector{Op}
    linear::Matrix{CNum}
    shift::Vector{CNum}
    relations::Vector{ParamRelation}
end

function AffineAction(
        structure::S, basis::Vector{Op}, linear::AbstractMatrix, shift::AbstractVector;
        relations::Vector{ParamRelation} = ParamRelation[],
    ) where {S <: AffineStructure}
    n = length(basis)
    size(linear) == (n, n) || unitary_error(
        "an affine action on $n generators needs a $n×$n linear map; got $(size(linear))",
    )
    length(shift) == n || unitary_error(
        "an affine action on $n generators needs $n shifts; got $(length(shift))",
    )
    length(Set(basis)) == n || unitary_error("an affine action basis cannot contain duplicates")

    coefficients = Matrix{CNum}(undef, n, n)
    offsets = Vector{CNum}(undef, n)
    for j in 1:n, i in 1:n
        coefficients[i, j] = to_cnum(linear[i, j])
    end
    for i in 1:n
        offsets[i] = to_cnum(shift[i])
    end
    return AffineAction{S}(
        structure, copy(basis), coefficients, offsets, copy(relations),
    )
end

AffineAction(basis::Vector{Op}, linear::AbstractMatrix, shift::AbstractVector; kwargs...) =
    AffineAction(GenericAffine(), basis, linear, shift; kwargs...)

function reduce_affine(c::CNum, relations::Vector{ParamRelation}, scratch::Vector{ParamRelation})
    isempty(relations) && return c
    return reduce_all(c, relations, true, scratch)
end

function dagger_linear(linear::Matrix{CNum})
    n, m = size(linear)
    out = Matrix{CNum}(undef, m, n)
    for j in 1:m, i in 1:n
        out[j, i] = conj_cnum(linear[i, j])
    end
    return out
end

function transpose_linear(linear::Matrix{CNum})
    n, m = size(linear)
    out = Matrix{CNum}(undef, m, n)
    for j in 1:m, i in 1:n
        out[j, i] = linear[i, j]
    end
    return out
end

function inverse_linear(linear::Matrix{CNum}, ::BosonicNambu, ::Vector{ParamRelation})
    n = size(linear, 1)
    iseven(n) || unitary_error("a bosonic Nambu action needs an even-dimensional basis")
    half = n ÷ 2
    out = dagger_linear(linear)
    for j in 1:n, i in 1:n
        left_sign = i <= half ? 1 : -1
        right_sign = j <= half ? 1 : -1
        left_sign == right_sign || (out[i, j] = neg_cnum(out[i, j]))
    end
    return out
end

function inverse_linear(
        linear::Matrix{CNum}, ::SymplecticPhaseSpace, ::Vector{ParamRelation},
    )
    n = size(linear, 1)
    iseven(n) || unitary_error("a phase-space action needs an even-dimensional basis")
    half = n ÷ 2
    transposed = transpose_linear(linear)
    out = Matrix{CNum}(undef, n, n)
    # Ω⁻¹ Aᵀ Ω for Ω = [0 I; -I 0].
    for j in 1:n, i in 1:n
        source_i = i <= half ? i + half : i - half
        source_j = j <= half ? j + half : j - half
        value = transposed[source_i, source_j]
        # Left multiplication by Ω⁻¹ and right multiplication by Ω contribute opposite
        # signs when exactly one index belongs to the momentum block.
        (i <= half) == (j <= half) || (value = neg_cnum(value))
        out[i, j] = value
    end
    return out
end

inverse_linear(
    linear::Matrix{CNum}, ::OrthogonalAction, ::Vector{ParamRelation},
) = transpose_linear(linear)

inverse_linear(
    linear::Matrix{CNum}, ::UnitaryLinearAction, ::Vector{ParamRelation},
) = dagger_linear(linear)

function inverse_linear(
        linear::Matrix{CNum}, ::GenericAffine, relations::Vector{ParamRelation},
    )
    n = size(linear, 1)
    augmented = Matrix{CNum}(undef, n, 2n)
    for j in 1:n, i in 1:n
        augmented[i, j] = linear[i, j]
        augmented[i, n + j] = i == j ? CNUM_ONE : CNUM_ZERO
    end

    scratch = ParamRelation[]
    for column in 1:n
        pivot_offset = findfirst(
            row -> !iszero_cnum(reduce_affine(augmented[row, column], relations, scratch)),
            column:n,
        )
        pivot_offset === nothing && unitary_error("affine linear map is singular")
        pivot = column - 1 + pivot_offset
        if pivot != column
            for j in 1:(2n)
                augmented[column, j], augmented[pivot, j] =
                    augmented[pivot, j], augmented[column, j]
            end
        end

        pivot_coefficient = reduce_affine(augmented[column, column], relations, scratch)
        pivot_inverse = inv(pivot_coefficient)
        for j in 1:(2n)
            augmented[column, j] = reduce_affine(
                mul_cnum(augmented[column, j], pivot_inverse), relations, scratch,
            )
        end

        for row in 1:n
            row == column && continue
            factor = reduce_affine(augmented[row, column], relations, scratch)
            iszero_cnum(factor) && continue
            for j in 1:(2n)
                augmented[row, j] = reduce_affine(
                    add_cnum(
                        augmented[row, j],
                        neg_cnum(mul_cnum(factor, augmented[column, j])),
                    ),
                    relations,
                    scratch,
                )
            end
        end
    end

    inverse = Matrix{CNum}(undef, n, n)
    for j in 1:n, i in 1:n
        inverse[i, j] = reduce_affine(augmented[i, n + j], relations, scratch)
    end
    return inverse
end

function Base.inv(action::AffineAction)
    linear = inverse_linear(action.linear, action.structure, action.relations)
    n = length(action.basis)
    shift = Vector{CNum}(undef, n)
    scratch = ParamRelation[]
    for i in 1:n
        value = CNUM_ZERO
        for j in 1:n
            value = add_cnum(value, mul_cnum(linear[i, j], action.shift[j]))
        end
        shift[i] = reduce_affine(neg_cnum(value), action.relations, scratch)
    end
    return AffineAction(
        action.structure, action.basis, linear, shift; relations = action.relations,
    )
end

function affine_rules(action::AffineAction)
    n = length(action.basis)
    rules = Dict{Op, QAdd}()
    sizehint!(rules, n)
    for i in 1:n
        pairs = Tuple{CNum, Vector{Op}}[]
        sizehint!(pairs, n + 1)
        for j in 1:n
            coefficient = action.linear[i, j]
            iszero_cnum(coefficient) || push!(pairs, (coefficient, Op[action.basis[j]]))
        end
        offset = action.shift[i]
        iszero_cnum(offset) || push!(pairs, (offset, Op[]))
        rules[action.basis[i]] = rule_qadd(pairs)
    end
    return rules
end

function static_transform(action::AffineAction)
    inverse_action = inv(action)
    return static_transform(
        affine_rules(action), affine_rules(inverse_action), action.relations,
    )
end
