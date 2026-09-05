# Internal affine representation for exact total canonical transformations.
#
# `AffineAction` is the construction IR. `UnitaryTransform` remains the compiled execution
# representation used by `conjugate` and `transform`.
struct AffineAction
    basis::Vector{Op}
    linear::Matrix{CNum}
    shift::Vector{CNum}
    relations::Vector{ParamRelation}
end

function AffineAction(
        basis::Vector{Op}, linear::AbstractMatrix, shift::AbstractVector;
        relations::Vector{ParamRelation} = ParamRelation[],
    )
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
    return AffineAction(copy(basis), coefficients, offsets, copy(relations))
end

function reduce_affine(c::CNum, relations::Vector{ParamRelation}, scratch::Vector{ParamRelation})
    isempty(relations) && return c
    return reduce_all(c, relations, true, scratch)
end

function inverse_linear(
        linear::Matrix{CNum}, relations::Vector{ParamRelation},
    )
    n = size(linear, 1)
    augmented = Matrix{CNum}(undef, n, 2n)
    for j in 1:n, i in 1:n
        augmented[i, j] = linear[i, j]
        augmented[i, n + j] = i == j ? CNUM_ONE : CNUM_ZERO
    end

    scratch = ParamRelation[]
    for column in 1:n
        pivot_row = findfirst(
            row -> !iszero_cnum(reduce_affine(augmented[row, column], relations, scratch)),
            column:n,
        )
        pivot_row === nothing && unitary_error("affine linear map is singular")
        pivot = column - 1 + pivot_row
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
    linear = inverse_linear(action.linear, action.relations)
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
    return AffineAction(action.basis, linear, shift; relations = action.relations)
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
