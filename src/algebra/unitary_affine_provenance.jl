# Proof-preserving helpers for affine actions that are canonical by construction.
#
# `static_transform(action)` remains the validation boundary for caller-supplied/raw affine
# data. Named constructors, exact generated flows, inversion, and composition arrive here
# only after their canonicality has already been established structurally.

function canonical_affine_inverse(action::AffineAction)
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

function canonical_transform(action::AffineAction)
    inverse_action = canonical_affine_inverse(action)
    return validated_transform(
        affine_rules(action), affine_rules(inverse_action), zero_qadd(), StaticTime(),
        action.relations, action,
    )
end

compiled_inverse_action_metadata(action::AffineAction) = canonical_affine_inverse(action)

function compile_composed_action_metadata(action::AffineAction)
    inverse_action = canonical_affine_inverse(action)
    return (affine_rules(action), affine_rules(inverse_action))
end
