# === Minimal conditional exact-transform layer ===
#
# A conditional transform is still compiled from the ordinary affine IR. The wrapper stores
# only zero-residual canonicality conditions that could not be established from the symbolic
# coefficients themselves. This is intentionally narrower than the diagnostics API tracked
# by #234.

"""
    ConditionalTransform

An exact affine transformation whose canonicality depends on symbolic zero constraints.
Use [`constraints`](@ref) to inspect the required residual equations. Application returns a
[`ConditionalExpression`](@ref), so those assumptions are not silently discarded.
"""
struct ConditionalTransform{U <: UnitaryTransform}
    transform::U
    residuals::Vector{CNum}

    function ConditionalTransform(transform::U, residuals::Vector{CNum}) where {U <: UnitaryTransform}
        isempty(residuals) && unitary_error(
            "an unconditional transform must be represented by `UnitaryTransform`",
        )
        return new{U}(transform, residuals)
    end
end

"""
    ConditionalExpression

A symbolic operator expression produced under unresolved exact constraints. The underlying
expression is returned by [`conditional_value`](@ref); every expression returned by
[`constraints`](@ref) must equal zero.
"""
struct ConditionalExpression
    value::QAdd
    residuals::Vector{CNum}
end

"""
    constraints(x)

Return the exact residual expressions that must equal zero for a conditional transformation
or expression to be valid. Ordinary `UnitaryTransform`s have no unresolved constraints.
"""
constraints(::UnitaryTransform) = Num[]
constraints(x::ConditionalTransform) = Num[to_num(c) for c in x.residuals]
constraints(x::ConditionalExpression) = Num[to_num(c) for c in x.residuals]

"""
    conditional_value(x::ConditionalExpression)

Return the symbolic expression carried by a conditional result. It is valid only on the
parameter stratum defined by `constraints(x) .== 0`.
"""
conditional_value(x::ConditionalExpression) = x.value

function reduce_constraint(c::CNum, relations::Vector{ParamRelation} = ParamRelation[])
    scratch = ParamRelation[]
    return reduce_all(c, relations, false, scratch)
end

function append_constraint!(out::Vector{CNum}, c::CNum)
    for existing in out
        (isequal(existing, c) || isequal(existing, neg_cnum(c))) && return out
    end
    push!(out, c)
    return out
end

function merge_constraints(first::Vector{CNum}, second::Vector{CNum})
    isempty(first) && return copy(second)
    isempty(second) && return copy(first)
    out = copy(first)
    for residual in second
        append_constraint!(out, residual)
    end
    return out
end

function unresolved_constraints(
        residuals::Vector{CNum}, relations::Vector{ParamRelation} = ParamRelation[];
        context::AbstractString = "conditional transform",
    )
    out = CNum[]
    for residual in residuals
        reduced = reduce_constraint(residual, relations)
        iszero_cnum(reduced) && continue
        is_native(reduced) && unitary_error(
            "$context violates the exact canonicality constraint `$(to_num(reduced)) = 0`",
        )
        append_constraint!(out, reduced)
    end
    return out
end

conditional_result(value::QAdd, residuals::Vector{CNum}) =
    isempty(residuals) ? value : ConditionalExpression(value, residuals)

# Applying a conditional map carries the conditions into the result rather than returning a
# bare QAdd whose domain of validity has been forgotten.
function conjugate(q::QAdd, U::ConditionalTransform)
    return ConditionalExpression(conjugate(q, U.transform), copy(U.residuals))
end

conjugate(o::QSym, U::ConditionalTransform) =
    conjugate(single_qadd(CNUM_ONE, Op[o]), U)

function transform(q::QAdd, U::ConditionalTransform)
    return ConditionalExpression(transform(q, U.transform), copy(U.residuals))
end

transform(o::QSym, U::ConditionalTransform) =
    transform(single_qadd(CNUM_ONE, Op[o]), U)

function conjugate(x::ConditionalExpression, U::UnitaryTransform)
    return ConditionalExpression(conjugate(x.value, U), copy(x.residuals))
end

function conjugate(x::ConditionalExpression, U::ConditionalTransform)
    residuals = merge_constraints(x.residuals, U.residuals)
    return ConditionalExpression(conjugate(x.value, U.transform), residuals)
end

function transform(x::ConditionalExpression, U::UnitaryTransform)
    return ConditionalExpression(transform(x.value, U), copy(x.residuals))
end

function transform(x::ConditionalExpression, U::ConditionalTransform)
    residuals = merge_constraints(x.residuals, U.residuals)
    return ConditionalExpression(transform(x.value, U.transform), residuals)
end

gauge_term(U::ConditionalTransform) = gauge_term(U.transform)
generators(U::ConditionalTransform) = generators(U.transform)

function Base.inv(U::ConditionalTransform)
    return ConditionalTransform(inv(U.transform), copy(U.residuals))
end

Base.adjoint(U::ConditionalTransform) = inv(U)

Base.:*(first::ConditionalTransform, second::ConditionalTransform) = ConditionalTransform(
    first.transform * second.transform,
    merge_constraints(first.residuals, second.residuals),
)
Base.:*(first::ConditionalTransform, second::UnitaryTransform) = ConditionalTransform(
    first.transform * second, copy(first.residuals),
)
Base.:*(first::UnitaryTransform, second::ConditionalTransform) = ConditionalTransform(
    first * second.transform, copy(second.residuals),
)

function substitute_relation(r::ParamRelation, rules::AbstractDict)
    hi = Symbolics.substitute(r.hi, rules)
    lo = Symbolics.substitute(r.lo, rules)
    return ParamRelation(hi, lo, r.sign)
end

function substitute_affine_action(action::AffineAction, rules::AbstractDict)
    n = length(action.basis)
    linear = Matrix{CNum}(undef, n, n)
    shift = Vector{CNum}(undef, n)
    for j in 1:n, i in 1:n
        linear[i, j] = substitute_cnum(action.linear[i, j], rules)
    end
    for i in 1:n
        shift[i] = substitute_cnum(action.shift[i], rules)
    end
    relations = ParamRelation[substitute_relation(r, rules) for r in action.relations]
    return AffineAction(
        action.structure, action.basis, linear, shift; relations = relations,
    )
end

substitute_time(time::StaticTime, ::AbstractDict) = time
function substitute_time(time::DynamicTime, rules::AbstractDict)
    substituted = Symbolics.substitute(time.variable, rules)
    isequal(substituted, time.variable) || unitary_error(
        "substitution of the differentiation variable `$(time.variable)` is not supported " *
            "for a conditional timed transform",
    )
    return time
end

function validate_transform_substitution_rules(rules::AbstractDict)
    for key in keys(rules)
        key isa QSym && unitary_error(
            "conditional-transform substitution resolves scalar parameters only; " *
                "operator substitutions belong on the transformed expression",
        )
    end
    return nothing
end

function conjugation_closed_substitutions(rules::AbstractDict)
    isempty(rules) && return rules
    closed = Dict{Any, Any}(rules)
    sizehint!(closed, 2 * length(rules))
    for (from, to) in rules
        raw_from = if from isa Num
            SymbolicUtils.unwrap(from)
        elseif from isa SymbolicUtils.BasicSymbolic
            from
        else
            continue
        end
        conjugated_from = raw_conj(raw_from)
        isequal(conjugated_from, raw_from) && continue
        conjugated_to = if to isa Number || to isa Num || to isa SymbolicUtils.BasicSymbolic
            conj(to)
        else
            continue
        end
        closed[conjugated_from] = conjugated_to
    end
    return closed
end

function substitute_unitary(
        U::UnitaryTransform{T, AffineAction}, rules::AbstractDict,
    ) where {T}
    validate_transform_substitution_rules(rules)
    isempty(rules) && return U
    time = substitute_time(U.time, rules)
    relations = ParamRelation[substitute_relation(r, rules) for r in U.relations]
    action = substitute_affine_action(U.action, rules)
    forward = Dict{Op, QAdd}()
    inverse = Dict{Op, QAdd}()
    sizehint!(forward, length(U.rules))
    sizehint!(inverse, length(U.inverse_rules))
    for (generator, image) in U.rules
        forward[generator] = substitute(image, rules)
    end
    for (generator, image) in U.inverse_rules
        inverse[generator] = substitute(image, rules)
    end
    gauge = substitute(U.gauge, rules)
    usable = all(is_usable_rel, relations) ? relations : filter(is_usable_rel, relations)
    return UnitaryTransform{T, AffineAction}(
        action, forward, inverse, U.generators, U.sites, gauge, time, usable,
        Val(:validated),
    )
end

function substitute_constraints(
        residuals::Vector{CNum}, rules::AbstractDict,
        relations::Vector{ParamRelation} = ParamRelation[];
        context::AbstractString = "conditional transform",
    )
    substituted = CNum[substitute_cnum(c, rules) for c in residuals]
    return unresolved_constraints(substituted, relations; context = context)
end

function SymbolicUtils.substitute(U::ConditionalTransform, rules::AbstractDict)
    isempty(rules) && return U
    closed_rules = conjugation_closed_substitutions(rules)
    transform = substitute_unitary(U.transform, closed_rules)
    residuals = substitute_constraints(
        U.residuals, closed_rules, transform.relations; context = "conditional transform",
    )
    return isempty(residuals) ? transform : ConditionalTransform(transform, residuals)
end

function SymbolicUtils.substitute(x::ConditionalExpression, rules::AbstractDict)
    isempty(rules) && return x
    closed_rules = conjugation_closed_substitutions(rules)
    value = substitute(x.value, closed_rules)
    residuals = substitute_constraints(
        x.residuals, closed_rules; context = "conditional expression",
    )
    return conditional_result(value, residuals)
end
