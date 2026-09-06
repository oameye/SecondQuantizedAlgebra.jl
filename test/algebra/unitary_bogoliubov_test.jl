using SecondQuantizedAlgebra
using LinearAlgebra: exp
using QuantumOpticsBase: FockBasis
using Test
using Symbolics: @variables
import SecondQuantizedAlgebra: expim

@testset "Exact bosonic Bogoliubov transformations" begin
    fock = FockSpace(:fock)
    a = Destroy(fock, :a)

    two_modes = FockSpace(:left) ⊗ FockSpace(:right)
    left = Destroy(two_modes, :left, 1)
    right = Destroy(two_modes, :right, 2)

    @variables r ϕ θ

    same_constraints(x, y) = begin
        left_constraints = constraints(x)
        right_constraints = constraints(y)
        length(left_constraints) == length(right_constraints) &&
            all(isequal(l, r) for (l, r) in zip(left_constraints, right_constraints))
    end

    @testset "single-mode squeezing is a Bogoliubov map" begin
        S = [
            cosh(r) expim(ϕ) * sinh(r)
            expim(-ϕ) * sinh(r) cosh(r)
        ]
        raw = Bogoliubov(a, S)
        named = Squeeze(a, r, ϕ)

        @test raw isa UnitaryTransform
        for op in (a, a')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
            @test iszero(
                simplify(conjugate(op, inv(raw)) - conjugate(op, inv(named))),
            )
        end
        @test iszero(simplify(conjugate(a, raw) - cosh(r) * a - expim(ϕ) * sinh(r) * a'))
    end

    @testset "passive two-mode mixing uses the same Nambu core" begin
        U = [cos(θ) sin(θ); -sin(θ) cos(θ)]
        V = [0 0; 0 0]
        raw = Bogoliubov((left, right), U, V)
        named = BeamSplitter(left, right, θ)

        @test raw isa UnitaryTransform
        for op in (left, right, left', right')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
        end
    end

    @testset "two-mode squeezing uses the same Nambu core" begin
        U = [cosh(r) 0; 0 cosh(r)]
        V = [0 sinh(r); sinh(r) 0]
        raw = Bogoliubov(Op[left, right], U, V)
        named = TwoModeSqueeze(left, right, r)

        @test raw isa UnitaryTransform
        for op in (left, right, left', right')
            @test iszero(simplify(conjugate(op, raw) - conjugate(op, named)))
        end
        @test iszero(
            simplify(
                conjugate(left, raw) - cosh(r) * left - sinh(r) * right',
            ),
        )
    end

    @testset "unresolved symbolic canonicality is explicit" begin
        @variables u::Number v::Number
        conditional = Bogoliubov(a, [u v; conj(v) conj(u)])
        @test conditional isa ConditionalTransform
        @test !isempty(constraints(conditional))
        required = simplify(u * conj(u) - v * conj(v) - 1)
        @test any(c -> isequal(simplify(c), required), constraints(conditional))

        applied = conjugate(a, conditional)
        @test applied isa ConditionalExpression
        @test same_constraints(applied, conditional)
        @test iszero(
            simplify(conditional_value(applied) - u * a - v * a'),
        )

        inverse = inv(conditional)
        @test inverse isa ConditionalTransform
        @test same_constraints(inverse, conditional)

        displacement = Displace(a, 1 // 3)
        mixed = conditional * displacement
        @test mixed isa ConditionalTransform
        mixed_applied = conjugate(a, mixed)
        sequential = conjugate(applied, displacement)
        @test mixed_applied isa ConditionalExpression
        @test sequential isa ConditionalExpression
        @test same_constraints(mixed_applied, conditional)
        @test same_constraints(sequential, conditional)
        @test iszero(
            simplify(
                conditional_value(mixed_applied) - conditional_value(sequential),
            ),
        )

        doubled = conditional * conditional
        @test doubled isa ConditionalTransform
        @test same_constraints(doubled, conditional)

        # Exact rational values discharge the canonicality condition without relying on
        # heuristic simplification of identities such as conj(cosh(r)) == cosh(r).
        resolved = substitute(conditional, Dict(u => 5 // 3, v => 4 // 3))
        @test resolved isa UnitaryTransform
        @test isempty(constraints(resolved))
        @test iszero(
            simplify(conjugate(a, resolved) - (5 // 3) * a - (4 // 3) * a'),
        )

        resolved_applied = substitute(applied, Dict(u => 5 // 3, v => 4 // 3))
        @test !(resolved_applied isa ConditionalExpression)
        @test iszero(
            simplify(resolved_applied - (5 // 3) * a - (4 // 3) * a'),
        )

        partial = substitute(conditional, Dict(u => 5 // 3))
        @test partial isa ConditionalTransform
        @test !isempty(constraints(partial))
        @test_throws ArgumentError substitute(conditional, Dict(u => 1, v => 1))
    end

    @testset "exact canonicality is enforced" begin
        @test_throws ArgumentError Bogoliubov(Op[], Matrix{Int}(undef, 0, 0))
        @test_throws ArgumentError Bogoliubov(a, [1 1; 1 1])
        @test_throws ArgumentError Bogoliubov(a, [1 0; 0 2])
        @test_throws ArgumentError Bogoliubov(a, [1 1; 0 1])
        @test_throws ArgumentError Bogoliubov(a, [1 0 0; 0 1 0])
        @test_throws ArgumentError Bogoliubov((left, left), [1 0; 0 1], [0 0; 0 0])
        @test_throws ArgumentError Bogoliubov((left, right), [1 0; 0 1], [1 0; 0 0])
        @test_throws ArgumentError Bogoliubov((left, right), [1 0; 0 1; 0 0], [0 0; 0 0])
        @test_throws ArgumentError Bogoliubov((left, right), [1 0; 0 1], [0 0 0; 0 0 0])

        identity_map = @inferred Bogoliubov(a, [1 0; 0 1])
        @test identity_map isa UnitaryTransform
        @test iszero(simplify(conjugate(a, identity_map) - a))
        @test iszero(simplify(conjugate(a, inv(identity_map)) - a))

        uv_identity = @inferred Bogoliubov(
            (left, right), [1 0; 0 1], [0 0; 0 0],
        )
        @test uv_identity isa UnitaryTransform
    end

    @testset "raw Bogoliubov map has an independent numeric oracle" begin
        dimension = 24
        lowering = zeros(ComplexF64, dimension, dimension)
        for n in 2:dimension
            lowering[n - 1, n] = sqrt(n - 1)
        end

        # Rational hyperbolic data are exactly canonical: c²-s² = 1. The matrix oracle
        # below is numerical and independent of the symbolic canonicality proof.
        c = 401 // 399
        s = 40 // 399
        S = [c s; s c]
        raw = Bogoliubov(a, S)
        actual = Matrix(to_numeric(conjugate(a, raw), FockBasis(dimension)).data)
        r0 = asinh(Float64(s))
        unitary = exp((r0 / 2) * ((lowering')^2 - lowering^2))
        expected = unitary' * lowering * unitary
        @test actual[1:6, 1:6] ≈ expected[1:6, 1:6] atol = 1.0e-10
    end
end
