using Test
using SecondQuantizedAlgebra
using Symbolics: @variables, Num
import SecondQuantizedAlgebra: Op, expim, to_num

@testset "Hamiltonian-derived displacement" begin
    fock = FockSpace(:automatic_fock)
    a = Destroy(fock, :a)

    phase = PhaseSpace(:automatic_phase)
    x = Position(phase, :x)
    p = Momentum(phase, :p)

    @variables ω Ω η g t dx dp ωd K
    @variables α::Number
    @variables envelope(t)

    @testset "static Fock equilibrium" begin
        reference = ω * a' * a + η * (a + a') + g
        U = Displace(a, reference)
        @test U isa UnitaryTransform
        @test isequal(conjugate(a, U), a - η / ω)
        @test isequal(transform(reference, U), ω * a' * a - η^2 / ω + g)

        complex_reference = ω * a' * a + α * a' + conj(α) * a
        complex_U = Displace(a, complex_reference)
        @test isequal(conjugate(a, complex_U), a - α / ω)
        @test isequal(conjugate(a', complex_U), conjugate(a, complex_U)')
    end

    @testset "bounded harmonic Fock response" begin
        reference = ω * a' * a - im * Ω * cos(ωd * t) * (a - a')
        U = Displace(a, reference, t)
        @test U isa UnitaryTransform
        transformed = transform(reference, U)
        operator_terms = [
            (term.ops, coefficient) for
                (term, coefficient) in transformed if !isempty(term.ops)
        ]
        @test length(operator_terms) == 1
        @test operator_terms[1][1] == Op[a', a]
        @test isequal(to_num(operator_terms[1][2]), Complex(Num(ω), Num(0)))

        constant_reference = ω * a' * a + η * (a + a') + g
        constant = Displace(a, constant_reference, t)
        @test isequal(conjugate(a, constant), a - η / ω)

        multitone_drive = η * cos(ωd * t) + g * sin(2ωd * t)
        multitone_reference = ω * a' * a + multitone_drive * (a + a')
        multitone = Displace(a, multitone_reference, t)
        multitone_terms = [
            term.ops for (term, _) in transform(multitone_reference, multitone) if
                !isempty(term.ops)
        ]
        @test multitone_terms == [Op[a', a]]
        @test isequal(conjugate(a', multitone), conjugate(a, multitone)')

        from_creation = Displace(a', reference, t)
        @test isequal(conjugate(a, from_creation), conjugate(a, U))

        Kerr = (K / 2) * a'^2 * a^2
        @test isequal(
            transform(reference + Kerr, U),
            transform(reference, U) + conjugate(Kerr, U),
        )
    end

    @testset "Fock validation" begin
        reference = ω * a' * a - im * Ω * cos(ωd * t) * (a - a')
        other = Destroy(FockSpace(:automatic_other), :b)
        @test_throws ArgumentError Displace(a, reference + K * a'^2 * a^2, t)
        @test_throws ArgumentError Displace(a, reference + g * other, t)
        @test_throws ArgumentError Displace(a, ω * a' * a + η * a', t)
        @test_throws ArgumentError Displace(
            a, (ω + g * cos(ωd * t)) * a' * a + η * (a + a'), t,
        )
        @test_throws ArgumentError Displace(a, ω * a' * a + envelope * (a + a'), t)

        nonlinear_phase = expim(t^2)
        @test_throws ArgumentError Displace(
            a, ω * a' * a + nonlinear_phase * a' + conj(nonlinear_phase) * a, t,
        )
        resonant = expim(-ω * t)
        @test_throws ArgumentError Displace(
            a, ω * a' * a + resonant * a' + conj(resonant) * a, t,
        )
        @test_throws ArgumentError Displace(a, η * (a + a'), t)
        @test_throws ArgumentError Displace(a, η * (a + a'))

        # A possible symbolic resonance remains an exact quotient by design.
        @test Displace(a, reference, t) isa UnitaryTransform
    end

    @testset "static quadrature equilibrium" begin
        reference =
            (ω / 2) * x^2 + (g / 2) * (x * p + p * x) +
            (Ω / 2) * p^2 + η * x + dx * p
        U = Displace(x, p, reference)
        @test U isa UnitaryTransform
        determinant = ω * Ω - g^2
        @test isequal(conjugate(x, U), x + (g * dx - Ω * η) / determinant)
        @test isequal(conjugate(p, U), p + (g * η - ω * dx) / determinant)
        operator_terms = [
            term.ops for (term, _) in simplify(transform(reference, U)) if
                !isempty(term.ops)
        ]
        @test Set(operator_terms) == Set((Op[x, x], Op[x, p], Op[p, p]))
    end

    @testset "bounded harmonic quadrature response" begin
        isotropic_reference = (ω / 2) * (x^2 + p^2) + η * cos(ωd * t) * x
        U = Displace(x, p, isotropic_reference, t)
        @test U isa UnitaryTransform
        operator_terms = [
            term.ops for (term, _) in transform(isotropic_reference, U) if
                !isempty(term.ops)
        ]
        @test Set(operator_terms) == Set((Op[x, x], Op[p, p]))

        multitone_reference =
            (ω / 2) * x^2 + (g / 2) * (x * p + p * x) +
            (Ω / 2) * p^2 +
            (η * cos(ωd * t) + dx * sin(2ωd * t)) * x +
            dp * cos(3ωd * t) * p
        multitone = Displace(x, p, multitone_reference, t)
        multitone_terms = [
            term.ops for (term, _) in transform(multitone_reference, multitone) if
                !isempty(term.ops)
        ]
        @test Set(multitone_terms) == Set((Op[x, x], Op[x, p], Op[p, p]))
    end

    @testset "quadrature validation" begin
        reference = (ω / 2) * (x^2 + p^2) + η * cos(ωd * t) * x
        other_phase = PhaseSpace(:automatic_selected) ⊗ PhaseSpace(:automatic_other_phase)
        other_x = Position(other_phase, :other_x, 2)
        nonlinear_phase = expim(t^2)

        @test_throws ArgumentError Displace(p, x, reference, t)
        @test_throws ArgumentError Displace(x, p, reference + other_x, t)
        @test_throws ArgumentError Displace(x, p, reference + K * x^3, t)
        @test_throws ArgumentError Displace(x, p, reference + im * x, t)
        @test_throws ArgumentError Displace(
            x, p,
            ((ω + g * cos(ωd * t)) / 2) * x^2 + (Ω / 2) * p^2 + η * x,
            t,
        )
        @test_throws ArgumentError Displace(
            x, p, (ω / 2) * (x^2 + p^2) + envelope * x, t,
        )
        @test_throws ArgumentError Displace(
            x, p,
            (ω / 2) * (x^2 + p^2) + (nonlinear_phase + conj(nonlinear_phase)) * x,
            t,
        )
        @test_throws ArgumentError Displace(
            x, p, (ω / 2) * (x^2 + p^2) + η * cos(ω * t) * x, t,
        )
        @test_throws ArgumentError Displace(x, p, (ω / 2) * x^2 + η * p)

        nonunique = Displace(x, p, (ω / 2) * x^2)
        @test isequal(conjugate(x + p, nonunique), x + p)

        # A possible symbolic determinant resonance remains an exact quotient.
        @test Displace(x, p, reference, t) isa UnitaryTransform
    end

    @testset "inference" begin
        fock_reference = ω * a' * a + η * (a + a')
        quadrature_reference = (ω / 2) * (x^2 + p^2) + η * x
        @test @inferred(Displace(a, fock_reference)) isa UnitaryTransform
        @test @inferred(Displace(a, fock_reference, t)) isa UnitaryTransform
        @test @inferred(Displace(x, p, quadrature_reference)) isa UnitaryTransform
        @test @inferred(Displace(x, p, quadrature_reference, t)) isa UnitaryTransform
    end
end
