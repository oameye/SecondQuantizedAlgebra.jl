# Changelog

All notable changes to [`SecondQuantizedAlgebra.jl`](https://github.com/qojulia/SecondQuantizedAlgebra.jl) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


##  [v0.11.1]

### Added

- Extend exact unitary transformations with a shared affine representation, general bosonic Bogoliubov maps, generator-derived exact transforms and rotating frames, and Hamiltonian-derived displacement frames.

### Fixed

- Make numeric conversion for QuantumToolbox.jl type-stable
- Follow the QuantumInterface 0.4.4 tensor-product ownership change by importing `tensor` and `⊗` directly from TensorCore.

##  [v0.11.0]

### Fixed

- Improve compact expim phases through simplification, substitution, differentiation, unitary transformations, printing, and averaging.

## [v0.10.1]

### Fixed

- `dagger(::QField)` now extends `QuantumInterface.dagger`, providing the shared QuantumOptics-ecosystem adjoint verb for operator expressions. It forwards to `qadjoint`; use `qadjoint` or `qconj` for scalar and symbolic adjoints.

## [v0.10.1]

### Added

- Exact named unitary transformations for Fock, phase-space, spin/Pauli, and ordinary N-level operators, with `conjugate`, time-aware `transform`, inversion, composition, and analytically derived gauges.
- `expim(x)` provides a canonical exact unit-phase algebra: products and integer powers add their real arguments, inverse phases cancel, and substitution, differentiation, projections, and explicit exponential/trigonometric conversion preserve the phase representation.
- Coefficient simplification reduces exact trigonometric and hyperbolic identities, including supported composite arguments and higher powers.
- A unitary-transformation guide and resonator-rotation example show how to prepare an exact moving-frame Hamiltonian for a later Floquet expansion.


## [v0.10.0]

Numeric conversion (`to_numeric`/`numeric_average`/`expect`) was redesigned to be extensible, type-stable, and multi-backend. This is a breaking release.

### Added

- A second numeric backend, [QuantumToolbox.jl](https://github.com/qutip/QuantumToolbox.jl), alongside QuantumOpticsBase.jl. Both are wired in through Julia package extensions and selected with the backend singletons `QuantumOpticsBackend()` / `QuantumToolboxBackend()`.
- A uniform Hilbert-space entry point `to_numeric(op, h::HilbertSpace, dims; backend, parameter, time_parameter, operators, adjoint_ops, op_type)`. It builds the backend basis from `dims` (Fock cutoff / spin number) and is the only form that works for both backends. The backend defaults to the single loaded one.
- Open backend hooks `numeric_operator`, `numeric_basis`, `numeric_subbasis`, `numeric_embed`, `numeric_identity`, `numeric_num_subsystems`, `numeric_assemble`, `numeric_assemble_td`, `numeric_materialize`, `numeric_expect`, and `numeric_backend` are exported. Downstream packages can implement another backend and add numeric support for existing operator roles; `OpKind` remains a closed symbolic-role set.
- Differentiable control: a `time_parameter` value may be a `(p, t) -> value` function of the solver parameter vector. On QuantumToolbox this yields a `QobjEvo` differentiable with respect to `p` (SciMLSensitivity with Enzyme/Mooncake), enabling gradient-based optimal control via `sesolve`/`mesolve`; QuantumOptics rejects the form. See the Kerr-resonator control example.

- Some TTFX motivated changes and more complicated precompile workload. [#223](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/223)

### Changed (breaking)

- `QuantumOpticsBase` moved from a hard dependency to a weak dependency. Using the numeric API now requires loading a backend: add `using QuantumOpticsBase` (or `using QuantumToolbox`) next to `using SecondQuantizedAlgebra`. The lightweight `QuantumInterface.jl` is a new hard dependency (it owns the `⊗`/`tensor`/`expect`/`basis` generics that the algebra extends).
- The time-dependent form (`time_parameter` non-empty) returns the backend's **native** time-dependent operator: a `TimeDependentSum` (QuantumOptics) or a `QobjEvo` (QuantumToolbox), both directly consumable by `mesolve`/`master_dynamic`/`sesolve`. It is no longer a `t -> op(t)` closure.
- Time-dependent conversion accepts only `op_type=nothing` or `identity`; eager `op_type` materializers are rejected because a time-varying operator cannot be materialized once during conversion.
- Product-space `dims` must have exactly one entry per symbolic subspace. Indexed numeric layouts validate that every physical slot is unique and in range instead of silently collapsing multiple sites onto a simple basis.

### Changed

- Averages of provably Hermitian operators (`adjoint(A) == A`) now carry `symtype === Real` instead of `Number`. This gives a faster `simplify` path and makes `conj(⟨a'a⟩)` fold to `⟨a'a⟩` rather than an inert `conj(...)` wrapper; indexed sums and lifted time-dependent variables inherit the typing, which survives `substitute`. Resolves [#171](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/171).

### Fixed

- An elementary function of a literal zero left unevaluated by Symbolics (e.g. the `exp(0)` factor produced when `exp(im*ω*t)` is Euler-expanded to `cos(ω t)*exp(0) + exp(0)*im*sin(ω t)`) is now folded to its exact value in the coefficient algebra, so it no longer leaks into printed equations. Only exact identities at argument `0` fold (`exp/cos/cosh → 1`, `sin/tan/sinh/tanh → 0`); non-zero or non-constant arguments (`exp(2)`, `sin(π)`) stay symbolic.

- `commutator` on two operator leaves no longer disagrees with `*`. Its fast path gated on a same-site test that ignored the operator name, so two distinct modes sharing one `FockSpace` (`a` and `b` from `Destroy(h, :a)`, `Destroy(h, :b)`) returned the `[a,a†] = 1` residual while `a*b' - b'*a` correctly gave `0`. The gate is now the `Equal` site comparison the per-operator hooks are documented to require.


## [v0.9.4]

### Fixed

- LaTeX rendering of an index whose name carries numeric per-slot suffixes no longer produces an invalid double subscript. As the pipeline transforms an index it accumulates suffixes through `(i::Index)(k)` (for example an index reaches `i_2_1` when a second distinct atom `i(2)` is later collapsed to its position-1 representative), and `_latex_index_suffix` emitted the name verbatim as `_{i_2_1}`. The unbraced `_2_1` is a double subscript that MathJax rejects with "Double subscripts: use braces to clarify", so the whole equation was left unrendered in Documenter pages and notebooks. The slot numbers are now joined into a single comma subscript (`i_2_1` renders as `i_{2,1}`); bare index names are unchanged.


### Added

- `CollectiveNLevelSpace` and `CollectiveTransition` provide an exact collective description of ``N`` identical multilevel systems in the permutation-symmetric subspace. A collective transition represents ``S^{ij} = \sum_k |i\rangle_k\langle j|`` and obeys the closed ``\mathfrak{su}(n)`` algebra ``[S^{ij}, S^{kl}] = \delta_{jk}S^{il} - \delta_{li}S^{kj}``, including the two-term case ``[S^{12}, S^{21}] = S^{11} - S^{22}``. Integer and symbolic levels, simple and product-space constructors, adjoints, fundamental-operator generation, Unicode/LaTeX printing, and the `is_collective_transition` predicate are supported. Collective transitions deliberately cannot be site-indexed and do not use the single-site completeness expansion: ``\sum_i S^{ii} = N I``, rather than ``I``.
- `to_numeric` converts a `CollectiveTransition` on a `QuantumOpticsBase.ManyBodyBasis` with an `NLevelBasis` one-body basis through `manybodyoperator`. This gives exact finite-``N`` dynamics in a symmetric space of dimension ``\binom{N+n-1}{n-1}`` and works unchanged inside composite bases. For two levels it reproduces the spin-``N/2`` representation, with ``S^{12}=S_+``, ``S^{21}=S_-``, and ``(S^{11}-S^{22})/2=S_z`` in the QuantumOpticsBase convention.
- `OP_COLLECTIVE_TRANSITION` extends the public `OpKind` tags with value `7`. It is appended after `OP_MOMENTUM`, preserving the values `0:6` and structural ordering keys of every existing operator role.
- A public constructor and body accessor for the moment-layer indexed-sum node, completing the read/write API alongside the existing `is_indexed_sum`, `has_sum_metadata`, `get_sum_indices`, and `get_sum_non_equal`. `SecondQuantizedAlgebra.indexed_sum(body, indices; non_equal = NonEqualPair[])` wraps an averaged `body` in a sum node over `indices` (mirroring what `average` builds for a summed `QAdd`), and `SecondQuantizedAlgebra.get_sum_body(x)` reads the summed body back out. This lets downstream code (for example cumulant expansion factorizing a summed moment into a product form) re-wrap a body while preserving its summation scope. Resolves [#209](https://github.com/qojulia/SecondQuantizedAlgebra.jl/issues/209).
- Public accessors and rebuilders for `Op`/`Index` internals, so downstream packages no longer reach into struct fields or the private constructors. `operator_index(op)` reads an operator's site `Index`; `acts_on(idx::Index)` returns an index's subspace as a one-element `Vector{Int}` (matching the return type of `acts_on(::QSym)`); `set_acts_on(op, aon)` rebinds an operator to a new subspace; and `SecondQuantizedAlgebra.rename(x, name)` renames an `Op` or an `Index` in place of the internal `Op`/`Index` constructor (preserving every other field). This is the SecondQuantizedAlgebra half of QuantumCumulants issue [#300](https://github.com/qojulia/QuantumCumulants.jl/issues/300).

### Changed

- The internal `_commute_pair` hook now carries up to two residual terms. Fock, phase-space, and spin operators continue to use one residual slot; collective transitions use both when both Kronecker deltas contribute. The eager commute pass and the leaf `commutator` fast path retain and combine both branches.
- Passing a `String` where a name `Symbol` is expected now throws a clear `ArgumentError` telling the user to use a `Symbol` (`:x`), instead of a cryptic `MethodError`. This covers Hilbert-space constructors (`FockSpace`, `NLevelSpace`, `PauliSpace`, `SpinSpace`, `PhaseSpace`), operator constructors (`Destroy`, `Create`, `Transition`, `Pauli`, `Spin`, `Position`, `Momentum`), and the index-related constructors `Index`, `IndexedVariable`, and `DoubleIndexedVariable`. Names remain `Symbol`s by design: a single canonical name type keeps name comparison and interning type-stable, so strings are rejected rather than silently converted.

### Documentation

- Added a guide comparing indexed sums with exact collective symmetric-subspace operators, developer documentation for the two-residual commute contract and the collective-transition/spin relationship, and a finite-``N`` Tavis–Cummings example using `ManyBodyBasis`

## [v0.9.3]

### Fixed

- The diagonal split in `_accumulate_with_diag!` no longer leaks the diagonal coefficient into the off-diagonal body when a summation index appears only in a term's coefficient. Multiplying a sum such as `Σ_k u(l,k)·σ_l` by an operator peels off the diagonal `k = l` slice, but the substitution dropped the split's `k ≠ l` constraint instead of rewriting it onto the collapsed index. A sibling scope then re-admitted `k = l`, so the off-diagonal body's coefficient became `u(l,k) + u(l,l)` and the diagonal was over-counted (for example `commutator(σ_l^{22}, Σ_{j,k} u(j,k)(σ_j^{21}+σ_j^{12}))` gained a spurious `u(l,l)` in every `Σ_k` term). The diagonal contribution now rewrites each `sum_idx ≠ β` constraint onto `ext_idx`, exactly as the `Σ` diagonal split already does, so the off-diagonal body keeps its bare coefficient and the diagonal is emitted once. This is the SecondQuantizedAlgebra half of QuantumCumulants issue [#198](https://github.com/qojulia/QuantumCumulants.jl/issues/198).

## [v0.9.2]

### Changed

- `to_numeric` now returns a concrete operator by default, and its representation depends only on `op_type`, not on the shape of the expression. Previously a bare operator translated to a `LazyTensor` and a sum to a `LazySum` (via the positional path), while the keyword path materialized `sparse` only in some branches, so composing subexpressions silently changed the return type and operators handed to solvers were lazy by default. The default `op_type` is now `sparse` (was `identity`), applied uniformly to bare operators, products, and sums. Pass `op_type=dense` for a dense operator, or `op_type=identity` to opt into the natural lazy representation (`LazyTensor`/`LazyProduct`/`LazySum`) for large tensor-product spaces where an operator is local to a few subsystems. The lazy embedding is kept internal, so a term is still assembled from lazy factors and materialized exactly once. Code relying on the previous lazy default must now pass `op_type=identity`.

## [v0.9.1]

### Fixed

- Printing a sum of averages whose operators have different concrete types (for example `⟨a'a⟩`, which wraps a `QAdd`, alongside `⟨σ⟩`, which wraps a `QSym`) no longer errors. When SymbolicUtils sorts such a sum for display it unwraps the operator constants and falls back to comparing the operator *types*, which previously threw `MethodError: no method matching isless(::Type{QAdd}, ::Type{...})`. `QField` types now compare by name, so heterogeneous operator sums render as before.

## [v0.9.0]

### Added

- `substitute` now replaces operators as well as scalar parameters in a single, simultaneous pass. A rule dictionary is split by key: `QSym` keys replace operator leaves, every other key substitutes into coefficients through SymbolicUtils. Operator replacements are applied once to the original leaves and the replacement expression is not searched again for operator keys, so a mode transformation such as `a => g*a + h*b` is well-defined even though its right-hand side contains the key `a`. Missing adjoint rules are generated automatically (`a => b` also implies `a' => b'`); pass `replace_adjoint=false` to match only the keys supplied explicitly.
- A keyword form of `to_numeric(op, basis; parameter, time_parameter, operators, adjoint_ops, op_type)`. Scalar `parameter`s are substituted first, then the expression is translated with custom numeric `operators` (missing adjoints are added when `adjoint_ops=true`), and each emitted operator is passed through `op_type`. When `time_parameter` is non-empty the result is a closure `t -> op(t)`; its values may be numbers or functions of time, and a key may be a bare variable `v` or `conj(v)`. A vector form `to_numeric(ops::AbstractVector, basis; kwargs...)` forwards the keywords to each element.
- `to_num(c::Coeff)` is now public and documented as the supported way to read the coefficient returned when iterating a `QAdd`. `Coeff` and `CNum` are now declared public, and the package's public-API declarations use [SciMLPublic.jl](https://github.com/SciML/SciMLPublic.jl) instead of a hand-rolled `@public` macro.

### Changed

- The ordering of the displayed expressions are now ordered first by Hilbert space and then by operator kind.

## [v0.8.3]

### Changed

- `average` of a many-term operator is now O(n) instead of folding with Base.:+? no ...
