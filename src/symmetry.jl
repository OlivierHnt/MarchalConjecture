import RadiiPolynomial: SymBaseSpace, indices, _findposition, _iscompatible

abstract type MarchalSymmetry <: SymBaseSpace end

struct EvenCos{T<:Real} <: MarchalSymmetry
    space :: Fourier{T}
    function EvenCos(space::Fourier{T}) where {T<:Real}
        ord = order(space)
        return new{T}(Fourier(ord-isodd(ord), frequency(space)))
    end
end
EvenCos(order::Int, frequency::Real) = EvenCos(Fourier(order, frequency))

struct OddCos{T<:Real} <: MarchalSymmetry
    space :: Fourier{T}
    function OddCos(space::Fourier{T}) where {T<:Real}
        ord = order(space)
        return new{T}(Fourier(ord-iseven(ord), frequency(space)))
    end
end
OddCos(order::Int, frequency::Real) = OddCos(Fourier(order, frequency))

struct EvenSin{T<:Real} <: MarchalSymmetry
    space :: Fourier{T}
    function EvenSin(space::Fourier{T}) where {T<:Real}
        ord = order(space)
        return new{T}(Fourier(ord-isodd(ord), frequency(space)))
    end
end
EvenSin(order::Int, frequency::Real) = EvenSin(Fourier(order, frequency))

struct OddSin{T<:Real} <: MarchalSymmetry
    space :: Fourier{T}
    function OddSin(space::Fourier{T}) where {T<:Real}
        ord = order(space)
        return new{T}(Fourier(ord-iseven(ord), frequency(space)))
    end
end
OddSin(order::Int, frequency::Real) = OddSin(Fourier(order, frequency))

struct EvenFourier{T<:Real} <: MarchalSymmetry
    space :: Fourier{T}
    function EvenFourier(space::Fourier{T}) where {T<:Real}
        ord = order(space)
        return new{T}(Fourier(ord-isodd(ord), frequency(space)))
    end
end
EvenFourier(order::Int, frequency::Real) = EvenFourier(Fourier(order, frequency))

#

Base.intersect(s₁::EvenCos, s₂::EvenCos) = EvenCos(intersect(desymmetrize(s₁), desymmetrize(s₂)))
Base.intersect(s₁::OddCos, s₂::OddCos) = OddCos(intersect(desymmetrize(s₁), desymmetrize(s₂)))
Base.intersect(s₁::EvenSin, s₂::EvenSin) = EvenSin(intersect(desymmetrize(s₁), desymmetrize(s₂)))
Base.intersect(s₁::OddSin, s₂::OddSin) = OddSin(intersect(desymmetrize(s₁), desymmetrize(s₂)))
Base.intersect(s₁::EvenFourier, s₂::EvenFourier) = EvenFourier(intersect(desymmetrize(s₁), desymmetrize(s₂)))

#

indices(s::EvenCos) = 0:2:order(s)
indices(s::OddCos) = 1:2:order(s)
indices(s::EvenSin) = 2:2:order(s)
indices(s::OddSin) = 1:2:order(s)
indices(s::EvenFourier) = -order(s):2:order(s)

_findposition(i::Int, ::EvenCos) = i ÷ 2 + 1
_findposition(i::Int, ::OddCos) = (i + 1) ÷ 2
_findposition(i::Int, ::EvenSin) = i ÷ 2
_findposition(i::Int, ::OddSin) = (i + 1) ÷ 2
_findposition(i::Int, s::EvenFourier) = (i + order(s)) ÷ 2 + 1

_iscompatible(s₁::EvenCos, s₂::EvenCos) = _iscompatible(desymmetrize(s₁), desymmetrize(s₂))
_iscompatible(s₁::OddCos, s₂::OddCos) = _iscompatible(desymmetrize(s₁), desymmetrize(s₂))
_iscompatible(s₁::EvenSin, s₂::EvenSin) = _iscompatible(desymmetrize(s₁), desymmetrize(s₂))
_iscompatible(s₁::OddSin, s₂::OddSin) = _iscompatible(desymmetrize(s₁), desymmetrize(s₂))
_iscompatible(s₁::EvenFourier, s₂::EvenFourier) = _iscompatible(desymmetrize(s₁), desymmetrize(s₂))





#- sym to no sym

function materialize_symmetry(s_sym::Union{EvenCos,OddCos}, s::Fourier)
    A = zeros(ExactReal{Int}, s_sym, s)
    for β ∈ indices(s_sym), α ∈ indices(s)
        if β == abs(α)
            A[α,β] = 1
        end
    end
    return coefficients(A)
end

function materialize_symmetry(s_sym::Union{EvenSin,OddSin}, s::Fourier)
    A = zeros(Complex{ExactReal{Int}}, s_sym, s)
    for β ∈ indices(s_sym), α ∈ indices(s)
        if β == abs(α)
            A[α,β] = -im*sign(α)
        end
    end
    return coefficients(A)
end

function materialize_symmetry(s_sym::EvenFourier, s::Fourier)
    A = zeros(ExactReal{Int}, s_sym, s)
    for β ∈ indices(s_sym), α ∈ indices(s)
        if β == α
            A[α,β] = 1
        end
    end
    return coefficients(A)
end

#- no sym to sym

function materialize_symmetry(s::Fourier, s_sym::Union{EvenCos,OddCos})
    A = zeros(ExactReal{Float64}, s, s_sym)
    for β ∈ indices(s), α ∈ indices(s_sym)
        if abs(β) == α
            A[α,β] = ifelse(β == 0, 1.0, 0.5)
        end
    end
    return coefficients(A)
end

function materialize_symmetry(s::Fourier, s_sym::Union{EvenSin,OddSin})
    A = zeros(Complex{ExactReal{Float64}}, s, s_sym)
    for β ∈ indices(s), α ∈ indices(s_sym)
        if abs(β) == α
            A[α,β] = 0.5im*sign(β)
        end
    end
    return coefficients(A)
end

function materialize_symmetry(s::Fourier, s_sym::EvenFourier)
    A = zeros(ExactReal{Int}, s, s_sym)
    for β ∈ indices(s), α ∈ indices(s_sym)
        if β == α
            A[α,β] = 1
        end
    end
    return coefficients(A)
end





#

function project_sym!(c::Sequence{<:CartesianSpace}, a::Sequence{<:CartesianSpace})
    space_c = space(c)
    # if space(a) == space_c
    #     coefficients(c) .= coefficients(a)
    # else
    for i ∈ 1:nspaces(space_c)
        project_sym!(component(c, i), component(a, i))
    end
    # end
    return c
end

function project_sym!(c::Sequence{ParameterSpace}, a::Sequence{ParameterSpace})
    c[1] = a[1]
    return c
end

function project_sym!(c::Sequence{<:MarchalSymmetry}, a::Sequence{<:Fourier})
    mul!(coefficients(c), materialize_symmetry(space(a), space(c)), coefficients(a))
    return c
end

function project_sym!(c::Sequence{<:Fourier}, a::Sequence{<:MarchalSymmetry})
    mul!(coefficients(c), materialize_symmetry(space(a), space(c)), coefficients(a))
    return c
end

#

function project_sym!(C::LinearOperator{<:CartesianSpace,<:CartesianSpace}, A::LinearOperator{<:CartesianSpace,<:CartesianSpace})
    domain_C = domain(C)
    codomain_C = codomain(C)
    # if domain(A) == domain_C && codomain(A) == codomain_C
    #     coefficients(C) .= coefficients(A)
    # else
    for j ∈ 1:nspaces(domain_C), i ∈ 1:nspaces(codomain_C)
        project_sym!(component(C, i, j), component(A, i, j))
    end
    # end
    return C
end

function project_sym!(C::LinearOperator{<:CartesianSpace,<:VectorSpace}, A::LinearOperator{<:CartesianSpace,<:VectorSpace})
    domain_C = domain(C)
    # codomain_C = codomain(C)
    # if domain(A) == domain_C && codomain(A) == codomain_C
    #     coefficients(C) .= coefficients(A)
    # else
    for j ∈ 1:nspaces(domain_C)
        project_sym!(component(C, j), component(A, j))
    end
    # end
    return C
end

function project_sym!(C::LinearOperator{<:VectorSpace,<:CartesianSpace}, A::LinearOperator{<:VectorSpace,<:CartesianSpace})
    # domain_C = domain(C)
    codomain_C = codomain(C)
    # if domain(A) == domain_C && codomain(A) == codomain_C
    #     coefficients(C) .= coefficients(A)
    # else
    for i ∈ 1:nspaces(codomain_C)
        project_sym!(component(C, i), component(A, i))
    end
    # end
    return C
end

function project_sym!(C::LinearOperator{ParameterSpace,ParameterSpace}, A::LinearOperator{ParameterSpace,ParameterSpace})
    C[1,1] = A[1,1]
    return C
end

function project_sym!(C::LinearOperator{<:MarchalSymmetry,<:MarchalSymmetry}, A::LinearOperator{<:Fourier,<:Fourier})
    mul!(coefficients(C), materialize_symmetry(codomain(A), codomain(C)), coefficients(A) * materialize_symmetry(domain(C), domain(A)))
    return C
end

function project_sym!(C::LinearOperator{<:MarchalSymmetry,ParameterSpace}, A::LinearOperator{<:Fourier,ParameterSpace})
    mul!(coefficients(C), coefficients(A), materialize_symmetry(domain(C), domain(A)))
    return C
end

function project_sym!(C::LinearOperator{ParameterSpace,<:MarchalSymmetry}, A::LinearOperator{ParameterSpace,<:Fourier})
    mul!(coefficients(C), materialize_symmetry(codomain(A), codomain(C)), coefficients(A))
    return C
end

function project_sym!(C::LinearOperator{<:Fourier,ParameterSpace}, A::LinearOperator{<:MarchalSymmetry,ParameterSpace})
    mul!(coefficients(C), LinearOperator(coefficients(A)), materialize_symmetry(domain(C), domain(A)))
    return C
end

function project_sym!(C::LinearOperator{ParameterSpace,<:Fourier}, A::LinearOperator{ParameterSpace,<:MarchalSymmetry})
    mul!(coefficients(C), materialize_symmetry(codomain(A), codomain(C)), coefficients(A))
    return C
end





#

function norm1_sym(c::Sequence{<:CartesianSpace}, ν)
    # @assert isreal(ν) && abs(ν) > 1
    z = abs(zero(eltype(c))) * zero(ν)
    for i ∈ 1:nspaces(space(c))
        z += norm1_sym(component(c, i), ν)
    end
    return z
end

function norm1_sym(c::Sequence{<:MarchalSymmetry}, ν)
    z = abs(zero(eltype(c))) * zero(ν)
    for k ∈ indices(space(c))
        z += ExactReal(ifelse(k == 0, 1, 2)) * abs(c[k]) * ν ^ ExactReal(abs(k))
    end
    return z
end

function norm1_sym(c::Sequence{<:EvenFourier}, ν)
    z = abs(zero(eltype(c))) * zero(ν)
    for k ∈ indices(space(c))
        z += abs(c[k]) * ν ^ ExactReal(abs(k))
    end
    return z
end

norm1_sym(c::Sequence{ParameterSpace}, ::Any) = abs(c[1])

function opnorm1_sym(C::LinearOperator{<:CartesianSpace,<:CartesianSpace}, ν)
    # @assert isreal(ν) && abs(ν) > 1
    z = zero(eltype(C)) * zero(ν)
    v = Vector{typeof(z)}(undef, size(coefficients(C), 2))
    for j ∈ axes(coefficients(C), 2)
        c = Sequence(codomain(C), view(C, :, j))
        v[j] = norm1_sym(c, ν)
    end
    return dual_norm1_sym(Sequence(domain(C), v), ν)
end


function dual_norm1_sym(c::Sequence{<:CartesianSpace}, ν)
    z = zero(eltype(c)) * zero(inv(ν))
    for i ∈ 1:nspaces(space(c))
        z = max(z, dual_norm1_sym(component(c, i), ν))
    end
    return z
end

function dual_norm1_sym(c::Sequence{<:MarchalSymmetry}, ν)
    ν⁻¹ = inv(ν)
    z = abs(zero(eltype(c))) * zero(ν⁻¹)
    for k ∈ indices(space(c))
        z = max(z, ExactReal(ifelse(k == 0, 1, 2)) \ abs(c[k]) * ν⁻¹ ^ ExactReal(abs(k)))
    end
    return z
end

function dual_norm1_sym(c::Sequence{<:EvenFourier}, ν)
    ν⁻¹ = inv(ν)
    z = abs(zero(eltype(c))) * zero(ν⁻¹)
    for k ∈ indices(space(c))
        z = max(z, abs(c[k]) * ν⁻¹ ^ ExactReal(abs(k)))
    end
    return z
end

dual_norm1_sym(c::Sequence{ParameterSpace}, ::Any) = abs(c[1])
