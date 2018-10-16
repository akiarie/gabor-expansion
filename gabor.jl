module Gabor

using LinearAlgebra

# Exceptions
struct InvalidFuncDimensions <: Exception
    funcSize::Int
    spaceSize::Int
end
Base.showerror(io::IO, e::InvalidFuncDimensions) =
    print(io, "The dimensions ($(e.funcSize)) of the Func values"
          * " do not match the size ($(e.spaceSize)) of the space!")

struct InvalidElemFuncDimensions <: Exception
    spaceSize::Int
    timeStep::Int
    freqStep::Int
end
Base.showerror(io::IO, e::InvalidElemFuncDimensions) =
    print(io, "Unable to fit time-step $(e.timeStep) and frequency-step"
         * " of $(e.freqStep) in space of size $(e.spaceSize)!")

struct FunctionsDoNotMatch <: Exception
end
Base.showerror(io::IO, e::FunctionsDoNotMatch) =
    print(io, "Cannot compare functions with different dimensions or"
         * " domains for biorthogonality")


Scalar = Complex{Float64}

struct Space
    domain::UnitRange{Int64}
end
size(space::Space)::Int = length(space.domain)

struct Func
    space::Space
    values::Vector{Scalar}
    function Func(space::Space, values::Vector{Scalar})
        if length(values) ≠ size(space)
            throw(InvalidFuncDimensions(length(values), size(space)))
        end 
        new(space, values)
    end
end
(f::Func)(k::Int)::Scalar = f.values[findfirst(isequal(k), f.space.domain)]

# the shift-modulation operator on funcs
function ψ(p::Int, q::Int, g::Func)::Func
    g_shift = circshift(g.values, p)
    L = size(g.space)
    Func(g.space, [exp(2π*im*q*k/L)*g_shift[k+1] for k in g.space.domain])
end

struct ElemFunc
    func::Func
    timeStep::Int
    freqStep::Int
    function ElemFunc(func::Func, timeStep::Int, freqStep::Int)
        sizeDivisors = [size(func.space) % i for i in [timeStep, freqStep]]
        nonDiv = filter(k -> k ≠ 0, sizeDivisors) # failed to divide evenly
        if (timeStep ≤ 0) || (freqStep ≤ 0) || (length(nonDiv) ≠ 0)
            throw(InvalidElemFuncDimensions(func.space.size, timeStep, freqStep))
        end
        new(func, timeStep, freqStep)
    end
end
(ψg::ElemFunc)(m::Int, n::Int)::Func = ψ(m*ψg.timeStep, n*ψg.freqStep, ψg.func)
dimensions(ψg::ElemFunc) = map(div -> Int(size(ψg.func.space)/div), [ψg.timeStep, ψg.freqStep])

function operator(ψg::ElemFunc)
    M, N = dimensions(ψg)
    shift_values = [ψg(m, n).values for m in 0:M-1, n in 0:N-1]
    sum([shift_g*transpose(shift_g) for shift_g in shift_values])
end

netΔ(A, B) = sum([abs(A[i]-B[i]) for i in 1:length(A)])

# returns the netΔ to the L×L identity of the outer product sum
function biorthogonal(ψg::ElemFunc, ψγ::ElemFunc)
    if (ψg.timeStep ≠ ψγ.timeStep) || (ψg.freqStep ≠ ψγ.timeStep) || (ψg.func.space ≠ ψγ.func.space)
        throw(FunctionsDoNotMatch())
    end
    L = size(ψg.func.space)
    M, N = dimensions(ψg)
    out_prod = sum([ψg(m, n).values*conj(transpose(ψγ(m, n).values)) for m in 0:M-1, n in 0:N-1])
    id_L = Matrix{Scalar}(LinearAlgebra.I, L, L)
    netΔ(id_L, out_prod)
end

function frame(ψg::ElemFunc)
    S = operator(ψg)
    γ = Func(ψg.func.space, inv(S)*ψg.func.values)
    ψγ = ElemFunc(γ, ψg.timeStep, ψg.freqStep)
    biorthogonal(ψg, ψγ)
end

# compute Wexler-Raz minimum energy dual
function wr_bio(ψg::ElemFunc)::ElemFunc
    L = size(ψg.func.space)
    M, N = dimensions(ψg)
    W = exp(2π*im/N)
    ρ = M*N/L
    μ = vcat(1, zeros(Scalar, L-1))
    lattice = [(m,n) for m in 0:M-1 for n in 0:N-1]
    G = [ρ*conj(ψg(m, n)(k)) for (m,n) in lattice, k in ψg.func.space.domain]
    γ = Func(ψg.func.space, G \ μ)
    ElemFunc(γ, ψg.timeStep, ψg.freqStep)
end

end


L = 64
periodic(k::Integer, period) = ((k%period)+period)%period
window(k::Integer)::Complex{Float64} = (2^(1/2)/8)^(1/2)*exp(-π*((periodic(k,L)-31.5)/8)^2)

space = Gabor.Space(0:63)
g = Gabor.Func(space, map(window, 0:63))
g_elem = Gabor.ElemFunc(g, 8, 8)
S = Gabor.operator(g_elem)
g̃ = Gabor.Func(space, inv(S)*g.values)
g̃_elem = Gabor.ElemFunc(g̃, 8, 8)

γ_elem = Gabor.wr_bio(g_elem)
γ = γ_elem.func
