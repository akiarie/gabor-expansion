module Gabor
import Base:*,+,-,/
import Base.size
using LinearAlgebra
using Plots
import Plots.plot
import Plots.surface

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

struct AnalyseMismatch <: Exception
end
Base.showerror(io::IO, e::AnalyseMismatch) =
print(io, "Domain of function being analysed does not match"
      * " elementary function domain")

struct SynthesiseMismatch <: Exception
end
Base.showerror(io::IO, e::SynthesiseMismatch) =
print(io, "Size of coefficients does not match elementary function dimensions")

struct RankDeficientMatrix <: Exception
end
Base.showerror(io::IO, e::RankDeficientMatrix) =
print(io, "The row rank of the G matrix is less than its number of rows")



Scalar = Complex{Float64}

struct Func
    domain::UnitRange{Int64}
    values::Vector{Scalar}
    function Func(domain::UnitRange{Int64}, values::Vector{Scalar})
        if length(values) ≠ length(domain)
            throw(InvalidFuncDimensions(length(values), length(domain)))
        end 
        new(domain, values)
    end
end
(f::Func)(k::Int)::Scalar = circshift(f.values, f.domain[1]-k)[1]
*(A, f::Func) = Func(f.domain, A*f.values)
*(f::Func, A) = Func(f.domain, f.values*A)
+(f::Func, g::Func)::Func = Func(f.domain, g.values+f.values)
-(f::Func, g::Func)::Func = Func(f.domain, g.values-f.values)
/(f::Func, div)::Func = Func(f.domain, f.values/div)
norm(f::Func)::Float64 = sqrt(sum([abs(f(k))^2 for k in f.domain]))
δ(f::Func, g::Func) = norm((f/norm(f)) - (g/norm(g)))
plot(f::Func) = plot(f.domain, map(real, f.values),
                           bar_width=0.2,
                           seriestype=:bar,
                           linestyle=:solid,
                           fillcolor=:black,
                           legend=false)
function compute_all(g::Func, M::Integer, N::Integer)
    L = length(g.domain)
    ψg = ElemFunc(g, Integer(L/M), Integer(L/N))
    S = operator(ψg)
    g̃ = inv(S)*g
    ψg̃ = ElemFunc(g̃, ψg.timeStep, ψg.freqStep)
    ψγ = wr_bio(ψg)
    γ = ψγ.func
    ψg, S, g̃, ψg̃, γ, ψγ
end


# the shift-modulation operator on funcs
function ψ(p::Int, q::Int, g::Func)::Func
    g_shift = Func(g.domain, circshift(g.values, p))
    L = length(g.domain)
    Func(g.domain, [exp(2π*im*q*k/L)*g_shift(k) for k in g.domain])
end

struct ElemFunc
    func::Func
    timeStep::Int
    freqStep::Int
    function ElemFunc(func::Func, timeStep::Int, freqStep::Int)
        sizeDivisors = [length(func.domain) % i for i in [timeStep, freqStep]]
        nonDiv = filter(k -> k ≠ 0, sizeDivisors) # failed to divide evenly
        if (timeStep ≤ 0) || (freqStep ≤ 0) || (length(nonDiv) ≠ 0)
            throw(InvalidElemFuncDimensions(length(func.domain), timeStep, freqStep))
        end
        new(func, timeStep, freqStep)
    end
    ElemFunc(f::Func, ψg::ElemFunc) = ElemFunc(f, ψg.timeStep, ψg.freqStep)
end
(ψg::ElemFunc)(m::Int, n::Int)::Func = ψ(m*ψg.timeStep, n*ψg.freqStep, ψg.func)
dimensions(ψg::ElemFunc) = vcat(map(div -> Int(length(ψg.func.domain)/div), [ψg.timeStep, ψg.freqStep]), length(ψg.func.domain))

function operator(ψg::ElemFunc)
    M, N, _ = dimensions(ψg)
    shift_values = [ψg(m, n).values for m in 0:M-1, n in 0:N-1]
    sum([v*conj(transpose(v)) for v in shift_values])
end

netΔ(A, B) = sum([abs(A[i]-B[i]) for i in 1:length(A)])

# returns the netΔ to the L×L identity of the outer product sum
function biorthogonal(ψg::ElemFunc, ψγ::ElemFunc)
    if (dimensions(ψg) ≠ dimensions(ψγ)) 
        throw(FunctionsDoNotMatch())
    end
    L = length(ψg.func.domain)
    M, N, _ = dimensions(ψg)
    out_prod = sum([ψg(m, n).values*conj(transpose(ψγ(m, n).values)) for m in 0:M-1, n in 0:N-1])
    id_L = Matrix{Scalar}(LinearAlgebra.I, L, L)
    netΔ(id_L, out_prod)
end

function frame(ψg::ElemFunc)
    S = operator(ψg)
    γ = inv(S)*ψg.func
    ψγ = ElemFunc(γ, ψg)
    biorthogonal(ψg, ψγ)
end

struct Lattice
    values::Array{Scalar, 2}
end
(c::Lattice)(m::Int, n::Int) = c.values[m+1, n+1]
size(c::Lattice) = size(c.values)
function surface(c::Lattice, filter=abs)
    M, N = size(c.values)
    mid_M, mid_N = Int(M/2), Int(N/2)
    m_range = circshift(0:N-1, mid_N)
    n_range = circshift(0:M-1, mid_M)

    # determine ticks
    freqs = [100, 10, 5, 2, 1]
    levels = [1000, 50, 25, 15, 1]
    tick_freq(S) = freqs[findfirst(k -> div(S, k) > 0, levels)]
    x_freq = tick_freq(N)
    y_freq = tick_freq(M)
    xticks = (0:x_freq:N-1, -mid_N:x_freq:(mid_N+N))
    yticks = (0:y_freq:M-1, -mid_M:y_freq:(mid_M+M))

    surface(m_range, n_range, map(filter, c.values),
                xticks=xticks, yticks=yticks,
                xlabel="n", ylabel="m",
                fill=(true, cgrad(:grays,[0.0,0.1,0.2,0.5,1.0])),
                legend=:none)
end

function analyse(ψg::ElemFunc, x::Func)
    M, N, _ = dimensions(ψg)
    if (x.domain ≠ ψg.func.domain)
        throw(AnalyseMismatch())
    end
    Σ(m,n) = sum([conj(ψg(m, n)(k))*x(k) for k in ψg.func.domain])
    Lattice([Σ(m,n) for m in 0:M-1, n in 0:N-1])
end

function synthesize(ψg::ElemFunc, c::Lattice)
    M, N, _ = dimensions(ψg)
    if (M,N) ≠ size(c)
        throw(SynthesiseMismatch())
    end
    x(k) = sum([c(m,n)*ψg(m,n)(k) for m in 0:M-1, n in 0:N-1])
    domain = ψg.func.domain
    Func(domain, map(x, domain))
end

# compute Wexler-Raz minimum energy dual
function wr_bio(ψg::ElemFunc)::ElemFunc
    M, N, L = dimensions(ψg)
    lattice = [(m,n) for m in 0:ψg.freqStep-1 for n in 0:ψg.timeStep-1]
    G = [conj(ψ(m*N, n*M, ψg.func)(k)) for (m,n) in lattice, k in ψg.func.domain]
    if rank(G) ≠ ψg.timeStep*ψg.freqStep
        throw(RankDeficientMatrix())
    end
    μ = vcat(L/(M*N), zeros(Scalar, ψg.timeStep*ψg.freqStep-1))
    γ = Func(ψg.func.domain, pinv(G)*μ)
    ElemFunc(γ, ψg)
end

function wr_bio(ψg::ElemFunc, w::Int)
    M, N, L = dimensions(ψg)
    M′ = floor(w/N) == (w/N) ? Int(w/N)-1 : Int(floor(w/N))
    lattice = [(m,n) for m in 0:M′ for n in 0:ψg.timeStep-1]
    G = [conj(ψ(m*N, n*M, ψg.func)(k)) for (m,n) in lattice, k in ψg.func.domain[1:w]]
    μ = vcat(L/(M*N), zeros(Scalar, (M′+1)*ψg.timeStep-1))
    γ_values = G \ μ
    γ = Func(ψg.func.domain, vcat(γ_values, zeros(Scalar, L-w)))
    ElemFunc(γ, ψg)
end

end
