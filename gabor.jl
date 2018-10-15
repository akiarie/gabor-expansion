module Gabor

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


Scalar = Complex{Float64}

struct Space
    domain::UnitRange{Int64}
end
size(space::Space) = length(space.domain)

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

# the shift-modulation operator on funcs
function ψ(p::Int, q::Int, g::Func)
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

function operator(g::ElemFunc)
    M, N = map(div -> Int(size(g.func.space)/div), [g.timeStep, g.freqStep])
    shift_values = [ψ(m*g.timeStep, n*g.freqStep, g.func).values for m in 0:M-1, n in 0:N-1]
    sum([ψg*transpose(ψg) for ψg in shift_values])
end

end
