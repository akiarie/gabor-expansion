module Gabor

# Exceptions
struct InvalidFunctionDimensions <: Exception
    funcSize::Int
    spaceSize::Int
end
Base.showerror(io::IO, e::InvalidFunctionDimensions) =
    print(io, "The dimensions ($(e.funcSize)) of the Function values"
          * " do not match the size ($(e.spaceSize)) of the space!")


scalar = Complex{Float64}

struct Space
    size::Int
end

struct Function
    space::Space
    values::Vector{scalar}
end
function Function(Space, values::Vector{scalar}) 
    if length(values) != Space.size
        throw(InvalidFunctionDimensions(length(values), Space.size))
    end
    Function(Space, values)
end

end
