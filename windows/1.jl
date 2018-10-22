# rectangular windows
domain = 0:63

g = Gabor.Func(domain, vcat(zeros(Gabor.Scalar, 28),
                            ones(Gabor.Scalar, 8),
                            zeros(Gabor.Scalar, 28)))

# efficient frames
ψg, S, g̃, ψg̃, γ, ψγ = Gabor.compute_all(g, 8, 8)
ψgᵣ, Sᵣ, g̃ᵣ, ψg̃ᵣ, γᵣ, ψγᵣ = Gabor.compute_all(g, 8, 16)
