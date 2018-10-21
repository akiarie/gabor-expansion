# gaussian windows
domain = 0:63

gaussian(k::Integer)::Gabor.Scalar = (2^(1/2)*8)^(1/2)*exp(-π*((k-31.5)/8)^2)

g = Gabor.Func(domain, map(gaussian, domain))

# efficient frames
ψg₁, S₁, g̃₁, ψg̃₁, γ₁, ψγ₁ = Gabor.compute_all(g, 8, 8)
ψg₂, S₂, g̃₂, ψg̃₂, γ₂, ψγ₂ = Gabor.compute_all(g, 4, 16)
ψg₃, S₃, g̃₃, ψg̃₃, γ₃, ψγ₃ = Gabor.compute_all(g, 16, 4)

# redundant
ψgᵣ₁, Sᵣ₁, g̃ᵣ₁, ψg̃ᵣ₁, γᵣ₁, ψγᵣ₁ = Gabor.compute_all(g, 16, 8)
ψgᵣ₂, Sᵣ₂, g̃ᵣ₂, ψg̃ᵣ₂, γᵣ₂, ψγᵣ₂ = Gabor.compute_all(g, 8, 16)
ψgᵣ₃, Sᵣ₃, g̃ᵣ₃, ψg̃ᵣ₃, γᵣ₃, ψγᵣ₃ = Gabor.compute_all(g, 16, 16)
ψgᵣ₄, Sᵣ₄, g̃ᵣ₄, ψg̃ᵣ₄, γᵣ₄, ψγᵣ₄ = Gabor.compute_all(g, 32, 32)
