# gaussian windows
domain = 0:63

x = Gabor.Func(domain, vcat(zeros(Gabor.Scalar, 23),
                            ones(Gabor.Scalar, 18),
                            zeros(Gabor.Scalar, 23)))

rect = Gabor.Func(domain, vcat(zeros(Gabor.Scalar, 28),
                            ones(Gabor.Scalar, 8),
                            zeros(Gabor.Scalar, 28)))

gaussian(k::Integer)::Gabor.Scalar = (2^(1/2)*8)^(1/2)*exp(-π*((k-31.5)/8)^2)

gauss = Gabor.Func(domain, map(gaussian, domain))

# efficient frames
ψg₁, S₁, g̃₁, ψg̃₁, γ₁, ψγ₁ = Gabor.compute_all(gauss, 8, 8)
ψg₂, S₂, g̃₂, ψg̃₂, γ₂, ψγ₂ = Gabor.compute_all(gauss, 4, 16)
ψg₃, S₃, g̃₃, ψg̃₃, γ₃, ψγ₃ = Gabor.compute_all(gauss, 16, 4)

# redundant
ψgᵣ₁, Sᵣ₁, g̃ᵣ₁, ψg̃ᵣ₁, γᵣ₁, ψγᵣ₁ = Gabor.compute_all(gauss, 16, 8)
ψgᵣ₂, Sᵣ₂, g̃ᵣ₂, ψg̃ᵣ₂, γᵣ₂, ψγᵣ₂ = Gabor.compute_all(gauss, 8, 16)

ψgᵣ₃, Sᵣ₃, g̃ᵣ₃, ψg̃ᵣ₃, γᵣ₃, ψγᵣ₃ = Gabor.compute_all(gauss, 16, 16)
ψgᵣ₄, Sᵣ₄, g̃ᵣ₄, ψg̃ᵣ₄, γᵣ₄, ψγᵣ₄ = Gabor.compute_all(gauss, 32, 32)
ψgᵣ₅, Sᵣ₅, g̃ᵣ₅, ψg̃ᵣ₅, γᵣ₅, ψγᵣ₅ = Gabor.compute_all(gauss, 64, 64)

ψgᵣ₆, Sᵣ₆, g̃ᵣ₆, ψg̃ᵣ₆, γᵣ₆, ψγᵣ₆ = Gabor.compute_all(gauss, 32, 16)
ψgᵣ₇, Sᵣ₇, g̃ᵣ₇, ψg̃ᵣ₇, γᵣ₇, ψγᵣ₇ = Gabor.compute_all(gauss, 16, 32)

ψgᵣ₈, Sᵣ₈, g̃ᵣ₈, ψg̃ᵣ₈, γᵣ₈, ψγᵣ₈ = Gabor.compute_all(gauss, 64, 16)
ψgᵣ₉, Sᵣ₉, g̃ᵣ₉, ψg̃ᵣ₉, γᵣ₉, ψγᵣ₉ = Gabor.compute_all(gauss, 16, 64)

ψgᵣ₁₀, Sᵣ₁₀, g̃ᵣ₁₀, ψg̃ᵣ₁₀, γᵣ₁₀, ψγᵣ₁₀ = Gabor.compute_all(gauss, 64, 32)
ψgᵣ₁₁, Sᵣ₁₁, g̃ᵣ₁₁, ψg̃ᵣ₁₁, γᵣ₁₁, ψγᵣ₁₁ = Gabor.compute_all(gauss, 32, 64)


ψgrect, Srect, g̃rect, ψg̃rect, γrect, ψγrect = Gabor.compute_all(rect, 8, 8)

c_rect = Gabor.analyse(ψγrect, x)

c_gauss₁ = Gabor.analyse(ψγ₁, x)
c_gauss₂ = Gabor.analyse(ψγ₂, x)
c_gauss₃ = Gabor.analyse(ψγ₃, x)

c_gaussᵣ₁ = Gabor.analyse(ψγᵣ₁, x)
c_gaussᵣ₂ = Gabor.analyse(ψγᵣ₂, x)
c_gaussᵣ₃ = Gabor.analyse(ψγᵣ₃, x)
c_gaussᵣ₄ = Gabor.analyse(ψγᵣ₄, x)
c_gaussᵣ₅ = Gabor.analyse(ψγᵣ₅, x)
c_gaussᵣ₆ = Gabor.analyse(ψγᵣ₆, x)
c_gaussᵣ₇ = Gabor.analyse(ψγᵣ₇, x)
c_gaussᵣ₈ = Gabor.analyse(ψγᵣ₈, x)
c_gaussᵣ₉ = Gabor.analyse(ψγᵣ₉, x)
c_gaussᵣ₁₀ = Gabor.analyse(ψγᵣ₁₀, x)
c_gaussᵣ₁₁ = Gabor.analyse(ψγᵣ₁₁, x)
