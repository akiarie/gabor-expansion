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


ψgrect, Srect, g̃rect, ψg̃rect, γrect, ψγrect = Gabor.compute_all(rect, 8, 8)

c_rect = Gabor.analyse(ψγrect, x)

c_gauss₁ = Gabor.analyse(ψγ₁, x)
c_gauss₂ = Gabor.analyse(ψγ₂, x)
c_gauss₃ = Gabor.analyse(ψγ₃, x)

c_gaussᵣ₁ = Gabor.analyse(ψγᵣ₁, x)
c_gaussᵣ₂ = Gabor.analyse(ψγᵣ₂, x)
c_gaussᵣ₃ = Gabor.analyse(ψγᵣ₃, x)
c_gaussᵣ₄ = Gabor.analyse(ψγᵣ₄, x)
