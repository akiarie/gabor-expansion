# Gabor Expansion in Julia

This repository contains a basic Julia implementation of the discrete Gabor expansion.

In it are defined structs and methods that provide for signal analysis and synthesis, as well as the
computation of biorthogonal functions, both using the inverse frame operator and with the Wexler-Raz
matrix.

All the main functionality is in [gabor.jl](gabor.jl), while examples of what can be generated with
the module are in [windows](windows/) and [signals](signals/).
