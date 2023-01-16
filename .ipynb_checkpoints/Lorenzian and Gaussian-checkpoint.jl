using SymPy

x = Sym("x") 
γ , σ = Sym("γ σ")

function lorentz(x, γ)
    # non-normnormalized Lorentzian
    return 1/(γ^2 + x^2)
end

print("Finding the norm of Lorentzian:")
L_int = sympy.Integral(lorentz(x,γ), (x, -oo, oo))
L_norm = L_int.doit()
latexify(Eq(L_int, L_norm ))

print("calculating the half width at half maximum (HWHM) of a Lorenzian:")
L_max = lorentz(0,γ)
L_HWHM = Eq(lorentz(x,γ),(1//2*L_max))
print(L_HWHM)
solve(L_HWHM, x)
print(solve(L_HWHM, x))


function gauss(x, σ)
# non-normalized Gaussian
    return exp(-x^2/2σ^2)
end

print("Finding the norm of Gaussian:")
G_int = sympy.Integral(gauss(x,σ), (x, -oo, oo))
G_norm = G_int.doit()
print(Eq(G_int, G_norm))

# calculating the half width at half maximum (HWHM) of a Gaussian
G_max = gauss(0, σ)
G_HWHM = Eq(gauss(x,σ),1//2*G_max)
print(G_HWHM)
solve(G_HWHM, x)
print(solve(G_HWHM, x))


# example: generating LaTeX code for the displayed equations:
# print(latex(Eq(L_int, L_norm )))