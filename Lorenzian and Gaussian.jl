using SymPy

x = Sym("x") 
gamma, sigma = Sym("gamma sigma")

function lorentz(x, gamma)
    # non-normnormalized Lorentzian
    return 1/(gamma^2 + x^2)
end

# Finding the norm of Lorentzian
L_int = sympy.Integral(lorentz(x,gamma), (x, -oo, oo))
L_norm = L_int.doit()
#display(Eq(L_int, L_norm ))

# calculating the half width at half maximum (HWHM) of a Lorenzian
L_max = lorentz(0,gamma)
L_HWHM = Eq(lorentz(x,gamma),(1//2*L_max))
#display(L_HWHM)
solve(L_HWHM, x)
#display(solve(L_HWHM, x))


function gauss(x,sigma)
# non-normalized Gaussian
    return exp(-x^2/2sigma^2)
end

# finding the norm of Gaussian
G_int = sympy.Integral(gauss(x,sigma), (x, -oo, oo))
G_norm = G_int.doit()
#display(Eq(G_int, G_norm))

# calculating the half width at half maximum (HWHM) of a Gaussian
G_max = gauss(0, sigma)
G_HWHM = Eq(gauss(x,sigma),1//2*G_max)
#display(G_HWHM)
solve(G_HWHM, x)
#display(solve(G_HWHM, x))


# example: generating LaTeX code for the displayed equations:
# print(latex(Eq(L_int, L_norm )))