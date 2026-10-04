# The two reflection methods, chosen by `reflections` in [model]: the fixed
# absence rules of the centering, and the full structure factor of the cell.
# Each method gives the reflections of a sample (`reflections`) and the
# scattering weight of each (`scattering_weights`). The form factor of the
# radiation (`form_factor`) and s = sin θ / λ (`scattering_s`) are methods on
# the mode.

"""
    reflections(model::ReflectionModel, sample::Sample, max_hkl_sq::Int)

The reflection families of `sample` with h²+k²+l² ≤ `max_hkl_sq`, as
`(indices, multiplicities)` (see `Miller_indices`).

With `AbsenceRules`, the systematic absences are those of the centering.
"""
function reflections(::AbsenceRules, sample::Sample, max_hkl_sq::Int)
    one_atom(sample)
    return Miller_indices(sample.centering, max_hkl_sq)
end


"""
    scattering_weights(model::ReflectionModel, mode::Radiation, sample::Sample,
                       indices, s) -> Vector{Float64}

The scattering weight of each reflection `indices[i]` at s[i] = sin θ / λ,
normalized to 1 at s = 0, so that the peaks keep their scale relative to the
background.

With `AbsenceRules`, (f(s)/f(0))² · exp(−2B s²) of the single atom, f the form
factor of the radiation (`form_factor`); it does not depend on hkl.
"""
function scattering_weights(::AbsenceRules, mode::Radiation, sample::Sample,
                            indices::AbstractVector{<:AbstractVector{<:Integer}},
                            s::AbstractVector{<:Real})
    atom = one_atom(sample)
    f₀ = form_factor(mode, atom.element, 0.0)
    return (form_factor.(Ref(mode), atom.element, s) ./ f₀) .^ 2 .* Debye_Waller.(s, atom.B)
end


# The single atom of a sample under the absence rules.
function one_atom(sample::Sample)::Atom
    length(sample.atoms) == 1 ||
        throw(ArgumentError("sample $(sample.name): the absence rules need a one-atom basis; use reflections = \"structure_factor\""))
    return only(sample.atoms)
end
