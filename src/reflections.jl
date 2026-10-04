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

With `StructureFactor`, every family of the cell is tried, and a family is
absent when F = 0 at every s: when for each kind of atom (element and B) the
sum Σⱼ exp(2πi h·rⱼ) over its sites vanishes, for every member of the family.
An accidental near-cancellation between different elements (the odd reflections
of KCl) is kept, weak.
"""
function reflections(::AbsenceRules, sample::Sample, max_hkl_sq::Int)
    one_atom(sample)
    return Miller_indices(sample.centering, max_hkl_sq)
end

function reflections(::StructureFactor, sample::Sample, max_hkl_sq::Int)
    groups = atom_groups(sample)
    indices, multiplicities = Miller_indices("SC", max_hkl_sq)
    present = [any(abs2(phase_sum(sites, hkl′)) > 1e-10 * length(sites)^2
                   for (_, _, sites) in groups for hkl′ in family_members(hkl))
               for hkl in indices]
    return indices[present], multiplicities[present]
end


"""
    scattering_weights(model::ReflectionModel, mode::Radiation, sample::Sample,
                       indices, s) -> Vector{Float64}

The scattering weight of each reflection `indices[i]` at s[i] = sin θ / λ,
normalized to 1 at s = 0, so that the peaks keep their scale relative to the
background.

With `AbsenceRules`, (f(s)/f(0))² · exp(−2B s²) of the single atom, f the form
factor of the radiation (`form_factor`); it does not depend on hkl.

With `StructureFactor`, |F(hkl)|² / F(000)², averaged over the members of the
family, with

    F(hkl) = Σⱼ fⱼ(s) exp(−Bⱼ s²) exp(2πi (h xⱼ + k yⱼ + l zⱼ))

over every atom of the unit cell (`unit_cell`), and F(000) = Σⱼ fⱼ(0). The
average is |F|² itself when the cell has the full cubic symmetry m-3m. For a
monatomic cell of n atoms, F = n f exp(−B s²) on the allowed reflections and
F(000) = n f(0), so the weight equals that of `AbsenceRules`.
"""
function scattering_weights(::AbsenceRules, mode::Radiation, sample::Sample,
                            indices::AbstractVector{<:AbstractVector{<:Integer}},
                            s::AbstractVector{<:Real})
    atom = one_atom(sample)
    f₀ = form_factor(mode, atom.element, 0.0)
    return (form_factor.(Ref(mode), atom.element, s) ./ f₀) .^ 2 .* Debye_Waller.(s, atom.B)
end

function scattering_weights(::StructureFactor, mode::Radiation, sample::Sample,
                            indices::AbstractVector{<:AbstractVector{<:Integer}},
                            s::AbstractVector{<:Real})
    groups = atom_groups(sample)
    F₀ = sum(length(sites) * form_factor(mode, element, 0.0) for (element, _, sites) in groups)
    weights = Vector{Float64}(undef, length(indices))
    for (i, (hkl, sᵢ)) in enumerate(zip(indices, s))
        # Scattering amplitude of one atom of each kind, with its thermal damping
        f = [form_factor(mode, element, sᵢ) * exp(-B * sᵢ^2) for (element, B, _) in groups]
        members = family_members(hkl)
        F² = sum(abs2(sum(f[g] * phase_sum(groups[g][3], hkl′) for g in eachindex(groups)))
                 for hkl′ in members)
        weights[i] = F² / (length(members) * F₀^2)
    end
    return weights
end


# The atoms of the unit cell, grouped by kind (element, B): a vector of
# (element, B, sites). Atoms of one kind scatter with the same f exp(−B s²).
function atom_groups(sample::Sample)
    groups = Tuple{String,Float64,Vector{NTuple{3,Float64}}}[]
    for atom in unit_cell(sample)
        g = findfirst(g -> g[1] == atom.element && g[2] == atom.B, groups)
        g === nothing ? push!(groups, (atom.element, atom.B, [atom.xyz])) : push!(groups[g][3], atom.xyz)
    end
    return groups
end

# Σⱼ exp(2πi (h xⱼ + k yⱼ + l zⱼ)) over the fractional positions `sites`.
phase_sum(sites, hkl) = sum(cispi(2 * (hkl[1] * x + hkl[2] * y + hkl[3] * z)) for (x, y, z) in sites)


# The single atom of a sample under the absence rules.
function one_atom(sample::Sample)::Atom
    length(sample.atoms) == 1 ||
        throw(ArgumentError("sample $(sample.name): the absence rules need a one-atom basis; use reflections = \"structure_factor\""))
    return only(sample.atoms)
end
