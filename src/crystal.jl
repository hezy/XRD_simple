# Cubic crystallography: Miller indices, multiplicities, d and g spacings.

"""
    cubic_multiplicity(h::Int, k::Int, l::Int)::Int

Compute reflection multiplicity for a cubic crystal from a canonical Miller index.

Assumes the input is in canonical form `h ≥ k ≥ l ≥ 0` and not `[0,0,0]`.
Returns the number of (sign, permutation) variants that share the same |G|² =
h²+k²+l², which equals the multiplicity of the {hkl} family under cubic (m-3m)
point-group symmetry.

# Examples
- {100} → 6, {110} → 12, {111} → 8, {210} → 24, {211} → 24, {321} → 48
"""
function cubic_multiplicity(h::Int, k::Int, l::Int)::Int
    nonzero = (h != 0) + (k != 0) + (l != 0)
    sign_variants = 2^nonzero

    # Given h ≥ k ≥ l ≥ 0, repeats appear only as h==k or k==l
    perms = if h == k == l
        1                 # {hhh}
    elseif h == k || k == l
        3                 # {hhl}, {hh0}, {h00}
    else
        6                 # all distinct
    end

    return perms * sign_variants
end


"""
    Miller_indices(cell_type::String, max_hkl_sq::Int)::Tuple{Vector{Vector{Int}}, Vector{Int}}

Generate canonical Miller indices and reflection multiplicities for cubic crystals.

Enumerates representatives `h ≥ k ≥ l ≥ 0` (excluding `[0,0,0]`) with
`h² + k² + l² ≤ max_hkl_sq`, applies the systematic absence rule for the given
centering, and returns each allowed reflection together with its multiplicity.
Callers sum one peak per representative, weighted by multiplicity — equivalent
to summing over every sign and permutation, at a fraction of the cost.

# Arguments
- `cell_type::String`: "SC", "BCC", or "FCC"
- `max_hkl_sq::Int`: Upper bound on h²+k²+l² (see `bragg_max_hkl_sq`)

# Returns
- `Vector{Vector{Int}}`: Canonical [h,k,l] representatives
- `Vector{Int}`: Multiplicity of each reflection family

# Throws
- `ArgumentError`: If `cell_type` is not "SC", "BCC", or "FCC"
- `ArgumentError`: If `max_hkl_sq < 1`
"""
function Miller_indices(cell_type::String,
                        max_hkl_sq::Int
                        )::Tuple{Vector{Vector{Int}}, Vector{Int}}

    cell_type in ("SC", "BCC", "FCC") || throw(ArgumentError("cell_type must be 'SC', 'BCC', or 'FCC', got '$cell_type'"))
    max_hkl_sq ≥ 1 || throw(ArgumentError("max_hkl_sq must be ≥ 1, got $max_hkl_sq"))

    max_idx = floor(Int, sqrt(max_hkl_sq))
    indices = Vector{Vector{Int}}()
    multiplicities = Vector{Int}()

    for h in 0:max_idx, k in 0:h, l in 0:k
        (h == 0 && k == 0 && l == 0) && continue
        h^2 + k^2 + l^2 > max_hkl_sq && continue

        allowed = if cell_type == "SC"
            true
        elseif cell_type == "BCC"
            iseven(h + k + l)
        else  # FCC
            (iseven(h) && iseven(k) && iseven(l)) || (isodd(h) && isodd(k) && isodd(l))
        end
        allowed || continue

        push!(indices, [h, k, l])
        push!(multiplicities, cubic_multiplicity(h, k, l))
    end

    return indices, multiplicities
end


"""
    d_list(indices::AbstractVector{<:AbstractVector{<:Integer}}, a::Real)::Vector{Float64}

Calculate the interplanar distances (d-spacing) for a cubic crystal structure given Miller indices
and lattice parameter.

# Arguments
- `indices::AbstractVector{<:AbstractVector{<:Integer}}`: Array of Miller indices, where each index is a vector of three 
   integers [h,k,l] representing crystallographic planes
- `a::Real`: Lattice parameter (unit cell edge length) in appropriate units

# Returns
- `Vector{Float64}`: Array of interplanar distances corresponding to each set of Miller indices

# Throws
- `DimensionMismatch`: If any Miller index vector doesn't contain exactly 3 components
- `DomainError`: If lattice parameter is not positive
"""
function d_list(indices::AbstractVector{<:AbstractVector{<:Integer}}, a::Real)::Vector{Float64}
    # Validate lattice parameter
    a > 0 || throw(DomainError(a, "Lattice parameter must be positive"))
    
    # Validate indices structure and dimensions
    for (idx, hkl) in enumerate(indices)
        length(hkl) == 3 || throw(DimensionMismatch(
            "Miller index at position $idx must have exactly 3 components"))
    end
    
    # Pre-allocate output array for better performance
    result = Vector{Float64}(undef, length(indices))
    
    # Calculate d-spacings using direct iteration instead of array comprehension
    # This avoids creating temporary arrays and is more memory efficient
    @inbounds for (i, (h, k, l)) in enumerate(indices)
        result[i] = a / sqrt(h^2 + k^2 + l^2)
    end
    
    return result
end


"""
    g_list(indices::AbstractVector{<:AbstractVector{<:Integer}}, a::Real)::Vector{Float64}

Scattering-vector magnitudes g = |G| = √(h²+k²+l²)/a (1/Å) for cubic Miller
indices. The reciprocal-space analogue of `d_list` (g = 1/d).
"""
function g_list(indices::AbstractVector{<:AbstractVector{<:Integer}}, a::Real)::Vector{Float64}
    a > 0 || throw(DomainError(a, "Lattice parameter must be positive"))

    result = Vector{Float64}(undef, length(indices))
    @inbounds for (i, hkl) in enumerate(indices)
        length(hkl) == 3 || throw(DimensionMismatch(
            "Miller index at position $i must have exactly 3 components"))
        h, k, l = hkl
        result[i] = sqrt(h^2 + k^2 + l^2) / a
    end
    return result
end
