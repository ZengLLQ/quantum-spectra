"""Normal mode calculation and vibrational frequency analysis."""

using LinearAlgebra
include("constants.jl")

struct NormalMode
    mode_index::Int
    wavenumber_cm::Float64 # in cm^-1
    is_imaginary::Bool
    eigenvector_mw::Vector{Float64} # mass-weighted displacements
    displacement_cart::Matrix{Float64} # N x 3 Cartesian displacement vectors
end

struct VibrationalAnalysis
    modes::Vector{NormalMode}
    wavenumbers::Vector{Float64}
    zero_point_energy_kcal_mol::Float64
end

function mass_weighted_hessian(mol::Molecule, hessian_cart_si::Matrix{Float64})::Matrix{Float64}
    N = natoms(mol)
    dim = 3 * N
    H_mw = zeros(dim, dim)

    # Convert masses to kg
    masses_kg = [a.mass_amu * AMU_TO_KG for a in mol.atoms]

    for i in 1:N
        for j in 1:N
            inv_sqrt_mi = 1.0 / sqrt(masses_kg[i])
            inv_sqrt_mj = 1.0 / sqrt(masses_kg[j])
            for a in 1:3
                for b in 1:3
                    idx_i = 3 * (i - 1) + a
                    idx_j = 3 * (j - 1) + b
                    H_mw[idx_i, idx_j] = hessian_cart_si[idx_i, idx_j] * inv_sqrt_mi * inv_sqrt_mj
                end
            end
        end
    end
    return H_mw
end

function calculate_normal_modes(
    mol::Molecule,
    hessian_cart_si::Matrix{Float64};
    project_external::Bool = true
)::VibrationalAnalysis
    N = natoms(mol)
    dim = 3 * N

    H_mw = mass_weighted_hessian(mol, hessian_cart_si)

    if project_external
        P = eckart_projection_matrix(mol)
        H_work = Symmetric(P * H_mw * P)
    else
        H_work = Symmetric(H_mw)
    end

    decomp = eigen(H_work)
    evals = decomp.values
    evecs = decomp.vectors

    modes = Vector{NormalMode}()
    wavenumbers = Vector{Float64}()

    masses_kg = [a.mass_amu * AMU_TO_KG for a in mol.atoms]

    mode_counter = 1
    for k in 1:dim
        λ = evals[k]
        ω = sqrt(abs(λ)) # rad/s
        ν_cm = ω / (2.0 * π * SPEED_OF_LIGHT_CM_S)

        # Ignore translation and rotation zero-frequency modes (< 15 cm^-1)
        if ν_cm > 15.0
            is_imag = λ < 0.0
            signed_nu = is_imag ? -ν_cm : ν_cm

            # Cartesian displacements: L_cart = M^(-1/2) * L_mw
            v_mw = evecs[:, k]
            disp = zeros(N, 3)
            for i in 1:N
                inv_sqrt_m = 1.0 / sqrt(masses_kg[i])
                disp[i, 1] = v_mw[3*(i-1) + 1] * inv_sqrt_m
                disp[i, 2] = v_mw[3*(i-1) + 2] * inv_sqrt_m
                disp[i, 3] = v_mw[3*(i-1) + 3] * inv_sqrt_m
            end
            # Normalize Cartesian displacement
            norm_disp = norm(disp)
            if norm_disp > 0
                disp ./= norm_disp
            end

            push!(modes, NormalMode(mode_counter, signed_nu, is_imag, v_mw, disp))
            push!(wavenumbers, signed_nu)
            mode_counter += 1
        end
    end

    # Calculate Zero-Point Vibrational Energy (ZPE) = 0.5 * sum(h * nu)
    zpe_joules = 0.0
    for nu in wavenumbers
        if nu > 0
            nu_hz = nu * SPEED_OF_LIGHT_CM_S
            zpe_joules += 0.5 * PLANCK_H * nu_hz
        end
    end
    zpe_kcal_mol = zpe_joules * J_TO_KCAL_PER_MOL

    return VibrationalAnalysis(modes, wavenumbers, zpe_kcal_mol)
end
