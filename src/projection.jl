"""Eckart translational and rotational projection operator."""

using LinearAlgebra

function eckart_projection_matrix(mol::Molecule)::Matrix{Float64}
    N = natoms(mol)
    c_mol = center_molecule(mol)
    dim = 3 * N

    # Build 6 external vectors
    vectors = Vector{Vector{Float64}}()

    # 3 Translations
    Tx = zeros(dim)
    Ty = zeros(dim)
    Tz = zeros(dim)
    for i in 1:N
        sqrt_m = sqrt(c_mol.atoms[i].mass_amu)
        Tx[3*(i-1) + 1] = sqrt_m
        Ty[3*(i-1) + 2] = sqrt_m
        Tz[3*(i-1) + 3] = sqrt_m
    end
    push!(vectors, Tx)
    push!(vectors, Ty)
    push!(vectors, Tz)

    # 3 Rotations: sqrt(m_i) * (e_alpha x r_i)
    Rx = zeros(dim)
    Ry = zeros(dim)
    Rz = zeros(dim)
    for i in 1:N
        sqrt_m = sqrt(c_mol.atoms[i].mass_amu)
        x, y, z = c_mol.atoms[i].position_angstrom
        # Rx: e_x x r = (0, -z, y)
        Rx[3*(i-1) + 2] = -sqrt_m * z
        Rx[3*(i-1) + 3] =  sqrt_m * y
        # Ry: e_y x r = (z, 0, -x)
        Ry[3*(i-1) + 1] =  sqrt_m * z
        Ry[3*(i-1) + 3] = -sqrt_m * x
        # Rz: e_z x r = (-y, x, 0)
        Rz[3*(i-1) + 1] = -sqrt_m * y
        Rz[3*(i-1) + 2] =  sqrt_m * x
    end
    push!(vectors, Rx)
    push!(vectors, Ry)
    push!(vectors, Rz)

    # Gram-Schmidt orthonormalization
    ortho_vectors = Vector{Vector{Float64}}()
    for v in vectors
        v_proj = copy(v)
        for u in ortho_vectors
            v_proj .-= dot(u, v_proj) .* u
        end
        norm_v = norm(v_proj)
        if norm_v > 1e-6
            push!(ortho_vectors, v_proj ./ norm_v)
        end
    end

    # Build projector P = I - sum(u_k * u_k^T)
    P = Matrix{Float64}(I, dim, dim)
    for u in ortho_vectors
        P .-= u * u'
    end

    return P
end
