"""Molecular geometry and physical properties."""

using LinearAlgebra

struct Atom
    symbol::String
    mass_amu::Float64
    position_angstrom::Vector{Float64} # [x, y, z] in Angstroms
end

struct Molecule
    atoms::Vector{Atom}
end

function natoms(mol::Molecule)::Int
    return length(mol.atoms)
end

function total_mass(mol::Molecule)::Float64
    return sum(a.mass_amu for a in mol.atoms)
end

function center_of_mass(mol::Molecule)::Vector{Float64}
    m_tot = total_mass(mol)
    com = zeros(3)
    for a in mol.atoms
        com .+= a.mass_amu .* a.position_angstrom
    end
    return com ./ m_tot
end

function center_molecule(mol::Molecule)::Molecule
    com = center_of_mass(mol)
    new_atoms = [
        Atom(a.symbol, a.mass_amu, a.position_angstrom .- com)
        for a in mol.atoms
    ]
    return Molecule(new_atoms)
end

function inertia_tensor(mol::Molecule)::Matrix{Float64}
    c_mol = center_molecule(mol)
    I = zeros(3, 3)
    for a in c_mol.atoms
        x, y, z = a.position_angstrom
        m = a.mass_amu
        I[1, 1] += m * (y^2 + z^2)
        I[2, 2] += m * (x^2 + z^2)
        I[3, 3] += m * (x^2 + y^2)
        I[1, 2] -= m * x * y
        I[2, 1] -= m * x * y
        I[1, 3] -= m * x * z
        I[3, 1] -= m * x * z
        I[2, 3] -= m * y * z
        I[3, 2] -= m * y * z
    end
    return I
end

function is_linear(mol::Molecule; tol::Float64=1e-3)::Bool
    if natoms(mol) < 3
        return true
    end
    I = inertia_tensor(mol)
    evals = eigvals(Symmetric(I))
    sort!(evals)
    return evals[1] / evals[3] < tol
end
