module QuantumSpectra

include("constants.jl")
include("molecule.jl")
include("projection.jl")
include("normal_modes.jl")
include("thermodynamics.jl")
include("spectrum.jl")

export Atom, Molecule, natoms, total_mass, center_of_mass, center_molecule, inertia_tensor, is_linear
export eckart_projection_matrix
export mass_weighted_hessian, calculate_normal_modes, NormalMode, VibrationalAnalysis
export calculate_vibrational_thermodynamics, VibrationalThermodynamics
export generate_ir_spectrum, IRSpectrum

end
