using Test
using LinearAlgebra
using QuantumSpectra

@testset "QuantumSpectra Test Suite" begin
    @testset "Molecular Geometry & Inertia" begin
        # Water molecule (H2O)
        # O at (0, 0, 0), H1 at (0, 0.757, 0.586), H2 at (0, -0.757, 0.586)
        h2o = Molecule([
            Atom("O", 15.9994, [0.0, 0.0, 0.0]),
            Atom("H", 1.00784, [0.0, 0.757, 0.586]),
            Atom("H", 1.00784, [0.0, -0.757, 0.586]),
        ])

        @test natoms(h2o) == 3
        @test total_mass(h2o) ≈ 18.01508 atol=1e-4

        com = center_of_mass(h2o)
        @test com[1] ≈ 0.0 atol=1e-6
        @test com[2] ≈ 0.0 atol=1e-6
        @test com[3] > 0.0 # Shifted towards oxygen

        centered = center_molecule(h2o)
        @test norm(center_of_mass(centered)) ≈ 0.0 atol=1e-6

        @test !is_linear(h2o)

        # Carbon monoxide (CO) - diatomic (linear)
        co = Molecule([
            Atom("C", 12.011, [0.0, 0.0, 0.0]),
            Atom("O", 15.999, [0.0, 0.0, 1.128]),
        ])
        @test is_linear(co)
    end

    @testset "Eckart Projection Operator" begin
        h2o = Molecule([
            Atom("O", 15.9994, [0.0, 0.0, 0.0]),
            Atom("H", 1.00784, [0.0, 0.757, 0.586]),
            Atom("H", 1.00784, [0.0, -0.757, 0.586]),
        ])

        P = eckart_projection_matrix(h2o)
        dim = 3 * natoms(h2o)

        # Check idempotent P^2 = P
        @test P * P ≈ P atol=1e-8
        # Check symmetric P^T = P
        @test P' ≈ P atol=1e-8
        # Trace of P should equal 3N - 6 = 9 - 6 = 3 (number of vibrational DOFs)
        @test tr(P) ≈ 3.0 atol=1e-4
    end

    @testset "Analytical Diatomic Harmonic Oscillator" begin
        # Diatomic molecule: two masses m1, m2 separated along z-axis
        m1_amu = 12.0
        m2_amu = 16.0
        r0 = 1.128 # Angstroms

        diatomic = Molecule([
            Atom("C", m1_amu, [0.0, 0.0, 0.0]),
            Atom("O", m2_amu, [0.0, 0.0, r0]),
        ])

        amu_kg = 1.66053906660e-27
        c_cm_s = 2.99792458e10
        mu_kg = (m1_amu * m2_amu / (m1_amu + m2_amu)) * amu_kg
        k_si = 1800.0 # N/m (typical strong bond)
        nu_exact_cm = (1.0 / (2.0 * π * c_cm_s)) * sqrt(k_si / mu_kg)

        H_cart = zeros(6, 6)
        H_cart[3, 3] = k_si
        H_cart[6, 6] = k_si
        H_cart[3, 6] = -k_si
        H_cart[6, 3] = -k_si

        result = calculate_normal_modes(diatomic, H_cart; project_external=true)

        @test length(result.modes) == 1
        @test result.modes[1].wavenumber_cm ≈ nu_exact_cm rtol=1e-3
        @test result.zero_point_energy_kcal_mol > 0.0
    end

    @testset "Water Molecule 3 Normal Modes" begin
        # O at (0, 0, 0), H1 at (0, 0.757, 0.586), H2 at (0, -0.757, 0.586)
        r_O = [0.0, 0.0, 0.0]
        r_H1 = [0.0, 0.757, 0.586]
        r_H2 = [0.0, -0.757, 0.586]

        h2o = Molecule([
            Atom("O", 15.9994, r_O),
            Atom("H", 1.00784, r_H1),
            Atom("H", 1.00784, r_H2),
        ])

        # Construct Urey-Bradley harmonic valence force field
        # Bond 1: O-H1 (k1 = 845 N/m)
        # Bond 2: O-H2 (k2 = 845 N/m)
        # Bond 3: H1-H2 (k3 = 120 N/m)
        H_cart = zeros(9, 9)

        function add_bond_hessian!(H, idx1, idx2, r1, r2, k)
            dr = r2 .- r1
            u = dr ./ norm(dr)
            block = k .* (u * u')
            i1 = 3 * (idx1 - 1) + 1
            i2 = 3 * (idx2 - 1) + 1
            H[i1:i1+2, i1:i1+2] .+= block
            H[i2:i2+2, i2:i2+2] .+= block
            H[i1:i1+2, i2:i2+2] .-= block
            H[i2:i2+2, i1:i1+2] .-= block
        end

        add_bond_hessian!(H_cart, 1, 2, r_O, r_H1, 845.0)
        add_bond_hessian!(H_cart, 1, 3, r_O, r_H2, 845.0)
        add_bond_hessian!(H_cart, 2, 3, r_H1, r_H2, 120.0)

        result = calculate_normal_modes(h2o, H_cart; project_external=true)

        # Must have exactly 3 internal vibrational normal modes (3N - 6)
        @test length(result.modes) == 3
        for m in result.modes
            @test m.wavenumber_cm > 100.0
            @test !m.is_imaginary
        end
        @test result.zero_point_energy_kcal_mol > 0.0

        # Test thermodynamics
        thermo_298 = calculate_vibrational_thermodynamics(result, 298.15)
        @test thermo_298.partition_function >= 1.0
        @test thermo_298.heat_capacity_cv_cal_mol_k >= 0.0
        @test thermo_298.entropy_cal_mol_k >= 0.0

        thermo_2000 = calculate_vibrational_thermodynamics(result, 2000.0)
        @test thermo_2000.heat_capacity_cv_cal_mol_k > thermo_298.heat_capacity_cv_cal_mol_k
    end

    @testset "Synthetic IR Spectrum" begin
        wavenumbers = [1595.0, 3657.0, 3756.0] # H2O experimental frequencies
        intensities = [20.0, 10.0, 50.0]

        spec_lor = generate_ir_spectrum(wavenumbers, intensities; lineshape=:lorentzian)
        @test length(spec_lor.absorbance) == 1000
        @test maximum(spec_lor.absorbance) > 0.0
        @test length(spec_lor.peaks) == 3

        spec_gauss = generate_ir_spectrum(wavenumbers, intensities; lineshape=:gaussian)
        @test maximum(spec_gauss.absorbance) > 0.0
    end
end
