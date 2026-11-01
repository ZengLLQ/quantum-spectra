"""Vibrational statistical mechanics and thermodynamic functions."""

include("constants.jl")

struct VibrationalThermodynamics
    temperature_k::Float64
    partition_function::Float64
    internal_energy_kcal_mol::Float64
    enthalpy_kcal_mol::Float64
    heat_capacity_cv_cal_mol_k::Float64
    entropy_cal_mol_k::Float64
end

function calculate_vibrational_thermodynamics(
    analysis::VibrationalAnalysis,
    temperature_k::Float64
)::VibrationalThermodynamics
    R = BOLTZMANN_KB * AVOGADRO_NA # J / (mol K) ~ 8.31446
    R_cal = R / 4.184 # cal / (mol K) ~ 1.9872

    qvib = 1.0
    u_vib_j_mol = 0.0
    cv_vib_j_mol_k = 0.0
    s_vib_j_mol_k = 0.0

    for nu in analysis.wavenumbers
        if nu > 10.0 # Only real vibrational modes
            nu_hz = nu * SPEED_OF_LIGHT_CM_S
            theta_v = (PLANCK_H * nu_hz) / BOLTZMANN_KB # Characteristic vibrational temp [K]
            x = theta_v / temperature_k

            # Q_vib
            if x < 100.0
                qvib *= 1.0 / (1.0 - exp(-x))

                exp_x = exp(x)
                # U_vib (including ZPE)
                u_vib_j_mol += AVOGADRO_NA * PLANCK_H * nu_hz * (0.5 + 1.0 / (exp_x - 1.0))

                # Cv_vib
                term_cv = (x^2 * exp_x) / ((exp_x - 1.0)^2)
                cv_vib_j_mol_k += R * term_cv

                # S_vib
                term_s = (x / (exp_x - 1.0)) - log(1.0 - exp(-x))
                s_vib_j_mol_k += R * term_s
            else
                # Low temperature limit: U = ZPE, Cv = 0, S = 0
                u_vib_j_mol += 0.5 * AVOGADRO_NA * PLANCK_H * nu_hz
            end
        end
    end

    u_kcal_mol = (u_vib_j_mol / 1000.0) / 4.184
    h_kcal_mol = u_kcal_mol # For vibrational part, H_vib = U_vib
    cv_cal_mol_k = cv_vib_j_mol_k / 4.184
    s_cal_mol_k = s_vib_j_mol_k / 4.184

    return VibrationalThermodynamics(
        temperature_k,
        qvib,
        u_kcal_mol,
        h_kcal_mol,
        cv_cal_mol_k,
        s_cal_mol_k
    )
end
