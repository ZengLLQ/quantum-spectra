"""Physical and chemical constants in SI and atomic units."""

# Planck constant [J s]
const PLANCK_H = 6.62607015e-34

# Reduced Planck constant [J s]
const HBAR = PLANCK_H / (2.0 * π)

# Speed of light in vacuum [m s^-1]
const SPEED_OF_LIGHT_M_S = 299792458.0

# Speed of light in vacuum [cm s^-1]
const SPEED_OF_LIGHT_CM_S = SPEED_OF_LIGHT_M_S * 100.0

# Avogadro constant [mol^-1]
const AVOGADRO_NA = 6.02214076e23

# Boltzmann constant [J K^-1]
const BOLTZMANN_KB = 1.380649e-23

# Unified atomic mass unit [kg]
const AMU_TO_KG = 1.66053906660e-27

# Angstrom to meter
const ANGSTROM_TO_M = 1.0e-10

# Joule to kcal/mol
const J_TO_KCAL_PER_MOL = 1.0e-3 * AVOGADRO_NA / 4184.0
