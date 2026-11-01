"""Synthetic Infrared (IR) absorption spectrum generator."""

struct IRSpectrum
    wavenumbers_grid::Vector{Float64}
    absorbance::Vector{Float64}
    peaks::Vector{Tuple{Float64, Float64}} # (wavenumber_cm, intensity)
end

function generate_ir_spectrum(
    wavenumbers::Vector{Float64},
    intensities::Vector{Float64};
    w_min::Float64 = 400.0,
    w_max::Float64 = 4000.0,
    n_points::Int = 1000,
    fwhm::Float64 = 20.0,
    lineshape::Symbol = :lorentzian
)::IRSpectrum
    grid = range(w_min, w_max, length=n_points)
    absorbance = zeros(n_points)

    gamma = fwhm / 2.0
    sigma = fwhm / (2.0 * sqrt(2.0 * log(2.0)))

    for (nu0, int0) in zip(wavenumbers, intensities)
        if nu0 >= w_min - 50.0 && nu0 <= w_max + 50.0
            for i in 1:n_points
                w = grid[i]
                if lineshape == :lorentzian
                    profile = (1.0 / π) * (gamma / ((w - nu0)^2 + gamma^2))
                else # :gaussian
                    profile = (1.0 / (sigma * sqrt(2.0 * π))) * exp(-((w - nu0)^2) / (2.0 * sigma^2))
                end
                absorbance[i] += int0 * profile
            end
        end
    end

    peaks = [
        (nu, int) for (nu, int) in zip(wavenumbers, intensities)
        if nu >= w_min && nu <= w_max
    ]

    return IRSpectrum(collect(grid), absorbance, peaks)
end
