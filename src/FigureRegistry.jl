const FIGURE_REGISTRY = Dict(
    "heatmaps" => :fig_maps,
    "pdfs" => :fig_pdf,
    "phase_diagram" => :fig_phase,
    "time_evolution" => :fig_time,
    "phase_magnetic_time" => :fig_phase_B_time,
    "magnetic_fit" => :fig_growth,
    "growth_rate_relations" => :fig_gamma_relations,
    "normalized_magnetic_relations" => :fig_normalized_B_relations,
    "normalized_magnetic_field" => :fig_logB,
    "magnetic_density" => :fig_bn,
    "hro" => :fig_hro,
    "hro_2d" => :fig_hro_2d,
    "hog" => :fig_hog,
    "energy_ratios" => :fig_energy,
    "energy_time" => :fig_energy_time,
    "vorticity" => :fig_vorticity,
    "vorticity_time" => :fig_vorticity_time,
    "enstrophy_density" => :fig_enstrophy_density,
    "power_spectra" => :fig_spectra,
    "density_spectra_time" => :fig_density_spectra_time,
    "velocity_spectra_time" => :fig_velocity_spectra_time,
    "vorticity_spectra_time" => :fig_vorticity_spectra_time,
    "enstrophy_spectra_time" => :fig_enstrophy_spectra_time,
    "magnetic_spectra_time" => :fig_magnetic_spectra_time,
    "structure_functions" => :fig_structure,
    "summary" => :fig_summary,
    "dust_polarization" => :fig_dust,
    "dust_structure" => :fig_dust_structure,
    "dust_pixel_spectrum" => :fig_dust_pixel_spectrum,
    "dust_statistics" => :fig_dust_statistics,
    "dust_p_column" => :fig_dust_p_column,
    "starlight_maps" => :fig_starlight_maps,
    "starlight_structure" => :fig_starlight_structure,
    "starlight_profiles" => :fig_starlight_profiles,
    "starlight_p_column" => :fig_starlight_p_column,
    "zeeman_maps" => :fig_zeeman_maps,
    "zeeman_structure" => :fig_zeeman_structure,
    "zeeman_spectra" => :fig_zeeman_spectra,
    "zeeman_p_column" => :fig_zeeman_p_column,
    "moose" => :fig_moose,
    "moose_structure" => :fig_moose_structure,
    "moose_power_spectra" => :fig_moose_power_spectra,
    "moose_tomography" => :fig_moose_tomography,
    "moose_rmsf" => :fig_moose_rmsf,
    "moose_p_column" => :fig_moose_p_column,
    "polarization_intensity" => :fig_polarization_intensity,
    "shine" => :fig_shine,
    "shine_structure" => :fig_shine_structure,
    "shine_power_spectra" => :fig_shine_power_spectra,
    "shine_rgb" => :fig_shine_rgb,
    "shine_spectrum" => :fig_shine_spectrum,
    "hi_faraday_hog" => :fig_hi_faraday_hog,
    "polarization_time" => :fig_polarization_time,
)

# These diagnostics use scalar or fully three-dimensional fields and therefore
# do not change when the viewing axis changes.  Every other registered figure
# is conservatively treated as LOS-dependent (projected maps, synthetic
# observations, 2D statistics, or mixed summary figures).
const LOS_INDEPENDENT_FIGURES = Set([
    "pdfs",
    "phase_diagram",
    "time_evolution",
    "phase_magnetic_time",
    "magnetic_fit",
    "growth_rate_relations",
    "normalized_magnetic_relations",
    "magnetic_density",
    "hro",
    "energy_ratios",
    "energy_time",
    "vorticity_time",
    "enstrophy_density",
    "power_spectra",
    "density_spectra_time",
    "velocity_spectra_time",
    "vorticity_spectra_time",
    "enstrophy_spectra_time",
    "magnetic_spectra_time",
    "structure_functions",
    "moose_rmsf",
])

function split_figures_by_los(figures)
    requested = unique(String.(figures))
    unknown = setdiff(requested, collect(keys(FIGURE_REGISTRY)))
    isempty(unknown) || error(
        "Unknown figures: $(join(unknown, ", ")).",
    )
    (
        independent = filter(name -> name in LOS_INDEPENDENT_FIGURES, requested),
        dependent = filter(name -> name ∉ LOS_INDEPENDENT_FIGURES, requested),
    )
end

const NOTEBOOK_FIGURES = Dict(
    "dynamo" => [
        "heatmaps",
        "pdfs",
        "time_evolution",
        "phase_magnetic_time",
        "magnetic_fit",
        "growth_rate_relations",
        "normalized_magnetic_relations",
        "normalized_magnetic_field",
        "magnetic_density",
        "hro",
        "hro_2d",
        "hog",
        "energy_ratios",
        "energy_time",
        "vorticity",
        "vorticity_time",
        "enstrophy_density",
        "power_spectra",
        "density_spectra_time",
        "velocity_spectra_time",
        "vorticity_spectra_time",
        "enstrophy_spectra_time",
        "magnetic_spectra_time",
        "structure_functions",
        "summary",
    ],
    "dust" => [
        "dust_polarization",
        "dust_structure",
        "dust_pixel_spectrum",
        "dust_statistics",
        "dust_p_column",
        "polarization_intensity",
    ],
    "starlightpol" => [
        "starlight_maps",
        "starlight_structure",
        "starlight_profiles",
        "starlight_p_column",
    ],
    "zeeman" => [
        "zeeman_maps",
        "zeeman_structure",
        "zeeman_spectra",
        "zeeman_p_column",
    ],
    "moose" => [
        "moose_tomography",
        "moose_rmsf",
        "hi_faraday_hog",
    ],
    "shine" => [
        "phase_diagram",
        "shine",
        "shine_structure",
        "shine_power_spectra",
        "shine_rgb",
        "shine_spectrum",
        "hi_faraday_hog",
        "polarization_time",
    ],
)

available_notebooks() = sort(collect(keys(NOTEBOOK_FIGURES)))

function figures_for_notebooks(notebook_names)
    normalized_names = lowercase.(strip.(String.(notebook_names)))
    unknown = setdiff(normalized_names, available_notebooks())
    isempty(unknown) || error(
        "Unknown notebooks: $(join(unknown, ", ")). Available notebooks: " *
        join(available_notebooks(), ", "),
    )
    unique(vcat((NOTEBOOK_FIGURES[name] for name in normalized_names)...))
end

all(
    figure_name in keys(FIGURE_REGISTRY)
    for figures in values(NOTEBOOK_FIGURES)
    for figure_name in figures
) || error("A notebook group references an unknown figure.")
