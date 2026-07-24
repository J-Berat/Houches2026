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
    "hog" => :fig_hog,
    "energy_ratios" => :fig_energy,
    "energy_time" => :fig_energy_time,
    "vorticity" => :fig_vorticity,
    "enstrophy_density" => :fig_enstrophy_density,
    "power_spectra" => :fig_spectra,
    "structure_functions" => :fig_structure,
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
    "moose_tomography" => :fig_moose_tomography,
    "moose_p_column" => :fig_moose_p_column,
    "shine" => :fig_shine,
    "shine_structure" => :fig_shine_structure,
    "shine_rgb" => :fig_shine_rgb,
    "shine_spectrum" => :fig_shine_spectrum,
)

const NOTEBOOK_FIGURES = Dict(
    "dynamo" => [
        "heatmaps",
        "pdfs",
        "phase_diagram",
        "time_evolution",
        "phase_magnetic_time",
        "magnetic_fit",
        "growth_rate_relations",
        "normalized_magnetic_relations",
        "normalized_magnetic_field",
        "magnetic_density",
        "hro",
        "hog",
        "energy_ratios",
        "energy_time",
        "vorticity",
        "enstrophy_density",
        "power_spectra",
        "structure_functions",
    ],
    "dust" => [
        "dust_polarization",
        "dust_structure",
        "dust_pixel_spectrum",
        "dust_statistics",
        "dust_p_column",
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
        "moose",
        "moose_structure",
        "moose_tomography",
        "moose_p_column",
    ],
    "shine" => [
        "shine",
        "shine_structure",
        "shine_rgb",
        "shine_spectrum",
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
