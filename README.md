# Houches2026 MHD Analysis Notebooks

Interactive [Pluto.jl](https://plutojl.org/) notebooks for exploring three-dimensional magnetohydrodynamic (MHD) simulations and producing synthetic observables. The project was prepared for the Les Houches 2026 summer school and is designed to run both on a personal computer and on a remote computing server through an SSH tunnel.

The notebooks read HDF5 or FITS simulation snapshots, expose the main physical and observational parameters as interactive controls, and generate publication-ready figures with CairoMakie.

## What is included?

| File | Purpose |
|---|---|
| `dynamo.jl` | General MHD diagnostics: maps, PDFs, phase diagrams, magnetic growth, energy ratios, vorticity, spectra, structure functions, comparative $B$--$n$ relations, HRO, and HOG. |
| `dust.jl` | Synthetic thermal-dust Stokes emission, polarization diagnostics, and structure functions for every projected observable. |
| `moose.jl` | MOOSE synchrotron emission, Faraday rotation, interferometric filtering, RM synthesis, and projected-observable structure functions. |
| `shine.jl` | SHINE synthetic H I 21-cm emission, phase-separated columns, velocity moments, spectra, RGB composites, and structure functions for all scalar maps. |
| `zeeman.jl` | Synthetic H I Stokes-I/V spectra, Zeeman magnetic-field recovery, and structure functions for all Zeeman maps. |
| `starlightpol.jl` | Dichroic starlight polarization using cell-by-cell Mueller propagation, with structure functions for all projected Stokes and derived maps. |

The five synthetic-observation notebooks compute axis-averaged, periodic, two-dimensional structure functions for every scalar sky map. Their common controls are in **Shared observational beam**: enable or disable the figures, select the order $p$, and choose the number of sampled separations. Polarization-angle increments use the shortest difference modulo $180^\circ$.

The Dynamo notebook also compares magnetic-field strength with number density, including a binned $B\propto n^\kappa$ fit. Its 3-D HRO measures the orientation of the magnetic field relative to isodensity structures as a function of density. Its projected HOG compares the gradients of column density and density-weighted magnetic-field strength and reports the normalized projected Rayleigh statistic.

The isotropic power-spectrum panels fit $E(k)=A k^\alpha$ over user-selected $k_{\min}$ and $k_{\max}$. Each panel reports $\alpha$ and $R^2$ and can overlay a $k^{-5/3}$ Kolmogorov reference normalized at the geometric center of the fitted interval.

All figures use LaTeX rendering consistently for x- and y-axis labels, numeric tick labels, and colorbar labels and ticks. Scientific-notation ticks are rendered as $a\times10^b$ rather than as plain-text Unicode.
| `run_pluto.jl` | Interactive launcher that asks which notebook to open. |
| `export_html.jl` | Executes one notebook from top to bottom and exports a self-contained HTML snapshot. |
| `dynamo_diagnostics.jl` | Master notebook containing every analysis section. This is the source of truth for shared code. |
| `split_notebooks.jl` | Regenerates the six focused notebooks from the master notebook and its dependency graph. |
| `Project.toml`, `Manifest-v1.11.toml`, `Manifest-v1.12.toml` | Reproducible environments with dependency versions resolved separately for Julia 1.11 and 1.12. |

## Requirements

- Julia 1.11 or Julia 1.12 (tested with Julia 1.11.7 and 1.12.6).
- A modern web browser.
- Access to HDF5 or FITS simulation snapshots.
- SSH access if Pluto is run on a remote server.

The Julia packages are declared in `Project.toml`; no manual package-by-package installation is required. Direct dependencies include Pluto, PlutoUI, CairoMakie, FFTW, FITSIO, HDF5, LaTeXStrings, StatsBase, and the required Julia standard libraries. `Manifest-v1.11.toml` and `Manifest-v1.12.toml` lock the complete dependency graph separately for the two supported Julia versions. Julia automatically selects the matching file.

Each generated notebook also embeds the same `Project.toml` and `Manifest.toml`, so opening a focused notebook directly in Pluto uses the reproducible environment rather than an unrelated global Julia environment.

## Quick start on a local computer

Clone the repository and enter it:

```bash
git clone https://github.com/J-Berat/Houches2026.git
cd Houches2026
```

Install the project dependencies once after cloning (and again when the manifest for your Julia version changes):

```bash
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

Start the interactive launcher:

```bash
julia --project=. run_pluto.jl
```

At startup, the launcher prints the Julia version and the selected manifest. For Julia 1.11 it must report `Manifest-v1.11.toml`; for Julia 1.12 it must report `Manifest-v1.12.toml`. It then displays a numbered menu and opens the selected notebook without repeating the package-installation step. For example:

```text
1. Dynamo
2. Dust
3. MOOSE
4. SHINE
5. ZEEMAN
6. StarlightPol
```

Enter a number from 1 to 6. On macOS or Windows, the browser opens automatically.

You may also bypass the menu:

```bash
julia --project=. run_pluto.jl dust
```

The argument may be `dynamo`, `dust`, `moose`, `shine`, `zeeman`, or `starlightpol`, with or without the `.jl` extension.

## Running on the PSMN or another remote server

Pluto runs on the server, while the web interface is displayed in your local browser. Keep the SSH connection open for the entire session.

On a shared server, **do not rely on Pluto's default port**. Choose an unused personal port for every session to avoid collisions with other users. The example below uses `15432`; replace it with another port if necessary. The same port number must be used in both the SSH tunnel and the `PLUTO_PORT` variable.

### 1. Open an SSH tunnel from your local computer

```bash
ssh -L 15432:127.0.0.1:15432 PSMN_sr650node230
```

### 2. Start Pluto on the server

The first time only, install the environment:

```bash
cd /Xnfs/Houches2026/DynSim/notebooks
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

After that one-time setup, start the launcher directly:

```bash
cd /Xnfs/Houches2026/DynSim/notebooks
PLUTO_PORT=15432 julia --project=. run_pluto.jl
```

Choose a notebook from the menu. Pluto prints a URL containing a session secret, similar to:

```text
http://localhost:15432/?secret=...
```

Open that exact URL in the browser on your local computer. The SSH tunnel securely forwards it to the server.

Press `Ctrl+C` in the server terminal when you want to stop Pluto.

## Selecting the data repository

The notebooks search for data in this order:

1. The path stored in the `DYNAMO_DATA_REPOSITORY` environment variable.
2. `/Xnfs/Houches2026/DynSim/cooling_freq_output`, when it exists.
3. The bundled `cooling/VaryingMach` directory, when present beside the notebooks.

To select a repository before starting Pluto:

```bash
export DYNAMO_DATA_REPOSITORY=/path/to/your/simulations
julia --project=. run_pluto.jl
```

You can also change the path interactively in the notebook and click **Load path**.

The loader only considers snapshot files located inside a directory named
`DataCubes`. The selected path may point to any of the following:

- a repository containing several simulation families;
- one family such as `VaryingMach`;
- one simulation directory;
- a `DataCubes` directory.

Directory names and nesting are discovered recursively, but directories that
contain snapshots directly without a `DataCubes` level are ignored. A typical
layout is:

```text
simulation-root/
└── VaryingMach/
    ├── run_turb_cooling_mhd_lo_mach/
    │   └── DataCubes/
    │       ├── info_00001.h5
    │       └── info_00002.h5
    ├── run_turb_cooling_mhd_mi_mach/
    │   └── DataCubes/
    └── run_turb_cooling_mhd_hi_mach/
        └── DataCubes/
```

## Supported snapshot formats

Recognized extensions are `.h5`, `.hdf5`, `.fits`, `.fit`, and `.fts`.

### HDF5

The loader searches recursively inside each HDF5 file. Common dataset aliases are recognized automatically:

| Physical quantity | Examples of accepted aliases |
|---|---|
| Density | `rho`, `density`, `massdensity` |
| Pressure | `p`, `pressure`, `press`, `thermalpressure`, `gaspressure` |
| Velocity | `vx`, `vy`, `vz`, `velx`, `vely`, `velz`, `velocityx`, `velocityy`, `velocityz` |
| Magnetic field | `bx`, `by`, `bz`, `bfieldx`, `bfieldy`, `bfieldz` |
| Box length | `l`, `boxlength`, `boxsize`, `domainlength` |
| Time | `t`, `time`, `snapshot_time`, `simtime` |

Centered magnetic components are preferred. Face-centered pairs such as `bx_l`/`bx_r` are also supported and are averaged to cell centers.

If zero or several datasets match a physical field, use the **HDF5 field mapping** table in the notebook to select the correct dataset explicitly. A manual selection overrides automatic alias matching.

### FITS

FITS snapshots may be either:

- one multi-extension FITS file with named extensions; or
- one directory per snapshot containing a separate FITS image for each field.

Useful extension names include `RHO`, `PRESSURE`, `VX`, `VY`, `VZ`, `BX`, `BY`, and `BZ`. Face-centered magnetic pairs such as `BX_L` and `BX_R` are also accepted. Box size and time may be supplied through `L`, `LBOX`, `BOXSIZE`, `TIME`, `T`, or `SIMTIME` metadata.

## Using a notebook

The notebooks start in lazy mode. Only the selected 3-D snapshot is loaded;
expensive time-series calculations do not run immediately.

1. Select a data repository and click **Load path**.
2. Select the active simulation, snapshot, and line of sight.
3. Check the physical-unit conversion factors.
4. Configure the diagnostic or synthetic-observation parameters.
5. Click a result cell and press `Shift+Enter`.

Pluto automatically evaluates every upstream dependency required by that result. This is usually faster and uses less memory than running the entire notebook.

Temporal figures are disabled initially. Enabling a temporal figure explicitly
reads the required snapshots sequentially, closes each HDF5 file immediately,
and retains only scalar summaries. The raw-cube cache contains at most one
snapshot; moving the snapshot slider evicts the previous raw cube before reading
the newly selected one.

Each simulation is capped at 40 snapshots. When a `DataCubes` directory contains
more than 40 files, the notebooks select 40 evenly spaced snapshots including
the first and last, preserving the full simulated time interval without opening
the omitted cubes.

For comparative diagnostics, select any number of simulations under **Simulations in comparative plots**. One simulation is selected initially to keep startup light. Comparative PDFs, distributions, and histograms use the chosen snapshot index for every selected run; if a run is shorter, its last available snapshot is used. Spatial maps continue to use the active **Run** and **Snapshot** selections.

## Exporting a complete notebook to HTML

The HTML exporter evaluates every cell and stops if any cell fails. Select a notebook with `DYNAMO_NOTEBOOK`:

```bash
DYNAMO_NOTEBOOK=dust.jl julia --project=. export_html.jl
```

This creates `dust.html` beside the notebook. To choose another output path:

```bash
DYNAMO_NOTEBOOK=dust.jl \
DYNAMO_HTML_PATH=/path/to/results/dust.html \
julia --project=. export_html.jl
```

The exported HTML is a read-only snapshot. Interactive controls require a live Pluto session.

## Regenerating the focused notebooks

Shared logic should be edited in `dynamo_diagnostics.jl`. Then regenerate all focused notebooks with:

```bash
julia --project=. split_notebooks.jl
```

The generator uses Pluto's dependency graph to include only the cells needed by each focused notebook.

> **Important:** do not make long-lived edits only in `dynamo.jl`, `dust.jl`, `moose.jl`, `shine.jl`, `zeeman.jl`, or `starlightpol.jl`. Regeneration can overwrite those edits. Apply shared or scientific changes to `dynamo_diagnostics.jl` first.

## Useful environment variables

| Variable | Meaning | Default |
|---|---|---|
| `DYNAMO_DATA_REPOSITORY` | Simulation-data root | PSMN path, then bundled data |
| `DYNAMO_NOTEBOOK` | Notebook used by the launcher or HTML exporter | Interactive menu for launcher; `dynamo.jl` for exporter |
| `DYNAMO_HTML_PATH` | HTML export destination | `<notebook-name>.html` |
| `PLUTO_HOST` | Pluto listening address | `127.0.0.1` |
| `PLUTO_PORT` | Pluto web-interface port | `1234` locally; set an unused personal port explicitly on a shared server |
| `PLUTO_LAUNCH_BROWSER` | Whether the launcher opens a browser | `true` on macOS/Windows, `false` on Linux |

Boolean values such as `1`, `true`, `yes`, and `on` are accepted for `PLUTO_LAUNCH_BROWSER`.

## Troubleshooting

### `julia` is not recognized

Julia is not available on your `PATH`. Load the Julia module provided by the server, or install Julia 1.11/1.12 and start a new terminal. Confirm the installation with `julia --version`.

### Packages are missing

Run the installation command once from the repository directory:

```bash
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

### The browser cannot reach Pluto on the server

Confirm that:

- the SSH session containing the port-forwarding option is still open;
- the local forwarding port, remote forwarding port, and `PLUTO_PORT` value are identical;
- you opened the complete URL printed by Pluto, including its `?secret=...` token.

If the selected port is already used, stop Pluto and choose a different number on both sides. For example:

```bash
ssh -L 16543:127.0.0.1:16543 PSMN_sr650node230
```

Then, on the server:

```bash
PLUTO_PORT=16543 julia --project=. run_pluto.jl
```

### No snapshots are found

Check that the repository contains a non-empty directory named `DataCubes`, that its snapshots use a supported extension, and that you have permission to read it. The path may point directly to a `DataCubes` directory. Set `DYNAMO_DATA_REPOSITORY` explicitly if the default PSMN location is not appropriate.

### An HDF5 field is missing or ambiguous

Read the accepted aliases shown in the mapping table, inspect the available-dataset dropdown, and choose the correct dataset manually. The chosen dataset must exist in every snapshot of the selected run.

### Calculations use too much memory

HDF5 access is serialized: only one file is open at a time, every dataset handle
is closed immediately after reading, the raw cache contains at most one cube,
and comparative diagnostics retain derived profiles or maps instead of a full
3-D cube for every simulation. Keep lazy
execution enabled, run only the needed result cells, and reduce the number of
selected simulations, maps, or spectral channels when memory is limited.

## Reproducibility notes

- Physical units and conversion factors are visible near the top of each notebook.
- Density-weighted projections, periodic-boundary operations, and observational conventions are documented beside their controls.
- `Manifest-v1.11.toml` and `Manifest-v1.12.toml` record dependency graphs resolved by the corresponding Julia versions. The generic `Manifest.toml` is kept only as a fallback for older tooling.
- The instrument-noise options use an explicit random seed.
- HTML export fails rather than silently publishing a notebook with errored cells.

## Repository structure

```text
Houches2026/
├── README.md
├── Project.toml
├── Manifest.toml
├── Manifest-v1.11.toml
├── Manifest-v1.12.toml
├── run_pluto.jl
├── export_html.jl
├── split_notebooks.jl
├── dynamo_diagnostics.jl
├── dynamo.jl
├── dust.jl
├── moose.jl
├── shine.jl
├── zeeman.jl
└── starlightpol.jl
```

## Getting help

When reporting a problem, include:

- the notebook name;
- the Julia version from `julia --version`;
- the data format and repository layout;
- the complete Pluto error message;
- the selected run, snapshot, and line of sight;
- any manual HDF5 field mappings or non-default unit conversions.
