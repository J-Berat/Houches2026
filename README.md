# Houches2026 — analyse MHD

Ce dépôt contient une seule source scientifique pour l'analyse de simulations
MHD et la production d'observables synthétiques :
`dynamo_diagnostics.jl`.

Le même code scientifique est utilisable de deux façons :

- avec Pluto, pour choisir manuellement un dossier et explorer les résultats ;
- avec `run_figures.jl`, pour produire des figures comparatives sans interface
  et sans initialiser Pluto.

## Structure du dépôt

```text
LesHouchesGit/
├── dynamo_diagnostics.jl       # source scientifique unique
├── run_pluto.jl                # lancement interactif
├── run_figures.jl              # calcul batch sans interface
├── export_html.jl              # export du notebook maître
├── Project.toml
├── Manifest.toml               # environnement de repli (toute version)
├── Manifest-v1.11.toml         # environnement résolu pour Julia 1.11
├── Manifest-v1.12.toml         # environnement résolu pour Julia 1.12
├── notebooks/                  # notebooks spécialisés générés
│   ├── dynamo.jl
│   ├── dust.jl
│   ├── moose.jl
│   ├── shine.jl
│   ├── starlightpol.jl
│   └── zeeman.jl
├── src/
│   ├── DynamoAnalysis.jl       # moteur batch Julia natif
│   ├── FigureRegistry.jl       # noms des figures disponibles
│   └── BatchCellIndex.jl       # dépendances générées, sans code scientifique
└── tools/
    ├── generate_batch_index.jl
    └── split_notebooks.jl
```

Les six notebooks spécialisés sont conservés pour l'utilisation interactive,
mais regroupés dans `notebooks/`. Ils sont générés depuis
`dynamo_diagnostics.jl` : il ne faut donc pas y faire de modifications
scientifiques manuelles. Cela évite plusieurs sources divergentes tout en
gardant les notebooks Dynamo, Dust, MOOSE, SHINE, StarlightPol et ZEEMAN.

## Installation

Julia 1.11 ou 1.12 est recommandé.

```bash
cd "/Users/jb270005/Desktop/LesHouchesGit"
julia --project=. -e 'import Pkg; Pkg.instantiate()'
```

Le dépôt versionne trois manifestes, ce qui est voulu : Julia sélectionne
automatiquement le manifeste correspondant à sa version
(`Manifest-v1.11.toml` sous Julia 1.11, `Manifest-v1.12.toml` sous Julia
1.12) et retombe sur `Manifest.toml` pour toute autre version. Chaque
manifeste est régénéré par `Pkg.instantiate()` / `Pkg.resolve()` sous la
version de Julia concernée ; il ne faut pas les éditer à la main.

## Utilisation interactive avec Pluto

```bash
cd "/Users/jb270005/Desktop/LesHouchesGit"
julia --project=. run_pluto.jl
```

Le lanceur demande s'il faut ouvrir le notebook complet ou l'un des six
notebooks spécialisés. Pour ouvrir directement Dust, par exemple :

```bash
julia --project=. run_pluto.jl dust
```

Chaque notebook démarre sans chemin de données et n'ouvre aucun cube. Dans
Pluto, renseigner **Data path**, puis cliquer sur **Load path**. Le chemin peut
pointer vers n'importe quel dossier local, disque externe ou système de
fichiers monté.

Il est aussi possible de fournir le chemin au lancement :

```bash
DYNAMO_DATA_REPOSITORY="/chemin/vers/les/simulations" \
julia --project=. run_pluto.jl
```

Le chemin peut désigner :

- un dossier contenant plusieurs familles de simulations ;
- une famille de simulations ;
- une simulation unique ;
- un dossier contenant directement des snapshots HDF5 ou FITS.

Les extensions reconnues sont `.h5`, `.hdf5`, `.fits`, `.fit` et `.fts`.

## Calcul des figures sans Pluto

`run_figures.jl` launches the three comparison folders below:

```text
/Xnfs/Houches2026/DynSim/cooling_freq_output/
├── VaryingMach/
├── VaryingRes/
└── VaryingRatio/
```

The run lists and folders are configured near the top of `run_figures.jl`.
`SELECTED_COMPARISONS` can contain `mach`, `resolution`, and/or `ratio`.
For each selected comparison and snapshot window, the script separates the
requested figures into:

- one `common` job for LOS-independent 3D and temporal diagnostics;
- one job in each of `los_x`, `los_y`, and `los_z` for projected maps,
  synthetic observations, 2D statistics, and mixed summary figures.

Dynamo still uses the first 20 snapshots. Dust, StarlightPol, ZEEMAN, MOOSE,
and SHINE still use the last 10 snapshots. The `common` products are not
recomputed for each viewing direction.

In batch mode, the comparative density PDFs, magnetic-field--density relation,
3D HRO, 2D HRO, and power spectra use every selected snapshot. Their solid
curves are snapshot medians and their shaded regions span the 16th--84th
percentiles. The 2D HRO reports the three lines of sight, bootstrap uncertainty
on the orientation parameter, and the number of contributing pixels. Power
spectra include compensated panels, the fitted interval, slope uncertainty,
the forcing wavenumber, and the Nyquist limit. A simulation keeps the same
color throughout a comparison group. The summary also includes the magnetic
field evolution separated into CNM, LNM, and WNM thermal phases.

The five time-resolved spectrum figures overlay every selected density,
velocity, vorticity, enstrophy, or magnetic snapshot in one panel per
simulation. Their horizontal coordinate is the dimensionless mode
\(kL/(2\pi)\). A shared physical-time colorbar, together with increasing
opacity and line width, shows temporal ordering while common axes preserve
direct comparison between simulations.

MOOSE and SHINE each export an additional four-panel spatial power-spectrum
figure. MOOSE covers Faraday depth, synchrotron intensity, polarized intensity,
and peak Faraday-spectrum amplitude. SHINE covers H I column density and the
zeroth, first, and second brightness-temperature moments. The spectra show
azimuthal-mode uncertainty, the fitted interval and slope, a forcing-scale
guide, and the projected Nyquist limit.

MOOSE uses a fixed physical line-of-sight depth of 100 pc for Faraday rotation
and synchrotron-emissivity integration. Its default LoTSS band is
120--168 MHz with 97.6 kHz channels. The MOOSE batch also exports the
rotation-measure spread function, including its measured and theoretical FWHM,
complex response, and strongest sidelobe level. The H I--Faraday HOG follows
the noise experiment of Berat et al. (2026): independent white Gaussian noise
is added to every native-resolution H I velocity channel with the default
\(\sigma_{\rm EBHIS}=90\,\mathrm{mK}\), without an EBHIS beam convolution.
Each simulation has a reproducible, independent realization; the Faraday cube
remains noise-free.

Check all 24 planned jobs (three comparisons × two snapshot windows ×
`common` plus three LOS folders) without computing:

```bash
DYNAMO_DRY_RUN=true \
julia --threads=auto --startup-file=no --project=. run_figures.jl
```

Puis lancer :

```bash
cd "/Users/jb270005/Desktop/LesHouchesGit"
julia --threads=auto --startup-file=no --project=. run_figures.jl
```

### Interactive terminal selection

To choose the comparison groups, figure families, individual plots, lines of
sight, snapshot window, output format, and destination without editing the
script, start the Julia terminal interface:

```bash
julia --threads=auto --startup-file=no --project=. \
    run_figures.jl --interactive
```

The equivalent Bash launcher is:

```bash
bash run_figures_interactive.sh
```

Paths below `/Xnfs` exist only on the cluster. To display the interactive menu
on the laptop while performing discovery, RAMSES conversion, and computation
remotely, use:

```bash
bash run_figures_interactive_cluster.sh
```

This allocates an SSH terminal on `PSMN_sr650node230`, runs the interactive
launcher in `~/Houches2026`, and downloads the figures after a successful run.
Set `DOWNLOAD_AFTER_RUN=false` to leave the results on the cluster, or override
the SSH alias with `CLUSTER_HOST`.

If `run_figures_interactive.sh` is started locally by mistake and an `/Xnfs`
path is entered, the Julia launcher now switches to this SSH workflow
automatically. The entered path is forwarded as the default data path in the
remote menu.

Selections use numbers: `1,3` selects two entries, `1-3` selects a range, and
`all` selects every displayed entry. Pressing Enter accepts the displayed
default. LOS-independent figures are always put in `common/` and computed only
once, even when several viewing axes are selected.

The first question selects the data layout. **Any directory** accepts either:

- the path of one simulation containing `DataCubes/` or snapshot files;
- a directory whose immediate subdirectories are simulations.

In the second case, the interface lists the simulations and lets the user
select one or several of them for comparative figures. No `VaryingMach`,
`VaryingRes`, or `VaryingRatio` suffix is added in this mode. The
**Preconfigured** option retains the three standard comparison folders.

Before asking which figures to compute, the interface scans every selected
simulation and prints a snapshot inventory: the resolved `DataCubes` directory,
the total number of HDF5/FITS or raw RAMSES snapshots, and the complete
numbered list of snapshot paths. A missing or empty cube directory is therefore
reported before any expensive calculation starts.

Raw RAMSES directories named `output_XXXXX` are supported through a persistent
yt conversion cache. Install the Python converter once on the cluster:

```bash
python3 -m pip install --user -r requirements-ramses.txt
```

When the selected simulation contains raw RAMSES outputs, the interface:

1. lists every detected `output_XXXXX`;
2. asks which outputs to cache;
3. asks for the uniform cube resolution;
4. converts them with yt into `.dynamo_cache/ramses_cubes/`;
5. runs the normal Julia analysis on those HDF5 cubes.

Each cache file contains density, pressure, three velocity components, three
cell-centred magnetic components, physical time, and box metadata. Existing
files are reused. A cache entry is rebuilt automatically when its RAMSES
`info_XXXXX.txt` fingerprint or requested resolution changes. Conversion uses
an atomic temporary file, so an interrupted job cannot leave a valid-looking
partial cube. Set `DYNAMO_PYTHON` if yt is installed in a different Python
environment.

The generated files are separated by comparison and snapshot window:

```text
figures/
├── varying_mach/
│   ├── dynamo_first20/
│   │   ├── common/
│   │   ├── los_x/
│   │   ├── los_y/
│   │   └── los_z/
│   └── observables_last10/
│       ├── common/
│       ├── los_x/
│       ├── los_y/
│       └── los_z/
├── varying_resolution/
│   └── ...
└── varying_ratio/
    └── ...
```

Le moteur batch :

- ne charge pas Pluto ;
- n'ouvre que les cubes nécessaires ;
- calcule l'union des dépendances des figures demandées une seule fois ;
- réutilise les cubes physiques et leurs champs dérivés avec des caches bornés
  en mémoire ;
- conserve dans `.dynamo_cache/` les réductions scientifiques compactes et un
  manifeste des snapshots pour les exécutions suivantes ;
- invalide automatiquement une entrée lorsque le fichier ou le dossier source
  est modifié ;
- limite séparément les workers de snapshots et les threads FFTW afin d'éviter
  la sur-souscription sur le système de fichiers partagé ;
- ferme les fichiers HDF5 après lecture ;
- affiche une barre de progression et la phase scientifique en cours ;
- écrit directement les figures PNG ou PDF dans `output_directory`.

Le cache persistant ne contient jamais les cubes bruts. Pour forcer un calcul
entièrement neuf, supprimer le dossier `.dynamo_cache/` avant le lancement.

À la fin, le script affiche la durée totale, le dossier de sortie absolu et le
chemin complet de chaque figure créée.

### Downloading the figures to the laptop

After the cluster computation finishes, run this command on the laptop from the
local repository:

```bash
cd "/Users/jb270005/Desktop/LesHouchesGit"
bash download_figures.sh
```

The script incrementally downloads:

```text
PSMN_sr650node230:~/Houches2026/figures/
```

into:

```text
/Users/jb270005/Desktop/LesHouchesGit/figures/
```

Existing unchanged figures are not transferred again. Custom locations can be
provided with `CLUSTER_HOST`, `REMOTE_FIGURE_DIRECTORY`, and
`LOCAL_FIGURE_DIRECTORY`.

For future runs, the complete cluster calculation and download can be started
from the laptop with one command:

```bash
bash run_cluster_and_download.sh
```

This command waits for a successful cluster computation before starting the
download.

La première simulation sert aux cartes non comparatives. Les diagnostics
comparatifs utilisent toutes les simulations listées.

Pour afficher la liste complète des noms de figures, laisser temporairement
`figures = String[]` dans `run_figures.jl`.

## Mise à jour du moteur batch

`src/BatchCellIndex.jl` contient seulement l'ordre et les dépendances des
cellules. Le code scientifique reste exclusivement dans
`dynamo_diagnostics.jl`.

Après une modification du notebook maître, régénérer l'index :

```bash
julia --startup-file=no --project=. tools/generate_batch_index.jl
```

Le batch vérifie l'empreinte du notebook et refuse d'utiliser un index périmé.

Pour mettre aussi à jour les six notebooks spécialisés après une modification
du notebook maître :

```bash
julia --startup-file=no --project=. tools/split_notebooks.jl
```

## Optimisations d'I/O

- Les répertoires et empreintes de snapshots sont mémorisés pendant la session.
- Un fichier HDF5 est parcouru une seule fois par chargement.
- Les composantes magnétiques centrées sont calculées en réutilisant un tampon.
- Le notebook interactif conserve par défaut un seul cube brut.
- Le batch peut conserver plusieurs simulations, dans la limite du nombre
  choisi et d'un plafond mémoire.
- Sur `/Xnfs`, le snapshot interactif peut être copié automatiquement sur le
  stockage temporaire local avant lecture.

Variables utiles :

| Variable | Rôle | Défaut |
|---|---|---|
| `DYNAMO_COMPARISON_REPOSITORY` | Parent de `VaryingMach`, `VaryingRes` et `VaryingRatio` | `/Xnfs/Houches2026/DynSim/cooling_freq_output` |
| `DYNAMO_DRY_RUN` | Affiche les jobs sans calculer | `false` |
| `DYNAMO_DATA_REPOSITORY` | Dossier des simulations | vide dans Pluto |
| `DYNAMO_LOCAL_HDF5_CACHE` | Cache local HDF5 : `auto`, `true`, `false` | `auto` |
| `DYNAMO_LOCAL_CACHE_DIRECTORY` | Parent du cache local | dossier temporaire |
| `DYNAMO_RAW_CUBE_CACHE_ENTRIES` | Nombre maximal de cubes en batch | nombre de simulations |
| `DYNAMO_RAW_CUBE_CACHE_MIB` | Plafond mémoire du cache | quart de la RAM |
| `DYNAMO_DERIVED_CACHE_MIB` | Plafond mémoire des champs dérivés | quart de la RAM |
| `DYNAMO_SNAPSHOT_WORKERS` | Nombre maximal de snapshots traités simultanément | `min(8, threads Julia)` |
| `DYNAMO_FFTW_THREADS` | Threads internes FFTW | `1` |
| `DYNAMO_CACHE_DIRECTORY` | Cache persistant des réductions et du manifeste | `.dynamo_cache` |
| `DYNAMO_PERSISTENT_CACHE_FILE` | Fichier de cache explicite | `scientific_cache_v1.bin` dans le dossier précédent |
| `PLUTO_HOST` | Adresse d'écoute | `127.0.0.1` |
| `PLUTO_PORT` | Port Pluto | `1234` |
| `PLUTO_LAUNCH_BROWSER` | Ouverture automatique du navigateur | oui sur macOS/Windows |

## Serveur distant

Exemple avec un tunnel SSH :

```bash
ssh -L 15432:127.0.0.1:15432 serveur
```

Puis sur le serveur :

```bash
cd /chemin/vers/LesHouchesGit
PLUTO_PORT=15432 julia --project=. run_pluto.jl
```

Ouvrir localement l'URL complète imprimée par Pluto, secret compris.

## Export HTML

```bash
julia --project=. export_html.jl
```

Pour choisir la destination :

```bash
DYNAMO_HTML_PATH="/chemin/resultat.html" \
julia --project=. export_html.jl
```

Pour exporter un notebook spécialisé :

```bash
DYNAMO_NOTEBOOK=dust.jl julia --project=. export_html.jl
```
