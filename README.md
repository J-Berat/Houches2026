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
For each selected comparison, the script creates two jobs:

- the first 20 snapshots for all 24 Dynamo figures, including separate 3D and
  projected 2D HRO outputs and a seven-panel summary figure;
- the last 10 snapshots for all 23 figures from Dust, StarlightPol, ZEEMAN,
  MOOSE, and SHINE.

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

Check all six planned jobs without computing:

```bash
DYNAMO_DRY_RUN=true \
julia --threads=auto --startup-file=no --project=. run_figures.jl
```

Puis lancer :

```bash
cd "/Users/jb270005/Desktop/LesHouchesGit"
julia --threads=auto --startup-file=no --project=. run_figures.jl
```

The generated files are separated by comparison and snapshot window:

```text
figures/
├── varying_mach/
│   ├── dynamo_first20/
│   └── observables_last10/
├── varying_resolution/
│   ├── dynamo_first20/
│   └── observables_last10/
└── varying_ratio/
    ├── dynamo_first20/
    └── observables_last10/
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
