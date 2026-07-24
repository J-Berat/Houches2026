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
├── Manifest-v1.11.toml
├── Manifest-v1.12.toml
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

Ouvrir `run_figures.jl` et modifier uniquement le bloc `CONFIG` :

```julia
const CONFIG = BatchConfig(
    data_repository = "/Xnfs/Houches2026/DynSim",
    simulations = [
        "simulation_1",
        "simulation_2",
        "simulation_3",
    ],
    snapshot = :last,
    line_of_sight = "z",
    figures = figures_for_notebooks(["dynamo", "dust"]),
    output_directory = joinpath(PROJECT_DIRECTORY, "figures", "dynamo_dust"),
    output_format = "png",
)
```

The notebook groups are `dynamo`, `dust`, `starlightpol`, `zeeman`, `moose`,
and `shine`. Selecting `dynamo` and `dust` computes all 23 figures from those
two notebooks.

Puis lancer :

```bash
cd "/Users/jb270005/Desktop/LesHouchesGit"
julia --threads=auto --startup-file=no --project=. run_figures.jl
```

Le moteur batch :

- ne charge pas Pluto ;
- n'ouvre que les cubes nécessaires ;
- calcule l'union des dépendances des figures demandées une seule fois ;
- réutilise les cubes avec un cache LRU borné en mémoire ;
- ferme les fichiers HDF5 après lecture ;
- affiche une barre de progression et la phase scientifique en cours ;
- écrit directement les figures PNG ou PDF dans `output_directory`.

À la fin, le script affiche la durée totale, le dossier de sortie absolu et le
chemin complet de chaque figure créée.

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
| `DYNAMO_DATA_REPOSITORY` | Dossier des simulations | vide dans Pluto |
| `DYNAMO_LOCAL_HDF5_CACHE` | Cache local HDF5 : `auto`, `true`, `false` | `auto` |
| `DYNAMO_LOCAL_CACHE_DIRECTORY` | Parent du cache local | dossier temporaire |
| `DYNAMO_RAW_CUBE_CACHE_ENTRIES` | Nombre maximal de cubes en batch | nombre de simulations |
| `DYNAMO_RAW_CUBE_CACHE_MIB` | Plafond mémoire du cache | quart de la RAM |
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
