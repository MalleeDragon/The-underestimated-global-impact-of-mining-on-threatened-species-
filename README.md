# The underestimated global impact of mining on threatened species

Code repository for:

> **The underestimated global impact of mining on threatened species**
> *Cell Reports Sustainability* (accepted 2026)

---

## Repository structure

```
Step 1_Make binary mines map.R          # Rasterize mine polygons
Step 2_Create Richness Maps from IUCN data.R  # Species richness rasters
Step 3_Create plots for publication.R   # Publication figures
Model and plot code.R                   # GLMMs: richness in mining vs pseudo-mining areas
model summaries update.R               # GLMs: mining threat probability by taxonomic group
mining_analysis_complete.R             # Combined script: runs all five steps in sequence

Data_S1.xlsx    # Per-species range overlap and threat status data (input for model summaries)
Data S2.xlsx    # Species richness counts per buffer polygon (input for GLMMs)
```

Scripts assume the following working directory structure (not all included — see **Data availability** below):

```
data_raw/
├── mines/74548 mine polygons/74548_projected.shp
├── IUCN/Range Maps/MAMMALS/MAMMALS.shp
├── IUCN/Range Maps/REPTILES/REPTILES_PART1.shp, REPTILES_PART2.shp
├── IUCN/Range Maps/AMPHIBIANS/AMPHIBIANS_PART1.shp, AMPHIBIANS_PART2.shp
├── IUCN/Range Maps/BIRDS/birds_threat_3.shp
├── Species_threatened_list.csv
└── mining_threatened_zero_overlap_species.csv

data_clean/          # Created by Steps 1–2; read by Step 3
outputs/             # Figures written here by Steps 2–3
```

---

## Scripts

### Step 1 — Make binary mines map (`Step 1_Make binary mines map.R`)

Rasterizes 74,548 mine polygons at three resolutions (50 km, 10 km, 5 km) using the `terra` and `sf` packages.

**Input:** `data_raw/mines/74548 mine polygons/74548_projected.shp`
**Output:** `data_clean/mines_raster_50km.tif`, `mines_raster_10km.tif`, `mines_raster_5km.tif`
**Packages:** `terra`, `sf`, `tidyverse`

---

### Step 2 — Create species richness maps (`Step 2_Create Richness Maps from IUCN data.R`)

Reads IUCN terrestrial vertebrate range maps (mammals, reptiles, amphibians, birds), transforms them to Mollweide equal-area projection (ESRI:54009), and rasterizes them to calculate:
- Total species richness globally
- Richness of mining-threatened species (overall and per taxonomic group)
- Proportion of species threatened by mining per cell
- Richness of species listed as mining-threatened but with no spatial overlap with mines

**Input:** IUCN range shapefiles, `Species_threatened_list.csv`, `mining_threatened_zero_overlap_species.csv`
**Output:** Multiple `.tif` rasters in `data_clean/`; `data_clean/mapped_species.csv`
**Packages:** `sf`, `terra`, `tidyverse`, `here`

---

### Step 3 — Create publication figures (`Step 3_Create plots for publication.R`)

Produces global map figures (proportion and richness) for all vertebrates combined and for each taxonomic group separately.

**Input:** Rasters from `data_clean/`; `data_clean/World Outline/ne_110m_land.shp`
**Output:** PDF figures in `outputs/`
**Packages:** `terra`, `sf`, `ggplot2`, `viridis`, `cowplot`

---

### Model and plot code (`Model and plot code.R`)

Fits negative binomial GLMMs (via `glmmTMB`) comparing threatened and total species richness inside mining areas vs. size-matched pseudo-mining areas, for each taxonomic group and combined. Generates coefficient plots and predicted-effect plots.

**Input:** `Data S2.xlsx`
**Packages:** `readxl`, `glmmTMB`, `ggplot2`, `ggeffects`, `dplyr`, `tidyr`, `patchwork`, `scales`, `gridExtra`, `stringr`

---

### Model summaries — threat probability (`model summaries update.R`)

Fits binomial GLMs to model the probability of a species having mining listed as a threat, as a function of the proportion of range overlapping with mines, separately for each taxonomic group.

**Input:** `Data_S1.xlsx` (columns used: `Group`, `Existing overlap proportion`, `Mining threat status (IUCN)`)
**Output:** `glm_prob_by_group_excluding_all_species.png`
**Packages:** `readxl`, `dplyr`, `ggplot2`

---

## Data availability

The supplementary data files (`Data_S1.xlsx` and `Data S2.xlsx`) required to reproduce the statistical models are included in this repository.

Large spatial datasets used in Steps 1–3 are not included due to size and licensing:

- **Mine polygons** (74,548 polygons): Tang, L. & Werner, T.T. (2023). Global mining footprint mapped from high-resolution satellite imagery. *Communications Earth & Environment*, 4, 1–12. https://doi.org/10.1038/s43247-023-00805-6. Freely available from the journal supplementary data.
- **IUCN Red List range maps**: available from the [IUCN Spatial Data Download portal](https://www.iucnredlist.org/resources/spatial-data-download) (free registration required).
- **World outline shapefile**: Natural Earth 1:110m land polygons, available from [naturalearthdata.com](https://www.naturalearthdata.com/downloads/110m-physical-vectors/).

---

## R session info

Scripts were developed and run in R 4.4–4.5. Key packages and versions used:

| Package   | Version |
|-----------|---------|
| terra     | ≥ 1.7   |
| sf        | ≥ 1.0   |
| glmmTMB   | ≥ 1.1   |
| ggplot2   | ≥ 3.4   |
| tidyverse | ≥ 2.0   |
| viridis   | ≥ 0.6   |
| cowplot   | ≥ 1.1   |
| patchwork | ≥ 1.2   |
| ggeffects | ≥ 1.5   |
| here      | ≥ 1.0   |

---

## Citation

Please cite the associated paper if you use this code.
