# ==============================================================================
# The underestimated global impact of mining on threatened species
# Complete analysis script
# ==============================================================================
#
# This script reproduces all analyses and figures in the paper.
# It combines five separate scripts in the order they must be run:
#
#   1. Rasterize mine polygons                    (Steps 1)
#   2. Create species richness rasters            (Step 2)
#   3. Create publication figures from rasters    (Step 3)
#   4. GLMMs: richness in mining vs pseudo-mining (Model and plot code)
#   5. GLMs: mining threat probability by taxon   (model summaries update)
#
# ------------------------------------------------------------------------------
# REQUIRED DATA
# ------------------------------------------------------------------------------
#
# Scripts 1-3 require large spatial datasets that are NOT included in this
# repository due to file size and licensing. These must be obtained separately
# and placed in the subdirectory structure described below.
#
# Mine polygons (Script 1 and 2):
#   Source: S&P Global Market Intelligence (commercial subscription required)
#   Place at: data_raw/mines/74548 mine polygons/74548_projected.shp
#
# IUCN terrestrial vertebrate range maps (Script 2):
#   Source: IUCN Red List Spatial Data Download (free, registration required)
#           https://www.iucnredlist.org/resources/spatial-data-download
#   Place at: data_raw/IUCN/Range Maps/MAMMALS/MAMMALS.shp
#             data_raw/IUCN/Range Maps/REPTILES/REPTILES_PART1.shp
#             data_raw/IUCN/Range Maps/REPTILES/REPTILES_PART2.shp
#             data_raw/IUCN/Range Maps/AMPHIBIANS/AMPHIBIANS_PART1.shp
#             data_raw/IUCN/Range Maps/AMPHIBIANS/AMPHIBIANS_PART2.shp
#             data_raw/IUCN/Range Maps/BIRDS/birds_threat_3.shp
#
# IUCN species list (Script 2):
#   Place at: data_raw/Species_threatened_list.csv
#   Columns required: internalTaxonId, scientificName, className
#
# Species with no spatial overlap with mines (Script 2):
#   Place at: data_raw/mining_threatened_zero_overlap_species.csv
#   Columns required: species
#
# World outline shapefile (Script 3):
#   Source: Natural Earth 1:110m land polygons
#           https://www.naturalearthdata.com/downloads/110m-physical-vectors/
#   Place at: data_clean/World Outline/ne_110m_land.shp
#
# Scripts 4-5 use the supplementary data files included in this repository:
#   Data S2.xlsx    — species richness counts per buffer polygon (Script 4)
#   Data_S1.xlsx    — per-species range overlap and threat status data (Script 5)
#
# ------------------------------------------------------------------------------
# OUTPUT DIRECTORIES
# ------------------------------------------------------------------------------
#
# Scripts 1-2 write intermediate rasters to: data_clean/
# Script 3 writes figures to: outputs/
# Script 5 writes a figure to the working directory.
#
# Create these directories before running if they do not exist:
#   dir.create("data_clean", showWarnings = FALSE)
#   dir.create("outputs", showWarnings = FALSE)
#
# ------------------------------------------------------------------------------
# PACKAGES
# ------------------------------------------------------------------------------
#
# Install all required packages with:
#   install.packages(c("terra", "sf", "tidyverse", "here",
#                      "ggplot2", "viridis", "cowplot",
#                      "glmmTMB", "ggeffects", "patchwork",
#                      "scales", "gridExtra", "stringr", "dplyr"))
#
# ==============================================================================


# ==============================================================================
# SCRIPT 1: Rasterize mine polygons
# ==============================================================================
#
# Reads the projected mine polygon shapefile and rasterizes it at three
# resolutions (50 km, 10 km, 5 km). Output rasters are binary (1 = mine
# present, 0 = no mine).
#
# Input:  data_raw/mines/74548 mine polygons/74548_projected.shp
# Output: data_clean/mines_raster_50km.tif
#         data_clean/mines_raster_10km.tif
#         data_clean/mines_raster_5km.tif

library(terra)
library(sf)
library(tidyverse)

mines = read_sf(here::here("data_raw", "mines", "74548 mine polygons", "74548_projected.shp"))

# Create template rasters at each resolution
res_km = 50
base = rast(extent = ext(mines),
            resolution = res_km*1000,
            vals = 1, crs = crs(mines))

res_10 = 10
base_10 = rast(extent = ext(mines),
               resolution = res_10*1000,
               vals = 1, crs = crs(mines))

res_5 = 5
base_5 = rast(extent = ext(mines),
              resolution = res_5*1000,
              vals = 1, crs = crs(mines))

# Rasterize: 1 where a mine polygon touches a cell, 0 otherwise
mines$dummy = 1
mines_raster_50 = rasterize(mines, base,    field = "dummy", background = 0, touches = TRUE, fun = 'max')
mines_raster_10 = rasterize(mines, base_10, field = "dummy", background = 0, touches = TRUE, fun = 'max')
mines_raster_5  = rasterize(mines, base_5,  field = "dummy", background = 0, touches = TRUE, fun = 'max')

terra::writeRaster(mines_raster_50, filename = "data_clean/mines_raster_50km.tif", overwrite = TRUE)
terra::writeRaster(mines_raster_10, filename = "data_clean/mines_raster_10km.tif", overwrite = TRUE)
terra::writeRaster(mines_raster_5,  filename = "data_clean/mines_raster_5km.tif",  overwrite = TRUE)


# ==============================================================================
# SCRIPT 2: Create species richness rasters from IUCN range maps
# ==============================================================================
#
# Reads IUCN terrestrial vertebrate range polygons for mammals, reptiles,
# amphibians, and birds. Transforms all to Mollweide equal-area projection
# (ESRI:54009) and rasterizes to the 50 km mine grid to calculate:
#   - Overall species richness
#   - Richness of mining-threatened species (total and per taxon)
#   - Proportion of species threatened by mining per cell
#   - Richness of species listed as mining-threatened with no spatial overlap
#
# Input:  IUCN range shapefiles (see header for paths)
#         data_raw/Species_threatened_list.csv
#         data_raw/mining_threatened_zero_overlap_species.csv
# Output: Multiple .tif rasters in data_clean/
#         data_clean/mapped_species.csv

library(sf)
library(terra)
library(tidyverse)

# Read IUCN range maps
mammals    = read_sf("data_raw/IUCN/Range Maps/MAMMALS/MAMMALS.shp")
reptiles1  = read_sf("data_raw/IUCN/Range Maps/REPTILES/REPTILES_PART1.shp")
reptiles2  = read_sf("data_raw/IUCN/Range Maps/REPTILES/REPTILES_PART2.shp")
amphibians1 = read_sf("data_raw/IUCN/Range Maps/AMPHIBIANS/AMPHIBIANS_PART1.shp")
amphibians2 = read_sf("data_raw/IUCN/Range Maps/AMPHIBIANS/AMPHIBIANS_PART2.shp")
birds      = read_sf("data_raw/IUCN/Range Maps/BIRDS/birds_threat_3.shp")

# Align column names and combine non-bird vertebrates
names(mammals) <- names(reptiles2)
species = rbind(mammals, reptiles1, reptiles2, amphibians1, amphibians2)

# Filter to terrestrial species and retain key columns
species = species %>% filter(terrestial == "true") %>% select(sci_name, class)

# Reproject to Mollweide equal-area projection
crs.set = "ESRI:54009"
species = species %>% st_transform(st_crs(crs.set))

# Add birds
birds = birds %>%
  select(SCI_NAME) %>% rename("sci_name" = "SCI_NAME") %>%
  mutate(class = "AVES") %>%
  st_transform(st_crs(species))
species = rbind(species, birds)

# Free memory
rm(amphibians1, amphibians2, reptiles1, reptiles2, mammals, birds)
gc()

# Create 50 km mine raster in Mollweide (used as the analysis grid)
mines = read_sf(here::here("data_raw", "mines", "74548 mine polygons", "74548_projected.shp")) %>%
  st_transform(st_crs(species))
base = rast(extent = ext(mines), resolution = 50000, vals = 1, crs = crs(mines))
mines$dummy = 1
mines_raster_50 = rasterize(mines, base, field = "dummy", background = 0, touches = TRUE, fun = 'max')
mines = mines_raster_50
writeRaster(mines, "data_clean/mines_raster_50km.tif", overwrite = TRUE)

# Read IUCN species list and identify mining-threatened species
spp.assessed = read.csv("data_raw/Species_threatened_list.csv")
spp.name     = unique(spp.assessed$scientificName)
spp.assessed.taxa = filter(spp.assessed, className %in% c("AVES", "REPTILIA", "AMPHIBIA", "MAMMALIA"))

# Filter range maps to mining-threatened species
species$dummy     = 1
species.mining    = species %>% filter(sci_name %in% spp.name)

# Record which species have spatial data vs. not
species.mining.spatial = data.frame(sci_name = species.mining$sci_name, mapped = 1) %>% distinct()
mappedvsnot = left_join(spp.assessed.taxa, species.mining.spatial, by = c("scientificName" = "sci_name"))
mappedvsnot$mapped = ifelse(is.na(mappedvsnot$mapped), 0, mappedvsnot$mapped)
xtabs(~mappedvsnot$mapped)
write.csv(mappedvsnot, "data_clean/mapped_species.csv")

# Calculate richness rasters
species.rich        = rasterize(species,        mines, field = "dummy", fun = 'sum')
species.mines.rich  = rasterize(species.mining, mines, field = "dummy", fun = 'sum')
names(species.rich)       = "SR_species"
names(species.mines.rich) = "SR_species_Mine"

# Proportion of all species that are mining-threatened
prop.species.mine           = species.mines.rich / species.rich
names(prop.species.mine)    = "Prop_species_Mine"

# Clip richness and proportion to mine cells only
species.rich.mines.clip     = species.mines.rich * mines
names(species.rich.mines.clip) = "SR_species_Mine_Clip"
prop.species.mines.clip     = prop.species.mine * mines
names(prop.species.mines.clip) = "Prop_species_Mine_Clip"

writeRaster(species.rich,            "data_clean/species_richness_global.tif",             overwrite = TRUE)
writeRaster(species.rich.mines.clip, "data_clean/species_richness_minethreat_global.tif",  overwrite = TRUE)
writeRaster(prop.species.mine,       "data_clean/prop_species_minethreat.tif",             overwrite = TRUE)
writeRaster(prop.species.mine,       "data_clean/prop_species_minethreat_clip.tif",        overwrite = TRUE)

# Per-taxon richness and proportions
mammals.rich    = rasterize(filter(species, class == "MAMMALIA"), mines, field = "dummy", fun = 'sum')
mammals.mines.rich = rasterize(filter(species.mining, class == "MAMMALIA"), mines, field = "dummy", fun = 'sum')
prop.mammals.mines = mammals.mines.rich / mammals.rich
names(prop.mammals.mines) = "Prop_Mammals_MineThreat"

birds.rich      = rasterize(filter(species, class == "AVES"), mines, field = "dummy", fun = 'sum')
birds.mines.rich = rasterize(filter(species.mining, class == "AVES"), mines, field = "dummy", fun = 'sum')
prop.birds.mines = birds.mines.rich / birds.rich
names(prop.birds.mines) = "Prop_Birds_MineThreat"

reptiles.rich   = rasterize(filter(species, class == "REPTILIA"), mines, field = "dummy", fun = 'sum')
reptiles.mines.rich = rasterize(filter(species.mining, class == "REPTILIA"), mines, field = "dummy", fun = 'sum')
prop.reptiles.mine = reptiles.mines.rich / reptiles.rich
names(prop.reptiles.mine) = "Prop_Reptiles_MineThreat"

amphibians.rich = rasterize(filter(species, class == "AMPHIBIA"), mines, field = "dummy", fun = 'sum')
amphibians.mines.rich = rasterize(filter(species.mining, class == "AMPHIBIA"), mines, field = "dummy", fun = 'sum')
prop.amphibians.mine = amphibians.mines.rich / amphibians.rich
names(prop.amphibians.mine) = "Prop_Amphibians_MineThreat"

writeRaster(prop.mammals.mines,   "data_clean/prop_mammals_mines.tif",   overwrite = TRUE)
writeRaster(prop.birds.mines,     "data_clean/prop_birds_mines.tif",     overwrite = TRUE)
writeRaster(prop.reptiles.mine,   "data_clean/prop_reptiles_mines.tif",  overwrite = TRUE)
writeRaster(prop.amphibians.mine, "data_clean/prop_amphibians_mines.tif", overwrite = TRUE)

writeRaster(mammals.mines.rich,   "data_clean/rich_mammals_mines.tif",   overwrite = TRUE)
writeRaster(birds.mines.rich,     "data_clean/rich_birds_mines.tif",     overwrite = TRUE)
writeRaster(reptiles.mines.rich,  "data_clean/rich_reptiles_mines.tif",  overwrite = TRUE)
writeRaster(amphibians.mines.rich,"data_clean/rich_amphibians_mines.tif", overwrite = TRUE)

# Clipped versions (mine cells only)
prop.mammals.mines.clip   = prop.mammals.mines   * mines
prop.birds.mines.clip     = prop.birds.mines     * mines
prop.reptiles.mines.clip  = prop.reptiles.mine   * mines
prop.amphibians.mines.clip = prop.amphibians.mine * mines

writeRaster(prop.mammals.mines.clip,   "data_clean/prop_mammals_mines_clip.tif",   overwrite = TRUE)
writeRaster(prop.birds.mines.clip,     "data_clean/prop_birds_mines_clip.tif",     overwrite = TRUE)
writeRaster(prop.reptiles.mines.clip,  "data_clean/prop_reptiles_mines_clip.tif",  overwrite = TRUE)
writeRaster(prop.amphibians.mines.clip,"data_clean/prop_amphibians_mines_clip.tif", overwrite = TRUE)

# Species listed as mining-threatened but with no spatial overlap with mines
species_nooverlap         = read.csv("data_raw/mining_threatened_zero_overlap_species.csv")
species_nooverlap_spatial = species.mining %>% filter(sci_name %in% species_nooverlap$species)
species.rich.nooverlap    = rasterize(species_nooverlap_spatial, mines, field = "dummy", fun = 'sum')
names(species.rich.nooverlap) = "SR_species_Mine_Overlap"
writeRaster(species.rich.nooverlap, "data_clean/species_rich_mines_nooverlap.tif", overwrite = TRUE)


# ==============================================================================
# SCRIPT 3: Create publication figures from rasters
# ==============================================================================
#
# Reads rasters from data_clean/ and produces global map figures for:
#   - Proportion and richness of mining-threatened species (all vertebrates)
#   - Per-taxon richness and proportion maps
#   - Richness of species listed as mining-threatened with no spatial overlap
#
# Input:  Rasters in data_clean/ (output of Script 2)
#         data_clean/World Outline/ne_110m_land.shp
# Output: PDF figures in outputs/

library(cowplot)
library(terra)
library(ggplot2)
library(sf)
library(viridis)

# Read world outline and 50 km mine raster
world = read_sf("data_clean/World Outline/ne_110m_land.shp") %>% st_transform("ESRI:54009")
mines = rast("data_clean/mines_raster_50km.tif")

# ---- All vertebrates: proportion threatened by mining ----
species.prop = rast("data_clean/prop_species_minethreat_clip.tif") * mines
species.rich = rast("data_clean/species_richness_minethreat_global.tif") * mines
species.prop[species.prop == 0] <- NA
species.rich[species.rich == 0] <- NA

breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 1)
labels <- c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", ">0.4")

species.prop.dat = as.data.frame(species.prop, xy = TRUE)
species.prop.dat$Prop_species_Mine_bin <- cut(species.prop.dat$Prop_species_Mine,
                                              breaks = breaks, labels = labels,
                                              include.lowest = TRUE)

global_prop_plot <- ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "lightgrey", alpha = 0.5) +
  geom_point(data = species.prop.dat, aes(x = x, y = y, colour = Prop_species_Mine_bin,
                                          size = Prop_species_Mine_bin)) +
  scale_colour_viridis_d(name = "Proportion", direction = 1, labels = labels, guide = "legend") +
  scale_size_manual(name = "Proportion", values = c(0.1, 1, 2.5, 4, 5),
                    labels = labels, guide = "legend") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(filename = "outputs/prop_alltaxa_plot.pdf", plot = global_prop_plot, height = 4, width = 8)

# ---- All vertebrates: richness ----
breaks <- c(0, 10, 20, 30, 40, 100)
labels <- c("0-10", "10-20", "20-30", "30-40", ">40")

species.rich.dat = as.data.frame(species.rich, xy = TRUE)
species.rich.dat$SR_species_Mine_Clip_bin <- cut(species.rich.dat$SR_species_Mine_Clip,
                                                  breaks = breaks, labels = labels,
                                                  include.lowest = TRUE)

global_rich_plot <- ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "lightgrey", alpha = 0.5) +
  geom_point(data = species.rich.dat, aes(x = x, y = y, colour = SR_species_Mine_Clip_bin,
                                           size = SR_species_Mine_Clip_bin)) +
  scale_colour_viridis_d(name = "Richness", direction = 1, labels = labels, guide = "legend") +
  scale_size_manual(name = "Richness", values = c(0.1, 1, 2.5, 4, 5),
                    labels = labels, guide = "legend") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(filename = "outputs/richness_alltaxa_plot.pdf", plot = global_rich_plot, height = 4, width = 8)

# ---- Per-taxon richness maps ----
mines = rast("data_clean/mines_raster_50km.tif")

mammals.rich   = rast("data_clean/rich_mammals_mines.tif")   * mines
birds.rich     = rast("data_clean/rich_birds_mines.tif")     * mines
reptiles.rich  = rast("data_clean/rich_reptiles_mines.tif")  * mines
amphibians.rich = rast("data_clean/rich_amphibians_mines.tif") * mines

mammals.rich[mammals.rich == 0]     <- NA
birds.rich[birds.rich == 0]         <- NA
reptiles.rich[reptiles.rich == 0]   <- NA
amphibians.rich[amphibians.rich == 0] <- NA

# Breaks/labels for richness (reused across taxa)
breaks <- c(0, 10, 20, 30, 40, 100)
labels <- c("0-10", "10-20", "20-30", "30-40", ">40")

make_richness_map <- function(raster_obj, layer_col, output_file) {
  dat <- as.data.frame(raster_obj, xy = TRUE)
  dat$bin <- cut(dat[[layer_col]], breaks = breaks, labels = labels, include.lowest = TRUE)
  p <- ggplot() +
    geom_sf(data = world, fill = "lightgrey", color = "lightgrey", alpha = 0.5) +
    geom_point(data = dat, aes(x = x, y = y, colour = bin, size = bin)) +
    scale_colour_viridis_d(name = "Richness", direction = 1, labels = labels, guide = "legend") +
    scale_size_manual(name = "Richness", values = c(0.1, 1, 2.5, 4, 5),
                      labels = labels, guide = "legend") +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  ggsave(filename = output_file, plot = p, height = 4, width = 8)
}

make_richness_map(mammals.rich,    "dummy", "outputs/richness_mammals_plot.pdf")
make_richness_map(birds.rich,      "dummy", "outputs/richness_birds_plot.pdf")
make_richness_map(amphibians.rich, "dummy", "outputs/richness_amphibians_plot.pdf")
make_richness_map(reptiles.rich,   "dummy", "outputs/richness_reptiles_plot.pdf")

# ---- Per-taxon proportion maps ----
mammals.prop   = rast("data_clean/prop_mammals_mines_clip.tif")
birds.prop     = rast("data_clean/prop_birds_mines_clip.tif")
reptiles.prop  = rast("data_clean/prop_reptiles_mines_clip.tif")
amphibians.prop = rast("data_clean/prop_amphibians_mines_clip.tif")

mammals.prop[mammals.prop == 0]     <- NA
birds.prop[birds.prop == 0]         <- NA
reptiles.prop[reptiles.prop == 0]   <- NA
amphibians.prop[amphibians.prop == 0] <- NA

breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 1)
labels <- c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", ">0.4")

make_prop_map <- function(raster_obj, layer_col, output_file) {
  dat <- as.data.frame(raster_obj, xy = TRUE)
  dat$bin <- cut(dat[[layer_col]], breaks = breaks, labels = labels, include.lowest = TRUE)
  p <- ggplot() +
    geom_sf(data = world, fill = "lightgrey", color = "lightgrey", alpha = 0.5) +
    geom_point(data = dat, aes(x = x, y = y, colour = bin, size = bin)) +
    scale_colour_viridis_d(name = "Proportion", direction = 1, labels = labels, guide = "legend") +
    scale_size_manual(name = "Proportion", values = c(0.1, 1, 2.5, 4, 5),
                      labels = labels, guide = "legend") +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  ggsave(filename = output_file, plot = p, height = 4, width = 8)
}

make_prop_map(mammals.prop,   "Prop_Mammals_MineThreat",   "outputs/prop_mammals_plot.pdf")
make_prop_map(birds.prop,     "Prop_Birds_MineThreat",     "outputs/prop_birds_plot.pdf")
make_prop_map(amphibians.prop,"Prop_Amphibians_MineThreat","outputs/prop_amphibian_plot.pdf")
make_prop_map(reptiles.prop,  "Prop_Reptiles_MineThreat",  "outputs/prop_reptiles_plot.pdf")

# ---- Species listed as mining-threatened with no spatial overlap ----
no.overlap = rast("data_clean/species_rich_mines_nooverlap.tif")
no.overlap[no.overlap == 0] <- NA

breaks <- c(0, 1, 2, 3, 4, 5)
labels <- c("1", "2", "3", "4", "5")

no.overlap.dat = as.data.frame(no.overlap, xy = TRUE)
no.overlap.dat$bin <- cut(no.overlap.dat$SR_species_Mine_Overlap,
                           breaks = breaks, labels = labels, include.lowest = TRUE)

global_rich_nooverlap <- ggplot() +
  geom_sf(data = world, fill = "lightgrey", color = "lightgrey", alpha = 0.5) +
  geom_tile(data = no.overlap.dat, aes(x = x, y = y, fill = SR_species_Mine_Overlap)) +
  scale_fill_viridis_c(name = "Richness", direction = 1, guide = "legend") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

ggsave(filename = "outputs/richness_allspecies_nominesmapped_plot.pdf",
       plot = global_rich_nooverlap, height = 4, width = 8)


# ==============================================================================
# SCRIPT 4: GLMMs — species richness in mining vs. pseudo-mining areas
# ==============================================================================
#
# Fits negative binomial GLMMs (glmmTMB, nbinom2 family) comparing threatened
# and total species richness inside mining areas against size-matched
# pseudo-mining control areas, for each taxonomic group and combined.
# Generates predicted-effect plots and coefficient (forest) plots.
#
# Input:  buffer_final_counts_with_all_birds.csv
#         Columns: FID, Treatment (mines_10km / pseudo_10km),
#                  amphibians_threatened_count, reptiles_threatened_count,
#                  mammals_threatened_count, birds_threatened_count,
#                  sum_threatened_count, [total count equivalents]
# Output: Plots displayed (not saved to file — save manually as needed)

rm(list = ls())

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(glmmTMB)
library(ggeffects)
library(gridExtra)
library(patchwork)
library(scales)
library(readxl)

polygons_species_richness <- read_excel("Data S2.xlsx")

# Set reference level: pseudo-mining areas as baseline
polygons_species_richness$Treatment = factor(polygons_species_richness$Treatment,
                                              levels = c("pseudo_10km", "mines_10km"))
polygons_species_richness$X <- paste0(polygons_species_richness$Treatment, "_",
                                       polygons_species_richness$FID)

# ---- GLMMs: threatened species richness ----
threatened_amphibians_model <- glmmTMB(amphibians_threatened_count ~ Treatment + (1|FID),
                                        data = polygons_species_richness, family = nbinom2())
threatened_reptiles_model   <- glmmTMB(reptiles_threatened_count   ~ Treatment + (1|FID),
                                        data = polygons_species_richness, family = nbinom2())
threatened_mammals_model    <- glmmTMB(mammals_threatened_count    ~ Treatment + (1|FID),
                                        data = polygons_species_richness, family = nbinom2())
threatened_birds_model      <- glmmTMB(birds_threatened_count      ~ Treatment + (1|FID),
                                        data = polygons_species_richness, family = nbinom2())
threatened_total_model      <- glmmTMB(sum_threatened_count        ~ Treatment + (1|FID),
                                        data = polygons_species_richness, family = nbinom2())

# Extract model summaries into a tidy table
extract_and_reshape_summary <- function(model, group) {
  coef_summary <- summary(model)$coef$cond
  intercept_row <- data.frame(
    Group = group, Parameter = "Intercept",
    Estimate = coef_summary["(Intercept)",       "Estimate"],
    SE       = coef_summary["(Intercept)",       "Std. Error"],
    z_value  = coef_summary["(Intercept)",       "z value"],
    Pr       = coef_summary["(Intercept)",       "Pr(>|z|)"]
  )
  treatment_row <- data.frame(
    Group = group, Parameter = "Treatment",
    Estimate = coef_summary["Treatmentmines_10km", "Estimate"],
    SE       = coef_summary["Treatmentmines_10km", "Std. Error"],
    z_value  = coef_summary["Treatmentmines_10km", "z value"],
    Pr       = coef_summary["Treatmentmines_10km", "Pr(>|z|)"]
  )
  rbind(intercept_row, treatment_row)
}

reshaped_results <- rbind(
  extract_and_reshape_summary(threatened_total_model,      "Total"),
  extract_and_reshape_summary(threatened_amphibians_model, "Amphibians"),
  extract_and_reshape_summary(threatened_reptiles_model,   "Reptiles"),
  extract_and_reshape_summary(threatened_mammals_model,    "Mammals"),
  extract_and_reshape_summary(threatened_birds_model,      "Birds")
)
print(reshaped_results)

# ---- Predicted-effect plots: threatened richness ----
predicted_amphibians <- ggpredict(threatened_amphibians_model, terms = "Treatment", ci.lvl = 0.95)
predicted_reptiles   <- ggpredict(threatened_reptiles_model,   terms = "Treatment", ci.lvl = 0.95)
predicted_mammals    <- ggpredict(threatened_mammals_model,    terms = "Treatment", ci.lvl = 0.95)
predicted_birds      <- ggpredict(threatened_birds_model,      terms = "Treatment", ci.lvl = 0.95)
predicted_total      <- ggpredict(threatened_total_model,      terms = "Treatment", ci.lvl = 0.95)

set_y_breaks <- function(data, n = 4) {
  seq(min(data$conf.low), max(data$conf.high), length.out = n)
}

make_predicted_plot <- function(predicted_values, dataset_name, mining_col = "#56B4E9", pseudo_col = "#D55E00") {
  predicted_data <- as.data.frame(predicted_values)
  predicted_data$x <- ifelse(predicted_data$x == "mines_10km", "Mining", "Pseudo")
  predicted_data$x <- factor(predicted_data$x, levels = c("Mining", "Pseudo"))
  y_breaks <- set_y_breaks(predicted_data, n = 4)
  ggplot(predicted_data, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = x)) +
    geom_errorbar(width = 0.1, linewidth = 1.5) +
    geom_point(size = 6, shape = 18) +
    scale_color_manual(values = c("Mining" = mining_col, "Pseudo" = pseudo_col)) +
    scale_y_continuous(breaks = y_breaks, labels = label_number(scale = 1, accuracy = 0.01)) +
    labs(title = dataset_name, x = "Area", y = "Species richness") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5), axis.line = element_line(linewidth = 1),
          axis.ticks = element_line(size = 1), axis.text = element_text(size = 14),
          axis.title = element_text(size = 14), panel.grid = element_blank(),
          legend.position = "none")
}

plot_total_t      <- make_predicted_plot(predicted_total,      "Total")
plot_amphibians_t <- make_predicted_plot(predicted_amphibians, "Amphibian") + theme(axis.title.x = element_blank())
plot_reptiles_t   <- make_predicted_plot(predicted_reptiles,   "Reptile")   + theme(axis.title.x = element_blank())
plot_mammals_t    <- make_predicted_plot(predicted_mammals,    "Mammal")
plot_birds_t      <- make_predicted_plot(predicted_birds,      "Bird")

# Coefficient (forest) plot for threatened richness
extract_coefficients <- function(model, model_name) {
  s <- summary(model)$coeff$cond %>% as.data.frame()
  s$term      <- rownames(s)
  s$conf.low  <- s$Estimate - 1.96 * s$`Std. Error`
  s$conf.high <- s$Estimate + 1.96 * s$`Std. Error`
  s$color     <- "black"
  s$color[s$term %in% c("(Intercept)", "Intercept (Outside)")] <- "gray"
  sel <- s$term == "Treatmentmines_10km" & s$`Pr(>|z|)` < 0.05
  s$color[sel & s$Estimate > 0] <- "#56B4E9"
  s$color[sel & s$Estimate < 0] <- "#FF7F7F"
  s$model <- model_name
  s
}

new_all_coefficients <- bind_rows(
  extract_coefficients(threatened_amphibians_model, "Amphibian"),
  extract_coefficients(threatened_reptiles_model,   "Reptile"),
  extract_coefficients(threatened_mammals_model,    "Mammal"),
  extract_coefficients(threatened_birds_model,      "Bird"),
  extract_coefficients(threatened_total_model,      "Total")
)
new_all_coefficients$model <- factor(new_all_coefficients$model,
                                      levels = c("Total", "Amphibian", "Reptile", "Mammal", "Bird"))

plot_pseudo_coefficient <- function(coefficients, x_color_pos = "#56B4E9", x_color_neg = "#FF7F7F") {
  pseudo <- coefficients %>%
    filter(term == "Treatmentmines_10km") %>%
    mutate(model = factor(model, levels = rev(levels(model))))
  ggplot(pseudo, aes(x = Estimate, y = model, color = color)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, linewidth = 1.5) +
    geom_point(size = 6, shape = 18) +
    scale_color_identity() +
    scale_x_continuous(labels = label_number(scale = 1, accuracy = 0.01)) +
    labs(x = "Mining area \u03b2", y = "Animal class") +
    theme_minimal() +
    theme(axis.line = element_line(linewidth = 1), axis.ticks = element_line(size = 1),
          axis.text = element_text(size = 14), axis.title = element_text(size = 14),
          panel.grid = element_blank(), legend.position = "none")
}

pseudo_plot_threatened <- plot_pseudo_coefficient(new_all_coefficients)

patchwork::wrap_plots(
  plot_total_t,      plot_amphibians_t,
  plot_reptiles_t,   plot_mammals_t,
  plot_birds_t,      pseudo_plot_threatened,
  ncol = 2, widths = c(1, 1), heights = c(1, 1, 1)
)

# ---- GLMMs: total (all) species richness ----
total_amphibians_model <- glmmTMB(amphibians_total_count ~ Treatment + (1|FID),
                                   data = polygons_species_richness, family = nbinom2())
total_reptiles_model   <- glmmTMB(reptiles_total_count   ~ Treatment + (1|FID),
                                   data = polygons_species_richness, family = nbinom2())
total_mammals_model    <- glmmTMB(mammals_total_count    ~ Treatment + (1|FID),
                                   data = polygons_species_richness, family = nbinom2())
total_birds_model      <- glmmTMB(birds_total_count      ~ Treatment + (1|FID),
                                   data = polygons_species_richness, family = nbinom2())
total_total_model      <- glmmTMB(sum_total_count        ~ Treatment + (1|FID),
                                   data = polygons_species_richness, family = nbinom2())

reshaped_results_total <- rbind(
  extract_and_reshape_summary(total_total_model,      "Total Sum"),
  extract_and_reshape_summary(total_amphibians_model, "Total Amphibians"),
  extract_and_reshape_summary(total_reptiles_model,   "Total Reptiles"),
  extract_and_reshape_summary(total_mammals_model,    "Total Mammals"),
  extract_and_reshape_summary(total_birds_model,      "Total Birds")
)
all_reshaped_results <- rbind(reshaped_results, reshaped_results_total)
print(all_reshaped_results)

# Predicted-effect plots for total richness (green/purple colour scheme)
predicted_amphibians <- ggpredict(total_amphibians_model, terms = "Treatment", ci.lvl = 0.95)
predicted_reptiles   <- ggpredict(total_reptiles_model,   terms = "Treatment", ci.lvl = 0.95)
predicted_mammals    <- ggpredict(total_mammals_model,    terms = "Treatment", ci.lvl = 0.95)
predicted_birds      <- ggpredict(total_birds_model,      terms = "Treatment", ci.lvl = 0.95)
predicted_total      <- ggpredict(total_total_model,      terms = "Treatment", ci.lvl = 0.95)

plot_amphibians_all <- make_predicted_plot(predicted_amphibians, "Amphibian", "#009E73", "#800080") + theme(axis.title.x = element_blank())
plot_reptiles_all   <- make_predicted_plot(predicted_reptiles,   "Reptile",   "#009E73", "#800080") + theme(axis.title.x = element_blank())
plot_mammals_all    <- make_predicted_plot(predicted_mammals,    "Mammal",    "#009E73", "#800080")
plot_birds_all      <- make_predicted_plot(predicted_birds,      "Bird",      "#009E73", "#800080")
plot_total_all      <- make_predicted_plot(predicted_total,      "Total",     "#009E73", "#800080")

extract_coefficients_total <- function(model, model_name) {
  df <- as.data.frame(summary(model)$coeff$cond)
  df$term      <- rownames(df)
  df$conf.low  <- df$Estimate - 1.96 * df$`Std. Error`
  df$conf.high <- df$Estimate + 1.96 * df$`Std. Error`
  df$color     <- "black"
  df$color[df$term %in% c("(Intercept)", "Intercept (Outside)")] <- "gray"
  sel <- df$term == "Treatmentmines_10km" & df$`Pr(>|z|)` < 0.05
  df$color[sel & df$Estimate > 0] <- "#009E73"
  df$color[sel & df$Estimate < 0] <- "#FF7F7F"
  df$model <- model_name
  df
}

new_all_coef_total <- bind_rows(
  extract_coefficients_total(total_amphibians_model, "Amphibian"),
  extract_coefficients_total(total_reptiles_model,   "Reptile"),
  extract_coefficients_total(total_mammals_model,    "Mammal"),
  extract_coefficients_total(total_birds_model,      "Bird"),
  extract_coefficients_total(total_total_model,      "Total")
)
new_all_coef_total$model <- factor(new_all_coef_total$model,
                                    levels = c("Total", "Amphibian", "Reptile", "Mammal", "Bird"))

pseudo_plot_total <- plot_pseudo_coefficient(new_all_coef_total)

patchwork::wrap_plots(
  plot_total_all,      plot_amphibians_all,
  plot_reptiles_all,   plot_mammals_all,
  plot_birds_all,      pseudo_plot_total,
  ncol = 2, widths = c(1, 1), heights = c(1, 1, 1)
)


# ==============================================================================
# SCRIPT 5: GLMs — probability of mining being listed as a threat
# ==============================================================================
#
# Fits binomial GLMs modelling the probability that a species has mining
# listed as a threat on the IUCN Red List, as a function of the proportion
# of its range overlapping with mines. Models are fitted separately for each
# taxonomic group (mammals, birds, reptiles, amphibians, invertebrates,
# plants, fungi).
#
# Input:  Data_S1.xlsx (supplementary data file included in this repository)
# Output: glm_prob_by_group_excluding_all_species.png

library(dplyr)
library(ggplot2)
library(readxl)

# Read from supplementary data file and align column names
raw <- read_excel("Data_S1.xlsx")
data <- raw %>%
  rename(
    group1                  = Group,
    prop                    = `Existing overlap proportion`,
    mining_listed_as_threat = `Mining threat status (IUCN)`
  ) %>%
  mutate(
    group1 = case_when(
      group1 %in% c("arthropod", "mollusca", "other") ~ "inverts",
      group1 == "plant" ~ "plants",
      TRUE ~ group1
    )
  )

data$group1 <- as.character(data$group1)
data$mining_listed_as_threat <- factor(data$mining_listed_as_threat, levels = c("0", "1"))

# Exclude the "all species" aggregate group and rename to display labels
plot_data <- data %>%
  filter(!tolower(group1) %in% c("all species", "all_species", "all")) %>%
  mutate(group1 = recode(group1,
                         bird      = "Birds",
                         amphibian = "Amphibians",
                         reptile   = "Reptiles",
                         mammal    = "Mammals",
                         inverts   = "Invertebrates",
                         plants    = "Plants",
                         fungi     = "Fungi"
  )) %>%
  mutate(group1 = factor(group1,
                         levels = c("Fungi", "Plants", "Invertebrates",
                                    "Amphibians", "Reptiles", "Birds", "Mammals")))

# Fit a binomial GLM per group and generate predicted probabilities
make_pred_df <- function(df, grp, n = 200) {
  m    <- glm(mining_listed_as_threat ~ prop, family = binomial, data = df)
  newd <- data.frame(prop = seq(0, 1, length.out = n))
  pr   <- predict(m, newdata = newd, type = "link", se.fit = TRUE)
  newd$fit   <- plogis(pr$fit)
  newd$lwr   <- plogis(pr$fit - 1.96 * pr$se.fit)
  newd$upr   <- plogis(pr$fit + 1.96 * pr$se.fit)
  newd$group1 <- grp
  df$group1   <- grp
  list(pred = newd, rug = df)
}

groups    <- levels(plot_data$group1)
pred_list <- lapply(groups, function(g) make_pred_df(filter(plot_data, group1 == g), g))
pred_df   <- bind_rows(lapply(pred_list, `[[`, "pred"))
rug_df    <- bind_rows(lapply(pred_list, `[[`, "rug"))

p <- ggplot() +
  geom_ribbon(data = pred_df, aes(x = prop, ymin = lwr, ymax = upr, fill = group1),
              alpha = 0.25, colour = NA) +
  geom_line(data = pred_df, aes(x = prop, y = fit, colour = group1), linewidth = 1) +
  geom_rug(data = rug_df,
           aes(x = prop, y = as.numeric(as.character(mining_listed_as_threat))),
           sides = "tb", alpha = 0.5) +
  facet_wrap(~ group1, ncol = 3) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Proportion of range impacted by mining",
       y = "Probability of mining being listed as a threat") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90"),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 2 * 11),
        axis.title.y = element_text(size = 2 * 11))

ggsave("glm_prob_by_group_excluding_all_species.png", p, width = 9, height = 8, dpi = 300)
p
