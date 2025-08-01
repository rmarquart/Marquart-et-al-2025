# Coupled Sea Ice Dynamics and Thermodynamics Simulations

This repository contains all necessary data, scripts, simulation, and solver files for the study presented in **Marquart et al. (2025)**. The study focuses on coupled (and decoupled) sea ice dynamics and thermodynamics simulations using [**OpenFOAM**](https://www.openfoam.com/) and **Python**.

---

## 📁 Folder Structure

### `ERA5_datasets`
Contains input variables for the thermodynamic model, sourced from ERA5:
- `specific_humidity_air.grib`
- `specific_humidity_surface.grib`
- `remaining_variables.grib`

These variables are also freely available at the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download).

---

### `Job_submission_scripts`
Scripts for submitting jobs to an HPC cluster.

#### Full-field Case (Wavelengths: 120m, 240m, 360m)
- **Dynamics-only model**: north, east, south, west directions
- **Thermodynamics-only model**: west direction
- **Fully coupled model**: west direction

#### Test Case (Wavelengths: 120m, 240m, 360m)
- **Dynamics-only model**: horizontal, vertical, diagonal, zigzag directions

> ⚠️ **Note**: Update file paths to match your directory structure before submitting.

---

### `Simulation_files`
Organized by simulation type and region.

#### Full-field Case
- **Dynamics**
  - `East`, `North`, `South`, `West`
    - `T120`, `T240`, `T360`: OpenFOAM simulation folders

- **Dynamics & Thermodynamics**
  - `West`
    - `Python_real`
      - `Previous index values`: generated index files for ice floes and grease ice
      - `Saved variables`: `.pkl` files for thermodynamics model input
      - `Storage alpha/csv/eta/h`: storage folders
      - `T120`, `T240`, `T360`: OpenFOAM model folders
      - `T120_cont`, `T240_cont`, `T360_cont`: thermodynamics model folders
    - `Python_storage_real`
      - `T120`, `T240`, `T360`: thickness file generation and updates
      - `T120_storage_real`, `T240_storage_real`, `T360_storage_real`: output storage
      - `Thermodynamics`: same structure as Dynamics & Thermodynamics  
        > Wave amplitude in `constant/transportProperties` is set to zero.  
        > Application in `system/controlDict` is set to `my_seaIce_06112024_v2`.

> ⚠️ **Note**: Update file paths to match your directory structure before submitting.

#### Test Case
- **Dynamics**
  - `Diagonal`, `Horizontal`, `Vertical`, `Zigzag`
    - `T120`, `T240`, `T360`: OpenFOAM simulation folders

---

### `Solver_files`
Contains OpenFOAM solvers used in the simulations:
- `my_seaIce_06112024`: for dynamics-only and fully coupled models
- `my_seaIce_06112024_v2`: for thermodynamics-only model

---

## 📌 Citation
If you use this repository, please cite:

Marquart, R., Alberello, A., Bogaers, A., De Santi, F., Vichi, M., 2025. WIce-FOAM 1.0: Coupled dynamic and thermodynamic modelling of heterogeneous sea ice and waves using OpenFOAM-v2306. EGUsphere 1–35. https://doi.org/10.5194/egusphere-2025-2184

---

## 📬 Contact
For questions or collaboration inquiries, please contact the corresponding author of the study.
