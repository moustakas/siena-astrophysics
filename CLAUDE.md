# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Purpose

This is a research and teaching repository for undergraduate students at Siena College doing astrophysics research with Prof. John Moustakas. It contains Python scripts, Jupyter notebooks, configuration files, and supporting materials for projects spanning galaxy clusters, large-scale structure, DESI, and N-body simulations.

## Environment Setup

This repo targets two primary machines: the Siena lab (Ubuntu workstations) and NERSC/nyx (HPC).

**DESI environment (lab machines):**
```bash
source bin/desi-env-lab
# or, if the alias is configured:
desi
```

**legacyhalos environment (nyx):**
```bash
source bin/legacyhalos-env-nyx
```

**Updating DESI packages from GitHub:**
```bash
bin/desi-update-packages   # requires $DESI_PRODUCT_ROOT to be set
```

Dot-files for bash profiles are in `etc/` (`dotbash_profile-lab`, `dotbash_profile-nyx`, `dotbashrc-hpcc`) and are meant to be symlinked or sourced from `~/.bash_profile`. The matplotlib config `etc/matplotlibrc` uses `TkAgg` backend by default.

## Running Code

Most work happens in Jupyter notebooks (`*.ipynb`). Launch via:
```bash
jupyter notebook
# or on NERSC:
# https://jupyter.nersc.gov
```

Standalone Python scripts are run directly:
```bash
python research/lss/cute-2ptcor.py --corrtype monopole --nrandom 1 --docute --qaplots
```

There is no project-wide test suite or build system.

## Repository Structure

```
bin/        Environment activation scripts for DESI and legacyhalos
etc/        Dot-files (bash profiles, matplotlibrc, emacs, jupyterconfig)
pro/        Legacy IDL/GDL scripts
research/   Student and faculty research projects (see below)
teaching/   Course and fellows teaching materials
```

**Active research directories:**

- `research/lss/` — Large-scale structure / BAO analysis using SDSS-III BOSS DR11. Uses the external `CUTE` C code (compiled separately, put in `$PATH`) for 2-point correlation functions. Requires `$LSS_CUTE` env var pointing to data directory.
- `research/csmf_redmapper/` — Conditional stellar mass function of brightest cluster galaxies (BCGs) from the redMaPPer catalog. Notebooks use `pyfits`, `esutil`, `healpy`.
- `research/redmapper/` — Stellar mass comparison between redMaPPer clusters and Prospector SED fits.
- `research/hff/` — Hubble Frontier Fields high-redshift galaxy SEDs and photometric catalogs.
- `research/desi/` — DESI target sky maps and related notebooks.
- `research/gadget2/` — Gadget-2 N-body simulation visualization via Blender.
- `research/advlab/desilens/` — Gravitational lensing simulation notebooks.
- `research/summer2024/` — Current summer research onboarding guide and resources.
- `research/summer-archive/` — Archived summer student work by year (2013–2023).

## Key External Dependencies

- **DESI software stack**: `desiutil`, `desimodel`, `desispec`, `desisim`, `desitarget`, `specter`, `redrock` — all from github.com/desihub, managed via `bin/desi-update-packages`
- **CUTE**: External C code for 2-point correlation functions (compiled separately)
- **Standard astro stack**: `numpy`, `matplotlib`, `astropy`, `healpy`, `scipy`
- **SED fitting**: `prospector`, `fsps` (for redmapper work)
- **Legacy**: `pyfits` (now `astropy.io.fits`), `esutil`

## Data Conventions

- FITS files are gitignored — data lives on NERSC (`/global/work/`) or nyx (`$IM_WORK_DIR`)
- Key env vars: `$DESI_ROOT`, `$DESI_PRODUCT_ROOT`, `$LEGACYHALOS_DIR`, `$LSS_CUTE`, `$DUST_DIR`, `$IM_WORK_DIR`
- Catalogs in `research/` (`.txt`, `.cat`, `.ecsv`) are small enough to commit; binary FITS data is not

## Computing Resources

- **Lab**: Dell Precision workstations in Roger Bacon Hall 113, Ubuntu, conda at `/home/desi/anaconda3`
- **nyx**: `nyx.siena.edu`, Anaconda at `/usr/local/anaconda3`
- **NERSC**: `https://jupyter.nersc.gov` for notebooks; DESI project directories at `/global/project/projectdirs/desi/`
- **Siena HPCC**: `https://www.siena.edu/centers-institutes/high-performance-computing-center/`
