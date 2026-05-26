Project 1: SGA-2025 Catalog Assembly and Validation
====================================================

Overview
--------

The [Siena Galaxy Atlas (SGA)](https://github.com/moustakas/SGA) is a
catalog of large, nearby galaxies built from the DESI Legacy Imaging Surveys.
The first version, SGA-2020, contains ~383,000 galaxies with multicolor images,
ellipse-fit photometry, and derived physical parameters. The 2025 edition
(SGA-2025) expands the sample to ~0.5 million objects and incorporates
improved photometry and new ancillary data.

This project focuses on validating and enriching the SGA-2025 catalog through
three complementary activities: assembling a comprehensive set of redshifts and
metadata, benchmarking a new ellipse-fitting code, and preparing an
MMA-ready version of the catalog.

Background Reading
------------------

*Required:*

* Moustakas et al. (2023) --- [The Siena Galaxy Atlas
  2020](https://ui.adsabs.harvard.edu/abs/2023ApJS..269....3M/abstract)
  (*ApJS*, 269, 3). Read this carefully; it describes the science case,
  methodology, and data products of SGA-2020 in detail.
* [DESI Data Release 1](https://data.desi.lbl.gov/doc/releases/dr1/) ---
  the primary new redshift source for SGA-2025.

*Recommended:*

* Toivonen et al. (2026) --- [SGA-2025 and Multi-Messenger
  Astrophysics](https://ui.adsabs.harvard.edu/abs/2026A%26A...706A.284T/abstract)
  (*A&A*, 706, A284). Describes how a nearby galaxy atlas can be used to
  prioritize targets for gravitational wave follow-up.
* [ISOSTER documentation](https://massiveseaotters.github.io/isoster/) ---
  the new isophotal fitting tool being benchmarked in Task 3.
* [NED (NASA/IPAC Extragalactic Database)](http://ned.ipac.caltech.edu) ---
  the primary resource for galaxy names, redshifts, and distances.

Key Resources
-------------

* [SGA-2020 web app](https://sga.legacysurvey.org/) --- browse the existing
  catalog, images, and photometric measurements interactively.
* [SGA code repository](https://github.com/moustakas/SGA) --- source code
  for catalog production and analysis.
* [DESI Legacy Imaging Surveys](https://legacysurvey.org) --- the underlying
  imaging data for the SGA.
* [DESI Trac/Wiki](https://desi.lbl.gov/trac/wiki) --- collaboration
  documentation, including data access instructions.
* [NOIRLab Astro Data Lab](https://datalab.noirlab.edu/) --- science
  platform with direct access to Legacy Survey and other catalogs via SQL.

Tasks
-----

### Task 1 — Metadata Assembly (Weeks 1--3)

Gather a standard set of ancillary data for all ~0.5M objects in the
SGA-2025. Moustakas will lead the pipeline infrastructure in Week 1.
Students will run, validate, and extend the pipeline in subsequent weeks.

* **Redshifts**: Cross-match the SGA-2025 against DESI/DR1 and other public
  redshift compilations. For each galaxy, identify the best available
  spectroscopic redshift and record its source.

* **Galaxy names**: Use the NED cross-identification tools to gather a
  standard set of names (NGC, IC, UGC, etc.) for each object. Write scripts
  to query the NED API in batch mode and parse the results.

* **Direct distances**: Compile published direct distance measurements
  (Cepheids, tip of the red giant branch, surface brightness fluctuations,
  etc.) from NED and the literature for the subset of SGA-2025 galaxies that
  have them. These are distinct from and more precise than redshift-based
  distances, and are essential for the MMA task.

**Shared activity**: Visual inspection of DESI spectra. All students will
spend time visually reviewing DESI spectra for a subset of SGA-2025 galaxies
to assess redshift quality and identify interesting or problematic cases.
This builds intuition for the data and directly feeds into the redshift
quality-assurance effort.

### Task 2 --- ISOSTER Benchmarking (Weeks 3--5)

Compare ellipse-fit photometry from
[ISOSTER](https://massiveseaotters.github.io/isoster/) against the existing
photutils-based measurements in SGA-2020 for a representative subset of
galaxies spanning a range of sizes, morphologies, and surface brightness
profiles.

* Select a benchmark sample of ~100--500 galaxies from the SGA-2020 with
  existing high-quality photutils measurements.
* Run ISOSTER on the same galaxies and images using comparable fitting
  parameters.
* Compare outputs: ellipse parameters (position angle, ellipticity, effective
  radius), curve-of-growth photometry, and runtime performance.
* Document discrepancies and identify cases where ISOSTER performs better or
  worse than photutils.
* Produce a benchmarking notebook with summary figures suitable for sharing
  with the ISOSTER developer.

### Task 3 --- MMA-Ready Catalog (Weeks 4--6)

Generate a version of the SGA-2025 optimized for Multi-Messenger Astrophysics
(MMA) searches, in collaboration with Moustakas's colleague at the University
of Arizona.

* Starting from the assembled metadata (Task 1), compute or compile the
  quantities needed for MMA prioritization: direct distances where available,
  redshift-based distances otherwise, stellar masses, and star formation rates.
* Follow the methodology described in Toivonen et al. (2026) to assess
  completeness and assign MMA utility scores to catalog objects.
* Produce the MMA catalog in a standard format (FITS and/or ECSV) with full
  provenance documentation.

Rough Timeline
--------------

+----------+-------------------------------------------------------------+
| Week 1   | Read SGA-2020 paper and explore the web app and code repo.  |
| (May 26) | Set up NERSC environment and data access. Moustakas leads   |
|          | metadata pipeline kickoff.                                  |
+----------+-------------------------------------------------------------+
| Week 2   | Begin NED name cross-matching and redshift assembly         |
| (June 1) | pipeline. Prepare for DESI observing run (read Support       |
|          | Observer instructions).                                     |
+----------+-------------------------------------------------------------+
| Week 3   | DESI observing shift (June 9--11). Continue metadata        |
| (June 8) | pipeline; begin visual inspection of DESI spectra. Start    |
|          | ISOSTER installation and test runs.                         |
+----------+-------------------------------------------------------------+
| Week 4   | Metadata pipeline wrap-up. ISOSTER benchmarking runs and    |
| (June 15)| comparisons. Begin MMA catalog assembly.                    |
+----------+-------------------------------------------------------------+
| Week 5   | Finalize ISOSTER benchmarking notebook. MMA catalog         |
| (June 22)| validation and completeness checks. Generate key figures.   |
+----------+-------------------------------------------------------------+
| Week 6   | Polish deliverables. Write symposium abstract. Prepare      |
| (June 29)| poster. Final catalog and notebook deposits in the repo.    |
+----------+-------------------------------------------------------------+

Deliverables
------------

* A set of Python scripts (committed to the SGA repo) for batch NED
  cross-matching, redshift assembly, and direct distance compilation.
* A Jupyter notebook documenting the ISOSTER benchmarking results, including
  summary figures and a performance comparison table.
* The MMA-ready SGA-2025 catalog in FITS/ECSV format with provenance headers.
* A research poster presented at the Siena University fall symposium.
