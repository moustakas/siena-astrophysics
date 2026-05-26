Project 3: DESI Steel Program Spectroscopic Validation
=======================================================

Overview
--------

[Steel](https://arxiv.org/abs/2603.10113) is a pilot program for DESI-II, the
planned next-generation extension of DESI, designed to calibrate the redshift
distributions of faint galaxies used as sources in upcoming weak lensing
surveys (LSST, Euclid, and Roman). Weak lensing uses the subtle distortion of
background galaxy shapes by foreground matter to map the large-scale
distribution of mass, but it requires precise knowledge of the redshifts of
those background galaxies. Steel addresses this by obtaining DESI spectra for
a carefully selected sample of faint galaxies (22 < *i* < 24), providing the
ground truth needed to calibrate photometric redshift estimates from imaging
alone. Our group is contributing to this collaboration on two fronts: visual
inspection (VI) of Steel spectra to confirm or reject automated redshift fits,
and using [FastSpecFit](https://github.com/desihub/fastspecfit) to help reduce
the VI workload by flagging objects with high-confidence measurements.

Background Reading
------------------

* DeRose et al. (2026) --- [Steeling Weak Lensing Source Galaxy Samples
  against Systematics using Wide Field
  Spectroscopy](https://arxiv.org/abs/2603.10113). The forecast paper
  motivating the Steel program. Pay particular attention to the science
  requirements on σ(⟨z⟩) and the role of representative spectroscopy.

* Ratajczak et al. (2026) --- [The Compilation and Validation of the
  Spectroscopic Redshift Catalogs for the DESI-COSMOS and DESI-XMMLSS
  Fields](https://arxiv.org/abs/2508.09286). Describes the VI-based
  validation pipeline for DESI deep fields, including the use of FastSpecFit
  metrics as part of the redshift quality classification framework. The
  COSMOS and XMM-LSS fields are also Steel fields.

* [FastSpecFit documentation and
  repository](https://github.com/desihub/fastspecfit) --- the spectral
  synthesis and emission-line fitting tool used in Task 2.
* Steel VI training slides --- to be shared by the collaboration before
  VI begins.

Key Resources
-------------

* [DESI Trac/Wiki](https://desi.lbl.gov/trac/wiki) --- data access,
  software documentation, and VI protocols.
* [FastSpecFit repository](https://github.com/desihub/fastspecfit) ---
  source code and documentation for the spectral fitting tool.
* [NOIRLab Astro Data Lab](https://datalab.noirlab.edu/) --- provides
  notebook environments and access to supporting photometric catalogs.
* NERSC --- reduced Steel spectra and VI output files are stored on NERSC
  (contact Moustakas for path and access instructions).

Tasks
-----

### Task 1 --- Visual Inspection of Steel Spectra (All Students, Weeks 1--4)

VI is the critical bottleneck for the Steel program. The automated redshift
fitter (Redrock) produces a figure of merit Δχ², but this statistic alone is
insufficient for confident redshifts: roughly 50% of Steel objects fall in the
range 2 < Δχ² < 100, where human inspection is essential to distinguish
secure redshifts from failures.

All students will be trained in the PROSPECT VI tool using training materials
provided by the collaboration. VI is organized into batches of 50 spectra,
sorted in order of decreasing Δχ² (so batches with larger page IDs contain the
most ambiguous spectra and are highest priority). A tracking spreadsheet on
NERSC records progress by field and reviewer.

**Authorship note**: Any collaborator who completes 60 or more VI batches ---
with at least 30 in the high-priority 2 < Δχ² < 100 range --- earns
co-authorship on the Steel data paper. This is a concrete and achievable goal
for a motivated student over 6 weeks.

Practical steps:

* Obtain NERSC access and locate the VI output directory for the reduced
  fields (0051, 0052; Moustakas will provide paths).
* Complete the VI training slides before inspecting any real data.
* Work through batches systematically, recording results in the tracking
  spreadsheet. Prioritize 2 < Δχ² < 100 batches.
* Flag unusual or interesting spectra (e.g., emission-line galaxies at
  unexpected redshifts, apparent blends, strong absorption systems) for
  group discussion.

### Task 2 --- FastSpecFit-Based VI Prioritization (Weeks 2--5)

The collaboration is actively exploring whether FastSpecFit outputs can be
used to predict which objects require VI and which can be auto-accepted,
substantially reducing the human VI burden. This is an open problem where our
group can make a direct contribution.

FastSpecFit fits the DESI spectrum and broadband photometry simultaneously,
producing emission-line fluxes, equivalent widths, continuum properties, and
associated quality metrics. For objects with secure emission-line detections
(e.g., strong [OII], Hα, or [OIII]), the FastSpecFit redshift and quality
flags may be sufficient to bypass VI entirely.

Practical steps:

* Install and configure FastSpecFit on NERSC following the repository
  documentation.
* Run FastSpecFit on the reduced Steel fields (0051, 0052) to generate
  emission-line catalogs and quality metrics.
* Cross-match FastSpecFit outputs against the completed VI results from Task 1
  to assess how well FastSpecFit quality flags predict VI outcomes.
* Develop a simple classification scheme (e.g., a threshold on emission-line
  signal-to-noise or FastSpecFit χ²) that identifies objects which can be
  safely auto-accepted.
* Quantify the reduction in required VI and the contamination rate (fraction
  of incorrectly auto-accepted failures).

### Task 3 --- Redshift Success Rate Characterization (Weeks 4--6)

After VI is substantially complete on the reduced fields, characterize the
Steel spectroscopic success rate and assess whether the program meets its
design requirements.

Key questions (from the collaboration's "next steps"):

* How does the success rate vary with *i*-band magnitude, color, and
  integration time?
* Does the *i* < 23.5 sample reach the required density and success rate
  without pushing to *i* < 24.0?
* How does the Δχ² distribution (and the FastSpecFit-predicted success rate)
  compare between fields?

Practical steps:

* Compile VI results for fields 0051 and 0052 into a summary catalog.
* Compute success rate as a function of magnitude (comparing *i* < 23.5 vs.
  *i* < 24.0), Δχ², and exposure time.
* Plot the redshift distribution *n(z)* for the confirmed Steel subsample and
  compare with photometric predictions.
* Produce a summary figure set suitable for sharing with the collaboration
  and for inclusion in the eventual data paper.

Rough Timeline
--------------

**Week 1 (May 26)**: Read required papers and Steel VI training slides. Set
up NERSC access and locate Steel data. Complete VI training. Begin VI on
high-Δχ² batches (easy end first).

**Week 2 (June 1)**: Continue VI. Begin FastSpecFit installation and test
runs on a small subset. Prepare for DESI observing run (read Support Observer
instructions).

**Week 3 (June 8)**: DESI observing shift (June 9--11). Continue VI, shifting
to 2 < Δχ² < 100 priority batches. FastSpecFit production runs on fields 0051
and 0052.

**Week 4 (June 15)**: Cross-match FastSpecFit outputs with VI results.
Develop and evaluate the VI prioritization classifier.

**Week 5 (June 22)**: Success rate characterization. Compare *i* < 23.5 vs.
*i* < 24.0. Generate *n(z)* plots and summary figures.

**Week 6 (June 29)**: Polish deliverables. Write symposium abstract. Prepare
poster. Share FastSpecFit results with collaboration.

Deliverables
------------

* Completed VI for as many Steel batches as possible, recorded in the
  collaboration tracking spreadsheet (target: ≥60 batches per student,
  with ≥30 in the 2 < Δχ² < 100 range for authorship eligibility).
* A Jupyter notebook documenting the FastSpecFit runs on Steel fields 0051
  and 0052, including the VI prioritization analysis and classifier evaluation.
* A summary figure set characterizing the Steel redshift success rate as a
  function of magnitude and Δχ², suitable for the collaboration and the
  data paper.
* A research poster presented at the Siena University fall symposium.
