Project 2: SGA-2025 Morphology Classification with Machine Learning
====================================================================

Overview
--------

The [Siena Galaxy Atlas (SGA-2025)](https://github.com/moustakas/SGA)
contains multiband images and photometry for ~0.5 million nearby
galaxies drawn from the [DESI Legacy Imaging
Surveys](https://legacysurvey.org). This project uses state-of-the-art
public machine-learning tools to organize and classify SGA-2025
objects by morphology and other visual properties.

Two complementary approaches will be explored. The first applies a
self-supervised representation learned directly from Legacy Survey
galaxy images, which can be used for similarity searches and
unsupervised clustering without any labeled data. The second uses
Zoobot, a deep-learning framework designed for supervised galaxy
morphology classification that can be fine-tuned to specific science
questions with relatively small training sets. Together these tools
can reveal structure in the SGA-2025 that is difficult to find through
catalog-based methods alone.

Background Reading
------------------

* Moustakas et al. (2023) --- [The Siena Galaxy Atlas
  2020](https://ui.adsabs.harvard.edu/abs/2023ApJS..269....3M/abstract)
  (*ApJS*, 269, 3). Provides essential context for the SGA dataset and its
  imaging data products.
* Walmsley et al. (2023) --- [Zoobot: Adaptable Deep Learning Models
  for Galaxy
  Morphology](https://ui.adsabs.harvard.edu/abs/2023JOSS....8.5312W/abstract)
  (*JOSS*, 8, 5312). This paper describes the Galaxy Zoo morphology
  classification framework and the pretrained models available for
  fine-tuning.
* Stein et al. --- [Self-Supervised Representation Learning for Astronomical
  Images](https://github.com/georgestein/ssl-legacysurvey). Read the
  repository README and any linked papers for background on the
  self-supervised approach applied to Legacy Survey data.

Key Resources
-------------

* [SGA-2020 web app](https://sga.legacysurvey.org/) --- explore galaxy
  images and morphologies interactively; useful for building visual intuition.
* [Legacy Survey dev sky viewer](https://www.legacysurvey.org/viewer-dev) ---
  an interactive sky viewer with dedicated SGA-2025 overlay layers:
  `sga2025-parent`, `sga2025-south`, and `sga2025-north`. Use these to
  visually explore the sample footprint and individual objects in context.
* [ssl-legacysurvey repository](https://github.com/georgestein/ssl-legacysurvey)
  --- pretrained model checkpoints, similarity search tools, and 20 TB of
  processed Legacy Survey galaxy cutouts (available at NERSC).
* [Zoobot documentation](https://zoobot.readthedocs.io/en/latest/index.html)
  --- tutorials and API reference for fine-tuning pretrained morphology
  classifiers.
* [DESI Legacy Imaging Surveys](https://legacysurvey.org) --- the underlying
  imaging dataset; image cutouts for SGA-2025 objects are accessible through
  the Legacy Survey cutout server.
* [NOIRLab Astro Data Lab](https://datalab.noirlab.edu/) --- provides
  notebook environments and SQL access to Legacy Survey catalogs.

Tasks
-----

### Task 1 --- Visual Inspection and Training Set Construction (Weeks 1--2)

Before running any models, build intuition by looking at the data.

**Shared activity**: All students will participate in structured visual
inspection of SGA-2025 galaxy images. The goals are to (1) develop a shared
vocabulary for describing galaxy morphology, (2) identify interesting
subpopulations (mergers, edge-on disks, low surface brightness galaxies,
etc.), and (3) construct a labeled training set for Zoobot fine-tuning.

* Explore the SGA-2020 web app and the Legacy Survey sky viewer to get
  a feel for the variety of objects in the sample.
* Use a simple annotation tool (e.g., a Jupyter notebook with ipywidgets, or
  a Google Form) to label several hundred galaxy images by morphological type.
* Document the labeling scheme and inter-rater agreement as part of the
  scientific record.

Tasks 2 and 3 are parallel tracks intended to be pursued simultaneously by
different students, both drawing on the labeled dataset from Task 1. Results
from both tracks will be compared and synthesized during weeks 5--6.

### Task 2 (Track A) --- Self-Supervised Representation Learning (Weeks 2--5)

Apply the pretrained ssl-legacysurvey model to SGA-2025 galaxy images to
extract latent representations, and use those representations to explore the
morphological diversity of the sample.

* Download or access the pretrained ResNet50/MoCov2 model checkpoints from
  the ssl-legacysurvey repository.
* Extract feature vectors for a representative subset of SGA-2025 galaxies
  using the pretrained encoder (no fine-tuning required at this stage).
* Apply dimensionality reduction (UMAP or t-SNE) to visualize the
  representation space and identify clusters.
* Perform similarity searches: given a query galaxy of interest (e.g., a
  known merger, a ring galaxy), retrieve the nearest neighbors in the
  representation space and assess whether the retrieved galaxies share
  visually similar properties.
* Produce interactive visualizations of the embedding space, annotated with
  galaxy thumbnails, for inclusion in the final poster.

### Task 3 (Track B) --- Morphology Classification with Zoobot (Weeks 2--5)

Use a pretrained Zoobot model to classify SGA-2025 galaxies by one or more
morphological properties of scientific interest, fine-tuning as needed using
the labeled training set from Task 1.

* Set up the Zoobot environment and work through the provided tutorials to
  understand the available pretrained models and classification framework.
* Select one or two specific classification questions well-suited to the
  SGA-2025 sample (e.g., disk vs. spheroid, presence of a bar, merger
  signatures).
* Apply pretrained Zoobot models to the SGA-2025 sample; fine-tune on the
  Task 1 training set if needed to improve performance on the target classes.
* Evaluate model performance on a held-out validation set and assess where
  the model is confident vs. uncertain.
* Run the classifier on the full subset of SGA-2025 galaxies with available
  cutouts and examine the distribution of predicted morphologies.

Rough Timeline
--------------

**Week 1 (May 26)**: Read required papers. Explore the SGA-2020 web app and
dev sky viewer. Begin visual inspection and image labeling. Set up NERSC
environment.

**Week 2 (June 1)**: Complete training set construction. *Track A*: begin
ssl-legacysurvey feature extraction. *Track B*: set up Zoobot environment and
run initial tests. Prepare for DESI observing run.

**Week 3 (June 8)**: DESI observing shift (June 9--11). *Track A*: run
dimensionality reduction and similarity searches. *Track B*: apply pretrained
Zoobot models to SGA-2025 subset; begin fine-tuning if needed.

**Week 4 (June 15)**: *Track A*: refine embeddings and produce interactive
visualizations. *Track B*: evaluate model performance and run classifier on
full subset.

**Week 5 (June 22)**: Both tracks generate summary figures. Begin comparison
and synthesis of Track A and Track B results.

**Week 6 (June 29)**: Polish deliverables. Write symposium abstract. Prepare
poster. Commit notebooks and results to the repo.

Deliverables
------------

* A labeled morphology training set for SGA-2025 (saved as an ECSV or FITS
  table in the repository), with documented labeling criteria.
* A Jupyter notebook demonstrating the ssl-legacysurvey feature extraction,
  UMAP/t-SNE visualization, and similarity search results.
* A Jupyter notebook documenting the Zoobot fine-tuning experiment, including
  training curves, validation metrics, and example classifications.
* A morphology classification table for the analyzed SGA-2025 subset,
  suitable for inclusion in a future catalog release.
* A research poster presented at the Siena University fall symposium.
