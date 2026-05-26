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

### Task 2 --- Self-Supervised Representation Learning (Weeks 2--4)

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

### Task 3 --- Zoobot Fine-Tuning (Weeks 3--5)

Fine-tune a pretrained Zoobot model on the labeled training set from Task 1
to classify SGA-2025 galaxies by one or more morphological properties of
scientific interest.

* Select one or two specific classification questions well-suited to the
  SGA-2025 sample (e.g., disk vs. spheroid, presence of a bar, merger
  signatures).
* Fine-tune the Zoobot model using the labeled training set, following the
  Zoobot tutorials for small-data fine-tuning.
* Evaluate model performance on a held-out validation set. Assess where the
  model is confident vs. uncertain.
* Run the trained classifier on the full subset of SGA-2025 galaxies with
  available cutouts and examine the distribution of predicted morphologies.
* Compare Zoobot classifications against visual inspection labels and against
  classifications from the ssl-legacysurvey approach (Task 2).

Rough Timeline
--------------

**Week 1 (May 26)**: Read required papers. Explore the SGA-2020 web app and
Legacy Survey cutout server. Begin visual inspection and image labeling. Set
up NERSC environment.

**Week 2 (June 1)**: Complete training set construction. Begin
ssl-legacysurvey feature extraction. Prepare for DESI observing run (read
Support Observer instructions).

**Week 3 (June 8)**: DESI observing shift (June 9--11). Continue feature
extraction; run dimensionality reduction and similarity searches. Begin Zoobot
environment setup and test runs.

**Week 4 (June 15)**: Zoobot fine-tuning runs. Evaluate model performance.
Compare ssl-legacysurvey and Zoobot outputs.

**Week 5 (June 22)**: Run classifiers on full SGA-2025 subset. Generate
interactive visualizations and summary figures.

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
