Moustakas Summer 2026 Research
==============================

Overview
--------

This document assembles a few notes and resources to help you get started with
summer research. Although research is very rewarding, it can also frequently be
extremely frustrating! My best advice is for you to try to keep your eye on the
*big picture* of what you're trying to accomplish and know that the struggle is
oftentimes just as important as the end result. Indeed, we often learn best by
struggling, and sometimes we make our greatest insights and discoveries when
things appear to be most bleak!

You should also not hesitate to help out one another. Research is
almost always a *collaborative* process; for example, I strongly
encourage you to [pair
code](https://stackify.com/pair-programming-advantages). We will also
use [Slack](https://slack.com) and, occasionally,
[Zoom](https://zoom.us) to stay in frequent contact. Like tackling a
hard quantum physics problem with a group of classmates, research is
best done in small groups and teams who communicate often.

Finally, at its heart, modern astrophysics is an applied data science:
Astronomers frequently mine large datasets---using one or more modern programming
languages like [Python](https://python.org)---in order to gain scientifically
interesting insight from those data. Throughout this journey, a good
understanding of statistics is crucial, as well as the ability to clearly
communicate your findings to your peers, both in writing and orally.

What do I need to get started?
-------------------------------

We will be using many of the same collaboration tools used by astronomers, as
well as some more specialized and proprietary access to data and computing
resources specific to astro research at Siena University. Therefore, before the
start of summer research, please be sure you have (or have completed) the
following:

* A [Github](https://github.com) account. Once you have your handle,
  please share it with Moustakas on Slack so you can be given access
  to the
  [siena-astrophysics](https://github.com/moustakas/siena-astrophysics)
  repository (and other relevant repositories).

* A completed [DESI Membership
  Form](https://desi.lbl.gov/desipub/app/MembershipForm/form) (get the
  login credentials from Moustakas). Once your membership is approved,
  Moustakas will share resources and getting-started guides for
  navigating the DESI collaboration tools and data.

* A [NERSC](https://www.nersc.gov/) account (contact Moustakas for
  instructions). NERSC, the National Energy Research Scientific
  Computing Center, is one of the world's most powerful supercomputing
  centers, and you will be using it frequently for your research,
  especially the [NERSC Jupyter server](https://jupyter.nersc.gov).

* A [Slack](https://slack.com) account so you can be added to the
  [DESI Slack workspace](https://desisurvey.slack.com).

* An [Overleaf](https://overleaf.com) account. We will use the
  typesetting platform [LaTeX](https://www.latex-project.org/) to
  write a summary research paper this summer.

Finally, I recommend a good lab book (or you might consider keeping an
electronic journal or set of notes, ideally on Github) and enthusiasm!

Responsible Conduct of Research
--------------------------------

Whether you are being supported this summer by CURCA or on a grant
(from either the National Science Foundation or the Department of
Energy), you must complete an online [Responsible Conduct of
Research](https://about.citiprogram.org/en/homepage) training program,
in which you will learn about collaborative research, responsible
authorship, data management, etc.

*Note: If you haven't completed the training in the previous three
years, it must be completed during the first week of your summer
research experience and no later than June 1.*

Computers and more!
--------------------

The astrophysics research lab, Roger Bacon Hall 113, has five
computers which you are free to use (ask Moustakas for the punch code
to the room and the login to the computers). Each machine is a Dell
Precision 3450 with 8 Intel i7 processors, 16 GB of RAM, and 500 GB of
SSD disk space, running [Ubuntu 24.04](https://ubuntu.com/). In
general, you will be using the Ubuntu/Linux operating system. Linux,
comparable to Windows and OSX, is the standard scientific computing
environment in astronomy.

In addition to these computers, during your summer research experience
you may utilize NERSC (particularly the [NERSC Jupyter
server](https://jupyter.nersc.gov)), [Google
Colaboratory](https://colab.research.google.com), the server
*nyx.siena.edu*, the [Siena University high-performance computing
cluster
(HPCC)](https://www.siena.edu/centers-institutes/high-performance-computing-center/),
and other computational resources.

Dark Energy Spectroscopic Instrument (DESI)
--------------------------------------------

[DESI](https://desi.lbl.gov) is a state-of-the-art, 5-year spectroscopic
redshift survey and [Stage IV dark energy
experiment](https://arxiv.org/pdf/1604.07626.pdf) which will place unprecedented
constraints on the expansion history of the universe and our fundamental
understanding of dark energy. The survey---and the grant which is paying you!---is
supported by the [Department of Energy (DOE) Office of
Science](https://www.energy.gov/science/office-science) as one of the core
experiments of the [Cosmic Frontier
Vision](https://science.osti.gov/hep/Research/Cosmic-Frontier).

Going forward,
[DESI-II](https://indico.in2p3.fr/event/31064/contributions/130380/attachments/81926/120926/DESI2_Future_Surveys.pdf)
(see also [this
talk](https://indico.physics.lbl.gov/event/2382/contributions/7561/attachments/3765/4996/Dawson_DESI-II_P5TownHall.pdf))
aims to push the cosmological constraints obtainable from
ground-based, highly multiplexed optical spectroscopy even further.

As time permits, I recommend you read the following overview papers
and data release documentation for DESI:

* [DESI Data Release 1 (DR1) Paper](https://arxiv.org/abs/2503.14745)
* [DESI/DR1 Documentation](https://data.desi.lbl.gov/doc/releases/dr1/)
* [DESI tutorials](https://github.com/desihub/tutorials/tree/main/01_getting_started).

DESI Legacy Imaging Surveys
-----------------------------

The [DESI Legacy Imaging Surveys](https://legacysurvey.org) are a set of
complementary ground-based imaging surveys that provided the photometric
targeting catalog for DESI. The surveys cover roughly 14,000 square degrees of
the extragalactic sky in three optical bands (*g*, *r*, *z*) and are combined
with mid-infrared photometry from the *WISE* satellite. The Legacy Surveys
website provides access to imaging data, photometric catalogs, and a variety of
tools for exploring the data.

Building your toolkit
----------------------

To start building your technical toolkit, you should work through the following
(or comparable) tutorials from [Software
Carpentry](https://software-carpentry.org/lessons) and elsewhere:

* Linux/Unix and Git:
  * [The Unix Shell](https://swcarpentry.github.io/shell-novice)
  * [Version Control with Git](https://swcarpentry.github.io/git-novice) (see also [git - the simple guide](https://rogerdudler.github.io/git-guide))
* Python:
  * [Python for Astronomers](https://prappleizer.github.io/)
  * [Programming with Python](https://swcarpentry.github.io/python-novice-inflammation)
  * [Google's Python Class](https://developers.google.com/edu/python)
  * [PEP 8 - Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/)
  * [Scipy Lecture Notes](https://scipy-lectures.org)
* Astronomical data and tools:
  * [DESI Getting Started Tutorials](https://github.com/desihub/tutorials/tree/main/01_getting_started)
  * [NOIRLab Astro Data Lab](https://datalab.noirlab.edu/) --- a science platform with access to large astronomical survey datasets and Jupyter notebook environments
* Miscellaneous:
  * [Learn LaTeX in 30 minutes](https://www.overleaf.com/learn/latex/Learn_LaTeX_in_30_minutes)
  * [Learning Markdown](https://daringfireball.net/projects/markdown/syntax)
  * [Installing Anaconda](https://docs.anaconda.com/anaconda/install) and [Getting Started with Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html)

Digging Deeper: Data Science and Machine Learning
--------------------------------------------------

Artificial intelligence (AI), machine learning (ML), and deep learning (DL) have
taken science---and our society---by storm the past couple of decades. You will
be gaining familiarity with the techniques underlying these enormously powerful
methods and how they are used to solve a variety of problems in astrophysics and
for DESI in particular. Additional resources and tutorials will be added here.
