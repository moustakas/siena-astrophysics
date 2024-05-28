Moustakas & Finn Summer 2024 Research
=====================================

Overview
--------

This document assembles a few notes and resources to help you get started with
summer research. Although research is very rewarding, it can also frequently be
extremely frustrating!  Our best advice is for you to try to keep your eye on the
*big picture* of what you're trying to accomplish and know that the struggle is
oftentimes just as important as the end result. Indeed, we often learn best by
struggling, and sometimes we make our greatest insights and discoveries when
things appear to be most bleak!

You should also not hesitate to help out one another. Research is almost always
a *collaborative* process; for example, we strongly encourage you to [pair
code](https://stackify.com/pair-programming-advantages).  We will also use
[Slack](https://slack.com) and, occasionally, [Zoom](https://zoom.us) to stay in
frequent contact. Like tackling a hard quantum physics problem with a group of
classmates, research is best done in small groups and teams who communicate
often.

Finally, at its heart, modern astrophysics is an applied data science:
Astronomers frequently mine large datasets--using one or more modern programming
languages like [Python](https://python.org)--in order to gain scientifically
interesting insight from those data. Throughout this journey, a good
understanding of statistics is crucial, as well as the ability to clearly
communicate your findings to your peers, both in writing and orally.

What do I need to get started?
------------------------------

We will be using many of the same collaboration tools used by astronomers, as
well as some more specialized and proprietary access to data and computing
resources specific to astro research at Siena. Therefore, before the start of
summer research, please be sure you have (or have completed) the following:

* A [Slack](https://slack.com) account so you can be added to the
  [siena-astrophysics Slack workspace](https://siena-astrophysics.slack.com).

* A [Github](https://github.com) account. Once you have your handle, please
  share it with Moustakas on Slack so you can be given access to the
  [siena-astrophysics](https://github.com/moustakas/siena-astrophysics)
  repository (and other relevant repositories).

* An [Overleaf](https://overleaf.com) account. We will use the typesetting
  platform [LaTeX](https://www.latex-project.org/) to write a summary research
  paper this summer.

* A completed [DESI membership
  form](https://desi.lbl.gov/desipub/app/MembershipForm/form) (get the login
  credentials from Moustakas). Once your application has been approved, you will
  receive instructions for how to access the [DESI Trac
  system](https://desi.lbl.gov/trac/wiki), which has a *ton* of useful and
  important information (see especially the [Getting
  Started](https://desi.lbl.gov/trac/wiki/GettingStarted) page).

* A [NERSC
  account](https://desi.lbl.gov/trac/wiki/Computing/AccessNersc). [NERSC](https://www.nersc.gov/),
  the National Energy Research Scientific Computing Center, is one of the
  world's most powerful supercomputing centers, and you will be using it
  frequently for your research, especially the [NERSC Jupyter
  server](https://jupyter.nersc.gov).

Finally, I recommend a good lab book (or you might consider keeping an
electronic journal or set of notes, ideally on Github) and enthusiasm!

Responsible Conduct of Research
-------------------------------

Whether you are being supported this summer by CURCA or on a grant (from either
the National Science Foundation or the Department of Energy), you must complete
an online [Responsible Conduct of
Research](https://about.citiprogram.org/en/homepage) training program (see
[these
instructions](https://www.siena.edu/files/resources/responsible-conduct-of-research-2016.pdf)),
in which you will learn about collaborative research, responsible authorship,
data management, etc.  

*Note: If you haven't completed the training in the previous three years, it
must be completed during the first week of your summer research experience.*

Computers and more!
-------------------

The astrophysics research lab, Roger Bacon Hall 113, has three computers which
you are free to use (ask either Moustakas or Finn for the punch code to the room
and the login to the computers). Each machine is a Dell Precision 3450 with 8
Intel i7 processors, 16 GB of RAM, 500 GB of SSD disk space, and a dual-boot
Windows 10 and [Linux/Ubuntu](https://ubuntu.com/) operating system
installed. In general, you will be using the Linux/Ubuntu operating system,
which is the standard scientific computing environment in astronomy.

In addition to these computers, during your summer research experience you may
utilize NERSC, particularly the [NERSC Jupyter
server](https://jupyter.nersc.gov), [Google
Colaboratory](https://colab.research.google.com), the server *nyx.siena.edu*,
the [Siena College high-performance computing cluster
(HPCC)](https://www.siena.edu/centers-institutes/high-performance-computing-center/),
and other computational tools.

Dark Energy Spectroscopic Instrument (DESI)
-------------------------------------------

[DESI](https://desi.lbl.gov) is a state-of-the-art, 5-year spectroscopic
redshift survey and [Stage IV dark energy
experiment](https://arxiv.org/pdf/1604.07626.pdf) which will place unprecedented
constraints on the expansion history of the universe and our fundamental
understanding of dark energy. The survey--and the grant which is paying you!--is
supported by the [Department of Energy (DOE) Office of
Science](https://www.energy.gov/science/office-science) as one of the core
experiments of the [Cosmic Frontier
Vision](https://science.osti.gov/hep/Research/Cosmic-Frontier). 

As time permits, I recommend you read the following overview papers on the
instrumentation and science case for DESI:
* [Overview of the Instrumentation for the Dark Energy Spectroscopic Instrument
](https://arxiv.org/abs/2205.10939)
* [The DESI Experiment Part I: Science,Targeting, and Survey Design](https://arxiv.org/abs/1611.00036) 

Building your toolkit
---------------------

To start building your technical toolkit, you should work through the following
(or comparable) tutorials from [Software
Carpentry](https://software-carpentry.org/lessons) and elsewhere:

* Linux/Unix, EMACS, and Git:
  * [The Unix Shell](http://swcarpentry.github.io/shell-novice)
  * [Absolute Beginner's Guide to EMACS](http://www.jesshamrick.com/2012/09/10/absolute-beginners-guide-to-emacs) (and you might as well learn a little bit of [vim](https://linuxconfig.org/vim-tutorial), too)
  * [Version Control with Git](http://swcarpentry.github.io/git-novice) (see also [Git Immersion](http://gitimmersion.com/) and [git - the simple guide](https://rogerdudler.github.io/git-guide)) 
* Python:
  * [Python for Astronomers](https://prappleizer.github.io/)
  * [Programming with Python](http://swcarpentry.github.io/python-novice-inflammation)
  * [Google's Python Class](https://developers.google.com/edu/python)
  * [PEP 8 - Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/)
  * [Scipy Lecture Notes](http://www.scipy-lectures.org/index.html)
* Miscellaneous:
  * [Learn LaTeX in 30 minutes](https://www.overleaf.com/learn/latex/Learn_LaTeX_in_30_minutes) 
  * [Learning Markdown](https://daringfireball.net/projects/markdown/syntax)
  * [Installing Anaconda](https://docs.anaconda.com/anaconda/install) and [Getting Started with Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html)

Digging Deeper: Data Science and Machine Learning
-------------------------------------------------

Artificial intelligence (AI), machine learning (ML), and deep learning (DL) have
taken science--and our society--by storm the past couple decades. (Here's a nice
recent article which talks about the enormous impact that machine learning is
having on physics and other lines of scientific inquiry: [Powerful "Machine
Scientists" Distill the Laws of Physics From Raw
Data](https://www.quantamagazine.org/machine-scientists-distill-the-laws-of-physics-from-raw-data-20220510/).)

You will be gaining familiarity with the techniques underlying these enormously
powerful methods and how they are used to solve a variety of problems in
astrophysics and for DESI in particular. With these ideas in mind, check out and
work through the following tutorials:
* [Python Data Science Handbook](https://github.com/jakevdp/PythonDataScienceHandbook)
* [Statistics, Data Mining, and Machine Learning in Astronomy](http://www.astroml.org/index.html)
* [TensorFlow Tutorials](https://www.tensorflow.org/tutorials)
* [scikit-learn: machine learning in Python](https://scipy-lectures.org/packages/scikit-learn/index.html#introducing-the-scikit-learn-estimator-object) 
