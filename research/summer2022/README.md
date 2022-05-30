Moustakas Summer 2022 Research Group
====================================

Overview
--------

This document assembles a few notes and resources to help you get started with
summer research.  Although research is very rewarding, it can also frequently be
extremely frustrating!  My best advice is for you to try to keep your eye on the
*big picture* of what you're trying to accomplish and know that the struggle is
oftentimes just as important as the end result.  Indeed, we often learn best by
struggling, and sometimes we make our greatest insights and discoveries when
things appear to be most bleak!

You should also not hesitate to help out one another.  Research is almost always
a *collaborative* process; for example, I strongly encourage you to [pair
code](https://stackify.com/pair-programming-advantages).  We will also use
[Slack](https://slack.com) and, occasionally, [Zoom](https://zoom.us) to stay in
frequent contact.  Like tackling a hard quantum physics problem with a group of
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
resources specific to my research group. Therefore, before the start of summer
research, please be sure you have (or have completed) the following:

1. A [Slack](https://slack.com) account so you can be added to the
   [siena-astrophysics Slack workspace](https://siena-astrophysics.slack.com).

2. A [Github](https://github.com) account. Once you have your handle, please
   share it with me on Slack so I can give you access to the
   [siena-astrophysics](https://github.com/moustakas/siena-astrophysics)
   repository (and other relevant repositories).

3. Fill out the [DESI membership
   form](https://desi.lbl.gov/trac/wiki/NewMembers#TheNewMemberForm) (get the
   login credentials from Moustakas). Once your application has been approved,
   you will receive instructions for how to access the [DESI Trac
   system](https://desi.lbl.gov/trac/wiki), which has a *ton* of useful and
   important information (see especially the [Getting
   Started](https://desi.lbl.gov/trac/wiki/GettingStarted) page).

4. [Get a NERSC
   account](https://desi.lbl.gov/trac/wiki/Computing/AccessNersc). [NERSC](https://www.nersc.gov/),
   the National Energy Research Scientific Computing Center, is one of the
   world's most powerful supercomputing centers, and you will be using it
   frequently for your research, especially the [NERSC Jupyter
   Hub server](https://jupyter.nersc.gov).


Finally, I recommend a good lab book (or you might consider keeping an
electronic journal or set of notes, ideally on Github) and enthusiasm!

Responsible Conduct of Research
-------------------------------

The grant which is supporting you this summer (whether from the National Science
Foundation or the Department of Energy) mandates that you complete an online
[Responsible Conduct of Research](https://about.citiprogram.org/en/homepage)
training program (see [these
instructions](https://www.siena.edu/files/resources/responsible-conduct-of-research-2016.pdf)),
in which you will learn about collaborative research, responsible authorship,
data management, etc.  *If you haven't completed the training in the previous
three years, it must be completed during the first week of your summer research
experience.*

Computers and more!
-------------------

The astrophysics research lab, Roger Bacon Hall 113, has three computers which
you are free to use (ask Moustakas for the punch code to the room and the login
to the computers). Each machine is a Dell Precision 3450 with 8 Intel i7
processors, 16 GB of RAM, 500 GB of SSD disk space, and a dual-boot Windows 10
and [Linux/Ubuntu](https://ubuntu.com/) operating system installed. In general,
you will be using the Linux/Ubuntu operating system, which is the standard
scientific computing environment in astronomy.

In addition to these computers, during your summer research experience you may
utilize NERSC, particularly the [NERSC Jupyter Hub
server](https://jupyter.nersc.gov), the [Google
Colaboratory](https://colab.research.google.com), my remote server
*nyx.siena.edu*, the Siena College high-performance computing cluster (HPCC),
and other computational tools.

Building your toolkit
---------------------

To start building your technical toolkit you will be working through several
lessons and tutorials from [Software Carpentry](https://software-carpentry.org/lessons) on other online resources.

* [The Unix Shell](http://swcarpentry.github.io/shell-novice)
* [Absolute Beginner's Guide to EMACS](http://www.jesshamrick.com/2012/09/10/absolute-beginners-guide-to-emacs)
* [Programming with Python](http://swcarpentry.github.io/python-novice-inflammation)
* [Version Control with Git](http://swcarpentry.github.io/git-novice) (see also [Git Immersion](http://gitimmersion.com/)) 

Finally, if you are interested in poking around at some of the
libraries and concepts we're going to be exploring this summer, you're
more than welcome to check these out--
  http://www.astroml.org/index.html
  https://www.tensorflow.org/tutorials
  https://github.com/desihub/timedomain/tree/master/desitrip


Additional resources
--------------------
Some additional resources you may find useful include:

* [Google's Python Class](https://developers.google.com/edu/python)
* [The Foundations of Data Science](https://ds8.gitbooks.io/textbook/content)
* [Scipy Lecture Notes](http://www.scipy-lectures.org/index.html)
* [Practical Python for Astronomers](https://python4astronomers.github.io)
* [Python Data Science Handbook](https://github.com/jakevdp/PythonDataScienceHandbook)
* [PEP 8 - Style Guide for Python Code](https://www.python.org/dev/peps/pep-0008/)
* [Learning Markdown](https://daringfireball.net/projects/markdown/syntax)

Papers and software packages
----------------------------

Write me!

This summer's research project focuses on understanding galaxy formation using
imaging and spectroscopy from the [Mapping Nearby Galaxies at APO
(MaNGA)](https://www.sdss.org/surveys/manga) survey. A non-exhaustive list of
relevant papers and software packages with which you will need to become
familiar include:

* [MaNGA Web portal](https://www.sdss.org/surveys/manga)
* [Marvin tool for accessing MaNGA data](https://www.sdss.org/dr15/manga/marvin)
* [Overview of the SDSS-IV MaNGA Survey: Mapping nearby Galaxies at Apache Point Observatory](https://ui.adsabs.harvard.edu/abs/2015ApJ...798....7B/abstract)

Need to add:
------------
* Logging into nyx.
* Jupyter notebook daily log.
* Anaconda/conda.
