#!/usr/bin/env python

"""
Try out BART.  See http://dan.iel.fm/bart/current for details.

"""

import bart
import numpy as np

# Set up a star with limb darkening.
star = bart.Star(ldp=bart.ld.QuadraticLimbDarkening(0.3, 0.1, 100))

# Set up a planet somewhat like the Earth.
planet = bart.Planet(r=0.01, period=365.0)

# Put it all together.
solar_system = bart.PlanetarySystem(star)
solar_system.add_planet(planet)

# "Observe" the light curve.
t = np.arange(364.0, 366.0, 0.02)
lc = solar_system.light_curve(t)
