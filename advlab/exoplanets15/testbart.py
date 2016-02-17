#!/usr/bin/env python

"""
Try out BART.  See http://dan.iel.fm/bart/current for details.

"""
#!/usr/bin/python
########################TESTING############################
import bart
import transit
from kplr.ld import get_quad_coeffs
import numpy as np
import matplotlib.pyplot as pl

##Set the Quadractic Limb Darkening

mu1,mu2 = get_quad_coeffs(5647.0, logg=4.236, feh = 0.34)

##Creates Stars and Planets

star = bart.Star(mass = 1.209, radius = 1.391, mu1 = mu1, mu2 = mu2)
planet = bart.Planet(r = 0.09829 * star.radius, period = 3.234723, b = 0.398)

##Creates Solar System

kepler6 = bart.PlanetarySystem(star)
kepler6.add_planet(planet)

##Try to plot the light curve(This is where it fails)

t = np.arange(120.0, 210.0, 0.5/ 24.)
lc = kepler6.light_curve(t)

pl.plot(t, lc, ".k")
pl.show()

