#!/usr/bin/env python

# Katie Hoag
# Code to test pyephem

import ephem.stars
import ephem

for star in ephem.stars.db.split("\n"):
  print star.split(",")[0]

u=ephem.Uranus('2016/5/24')
print("%s %s %s" % (u.ra, u.dec, u.mag))
