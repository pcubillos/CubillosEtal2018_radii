#! /usr/bin/env python

import sys, os
import numpy as np
import scipy.interpolate as si

sys.path.append("../code")
import sma   as sma

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.wine      as pw


pyrat = pb.pyrat.init("MRTz_grid.cfg")
pyrat.od.path = "transit"  # Set geometry
pyrat.verb    = 0          # Mute it
press = pyrat.atm.press

metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
nz = len(metal)

# Now that we know p0, run transmission RT to get the transit radius:
for z in np.arange(nz):
  # Load grid:
  grid = np.load("./MRT_{:s}.npz".format(metal[z]))
  Rp, Mp, Teq     = grid["Rp"], grid["Mp"],   grid["Teq"]
  p0, flag, Hhill = grid["p0"], grid["flag"], grid["Hhill"]
  nr = len(Rp)
  nm = len(Mp)
  nt = len(Teq)
  nf = pyrat.obs.nfilters

  # Planet transit radius and pressure:
  Rtransit = np.zeros(np.shape(p0))
  ptransit = np.zeros(np.shape(p0))
  smaxis = sma.sma(Teq, pyrat.phy.mstar/pc.msun) * pc.au
  for itemp in np.arange(nt):
    spec, p, temp, q = pb.atmosphere.readatm("../run01_atm/{:s}_{:04.0f}K.atm".
                                           format(metal[z], Teq[itemp]))
    pyrat.phy.smaxis  = smaxis[itemp]
    mm = np.sum(q*pyrat.mol.mass, axis=1)
    for imass in np.arange(nm):
      pyrat.phy.mplanet = Mp[imass] * pc.mearth
      pyrat.phy.rhill   = pyrat.phy.smaxis * (pyrat.phy.mplanet /
                                              (3*pyrat.phy.mstar))**(1.0/3.0)
      for irad in np.arange(nr):
        if not flag[irad,imass,itemp] or Hhill[irad,imass,itemp]<2.0:
          continue
        print("\n::::::::::::::::::::::::::::::::::::\n"
              "r={:02d}/{:02d}  m={:03d}/{:03d}  t={:02d}/{:02d}  {:s}\n"
              "::::::::::::::::::::::::::::::::::::".
              format(irad, nr-1, imass, nm-1, itemp, nt-1, metal[z]))
        pyrat.phy.rplanet = Rp[irad] * pc.rearth
        # pressure reference level:
        pyrat.refpressure = p0[irad,imass,itemp,0]
        radius = pyrat.hydro(press, temp, mm, pyrat.phy.gplanet,
                      pyrat.phy.mplanet, pyrat.refpressure, pyrat.phy.rplanet)
        radius[radius<0] = np.inf
        # Compute spectra:
        pyrat = pb.pyrat.run(pyrat, [temp, q, radius])
        # Band-integrated Rp/Rs:
        rprs = np.sqrt(pw.bandintegrate(pyrat=pyrat))
        # Save results:
        Rtransit[irad,imass,itemp] = rprs*pyrat.phy.rstar
        pinterp = si.interp1d(radius, press)
        ptransit[irad,imass,itemp] = pinterp(Rtransit[irad,imass,itemp])

      np.savez("rtransit_{:s}.npz".format(metal[z]), Rt=Rtransit, pt=ptransit)
