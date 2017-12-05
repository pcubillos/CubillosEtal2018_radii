#! /usr/bin/env python

import sys
import os
import ctypes
import numpy as np
import scipy.interpolate as si
import multiprocessing as mp

sys.path.append("../code")
import sma   as sma

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.wine       as pw
import pyratbay.atmosphere as pa


def worker(pyrat, ptransit, rtransit, sm_flag, itemp, metal, Rp, Mp, Teq, ID):
  """
  Compute the photospheric pressure (p0) corresponding to Rp for each
  model in the grid.

  This function finds a [mass,temp] cell that has not been previously
  computed (flag[m,t]==0) and then computes p0 for the whole range of Rp.

  Parameters
  ----------
  pyrat: Pyrat object
     The atmospheric model object
  p0: 4D shared-memory float ndarray
     Output photospheric pressure array [nrad,nmass,ntemp,nfilter].
  sm_flag: multipprocessing shared-memory integer Array
     Grid's completion flag [nmass,ntemp].
  itemp: Integer
     Initial temperature index in the grid.
  metal: String
  Rp: 1D float ndarray
     The grid's radius array.
  Mp: 1D float ndarray
     The grid's mass array.
  Teq: 1D float ndarray
     The grid's temperature array.
  ID: Integer
     The subprocess' ID index.
  """
  # Load grid:
  grid = np.load("./MRT_{:s}.npz".format(metal))
  p0, Hhill = grid["p0"], grid["Hhill"]
  nr, nm, nt, nf = np.shape(p0)

  # Use this flag to synchronize processes:
  flag = np.ctypeslib.as_array(sm_flag.get_obj()).reshape((nm,nt))

  # Compute orbital semi-major axis from Teq and stellar mass:
  smaxis = sma.sma(Teq, pyrat.phy.mstar/pc.msun) * pc.au

  # Starting point:
  imass  =  0
  itemp0 = -1

  while True:
    # Search next available [itemp,imass] to run (use lock):
    with sm_flag.get_lock():
      if np.sum(flag) == nt*nm:
        break
      # Search for next imass in current itemp:
      if np.any(flag[:,itemp] == 0):
        while flag[imass,itemp] == 1:
         imass = (imass+1)%nm
      # Find another itemp:
      else:
        if np.any(np.all(flag==0, axis=0)):
          while np.any(flag[:,itemp] != 0): # First search empty itemp
            itemp = (itemp-1)%nt
        else:
          while np.all(flag[:,itemp] == 1): # Else, just get next available
            itemp = (itemp-1)%nt
        # Now look for imass:
        while flag[imass,itemp] == 1:
          imass = (imass+1)%nm
      # Update flag:
      flag[imass,itemp] = 1

    # Update atmosphere:
    if itemp != itemp0:
      spec, p, temp, q = pa.readatm("../run01_atm/{:s}_{:04.0f}K.atm".
                                    format(metal, Teq[itemp]))
      pyrat.phy.smaxis = smaxis[itemp]
      mm = np.sum(q*pyrat.mol.mass, axis=1)
      itemp0 = itemp

    pyrat.phy.mplanet = Mp[imass] * pc.mearth
    pyrat.phy.rhill   = pyrat.phy.smaxis * (pyrat.phy.mplanet /
                                           (3*pyrat.phy.mstar))**(1.0/3.0)

    # Compute ptransit/rtransit for all radii in this [imass,itemp] cell:
    for irad in np.arange(nr):
      if Hhill[irad,imass,itemp] < 2.0: # or not flag[irad,imass,itemp]:
        continue
      if ID == 0:
        print("\n::::::::::::::::::::::::::::::::::::\n"
              "r={:02d}/{:02d}  m={:03d}/{:03d}  t={:02d}/{:02d}  {:s}\n"
              "::::::::::::::::::::::::::::::::::::".
              format(irad, nr-1, imass, nm-1, itemp, nt-1, metal))
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
      rtransit[irad,imass,itemp] = rprs*pyrat.phy.rstar
      pinterp = si.interp1d(radius, press)
      ptransit[irad,imass,itemp] = pinterp(rtransit[irad,imass,itemp])
    if ID == 0:
      np.savez("rtransit_{:s}.npz".format(metal), Rt=rtransit, pt=ptransit)
    print("Completed {}/{} [ID={}, {}/{}]".format(np.sum(flag), nm*nt,
                                                    ID, itemp, imass))

# Load Pyrat object:
pyrat = pb.pyrat.init("MRTz_grid.cfg")
press = pyrat.atm.press
pyrat.od.path = "transit"  # Reset geometry
pyrat.verb    = 0          # Mute it

# Metallicity:
metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
nz = len(metal)

# Input grid arrays and shapes:
grid = np.load("./MRT_{:s}.npz".format(metal[0]))
Rp, Mp, Teq = grid["Rp"], grid["Mp"], grid["Teq"]
nr, nm, nt, nf = np.shape(grid["p0"])

# Number of parallel CPUs:
ncpu = 23

# Shared memory arrays:
sm_ptransit = mp.Array(ctypes.c_double, nr*nm*nt*nf)
sm_rtransit = mp.Array(ctypes.c_double, nr*nm*nt*nf)
sm_flag     = mp.Array(ctypes.c_int,       nm*nt)
# See description in the worker() docstring:
ptransit = np.ctypeslib.as_array(sm_ptransit.get_obj()).reshape((nr,nm,nt,nf))
rtransit = np.ctypeslib.as_array(sm_rtransit.get_obj()).reshape((nr,nm,nt,nf))
flag     = np.ctypeslib.as_array(sm_flag.get_obj()    ).reshape((   nm,nt))

# Now with known p0, run transmission to get the transit radius and pressure:
for z in np.arange(nz):
  # Reset outputs:
  rtransit[:] = 0.0
  ptransit[:] = 0.0
  flag    [:] = 0

  # Start subprocesses:
  procs = []
  for i in np.arange(ncpu):
    i0 = int(i*nt/ncpu)
    p = mp.Process(target=worker, args=(pyrat, ptransit, rtransit, sm_flag,
                   i0, metal[z], Rp, Mp, Teq, i))
    p.start()
    procs.append(p)
  # Wait until they probe the whole grid:
  for i in np.arange(ncpu):
    procs[i].join()
  # Final save:
  np.savez("rtransit_{:s}.npz".format(metal[z]), Rt=rtransit, pt=ptransit)
