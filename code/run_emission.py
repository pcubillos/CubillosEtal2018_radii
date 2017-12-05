#! /usr/bin/env python
import sys
import os
import ctypes
import numpy as np
import multiprocessing as mp

sys.path.append("../code")
import sma as sma

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.tools      as pt
import pyratbay.atmosphere as pa


def worker(pyrat, p0, sm_flag, niter, Hhill, itemp, metal, Rp, Mp, Teq, ID):
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
  niter: 3D shared-memory integer ndarray
     Number of iterations needed for each grid model [nrad,nmass,ntemp].
  Hhill: 3D shared-memory float ndarray
     Height of Hill radius over planetary radius in scale heights.
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
  # Grid shape:
  nr, nm, nt, nf = np.shape(p0)
  flag = np.ctypeslib.as_array(sm_flag.get_obj()).reshape((nm,nt))

  # Max number of iterations per model:
  nmax = 7
  ph = np.zeros(nmax,   np.double) # Hill radius pressure
  rp = np.zeros(nmax+1, np.double) # Reference pressure at each iteration

  # Density:
  dens = Mp*pc.mearth/(4./3*np.pi*np.expand_dims(Rp*pc.rearth,1)**3)
  # Density flag:
  dflag = np.ones((len(Rp), len(Mp)), bool)
  dflag[np.where(dens <  0.03)] = 0
  dflag[np.where(dens > 30.00)] = 0

  # Compute orbital semi-major axis from Teq and stellar mass:
  smaxis = sma.sma(Teq, pyrat.phy.mstar/pc.msun) * pc.au

  # Starting point:
  imass = nm-1
  itemp0 = -1

  while True:
    # Search next available [itemp,imass] to run (use lock):
    with sm_flag.get_lock():
      if np.sum(flag) == nt*nm:
        break
      # Search for next imass in current itemp:
      if np.any(flag[0:imass,itemp] == 0):
        while flag[imass,itemp] == 1:
         imass = (imass-1)%nm
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
          imass = (imass-1)%nm
      # Update flag:
      flag[imass,itemp] = 1

    # Update atmosphere:
    if itemp != itemp0:
      spec, p, temp, q = pa.readatm("../run01_atm/{:s}_{:04.0f}K.atm".
                                    format(metal, Teq[itemp]))
      pyrat.phy.smaxis = smaxis[itemp]
      itemp0 = itemp

    pyrat.phy.mplanet = Mp[imass] * pc.mearth
    pyrat.phy.rhill   = pyrat.phy.smaxis * (pyrat.phy.mplanet /
                                           (3*pyrat.phy.mstar))**(1.0/3.0)

    # Compute p0 for all radii in this [imass,itemp] cell:
    for irad in np.arange(nr):
      # Invalid density:
      if not dflag[irad, imass]:
        continue
      if ID == 0:
        print("\n::::::::::::::::::::::::::::::::::::\n"
              "r={:02d}/{:02d}  m={:03d}/{:03d}  t={:02d}/{:02d}  {:s}\n"
              "::::::::::::::::::::::::::::::::::::".
              format(irad, nr-1, imass, nm-1, itemp, nt-1, metal))
      pyrat.phy.rplanet = Rp[irad] * pc.rearth
      if pyrat.phy.rplanet > pyrat.phy.rhill:
        break
      # Initial guess for the pressure reference level:
      rp[:] = 0.1*pc.bar  # Default
      if   irad != 0     and p0[irad-1,imass,itemp,0] != 0.0:  # From prev Rp
        rp[:] = p0[irad-1,imass,itemp,0]
      elif imass != nm-1 and p0[irad,imass+1,itemp,0] != 0.0:  # From prev Mp
        rp[:] = p0[irad,imass+1,itemp,0]
      # Iterate (nmax times) untill ref. pressure converges:
      for i in np.arange(nmax):
        # Compute spectra:
        pyrat.refpressure = rp[i]
        pyrat = pb.pyrat.run(pyrat, [temp, q])
        cf  = pt.cf(pyrat.od.depth, pyrat.atm.press, pyrat.od.B)
        bcf = pt.bandcf(cf, pyrat.obs.bandtrans, pyrat.spec.wn,
                        pyrat.obs.bandidx)
        rp[i+1] = 10**np.average(np.log10(press), weights=bcf[0])
        ph[i]   = press[pyrat.atm.rtop]
        if ID == 0:
          print("Photosphere: {:.4f} bar  --  dH={:6.3f}  --  update={:6.2f}%".
              format(rp[i+1]/pc.bar, np.log(rp[i+1]/ph[i]),
                     100.0*np.abs(rp[i+1]-rp[i])/rp[i+1]))
        # Ref. pressure converged (improvement smaller than 2% in pressure):
        if i>=2 and np.abs(rp[i+1]-rp[i])/rp[i+1] < 0.02:
          break
      # Save results:
      for j in np.arange(pyrat.obs.nfilters):
        p0 [irad,imass,itemp,j] = 10**np.average(np.log10(press),
                                                 weights=bcf[j])
      Hhill[irad,imass,itemp] = np.log(rp[i+1]/ph[i])
      niter[irad,imass,itemp] = i
      if np.log(rp[i+1]/ph[i]) < 2.0: # Break radii loop if < 2 scale heights
        break
    # Final save:
    if ID == 0:
      np.savez("MRT_{:s}.npz".format(metal), Rp=Rp, Mp=Mp, Teq=Teq,
             p0=p0, niter=niter, Hhill=Hhill)
    print("Completed {}/{} [ID={}, {}/{}]".format(np.sum(flag), nm*nt,
                                                  ID, itemp, imass))


# Load Pyrat object:
pyrat   = pb.pyrat.init("MRTz_grid.cfg")
press   = pyrat.atm.press    # Atmospheric pressure array (barye)
nlayers = pyrat.atm.nlayers
pyrat.verb = 0  # Mute it

# Metallicity:
metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
nz = len(metal)

# Create new grid:
nf = pyrat.obs.nfilters  # Number of observing filter
nt = 28  # Number of temperature values
nr = 79  # Number of planetary radius values
nm = 45  # Number of planetary mass values
# The arrays:
Teq = np.linspace(300, 3000, nt)  # K
Rp  = np.linspace(1.0, 40.0, nr)  # Rearth
Mp  = np.logspace(np.log10(1.0), np.log10(630), 45)  # Mearth
Mp[-1]  = 600.0

# Number of parallel CPUs:
ncpu = 23

# Shared memory arrays:
sm_p0    = mp.Array(ctypes.c_double, nr*nm*nt*nf)
sm_flag  = mp.Array(ctypes.c_int,       nm*nt)  # ctypeslib fails for c_bool
sm_niter = mp.Array(ctypes.c_int,    nr*nm*nt)
sm_Hhill = mp.Array(ctypes.c_double, nr*nm*nt)
# See description in the worker() docstring:
p0    = np.ctypeslib.as_array(sm_p0.get_obj()   ).reshape((nr,nm,nt,nf))
flag  = np.ctypeslib.as_array(sm_flag.get_obj() ).reshape((   nm,nt))
niter = np.ctypeslib.as_array(sm_niter.get_obj()).reshape((nr,nm,nt))
Hhill = np.ctypeslib.as_array(sm_Hhill.get_obj()).reshape((nr,nm,nt))

# Main loop:
for z in np.arange(nz):
  # Reset arrays for each metallicity
  p0   [:] = 0.0
  flag [:] = 0
  niter[:] = 0
  Hhill[:] = 0.0

  # Start subprocesses:
  procs = []
  for i in np.arange(ncpu):
    i0 = int(i*nt/ncpu)
    p = mp.Process(target=worker, args=(pyrat, p0, sm_flag, niter, Hhill, i0,
                metal[z], Rp, Mp, Teq, i))
    p.start()
    procs.append(p)
  # Wait until they probe the whole grid:
  for i in np.arange(ncpu):
    procs[i].join()
  # Final save:
  np.savez("MRT_{:s}.npz".format(metal[z]), Rp=Rp, Mp=Mp, Teq=Teq,
           p0=p0, niter=niter, Hhill=Hhill)
