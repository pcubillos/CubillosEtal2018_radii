#! /usr/bin/env python
import sys
import numpy as np
import scipy.constants as sc
import scipy.optimize  as so
import matplotlib.pyplot as plt

sys.path.append("./pyratbay")
import pyratbay.constants as pc


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

grid = np.load("run02_grid/MRT_solar.npz")
Rp, Mp, Teq = grid["Rp"], grid["Mp"], grid["Teq"]

metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
nz = len(metal)

# Restricted Jeans escape parameter:
mh = sc.m_p + sc.m_e
R = Rp [:,np.newaxis,np.newaxis]
M = Mp [np.newaxis,:,np.newaxis]
T = Teq[np.newaxis,np.newaxis,:]
Lambda = 0.1*pc.mearth/pc.rearth * sc.G*M*mh/(sc.k*T*R)


def powlaw(p, Rp, Mp, Teq):
  mesh = np.meshgrid(Mp, Rp, Teq)
  return p[0] * mesh[1]**p[1] * mesh[0]**p[2] * mesh[2]**p[3]

def error(p, Rp, Mp, Teq, y, Lambda0=0.0):
  diff = y - powlaw(p, Rp, Mp, Teq/300.0)
  # Ignore points where y=0
  diff[~np.isfinite(diff)] = 0.0
  diff[y==0] = 0.0
  diff[Lambda < Lambda0] = 0.0
  return diff.flatten()
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Lambda0 = 10.0
pfit = np.zeros((nz,4))
tfit = np.zeros((nz,4))

for z in np.arange(len(metal)):
  # Photosphere pressures:
  grid = np.load("run02_grid/MRT_{:s}.npz".format(metal[z]))
  p0 = grid["p0"][:,:,:,0]/pc.bar
  guess = [0.1, 0.0, 0.0, 0.0]
  fit = so.leastsq(error, guess, args=(Rp,Mp,Teq,p0,Lambda0), full_output=True)
  pfit[z] = fit[0]
  # Transit pressure:
  grid = np.load("run02_grid/rtransit_{:s}.npz".format(metal[z]))
  pt = grid["pt"][:,:,:,0]/pc.bar
  guess = [1.0, 0.0, 0.0, 0.0]
  fit = so.leastsq(error, guess, args=(Rp,Mp,Teq,pt,Lambda0), full_output=True)
  tfit[z] = fit[0]

np.set_printoptions(precision=2)
print("Photosphere:\n{}".format(pfit))
print("\nTransit:\n{}".format(tfit))
np.set_printoptions(precision=8)

