import numpy as np
import scipy.interpolate as si


def sma(temp, ms):
  """
  Compute semi-major axis for requested equilibrium temperatures and
  stellar masses, extrapolating from main-sequence stellar values.

  Parameters
  ----------
  temp: Float or 1D float ndarray
     Planetary equilibrium temperature (K).
  ms: Float or 1D float ndarray
     Stellar mass (Msun).

  Return
  ------
  a: 1D float ndarray
     Semi-major axis (AU).
  """
  # Stellar mass array in solar masses:
  Ms = np.linspace(0.4, 1.3, 10)
  # Semi-major axis at 250 K for stellar masses from 0.4--1.3 (dm=0.1):
  sma250 = np.array([0.13340, 0.22680, 0.34865, 0.49505, 0.69795,
                     1.01260, 1.45355, 1.66855, 1.91340, 2.18890])
  # SMA interpolator at 250 K:
  a250 = si.interp1d(Ms, sma250)

  # Check data types:
  if       isinstance(temp, (int,float)) and     isinstance(ms, (int,float)):
    temp = np.asarray(temp)
    ms   = np.asarray(ms)
  elif     isinstance(temp, (int,float)) and not isinstance(ms, (int,float)):
    temp = np.tile(temp, len(ms))
  elif not isinstance(temp, (int,float)) and     isinstance(ms, (int,float)):
    ms   = np.tile(ms, len(temp))
  elif len(temp) != len(ms):
    print("Size of input arrays do not match!")
    return None

  # Out of bounds error:
  if np.any(ms < Ms[0]) or np.any(ms > Ms[-1]):
    print("Out of bound inputs!")
    return None

  return a250(ms) * (temp/250.0)**(-2)
