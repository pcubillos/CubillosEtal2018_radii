#! /usr/bin/env python 

import sys
import numpy as np

sys.path.append("../pyratbay")
import pyratbay as pb

for temp in np.arange(300, 3001, 100):
  pb.pbay.run("atm_solar_{:04d}K.cfg".format(temp))
