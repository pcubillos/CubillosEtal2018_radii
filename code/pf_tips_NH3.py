#! /usr/bin/env python

import sys
import numpy as np

sys.path.append("../pyratbay/modules/pytips")
import pytips as p

molname = "NH3"
temp    = np.linspace(70, 3000, 294)
p.to_file("./PF_tips_NH3.dat", molname, temp)
