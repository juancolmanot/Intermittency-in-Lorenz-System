import sys
import numpy as np

evol1_file = sys.argv[1]
evol2_file = sys.argv[2]

evol_dat1 = np.loadtxt(evol1_file)
evol_dat2 = np.loadtxt(evol2_file)

for i in range(2000):
	print(evol_dat1[i, 1] - evol_dat2[i, 1])