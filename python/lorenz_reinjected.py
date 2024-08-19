import sys
import numpy as np

mapped_file = sys.argv[1]
reinjection_file = sys.argv[2]

x_mapped = np.loadtxt(mapped_file)

def reinjection_1d(x, x1, xf, c):

	if (x < xf - c or x > xf + c):
		if (xf - c < x1 < xf + c):
			return True

	return False

xreinj = []

yf = 41.2861
clam = 1.85

for i in range(len(x_mapped) - 1):

	yr, _ = x_mapped[i]
	yr1, _ = x_mapped[i + 1]

	if (reinjection_1d(yr, yr1, yf, clam)):
		xreinj.append([yr, yr1])

with open(reinjection_file, 'w') as f:

	for i, yri in xreinj:
		f.write(f'{i} {yri[0]} {yri[1]}\n')