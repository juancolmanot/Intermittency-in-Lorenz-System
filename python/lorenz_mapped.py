import sys
import numpy as np

trajectory_file = sys.argv[1]
mapped_file = sys.argv[2]

x_trajectory = np.loadtxt(trajectory_file)

xi = np.zeros(3)
yi = np.zeros(3)
zi = np.zeros(3)

xfit, yfit, zfit = np.zeros(3), np.zeros(3), np.zeros(3)

yreg, zreg = [], []

target = 50000

rcount = 0

xp = 0

for i, xti in enumerate(x_trajectory):

	xi[0:2], yi[0:2], zi[0:2] = xi[1:3], yi[1:3], zi[1:3]

	xi[2], yi[2], zi[2] = xti[1], xti[2], xti[3]

	if (xi[0] < xp and xi[1] > xp):
		
		rcount += 1

		a1, b1, c1 = np.polyfit(xi, yi, 2)
		a2, b2, c2 = np.polyfit(xi, zi, 2)

		yreg.append(a1 * xp * xp + b1 * xp + c1)
		zreg.append(a2 * xp * xp + b2 * xp + c2)

		if (rcount >= target):
			break

print(rcount)

with open(mapped_file, 'w') as f:
	for i in range(len(yreg)):

		f.write(f'{yreg[i]} {zreg[i]}\n')