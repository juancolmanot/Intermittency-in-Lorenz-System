import sys
import time
import numpy as np
from scipy.integrate import odeint, solve_ivp

def lorenz(t, x, p):
	sigma, rho, beta = p
	dfdx = np.zeros(3)

	dfdx[0] = sigma * (x[1] - x[0])
	dfdx[1] = ((rho - x[2]) * x[0]) - x[1]
	dfdx[2] = (x[0] * x[1]) - beta * x[2]

	return dfdx

def jac_lorenz(t, x, p):
	sigma, rho, beta = p
	
	jac = np.zeros((3, 3))
	
	jac[0, 0] = -sigma
	jac[0, 1] = sigma
	jac[0, 2] = 0.0
	jac[1, 0] = rho - x[2]
	jac[1, 1] = -1.0
	jac[1, 2] = - x[0]
	jac[2, 0] = x[1]
	jac[2, 1] = x[0]
	jac[2, 2] = - beta

	return jac

def reinjection_1d(x, x1, xf, c):

	if (x < xf - c or x > xf + c):
		if (xf - c < x1 < xf + c):
			return True

	return False

sigma, rho, beta = 10.0, 166.07, 8.0 / 3.0
params = [sigma, rho, beta]

# Evol transient
t_start = 0.0
t_end = 500.0
dt = 0.01

t_eval = np.arange(t_start, t_end, dt)

x0 = [2.0, -1.0, 150.0]
sol = solve_ivp(lorenz, (t_start, t_end), x0, args=(params,), t_eval=t_eval, jac=jac_lorenz, method='BDF')

# Integrate stationary state
t_start = 0.0
t_end = 50000.0
dt = 0.01

t_eval = np.arange(t_start, t_end, dt)

x0 = sol.y[:, -1]
trajectory = np.empty((len(t_eval) + 2, 3))
trajectory[0:3] = sol.y[:, -4:-1]

xreg = [[0, 0]]
target = 1000
mapcount = 0
xp = 0

xi, yi, zi = trajectory[0:3]
xfit, yfit, zfit = np.zeros(3), np.zeros(3), np.zeros(3)

tprev = time.time()

xreinj = []
yf = 41.2861
clam = 1.85
reinjcount = 0

for i in range(1, len(t_eval)):
	t_span = [t_eval[i 	- 1], t_eval[i]]
	sol = solve_ivp(lorenz, t_span, trajectory[i - 1], args=(params,), jac=jac_lorenz, method='BDF')
	trajectory[i] = sol.y[:, -1]

	xi[0:2], yi[0:2], zi[0:2] = xi[1:3], yi[1:3], zi[1:3]
	xi[2], yi[2], zi[2] = trajectory[i]

	if (xi[0] < xp and xi[1] > xp):

		a1, b1, c1 = np.polyfit(xi, yi, 2)
		a2, b2, c2 = np.polyfit(xi, zi, 2)
		xreg.append([a1 * xp * xp + b1 * xp + c1, a2 * xp * xp + b2 * xp + c2])
		yr, _ = xreg[mapcount]
		yr1, _ = xreg[mapcount + 1]
		mapcount += 1

		if (reinjection_1d(yr, yr1, yf, clam)):
			xreinj.append([yr, yr1])
			reinjcount += 1

			if (reinjcount >= target):
				break

	tnow = time.time()
	if (tnow - tprev > 2):
		print(f'Progress reinj: {reinjcount * 100 / target} %')
		print(f'Progress loop: {i * 100 / len(t_eval)} %')
		tprev = tnow


reinj_file = sys.argv[1]

print(reinjcount)

with open(reinj_file, 'w') as f:

	for i, xri in enumerate(xreinj):
		f.write(f'{i} {xri[0]} {xri[1]}\n')