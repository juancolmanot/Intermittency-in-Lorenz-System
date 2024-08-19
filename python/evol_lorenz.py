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


sigma, rho, beta = 10.0, 166.07, 8.0 / 3.0
params = [sigma, rho, beta]

x0 = [2.0, -1.0, 150.0]

t_start = 0.0
t_end = 50000.0
dt = 0.1
t_eval = np.arange(t_start, t_end, dt)
trajectory = np.empty((len(t_eval) + 2, 3))
trajectory[0] = x0

trajectory_file = sys.argv[1]

tprev = time.time()

with open(trajectory_file, 'w') as f:

	for i in range(1, len(t_eval)):
		t_span = [t_eval[i 	- 1], t_eval[i]]
		sol = solve_ivp(lorenz, t_span, trajectory[i - 1], args=(params,), jac=jac_lorenz, method='BDF')
		trajectory[i] = sol.y[:, -1]

		f.write(f'{t_eval[i]} {trajectory[i, 0]} {trajectory[i, 1]} {trajectory[i, 2]}\n')

		tnow = time.time()
		if (tnow - tprev > 2):
			print(f'Progress loop: {i * 100 / len(t_eval)} %')
			tprev = tnow