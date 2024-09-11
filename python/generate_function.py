import numpy as np
import sys

def read_parameters(param_file):
	with open(param_file, 'r') as file:
		params = {}
		for line in file:
			key, value = line.split(':')
			params[key.strip()] = float(value.strip())

		return params['xmin'], params['xmax'], int(params['n_points'])

def generate_function_data(
	xmin,
	xmax,
	n_points,
	function,
	output_file
):

	x = np.linspace(xmin, xmax, n_points)

	y = function(x)

	data = np.column_stack((x, y))

	np.savetxt(output_file, data, fmt='%.6f', delimiter='\t', header='x\tf(x)')
	print(f"Data written to {output_file}")

def function_computed(x):

	a, b = 4.1, -107.5
	#a, b = 4.1, -108.3

	return a * x + b

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: python generate_function.py parameters.txt output.dat")
		sys.exit(1)

	param_file = sys.argv[1]
	output_file = sys.argv[2]

	xmin, xmax, n_points = read_parameters(param_file)

	generate_function_data(xmin, xmax, n_points, function_computed, output_file)
	