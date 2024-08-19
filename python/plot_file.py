import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator, ScalarFormatter

def read_datafile(filename, x_col, y_col):
	data = np.loadtxt(filename)
	x = data[:, x_col]
	y = data[:, y_col]
	return x, y

def plot_from_instructions(instructions_file):
	with open(instructions_file, 'r') as file:
		lines = file.readlines()

	total_axes = 0
	layout = (1, 1)
	figsize = (5, 5)
	axes_instructions = []
	save_file = 'output_plot.pdf'

	current_axes = None
	for line in lines:
		line = line.strip()
		if line.startswith("total_axes:"):
			total_axes = int(line.split(":")[1].strip())
		elif line.startswith("layout:"):
			layout = tuple(map(int, line.split(":")[1].strip().split(',')))
		elif line.startswith("figsize:"):
			figsize = tuple(map(float, line.split(":")[1].strip().split(',')))
		elif line.startswith("axes"):
			if current_axes is not None:
				axes_instructions.append(current_axes)
			current_axes = {"plots": 0, "data": [], "settings": {}}
		elif line.startswith("plots:"):
			current_axes["plots"] = int(line.split(":")[1].strip())
		elif line.startswith("file:"):
			data_info = {"file": line.split(":")[1].strip()}
			current_axes["data"].append(data_info)
		elif line.startswith("method:"):
			if current_axes["data"]:
				current_axes["data"][-1]["method"] = line.split(":")[1].strip()
			else:
				print("Error: 'x_col' specified withou")
		elif line.startswith("x_col:"):
			if current_axes["data"]:
				current_axes["data"][-1]["x_col"] = int(line.split(":")[1].strip())
			else:
				print("Error: 'x_col' specified without a preceding 'file'")
		elif line.startswith("y_col:"):
			if current_axes["data"]:
				current_axes["data"][-1]["y_col"] = int(line.split(":")[1].strip())
			else:
				print("Error: 'y_col' specified without a preceding 'file'")
		elif line.startswith("label:"):
			if current_axes["data"]:
				current_axes["data"][-1]["label"] = line.split(":")[1].strip()
			else:
				print("Error: 'label' specified without a preceding 'file'")
		elif line.startswith("marker:"):
			if current_axes["data"]:
				current_axes["data"][-1]["marker"] = line.split(":")[1].strip()
			else:
				print("Error: 'marker' specified without a preceding 'file'")
		elif line.startswith("every_n"):
			if current_axes["data"]:
				current_axes["data"][-1]["every_n"] = int(line.split(":")[1].strip())
		elif line.startswith("title:"):
			current_axes["settings"]["title"] = line.split(":")[1].strip()
		elif line.startswith("xlabel:"):
			current_axes["settings"]["xlabel"] = line.split(":")[1].strip()
		elif line.startswith("ylabel:"):
			current_axes["settings"]["ylabel"] = line.split(":")[1].strip()
		elif line.startswith("color:"):
			if current_axes["data"]:
				current_axes["data"][-1]["color"] = line.split(":")[1].strip()
			else:
				print("Error: 'color' specified without a preceding 'file'")
		elif line.startswith("dot_size:"):
			if current_axes["data"]:
				current_axes["data"][-1]["dot_size"] = float(line.split(":")[1].strip())
			else:
				print("Error: 'dot_size' specified without a preceding 'file'")
		elif line.startswith("linewidth:"):
			if current_axes["data"]:
				current_axes["data"][-1]["linewidth"] = float(line.split(":")[1].strip())
			else:
				print("Error: 'linewidth' specified without a preceding 'file'")
		elif line.startswith("xlim:"):
			current_axes["settings"]["xlim"] = eval(line.split(":")[1].strip())
		elif line.startswith("ylim:"):
			current_axes["settings"]["ylim"] = eval(line.split(":")[1].strip())
		elif line.startswith("xscale:"):
			current_axes["settings"]["xscale"] = line.split(":")[1].strip()
		elif line.startswith("yscale"):
			current_axes["settings"]["yscale"] = line.split(":")[1].strip()
		elif line.startswith("save_file:"):
			save_file = line.split(":")[1].strip()

	if current_axes is not None:
		axes_instructions.append(current_axes)

	fig, axes = plt.subplots(layout[0], layout[1], figsize=figsize)
	if layout[0] == 1 and layout[1] == 1:
		axes = [axes]
	elif layout[0] == 1 or layout[1] == 1:
		axes = axes.flatten()
	else:
		axes = axes.ravel()

	for ax_index, ax_info in enumerate(axes_instructions):
		ax = axes[ax_index]
		for plot_info in ax_info['data']:
			x, y = read_datafile(plot_info['file'], plot_info['x_col'], plot_info['y_col'])

			every_n = plot_info.get('every_n', 1)
			x, y = x[::every_n], y[::every_n]

			method = plot_info.get('method', 'scatter')
			marker = plot_info.get('marker', 'o')
			if method == 'plot':
				ax.plot(x, y, color=plot_info['color'], lw=plot_info['linewidth'], label=plot_info.get('label', ''))
			elif method == 'scatter':
				ax.scatter(x, y, color=plot_info['color'], s=plot_info['dot_size'], label=plot_info.get('label', ''), marker=marker)


		font_properties = {'family': 'serif', 'size': 14}

		'''ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
		ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
		ax.xaxis.get_major_formatter().set_scientific(True)
		ax.xaxis.get_major_formatter().set_powerlimits((0, 0))

		ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
		ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
		ax.yaxis.get_major_formatter().set_scientific(True)
		ax.yaxis.get_major_formatter().set_powerlimits((0, 0))'''

		# Set x-scale and y-scale if provided:
		if 'xscale' in ax_info['settings']:
			ax.set_xscale(ax_info['settings']['xscale'])
		if 'yscale' in ax_info['settings']:
			ax.set_yscale(ax_info['settings']['yscale'])

		# Set xlim and ylim if provided
		if 'xlim' in ax_info['settings']:
			ax.set_xlim(ax_info['settings']['xlim'])
		if 'ylim' in ax_info['settings']:
			ax.set_ylim(ax_info['settings']['ylim'])

		ax.set_title(ax_info['settings'].get('title', ''), **font_properties)
		ax.set_xlabel(ax_info['settings'].get('xlabel', ''), **font_properties)
		ax.set_ylabel(ax_info['settings'].get('ylabel', ''), **font_properties)
		ax.legend()

	plt.tight_layout()
	plt.savefig(save_file)

if __name__ == "__main__":
	inst_file = sys.argv[1]
	plot_from_instructions(inst_file)
