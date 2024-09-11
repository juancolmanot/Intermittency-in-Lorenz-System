import sys
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def read_datafile(filename, x_col, y_col, z_col=None):
    data = np.loadtxt(filename)
    x = data[:, x_col]
    y = data[:, y_col]
    if z_col is not None:
        z = data[:, z_col]
        return x, y, z
    return x, y

def plot_from_instructions(instructions_file):
    with open(instructions_file, 'r') as file:
        lines = file.readlines()

    total_axes = 0
    layout = (1, 1)
    figsize = (5, 5)
    axes_instructions = []
    save_file = 'output_plot.pdf'
    elevation = 30  # Default elevation angle
    azimuth = 45    # Default azimuth angle

    current_axes = None
    for line in lines:
        line = line.strip()
        if line.startswith("total_axes:"):
            total_axes = int(line.split(":")[1].strip())
        elif line.startswith("layout:"):
            layout = tuple(map(int, line.split(":")[1].strip().split(',')))
        elif line.startswith("figsize:"):
            figsize = tuple(map(float, line.split(":")[1].strip().split(',')))
        elif line.startswith("elevation:"):
            elevation = float(line.split(":")[1].strip())
        elif line.startswith("azimuth:"):
            azimuth = float(line.split(":")[1].strip())
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
        elif line.startswith("x_col:"):
            if current_axes["data"]:
                current_axes["data"][-1]["x_col"] = int(line.split(":")[1].strip())
        elif line.startswith("y_col:"):
            if current_axes["data"]:
                current_axes["data"][-1]["y_col"] = int(line.split(":")[1].strip())
        elif line.startswith("z_col:"):
            if current_axes["data"]:
                current_axes["data"][-1]["z_col"] = int(line.split(":")[1].strip())
        elif line.startswith("label:"):
            if current_axes["data"]:
                current_axes["data"][-1]["label"] = line.split(":")[1].strip()
        elif line.startswith("color:"):
            if current_axes["data"]:
                current_axes["data"][-1]["color"] = line.split(":")[1].strip()
        elif line.startswith("dot_size:"):
            if current_axes["data"]:
                current_axes["data"][-1]["dot_size"] = float(line.split(":")[1].strip())
        elif line.startswith("linewidth:"):
            if current_axes["data"]:
                current_axes["data"][-1]["linewidth"] = float(line.split(":")[1].strip())
        elif line.startswith("title:"):
            current_axes["settings"]["title"] = line.split(":")[1].strip()
        elif line.startswith("xlabel:"):
            current_axes["settings"]["xlabel"] = line.split(":")[1].strip()
        elif line.startswith("ylabel:"):
            current_axes["settings"]["ylabel"] = line.split(":")[1].strip()
        elif line.startswith("zlabel:"):
            current_axes["settings"]["zlabel"] = line.split(":")[1].strip()
        elif line.startswith("xlim:"):
            current_axes["settings"]["xlim"] = eval(line.split(":")[1].strip())
        elif line.startswith("ylim:"):
            current_axes["settings"]["ylim"] = eval(line.split(":")[1].strip())
        elif line.startswith("zlim:"):
            current_axes["settings"]["zlim"] = eval(line.split(":")[1].strip())
        elif line.startswith("xscale:"):
            current_axes["settings"]["xscale"] = line.split(":")[1].strip()
        elif line.startswith("yscale:"):
            current_axes["settings"]["yscale"] = line.split(":")[1].strip()
        elif line.startswith("zscale:"):
            current_axes["settings"]["zscale"] = line.split(":")[1].strip()
        elif line.startswith("save_file:"):
            save_file = line.split(":")[1].strip()

    if current_axes is not None:
        axes_instructions.append(current_axes)

    fig = plt.figure(figsize=figsize)
    axes = []

    for i in range(layout[0] * layout[1]):
        ax = fig.add_subplot(layout[0], layout[1], i + 1, projection='3d')
        axes.append(ax)

    for ax_index, ax_info in enumerate(axes_instructions):
        ax = axes[ax_index]
        for plot_info in ax_info['data']:
            if 'z_col' in plot_info:
                x, y, z = read_datafile(plot_info['file'], plot_info['x_col'], plot_info['y_col'], plot_info['z_col'])
            else:
                x, y = read_datafile(plot_info['file'], plot_info['x_col'], plot_info['y_col'])
                z = None

            method = plot_info.get('method', 'scatter')

            if method == 'plot':
                ax.plot(x, y, z, color=plot_info['color'], lw=plot_info['linewidth'], label=plot_info.get('label', ''))
            elif method == 'scatter':
                ax.scatter(x, y, z, color=plot_info['color'], s=plot_info['dot_size'], label=plot_info.get('label', ''))
            elif method == 'surface':
                X, Y = np.meshgrid(x, y)
                Z = z.reshape(X.shape)
                ax.plot_surface(X, Y, Z, color=plot_info['color'], label=plot_info.get('label', ''))
            elif method == 'wireframe':
                X, Y = np.meshgrid(x, y)
                Z = z.reshape(X.shape)
                ax.plot_wireframe(X, Y, Z, color=plot_info['color'], label=plot_info.get('label', ''))

        ax.set_title(ax_info['settings'].get('title', ''))
        ax.set_xlabel(ax_info['settings'].get('xlabel', ''))
        ax.set_ylabel(ax_info['settings'].get('ylabel', ''))
        ax.set_zlabel(ax_info['settings'].get('zlabel', ''))

        if 'xlim' in ax_info['settings']:
            ax.set_xlim(ax_info['settings']['xlim'])
        if 'ylim' in ax_info['settings']:
            ax.set_ylim(ax_info['settings']['ylim'])
        if 'zlim' in ax_info['settings']:
            ax.set_zlim(ax_info['settings']['zlim'])

        if 'xscale' in ax_info['settings']:
            ax.set_xscale(ax_info['settings']['xscale'])
        if 'yscale' in ax_info['settings']:
            ax.set_yscale(ax_info['settings']['yscale'])
        if 'zscale' in ax_info['settings']:
            ax.set_zscale(ax_info['settings']['zscale'])

        ax.legend()

    ax.view_init(elev=elevation, azim=azimuth)
    plt.tight_layout()
    plt.savefig(save_file)

if __name__ == "__main__":
    inst_file = sys.argv[1]
    plot_from_instructions(inst_file)
