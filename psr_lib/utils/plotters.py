from matplotlib import pyplot as plt
import numpy as np

def polar_plot(func, xs=None, phis=None, use_shading=False, file_name=None, label=None, contour_level=None, cmap='viridis'):
    fig = plt.figure(figsize=[5,5])
    ax = fig.add_axes([0.1,0.1,0.8,0.8], polar=True)
    ax.grid(False)
    if callable(func):
        if xs is None:
            xs = np.linspace(0.01, 0.99, 100)
        if phis is None:
            phis = np.linspace(0, 2*np.pi, 100, endpoint=False)
        xs_grid, phis_grid = np.meshgrid(xs, phis, indexing='ij')
        funcs = func(xs_grid, phis_grid)
    else: 
        if xs is None or phis is None:
            raise ValueError("xs ans phis must be passed if func is not callable")
        funcs = func
    if use_shading == True:
        bar = ax.pcolormesh(phis, xs, funcs, edgecolors='face', shading='gouraud', cmap=cmap)
    else:
        bar = ax.pcolormesh(phis, xs, funcs, edgecolors='face', cmap=cmap)
    if contour_level != None:
        ax.contour(phis, xs, funcs, contour_level, colors=['black'])
    if label is not None:
        fig.colorbar(bar, ax=ax, label=label)
    else: 
        fig.colorbar(bar, ax=ax)
    if file_name is not None:
        fig.savefig(file_name + '.png', dpi=400, bbox_inches='tight')
    fig.show()
    return fig, ax