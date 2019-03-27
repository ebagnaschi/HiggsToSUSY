# Matplotlib
import numpy as np
import matplotlib as mpl
import matplotlib.mlab as mlab
mpl.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}',
    r'\usepackage{slashed}',]
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def create_fig_and_ax(figsize = None):
    print("Creating figure and axes")
    if (figsize == None):
        fig, axes = plt.subplots(1,1)
    else:
        fig, axes = plt.subplots(1,1, figsize = figsize)
    return fig, axes

def plot_line_datasets(fig, ax, datasets, alpha = 1, linewidth = 1, linestyle = 'solid'):
    for idata in datasets:
        plot_line(fig, ax, idata[0], idata[1], alpha = alpha, linewidth = linewidth, linestyle = linestyle)

def plot_line(fig, ax, x = None, y = None, alpha = 1, linewidth = 1, linestyle = 'solid'):
    if ((x == None) and (y != None)):
        ax.plot(y, alpha = alpha, linewidth = linewidth, linestyle = linestyle)
    else:
        ax.plot(x, y, alpha = alpha, linewidth = linewidth, linestyle = linestyle)

def plot_scatter(fig, ax, x, y, colors, color_range = "dynamic", marker='s', marker_size = 0.1, rasterized = False, alpha = 1., cmap = None, cbar_flag = False, plot_lims = None, cbar_label = ""):
    print("Plot scatterplot")
    if cmap == None:
        cmap_obj = mpl.cm.get_cmap("viridis")
    else:
        cmap_obj = mpl.cm.get_cmap(cmap)
    if (color_range == "dynamic"):
        scatterplot = ax.scatter(x ,y ,c=colors, s=marker_size, marker=marker, rasterized = rasterized, alpha = alpha, cmap = cmap_obj)
    else:
        scatterplot = ax.scatter(x ,y ,c=colors, s=marker_size, marker=marker, vmin = color_range[0], vmax = color_range[1], rasterized = rasterized, alpha = alpha, cmap = cmap_obj)
    if (plot_lims == None):
        ax.set_xlim(min(x),max(x))
        ax.set_ylim(min(y),max(y))
    else:
        ax.set_xlim(*plot_lims[0])
        ax.set_ylim(*plot_lims[1])
    if (cbar_flag == True):
        cb = fig.colorbar(scatterplot)
        cb.set_label(cbar_label, rotation=-90, labelpad = 15)
    return scatterplot

def set_labels_title(ax, suptitle, x_ax_label, y_ax_label, tick_fs = 13, label_fs = 15, suptitle_fs = 15, y_sup = 0.95):
    ax.set_xlabel(x_ax_label, fontsize = label_fs)
    ax.set_ylabel(y_ax_label, fontsize = label_fs)
    plt.suptitle(suptitle, fontsize = suptitle_fs, y = y_sup)
    ax.tick_params(axis='both', which='major', labelsize = tick_fs, bottom = True, top = True, left = True, right = True)
    ax.tick_params(axis='both', which='minor', labelsize = tick_fs, bottom = True, top = True, left = True, right = True)

def set_ranges(fig, ax, xlim, ylim):
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def set_x_major_ticks(fig, ax, interval):
    majorLocator = MultipleLocator(interval)
    ax.xaxis.set_major_locator(majorLocator)

def set_y_major_ticks(fig, ax, interval):
    majorLocator = MultipleLocator(interval)
    ax.yaxis.set_major_locator(majorLocator)

def set_x_minor_ticks(fig, ax, interval):
    minorLocator = MultipleLocator(interval)
    ax.xaxis.set_minor_locator(minorLocator)

def set_y_minor_ticks(fig, ax, interval):
    minorLocator = MultipleLocator(interval)
    ax.yaxis.set_minor_locator(minorLocator)


def plot_tricontourf(fig, ax, x, y, colors, color_range = "dynamic", level_contourf = None,
                  cbar_options = {}, cmap = None, vrange = None):
    print("Plotting contourf")
    if cmap == None:
        cmap_obj = None
    else:
        cmap_obj = mpl.cm.get_cmap(cmap)
    if (level_contourf == None):
        cs = ax.tricontourf(x, y, colors, cmap = cmap_obj)
    else:
        cs = ax.tricontourf(x, y, colors, level_contourf, cmap = cmap_obj)
    cbar = fig.colorbar(cs)
    if "label" in cbar_options:
        cbar.set_label(cbar_options["label"], fontsize=15)
    cbar.ax.tick_params(labelsize=14)
    return cs

def plot_tricontour(fig, ax, x, y, colors, color_range = "dynamic", level_contour = None,
                  fmt = '%2.1f', linestyles = "solid", cmap = None, legend_label = None, vrange = None):
    print("Plotting tricontour")
    if cmap == None:
        cmap_obj = mpl.cm.get_cmap("viridis")
    else:
        cmap_obj = mpl.cm.get_cmap(cmap)
    print(colors)
    if (level_contour == None):
        co = ax.tricontour(x, y, colors, linestyles = linestyles, cmap = cmap_obj)
    else:
        co = ax.tricontour(x, y, colors, level_contour, linestyles = linestyles, cmap = cmap_obj)
#    ax.clabel(co, fmt=fmt, colors='black', fontsize=14)
    ax.clabel(co, fmt=fmt,  fontsize=14)
    if legend_label != None:
        line_custom_legend = Line2D([0,1], [0,1], color=cmap_obj(0.), lw=2, linestyle=linestyles)
        line_custom_legend.set_dashes([0,2])
        proxy_artist = [line_custom_legend]
        ax.legend(proxy_artist, [legend_label], loc = 3, fontsize = 12)
    return co

def plot_contour(fig, ax, x, y, z, color_range = "dynamic", level_contour = None,
                  fmt = '%2.1f', linestyles = "solid", cmap = None, colors = None, legend_label = None, vrange = None,
                 interp_grid = True, alpha = 1, clabel_color = None):
    print("Plotting contour")
    if cmap == None:
        cmap_obj = None
    else:
        cmap_obj = mpl.cm.get_cmap(cmap)

    if (interp_grid == True):
        xi = np.linspace(min(x), max(x))
        yi = np.linspace(min(y), max(y))
        zi = mlab.griddata(x, y, z, xi, yi, interp='nn')
    else:
        xi = x
        yi = y
        zi = z

    if (level_contour == None):
        co = ax.contour(xi, yi, zi, linestyles = linestyles, cmap = cmap_obj, colors = colors, alpha = alpha)
    else:
        co = ax.contour(xi, yi, zi, level_contour, linestyles = linestyles, cmap = cmap_obj, colors = colors, alpha = alpha)
    ax.clabel(co, fmt=fmt, colors = clabel_color, fontsize=14)
    if legend_label != None:
        line_custom_legend = Line2D([0,1], [0,1], color=cmap_obj(0.), lw=2, linestyle=linestyles)
        line_custom_legend.set_dashes([0,2])
        proxy_artist = [line_custom_legend]
        ax.legend(proxy_artist, [legend_label], loc = 3, fontsize = 12)
    return co

def plot_contourf(fig, ax, x, y, colors, color_range = "dynamic", level_contour = None,
                  fmt = '%2.1f', linestyles = "solid", cmap = None, legend_label = None, vrange = None,
                  interp_grid = True):
    print("Plotting contour")
    if cmap == None:
        cmap_obj = mpl.cm.get_cmap("viridis")
    else:
        cmap_obj = mpl.cm.get_cmap(cmap)

    if (interp_grid == True):
        xi = np.linspace(min(x), max(x))
        yi = np.linspace(min(y), max(y))
        zi = mlab.griddata(x, y, colors, xi, yi, interp='nn')
    else:
        xi = x
        yi = y
        zi = colors

    if (level_contour == None):
        co = ax.contourf(xi, yi, zi, linestyles = linestyles, cmap = cmap_obj)
    else:
        co = ax.contourf(xi, yi, zi, level_contour, linestyles = linestyles, cmap = cmap_obj)
    return co


def plot_imshow(fig, ax, x, y, z, vrange = [None, None], interp_grid = True, cbar_flag = False, cbar_label = None):
    if (interp_grid == True):
        xi = np.linspace(min(x), max(x))
        yi = np.linspace(min(y), max(y))
        zi = mlab.griddata(x, y, colors, xi, yi, interp='nn')
    else:
        xi = x
        yi = y
        zi = colors

    imshow = ax.imshow(xi,yi,zi, vmin = vrange[0], vmax = vrange[1])
    if (cbar_flag == True):
        cb = fig.colorbar(imshow)
        cb.set_label(cbar_label, rotation=-90, labelpad = 15)
    return imshow

def plot_pcolormesh(fig, ax, x, y, z, cmap = None, vrange = [None, None], interp_grid = True, cbar_flag = False, cbar_label = None):
    if cmap == None:
        cmap_obj = mpl.cm.get_cmap("viridis")
    else:
        cmap_obj = mpl.cm.get_cmap(cmap)
    if (interp_grid == True):
        xi = np.linspace(min(x), max(x))
        yi = np.linspace(min(y), max(y))
        zi = mlab.griddata(x, y, z, xi, yi, interp='nn')
    else:
        xi = x
        yi = y
        zi = colors
    pcol = ax.pcolormesh(xi,yi,zi, cmap = cmap_obj, vmin = vrange[0], vmax = vrange[1])
    if (cbar_flag == True):
        cb = fig.colorbar(pcol)
        cb.set_label(cbar_label, rotation=-90, labelpad = 15)
    return pcol


def save_fig(fig, ax, outfilename = "plot.pdf", dpi=None):
    print("Saving plot to {}".format(outfilename))
    if (dpi == None):
        plt.savefig(outfilename)
    else:
        plt.savefig(outfilename, dpi = dpi)

def close_fig(fig):
    plt.close(fig)
