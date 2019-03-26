#!/usr/bin/python
import os
import sys
import math
import argparse
# Matplotlib
import matplotlib as mpl
mpl.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{amsmath}',
    r'\usepackage{amssymb}',
    r'\usepackage{slashed}',]


def parse_args():
    """Function parses command line options"""
    parser = argparse.ArgumentParser(
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            fromfile_prefix_chars='@')
    parser.add_argument('--dir-in-data', default=["SLHA"], required=False,
                        help='Directories with SLHA input files', nargs="+")
    parser.add_argument('--grid-file', default="withNorm_mh125_char2_13_HH.txt", required=False,
                        help='Grid file with XS from SusHi')
    return parser.parse_args()

def find_data_directories(dirlist, filename_template=".slha"):
    counter=0
    data_file_list = []
    for iroot in dirlist:
        print("Searching for data in {}".format(iroot))
        for root, dirs, files in os.walk(iroot, followlinks=True):
            for ifile in files:
                if filename_template in ifile:
                    data_file_list[len(data_file_list):] = ["{}/{}".format(root, ifile)]
    if len(data_file_list)>0:
        print("Found a total of {} data files".format(len(data_file_list)))
    else:
        print("No data found in {}".format(dirlist))
        sys.exit(1)
    return data_file_list

def load_data(data_file_list, data_description, decay_list = []):
    data_dict = {}
    decay_found_flag = {}
    for idata in data_description:
        data_dict[idata[3]] = []
    for idecay in decay_list:
        data_dict[idecay[3]] = []
    # Loop aver the files
    for ifile in data_file_list:
        with open(ifile, 'r') as fh:
            filetext = fh.read()
        current_block = None
        decay_block_flg = False
        for idecay in decay_list:
            decay_found_flag[idecay[3]] = False
        for iline in filetext.splitlines():
            isplit = iline.split()
            if ((isplit[0].lower() == "block") or (isplit[0].upper() == "DECAY")):
                current_block = isplit[1]
                current_block_data = isplit[2:]
            slha_comment = isplit[-1]
            value_to_extract = list(filter(lambda x: (x[0] == current_block and x[1] == slha_comment) ,
                                                      data_description))
            if len(value_to_extract) > 0:
                if len(value_to_extract) > 1:
                    print("Error: more than matching for line {}".format(isplit))
                    print("Matching data definitions are {}".format(value_to_extract))
                data_key = value_to_extract[0][3]
                data_slha_idx = value_to_extract[0][2]
                data_dict[data_key].append(float(isplit[data_slha_idx]))
            if isplit[0].upper() == "DECAY":
                decay_block_flg = False
                # Is it a decay block we are searching for?
#                print(isplit)
                # Match if the particle decaying is specific in the DECAY entry
                decay_to_extract_list = list(filter(lambda x: (x[0] == float(isplit[1])), decay_list))
                if len(decay_to_extract_list) > 0:
                    decay_block_flg = True
            if (decay_block_flg == True):
                try:
                    # 2->2 decays, check the two decay products
                    matched_decay = list(filter(lambda x: (x[1] == float(isplit[2]) and x[2] == float(isplit[3])),
                                            decay_to_extract_list))
                except ValueError:
                    pass
                else:
                    if (len(matched_decay)>0):
 #                       print("matched decay", matched_decay)
 #                       print(isplit)
                        # The BR is always the first entry
                        data_dict[matched_decay[0][3]].append(float(isplit[0]))
                        decay_found_flag[matched_decay[0][3]] = True
        for idecay in decay_list:
            if decay_found_flag[idecay[3]] == False:
                print("BR {} set to zero".format(idecay))
                data_dict[idecay[3]].append(0.)
#        sys.exit(0)
    print("Sanity check")
    nfiles = len(data_file_list)
    print("Read a total of {} files".format(nfiles))
    for ikey in data_dict:
        print("Len data_dict[{}] = {}".format(ikey, len(data_dict[ikey])))
        assert(len(data_dict[ikey]) == nfiles)
    return data_dict

def load_grid(grid_file):
    tb_list = []
    ma_list = []
    with open(grid_file, 'r') as fh:
        for iline in fh:
            isplit = iline.split()
            tb = float(isplit[0])
            ma = float(isplit[1])
            tb_list[len(tb_list):] = [tb]
            ma_list[len(tb_list):] = [ma]
    return ma_list, tb_list

def create_fig_and_ax():
    print("Creating figure and axes")
    fig, axes = plt.subplots(1,1)
    return fig, axes

def plot_line(fig, ax):
    pass

def plot_scatter(fig, ax, x, y, colors, color_range = "dynamic", marker='s', marker_size = 0.1, rasterized = False, alpha = 1.):
    print("Plot scatterplot")
    if (color_range == "dynamic"):
        scatterplot = ax.scatter(x ,y ,c=colors, s=marker_size, marker=marker, rasterized = rasterized, alpha = alpha)
    else:
        scatterplot = ax.scatter(x ,y ,c=colors, s=marker_size, marker=marker, vmin = color_range[0], vmax = color_range[1], rasterized = rasterized, alpha = alpha)
    ax.set_xlim(min(x),max(x))
    ax.set_ylim(min(y),max(y))
#    fig.colorbar(scatterplot)
    return scatterplot

def set_labels_title(ax, suptitle, x_ax_label, y_ax_label, label_fs = 15, suptitle_fs = 15):
    ax.set_xlabel(x_ax_label, fontsize = label_fs)
    ax.set_ylabel(y_ax_label, fontsize = label_fs)
    plt.suptitle(suptitle, fontsize = suptitle_fs)

def plot_imshow(fig, ax):
    pass

def save_fig(fig, ax, outfilename = "plot.pdf", dpi=None):
    print("Saving plot to {}".format(outfilename))
    if (dpi == None):
        plt.savefig(outfilename)
    else:
        plt.savefig(outfilename, dpi = dpi)

if __name__ == "__main__":
    cmnd_args = vars(parse_args())
    ma_list, tb_list = load_grid(cmnd_args["grid_file"])
    slha_file_list = find_data_directories(cmnd_args["dir_in_data"])
    slha_file_list_ewkinos = find_data_directories(["run-03-2019"])
    data_set = load_data(slha_file_list,
                         [ ["MINPAR", "TB", 1, "tb"], ["EXTPAR", "MA0", 1, "MA"]]
    )
    data_set_ewkinos = load_data(slha_file_list_ewkinos,
                         [ ["MINPAR", "TB", 1, "tb"], ["EXTPAR", "MA0", 1, "MA"]]
    )
    fig, ax = create_fig_and_ax()
    plot_scatter(fig, ax, data_set["MA"], data_set["tb"], len(data_set["tb"])*[1], marker_size = 0.1, rasterized = True)
    plot_scatter(fig, ax, data_set_ewkinos["MA"], data_set_ewkinos["tb"], len(data_set_ewkinos["tb"])*["blue"], marker_size = 3, rasterized = True)
    plot_scatter(fig, ax, ma_list, tb_list, len(ma_list)*["red"], marker="x", marker_size = 0.05, rasterized = True, alpha=0.5)
    set_labels_title(ax, r"Light electroweakinos, grids",
                     r"$M_A$~[GeV]", r"$\tan\beta$")
    save_fig(fig, ax, outfilename = "light-ewkinos-ma-tb-plane.pdf", dpi = 300)
