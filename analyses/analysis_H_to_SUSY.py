#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import re

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
    parser.add_argument('--dir-in-data', default=["lightewkinos"], required=False,
                        help='Directories with SLHA input files', nargs="+")
    return parser.parse_args()

def find_data_directories(dirlist, filename_template=".slha"):
    counter=0
    data_dir_dict = {}
    result_re = re.compile("result-[0-9]*.txt")
    for iroot in dirlist:
        print("Searching for data in {}".format(iroot))
        for root, dirs, files in os.walk(iroot, followlinks=True):
            for ifile in files:
                if (filename_template in ifile):
                    if root in data_dir_dict:
                        data_dir_dict[root]["slha"] = "{}/{}".format(root, ifile)
                    else:
                        data_dir_dict[root] = {}
                        data_dir_dict[root]["slha"] = "{}/{}".format(root, ifile)
                if (result_re.match(ifile) != None):
                    if root in data_dir_dict:
                        if "results" in data_dir_dict[root]:
                            data_dir_dict[root]["results"].append("{}/{}".format(root, ifile))
                        else:
                            data_dir_dict[root]["results"] = []
                            data_dir_dict[root]["results"].append("{}/{}".format(root, ifile))
                    else:
                        data_dir_dict[root] = {}
                        data_dir_dict[root]["results"] = []
                        data_dir_dict[root]["results"].append("{}/{}".format(root, ifile))
    if len(data_dir_dict)>0:
        print("Found a total of {} data dir".format(len(data_dir_dict)))
    else:
        print("No data found in {}".format(dirlist))
        sys.exit(1)
    return data_dir_dict

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

def create_fig_and_ax(figsize = None):
    print("Creating figure and axes")
    if (figsize == None):
        fig, axes = plt.subplots(1,1)
    else:
        fig, axes = plt.subplots(1,1, figsize = figsize)
    return fig, axes

def plot_line(fig, ax, x, y, alpha = 1, linewidth = 1, linestyle = 'solid'):
    ax.plot(y, alpha = alpha, linewidth = linewidth, linestyle = linestyle)

def plot_scatter(fig, ax, x, y, colors, color_range = "dynamic", marker='s', marker_size = 0.1, rasterized = False, alpha = 1., cmap = None):
    print("Plot scatterplot")
    if (color_range == "dynamic"):
        scatterplot = ax.scatter(x ,y ,c=colors, s=marker_size, marker=marker, rasterized = rasterized, alpha = alpha, cmap = cmap)
    else:
        scatterplot = ax.scatter(x ,y ,c=colors, s=marker_size, marker=marker, vmin = color_range[0], vmax = color_range[1], rasterized = rasterized, alpha = alpha, cmap = cmap)
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

def close_fig(fig):
    plt.close(fig)

def get_r_vs_iteration_data(data_dir_dict):
    r_all_points = []
    for idir in sorted(data_dir_dict):
        sorted_file_list = sorted(data_dir_dict[idir]["results"],
                                  key = lambda x: int(os.path.basename(x).replace(".txt", "").replace("result-", "")))
        r_values = []
        analysis_list = []
        SR_list = []
        for ires in sorted_file_list:
            with open(ires, 'r') as fh:
                for iline in fh:
                    isplit = iline.split()
                    if ((isplit[0] == "Result") and (isplit[1] == "for")):
                        ir = float(isplit[-1])
                        r_values[len(r_values):] = [ir]
                    if (isplit[0] == "Analysis:"):
                        analysis_list[len(analysis_list):] = [isplit[-1]]
                    if (isplit[0] == "SR:"):
                        SR_list[len(SR_list):] = [isplit[-1]]
        data_dir_dict[idir]["r-value"] = r_values[-1]
        data_dir_dict[idir]["analysis"] = analysis_list[-1]
        data_dir_dict[idir]["SR"] = SR_list[-1]
        r_all_points[len(r_all_points):] = [r_values]
    return r_all_points

def load_signal_regions(data_dir_dict, fn):
    best_signal_dict = {}
    dict_analysis = {}
    for idir in data_dir_dict:
        best_signal_dict[idir] = {}
        best_signal_fn = "{}/evaluation/{}".format(idir, fn)
        with open(best_signal_fn, 'r') as fh:
            for iline in fh:
                isplit = iline.split()
                ## This is always the first line
                if (isplit[0] == "analysis"):
                    indx_r = isplit.index("robscons")
                else:
                    analysis = isplit[0]
                    sr = isplit[1]
                    if not (analysis in dict_analysis):
                        dict_analysis[analysis] = set(["{}".format(sr)])
                    else:
                        dict_analysis[analysis].add("{}".format(sr))
                    if not (analysis in best_signal_dict[idir]):
                        best_signal_dict[idir][analysis] = {}
                    best_signal_dict[idir][analysis][sr] = isplit[indx_r]
    return best_signal_dict, dict_analysis


def get_max_r_analysis(analysis, total_results_dict):
    # NOTE: SORTED for consistency
    r_values = []
    for idir in sorted(total_results_dict.keys()):
        if not (analysis in total_results_dict[idir]):
            r_values[len(r_values):] = [0]
        else:
            all_r = [float(total_results_dict[idir][analysis][sr]) for sr in total_results_dict[idir][analysis]]
            r_values[len(r_values):] = [max(all_r)]
    return r_values

def get_max_r_analysis_sr(analysis, sr, total_results_dict):
    # NOTE: SORTED for consistency
    r_values = []
    for idir in sorted(total_results_dict.keys()):
        if not (analysis in total_results_dict[idir]):
            r_values[len(r_values):] = [0]
        else:
            if not (sr in total_results_dict[idir][analysis]):
                r_values[len(r_values):] = [0]
            else:
                r_values[len(r_values):] = [float(total_results_dict[idir][analysis][sr])]
    return r_values


def load_process_detail_analysis_sr(data_dir_dict):
    process_list_fn = ["GG_H_MSSM_processResults.txt", "GG_A_MSSM_processResults.txt", "BB_H_MSSM_processResults.txt", "BB_A_MSSM_processResults.txt",
                    "other_proc_H_A_MSSM_processResults.txt",
                    "directewkinos_processResults.txt"]

    process_name = ["ggH", "ggA", "bbH", "bbA", "otherHA", "directEWKinos"]

    dict_analysis = {}

    dict_process_detail = {}
    for idir in sorted(data_dir_dict):
        dict_process_detail[idir] = {}
        for iproc in process_name:
            dict_process_detail[idir][iproc] = {}
        for i, iprocess in enumerate(process_list_fn):
            process_key = process_name[i]
            process_fn = "{}/evaluation/{}".format(idir, iprocess)
            with open(process_fn, 'r') as fh:
                for iline in fh:
                    isplit = iline.split()
                    if isplit[0] == "analysis":
                        indx_norm_events = isplit.index("signalnormevents")
                    else:
                        ianalysis = isplit[0]
                        isr = isplit[1]
                        inorm_events = isplit[indx_norm_events]
                        if not (ianalysis in dict_process_detail[idir][process_key]):
                            dict_process_detail[idir][process_key][ianalysis] = {}
                            dict_analysis[ianalysis] = set()
                        dict_analysis[ianalysis].add(isr)
                        dict_process_detail[idir][process_key][ianalysis][isr] = inorm_events

    return dict_process_detail, dict_analysis

def get_dominant_process_sr_analysis(analysis, sr, process_detail_dict):
    process_name = ["ggH", "ggA", "bbH", "bbA", "otherHA", "directEWKinos"]
    color_process = ["red", "blue", "green", "magenta", "pink", "black"]
    # NOTE: SORTED for consistency
    dominant_proc = []
    for idir in sorted(process_detail_dict.keys()):
        iproc_list = []
        for iproc in process_name:
            if not (analysis in process_detail_dict[idir][iproc]):
                iproc_list[len(iproc_list):] = [0]
            else:
                if not (sr in process_detail_dict[idir][iproc][analysis]):
                    iproc_list[len(iproc_list):] = [0]
                else:
                    iproc_list[len(iproc_list):] = [float(process_detail_dict[idir][iproc][analysis][sr])]
        max_events = max(iproc_list)
        if (max_events > 0):
            idx_max = iproc_list.index(max_events)
            max_proc_key = process_name[idx_max]
            icolor = mpl.colors.to_rgba(color_process[idx_max])
            # change alpha
            icolor = (icolor[0], icolor[1], icolor[2], max_events/sum(iproc_list))
        else:
            icolor = "white"
        dominant_proc[len(dominant_proc):] = [icolor]
    return dominant_proc

def get_dominant_process_analysis(analysis, process_detail_dict, best_analysis_dict):
    process_name = ["ggH", "ggA", "bbH", "bbA", "otherHA", "directEWKinos"]
    color_process = ["red", "blue", "green", "magenta", "pink", "black"]
    # NOTE: SORTED for consistency
    dominant_proc = []
    for idir in sorted(process_detail_dict.keys()):
        iproc_list = []
        if not (analysis in best_analysis_dict[idir]):
            max_events = 0
        else:
            best_sr_analysis_list = list(best_analysis_dict[idir][analysis].keys())
            if (len(best_sr_analysis_list) > 1):
                print("Error: incosistent best_signal_regions.txt")
                sys.exit(1)
            best_sr_analysis = best_sr_analysis_list[0]
            for iproc in process_name:
                if not (analysis in process_detail_dict[idir][iproc]):
                    iproc_list[len(iproc_list):] = [0]
                else:
                    if not (best_sr_analysis in process_detail_dict[idir][iproc][analysis]):
                        iproc_list[len(iproc_list):] = [0]
                    else:
                        iproc_list[len(iproc_list):] = [float(process_detail_dict[idir][iproc][analysis][best_sr_analysis])]
            max_events = max(iproc_list)
        if (max_events > 0):
            idx_max = iproc_list.index(max_events)
            max_proc_key = process_name[idx_max]
            icolor = mpl.colors.to_rgba(color_process[idx_max])
            # change alpha
            icolor = (icolor[0], icolor[1], icolor[2], max_events/sum(iproc_list))
        else:
            icolor = (1,1,1,1)
        dominant_proc[len(dominant_proc):] = [icolor]
    return dominant_proc

def get_dominant_proc(data_dir_dict, process_detail_dict):
    process_name = ["ggH", "ggA", "bbH", "bbA", "otherHA", "directEWKinos"]
    color_process = ["red", "blue", "green", "magenta", "pink", "black"]
    dominant_proc_list = []
    # NOTE: SORTED for consistency
    for idir in sorted(data_dir_dict):
        analysis = data_dir_dict[idir]["analysis"]
        best_sr_analysis = data_dir_dict[idir]["SR"]
        iproc_list = []
        # Loop all over the processes
#        print(idir)
        for iproc in process_name:
            iproc_list[len(iproc_list):] = [float(process_detail_dict[idir][iproc][analysis][best_sr_analysis])]
#            print(iproc, analysis)
#            print(process_detail_dict[idir][iproc][analysis])
#            sys.exit(1)
        # Find the largest one
#        print(idir, analysis, best_sr_analysis, iproc_list, data_dir_dict[idir]["r-value"])
        max_events = max(iproc_list)
        idx_max = iproc_list.index(max_events)
        max_proc_key = process_name[idx_max]
        icolor = mpl.colors.to_rgba(color_process[idx_max])
        # change alpha
        icolor = (icolor[0], icolor[1], icolor[2], max_events/sum(iproc_list))
        dominant_proc_list[len(dominant_proc_list):] = [icolor]
    return dominant_proc_list

if __name__ == "__main__":
    cmnd_args = vars(parse_args())
    data_dir_dict = find_data_directories(cmnd_args["dir_in_data"])
    # NOTE: SORTED for consistency
    slha_file_list = [ data_dir_dict[idir]["slha"] for idir in sorted(data_dir_dict)]
    best_signal_regions_dict, best_analysis_dict = load_signal_regions(data_dir_dict, "best_signal_regions.txt")
    total_results_dict, total_analysis_dict = load_signal_regions(data_dir_dict, "total_results.txt")
    process_detail_dict, process_analysis_dict = load_process_detail_analysis_sr(data_dir_dict)
    # r vs iteration plots, but also adds to the dict the r-value, and the best analysis-SR
    r_vs_iteration_list = get_r_vs_iteration_data(data_dir_dict)

    # best-r in the MA-tanb plane, split by analysis
    ma_tb_data_set = load_data(slha_file_list,
                         [ ["MINPAR", "TB", 1, "tb"], ["EXTPAR", "MA0", 1, "MA"]]
    )
    ma_array = np.array(ma_tb_data_set["MA"])
    tb_array = np.array(ma_tb_data_set["tb"])
    cmap = plt.get_cmap("coolwarm")

    global_proc_flag = False
    # Plot the mA-tanb plane, showing for each point the dominant process
    if (global_proc_flag == True):
        global_proc = get_dominant_proc(data_dir_dict, process_detail_dict)
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_tb_data_set["MA"], ma_tb_data_set["tb"], global_proc, marker_size = 10)
        save_fig(fig, ax, outfilename = "light-ewkinos-ma-tanb-dominant-proc.pdf")
        close_fig(fig)

    global_analysis_flag = True
    if (global_analysis_flag == True):
        global_proc = get_dominant_analysis(data_dir_dict)
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_tb_data_set["MA"], ma_tb_data_set["tb"], global_proc, marker_size = 10)
        save_fig(fig, ax, outfilename = "light-ewkinos-ma-tanb-dominant-analysis.pdf")
        close_fig(fig)
    # Plot the mA-tanb plane, showing for each point the dominant process

    r_study_flag = False
    ## Plot r vs iteration for a random selection of 100 points
    if (r_study_flag == True):
        # the r-value item is added by the call above to get_r_vs_iteration_data
        r_values = [ data_dir_dict[idir]["r-value"] for idir in sorted(data_dir_dict)]
        fig, ax = create_fig_and_ax()
        for iridx in set(np.random.randint(0, len(r_vs_iteration_list), 100)):
            plot_line(fig, ax, None, r_vs_iteration_list[iridx], alpha = 0.8)
        set_labels_title(ax, "r vs event number", "iteration", "r")
        save_fig(fig, ax, outfilename = "light-ewkinos-r-evolution.pdf")
        close_fig(fig)

        # global best-R in the MA-tanb plane
        fig, ax = create_fig_and_ax()
        cmap = plt.get_cmap("coolwarm")
        ma_array = np.array(ma_tb_data_set["MA"])
        tb_array = np.array(ma_tb_data_set["tb"])
        r_values_array = np.array(r_values)
        condition_r = np.where(r_values_array > 1)
        plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [-3,5], cmap = cmap)
        condition_r = np.where(r_values_array <= 1)
        plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [0,2], cmap = cmap)
        save_fig(fig, ax, outfilename = "light-ewkinos-r-ma-tb-plane.pdf")
        close_fig(fig)

#### Detailed per-analysis/per-analysis-SR planes

    analysis_flag = False
    ## Plot a mA-tanb plane for each analysis (of which the best SR is chosen), showing the r-value with a color-coded scale
    if (analysis_flag == True):
        analysis_plot_dir = "per-analysis"
        if not os.path.isdir(analysis_plot_dir):
            os.makedirs(analysis_plot_dir)

        for ianalysis in sorted(total_analysis_dict.keys()):
            print("Plotting analysis {}".format(ianalysis))
            r_values = np.array(get_max_r_analysis(ianalysis, total_results_dict))
            r_values_array = np.array(r_values)
            fig, ax = create_fig_and_ax()
            set_labels_title(ax, r"{}".format(ianalysis.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
            condition_r = np.where(r_values_array > 1)
            if len(r_values_array[condition_r]) > 0:
                plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [-3,5], cmap = cmap)
            condition_r = np.where(r_values_array <= 1)
            if len(r_values_array[condition_r]) > 0:
                plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [0,2], cmap = cmap)
            save_fig(fig, ax, outfilename = "{}/light-ewkinos-r-ma-tb-plane-{}.pdf".format(analysis_plot_dir, ianalysis))
            close_fig(fig)

    sr_analysis_flag = False
    ## Plot a mA-tanb plane for each analysis and SR, showing the r-value with a color-coded scale
    if (sr_analysis_flag == True):
        for ianalysis in sorted(total_analysis_dict.keys()):
            for isr in sorted(total_analysis_dict[ianalysis]):
                r_values = np.array(get_max_r_analysis_sr(ianalysis, isr, total_results_dict))
                r_values_array = np.array(r_values)
                fig, ax = create_fig_and_ax()
                set_labels_title(ax, r"{}~~{}".format(ianalysis.replace("_","\_"), isr.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                condition_r = np.where(r_values_array > 1)
                if len(r_values_array[condition_r]) > 0:
                    plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [-3,5], cmap = cmap)
                condition_r = np.where(r_values_array <= 1)
                if len(r_values_array[condition_r]) > 0:
                    plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [0,2], cmap = cmap)

                sr_analysis_plot_dir = "per-analysis-sr/{}".format(ianalysis)
                if not os.path.isdir(sr_analysis_plot_dir):
                    os.makedirs(sr_analysis_plot_dir)
                save_fig(fig, ax, outfilename = "{}/light-ewkinos-r-ma-tb-plane-{}-sr-{}.pdf".format(sr_analysis_plot_dir, ianalysis, isr))
                close_fig(fig)

    dominant_process_analysis_sr_flag = False
    ## Plot a mA-tanb plane for each analysis and SR, showing the dominant process for each point using a color coding
    if (dominant_process_analysis_sr_flag == True):
        for ianalysis in process_analysis_dict:
            for isr in process_analysis_dict[ianalysis]:
                dominant_process_list = get_dominant_process_sr_analysis(ianalysis, process_detail_dict)

                fig, ax = create_fig_and_ax()

                set_labels_title(ax, r"dominant proc.~{}~~{}".format(ianalysis.replace("_","\_"), isr.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                plot_scatter(fig, ax, ma_tb_data_set["MA"], ma_tb_data_set["tb"], dominant_process_list, marker_size = 10)

                sr_analysis_plot_dir = "dominant-process-per-analysis-sr/{}".format(ianalysis)
                if not os.path.isdir(sr_analysis_plot_dir):
                    os.makedirs(sr_analysis_plot_dir)
                save_fig(fig, ax, outfilename = "{}/light-ewkinos-dominant-process-ma-tb-plane-{}-sr-{}.pdf".format(sr_analysis_plot_dir, ianalysis, isr))
                close_fig(fig)

    dominant_process_analysis_flag = False
    ## Plot a mA-tanb plane for each analysis (chosing the best SR on a point-by-point basis), showing the dominant process for each point using a color coding
    if (dominant_process_analysis_flag == True):
        for ianalysis in process_analysis_dict:
            print("Plotting dominant process for analysis {}".format(ianalysis))
            # atlas_1708_07875
            if (ianalysis != "atlas_1708_07875"):
                continue
            dominant_process_list = get_dominant_process_analysis(ianalysis, process_detail_dict, best_signal_regions_dict)

            fig, ax = create_fig_and_ax()

            set_labels_title(ax, r"dominant proc.~{}".format(ianalysis.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
            print(dominant_process_list)
            plot_scatter(fig, ax, ma_tb_data_set["MA"], ma_tb_data_set["tb"], dominant_process_list, marker_size = 10)

            proc_analysis_plot_dir = "dominant-process-per-analysis/{}".format(ianalysis)
            if not os.path.isdir(proc_analysis_plot_dir):
                os.makedirs(proc_analysis_plot_dir)
            save_fig(fig, ax, outfilename = "{}/light-ewkinos-dominant-process-ma-tb-plane-{}.pdf".format(proc_analysis_plot_dir, ianalysis))
            close_fig(fig)
