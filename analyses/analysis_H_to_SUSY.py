#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import re
import itertools
from plotting import *
from slha_utils import *
from CM_analysis import *
from scenario_spectrum import *

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
    # The format for the standard blocks is ["BLOCK", "LAST COLUMN", field of data, key used in the dictionary]
    ma_tb_data_set = load_data(slha_file_list,
                         [ ["MINPAR", "TB", 1, "tb"], ["EXTPAR", "MA0", 1, "MA"],
                           ["MASS", "MHH", 1, "MH"], ["MASS", "Mh0", 1, "Mh"],
                           ["MASS", "MNeu(1)", 1, "mneu1"], ["MASS", "MNeu(2)", 1, "mneu2"],
                           ["MASS", "MNeu(3)", 1, "mneu3"], ["MASS", "MNeu(4)", 1, "mneu4"],
                           ["MASS", "MCha(1)", 1, "mcha1"], ["MASS", "MCha(2)", 1, "mcha2"],
                           ["EXTPAR", "M1", 1, "M1"], ["EXTPAR", "M2", 1, "M2"], ["EXTPAR", "MUE", 1, "mu"],
                         ]
    )
    ma_array = np.array(ma_tb_data_set["MA"])
    tb_array = np.array(ma_tb_data_set["tb"])
    M1_array = np.array(ma_tb_data_set["M1"])
    #scenario_selection = np.where(M1_array < 100)
    scenario_selection = np.where((M1_array > 100) & (M1_array < 200))
    cmap = plt.get_cmap("coolwarm")

    ## Full dataset for the mass planes and the decays
    mass_planes_flag = True
    if (mass_planes_flag == True):
        spectrum_analysis()

    global_proc_flag = False
    # Plot the mA-tanb plane, showing for each point the dominant process
    if (global_proc_flag == True):
        global_proc = get_dominant_proc(data_dir_dict, process_detail_dict)
        global_proc_arr = np.array(global_proc)
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], global_proc_arr[scenario_selection], marker_size = 10)
        save_fig(fig, ax, outfilename = "light-ewkinos-ma-tanb-dominant-proc.pdf")
        close_fig(fig)

    global_analysis_flag = False
    # Plot the mA-tanb plane, showing for each point the dominant process
    if (global_analysis_flag == True):
        global_proc = get_dominant_analysis(data_dir_dict)
        global_proc_arr = np.array(global_proc)
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], global_proc_arr[scenario_selection], marker_size = 10)
        save_fig(fig, ax, outfilename = "light-ewkinos-ma-tanb-dominant-analysis.pdf")
        close_fig(fig)

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
        r_values_array = np.array(r_values)
        condition_r = np.where((r_values_array > 1) & (M1_array > 100) & (M1_array < 100))
        cmap = plt.get_cmap("Reds")
        plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [-3,1], cmap = cmap)
        condition_r = np.where((r_values_array <= 1) & (M1_array > 100) & (M1_array < 100))
        cmap = plt.get_cmap("winter")
        plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [1,2], cmap = cmap)
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
            condition_r = np.where((r_values_array > 1) & (M1_array > 100) & (M1_array < 100))
            if len(r_values_array[condition_r]) > 0:
                plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [-3,5], cmap = cmap)
            condition_r = np.where(r_values_array <= 1 & (M1_array > 100) & (M1_array < 100))
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
                condition_r = np.where(r_values_array > 1 & (M1_array > 100) & (M1_array < 100))
                if len(r_values_array[condition_r]) > 0:
                    plot_scatter(fig, ax, ma_array[condition_r], tb_array[condition_r], r_values_array[condition_r], marker_size = 10, color_range = [-3,5], cmap = cmap)
                condition_r = np.where(r_values_array <= 1 & (M1_array > 100) & (M1_array < 100))
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
                dominant_process_arr = np.array(dominant_process_list)
                fig, ax = create_fig_and_ax()

                set_labels_title(ax, r"dominant proc.~{}~~{}".format(ianalysis.replace("_","\_"), isr.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], dominant_process_arr[scenario_selection], marker_size = 10)

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
            dominant_process_arr = np.array(dominant_process_list)
            fig, ax = create_fig_and_ax()

            set_labels_title(ax, r"dominant proc.~{}".format(ianalysis.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
            print(dominant_process_list)
            plot_scatter(fig, ax, ma_array[scenario_selection], ma_tb_data_set[scenario_selection], dominant_process_arr[scenario_selection], marker_size = 10)

            proc_analysis_plot_dir = "dominant-process-per-analysis/{}".format(ianalysis)
            if not os.path.isdir(proc_analysis_plot_dir):
                os.makedirs(proc_analysis_plot_dir)
            save_fig(fig, ax, outfilename = "{}/light-ewkinos-dominant-process-ma-tb-plane-{}.pdf".format(proc_analysis_plot_dir, ianalysis))
            close_fig(fig)
