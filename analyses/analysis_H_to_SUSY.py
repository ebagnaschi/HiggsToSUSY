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
from scenario_ewkino_xs import *
from extract_xs_ewkino import *
import pickle

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
    result_MG_re = re.compile("result-MG-[0-9]*.html")
    for iroot in dirlist:
        print("Searching for data in {}".format(iroot))
        for root, dirs, files in os.walk(iroot, followlinks=True):
#            if (counter >= 20):
#                break
            for ifile in files:
                if (filename_template in ifile):
                    if root in data_dir_dict:
#                        counter = counter+1
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
                if (result_MG_re.match(ifile) != None):
                    if root in data_dir_dict:
                        if "resultsMG" in data_dir_dict[root]:
                            data_dir_dict[root]["resultsMG"].append("{}/{}".format(root, ifile))
                        else:
                            data_dir_dict[root]["resultsMG"] = []
                            data_dir_dict[root]["resultsMG"].append("{}/{}".format(root, ifile))
                    else:
                        data_dir_dict[root] = {}
                        data_dir_dict[root]["resultsMG"] = []
                        data_dir_dict[root]["resultsMG"].append("{}/{}".format(root, ifile))
    if len(data_dir_dict)>0:
        print("Found a total of {} data dir".format(len(data_dir_dict)))
    else:
        print("No data found in {}".format(dirlist))
        sys.exit(1)
    return data_dir_dict

if __name__ == "__main__":
    cmnd_args = vars(parse_args())

    res_info_pkl = "res_loc.pkl"
    if os.path.isfile(res_info_pkl):
        print("Loading result location from {}".format(res_info_pkl))
        with open(res_info_pkl, 'rb') as fh:
            data_dir_dict, slha_file_list = pickle.load(fh)
    else:
        data_dir_dict = find_data_directories(cmnd_args["dir_in_data"])
        # NOTE: SORTED for consistency
        slha_file_list = [ data_dir_dict[idir]["slha"] for idir in sorted(data_dir_dict)]
        with open(res_info_pkl, 'wb') as fh:
            pickle.dump([data_dir_dict, slha_file_list], fh)

    ewkino_pkl_file = "ewkino_xs.pkl"
    if os.path.isfile(ewkino_pkl_file):
        print("Loading ewkino xs from pickl'ed file")
        with open(ewkino_pkl_file, 'rb') as fh:
            ewkino_xs = pickle.load(fh)
    else:
        ewkino_xs = extract_ewkino_xs(data_dir_dict)
        with open(ewkino_pkl_file, 'wb') as fh:
            pickle.dump(ewkino_xs, fh)

    cm_pkl_file = "CM-result.pkl"
    if os.path.isfile(cm_pkl_file):
        print("Loading CM results from {}".format(cm_pkl_file))
        with open(cm_pkl_file, 'rb') as fh:
            best_signal_regions_dict, best_analysis_dict, total_results_dict, total_analysis_dict, process_detail_dict, process_analysis_dict, r_vs_iteration_list, r_values, data_dir_dict = pickle.load(fh)
    else:
        # best-r in the MA-tanb plane, split by analysis
        best_signal_regions_dict, best_analysis_dict = load_signal_regions(data_dir_dict, "best_signal_regions.txt")
        total_results_dict, total_analysis_dict = load_signal_regions(data_dir_dict, "total_results.txt")
        process_detail_dict, process_analysis_dict = load_process_detail_analysis_sr(data_dir_dict)
        # r vs iteration plots, but also adds to the dict the r-value, and the best analysis-SR. It loads the data from results-*.txt.
        r_vs_iteration_list = get_r_vs_iteration_data(data_dir_dict)
        # the r-value item is added by the call above to get_r_vs_iteration_data
        r_values = [ data_dir_dict[idir]["r-value"] for idir in sorted(data_dir_dict)]
        pickle_list = [best_signal_regions_dict, best_analysis_dict,
                       total_results_dict, total_analysis_dict,
                       process_detail_dict, process_analysis_dict,
                       r_vs_iteration_list,
                       r_values, data_dir_dict]
        with open(cm_pkl_file, 'wb') as fh:
            pickle.dump(pickle_list, fh)

    r_vs_iteration_arr = np.array(r_vs_iteration_list)
    r_values_array = np.array(r_values)

    # The format for the standard blocks is ["BLOCK", "LAST COLUMN", field of data, key used in the dictionary]
    ma_tb_data_set_pkl_fn = "ma_tb_data_set.pkl"
    if os.path.isfile(ma_tb_data_set_pkl_fn):
        print("Loading data of the main set SLHA files from pkl file {}".format(ma_tb_data_set_pkl_fn))
        with open(ma_tb_data_set_pkl_fn, 'rb') as fh:
            ma_tb_data_set = pickle.load(fh)
    else:
        ma_tb_data_set = load_data(slha_file_list,
                         [ ["MINPAR", "TB", 1, "tb"], ["EXTPAR", "MA0", 1, "MA"],
                           ["MASS", "MHH", 1, "MH"], ["MASS", "Mh0", 1, "Mh"],
                           ["MASS", "MNeu(1)", 1, "mneu1"], ["MASS", "MNeu(2)", 1, "mneu2"],
                           ["MASS", "MNeu(3)", 1, "mneu3"], ["MASS", "MNeu(4)", 1, "mneu4"],
                           ["MASS", "MCha(1)", 1, "mcha1"], ["MASS", "MCha(2)", 1, "mcha2"],
                           ["EXTPAR", "M1", 1, "M1"], ["EXTPAR", "M2", 1, "M2"], ["EXTPAR", "MUE", 1, "mu"],
                         ]
        )
        with open(ma_tb_data_set_pkl_fn, 'wb') as fh:
            pickle.dump(ma_tb_data_set, fh)

    ma_array = np.array(ma_tb_data_set["MA"])
    tb_array = np.array(ma_tb_data_set["tb"])
    M1_array = np.array(ma_tb_data_set["M1"])
    M2_array = np.array(ma_tb_data_set["M2"])
    mu_array = np.array(ma_tb_data_set["mu"])

    # Load external curves
    higgs_prop = np.loadtxt("higgs_properties.txt")
    higgs_prop_x = [idata[0] for idata in higgs_prop]
    higgs_prop_y = [idata[1] for idata in higgs_prop]
    higgs_searches = np.loadtxt("higgs_searches.txt")
    higgs_searches_x = [idata[0] for idata in higgs_searches]
    higgs_searches_y = [idata[1] for idata in higgs_searches]

    # Scenario selection
    scenario_selection = np.where((M1_array == 160.) & (M2_array == 180.) & (mu_array == 180.))
    print("Found a total of {} points in the given scenario".format(len(ma_array[scenario_selection])))
    ## Analysis selection
    # Scenario properties
    mass_planes_flag = False
    xs_planes_flag = False
    # CheckMATE global analyses
    global_proc_flag = False
    global_analysis_flag = False
    global_SR_flag = False
    r_study_flag = True
    # Per analysis plots
    analysis_flag = True
    dominant_process_analysis_flag = True
    # Per SR plots
    sr_analysis_flag = True
    dominant_process_analysis_sr_flag = False
    ## End analysis selection

    # Colormap management
    cmap = plt.get_cmap("coolwarm")
    max_r = 2
    shiftedColorMap(cmap, start=0, midpoint=1, stop=max_r, name='coolwarmshifted')
    cmap = plt.get_cmap("coolwarmshifted")
    cmap_viridis = plt.get_cmap("viridis")

    ## Full dataset for the mass planes and the decays

    if (mass_planes_flag == True):
        spectrum_analysis(ma_tb_data_set, scenario_selection)

    if (xs_planes_flag == True):
        ewkino_xs_analysis(ma_tb_data_set, scenario_selection, ewkino_xs)

    global_plot_dir = "global-analysis"
    if not os.path.isdir(global_plot_dir):
        os.makedirs(global_plot_dir)

    # Plot the mA-tanb plane, showing for each point the dominant process
    if (global_proc_flag == True):
        global_proc = get_dominant_proc(data_dir_dict, process_detail_dict)
        global_proc_arr = np.array(global_proc)
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], global_proc_arr[scenario_selection], marker_size = 10)
        set_labels_title(ax, r"Dominant process", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/light-ewkinos-ma-tanb-dominant-proc.pdf".format(global_plot_dir))
        close_fig(fig)

    # Plot the mA-tanb plane, showing for each point the dominant process
    if (global_analysis_flag == True):
        global_proc = get_dominant_analysis(data_dir_dict)
        global_proc_arr = np.array(global_proc)
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], global_proc_arr[scenario_selection], marker_size = 10)
        set_labels_title(ax, r"Best analysis", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/light-ewkinos-ma-tanb-dominant-analysis.pdf".format(global_plot_dir))
        close_fig(fig)

    # Plot the mA-tanb plane, showing for each point the dominant process
    if (global_SR_flag == True):
        global_proc = get_dominant_SR(data_dir_dict)
        global_proc_arr = np.array(global_proc)
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], global_proc_arr[scenario_selection], marker_size = 10)
        set_labels_title(ax, r"Best SR-analysis", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/light-ewkinos-ma-tanb-dominant-analysis-SR.pdf".format(global_plot_dir))
        close_fig(fig)


    ## Plot r vs iteration for a random selection of 100 points
    if (r_study_flag == True):
        fig, ax = create_fig_and_ax()
        max_iter = 0
        for iridx in r_vs_iteration_arr[scenario_selection]:
            x_it = [x+1 for x in range(len(iridx))]
            plot_line(fig, ax, x_it, iridx, alpha = 0.5)
            max_iter = max(len(iridx), max_iter)
        plot_axhline(fig, ax, y = 1, linewidth = 2)
        set_ranges(fig, ax, [1, 30], [0, 2])
        set_labels_title(ax, "r vs iteration (cumulative event number)", "Iteration", "r")
        save_fig(fig, ax, outfilename = "light-ewkinos-r-evolution.pdf")
        close_fig(fig)

        # global best-R in the MA-tanb plane as a function of iteration
        iter_fn_list = []
        for iiter in range(0, max([len(x) for x in r_vs_iteration_arr])):
            ir = []
            for i_iter_point in r_vs_iteration_arr:
                if iiter > len(i_iter_point)-1:
                    ir[len(ir):] = [i_iter_point[-1]]
                else:
                    ir[len(ir):] = [i_iter_point[iiter]]
            ir_array = np.array(ir)
            fig, ax = create_fig_and_ax()
            r_values_array = np.array(r_values)
            plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], ir_array[scenario_selection],
                         marker_size = 10, color_range = [0, max_r], cmap = cmap, cbar_flag = True, cbar_label = "r")
            set_labels_title(ax, r"{{\tt CheckMATE}} analysis, iteration = {}".format(iiter+1), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
            plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
            plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
            set_x_major_ticks(fig, ax, 500)
            set_x_minor_ticks(fig, ax, 100)
            set_y_minor_ticks(fig, ax, 2)
            set_ranges(fig, ax, [90, 2000], [1, 60])
            iter_filename = "{}/light-ewkinos-r-ma-tb-plane-{}.pdf".format(global_plot_dir, iiter)
            iter_fn_list[len(iter_filename):] = [iter_filename]
            save_fig(fig, ax, outfilename = iter_filename)
            close_fig(fig)
        os.system("pdftk {} cat output {}/light-ewkinos-r-ma-tb-plane-iteration-combined.pdf".format(" ".join(iter_fn_list), global_plot_dir))

        # global best-R in the MA-tanb plane
        fig, ax = create_fig_and_ax()
        r_values_array = np.array(r_values)
        plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], r_values_array[scenario_selection],
                     marker_size = 10, color_range = [0, max_r], cmap = cmap, cbar_flag = True, cbar_label = "r")
        set_labels_title(ax, r"{\tt CheckMATE} analysis", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/light-ewkinos-r-ma-tb-plane.pdf".format(global_plot_dir))
        close_fig(fig)

        # Dump the r-values to file
        slha_r = list(zip(slha_file_list, r_values))
        slha_r_arr = np.array(slha_r)
        with open("{}/r-from-results.txt".format(global_plot_dir), "w") as fh:
            fh.write("#SLHA-file\t\tr\n")
            for ipoint in slha_r_arr[scenario_selection]:
                fh.write("{}\t\t{}\n".format(ipoint[0], ipoint[1]))

#### Detailed per-analysis/per-analysis-SR planes


    ## Plot a mA-tanb plane for each analysis (of which the best SR is chosen), showing the r-value with a color-coded scale
    if (analysis_flag == True):
        analysis_plot_dir = "per-analysis"
        if not os.path.isdir(analysis_plot_dir):
            os.makedirs(analysis_plot_dir)

        plot_list = ["totalmcevents", "totalnormevents", "totalsumofweights", "totalsumofweights2", "signalsumofweights", "signalsumofweights2",
                     "signalnormevents", "signal_err_stat", "signal_err_sys", "signal_err_tot", "obs", "bkg", "bkgerr",  "eff", "eff_err_stat",
                     "eff_err_sys", "eff_err_tot", "s95obs", "s95exp", "robs", "robscons", "robsconssysonly", "rexp", "rexpcons", "rexpconssysonly",
                     "SR"]
        for iplot in plot_list:
            for ianalysis in sorted(total_analysis_dict.keys()):
                print("Plotting analysis {} for {}".format(iplot, ianalysis))
                plot_dir = "{}/{}".format(analysis_plot_dir, ianalysis)
                if not os.path.isdir(plot_dir):
                    os.makedirs(plot_dir)
                if (iplot == "SR"):
                    r_values_analysis = np.array(get_max_SR_analysis(ianalysis, total_results_dict, ikey = iplot))
                    r_values_analysis_array = np.array(r_values_analysis)
                else:
                    r_values_analysis = np.array(get_max_r_analysis(ianalysis, total_results_dict, ikey = iplot))
                    r_values_analysis_array = np.array(r_values_analysis)
                fig, ax = create_fig_and_ax()
                set_labels_title(ax, r"{}".format(ianalysis.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
    #            print(r_values_analysis_array[scenario_selection])
                if iplot in ["robs", "robscons", "robsconssysonly", "rexp", "rexpcons", "rexpconssysonly"]:
                    plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], r_values_analysis_array[scenario_selection],
                                 marker_size = 10, color_range = [0,max_r], cmap = cmap, cbar_flag = True, cbar_label = "{}".format(iplot.replace("_","\_")))
                elif (iplot == "SR"):
                    plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], r_values_analysis_array[scenario_selection],
                                 marker_size = 10)
                else:
                    plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], r_values_analysis_array[scenario_selection],
                             marker_size = 10, cmap = cmap_viridis, cbar_flag = True, cbar_label = "{}".format(iplot.replace("_","\_")))
                plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
                plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
                set_x_major_ticks(fig, ax, 500)
                set_x_minor_ticks(fig, ax, 100)
                set_y_minor_ticks(fig, ax, 2)
                set_ranges(fig, ax, [90, 2000], [1, 60])
                save_fig(fig, ax, outfilename = "{}/{}/light-ewkinos-{}-ma-tb-plane-{}.pdf".format(analysis_plot_dir, ianalysis, iplot, ianalysis))
                close_fig(fig)

                # Dump the r-values to file
                slha_r = list(zip(slha_file_list, r_values_analysis))
                slha_r_arr = np.array(slha_r)
                with open("{}/{}/{}-from-{}.txt".format(analysis_plot_dir, ianalysis, iplot, ianalysis), "w") as fh:
                    fh.write("#SLHA-file\t\tr\n")
                    for ipoint in slha_r_arr[scenario_selection]:
                        fh.write("{}\t\t{}\n".format(ipoint[0], ipoint[1]))

    ## Plot a mA-tanb plane for each analysis and SR, showing the r-value with a color-coded scale
    if (sr_analysis_flag == True):
        for ianalysis in sorted(total_analysis_dict.keys()):
            for isr in sorted(total_analysis_dict[ianalysis]):
                r_values_an_sr = np.array(get_max_r_analysis_sr(ianalysis, isr, total_results_dict))
                r_values_an_sr_array = np.array(r_values_an_sr)
                fig, ax = create_fig_and_ax()

                set_labels_title(ax, r"{}~~{}".format(ianalysis.replace("_","\_"), isr.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], r_values_an_sr_array[scenario_selection],
                             marker_size = 10, color_range = [0, max_r], cmap = cmap, cbar_flag = True, cbar_label = "r")
                plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
                plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
                set_x_major_ticks(fig, ax, 500)
                set_x_minor_ticks(fig, ax, 100)
                set_y_minor_ticks(fig, ax, 2)
                set_ranges(fig, ax, [90, 2000], [1, 60])

                sr_analysis_plot_dir = "per-analysis-sr/{}".format(ianalysis)
                if not os.path.isdir(sr_analysis_plot_dir):
                    os.makedirs(sr_analysis_plot_dir)
                save_fig(fig, ax, outfilename = "{}/light-ewkinos-r-ma-tb-plane-{}-sr-{}.pdf".format(sr_analysis_plot_dir, ianalysis, isr))
                close_fig(fig)

                # Dump the r-values to file
                slha_r = list(zip(slha_file_list, r_values_an_sr))
                slha_r_arr = np.array(slha_r)
                with open("{}/r-from-{}-{}.txt".format(sr_analysis_plot_dir, ianalysis, isr), "w") as fh:
                    fh.write("#SLHA-file\t\tr\n")
                    for ipoint in slha_r_arr[scenario_selection]:
                        fh.write("{}\t\t{}\n".format(ipoint[0], ipoint[1]))

    ## Plot a mA-tanb plane for each analysis and SR, showing the dominant process for each point using a color coding
    if (dominant_process_analysis_sr_flag == True):
        for ianalysis in process_analysis_dict:
            for isr in process_analysis_dict[ianalysis]:
                dominant_process_list, numeric_dominance = get_dominant_process_sr_analysis(ianalysis, isr, process_detail_dict)
                dominant_process_arr = np.array(dominant_process_list)
                numeric_dominance_arr = np.array(numeric_dominance)

                fig, ax = create_fig_and_ax()

                set_labels_title(ax, r"dominant proc.~{}~~{}".format(ianalysis.replace("_","\_"), isr.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                print(dominant_process_arr[scenario_selection])
                plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], dominant_process_arr[scenario_selection], marker_size = 10)
                plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
                plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
                set_x_major_ticks(fig, ax, 500)
                set_x_minor_ticks(fig, ax, 100)
                set_y_minor_ticks(fig, ax, 2)
                set_ranges(fig, ax, [90, 2000], [1, 60])

                sr_analysis_plot_dir = "dominant-process-per-analysis-sr/{}".format(ianalysis)
                if not os.path.isdir(sr_analysis_plot_dir):
                    os.makedirs(sr_analysis_plot_dir)
                save_fig(fig, ax, outfilename = "{}/light-ewkinos-dominant-process-ma-tb-plane-{}-sr-{}.pdf".format(sr_analysis_plot_dir, ianalysis, isr))
                close_fig(fig)

                fig, ax = create_fig_and_ax()

                set_labels_title(ax, r"dominant proc.~{}~~{}".format(ianalysis.replace("_","\_"), isr.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                print(dominant_process_arr[scenario_selection])
                plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], numeric_dominance_arr[scenario_selection], marker_size = 10,
                             cbar_flag = True, cbar_label = "Relative importance")
                plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
                plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
                set_x_major_ticks(fig, ax, 500)
                set_x_minor_ticks(fig, ax, 100)
                set_y_minor_ticks(fig, ax, 2)
                set_ranges(fig, ax, [90, 2000], [1, 60])

                save_fig(fig, ax, outfilename = "{}/light-ewkinos-dominant-process-ma-tb-plane-{}-sr-{}-relative-importance.pdf".format(sr_analysis_plot_dir, ianalysis, isr))
                close_fig(fig)



    ## Plot a mA-tanb plane for each analysis (chosing the best SR on a point-by-point basis), showing the dominant process for each point using a color coding
    if (dominant_process_analysis_flag == True):
        for ianalysis in process_analysis_dict:
            print("Plotting dominant process for analysis {}".format(ianalysis))
            # atlas_1708_07875
            #if (ianalysis != "cms_sus_16_039"):
            #    continue
            dominant_process_list, numeric_dominance = get_dominant_process_analysis(ianalysis, process_detail_dict, best_signal_regions_dict)

            dominant_process_arr = np.array(dominant_process_list)
            numeric_dominance_arr = np.array(numeric_dominance)

            fig, ax = create_fig_and_ax()
            set_labels_title(ax, r"dominant proc.~{}".format(ianalysis.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")

            plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], dominant_process_arr[scenario_selection], marker_size = 10)

            plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
            plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
            set_x_major_ticks(fig, ax, 500)
            set_x_minor_ticks(fig, ax, 100)
            set_y_minor_ticks(fig, ax, 2)
            set_ranges(fig, ax, [90, 2000], [1, 60])

            proc_analysis_plot_dir = "dominant-process-per-analysis/{}".format(ianalysis)
            if not os.path.isdir(proc_analysis_plot_dir):
                os.makedirs(proc_analysis_plot_dir)
            save_fig(fig, ax, outfilename = "{}/light-ewkinos-dominant-process-ma-tb-plane-{}.pdf".format(proc_analysis_plot_dir, ianalysis))
            close_fig(fig)

            fig, ax = create_fig_and_ax()
            set_labels_title(ax, r"dominant proc.~{}".format(ianalysis.replace("_","\_")), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")

            plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection], numeric_dominance_arr[scenario_selection], marker_size = 10, cbar_flag = True, cbar_label = "Relative importance")

            plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
            plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
            set_x_major_ticks(fig, ax, 500)
            set_x_minor_ticks(fig, ax, 100)
            set_y_minor_ticks(fig, ax, 2)
            set_ranges(fig, ax, [90, 2000], [1, 60])

            save_fig(fig, ax, outfilename = "{}/light-ewkinos-dominant-process-ma-tb-plane-{}-relative_importance.pdf".format(proc_analysis_plot_dir, ianalysis))
            close_fig(fig)
