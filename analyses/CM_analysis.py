import os
import sys
import matplotlib as mpl
# used for the legend
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

def get_r_vs_iteration_data(data_dir_dict):
    """Return an array of list, where each list is r computed by CM as a function of the iteration"""
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
    key_list = []
    for idir in data_dir_dict:
        best_signal_dict[idir] = {}
        best_signal_fn = "{}/evaluation/{}".format(idir, fn)
        with open(best_signal_fn, 'r') as fh:
            for iline in fh:
                isplit = iline.split()
                ## This is always the first line
                if (isplit[0] == "analysis"):
#                    indx_r = isplit.index("robscons")
                    key_list = isplit
                else:
                    analysis = isplit[0]
                    sr = isplit[1]
                    if not (analysis in dict_analysis):
                        dict_analysis[analysis] = set(["{}".format(sr)])
                    else:
                        dict_analysis[analysis].add("{}".format(sr))

                    if not (analysis in best_signal_dict[idir]):
                        best_signal_dict[idir][analysis] = {}
                    best_signal_dict[idir][analysis][sr] = {}
                    for ipos, ikey in enumerate(key_list[2:]):
                        best_signal_dict[idir][analysis][sr][ikey] = isplit[ipos+2]
    return best_signal_dict, dict_analysis


def get_max_r_analysis(analysis, total_results_dict, ikey = "robscons"):
    r_key = "rexpcons"
    ivalues = []
    for idir in sorted(total_results_dict.keys()):
        if not (analysis in total_results_dict[idir]):
            ivalues[len(ivalues):] = [0]
        else:
            all_r = [float(total_results_dict[idir][analysis][sr][r_key]) for sr in total_results_dict[idir][analysis]]
            max_r = max(all_r)
            all_key = [float(total_results_dict[idir][analysis][sr][ikey]) for sr in total_results_dict[idir][analysis]]
            ivalues[len(ivalues):] = [all_key[all_r.index(max_r)]]
    return ivalues

def get_max_SR_analysis(analysis, total_results_dict, ikey = "robscons"):
    r_key = "rexpcons"
    ivalues = []
    for idir in sorted(total_results_dict.keys()):
        if not (analysis in total_results_dict[idir]):
            ivalues[len(ivalues):] = [0]
        else:
            sr_list = sorted(total_results_dict[idir][analysis].keys())
            all_r = [(sr, float(total_results_dict[idir][analysis][sr][r_key])) for sr in sr_list]
            max_r = max(all_r, key=lambda x: x[1])
            ivalues[len(ivalues):] = [sr_list[all_r.index(max_r)]]

    single_SRs = set(ivalues)
    colors = [np.random.rand(3,) for i in range(len(single_SRs))]
    color_dict = {}
    for ii, isr in enumerate(single_SRs):
        color_dict[isr] = colors[ii]
    ivalues_colors = [color_dict[isr] for isr in ivalues]

    legend_files = "per-analysis/{}/analysis_SR_legend.pdf".format(analysis)
    print("Saving {} SR legend to {}".format(analysis, legend_files))
    patch_list = []
    for isr in single_SRs:
        ipatch = mpatches.Patch(color=color_dict[isr], label=isr.replace("_","\_"))
        patch_list.append(ipatch)
    plt.axis('off')
    plt.legend(handles=patch_list, ncol=2)
    plt.axis('off')
    plt.savefig(legend_files)
    os.system("pdfcrop {} {}".format(legend_files, legend_files))

    return ivalues_colors

def get_max_r_analysis_sr(analysis, sr, total_results_dict, r_key = "robscons"):
    # NOTE: SORTED for consistency
    r_values = []
    for idir in sorted(total_results_dict.keys()):
        if not (analysis in total_results_dict[idir]):
            r_values[len(r_values):] = [0]
        else:
            if not (sr in total_results_dict[idir][analysis]):
                r_values[len(r_values):] = [0]
            else:
                r_values[len(r_values):] = [float(total_results_dict[idir][analysis][sr][r_key])]
    return r_values


def load_process_detail_analysis_sr(data_dir_dict):
#    process_list_fn = ["GG_H_MSSM_processResults.txt", "GG_A_MSSM_processResults.txt", "BB_H_MSSM_processResults.txt", "BB_A_MSSM_processResults.txt",                    "other_proc_H_A_MSSM_processResults.txt"#,                    "directewkinos_processResults.txt"]

#    process_name = ["ggH", "ggA", "bbH", "bbA", "otherHA", "directEWKinos"]

    process_list_fn = ["directewkinos_processResults.txt"]

    process_name = ["directEWKinos"]

    dict_analysis = {}

    dict_process_detail = {}

    key_list = []
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
#                        indx_norm_events = isplit.index("signalnormevents")
                        key_list = isplit
                    else:
                        ianalysis = isplit[0]
                        isr = isplit[1]
                        if not (ianalysis in dict_process_detail[idir][process_key]):
                            dict_process_detail[idir][process_key][ianalysis] = {}
                            dict_analysis[ianalysis] = set()
                        dict_analysis[ianalysis].add(isr)
                        dict_process_detail[idir][process_key][ianalysis][isr] = {}
                        for ipos, ikey in enumerate(key_list[2:]):
                            dict_process_detail[idir][process_key][ianalysis][isr][ikey]= isplit[ipos+2]

    return dict_process_detail, dict_analysis

def get_dominant_process_sr_analysis(analysis, sr, process_detail_dict, sig_key = "signalnormevents"):
    process_name = ["ggH", "ggA", "bbH", "bbA", "otherHA", "directEWKinos"]
    color_process = ["red", "blue", "green", "magenta", "pink", "black"]
    # NOTE: SORTED for consistency
    dominant_proc = []
    numeric_dominance = []
    for idir in sorted(process_detail_dict.keys()):
        iproc_list = []
        for iproc in process_name:
            if not (analysis in process_detail_dict[idir][iproc]):
                iproc_list[len(iproc_list):] = [0]
            else:
                if not (sr in process_detail_dict[idir][iproc][analysis]):
                    iproc_list[len(iproc_list):] = [0]
                else:
                    iproc_list[len(iproc_list):] = [float(process_detail_dict[idir][iproc][analysis][sr][sig_key])]
        max_events = max(iproc_list)
        if (max_events > 0):
            idx_max = iproc_list.index(max_events)
            max_proc_key = process_name[idx_max]
            icolor = mpl.colors.to_rgba(color_process[idx_max])
            # change alpha
            relative_importance = max_events/sum(iproc_list)
            icolor = (icolor[0], icolor[1], icolor[2], relative_importance)
        else:
            icolor = (1,1,1,1)
            relative_importance = 0
        dominant_proc[len(dominant_proc):] = [icolor]
        numeric_dominance[len(numeric_dominance):] = [relative_importance]
    return dominant_proc, numeric_dominance

def get_dominant_process_analysis(analysis, process_detail_dict, best_analysis_dict, sig_key = "signalnormevents"):
    process_name = ["ggH", "ggA", "bbH", "bbA", "otherHA", "directEWKinos"]
    color_process = ["red", "blue", "green", "magenta", "pink", "black"]
    # NOTE: SORTED for consistency
    dominant_proc = []
    numeric_dominance = []
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
                        iproc_list[len(iproc_list):] = [float(process_detail_dict[idir][iproc][analysis][best_sr_analysis][sig_key])]
            max_events = max(iproc_list)
        if (max_events > 0):
            idx_max = iproc_list.index(max_events)
            max_proc_key = process_name[idx_max]
            icolor = mpl.colors.to_rgba(color_process[idx_max])
            # change alpha
            relative_importance = max_events/sum(iproc_list)
            icolor = (icolor[0], icolor[1], icolor[2], relative_importance)
        else:
            icolor = (1,1,1,1)
            relative_importance = 0.
        dominant_proc[len(dominant_proc):] = [icolor]
        numeric_dominance[len(numeric_dominance):] = [relative_importance]
    return dominant_proc, numeric_dominance

def get_dominant_proc(data_dir_dict, process_detail_dict, sig_key = "signalnormevents"):
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
            iproc_list[len(iproc_list):] = [float(process_detail_dict[idir][iproc][analysis][best_sr_analysis][sig_key])]
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

def get_dominant_analysis(data_dir_dict):
    color_analysis = {"cms_sus_16_039": "darkred", "cms_sus_16_025": "pink",  "cms_sus_16_025":"lightgreen", "cms_sus_16_048":"blue", "atlas_1712_08119":"black",
                      "atlas_conf_2016_096":"orange"}

    dominant_analysis_list = []
    analysis_set = set()
    # NOTE: SORTED for consistency
    for idir in sorted(data_dir_dict):
        analysis = data_dir_dict[idir]["analysis"]
        analysis_set.add(analysis)
        icolor = mpl.colors.to_rgba(color_analysis[analysis])
        dominant_analysis_list[len(dominant_analysis_list):] = [icolor]

    legend_files = "legends/dominant_analysis_legend.pdf"
    if not os.path.isfile(legend_files):
        print("Saving dominant-analysis legend")
        patch_list = []
        for ianalysis in analysis_set:
            ipatch = mpatches.Patch(color=color_analysis[ianalysis], label=ianalysis.replace("_","\_"))
            patch_list.append(ipatch)
        plt.axis('off')
        plt.legend(handles=patch_list, ncol=2)
        plt.axis('off')
        plt.savefig(legend_files)
        os.system("pdfcrop {} {}".format(legend_files, legend_files))

    return dominant_analysis_list

def get_dominant_SR(data_dir_dict):
    color_analysis = {"SR_G05": "darkred", "SS14": "blue", "SS15": "orange", "SS12": "red", "SR_A14":"pink", "SR_A30":"lightgreen", "SR_A12":"magenta",
                      "SS30": "lightblue", "SR_B02": "orange",  "SR2_stop_1low_pt_1": "cyan", "SR_G03": "darkgreen", "3LI": "darkblue", "SR_G04":"0.75",
                      "SR2_stop_1low_pt_2":"navy", "SR2_stop_3high_pt_1":"darksalmon"}

    dominant_analysis_list = []
    analysis_SR_set = set()
    # NOTE: SORTED for consistency
    for idir in sorted(data_dir_dict):
        analysis = data_dir_dict[idir]["SR"]
        analysis_SR_set.add(analysis)
        icolor = mpl.colors.to_rgba(color_analysis[analysis])
        dominant_analysis_list[len(dominant_analysis_list):] = [icolor]

    legend_files = "legends/dominant_analysis_SR_legend.pdf"
    if not os.path.isfile(legend_files):
        print("Saving dominant-analysis-SR legend")
        patch_list = []
        for ianalysis in analysis_SR_set:
            ipatch = mpatches.Patch(color=color_analysis[ianalysis], label=ianalysis.replace("_","\_"))
            patch_list.append(ipatch)
        plt.axis('off')
        plt.legend(handles=patch_list, ncol=2)
        plt.axis('off')
        plt.savefig(legend_files)
        os.system("pdfcrop {} {}".format(legend_files, legend_files))


    return dominant_analysis_list
