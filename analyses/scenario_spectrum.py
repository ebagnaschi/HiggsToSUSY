import itertools
import os
from plotting import *
from slha_utils import *
from ewkinos_decays import *
import pickle

decay_id_to_str = { 1000022: "neu1", 1000023: "neu2", 1000025: "neu3", 1000035: "neu4",
                    1000024: "cha1p", -1000024: "cha1m", 1000037: "cha2p", -1000037: "cha2m",
                    22: "gamma", 25: "h", 35: "H", 36:"A", 37:"H+", -37:"H-",
                    1000021: "gluino",
                    -1:"dbar", 1:"d", -2:"ubar", 2:"u",
                    -3:"sbar", 3:"s", -4:"cbar", 4:"c",
                    -5:"bbar", 5:"b", -6:"tbar", 6:"t",
                    -11:"e+", 11:"e-", -13:"mu+", 13:"mu", -15:"tau+", 15:"tau-",
                    -12:"nue", 12:"nue", -14:"num", 14:"mum", -16:"ntau", 16:"ntau",
}

def calculate_components_neutralino(NN_row):
    bino_component = NN_row[0]**2
    wino_component = NN_row[1]**2
    higgsino_component = NN_row[2]**2+NN_row[3]**2
    return bino_component, wino_component, higgsino_component

def calculate_components_chargino(CC_row):
    higgsino_component = CC_row[0]**2
    wino_component = CC_row[1]**2
    return wino_component, higgsino_component

def spectrum_analysis(ma_tb_data_set, scenario_selection_CM):

    slha_file_list_full = find_slha_directories(["SLHA"])
    #slha_file_list_rand = [ slha_file_list_full[irand] for irand in set(np.random.randint(0, len(slha_file_list_full), 3000)) ]
    #print("Randomly selected {} files".format(len(slha_file_list_rand)))

    # Load external curves
    higgs_prop = np.loadtxt("higgs_properties.txt")
    higgs_prop_x = [idata[0] for idata in higgs_prop]
    higgs_prop_y = [idata[1] for idata in higgs_prop]
    higgs_searches = np.loadtxt("higgs_searches.txt")
    higgs_searches_x = [idata[0] for idata in higgs_searches]
    higgs_searches_y = [idata[1] for idata in higgs_searches]

    tex_particle= {"neu1": r"\tilde{\chi}^0_1",  "neu2": r"\tilde{\chi}^0_2",  "neu3": r"\tilde{\chi}^0_3",  "neu4": r"\tilde{\chi}^0_4",
                   "cha1": r"\tilde{\chi}^{\pm}_1", "cha2":r"\tilde{\chi}^{\pm}_2",
                   "cha1p": r"\tilde{\chi}^{+}_1", "cha1m":r"\tilde{\chi}^{-}_2",
                   "cha2p": r"\tilde{\chi}^{+}_2", "cha2m":r"\tilde{\chi}^{-}_2",
    }


    hh_decay_to_charginos = [[1000024, -1000024], [1000024, -1000037], [1000037, -1000024], [1000037, -1000037]]
    hh_decay_to_neutralinos = [[1000022, 1000022], [1000022, 1000023], [1000022, 1000025], [1000022, 1000035],
        [1000023, 1000023], [1000023, 1000025], [1000023, 1000035], [1000025, 1000025],
        [1000025, 1000035], [1000035, 1000035]]

    hh_decay_to_ewkinos = [
        [1000024, -1000024], [1000024, -1000037], [1000037, -1000024], [1000037, -1000037],
        [1000022, 1000022], [1000022, 1000023], [1000022, 1000025], [1000022, 1000035],
        [1000023, 1000023], [1000023, 1000025], [1000023, 1000035], [1000025, 1000025],
        [1000025, 1000035], [1000035, 1000035],
        ]
    H_decay_list = []
    A_decay_list = []

    for idecay in hh_decay_to_ewkinos:
        H_decay_dict_key = "H_to_{}_{}".format(decay_id_to_str[idecay[0]], decay_id_to_str[idecay[1]])
        H_decay_list[len(H_decay_list):] = [ [35, idecay[0], idecay[1], H_decay_dict_key] ]
        A_decay_dict_key = "A_to_{}_{}".format(decay_id_to_str[idecay[0]], decay_id_to_str[idecay[1]])
        A_decay_list[len(A_decay_list):] = [ [36, idecay[0], idecay[1], A_decay_dict_key] ]

    neu2_decay_list, neu2_decay_keys = neu2_decay()
    neu3_decay_list, neu3_decay_keys = neu3_decay()
    neu4_decay_list, neu4_decay_keys = neu4_decay()
    cha1_decay_list, cha1_decay_keys = cha1_decay()
    cha2_decay_list, cha2_decay_keys = cha2_decay()

    neu_mix_matrix = []
    for i in range(1,5):
        for j in range(1,5):
            neu_mix_matrix[len(neu_mix_matrix):] = [ ["NMIX", "ZNeu({},{})".format(i,j), 2, "N{}{}".format(i,j)] ]

    chap_mix_matrix = []
    for i in range(1,3):
        for j in range(1,3):
            chap_mix_matrix[len(chap_mix_matrix):] = [ ["UMIX", "UCha({},{})".format(i,j), 2, "U{}{}".format(i,j)] ]

    cham_mix_matrix = []
    for i in range(1,3):
        for j in range(1,3):
            cham_mix_matrix[len(cham_mix_matrix):] = [ ["VMIX", "VCha({},{})".format(i,j), 2, "V{}{}".format(i,j)] ]

    pkl_fn = "ma_tb_slha_full.pkl"
    if os.path.isfile(pkl_fn):
        with open(pkl_fn, 'rb') as fh:
            ma_tb_data_set_full = pickle.load(fh)
    else:
        ma_tb_data_set_full = load_data(slha_file_list_full,
                     [ ["MINPAR", "TB", 1, "tb"], ["EXTPAR", "MA0", 1, "MA"],
                       ["MASS", "MHH", 1, "MH"], ["MASS", "Mh0", 1, "Mh"],
                       ["MASS", "MNeu(1)", 1, "mneu1"], ["MASS", "MNeu(2)", 1, "mneu2"],
                       ["MASS", "MNeu(3)", 1, "mneu3"], ["MASS", "MNeu(4)", 1, "mneu4"],
                       ["MASS", "MCha(1)", 1, "mcha1"], ["MASS", "MCha(2)", 1, "mcha2"],
                       ["EXTPAR", "M1", 1, "M1"], ["EXTPAR", "M2", 1, "M2"], ["EXTPAR", "MUE", 1, "mu"],
                     ] + neu_mix_matrix + chap_mix_matrix + cham_mix_matrix,
                    H_decay_list + A_decay_list + neu2_decay_list + neu3_decay_list + neu4_decay_list + cha1_decay_list + cha2_decay_list
        )
        with open(pkl_fn, 'wb') as fh:
            pickle.dump(ma_tb_data_set_full, fh)

    print("Scenario summary")
    print("M1: {}".format(set(ma_tb_data_set_full["M1"])))
    print("M2: {}".format(set(ma_tb_data_set_full["M2"])))
    print("mu: {}".format(set(ma_tb_data_set_full["mu"])))

    # sum all over the decay channels
    BR_H_to_ewkinos = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
    H_decay_list_key = [x[3] for x in H_decay_list]

    for idecay in H_decay_list_key:
        for jj, ipoint_br in enumerate(ma_tb_data_set_full[idecay]):
            BR_H_to_ewkinos[jj] = BR_H_to_ewkinos[jj] + ipoint_br

    BR_H_to_neutralinos = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
    # decay to two-neutralinos
    for idecay in hh_decay_to_neutralinos:
        idecay_key = "H_to_{}_{}".format(decay_id_to_str[idecay[0]], decay_id_to_str[idecay[1]])
        for jj, ipoint_br in enumerate(ma_tb_data_set_full[idecay_key]):
            BR_H_to_neutralinos[jj] = BR_H_to_neutralinos[jj] + ipoint_br

    BR_H_to_charginos = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
    # decay to two-charginos
    for idecay in hh_decay_to_charginos:
        idecay_key = "H_to_{}_{}".format(decay_id_to_str[idecay[0]], decay_id_to_str[idecay[1]])
        for jj, ipoint_br in enumerate(ma_tb_data_set_full[idecay_key]):
            BR_H_to_charginos[jj] = BR_H_to_charginos[jj] + ipoint_br

    BR_A_to_ewkinos = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
    A_decay_list_key = [x[3] for x in A_decay_list]

    for idecay in A_decay_list_key:
        for jj, ipoint_br in enumerate(ma_tb_data_set_full[idecay]):
            BR_A_to_ewkinos[jj] = BR_A_to_ewkinos[jj] + ipoint_br

    BR_A_to_neutralinos = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
    # decay to two-neutralinos
    for idecay in hh_decay_to_neutralinos:
        idecay_key = "A_to_{}_{}".format(decay_id_to_str[idecay[0]], decay_id_to_str[idecay[1]])
        for jj, ipoint_br in enumerate(ma_tb_data_set_full[idecay_key]):
            BR_A_to_neutralinos[jj] = BR_A_to_neutralinos[jj] + ipoint_br

    BR_A_to_charginos = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
    # decay to two-charginos
    for idecay in hh_decay_to_charginos:
        idecay_key = "A_to_{}_{}".format(decay_id_to_str[idecay[0]], decay_id_to_str[idecay[1]])
        for jj, ipoint_br in enumerate(ma_tb_data_set_full[idecay_key]):
            BR_A_to_charginos[jj] = BR_A_to_charginos[jj] + ipoint_br

    ma_array_CM_points = np.array(ma_tb_data_set["MA"])
    tb_array_CM_points = np.array(ma_tb_data_set["tb"])
    CM_point_array = np.array(len(ma_array_CM_points)*[1])

    ma_array = np.array(ma_tb_data_set_full["MA"])
    tb_array = np.array(ma_tb_data_set_full["tb"])

    M1_array = np.array(ma_tb_data_set_full["M1"])
    M2_array = np.array(ma_tb_data_set_full["M2"])
    mu_array = np.array(ma_tb_data_set_full["mu"])

    scenario_selection = np.where( (M1_array == 160.) & (M2_array == 180.) & (mu_array == 180.)  )
    #scenario_selection = np.where( (M1_array == 68.) & (M2_array == 180.) & (mu_array == 180.)  )
    #scenario_selection = np.where( (M1_array == 1000.) & (M2_array == 1000.) & (mu_array == 1000.)  )
    print("Found a total of {} points in the given scenario".format(len(ma_array[scenario_selection])))
    scenario_str = "light-ewkinos-paper"

    BR_H_to_ewkinos_array = np.array(BR_H_to_ewkinos)
    BR_H_to_neutralinos_array = np.array(BR_H_to_neutralinos)
    BR_H_to_charginos_array = np.array(BR_H_to_charginos)

    BR_A_to_ewkinos_array = np.array(BR_A_to_ewkinos)
    BR_A_to_neutralinos_array = np.array(BR_A_to_neutralinos)
    BR_A_to_charginos_array = np.array(BR_A_to_charginos)


    higgs_decays_flag = True

    if (higgs_decays_flag == True):

        H_decay_dir = "H-decays"
        if not os.path.isdir(H_decay_dir):
            os.makedirs(H_decay_dir)

        # H to chi chi -- scatterplot
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_array, tb_array, BR_H_to_ewkinos_array, color_range = [0,1],
                     cbar_flag = True, rasterized = True, cbar_label = r"$BR(H \to \tilde{\chi} \tilde{\chi})$")
        set_labels_title(ax, "", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        save_fig(fig, ax, outfilename = "{}/{}-H-BR-to-ekwinos-ma-tb-scatter.pdf".format(H_decay_dir, scenario_str))
        close_fig(fig)

        # H to chi chi -- contourplot
        fig, ax = create_fig_and_ax()
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     BR_H_to_ewkinos_array[scenario_selection],
                     level_contour = [x for x in np.arange(0,1.1,0.1)])
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
        set_labels_title(ax, r"$BR(H \to \tilde{\chi} \tilde{\chi})$",
                         r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/{}-H-BR-to-ekwinos-ma-tb-contour.pdf".format(H_decay_dir, scenario_str))
        close_fig(fig)

        # # H to chi0 chi0 -- contourplot
        fig, ax = create_fig_and_ax()
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     BR_H_to_neutralinos_array[scenario_selection],
                     level_contour = [x for x in np.arange(0,1.1,0.1)])
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
        set_labels_title(ax, r"$BR(H \to \tilde{\chi}^0_i \tilde{\chi}^0_j)$",
                         r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/{}-H-BR-to-neutralinos-ma-tb-contour.pdf".format(H_decay_dir, scenario_str))
        close_fig(fig)

        # H to chi^pm chi^pm -- contourplot
        fig, ax = create_fig_and_ax()
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     BR_H_to_charginos_array[scenario_selection],
                     level_contour = [x for x in np.arange(0,1.1,0.1)])
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
        set_labels_title(ax, r"$BR(H \to \tilde{\chi}^{\pm}_i \tilde{\chi}^{\pm}_j)$",
                         r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/{}-H-BR-to-charginos-ma-tb-contour.pdf".format(H_decay_dir, scenario_str))
        close_fig(fig)

        # Single-decay channel plot
        for idecay in H_decay_list_key:
            fig, ax = create_fig_and_ax()
            BR_array = np.array(ma_tb_data_set_full[idecay])
            plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                         BR_array[scenario_selection],
                         level_contour = [x for x in np.arange(0,1.1,0.1)])
            plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
            plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ] , alpha = 1, linewidth = 1, linestyle = 'dashed')
            plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
            idecay_split = idecay.split('_')
            set_labels_title(ax, r"$\mathrm{{BR}}(H \to {p1} {p2})$".format(p1 = tex_particle[idecay_split[-2]], p2 = tex_particle[idecay_split[-1]]),
                             r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
            set_x_major_ticks(fig, ax, 500)
            set_x_minor_ticks(fig, ax, 100)
            set_y_minor_ticks(fig, ax, 2)
            set_ranges(fig, ax, [90, 2000], [1, 60])
            save_fig(fig, ax, outfilename = "{}/{}-{}-ma-tb-contour.pdf".format(H_decay_dir, scenario_str, idecay))
            close_fig(fig)

        A_decay_dir = "A-decays"
        if not os.path.isdir(A_decay_dir):
            os.makedirs(A_decay_dir)

        # A to chi chi -- scatterplot
        fig, ax = create_fig_and_ax()
        plot_scatter(fig, ax, ma_array, tb_array, BR_A_to_ewkinos_array, color_range = [0,1],
                     cbar_flag = True, rasterized = True, cbar_label = r"$BR(A \to \tilde{\chi} \tilde{\chi})$")
        set_labels_title(ax, "", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        save_fig(fig, ax, outfilename = "{}/{}-A-BR-to-ekwinos-ma-tb-scatter.pdf".format(A_decay_dir, scenario_str))
        close_fig(fig)

        # A to chi chi -- contourplot
        fig, ax = create_fig_and_ax()
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     BR_A_to_ewkinos_array[scenario_selection],
                     level_contour = [x for x in np.arange(0,1.1,0.1)])
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
        set_labels_title(ax, r"$BR(A \to \tilde{\chi} \tilde{\chi})$",
                         r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/{}-A-BR-to-ekwinos-ma-tb-contour.pdf".format(A_decay_dir, scenario_str))
        close_fig(fig)

        # A to chi0 chi0 -- contourplot
        fig, ax = create_fig_and_ax()
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     BR_A_to_neutralinos_array[scenario_selection],
                     level_contour = [x for x in np.arange(0,1.1,0.1)])
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
        set_labels_title(ax, r"$BR(A \to \tilde{\chi}^0_i \tilde{\chi}^0_j)$",
                         r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/{}-A-BR-to-neutralinos-ma-tb-contour.pdf".format(A_decay_dir, scenario_str))
        close_fig(fig)

        # A to chi^pm chi^pm -- contourplot
        fig, ax = create_fig_and_ax()
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     BR_A_to_charginos_array[scenario_selection],
                     level_contour = [x for x in np.arange(0,1.1,0.1)])
        plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
        plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
        set_labels_title(ax, r"$BR(A \to \tilde{\chi}^{\pm}_i \tilde{\chi}^{\pm}_j)$",
                         r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/{}-A-BR-to-charginos-ma-tb-contour.pdf".format(A_decay_dir, scenario_str))
        close_fig(fig)

        # Single-decay channel plot
        for idecay in A_decay_list_key:
            fig, ax = create_fig_and_ax()
            BR_array = np.array(ma_tb_data_set_full[idecay])
            plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                         BR_array[scenario_selection],
                         level_contour = [x for x in np.arange(0,1.1,0.1)])
            plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
            plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
            plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
            idecay_split = idecay.split('_')
            set_labels_title(ax, r"$\mathrm{{BR}}(A \to {p1} {p2})$".format(p1 = tex_particle[idecay_split[-2]], p2 = tex_particle[idecay_split[-1]]),
                             r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
            set_x_major_ticks(fig, ax, 500)
            set_x_minor_ticks(fig, ax, 100)
            set_y_minor_ticks(fig, ax, 2)
            set_ranges(fig, ax, [90, 2000], [1, 60])
            save_fig(fig, ax, outfilename = "{}/{}-{}-ma-tb-contour.pdf".format(A_decay_dir, scenario_str, idecay))
            close_fig(fig)

    # M1 plane
    # fig, ax = create_fig_and_ax()
    # plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
    #              M1_array[scenario_selection],  rasterized = True,
    #              marker_size = 1,  alpha = 0.7, cbar_flag = True)
    # set_labels_title(ax, r"$M_1~\mathrm{[GeV]}$", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
    # save_fig(fig, ax, outfilename = "{}-M1-plane.pdf".format(scenario_str))
    # close_fig(fig)

    ewkino_masses_flag = True

    if (ewkino_masses_flag == True):

        # Mass planes
        mass_dir = "ewkinos-masses"
        if not os.path.isdir(mass_dir):
             os.makedirs(mass_dir)
        ewkino_masses = ["mneu1", "mneu2", "mneu3", "mneu4", "mcha1", "mcha2"]
        for imass in ewkino_masses:
            fig, ax = create_fig_and_ax()
            mass_arr = np.array(ma_tb_data_set_full[imass])
            plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                          mass_arr[scenario_selection], level_contour = 10)
            plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
            plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
            plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
            set_labels_title(ax, r"$M_{{ {} }}~\mathrm{{[GeV]}}$".format(tex_particle[imass[1:]]),r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
            set_x_major_ticks(fig, ax, 500)
            set_x_minor_ticks(fig, ax, 100)
            set_y_minor_ticks(fig, ax, 2)
            set_ranges(fig, ax, [90, 2000], [1, 60])
            save_fig(fig, ax, outfilename = "{}/{}-{}-ma-tb-contour.pdf".format(mass_dir, scenario_str, imass))
            close_fig(fig)

        # Mass difference planes
        for imass1, imass2 in itertools.combinations(ewkino_masses, 2):
            delta_arr = np.array([abs(abs(x[0])-abs(x[1])) for x in zip(ma_tb_data_set_full[imass1], ma_tb_data_set_full[imass2])])
            fig, ax = create_fig_and_ax()
            plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                          delta_arr[scenario_selection], level_contour = 10)
            plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
            plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
            plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
            set_labels_title(ax, r"$|M_{{ {} }}-M_{{ {} }}|~\mathrm{{[GeV]}}$".format(tex_particle[imass1[1:]], tex_particle[imass2[1:]]),r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
            set_x_major_ticks(fig, ax, 500)
            set_x_minor_ticks(fig, ax, 100)
            set_y_minor_ticks(fig, ax, 2)
            set_ranges(fig, ax, [90, 2000], [1, 60])
            save_fig(fig, ax, outfilename = "{}/{}-delta-{}-{}-ma-tb-contour.pdf".format(mass_dir, scenario_str, imass1, imass2))
            close_fig(fig)


    decay_product = { "gluinoqqbar": r"\tilde{g} q \bar{q}",
                      "neu1ll": r"\tilde{\chi}^0_1 l \bar{l}",
                      "neu1qqbar": r"\tilde{\chi}^0_1 q \bar{q}",
                      "char1nul": r"\tilde{\chi}^{\pm}_1 \nu l^{\mp}",
                      "char2nul": r"\tilde{\chi}^{\pm}_2 \nu l^{\mp}",
                      "char1qbarq": r"\tilde{\chi}^{\pm}_1 q \bar{q}",
                      "char2qbarq": r"\tilde{\chi}^{\pm}_1 q \bar{q}",
                      "char1qqp": r"\tilde{\chi}^{\pm}_1 q q'",
                      "char2qqp": r"\tilde{\chi}^{\pm}_1 q q'",
                      "squarkquark": r"\tilde{q} q'",
                      "sleptonneu": r"\tilde{l} \nu",
                      "neu1lnu": r"\tilde{\chi}^0_1 l \nu",
                      "neu1qqp": r"\tilde{\chi}^0_1 q q'",
                      "neu1A": r"\tilde{\chi}^0_1 A",
                      "neu1H": r"\tilde{\chi}^0_1 H",
                      "neu1h": r"\tilde{\chi}^0_1 h",
                      "neu1gamma": r"\tilde{\chi}^0_1 \gamma",
                      "neu1Z": r"\tilde{\chi}^0_1 Z",
                      "neu1W-": r"\tilde{\chi}^0_1 W^{-}",
                      "neu1W+": r"\tilde{\chi}^0_1 W^{+}",
                      "neu1H+": r"\tilde{\chi}^0_1 H^{+}",
                      "neu1H-": r"\tilde{\chi}^0_1 H^{-}",
                      "neu2gamma": r"\tilde{\chi}^0_2 \gamma",
                      "neu1nunu": r"\tilde{\chi}^0_1 \nu \nu",
                      "neu2ll": r"\tilde{\chi}^0_2 l^{\pm} l^{\mp}",
                      "neu2A": r"\tilde{\chi}^0_2 A",
                      "neu2H": r"\tilde{\chi}^0_2 H",
                      "neu2h": r"\tilde{\chi}^0_2 h",
                      "neu2Z": r"\tilde{\chi}^0_2 Z",
                      "neu2H+": r"\tilde{\chi}^0_2 H^{+}",
                      "neu2H-": r"\tilde{\chi}^0_2 H^{-}",
                      "neu2W+": r"\tilde{\chi}^0_2 W^{+}",
                      "neu2W-": r"\tilde{\chi}^0_2 W^{-}",
                      "neu2nunu": r"\tilde{\chi}^0_2 \nu \nu",
                      "neu2qqbar": r"\tilde{\chi}^0_2 q \bar{q}",
                      "neu3A": r"\tilde{\chi}^0_3 A",
                      "neu3H": r"\tilde{\chi}^0_3 H",
                      "neu3h": r"\tilde{\chi}^0_3 h",
                      "neu3Z": r"\tilde{\chi}^0_3 Z",
                      "neu3H+": r"\tilde{\chi}^0_3 H^{+}",
                      "neu3H-": r"\tilde{\chi}^0_3 H^{-}",
                      "neu3W+": r"\tilde{\chi}^0_3 W^{+}",
                      "neu3W-": r"\tilde{\chi}^0_3 W^{-}",
                      "neu4A": r"\tilde{\chi}^0_4 A",
                      "neu4H": r"\tilde{\chi}^0_4 H",
                      "neu4h": r"\tilde{\chi}^0_4 h",
                      "neu4Z": r"\tilde{\chi}^0_4 Z",
                      "neu4H+": r"\tilde{\chi}^0_4 H^{+}",
                      "neu4H-": r"\tilde{\chi}^0_4 H^{-}",
                      "neu4W+": r"\tilde{\chi}^0_4 W^{+}",
                      "neu4W-": r"\tilde{\chi}^0_4 W^{-}",
                      "cha1W": r"\tilde{\chi}^{\pm}_1 W^{\mp}",
                      "cha1Hpm": r"\tilde{\chi}^{\pm}_1 H^{\mp}",
                      "cha1pA": r"\tilde{\chi}^{+}_1 A",
                      "cha1pZ": r"\tilde{\chi}^{+}_1 Z",
                      "cha1ph": r"\tilde{\chi}^{+}_1 h",
                      "cha1pH": r"\tilde{\chi}^{+}_1 H",
    }


    ewkino_decays_flag = False

    if (ewkino_decays_flag == True):

        # EWKino decay planes
        fh_tmp = open("temp.dat", 'w')
        ewkinos_decay_dir = "ewkino-decay-dir"
        if not os.path.isdir(ewkinos_decay_dir):
            os.makedirs(ewkinos_decay_dir)
        #decay_list = [ ("neu2", neu2_decay_keys), ("neu3", neu3_decay_keys), ("neu4", neu4_decay_keys), ("cha1p", cha1_decay_keys), ("cha2p", cha2_decay_keys) ]
        decay_list = [  ("neu3", neu3_decay_keys) ]
        for iewkino, idecay_groups in decay_list:
            tot_particle = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
            for idecay in sorted(idecay_groups):
                fh_tmp.write("{}\t{}\n".format(iewkino, idecay))
                print("Plotting {} to {}".format(iewkino, idecay))
                # Sum up all the relevant decays
                BR_ewkino = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
                for idecay_key in idecay_groups[idecay]:
                    for jj, ipoint_br in enumerate(ma_tb_data_set_full[idecay_key]):
                        BR_ewkino[jj] = BR_ewkino[jj] + ipoint_br
                        tot_particle[jj] = tot_particle[jj] + ipoint_br
                BR_ewkino_array = np.array(BR_ewkino)
                # plot the BR in the mA-tanb plane
                fig, ax = create_fig_and_ax()
    #            plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
    #                         BR_ewkino_array[scenario_selection], level_contour = 20, fmt = '%1.3e')
                mother_particle = idecay.split('_')[0]
                children_particle = "".join(idecay.split('_')[2:])
                tex_name = r"{} \to {}".format(tex_particle[mother_particle], decay_product[children_particle])
                plot_pcolormesh(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                                BR_ewkino_array[scenario_selection], cbar_flag = True, cbar_label = r"$BR({})$".format(tex_name), num = 200, rasterized = True)
                plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
                plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
                #plot_scatter(fig, ax, ma_array_CM_points[scenario_selection_CM], tb_array_CM_points[scenario_selection_CM], CM_point_array[scenario_selection_CM], marker_size = 1)
                set_labels_title(ax, r"$BR({})$".format(tex_name), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                set_x_major_ticks(fig, ax, 500)
                set_x_minor_ticks(fig, ax, 100)
                set_y_minor_ticks(fig, ax, 2)
                set_ranges(fig, ax, [90, 2000], [1, 60])
                save_fig(fig, ax, outfilename = "{}/{}-BR-to-{}-ma-tb-contour.pdf".format(ewkinos_decay_dir, idecay, scenario_str))
                close_fig(fig)
    #        print(tot_particle[0:100])
    #        print(slha_file_list_full[0:100])
            check_array = all([abs(x-1) < 1e-4 for x in tot_particle[0:100]])
    #        check_array = [abs(x-1) < 1e-4 for x in tot_particle[0:100]]
    #        print(slha_file_list_full[check_array.index(False)])
            assert(check_array == True)
        fh_tmp.close()


    ewkino_component_flag = True

    if (ewkino_component_flag == True):

        # cycle over the neturalinos
        ewkinos_component_dir = "ewkino-component-dir"
        if not os.path.isdir(ewkinos_component_dir):
            os.makedirs(ewkinos_component_dir)

        # neutralino components
        for i in range(1,5):
            # retrieve the elements
            NN_row = [ np.array(ma_tb_data_set_full["N{}{}".format(i, j)]) for  j in range(1,5) ]
            bino_component, wino_component, higgsino_component = calculate_components_neutralino(NN_row)
            assert(all(abs(bino_component + wino_component + higgsino_component - 1) < 1e-4))
            neutralino_component_label = ["$N_{{ {} 1 }}^2$".format(i), "$N_{{ {} 2}}^2$".format(i), "$N_{{ {} 3 }}^2 +  N_{{ {} 4 }}^2$".format(i,i)]
            string_component = ["bino", "wino", "higgsino"]
            for ic, icomponent in enumerate([bino_component, wino_component, higgsino_component]):
                fig, ax = create_fig_and_ax()
                plot_pcolormesh(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                                icomponent[scenario_selection], cbar_flag = True, cbar_label = neutralino_component_label[ic], num = 200, rasterized = True)
                plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
                plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
                print(ic)
                set_labels_title(ax, r"$\tilde{{\chi}}^0_{}$ {} component".format(i, string_component[ic]), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                set_x_major_ticks(fig, ax, 500)
                set_x_minor_ticks(fig, ax, 100)
                set_y_minor_ticks(fig, ax, 2)
                set_ranges(fig, ax, [90, 2000], [1, 60])
                save_fig(fig, ax, outfilename = "{}/{}-{}-ma-tb-contour.pdf".format(ewkinos_component_dir, "neu{}".format(i), string_component[ic], scenario_str))
                close_fig(fig)

        # chargino+ components
        for i in range(1,3):
            # retrieve the elements
            U_row = [ np.array(ma_tb_data_set_full["U{}{}".format(i, j)]) for  j in range(1,3) ]
            wino_component, higgsino_component = calculate_components_chargino(U_row)
            assert(all(abs(wino_component + higgsino_component - 1) < 1e-4))
            neutralino_component_label = ["$U_{{ {} 1}}^2$".format(i), "$U_{{ {} 2 }}^2$".format(i)]
            string_component = ["wino", "higgsino"]
            for ic, icomponent in enumerate([wino_component, higgsino_component]):
                fig, ax = create_fig_and_ax()
                plot_pcolormesh(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                                icomponent[scenario_selection], cbar_flag = True, cbar_label = neutralino_component_label[ic], num = 200, rasterized = True)
                plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
                plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
                set_labels_title(ax, r"$\tilde{{\chi}}^+_{}$ {} component".format(i, string_component[ic]), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                set_x_major_ticks(fig, ax, 500)
                set_x_minor_ticks(fig, ax, 100)
                set_y_minor_ticks(fig, ax, 2)
                set_ranges(fig, ax, [90, 2000], [1, 60])
                save_fig(fig, ax, outfilename = "{}/{}-{}-ma-tb-contour.pdf".format(ewkinos_component_dir, "cha{}p".format(i), string_component[ic], scenario_str))
                close_fig(fig)

        # chargino- components
        for i in range(1,3):
            # retrieve the elements
            V_row = [ np.array(ma_tb_data_set_full["V{}{}".format(i, j)]) for  j in range(1,3) ]
            wino_component, higgsino_component = calculate_components_chargino(V_row)
            assert(all(abs(wino_component + higgsino_component - 1) < 1e-4))
            neutralino_component_label = ["$V_{{ {} 1}}^2$".format(i), "$V_{{ {} 2 }}^2$".format(i)]
            string_component = ["wino", "higgsino"]
            for ic, icomponent in enumerate([wino_component, higgsino_component]):
                fig, ax = create_fig_and_ax()
                plot_pcolormesh(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                                icomponent[scenario_selection], cbar_flag = True, cbar_label = neutralino_component_label[ic], num = 200, rasterized = True)
                plot_line_datasets(fig, ax, [ [higgs_prop_x, higgs_prop_y] ], alpha = 1, linewidth = 1, linestyle = 'dotted')
                plot_line_datasets(fig, ax, [ [higgs_searches_x, higgs_searches_y] ], alpha = 1, linewidth = 1, linestyle = 'dashed')
                set_labels_title(ax, r"$\tilde{{\chi}}^-_{}$ {} component".format(i, string_component[ic]), r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
                set_x_major_ticks(fig, ax, 500)
                set_x_minor_ticks(fig, ax, 100)
                set_y_minor_ticks(fig, ax, 2)
                set_ranges(fig, ax, [90, 2000], [1, 60])
                save_fig(fig, ax, outfilename = "{}/{}-{}-ma-tb-contour.pdf".format(ewkinos_component_dir, "cha{}m".format(i), string_component[ic], scenario_str))
                close_fig(fig)
