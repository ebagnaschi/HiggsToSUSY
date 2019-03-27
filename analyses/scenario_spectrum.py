from plotting import *
from slha_utils import *

def spectrum_analysis():
    slha_file_list_full = find_slha_directories(["SLHA"])
    #slha_file_list_rand = [ slha_file_list_full[irand] for irand in set(np.random.randint(0, len(slha_file_list_full), 3000)) ]
    #print("Randomly selected {} files".format(len(slha_file_list_rand)))

    # Load external curves
    higgs_prop = np.loadtxt("higgs_properties.txt")
    higgs_prop_x = [idata[0] for idata in higgs_prop]
    higgs_prop_y = [idata[1] for idata in higgs_prop]
    higgs_searches = np.loadtxt("higgs_searches.txt")
    higgs_searches_x = [idata[0] for idata in higgs_prop]
    higgs_searches_y = [idata[1] for idata in higgs_prop]

    decay_id_to_str = { 1000022: "neu1", 1000023: "neu2", 1000025: "neu3", 1000035: "neu4",
                        1000024: "cha1p", -1000024: "cha1m", 1000037: "cha2p", -1000037: "cha2m" }

    hh_decay_to_ewkinos = [
        [-1000024, 1000024], [-1000037, 1000024], [-1000024, 1000037], [-1000037, 1000037],
        [1000022, 1000022], [1000022, 1000023], [1000022, 1000025], [1000022, 1000035],
        [1000023, 1000023], [1000023, 1000025], [1000023, 1000035], [1000025, 1000025],
        [1000025, 1000035], [1000035, 1000035],
        ]
    H_decay_list = []
    A_decay_list = []

    for idecay in hh_decay_to_ewkinos:
        H_decay_dict_key = "Hto-{}-{}".format(decay_id_to_str[idecay[0]], decay_id_to_str[idecay[1]])
        H_decay_list[len(H_decay_list):] = [ [35, idecay[0], idecay[1], H_decay_dict_key] ]
        A_decay_dict_key = "Ato-{}-{}".format(decay_id_to_str[idecay[0]], decay_id_to_str[idecay[1]])
        A_decay_list[len(A_decay_list):] = [ [36, idecay[0], idecay[1], A_decay_dict_key] ]

    ma_tb_data_set_full = load_data(slha_file_list_full,
                     [ ["MINPAR", "TB", 1, "tb"], ["EXTPAR", "MA0", 1, "MA"],
                       ["MASS", "MHH", 1, "MH"], ["MASS", "Mh0", 1, "Mh"],
                       ["MASS", "MNeu(1)", 1, "mneu1"], ["MASS", "MNeu(2)", 1, "mneu2"],
                       ["MASS", "MNeu(3)", 1, "mneu3"], ["MASS", "MNeu(4)", 1, "mneu4"],
                       ["MASS", "MCha(1)", 1, "mcha1"], ["MASS", "MCha(2)", 1, "mcha2"],
                       ["EXTPAR", "M1", 1, "M1"], ["EXTPAR", "M2", 1, "M2"], ["EXTPAR", "MUE", 1, "mu"],
                     ],
                    H_decay_list + A_decay_list
    )

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

    BR_A_to_ewkinos = [0 for x in range(0, len(ma_tb_data_set_full["MA"])) ]
    A_decay_list_key = [x[3] for x in A_decay_list]

    for idecay in A_decay_list_key:
        for jj, ipoint_br in enumerate(ma_tb_data_set_full[idecay]):
            BR_A_to_ewkinos[jj] = BR_A_to_ewkinos[jj] + ipoint_br

    ma_array = np.array(ma_tb_data_set_full["MA"])
    tb_array = np.array(ma_tb_data_set_full["tb"])

    M1_array = np.array(ma_tb_data_set_full["M1"])
    M2_array = np.array(ma_tb_data_set_full["M2"])
    mu_array = np.array(ma_tb_data_set_full["mu"])

    scenario_selection = np.where( (M1_array == 160.) & (M2_array == 180.) & (mu_array == 180.)  )
    scenario_str = "light-ewkinos-paper"

    BR_H_to_ewkinos_array = np.array(BR_H_to_ewkinos)
    BR_A_to_ewkinos_array = np.array(BR_A_to_ewkinos)

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
    plot_line_datasets(fig, ax, [higgs_prop_x, higgs_prop_x], alpha = 1, linewidth = 1, linestyle = 'dotted')
    plot_line_datasets(fig, ax, [higgs_searches_x, higgs_searches_x], alpha = 1, linewidth = 1, linestyle = 'dashed')
    # plot_contourf(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
    #             BR_H_to_ewkinos_array[scenario_selection],
    #             level_contour = [x for x in np.arange(0,1.1,0.1) ])
    # plot_pcolormesh(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
    #                 BR_H_to_ewkinos_array[scenario_selection], cbar_flag = True, cbar_label = r"$BR(H \to \tilde{\chi} \tilde{\chi})$")
    set_labels_title(ax, r"$BR(H \to \tilde{\chi} \tilde{\chi})$",
                     r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
    set_x_major_ticks(fig, ax, 500)
    set_x_minor_ticks(fig, ax, 100)
    set_y_minor_ticks(fig, ax, 2)
    set_ranges(fig, ax, [90, 2000], [1, 60])
    save_fig(fig, ax, outfilename = "{}/{}-H-BR-to-ekwinos-ma-tb-contour.pdf".format(H_decay_dir, scenario_str))
    close_fig(fig)

    # Single-decay channel plot
    for idecay in H_decay_list_key:
        fig, ax = create_fig_and_ax()
        BR_array = np.array(ma_tb_data_set_full[idecay])
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     BR_array[scenario_selection],
                     level_contour = [x for x in np.arange(0,1.1,0.1)])
        plot_line_datasets(fig, ax, [higgs_prop_x, higgs_prop_x], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [higgs_searches_x, higgs_searches_x], alpha = 1, linewidth = 1, linestyle = 'dashed')
        set_labels_title(ax, r"{}".format(idecay),
                         r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/{}-H-{}-ma-tb-contour.pdf".format(H_decay_dir, scenario_str, idecay))
        close_fig(fig)


    A_decay_dir = "A-decays"
    if not os.path.isdir(A_decay_dir):
        os.makedirs(A_decay_dir)

    # A to chi chi -- scatterplot
    fig, ax = create_fig_and_ax()
    plot_scatter(fig, ax, ma_array, tb_array, BR_A_to_ewkinos_array, color_range = [0,1],
                 cbar_flag = True, rasterized = True, cbar_label = r"$BR(H \to \tilde{\chi} \tilde{\chi})$")
    set_labels_title(ax, "", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
    save_fig(fig, ax, outfilename = "{}/{}-A-BR-to-ekwinos-ma-tb-scatter.pdf".format(A_decay_dir, scenario_str))
    close_fig(fig)

    # A to chi chi -- contourplot
    fig, ax = create_fig_and_ax()
    plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                 BR_A_to_ewkinos_array[scenario_selection],
                 level_contour = [x for x in np.arange(0,1.1,0.1)])
    set_labels_title(ax, r"$BR(A \to \tilde{\chi} \tilde{\chi})$",
                     r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
    set_x_major_ticks(fig, ax, 500)
    set_x_minor_ticks(fig, ax, 100)
    set_y_minor_ticks(fig, ax, 2)
    set_ranges(fig, ax, [90, 2000], [1, 60])
    save_fig(fig, ax, outfilename = "{}/{}-A-BR-to-ekwinos-ma-tb-contour.pdf".format(A_decay_dir, scenario_str))
    close_fig(fig)

    # Single-decay channel plot
    for idecay in A_decay_list_key:
        fig, ax = create_fig_and_ax()
        BR_array = np.array(ma_tb_data_set_full[idecay])
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     BR_array[scenario_selection],
                     level_contour = [x for x in np.arange(0,1.1,0.1)])
        plot_line_datasets(fig, ax, [higgs_prop_x, higgs_prop_x], alpha = 1, linewidth = 1, linestyle = 'dotted')
        plot_line_datasets(fig, ax, [higgs_searches_x, higgs_searches_x], alpha = 1, linewidth = 1, linestyle = 'dashed')
        set_labels_title(ax, r"{}".format(idecay),
                         r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}/{}-A-{}-ma-tb-contour.pdf".format(A_decay_dir, scenario_str, idecay))
        close_fig(fig)

    # M1 plane
    fig, ax = create_fig_and_ax()
    plot_scatter(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                 M1_array[scenario_selection],  rasterized = True,
                 marker_size = 1,  alpha = 0.7, cbar_flag = True)
    set_labels_title(ax, r"$M_1~\mathrm{[GeV]}$", r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
    save_fig(fig, ax, outfilename = "{}-M1-plane.pdf".format(scenario_str))
    close_fig(fig)

    # Mass planes
    tex_particle= {"neu1": r"\tilde{\chi}^0_1$",  "neu2": r"\tilde{\chi}^0_2$",  "neu3": r"\tilde{\chi}^0_3$",  "neu4": r"\tilde{\chi}^0_4$",
                   "cha1": r"\tilde{\chi}^{\pm}_1$", "cha2":r"\tilde{\chi}^{\pm}_2$"}

    for imass in ["mneu1", "mneu2", "mneu3", "mneu4", "mcha1", "mcha2"]:
        fig, ax = create_fig_and_ax()
        mass_arr = np.array(ma_tb_data_set_full[imass])
        plot_contour(fig, ax, ma_array[scenario_selection], tb_array[scenario_selection],
                     mass_arr[scenario_selection], level_contour = 10)
        print(imass[1:])
        set_labels_title(ax, r"$M_{{}}~\mathrm{[GeV]}$".format(tex_particle[imass[1:]]),r"$M_A~\mathrm{[GeV]}$", r"$\tan\beta$")
        set_x_major_ticks(fig, ax, 500)
        set_x_minor_ticks(fig, ax, 100)
        set_y_minor_ticks(fig, ax, 2)
        set_ranges(fig, ax, [90, 2000], [1, 60])
        save_fig(fig, ax, outfilename = "{}-{}-ma-tb-contour.pdf".format(scenario_str, imass))
        close_fig(fig)
