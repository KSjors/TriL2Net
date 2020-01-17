import logging

logging_level = logging.INFO
# set up logging to file - see previous section for more details
logging.basicConfig(level=logging_level,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='main.log',
                    filemode='w')
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)

fh = logging.FileHandler(filename='file_log.log')
fh.setLevel(logging_level)
# set a format which is simpler for console use
formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
console.setFormatter(formatter)
fh.setFormatter(formatter)
# add the handler to the root logger
logging.getLogger('').addHandler(console)
# logging.getLogger('network').addHandler(fh)

from data.reference_networks import get_standard_binets_trinets, regenerate_standard_binets_trinets, TRINET_LIST
from data.generators import GENERATOR_DICT
from utils.help_functions import *
from tqdm import tqdm
import matplotlib.pyplot as plt
from experiments.experiments import TrinetSetGenerator, SolverGenerator, ExperimentGenerator, Experiment, ExperimentSet, \
    NetworkGeneratorSubnets
from config import settings
from datastructures.rooted_level_k_network import NetworkSet, RootedLevelKNetwork
import pandas as pd
import winsound
from scipy import optimize

if __name__ == '__main__':
    print(TRINET_LIST.volume())
    kk
    # regenerate_standard_binets_trinets()
    # kk

    # Processing power
    max_processes = 16

    # Experiment settings
    experiment_name = 'consistency_scores\\tail_move\\level_2\\FULL'
    repeat = 5

    # Network Generator
    generators = [GENERATOR_DICT['2a'], GENERATOR_DICT['2b'], GENERATOR_DICT['2c'], GENERATOR_DICT['2d']]
    # generators = [GENERATOR_DICT['1']]
    num_leaves = [15, 20, 25]
    subnets = NetworkSet.subnets_from_generators(generators, 2)  # TODO important change for level 1 vs 2
    ng = NetworkGeneratorSubnets(num_leaves, subnets)

    # Trinet Set Generator
    tail_move_prob = [0.01, 0.05, 0.10, 0.20]
    uni_prob = [0]
    del_prob = [0]
    max_replacement_level = 2
    tsg = TrinetSetGenerator(tail_move_prob, uni_prob, del_prob, max_replacement_level)

    # Solver Generator
    parameters = {
        'cut_arc_set_count_method'          : [settings.MAXIMUM_MULTIPLICITY]
        , 'minimal_sink_set_method'         : [settings.EXPAND_FIRST_SCC]
        , 'leaf_locator_method'             : [settings.DEFAULT]
        , 'level_threshold_method'          : [settings.DEFAULT]
        , 'level_count_method'              : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'generator_count_method'          : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'symmetric_sides_set_count_method': [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'leaf_order_count_method'         : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'leaf_order_method'               : [settings.DEFAULT]
        , 'fill_gaps'                       : [0]
    }
    # [0, 0, 0, 0, 1, 1, 1, 1, 0, 0]
    sg = SolverGenerator(**parameters)
    #
    """ ---------------------------------------------------------------------------------------------------
                                            EXPERIMENT GENERATOR
        # ---------------------------------------------------------------------------------------------------"""
    # eg = ExperimentGenerator(ng, tsg, sg, name=experiment_name, max_processes=max_processes, consistency_method=settings.WEIGHTED_AVERAGE)
    # eg.run_times_2(times=repeat, reuse_network=True, reuse_trinet_set=True, compute=[1, 0, 0, 0, 0, 0])

    """ ---------------------------------------------------------------------------------------------------
                                            EXPERIMENT SET
        ---------------------------------------------------------------------------------------------------"""
    es = ExperimentSet.load(f'C:\\Users\\sjors\\PycharmProjects\\TriLoNet-2\\experiments\\{experiment_name}', experiment_name)
    data = es.to_df()
    print(data.keys())
    # x_axis = ['parameters solver level_count_method', 'parameters distortion tail_move']
    # y_axis = ['consistency scores trinet']
    # data_2 = data[y_axis + x_axis]
    # data_2 = data_2.groupby(x_axis).std()
    # print(data_2)
    # data_2.unstack(0).plot.bar(legend=False)
    # plt.ylim([0, 1])
    # plt.xlabel('tail move (%)')
    # plt.ylabel('trinet set consistency')
    # plt.xticks([0, 1, 2, 3], ['1', '5', '10', '20'], rotation='horizontal')
    # plt.show()
    # kk

    # data = data[data['parameters input network n'] != 20]

    #
    " --- Plots -----------"


    def func(n, *args):
        return sum(p * n ** i for i, p in enumerate(args))


    def func_2(n, a1, a4, a5, a6, a7):
        return a1 * n + a4 ** 4 + a5 * n ** 5 + a6 * n ** 6 + a7 * n ** 7


    # x_axis = 'parameters input network n'
    # x_axis = 'parameters solver level_count_method'
    x_axis = 'parameters distortion tail_move'
    # y_axes = ['run times solving']
    y_axes = ['trinet', 'triplet', 'cluster']
    # data_std = data.groupby([x_axis], as_index=False).std()
    data_mean = data.groupby([x_axis], as_index=False).mean()

    # x_ticks = [0.01, 0.02, 0.05, 0.10, 0.15, 0.2]
    x_ticks = [0.01, 0.05, 0.10, 0.2, 0.3, 0.4, 0.5]

    cmap = plt.get_cmap("tab10")
    for i, y_axis in enumerate(y_axes):
        x, y = data[x_axis].values, data_mean['consistency scores ' + y_axis].values
        plt.plot(x_ticks, y, 'o', label=y_axis, color=cmap(i))
        p0 = [1, -1, -1, -1, -1]
        # p0 = [1, 1, 1, 1, 1]
        popt, pcov = optimize.curve_fit(func, xdata=x_ticks, ydata=y, bounds=([0.9999, -np.inf, -np.inf, -np.inf, -np.inf], [1, 0, 0, 0, 0]), p0=p0)
        # popt, pcov = optimize.curve_fit(func, xdata=[0.01, 0.02, 0.05, 0.10, 0.15, 0.2-], ydata=y, bounds=(0, np.inf), p0=p0)

        print()
        print(popt)
        x_data = np.linspace(0, 0.55, 100)
        y_data = [func(xi, *popt) for xi in x_data]
        plt.plot(x_data, y_data, color=cmap(i))
        plt.hlines(np.arange(0.1, 1.1, 0.1), 0, 1, linestyles='dotted', colors='grey')
        # plt.vlines([15, 20, 25, 30, 35, 40], 0, 800, linestyles='dotted', colors='grey')
    xlim = [0, 0.55]
    ylim = [0, 1]
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
               ncol=3, mode="expand", borderaxespad=0.)
    # plt.title()
    plt.ylabel('consitency score')
    plt.xlabel('tail move (%)')
    # plt.xticks(x_ticks, ['1', '2', '5', '10', '15', '20'], rotation='horizontal')
    plt.xticks(x_ticks, ['1', '5', '10', '20', '30', '40', '50'], rotation='horizontal')

    # plt.xticks([15, 20, 25, 30, 35, 40], ['15', '20', '25', '30', '35', '40'], rotation='horizontal')
    plt.show()

    #
    # # popt = es.plot_(x_axis, 'consistency scores distorted trinet set'
    # #                 , fit=fit, fit_order=2, data=data)
    # # popt = es.plot_(x_axis, 'summaries distorted trinet set density-0'
    # #                 , fit=False, fit_order=2, bounds=([0.99, -np.inf], [1, 0]), p0=[1, -1], data=data)
    # # popt = es.plot_(x_axis, 'run times solving'
    # #                 , fit=fit, fit_order=5, bounds=(0, np.inf), data=data)
    # # print(popt)
    # popt = es.plot_(x_axis, 'consistency scores trinet'
    #                 , fit=fit, fit_order=5, data=data, bounds=([0.99, -np.inf, -np.inf, -np.inf, -np.inf], [1, 0, 0, 0, 0]), p0=[1, -1, -1, -1, -1])
    # print(popt)
    # # popt = es.plot_(x_axis, 'consistency scores triplet'
    # #                 , fit=fit, fit_order=3, data=data, bounds=([0.99, -np.inf, -np.inf, -np.inf], [1, 0, 0, 0]), p0=[1, -1, -1, -1])
    # # popt = es.plot_(x_axis, 'consistency scores cluster'
    # #                 , fit=fit, fit_order=3, data=data, bounds=([0.99, -np.inf, -np.inf, -np.inf], [1, 0, 0, 0]), p0=[1, -1, -1, -1])
    # # popt = es.plot_(x_axis, 'consistency scores cut-arc set'
    # #                 , fit=fit, fit_order=2, data=data)
    #
    # # plot ...
    # # popt = es.plot_('parameters input network recombination', 'summaries input network total reticulations', fit=True, fit_order=3)
    # # popt = es.plot_('summaries input network total reticulations', 'run times solving', fit=False, fit_order=3)
    #
    # """ ---------------------------------------------------------------------------------------------------
    #                                         EXPERIMENT
    #     ---------------------------------------------------------------------------------------------------"""
    # e = es[16]
    # e = e.output_network

    # f = RootedLevelKNetwork.from_enewick(ene)
    # f.visualize()

    # """--  Load sets -------------------  """
    # # itns, otns = e.recreate_network_set('input trinet set'), e.recreate_network_set('output trinet set')
    # # itps, otps = e.recreate_network_set('input triplet set'), e.recreate_network_set('output triplet set')
    # # ics, ocs = e.recreate_network_set('input cluster set'), e.recreate_network_set('output cluster set')
    #
    # """--  Recreate sets -------------------  """
    # # input_network, output_network = e.input_network, e.output_network
    # # itns, otns = NetworkSet.induced_strict_network_set(input_network, 3), NetworkSet.induced_strict_network_set(output_network, 3)
    # # itps, otps = NetworkSet.induced_strict_tree_set_of_network_set(itns, 3), NetworkSet.induced_strict_tree_set_of_network_set(otns, 3)
    # # ics, ocs = NetworkSet.induced_cluster_set(input_network), NetworkSet.induced_cluster_set(output_network)
    #
    # """--  Compare sets -------------------  """
    # # tn_cs = itns.consistency_score(otns, settings.WEIGHTED_AVERAGE)
    # # tp_cs = itps.consistency_score(otps, settings.WEIGHTED_AVERAGE)
    # # c_cs = ics.consistency_score(ocs, settings.WEIGHTED_AVERAGE)
    # # print(tn_cs, tp_cs, c_cs)
    #
    # """--  Manual Compare sets -------------------  """
    # # input_network.visualize()
    # # output_network.visualize()
    # # print(itps.volume())
    # # print(itns.volume())
    #
    # """--  Rerun -------------------  """
    # # e.rerun(all_trinet_list, settings.ITERATIVE, 1)
    # # e.plot_consistency_scores()
    # # e.visualize_input_output()
    #
    # """ ---------------------------------------------------------------------------------------------------
    #                                         TEST
    #     ---------------------------------------------------------------------------------------------------"""
    #
    # # frequency = 1000  # Set Frequency To 2500 Hertz
    # # duration = 1000  # Set Duration To 1000 ms == 1 second
    # # for _ in range(10):
    # #     winsound.Beep(frequency, duration)
