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
logging.getLogger('network').addHandler(fh)

from data.all_trinets import get_standard_networks, regenerate_standard_networks
from utils.help_functions import *
from tqdm import tqdm
import matplotlib.pyplot as plt
from experiments.experiments import NetworkGenerator, TrinetSetGenerator, SolverGenerator, ExperimentGenerator, Experiment, ExperimentSet
from config import settings
from datastructures.rooted_level_k_network import NetworkSet, RootedLevelKNetwork

if __name__ == '__main__':
    regenerate_standard_networks()
    _, biconnected_trinet_binet_list, all_trinet_list = get_standard_networks()

    # Network Generator
    num_leaves = [10, 15, 20, 25]
    recomb_prob = [0.05]
    term_prob = [0.05]
    ng = NetworkGenerator(num_leaves, recomb_prob, term_prob)

    # Trinet Set Generator
    tail_move_prob = [0.05]
    uni_prob = [0]  # TODO -> not working atm
    del_prob = [0]
    tsg = TrinetSetGenerator(tail_move_prob, uni_prob, del_prob, all_trinet_list)

    # Solver Generator
    parameters = {
        'cut_arc_set_count_method'          : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'minimal_sink_set_method'         : [settings.DEFAULT_SINK_SET]
        , 'leaf_locator_method'             : [settings.GREEDY]  # , settings.ILP]
        , 'level_threshold_method'          : [settings.DEFAULT_THRESHOLD]
        , 'level_count_method'              : [settings.WEIGHTED_AVERAGE]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'generator_count_method'          : [settings.WEIGHTED_AVERAGE]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'symmetric_sides_set_count_method': [settings.WEIGHTED_AVERAGE]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'leaf_order_count_method'         : [settings.WEIGHTED_AVERAGE]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'leaf_order_method'               : [settings.DEFAULT_ORDER]
    }
    sg = SolverGenerator(all_trinet_list, **parameters)

    """ ---------------------------------------------------------------------------------------------------
                                            EXPERIMENT GENERATOR
        ---------------------------------------------------------------------------------------------------"""
    eg = ExperimentGenerator(ng, tsg, sg, name="4", trinet_set_method=settings.ITERATIVE, max_processes=8)
    es = eg.run_times(times=1)

    """ ---------------------------------------------------------------------------------------------------
                                            EXPERIMENT SET
        ---------------------------------------------------------------------------------------------------"""
    es = ExperimentSet.load(r'C:\Users\sjors\PycharmProjects\TriLoNet-2\experiments\4', 'Rerun')
    es.plot_consistency_scores()
    # df = es.to_df()
    # print(df.keys())
    # df = df.groupby(['number'])

    """ ---------------------------------------------------------------------------------------------------
                                            EXPERIMENT
        ---------------------------------------------------------------------------------------------------"""
    # e = Experiment.load(r'C:\Users\sjors\PycharmProjects\TriLoNet-2\experiments\3\2019-12-20 153130-233317.txt', '1')
    # e.plot_consistency_scores()
    # e.visualize_input_output()

    """--  Load setse -------------------  """
    # itns, otns = e.recreate_network_set('input trinet set'), e.recreate_network_set('output trinet set')
    # itps, otps = e.recreate_network_set('input triplet set'), e.recreate_network_set('output triplet set')
    # ics, ocs = e.recreate_network_set('input cluster set'), e.recreate_network_set('output cluster set')

    """--  Recreate sets -------------------  """
    # input_network, output_network = e.input_network, e.output_network
    # itns, otns = NetworkSet.induced_strict_network_set(input_network, 3), NetworkSet.induced_strict_network_set(output_network, 3)
    # itps, otps = NetworkSet.induced_strict_tree_set_of_network_set(itns, 3), NetworkSet.induced_strict_tree_set_of_network_set(otns, 3)
    # ics, ocs = NetworkSet.induced_cluster_set(input_network), NetworkSet.induced_cluster_set(output_network)

    """--  Compare sets -------------------  """
    # tn_cs = itns.consistency_score(otns)
    # tp_cs = itps.consistency_score(otps)
    # c_cs = ics.consistency_score(ocs)
    # print(tn_cs, tp_cs, c_cs)

    """--  Rerun -------------------  """
    # e.rerun(all_trinet_list, settings.ITERATIVE, 1)
    # e.plot_consistency_scores()
    # e.visualize_input_output()

    # Analyse
    # es = ExperimentSet.load(r'C:\Users\sjors\PycharmProjects\TriLoNet-2\experiments\1', '1-Analysis')
    # es = ExperimentSet.filter(es, 'trinet', 1, False)
    # es.plot_consistency_scores()
    # e = es[0]
    # e.rerun(all_trinet_list, settings.ITERATIVE, 8)
    # e.plot_consistency_scores()
    # ni, no = e.input_network, e.output_network
    # itns ,otns = NetworkSet.induced_strict_network_set(ni, 3), NetworkSet.induced_strict_network_set(no, 3)
    # itps, otps = NetworkSet.induced_strict_tree_set_of_network_set(itns, 3), NetworkSet.induced_strict_tree_set_of_network_set(otns, 3)
    # itps, otps = e.recreate_network_set('input triplet set'), e.recreate_network_set('output triplet set')
    # print(itps.consistency_score(otps, method=settings.WEIGHTED_SUM))
    # e.plot_consistency_scores()
    # # ssi, sso = itps[('a', 'd', 'j')], otps[('a', 'd', 'j')]
    # # print(ssi)
    # # print([network_info.multiplicity for network_info in ssi['network_info_set']])
    # # print(sso)
    # # print([network_info.multiplicity for network_info in sso['network_info_set']])
    # # for network_info in ssi['network_info_set']:
    # #     network_info.network.visualize()
    # # for network_info in sso['network_info_set']:
    # #     network_info.network.visualize()

    # e = Experiment.load(r'C:\Users\sjors\PycharmProjects\TriLoNet-2\experiments\failures\1\2019-12-19 165627-489870.txt', 'Rerun-2')
    # e.rerun(all_trinet_list, settings.ITERATIVE, 1)
    # itps, otps = e.recreate_network_set('input triplet set'), e.recreate_network_set('output triplet set')
    # print(itps.consistency_score(otps))
