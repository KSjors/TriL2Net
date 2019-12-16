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

from data.all_trinets import get_standard_networks
from utils.help_functions import *
from tqdm import tqdm
import matplotlib.pyplot as plt
from experiments.experiments import NetworkGenerator, TrinetSetGenerator, SolverGenerator, ExperimentGenerator, Experiment, ExperimentSet
from config import settings
from datastructures.rooted_level_k_network import NetworkSet, RootedLevelKNetwork

if __name__ == '__main__':
    _, biconnected_trinet_binet_list, all_trinet_list = get_standard_networks()

    # Network Generator
    num_leaves = [10, 15, 20]
    recomb_prob = [0.1]
    term_prob = [0]
    ng = NetworkGenerator(num_leaves, recomb_prob, term_prob)

    # Trinet Set Generator
    tail_move_prob = [0.15]
    uni_prob = [0]
    del_prob = [0]
    tsg = TrinetSetGenerator(tail_move_prob, uni_prob, del_prob, all_trinet_list)

    # Solver Generator
    parameters = {
        'cut_arc_set_count_method'          : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'minimal_sink_set_method'         : [settings.DEFAULT_SINK_SET]
        , 'leaf_locator_method'             : [settings.DEFAULT_ORDER]
        , 'level_threshold_method'          : [settings.DEFAULT_THRESHOLD]
        , 'level_count_method'              : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'generator_count_method'          : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'symmetric_sides_set_count_method': [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'leaf_order_count_method'         : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
        , 'leaf_order_method'               : [settings.MAXIMUM_MULTIPLICITY]  # , settings.WEIGHTED_SUM, settings.WEIGHTED_AVERAGE]
    }
    sg = SolverGenerator(all_trinet_list, **parameters)

    # Experiment
    # eg = ExperimentGenerator(ng, tsg, sg, name="Test", trinet_set_method=settings.ITERATIVE, max_processes=4)
    # es = eg.run_times(times=2)

    # # Analyse
    # es = ExperimentSet.load(r'C:\Users\sjors\PycharmProjects\TriLoNet-2\experiments\Test', 'Test-2')
    # print(es.experiments)
    # es.plot_consistency_scores()
    # es.rerun(all_trinet_list, settings.ITERATIVE, 1)
    # es[0].plot_run_times()
    # print(es[0].run_time)
    # e = Experiment.load(r'C:\Users\sjors\PycharmProjects\TriLoNet-2\experiments\Test\2019-12-14 161943-007493.txt', name="Check")
    # e.input_network.visualize()
    #
    exp = Experiment.load(r'C:\Users\sjors\PycharmProjects\TriLoNet-2\experiments\Test\2019-12-16 122342-883682.txt', name='Rerun')
    # # ene = exp['input network']
    exp.rerun(all_trinet_list, settings.ITERATIVE, 1)
    # print(exp.output_network.equal_structure(exp.input_network, equal_naming=True))
    # exp.visualize_input_output()