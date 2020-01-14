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

from data.reference_networks import get_standard_binets_trinets, regenerate_standard_binets_trinets
from data.generators import GENERATOR_DICT
from utils.help_functions import *
from tqdm import tqdm
import matplotlib.pyplot as plt
from experiments.experiments import TrinetSetGenerator, SolverGenerator, ExperimentGenerator, ExperimentGenerator2, Experiment, ExperimentSet, \
    NetworkGeneratorSubnets
from config import settings
from datastructures.rooted_level_k_network import NetworkSet, RootedLevelKNetwork
import pandas as pd
import winsound

if __name__ == '__main__':
    # regenerate_standard_binets_trinets()
    # kk

    # Processing power
    max_processes = 16

    # ---------------------------------------------------------------------- 2
    # Experiment settings
    experiment_name = 'runtimes\\level_2\\FULL'
    repeat = 4

    generators = [GENERATOR_DICT['2a'], GENERATOR_DICT['2b'], GENERATOR_DICT['2c'], GENERATOR_DICT['2d']]
    num_leaves = [15, 20, 25, 30, 35, 40]
    subnets = NetworkSet.subnets_from_generators(generators, 2)
    ng = NetworkGeneratorSubnets(num_leaves, subnets)

    # Trinet Set Generator
    tail_move_prob = [0.025]
    uni_prob = [0.025]
    del_prob = [0]
    max_replacement_level = 2
    tsg = TrinetSetGenerator(tail_move_prob, uni_prob, del_prob, max_replacement_level)

    # ------------------------- MM
    parameters = {
        'cut_arc_set_count_method'          : [settings.MAXIMUM_MULTIPLICITY]
        , 'minimal_sink_set_method'         : [settings.EXPAND_FIRST_SCC]
        , 'leaf_locator_method'             : [settings.DEFAULT]
        , 'level_threshold_method'          : [settings.DEFAULT]
        , 'level_count_method'              : [settings.WEIGHTED_AVERAGE]
        , 'generator_count_method'          : [settings.WEIGHTED_AVERAGE]
        , 'symmetric_sides_set_count_method': [settings.WEIGHTED_AVERAGE]
        , 'leaf_order_count_method'         : [settings.WEIGHTED_AVERAGE]
        , 'leaf_order_method'               : [settings.DEFAULT]
        , 'fill_gaps'                       : [0]
    }
    sg = SolverGenerator(**parameters)
    eg = ExperimentGenerator(ng, tsg, sg, name=experiment_name, max_processes=max_processes, consistency_method=settings.WEIGHTED_AVERAGE)
    eg.run_times_2(times=repeat, reuse_network=False, reuse_trinet_set=False, compute=[0, 0, 0, 0, 0, 0])
