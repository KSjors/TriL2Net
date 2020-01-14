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
from config import settings
from datastructures.rooted_level_k_network import NetworkSet, RootedLevelKNetwork
from datastructures.solver import Solver
import pandas as pd
import winsound

if __name__ == '__main__':
    # Processing power
    max_processes = 16

    # Experiment settings
    experiment_name = 'cooper_test_1'
    text_file = r'C:\Users\sjors\PycharmProjects\TriLoNet-2\data\TriLoNet Data\Bollyky.txt'

    # solver settings
    parameters = {
        'cut_arc_set_count_method'          : settings.MAXIMUM_MULTIPLICITY
        , 'minimal_sink_set_method'         : settings.EXPAND_FIRST_SCC
        , 'leaf_locator_method'             : settings.DEFAULT
        , 'level_threshold_method'          : settings.DEFAULT
        , 'level_count_method'              : settings.MAXIMUM_MULTIPLICITY
        , 'generator_count_method'          : settings.MAXIMUM_MULTIPLICITY
        , 'symmetric_sides_set_count_method': settings.MAXIMUM_MULTIPLICITY
        , 'leaf_order_count_method'         : settings.MAXIMUM_MULTIPLICITY
        , 'leaf_order_method'               : settings.DEFAULT
        , 'fill_gaps'                       : settings.FALSE
    }

    # Get input
    TS_input = NetworkSet.from_text_trinet_format(text_file)
    # print(TS_input.summary())
    # TP_input = NetworkSet.induced_tree_set_of_network_set(TS_input, progress_bar=True,  max_processes=1)

    # Create solver
    solver = Solver(TRINET_LIST, TS_input, **parameters)
    network_2a, _, _ = solver.solve()
    print(network_2a.enewick())
    network_2b, _, _ = solver.solve()
    print(network_2b.enewick())
    network_2c, _, _ = solver.solve()
    print(network_2c.enewick())
    network_2d, _, _ = solver.solve()
    print(network_2d.enewick())
    ab = network_2a.equal_structure(network_2b, equal_naming=True)[0]
    bc = network_2b.equal_structure(network_2c, equal_naming=True)[0]
    cd = network_2c.equal_structure(network_2d, equal_naming=True)[0]
    print(ab, bc, cd)

    # # load network 1
    # ene1 = '((HBVADW4A, (1#H1, ((HPBETNC, ((HPBADRC, (HPBCG, (HPBADRA, (HEHBVAYR, (HPBADR1CG)2#H2)25)24)23)31, (HBVADRM, HPBCGADR)22)21)30, (HBVADR4, 2#H2)20)19)18)29, (((HPBMUT, (HPBADW1, (HBVAYWMCG, (HPBHBVAA, (XXHEPAV, (XXHEPA, (HBVDNA)3#H3)17)16)15)14)13)28, (HPBADW2, ((HPBADWZCG, (3#H3, ((HVHEPB, (HUMPRECX, (HBVADW)4#H4)12)27, (HBVADW2, 4#H4)11)10)9)26, (HPBADW3, HPBADWZ)8)7)6)5)1#H1)0'
    # network_1 = RootedLevelKNetwork.from_enewick(ene1)
    # network_1.visualize(internal_node_labels=False, rankdir='LR', file_path=None, format='pdf')
    # TS_1 = NetworkSet.induced_strict_network_set(network_1, 3, progress_bar=True, max_processes=16)
    # TP_1 = NetworkSet.induced_tree_set_of_network_set(TS_1, progress_bar=True, max_processes=1)
    #
    # # load network 2
    # # ene2 = '((((((HPBADWZ, HPBADW3)31, (HPBADWZCG, ((((HVHEPB, (HBVADW)26#H3)27, (HUMPRECX, (HBVADW2, 26#H3)29)28)25, (HBVDNA)5#H1)15, HPBADW1)14)13)12, HPBADW2)11, (((((5#H1, XXHEPA)10, XXHEPAV)9, HPBHBVAA)8, HBVAYWMCG)7, HPBMUT)6)4)1#H0, (HBVADW4A, (1#H0, ((HPBADR1CG)17#H2, (((HBVADRM, HPBCGADR)30, (HPBADRC, (((HBVADR4, (HEHBVAYR, 17#H2)24)23, HPBADRA)22, HPBCG)21)20)19, HPBETNC)18)16)3)2)0'
    # # network_2 = RootedLevelKNetwork.from_enewick(ene2)
    # # network_2.visualize(internal_node_labels=False, rankdir='LR', file_path='HBP_tril2net', format='pdf')
    #
    # # Compute network sets
    # TS_2 = NetworkSet.induced_strict_network_set(network_2, 3, progress_bar=True, max_processes=16)
    # TP_2 = NetworkSet.induced_tree_set_of_network_set(TS_2, progress_bar=True, max_processes=1)
    # ene = network_2.enewick()
    # # print(ene)
    #
    # # Scores
    # c1 = TS_input.consistency_score(TS_1, method=settings.WEIGHTED_AVERAGE)
    # c3 = TP_input.consistency_score(TP_1, method=settings.WEIGHTED_AVERAGE)
    # c2 = TS_input.consistency_score(TS_2, method=settings.WEIGHTED_AVERAGE)
    # c4 = TP_input.consistency_score(TP_2, method=settings.WEIGHTED_AVERAGE)
    # print(c1, c2, c3, c4)
    #
    #









