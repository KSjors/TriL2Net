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

from data.reference_networks import TRINET_LIST
from data.generators import GENERATOR_LEVEL_LIST
from config import settings
from datastructures.rooted_level_k_network import NetworkSet, RootedLevelKNetwork
from datastructures.solver import Solver

if __name__ == '__main__':
    # SUBNET SETUP
    subnet_level = 1
    number_of_subnets = 10
    number_of_leaves_per_internal_arc = 4

    # NETWORK SETUP
    number_of_species = 20

    # DISTORTION SETUP
    tail_move_fraction = 0
    uniform_noise_fraction = 0.1

    # Create network
    generators = GENERATOR_LEVEL_LIST[subnet_level]
    subnets = NetworkSet.subnets_from_generators(generators, number_of_leaves_per_internal_arc)
    input_network = RootedLevelKNetwork.from_subnets(subnets, number_of_subnets, number_of_species)
    input_network.visualize()

    # Create and distort trinet set
    trinet_set = NetworkSet.induced_strict_network_set(input_network, 3, 8)
    distorted_trinet_set = NetworkSet.distort(trinet_set, uniform_noise_fraction, tail_move_fraction, 0, TRINET_LIST, subnet_level)

    # Create solver and solve
    solver = Solver(TRINET_LIST, distorted_trinet_set)
    output_network = solver.solve()[0]
    output_network.visualize()














