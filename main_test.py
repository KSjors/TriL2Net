import logging

# TODO
#  -- CORRECT INPUT --
#  LEVEL-1 input
#  -- TEST ALTERED INPUT --
#  Swap leaves in some trinets
#  Swap trinet for other trinet
#  -- PROOF --
#  Pseudocode interesting functions
#  Proof that for correct input output also correct
#  Proof auxilliary graph theorem
#  -- OPTIMIZE --
#  Pruning
#  Isomorphism checking --> separate code for trees?
#  -- RANDOM --
#  Do not add leaf which increases level
#  -- WRITE --


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

import os
from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkSet, Omega
from data.reference_networks import get_standard_binets_trinets, regenerate_standard_binets_trinets, pickle_save, pickle_read
from utils.help_functions import *
from data.generators import *
from config import settings

if __name__ == '__main__':

    n = RootedLevelKNetwork.random(10, 0.4, 1)
    n.visualize()
    ts = NetworkSet.induced_strict_network_set(n, 3)
    ts.save_to_file('jojoj', frmt=settings.FORMAT_tnets)
    ts2 = NetworkSet.from_named_trinet_format('jojoj.tnet')

    print(ts.consistency_score(ts2, method=settings.WEIGHTED_AVERAGE))
