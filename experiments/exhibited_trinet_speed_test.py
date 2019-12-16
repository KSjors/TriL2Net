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

from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkSet

if __name__ == "__main__":
    network = RootedLevelKNetwork.random(50, 0, 2)

    trinets = NetworkSet.induced_strict_network_set(network, 3, 4, True, method='Recursive')
    trinets_2 = NetworkSet.induced_strict_network_set(network, 3, 4, True, method='Iterative')
