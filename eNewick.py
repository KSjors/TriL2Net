import logging
import sys
from data.reference_networks import TRINET_LIST
from config import settings
from datastructures.rooted_level_k_network import NetworkSet
from datastructures.solver import Solver
import json

if __name__ == '__main__':
    a = sys.argv
    input_file = a[1]
    try:
        output_file = a[2]
    except IndexError:
        output_file = 'output'
    try:
        parameter_file = a[3]
        with open(parameter_file) as file:
            parameters = json.load(file)
    except IndexError:
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

    root = logging.getLogger('')
    root.setLevel(logging.INFO)

    # sh = logging.StreamHandler(sys.stdout)
    # sh.setLevel(logging.INFO)
    # sformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # sh.setFormatter(sformatter)
    # root.addHandler(sh)

    with open(f"{output_file}.log", 'w+') as _:
        pass

    fh = logging.FileHandler(f"{output_file}.log")
    fh.setLevel(logging.INFO)
    fformatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(fformatter)
    root.addHandler(fh)

    TS_input = NetworkSet.from_enewick_format(input_file)
    solver = Solver(TRINET_LIST, TS_input, **parameters)
    network, _, _ = solver.solve()

    with open(f"{output_file}.eNewick", 'w+') as f:
        f.write(network.enewick())
