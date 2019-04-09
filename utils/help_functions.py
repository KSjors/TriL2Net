import itertools
import copy
import time
import random
import socket
import hashlib
import logging
logging.basicConfig(level=logging.DEBUG)

def list_to_name(lst, name_list):
    return [name_list[x] for x in lst]


def list_of_list_to_name(lst_of_lst, name_list):
    return [list_to_name(lst, name_list) for lst in lst_of_lst]


def all_combinations(any_list, min_len, max_len, direction=1):
    if direction == 1:
        return itertools.chain.from_iterable(
            itertools.combinations_with_replacement(any_list, i)
            for i in range(min_len, max_len))
    else:
        itertools.chain.from_iterable(
            itertools.combinations_with_replacement(any_list, i)
            for i in range(min_len, max_len).__reversed__())


def leaf_name_iterator(min_len, max_len):
    alphabet = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    alphabet_iterator = all_combinations(alphabet, min_len, max_len + 1)
    return alphabet_iterator


def trinets_of_generator(generator):
    edges = generator.get_edges(leafless=True)
    number_of_generator_leaves = generator.number_of_leaves

    extra_leaves_iterator = itertools.combinations(edges, 3 - number_of_generator_leaves)
    trinets_gen_sides = []
    for extra_leaves in extra_leaves_iterator:
        current_trinet = copy.deepcopy(generator)
        for extra_leaf in extra_leaves:
            current_trinet.add_leaf_to_edge(extra_leaf[0], extra_leaf[1])
        trinets_gen_sides.append([current_trinet, generator, extra_leaves])

    if number_of_generator_leaves == 1:
        extra_leaves_iterator = itertools.combinations(edges, 1)
        for extra_leaves in extra_leaves_iterator:
            current_trinet = copy.deepcopy(generator)
            new_node_name, _ = current_trinet.add_leaf_to_edge(extra_leaves[0][0], extra_leaves[0][1])
            current_trinet.add_leaf_to_edge(new_node_name, extra_leaves[0][1])
            trinets_gen_sides.append([current_trinet, generator, (extra_leaves[0], copy.deepcopy(extra_leaves[0]))])

    return trinets_gen_sides


def guid():
    """Generate a universally unique ID."""
    logging.debug("Generating unique ID.")
    t = time.time() * 10 ** 3
    r = random.random() * 10 ** 17
    try:
        a = socket.gethostbyname(socket.gethostname())
    except:
        # if we can't get a network address, just imagine one
        a = random.random() * 10 ** 17
    data = str(t) + str(r) + str(a)
    uid = hashlib.md5(data.encode('utf-8')).hexdigest()
    return uid
