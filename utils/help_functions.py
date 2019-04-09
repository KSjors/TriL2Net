import itertools
import copy
import time
import random
import socket
import hashlib
import logging
from bidict import bidict
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


def is_symmetry(symmetrical_nodes, symmetry_check_list, list_of_lists):
    """Check if one of list of lists isomorphisms (def by symmetrical nodes) is in symmetry_check_list"""
    logging.debug("Checking if {} is a symmetry of one of {} through {}.".format(list_of_lists, symmetry_check_list, symmetrical_nodes))
    number_of_symmetries = len(symmetrical_nodes)
    symmetry_iterator = all_combinations(list(symmetrical_nodes), 1, number_of_symmetries+1)
    for symmetries in symmetry_iterator:
        current_symmetry_bidict = bidict()
        for node in symmetries:
            current_symmetry_bidict.put(node, symmetrical_nodes[node])
            current_symmetry_bidict.put(symmetrical_nodes[node], node)
        a = substitute(current_symmetry_bidict, list_of_lists)
        if a in symmetry_check_list:
            return True
    return False


def substitute(symmetries: bidict, list_of_lists):
    result = [[symmetries.get(i, i) for i in lst] for lst in list_of_lists]
    result.sort()
    return result


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
