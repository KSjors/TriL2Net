import logging
import itertools
import operator as op
from functools import reduce
import copy
import time
import random
import socket
import hashlib
import numpy as np
from mip import model
from bidict import bidict


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


def leaf_name_iterator(min_len, max_len, char_type='alph'):
    assert char_type in ('alph', 'ALPH', 'uid')
    alphabet = ""
    if char_type == 'alph':
        alphabet = list('abcdefghijklmnopqrstuvwxyz')
    elif char_type == 'ALPH':
        alphabet = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    elif char_type == 'uid':
        return guid_generator()
    alphabet_iterator = all_combinations(alphabet, min_len, max_len + 1)
    return alphabet_iterator


def is_symmetry(symmetrical_nodes, symmetry_check_list, list_of_lists):
    """Check if one of list of lists isomorphisms (def by symmetrical nodes) is in symmetry_check_list"""
    logging.debug("Checking if {} is a symmetry of one of {} through {}.".format(list_of_lists, symmetry_check_list, symmetrical_nodes))
    number_of_symmetries = len(symmetrical_nodes)
    symmetry_iterator = all_combinations(list(symmetrical_nodes), 1, number_of_symmetries + 1)
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
    logging.debug("Generated unique ID: {}.".format(uid))
    return uid

def guid_generator():
    while True:
        yield guid()


def leaf_names_to_identifier(leaf_names: list) -> tuple:
    alphabetical_leaf_names = sorted(copy.copy(leaf_names))
    return tuple(alphabetical_leaf_names)


def mss_leaf_name(mss):
    return '(' + "-".join(sorted(mss)) + ')'


def coalesce(value, default=None):
    return value if value is not None else default


def determine_ranking(score_matrix, ranking=None):
    score_matrix = copy.copy(score_matrix)
    if ranking is None:
        ranking = []
    if sum(sum(score_matrix)) == 0:
        all_players = list(range(len(score_matrix)))
        for player in ranking:
            all_players.remove(player)
        return ranking + all_players
    player_wins = sum(score_matrix.T)
    best_players = np.where(player_wins == player_wins.max())[0]
    if len(best_players) == 1 or len(best_players) == len(score_matrix) - len(ranking):
        # TODO all same ranking in one go?
        best_player = best_players[0]
    else:
        intermediate_ranking = determine_ranking(score_matrix[best_players][:, best_players])
        best_player = best_players[intermediate_ranking[0]]
    score_matrix[best_player, :] = 0
    score_matrix[:, best_player] = 0
    ranking.append(best_player)
    return determine_ranking(score_matrix, ranking)


def enewick(string):
    result = []
    while len(string) > 0:
        try:
            left_bracket_index = string.index("(")
        except ValueError:
            left_bracket_index = len(string) - 1
        try:
            right_bracket_index = string.index(")")
        except ValueError:
            right_bracket_index = len(string) - 1
        try:
            comma_index = string.index(",")
        except ValueError:
            comma_index = len(string) - 1
        try:
            space_index = string.index(" ")
        except ValueError:
            space_index = len(string) - 1
        m = min([left_bracket_index, right_bracket_index, comma_index, space_index])
        if len(string[:m]) > 0:
            result.append(string[:m])
        result.append(string[m])
        string = string[m + 1:]

    adjacency_dict = dict()
    leaf_names = set()
    enewick_helper(np.array(result), adjacency_dict, leaf_names)
    return adjacency_dict, leaf_names


def enewick_helper(string, adjacency_dict, leaf_names):
    if len(string) == 1:
        leaf_names.add(string[0])
        return

    root = string[-1]
    string = string[1:-2]

    left_bracket_count = (string == '(').astype(int)
    right_bracket_count = (string == ')').astype(int)
    cumsum = np.cumsum(left_bracket_count - right_bracket_count)
    spaces = set(np.where(string == " ")[0])
    zeroes = set(np.where(cumsum == 0)[0])
    middle = spaces.intersection(zeroes)

    if len(middle) != 0:
        middle = middle.pop()
        left_side = string[:middle - 1]
        right_side = string[middle + 1:]
        try:
            adjacency_dict[root].extend([left_side[-1], right_side[-1]])
        except KeyError:
            adjacency_dict[root] = [left_side[-1], right_side[-1]]
        enewick_helper(left_side, adjacency_dict, leaf_names)
        enewick_helper(right_side, adjacency_dict, leaf_names)
    else:
        try:
            adjacency_dict[root].extend([string[-1]])
        except KeyError:
            adjacency_dict[root] = [string[-1]]
        enewick_helper(string, adjacency_dict, leaf_names)


def check_bidict(bdict: bidict, first, second) -> bool:
    try:
        if bdict[first] != second:
            return False
    except KeyError:
        pass
    try:
        if bdict.inverse[second] != first:
            return False
    except KeyError:
        pass
    return True


def ncr(n, r):
    r = min(r, n - r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer / denom


def simplex_to_ILP(c, A_eq, b_eq, sense="maximize"):
    if sense == "maximize":
        m = model.Model(sense=model.MAXIMIZE, solver_name=model.CBC)
    else:
        m = model.Model(sense=model.MINIMIZE, solver_name=model.CBC)
    m.verbose = 0
    number_of_variables = A_eq.shape[1]
    y = [m.add_var(var_type=model.BINARY) for i in range(number_of_variables)]
    for row, b in zip(A_eq, b_eq):
        m += model.xsum(row[i] * y[i] for i in range(number_of_variables)) == b
    m.objective = model.xsum(c[i] * y[i] for i in range(number_of_variables))
    m.max_mip_gap_abs = 0.05
    status = m.optimize(max_seconds=300)
    assert status == model.OptimizationStatus.OPTIMAL or status == model.OptimizationStatus.FEASIBLE, "Could not find feasible solution"
    solution = [v.x for v in m.vars]
    return solution, m.objective_value, status


def shifted_node_number(node_number, other_node_numbers):
    return node_number - sum([other_node_number < node_number for other_node_number in other_node_numbers])



