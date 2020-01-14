import logging
import itertools
import operator as op
from functools import reduce
import collections
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
    assert char_type in ('alph', 'ALPH', 'uid'), "Unknown char type"
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
    # TODO, if default value takes long to calculate --> ....
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
        if left_bracket_index == right_bracket_index == comma_index == space_index:
            result.append(string)
            string = []
        else:
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
        leaf_name = string[0]
        if '#' not in leaf_name:
            leaf_names.add(leaf_name)
        return

    root = string[-1]
    root_name = string[-1] if '#' not in string[-1] else string[-1][: string[-1].index('#')]
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
        left_side_name = left_side[-1] if '#' not in left_side[-1] else left_side[-1][: left_side[-1].index('#')]
        right_side_name = right_side[-1] if '#' not in right_side[-1] else right_side[-1][: right_side[-1].index('#')]
        try:
            adjacency_dict[root_name].extend([left_side_name, right_side_name])
        except KeyError:
            adjacency_dict[root_name] = [left_side_name, right_side_name]
        enewick_helper(left_side, adjacency_dict, leaf_names)
        enewick_helper(right_side, adjacency_dict, leaf_names)
    else:
        string_name = string[-1] if '#' not in string[-1] else string[-1][: string[-1].index('#')]
        try:
            adjacency_dict[root_name].appenmd(string_name)
        except KeyError:
            adjacency_dict[root_name] = [string_name]
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


def original_node_number(node_number, other_node_numbers):
    while True:
        node_numbers_below = [other_node_number for other_node_number in other_node_numbers if other_node_number <= node_number]
        other_node_numbers = [other_node_number for other_node_number in other_node_numbers if other_node_number > node_number]
        if len(node_numbers_below) >= 1:
            node_number += len(node_numbers_below)
        else:
            break
    return node_number


def original_node_numbers(node_numbers):
    result = []
    for index, node_number in enumerate(node_numbers):
        result.append(original_node_number(node_number, node_numbers[:index]))
    return result


def power_set(iterable):
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s) + 1))


def flatten_dictionary(d, parent_key='', sep=' '):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten_dictionary(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def chebyshev_interpolation_points(start, stop, n, round_method=None):
    res = [0.5 * (start + stop) + 0.5 * (stop - start) * np.cos((2 * k + 1) / (2 * n) * np.pi) for k in range(n)]
    if round_method:
        res = [round_method(x) for x in res]
    return sorted(res)
