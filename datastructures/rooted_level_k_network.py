# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 10:32:22 2019

@author: Sjors
"""
import logging
import numpy as np
import pandas as pd
import copy
import itertools
from utils.help_functions import coalesce
from tqdm import tqdm
from bidict import bidict
from graphviz import Digraph
from utils.help_functions import *
from tarjan import tarjan
import math
import time
import os
import pprint
import random
import scipy
import collections
import multiprocessing

pp = pprint.PrettyPrinter(indent=4)

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'


class RootedLevelKNetwork:
    def __init__(self, adj_matrix: np.ndarray, node_name_map: bidict, leaf_names: list = None, leaf_numbers: list = None, level: int = 2, dimension: int = 3):
        self.uid = guid()
        self.logger = logging.getLogger('network.{}'.format(self.uid))
        self.logger.debug("Creating new rooted level {} network of dimension {}.".format(level, dimension))
        shape = adj_matrix.shape
        assert len(node_name_map) == shape[0], "Number of names does not match number of nodes."
        assert shape[0] == shape[1], "Adjacency matrix is not square."

        self.adj_matrix = adj_matrix
        self.node_name_map = node_name_map

        if leaf_numbers is not None:
            self.leaf_numbers = leaf_numbers
        else:
            self.leaf_numbers = [node_name_map[leaf_name] for leaf_name in leaf_names]

        # Set shape numbers
        self.number_of_nodes = shape[1]

        # Set network properties
        self.level = level
        self.dimension = dimension

        # optimization variables
        self._o_cut_arc_matrix = None
        self._o_cut_arc_sets = None
        self._o_biconnected_components = None
        self._o_partial_ordering = None
        self._o_leaf_ordering = None

    # ---------------------------------------------------------------- CLASS METHODS ---------------------------------------------------------------- #
    def __copy__(self):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.adj_matrix = copy.copy(self.adj_matrix)
        cp.node_name_map = copy.copy(self.node_name_map)
        cp.leaf_numbers = copy.copy(self.leaf_numbers)
        cp.number_of_nodes = copy.copy(self.number_of_nodes)
        cp.level = copy.copy(self.level)
        cp.dimension = copy.copy(self.dimension)
        cp.uid = guid()
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp._o_cut_arc_matrix = copy.copy(self._o_cut_arc_matrix)
        cp._o_cut_arc_sets = copy.copy(self._o_cut_arc_sets)
        cp._o_biconnected_components = copy.copy(self._o_biconnected_components)
        cp._o_partial_ordering = copy.copy(self._o_partial_ordering)
        cp._o_leaf_ordering = copy.copy(self._o_leaf_ordering)
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.adj_matrix = copy.copy(self.adj_matrix)
        cp.node_name_map = copy.copy(self.node_name_map)
        cp.leaf_numbers = copy.copy(self.leaf_numbers)
        cp.number_of_nodes = copy.copy(self.number_of_nodes)
        cp.level = copy.copy(self.level)
        cp.dimension = copy.copy(self.dimension)
        cp.uid = guid()
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp._o_cut_arc_matrix = copy.deepcopy(self._o_cut_arc_matrix)
        cp._o_cut_arc_sets = copy.deepcopy(self._o_cut_arc_sets)
        cp._o_biconnected_components = copy.deepcopy(self._o_biconnected_components)
        cp._o_partial_ordering = copy.deepcopy(self._o_partial_ordering)
        cp._o_leaf_ordering = copy.deepcopy(self._o_leaf_ordering)
        return cp

    @classmethod
    def from_enewick(cls, string, check_valid=True):
        adjacency_dict, leaf_names = enewick(string)
        result = cls.from_connections_dict(adjacency_dict, leaf_names=list(leaf_names), check_valid=check_valid)
        return result

    @classmethod
    def restrict(cls, network, leaf_names: list, suppress_redundant='all', suppress_parallel=True):
        return cls._restrict(network, network.get_node_numbers_of(leaf_names), suppress_redundant, suppress_parallel)

    @classmethod
    def _restrict(cls, network, leaf_numbers: list, suppress_redundant='all', suppress_parallel=True):
        assert set(leaf_numbers).issubset(
            set(network.leaf_numbers)), "Can not create sub-network of network {} using {} as they are not leaves of network.".format(
            network.uid,
            leaf_numbers)
        new_network = copy.deepcopy(network)
        print(network.get_node_names_of(leaf_numbers))
        new_network._terminate_leaves(leaf_numbers_to_keep=leaf_numbers)
        new_network.prune(suppress_redundant=suppress_redundant, suppress_parallel=suppress_parallel)
        return new_network

    @classmethod
    def from_dir_adj_matrix(cls, dir_adj_matrix: np.ndarray, level=2, dimension=2, check_valid=True, char_type='alph'):
        shape = dir_adj_matrix.shape
        adj_matrix = np.zeros((shape[1], shape[1])).astype(float)
        adj_matrix[:shape[0], :shape[1]] += dir_adj_matrix
        adj_matrix[:shape[1], :shape[0]] -= dir_adj_matrix.T
        node_name_map = bidict()
        leaf_names = list()
        ln_iterator = leaf_name_iterator(1, 100, char_type=char_type)
        for i in range(shape[1]):
            if i < shape[0]:
                node_name_map.put(str(i), i)
            else:
                leaf_name = "".join(next(ln_iterator))
                node_name_map.put(leaf_name, i)
                leaf_names.append(leaf_name)
        network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_names=leaf_names, level=level, dimension=dimension)
        assert (not check_valid) or network.is_valid(), "Connection dictionary results in a invalid network."
        return network

    @classmethod
    def from_connections_dict(cls, connections_dict: dict, leaf_names: list = None, level=2, dimension=2, check_valid=True):
        node_name_map = bidict()
        leaf_names = coalesce(leaf_names, [])

        # Empty adjacency matrix
        adj_matrix = np.zeros((0, 0)).astype(int)
        network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_names=[], level=level, dimension=dimension)

        for from_node_name in connections_dict.keys():
            network._add_node_to_network(str(from_node_name))

        for to_node_names in connections_dict.values():
            for to_node_name in to_node_names:
                if str(to_node_name) not in network.node_name_map:
                    network._add_node_to_network(str(to_node_name), leaf=True)

        for from_node_name, to_node_names in connections_dict.items():
            from_node_number = network.get_node_number_of(str(from_node_name))
            for to_node_name in to_node_names:
                to_node_number = network.get_node_number_of(str(to_node_name))
                network._add_arc(from_node_number, to_node_number)

        for leaf_name in leaf_names:
            if leaf_name not in network.node_name_map:
                network._add_node_to_network(str(leaf_name), leaf=True)
        assert (not check_valid) or network.is_valid(), "Connection dictionary results in an invalid network."
        return network

    @classmethod
    def get_network_below_node(cls, network, node_name: str):
        return RootedLevelKNetwork._get_network_below_node(network, network.get_node_number_of(node_name))

    @classmethod
    def _get_network_below_node(cls, network, node_number: int):
        node_numbers_below_node = sorted(list(network._nodes_below_nodes({node_number})))
        node_names_below_node = network.get_node_names_of(node_numbers_below_node)
        node_name_map = bidict({name: number for number, name in enumerate(node_names_below_node)})
        adj_matrix = network.adj_matrix[:, node_numbers_below_node][node_numbers_below_node, :]
        leaf_names = [name for name in node_names_below_node if network.is_leaf_node(name)]
        new_network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_names=leaf_names, level=network.level, dimension=network.dimension)
        return new_network

    @classmethod
    def random(cls, number_of_leaves: int, recombination_chance: int, level: int):
        enewick_string = '(a, b)0'
        network = cls.from_enewick(enewick_string)
        network.level = level
        network.evolve_times(number_of_leaves-2, recombination_chance, progress_bar=False)
        network.standardize_node_names()
        return network

    def is_valid(self) -> bool:
        """Check if network is valid: has right degrees and number of roots and reticulations"""
        self.logger.debug("Checking validity.")
        in_degrees = np.array(self.get_in_degrees())
        out_degrees = np.array(self.get_out_degrees())

        number_of_roots = sum(in_degrees == 0)
        assert number_of_roots == 1, "Network {} is not valid, it has more than one root ({}).".format(self.uid, number_of_roots)
        assert sum(in_degrees <= self.dimension - 1), "Network {} is not valid, it has a node with more arcs entering than allowed.".format(self.uid)
        assert sum(out_degrees <= self.dimension - 1), "Network {} is not valid, it has a node with more arcs exiting than allowed.".format(self.uid)
        assert sum(
            (in_degrees + out_degrees) <= self.dimension), "Network {} is not valid, it has a node with more arcs entering and exiting than allowed.".format(
            self.uid)
        # TODO: check if leaves in the end of the matrix
        # TODO: acyclic
        # TODO: root can reach everything

        total_number_of_reticulations = sum(in_degrees > 1)
        if total_number_of_reticulations <= self.level:
            return True
        biconnected_components = self._biconnected_components
        for bc in biconnected_components:
            in_sum_bc = in_degrees[bc]
            number_of_reticulations = sum(in_sum_bc > 1)
            assert number_of_reticulations <= self.level, "Network {} is not valid, biconnected component {} has to many reticulations ({} > {}).".format(
                self.uid, bc, number_of_reticulations, self.level)
        return True

    # --------------------------------------------------------- NODE NAME <--> NODE NUMBER METHODS --------------------------------------------------------- #
    def get_node_names_of(self, node_number_list: list) -> list:
        """Retrieve names corresponding to node numbers."""
        # return [self.get_node_name_of(node_number) for node_number in node_number_list]
        return list(map(self.get_node_name_of, node_number_list))

    def get_node_numbers_of(self, node_name_list: list) -> list:
        """Retrieve numbers corresponding to node names."""
        # return [self.get_node_number_of(node_name) for node_name in node_name_list]
        return list(map(self.get_node_number_of, node_name_list))

    def get_node_name_of(self, node_number: int) -> str:
        """Retrieve name corresponding to node number"""
        return self.node_name_map.inverse[node_number]

    def get_node_number_of(self, node_name: str) -> int:
        """Retrieve number corresponding to node name."""
        return self.node_name_map[node_name]

    # --------------------------------------------------------- PROPERTIES --------------------------------------------------------- #
    @property
    def node_names(self) -> list:
        """Retrieve list of node names"""
        return self.get_node_names_of(list(self.node_name_map.values()))

    @property
    def leaf_names(self) -> list:
        """Retrieve list of leaf names"""
        return self.get_node_names_of(self.leaf_numbers)

    @property
    def internal_node_numbers(self) -> list:
        """Retrieve list of leaf names ordered by their number"""
        return [node_number for node_number in self.node_name_map.values() if not self._is_leaf_node(node_number)]

    @property
    def internal_node_names(self) -> list:
        """Retrieve list of leaf names ordered by their number"""
        return self.get_node_names_of([node_number for node_number in self.node_name_map.values() if not self._is_leaf_node(node_number)])

    @property
    def leaf_mask(self):
        mask = np.ones(self.number_of_nodes)
        mask[self.leaf_numbers] = 0
        return mask

    @property
    def root(self) -> str:
        """Retrieve the name of the root."""
        return self.get_node_name_of(self._root)

    @property
    def _root(self) -> int:
        return np.where(sum(self.adj_matrix > 0) == 0)[0][0]

    @property
    def reticulations(self) -> list:
        """Retrieve node names of all reticulations."""
        return self.get_node_names_of(self._reticulations)

    @property
    def _reticulations(self) -> list:
        return list(np.where(np.array(self.get_in_degrees()) > 1)[0])

    @property
    def biconnected_components(self) -> list:
        """Compute biconnected components."""
        return [self.get_node_names_of(biconnected_component) for biconnected_component in self._biconnected_components]

    @property
    def _biconnected_components(self) -> list:
        """Compute biconnected components."""
        if self._o_biconnected_components is None:
            self._compute_biconnected_components()

        return [[node_number for node_number in component if not self._is_leaf_node(node_number)]
                for component in self._o_biconnected_components if not (len(component) == 1 and self._is_leaf_node(component[0]))]

    @property
    def biconnected(self) -> bool:
        """Check if network is biconnected"""
        return len(self._biconnected_components) == 1

    @property
    def cut_arc_sets(self) -> list:
        """Retrieve cut-arc sets of network."""
        return [self.get_node_names_of(cut_arc_set) for cut_arc_set in self._cut_arc_sets]

    @property
    def _cut_arc_sets(self) -> list:
        if self._o_cut_arc_sets is None:
            self._compute_cut_arc_sets()
            cut_arc_matrix = self.cut_arc_matrix
            _, to_nodes = np.where(cut_arc_matrix == 1)
            self._o_cut_arc_sets = []
            for to_node in to_nodes:
                ca_set = list(self._nodes_below_nodes({to_node}))
                self._o_cut_arc_sets.append(ca_set)
        return [[node_number for node_number in cut_arc_set if self._is_leaf_node(node_number)] for cut_arc_set in self._o_cut_arc_sets]

    @property
    def cut_arc_matrix(self) -> np.ndarray:
        """Compute indicator matrix for arcs which are cut-arcs."""
        if self._o_cut_arc_matrix is None:
            self._compute_cut_arc_matrix()

        return self._o_cut_arc_matrix

    @property
    def partial_ordering(self) -> list:
        return [self.get_node_names_of(partial_ordering) for partial_ordering in self._partial_ordering]

    @property
    def _partial_ordering(self) -> list:
        if self._o_partial_ordering is None:
            self._compute_partial_ordering()
        return self._o_partial_ordering

    @property
    def leaf_ordering(self) -> list:
        return [self.get_node_names_of(leaf_ordering) for leaf_ordering in self._leaf_ordering]

    @property
    def _leaf_ordering(self):
        if self._o_leaf_ordering is None:
            self._compute_leaf_ordering()
        return self._o_leaf_ordering

    @property
    def number_of_edges(self):
        mask = self.adj_matrix > 0
        return np.sum(np.multiply(self.adj_matrix, mask))

    @property
    def edges(self) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        return [(self.get_node_name_of(from_node), self.get_node_name_of(to_node)) for from_node, to_node in self._edges]

    @property
    def _edges(self) -> list:
        i = 1
        from_nodes, to_nodes = [], []
        while np.any(self.adj_matrix >= i):
            extra_from_nodes, extra_to_nodes = np.where(self.adj_matrix >= i)
            to_nodes.extend(extra_to_nodes)
            from_nodes.extend(extra_from_nodes)
            i += 1
        return [(from_node, to_node) for from_node, to_node in zip(np.array(from_nodes), np.array(to_nodes))]

    @property
    def internal_edges(self) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        return [(self.get_node_name_of(from_node), self.get_node_name_of(to_node)) for from_node, to_node in self._internal_edges]

    @property
    def _internal_edges(self) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        mask = self.leaf_mask
        i = 1
        from_nodes, to_nodes = [], []
        while np.any(np.multiply(self.adj_matrix, mask) >= i):
            extra_from_nodes, extra_to_nodes = np.where(np.multiply(self.adj_matrix, mask) >= i)
            to_nodes.extend(extra_to_nodes)
            from_nodes.extend(extra_from_nodes)
            i += 1
        return [(from_node, to_node) for from_node, to_node in zip(np.array(from_nodes), np.array(to_nodes))]

    @property
    def number_of_reticulations(self):
        col_sum = sum(self.adj_matrix > 0)
        if type(col_sum) == int:
            col_sum = np.array([col_sum])
        number_of_reticulations = sum(col_sum > 1)
        return number_of_reticulations

    @property
    def number_of_leaves(self):
        return len(self.leaf_numbers)

    def strict_level(self) -> int:
        in_degrees = np.array(self.get_in_degrees())

        level = 0
        biconnected_components = self._biconnected_components
        for bc in biconnected_components:
            in_sum_bc = in_degrees[bc]
            level = max(level, sum(in_sum_bc > 1))

        return level

    # --------------------------------------------------------- NODE TYPES --------------------------------------------------------- #
    def is_connected_node(self, node_name: str) -> bool:
        """Check whether node is connected to the graph."""
        return self._is_connected_node(self.get_node_number_of(node_name))

    def _is_connected_node(self, node_number: int) -> bool:
        """Check whether node is connected to the graph."""
        return sum(self.adj_matrix[node_number, :] != 0) > 0

    # --------------------------------------------------------- ALTER NODE NAME METHODS --------------------------------------------------------- #
    def rename_node(self, old_name: str, new_name: str) -> None:
        """Replace name of node with name old_name by new_name."""
        assert new_name not in self.node_names, f"There already exists a node with node name {new_name}."

        # Rename node in node_names
        node_number = self.node_name_map.pop(old_name)
        self.node_name_map.put(new_name, node_number)

    def _rename_node(self, node_number: int, new_name: str) -> None:
        assert new_name not in self.node_names, f"There already exists a node with node name {new_name}."

        # Rename node in node_names
        old_name = self.node_name_map.inverse[node_number]
        self.node_name_map.pop(old_name)
        self.node_name_map.put(new_name, node_number)

    def first_unused_internal_name(self) -> str:
        """Find first unused leaf name in numerical order."""
        i = 0
        while str(i) in self.node_name_map:
            i += 1
        return str(i)

    def first_unused_leaf_name(self, char_type='alph') -> str:
        """Find first unused leaf name in alphabetical order."""
        leaf_name_length = 1
        found = False
        leaf_name = None
        while not found:
            leaf_names = leaf_name_iterator(leaf_name_length, leaf_name_length, char_type)
            while True:
                try:
                    leaf_name = "".join(next(leaf_names))
                except StopIteration:
                    leaf_name_length += 1
                    break
                try:
                    if not self.is_leaf_node(leaf_name):
                        found = True
                        break
                except KeyError:
                    found = True
                    break
        assert leaf_name is not None, "Could not find a leaf name."
        return leaf_name

    def standardize_internal_node_names(self) -> None:
        for internal_node in self.internal_node_numbers:
            self._rename_node(internal_node, self.get_node_name_of(internal_node) + "^")

        for index, internal_node in enumerate(self.internal_node_numbers):
            self._rename_node(internal_node, str(index))

    def standardize_leaf_node_names(self, char_type='alph') -> None:
        for leaf_number in self.leaf_numbers:
            self._rename_node(leaf_number, self.get_node_name_of(leaf_number) + "^")

        for leaf_node in self.leaf_numbers:
            self._rename_node(leaf_node, self.first_unused_leaf_name(char_type))

    def standardize_node_names(self, char_type='alph') -> None:
        for node_number in range(self.number_of_nodes):
            self._rename_node(node_number, self.get_node_name_of(node_number) + "^")

        for index, internal_node in enumerate(self.internal_node_numbers):
            self._rename_node(internal_node, str(index))

        for leaf_node in self.leaf_numbers:
            self._rename_node(leaf_node, self.first_unused_leaf_name(char_type))

    # TODO: to networkinfolist
    def collapse(self, leaf_names_set) -> list:
        """Returns list of shrunken networks"""
        intersection = set(self.leaf_names).intersection(leaf_names_set)
        if len(intersection) == len(self.leaf_names):
            return []
        if len(intersection) == 0:
            return [copy.deepcopy(self)]
        leaf_numbers_set = set(self.get_node_numbers_of(list(intersection)))
        mss_name = mss_leaf_name(leaf_names_set)
        return self._collapse(leaf_numbers_set, mss_name)

    def _collapse(self, leaf_numbers_set: set, new_name: str) -> list:
        result = []
        for leaf_to_keep in leaf_numbers_set:
            new_network = copy.deepcopy(self)
            new_network._rename_node(leaf_to_keep, new_name)
            for leaf_number in leaf_numbers_set.difference([leaf_to_keep]):
                new_network._terminate_leaf(leaf_number)
            new_network.prune()

            result.append(new_network)
        return result

    # --------------------------------------------------------- ALTER EDGE METHODS --------------------------------------------------------- #
    def is_arc(self, from_node_name: str, to_node_name: str):
        return self._is_arc(self.get_node_number_of(from_node_name), self.get_node_number_of(to_node_name))

    def _is_arc(self, from_node_number: int, to_node_number: int):
        """Add connection between from_node_name and to_node_name."""
        return self.adj_matrix[from_node_number][to_node_number] >= 1

    def _add_arc(self, from_node_number: int, to_node_number: int):
        """Add connection between from_node_name and to_node_name."""
        self.adj_matrix[from_node_number][to_node_number] += 1
        self.adj_matrix[to_node_number][from_node_number] -= 1

        self.reset_optimization_variables()

    def _remove_arc(self, from_node_number: int, to_node_number: int):
        """Remove connection between from_node_name and to_node_name."""
        assert self.adj_matrix[from_node_number][to_node_number] > 0 > self.adj_matrix[to_node_number][
            from_node_number], "Cannot remove non-existent connection"

        self.adj_matrix[from_node_number][to_node_number] -= 1
        self.adj_matrix[to_node_number][from_node_number] += 1

        self.reset_optimization_variables()

    # --------------------------------------------------------- ALTER NODE METHODS --------------------------------------------------------- #
    def _remove_leaf_status_node(self, node_number: int):
        """Remove leaf with name node_name"""
        self.leaf_numbers.remove(node_number)

    def _add_leaf_status_node(self, node_number: int):
        """Remove leaf with name node_name"""
        self.leaf_numbers.append(node_number)

    def _remove_node_from_dict(self, node_number: int) -> None:
        """Remove node with name node_name from node_names dictionary."""
        self.node_name_map.inverse.pop(node_number)

        # Decrease node number of nodes with number higher than node_name's
        for y in range(node_number + 1, self.number_of_nodes):
            self.node_name_map[self.node_name_map.inverse[y]] -= 1
            if self._is_leaf_node(y):
                self.leaf_numbers.remove(y)
                self.leaf_numbers.append(y - 1)

    def _remove_node_from_adj_matrix(self, node_number: int):
        """Remove node from adj_matrix."""
        mask = np.ones(self.number_of_nodes, dtype=bool)
        mask[node_number] = False
        self.adj_matrix = self.adj_matrix[mask, :][:, mask]

    def _remove_node_from_network(self, node_number: int):
        """Remove node fully from network."""
        if self._is_leaf_node(node_number):
            self._remove_leaf_status_node(node_number)

        self._remove_node_from_adj_matrix(node_number)
        self._remove_node_from_dict(node_number)
        self.number_of_nodes -= 1

        self.reset_optimization_variables()

    def _add_node_to_dict(self, node_name: str = None, leaf=False, char_type='alph') -> (str, int):
        """Add internal node to node_names dictionary. Returns its number."""

        # Get next number, and get a node_name if no name is given
        node_number = self.number_of_nodes
        if leaf:
            node_name = self.first_unused_leaf_name(char_type) if node_name is None else node_name
        else:
            node_name = self.first_unused_internal_name() if node_name is None else node_name

        self.node_name_map.put(node_name, node_number)
        return node_name, node_number

    def _add_node_to_adj_matrix(self) -> None:
        """Add internal node to adj matrix """

        v = np.zeros((1, self.number_of_nodes))
        self.adj_matrix = np.concatenate((self.adj_matrix, v), axis=0)

        h = np.zeros((self.number_of_nodes + 1, 1))
        self.adj_matrix = np.concatenate((self.adj_matrix, h), axis=1)

    def _add_node_to_network(self, node_name: str = None, leaf: bool = False, char_type='alph') -> (str, int):
        """Add internal node to network."""
        self.logger.debug("Adding internal node to network.")

        node_name, node_number = self._add_node_to_dict(node_name, leaf=leaf, char_type=char_type)
        self._add_node_to_adj_matrix()

        # Update numbers
        if leaf:
            self.leaf_numbers.append(node_number)
        self.number_of_nodes += 1
        self.reset_optimization_variables()
        return node_name, node_number

    def _remove_component(self, component: list) -> None:
        """Remove all nodes in component."""
        for node_number in sorted(component, reverse=True):
            self._remove_node_from_network(node_number=node_number)

    def _suppress_node(self, node_number: int):
        """Remove node and connect its parent to its child if they both exist."""
        self._suppress_component([node_number])

    def _suppress_component(self, component: list):
        """Remove all nodes in component and connect its parent to its child if they both exist."""
        in_nodes = self._in_nodes_of_nodes(component)
        out_nodes = self._out_nodes_of_nodes(component)
        self._remove_component(component)

        if len(in_nodes) == 1 and len(out_nodes) == 1:
            out_node = out_nodes[0]
            in_node = in_nodes[0]
            out_node = shifted_node_number(out_node, component)
            in_node = shifted_node_number(in_node, component)
            if not self._is_arc(in_node, out_node):
                self._add_arc(in_node, out_node)

    # --------------------------------------------------------- EVOLVE/TERMINATE METHODS --------------------------------------------------------- #
    def terminate_leaves(self, leaf_names_to_terminate: list = None, leaf_names_to_keep: list = None):
        assert bool(leaf_names_to_terminate is not None) or bool(leaf_names_to_keep is not None), "Exactly one of these parameters has to be set"
        leaf_numbers_to_terminate = self.get_node_numbers_of(coalesce(leaf_names_to_terminate, []))
        leaf_numbers_to_keep = self.get_node_numbers_of(coalesce(leaf_names_to_keep, []))
        self._terminate_leaves(leaf_numbers_to_terminate=leaf_numbers_to_terminate, leaf_numbers_to_keep=leaf_numbers_to_keep)

    def _terminate_leaves(self, leaf_numbers_to_terminate: list = None, leaf_numbers_to_keep: list = None):
        leaf_numbers_to_terminate = coalesce(leaf_numbers_to_terminate, [])
        leaf_numbers_to_keep = coalesce(leaf_numbers_to_keep, [])
        if len(leaf_numbers_to_terminate) == 0:
            leaf_numbers_to_terminate = sorted(list(set(self.leaf_numbers).difference(leaf_numbers_to_keep)), reverse=True)
        all_removed_node_numbers = []
        print(self.get_node_names_of(leaf_numbers_to_terminate))
        for current_leaf in leaf_numbers_to_terminate:
            print()
            self._terminate_node(current_leaf, all_removed_node_numbers)


    def _terminate_leaf(self, leaf_number: int):
        """Remove leaf and prune to make a proper network again"""
        self._terminate_node(leaf_number)

    def terminate_leaf(self, leaf_name: str):
        """Remove leaf and prune to make a proper network again"""
        assert self.is_leaf_node(leaf_name), "Node must be a leaf"
        self._terminate_leaf(self.get_node_number_of(leaf_name))

    def add_leaf_to_edge(self, edge: list, leaf_name: str = None, char_type='alph') -> (str, str):
        """Add node between parent and child and attach leaf with leaf_name to it."""
        return self.get_node_names_of(list(self._add_leaf_to_edge(self.get_node_numbers_of(edge), leaf_name, char_type)))

    def _add_leaf_to_edge(self, edge, leaf_name: str = None, char_type='alph') -> (int, int):
        parent = edge[0]
        child = edge[1]

        internal_name, internal_number = self._add_node_to_network()
        leaf_name, leaf_number = self._add_node_to_network(leaf_name, leaf=True, char_type=char_type)

        self._remove_arc(parent, child)

        self._add_arc(parent, internal_number)
        self._add_arc(internal_number, child)
        self._add_arc(internal_number, leaf_number)
        return internal_number, leaf_number

    # TODO: recombine multiple
    def recombine_leaves(self, leaf_name_1: str, leaf_name_2: str) -> None:
        assert self.is_leaf_node(leaf_name_1) and self.is_leaf_node(leaf_name_2), "Can not recombine non-leaf nodes"
        self._recombine_leaves(self.get_node_number_of(leaf_name_1), self.get_node_number_of(leaf_name_2))

    def _recombine_leaves(self, leaf_number_1: int, leaf_number_2: int) -> None:
        assert self._is_leaf_node(leaf_number_1) and self._is_leaf_node(leaf_number_2), "Can not recombine non-leaf nodes"
        self._remove_leaf_status_node(leaf_number_1)
        self._remove_leaf_status_node(leaf_number_2)

        leaf_name_1 = self.get_node_name_of(leaf_number_1)
        leaf_name_2 = self.get_node_name_of(leaf_number_2)

        self._rename_node(leaf_number_1, self.first_unused_internal_name())
        self._rename_node(leaf_number_2, self.first_unused_internal_name())

        new_leaf_name_1 = leaf_name_1
        new_leaf_name_2 = leaf_name_2
        recombined_leaf_name = "(" + leaf_name_1 + "+" + leaf_name_2 + ")"

        new_leaf_name_1, new_leaf_number_1 = self._add_node_to_network("(" + new_leaf_name_1 + "-R)", leaf=True)
        new_leaf_name_2, new_leaf_number_2 = self._add_node_to_network("(" + new_leaf_name_2 + "-R)", leaf=True)
        internal_leaf_name, internal_leaf_number = self._add_node_to_network(leaf=False)
        recombined_leaf_name, recombined_leaf_number = self._add_node_to_network(recombined_leaf_name, leaf=True)

        self._add_arc(leaf_number_1, new_leaf_number_1)
        self._add_arc(leaf_number_2, new_leaf_number_2)
        self._add_arc(leaf_number_1, internal_leaf_number)
        self._add_arc(leaf_number_2, internal_leaf_number)
        self._add_arc(internal_leaf_number, recombined_leaf_number)

    # TODO --> into n leaves
    def split_leaf(self, leaf_name_to_split: str):
        """ Split leaf into two leaves"""
        assert self.is_leaf_node(leaf_name_to_split), "Can not split non-leaf-node"
        self._split_leaf(self.get_node_number_of(leaf_name_to_split))

    def _split_leaf(self, leaf_number_to_split: int):
        self._remove_leaf_status_node(leaf_number_to_split)
        leaf_name_to_split = self.get_node_name_of(leaf_number_to_split)
        _, left_leaf_number = self._add_node_to_network(node_name='(' + leaf_name_to_split + "-I)", leaf=True)
        _, right_leaf_number = self._add_node_to_network(node_name='(' + leaf_name_to_split + "-II)", leaf=True)
        self._add_arc(leaf_number_to_split, left_leaf_number)
        self._add_arc(leaf_number_to_split, right_leaf_number)

    def replace_leaf_with_network(self, leaf_name_to_replace: str, replacement_network, replace_names: bool = False, char_type='alph'):
        assert self.is_leaf_node(leaf_name_to_replace), f"Cannot replace leaf {leaf_name_to_replace} with network as it is not a leaf."
        assert replace_names or set(self.leaf_names).isdisjoint(
            replacement_network.leaf_names), f"Cannot replace leaf {leaf_name_to_replace} with network \n \n {replacement_network}  \n \n as it has some leafs same as \n \n {self}"
        leaf_number_to_replace = self.get_node_number_of(leaf_name_to_replace)
        return self._replace_leaf_with_network(leaf_number_to_replace, replacement_network, replace_names, char_type)

    def _replace_leaf_with_network(self, leaf_number_to_replace: int, replacement_network, replace_names: bool = False, char_type='alph'):
        replacement_network = copy.deepcopy(replacement_network)
        replacement_network_internal_node_names = set(replacement_network.internal_node_names)
        replacement_network_root_name = replacement_network.root

        # Add internal nodes from replacement network to current network with primed names
        for node_name in replacement_network.internal_node_names:
            self._add_node_to_network(node_name + "*", leaf=False, char_type=char_type)

        # Add leaves from replacement network to current network
        if replace_names:
            for leaf_name in replacement_network.leaf_names:
                self._add_node_to_network(leaf_name + "*", leaf=True, char_type=char_type)
        else:
            for leaf_name in replacement_network.leaf_names:
                self._add_node_to_network(leaf_name, leaf=True, char_type=char_type)

        # Replace leaf_name with root of replacement network
        parent_of_leaf = self._in_nodes_of_node(leaf_number_to_replace)[0]
        self._remove_node_from_network(leaf_number_to_replace)
        self._add_arc(parent_of_leaf, self.get_node_number_of(replacement_network_root_name + "*"))

        # Add all connections from replacement network to current network
        for internal_node_name in replacement_network_internal_node_names:
            to_node_names = replacement_network.out_nodes_of_node(internal_node_name)
            internal_node_number = self.get_node_number_of(internal_node_name + "*")
            for to_node_name in to_node_names:
                if to_node_name in replacement_network_internal_node_names or replace_names:
                    to_node_name += "*"
                self._add_arc(internal_node_number, self.get_node_number_of(to_node_name))

        if replace_names:
            self.standardize_node_names(char_type)

    # ---------------------------------------------------------------- DEGREE METHODS ---------------------------------------------------------------- #
    def in_degree_node(self, node_name: str) -> int:
        """Retrieve number of arcs entering node_name."""
        return self._in_degree_node(self.get_node_number_of(node_name))

    def _in_degree_node(self, node_number: int) -> int:
        """Retrieve number of arcs entering node_name."""
        return sum(np.multiply(self.adj_matrix.T[node_number], self.adj_matrix.T[node_number] > 0))

    def out_degree_node(self, node_name: str) -> int:
        """Retrieve number of arcs exiting node_name."""
        return self._out_degree_node(self.get_node_number_of(node_name))

    def _out_degree_node(self, node_number: int) -> int:
        return sum(np.multiply(self.adj_matrix[node_number], self.adj_matrix[node_number] > 0))

    def get_in_degrees(self, leafless: bool = False) -> list:
        """Retrieve number of arcs entering each node."""
        mask = self.leaf_mask if leafless else np.ones(self.number_of_nodes)
        in_degree = np.sum(np.multiply(np.multiply(self.adj_matrix, self.adj_matrix > 0), mask), axis=0)
        try:
            return list(in_degree)
        except TypeError:
            return [in_degree]

    def get_out_degrees(self, leafless: bool = False) -> list:
        """Retrieve number of arcs exiting each node."""
        mask = self.leaf_mask if leafless else np.ones(self.number_of_nodes)
        out_degree = abs(np.sum(np.multiply(np.multiply(self.adj_matrix.T, self.adj_matrix.T > 0), mask), axis=0))
        try:
            return list(out_degree)
        except TypeError:
            return [out_degree]

    def in_out_degree_node(self, node_name: str):
        return self._in_out_degree_node(self.get_node_number_of(node_name))

    def _in_out_degree_node(self, node_number: int):
        return self._in_degree_node(node_number), self._out_degree_node(node_number)

    # ---------------------------------------------------------------- CONNECTION METHODS ---------------------------------------------------------------- #
    def in_nodes_of_node(self, node_name: str) -> list:
        """Retrieve all nodes which enter node node_name."""
        return self.get_node_names_of(self._in_nodes_of_node(self.get_node_number_of(node_name)))

    def _in_nodes_of_node(self, node_number: int) -> list:
        return list(np.where(self.adj_matrix[node_number, :] == -1)[0]) + list(np.where(self.adj_matrix[node_number, :] == -2)[0]) * 2

    def out_nodes_of_node(self, node_name: str) -> list:
        """Retrieve all nodes which exit node node_name."""
        return self.get_node_names_of(self._out_nodes_of_node(self.get_node_number_of(node_name)))

    def _out_nodes_of_node(self, node_number: int) -> list:
        return list(np.where(self.adj_matrix[node_number, :] == 1)[0]) + list(np.where(self.adj_matrix[node_number, :] == 2)[0]) * 2

    def in_nodes_of_nodes(self, component: list) -> list:
        """Retrieve all nodes which enter component."""
        return self.get_node_names_of(self._in_nodes_of_nodes(self.get_node_numbers_of(component)))

    def _in_nodes_of_nodes(self, component: list) -> list:
        in_nodes = np.array([])
        for node_number in component:
            in_nodes = np.array(list(itertools.chain(in_nodes, self._in_nodes_of_node(node_number))))
        return [node for node in in_nodes if node not in component]

    def out_node_of_nodes(self, component: list) -> list:
        """Retrieve all nodes which exit component."""
        return self.get_node_names_of(self._out_nodes_of_nodes(self.get_node_numbers_of(component)))

    def _out_nodes_of_nodes(self, component: list) -> list:
        in_nodes = np.array([])
        for node_number in component:
            in_nodes = np.array(list(itertools.chain(in_nodes, self._out_nodes_of_node(node_number))))
        return [node for node in in_nodes if node not in component]

    def connections_of_node(self, node_name: str) -> list:
        """Retrieve all nodes which enter and exit node."""
        return self.get_node_names_of(self._connections_of_node(self.get_node_number_of(node_name)))

    def _connections_of_node(self, node_number: int) -> list:
        return self._in_nodes_of_node(node_number) + self._out_nodes_of_node(node_number)

    def leaves_below_nodes(self, parent_names: set, max_depth: int = -1) -> set:
        """Retrieve all leaves max_depth below parents."""
        parent_numbers = set(self.get_node_numbers_of(list(parent_names)))
        leaf_numbers = list(self._leaves_below_nodes(parent_numbers, max_depth))
        return set(self.get_node_names_of(leaf_numbers))

    def _leaves_below_nodes(self, parent_names: set, max_depth: int = -1) -> set:
        """Retrieve all leaves max_depth below parents."""
        return set([child for child in self._nodes_below_nodes(parent_names, max_depth) if self._is_leaf_node(child)])

    def nodes_below_nodes(self, parent_names: set, max_depth: int = -1) -> set:
        """Retrieve all nodes max_depth below parents."""
        parent_numbers = set(self.get_node_numbers_of(list(parent_names)))
        children_numbers = self._nodes_below_nodes(parent_numbers, max_depth)
        return set(self.get_node_names_of(list(children_numbers)))

    def _nodes_below_nodes(self, parent_numbers: set, max_depth: int = -1) -> set:
        """Retrieve all nodes max_depth below parents."""
        if max_depth == 0:
            return parent_numbers
        next_generation = set()
        for parent_number in parent_numbers:
            next_generation.update(set(self._out_nodes_of_node(parent_number)))
        if next_generation:
            parent_numbers.update(self._nodes_below_nodes(next_generation, max_depth - 1))
            return parent_numbers
        return parent_numbers

    def nodes_above_nodes(self, children_names: set, max_height: int = -1) -> set:
        """Retrieve all nodes max_depth above children."""
        children_numbers = set(self.get_node_numbers_of(list(children_names)))
        parent_numbers = self._nodes_above_nodes(children_numbers, max_height)
        return set(self.get_node_names_of(list(parent_numbers)))

    def _nodes_above_nodes(self, children_numbers: set, max_height: int = -1) -> set:
        if max_height == 0:
            return children_numbers
        prev_generation = set()
        for child_number in children_numbers:
            prev_generation.update(set(self._in_nodes_of_node(child_number)))
        if prev_generation:
            children_numbers.update(self._nodes_above_nodes(prev_generation, max_height - 1))
            return children_numbers
        return children_numbers

    def connected_components(self) -> list:
        connected_components = self._connected_components()
        return [self.get_node_names_of(connected_component) for connected_component in connected_components]

    def _connected_components(self) -> list:
        unchecked_nodes = set(self.node_name_map.values())
        components = []
        while unchecked_nodes:
            node_number = unchecked_nodes.pop()
            component = [node_number]
            unchecked_connections = set(self._connections_of_node(node_number)).intersection(unchecked_nodes)
            while unchecked_connections:
                connection_number = unchecked_connections.pop()
                unchecked_nodes.remove(connection_number)
                new_connections = set(self._connections_of_node(connection_number))
                unchecked_connections = (unchecked_connections.union(new_connections)).intersection(unchecked_nodes)
                component.append(connection_number)
            if not (len(component) == 1 and self._is_leaf_node(component[0])):
                components.append(component)
        return components

    # ---------------------------------------------------------------- OPTIMIZATION METHODS ---------------------------------------------------------------- #
    def optimization_available(self) -> bool:
        return self._biconnected_components is not None and self._o_partial_ordering is not None and self._o_cut_arc_matrix is not None and self._o_cut_arc_sets is not None

    def reset_optimization_variables(self) -> None:
        self._o_cut_arc_matrix = None
        self._o_cut_arc_sets = None
        self._o_biconnected_components = None
        self._o_partial_ordering = None
        self._o_leaf_ordering = None

    def calculate_optimization_variables(self) -> None:
        self._compute_cut_arc_matrix()
        self._compute_cut_arc_sets()
        self._compute_biconnected_components()
        self._compute_partial_ordering()
        self._compute_leaf_ordering()

    def _compute_cut_arc_matrix(self) -> None:
        visited = [False] * self.number_of_nodes
        disc = [-1] * self.number_of_nodes
        low = [-1] * self.number_of_nodes
        parent = [None] * self.number_of_nodes
        self._o_cut_arc_matrix = np.zeros((self.number_of_nodes, self.number_of_nodes))

        for i in range(self.number_of_nodes):
            if not visited[i]:
                self._cut_arc_helper(i, visited, disc, low, parent, self._o_cut_arc_matrix)

    def _cut_arc_helper(self, u_number, visited, disc, low, parent, cut_arc_matrix, t=0):
        # TODO source: https://www.geeksforgeeks.org/bridge-in-a-graph/
        visited[u_number] = True
        disc[u_number] = t
        low[u_number] = t
        t += 1

        for v_number in self._connections_of_node(u_number):
            if not visited[v_number]:
                parent[v_number] = u_number
                self._cut_arc_helper(v_number, visited, disc, low, parent, cut_arc_matrix, t)

                low[u_number] = min(low[u_number], low[v_number])

                if low[v_number] > disc[u_number]:
                    cut_arc_matrix[u_number][v_number] = 1
            elif v_number != parent[u_number]:
                low[u_number] = min(low[u_number], disc[v_number])

    def _compute_cut_arc_sets(self) -> None:
        """Retrieve cut-arc sets of network."""
        cut_arc_matrix = self.cut_arc_matrix
        _, to_nodes = np.where(cut_arc_matrix == 1)
        self._o_cut_arc_sets = []
        for to_node in to_nodes:
            ca_set = list(self._nodes_below_nodes({to_node}))
            self._o_cut_arc_sets.append(ca_set)

    def _compute_biconnected_components(self) -> None:
        cut_arc_matrix = self.cut_arc_matrix
        self.adj_matrix -= cut_arc_matrix - cut_arc_matrix.T
        self._o_biconnected_components = self._connected_components()
        self.adj_matrix += cut_arc_matrix - cut_arc_matrix.T

    def _compute_partial_ordering(self) -> None:
        root = self._root
        if self._is_leaf_node(root):
            self._o_partial_ordering = [[]]
        new_result = [[root]]
        self._o_partial_ordering = []
        changes = True
        while changes:
            changes = False
            self._o_partial_ordering = new_result
            new_result = []
            for track in self._o_partial_ordering:
                new_tracks = []
                children = self._out_nodes_of_node(track[-1])
                for child in children:
                    if not self._is_leaf_node(child):
                        changes = True
                        new_tracks.append(track + [child])
                for new_track in new_tracks:
                    new_result.append(new_track)
                if len(new_tracks) == 0:
                    new_result.append(track)

    def _compute_leaf_ordering(self) -> None:
        partial_node_ordering = self._partial_ordering
        self._o_leaf_ordering = []
        for track in partial_node_ordering:
            leaf_track = []
            for node in track:
                children = self._out_nodes_of_node(node)
                if len(children) == 1:
                    if self._is_leaf_node(children[0]):
                        leaf_track.append(children[0])
                        self._o_leaf_ordering.append(leaf_track)
                else:
                    count = 0
                    leaf_track_1, leaf_track_2 = [], []
                    if self._is_leaf_node(children[0]):
                        leaf_track_1 = copy.copy(leaf_track)
                        leaf_track_1.append(children[0])
                        count += 1
                    if self._is_leaf_node(children[1]):
                        leaf_track_2 = copy.copy(leaf_track)
                        leaf_track_2.append(children[1])
                        count += 2
                    if count == 1:
                        leaf_track = leaf_track_1
                    elif count == 2:
                        leaf_track = leaf_track_2
                    elif count == 3:
                        self._o_leaf_ordering.append(leaf_track_1)
                        self._o_leaf_ordering.append(leaf_track_2)

    # @property
    # def cut_arc_matrix_old(self) -> np.ndarray:
    #     if self.__cut_arc_matrix is None:
    #         self.__cut_arc_matrix = np.zeros((self.number_of_nodes, self.number_of_nodes))
    #         edges = self._edges()
    #         for edge in edges:
    #             from_node, to_node = edge
    #             self.adj_matrix[from_node][to_node] -= 1
    #             self.adj_matrix[to_node][from_node] += 1
    #             components = self.connected_components()
    #             if len(components) > 1:
    #                 self.__cut_arc_matrix[from_node][to_node] = 1
    #                 self.__cut_arc_matrix[to_node][from_node] = 1
    #             self.adj_matrix[from_node][to_node] += 1
    #             self.adj_matrix[to_node][from_node] -= 1
    #     return self.__cut_arc_matrix

    # ---------------------------------------------------------------- NODE TYPE METHODS ---------------------------------------------------------------- #
    def is_leaf_node(self, node_name: str) -> bool:
        return self._is_leaf_node(self.get_node_number_of(node_name))

    def _is_leaf_node(self, node_number: int) -> bool:
        return node_number in self.leaf_numbers

    def is_reticulation_node(self, node_name: str) -> bool:
        return self._is_reticulation_node(self.get_node_number_of(node_name))

    def _is_reticulation_node(self, node_number: int) -> bool:
        return self._in_degree_node(node_number) > 1

    def is_root_node(self, node_name: str) -> bool:
        return self._is_root_node(self.get_node_number_of(node_name))

    def _is_root_node(self, node_number: int) -> bool:
        return self._in_degree_node(node_number) == 0

    def is_tree_node(self, node_name: str) -> bool:
        return self._is_tree_node(self.get_node_number_of(node_name))

    def _is_tree_node(self, node_number: int) -> bool:
        return not self._is_root_node(node_number) and self._out_degree_node(node_number) > 1

    def get_node_type(self, node_name: str) -> str:
        return self._get_node_type(self.get_node_number_of(node_name))

    def _get_node_type(self, node_number: int) -> str:
        """Find out what type of node node_name is."""
        if self._is_leaf_node(node_number):
            return "leaf"
        if self._is_reticulation_node(node_number):
            return "reticulation"
        if self._is_root_node(node_number):
            return "root"
        if self._is_tree_node(node_number):
            return "tree"

    # ---------------------------------------------------------------- COMPONENT METHODS ---------------------------------------------------------------- #
    def _is_redundant_component(self, component: list) -> bool:
        """Check if component is redundant."""
        self.logger.debug("Checking if component {} is redundant.".format(component))
        out_nodes = self._out_nodes_of_nodes(component)
        if len(out_nodes) != 1:
            return False
        return True

    def _is_strongly_redundant_component(self, component: list) -> bool:
        """Check if component is strongly redundant."""
        self.logger.debug("Checking if component {} is strongly redundant.".format(component))
        out_nodes = self.out_node_of_nodes(component)
        if len(out_nodes) != 1:
            return False
        temp_network = copy.deepcopy(self)
        temp_network._remove_component(component)
        connected_comp = temp_network.connected_components()
        if len(connected_comp) > 1:
            return False
        return True

    def _get_component_list_names(self, component_list: list) -> list:
        """Retrieve names of list of lists of node numbers."""
        self.logger.debug("Retrieve node names of {}.".format(component_list))
        return [self.get_node_name_of(comp) for comp in component_list]

    def level_of_component(self, component):
        return self._level_of_component(self.get_node_numbers_of(component))

    def _level_of_component(self, component):
        return len(self._reticulations_of_component(component))

    def reticulations_of_component(self, component):
        return self._reticulations_of_component(self.get_node_numbers_of(component))

    def _reticulations_of_component(self, component):
        return [self._is_reticulation_node(node_number) for node_number in component]

    # ---------------------------------------------------------------- EVOLUTION/TERMINATE METHODS ---------------------------------------------------------- #
    def _terminate_node(self, node_number: int, removed_node_numbers: list = None) -> list:
        removed_node_numbers = coalesce(removed_node_numbers, [])
        node_number = shifted_node_number(node_number, removed_node_numbers)
        print(f"Terminating node {self.get_node_name_of(node_number)}")

        children_numbers = self._nodes_below_nodes({node_number}, max_depth=1).difference([node_number])
        parent_numbers = self._nodes_above_nodes({node_number}, max_height=1).difference([node_number])
        if len(children_numbers) == 0:
            print("Node is a dead end")
            self._remove_node_from_network(node_number)
            removed_node_numbers.append(node_number)
            for p_node_number in parent_numbers:
                self._terminate_node(p_node_number, removed_node_numbers)

        elif len(parent_numbers) == 0 and len(children_numbers) == 1:
            print("Node is a dead root")
            self._remove_node_from_network(node_number)
            removed_node_numbers.append(node_number)
        elif len(parent_numbers) == 1 and len(children_numbers) == 1:
            print("Node is suppressable")
            self._suppress_node(node_number)
            removed_node_numbers.append(node_number)

        return original_node_numbers(removed_node_numbers)

    def prune(self, suppress_redundant: str = 'all', suppress_parallel: bool = True):
        """Suppress al unnecessary/redundant/parallel nodes/edges in network."""
        changes = True
        while changes:
            changes = False

            # Check for parallel arcs
            if suppress_parallel:
                parallel_arcs = np.where(self.adj_matrix >= 2)
                parallel_arcs = zip(parallel_arcs[0], parallel_arcs[1])
                for from_node, to_node in parallel_arcs:
                    changes = True
                    self._remove_arc(from_node, to_node)

            # Check for redundant components:
            if suppress_redundant != 'none':
                bcs_to_remove = []
                bcs = self._biconnected_components
                for bc in bcs:
                    if suppress_redundant == 'strongly' and self._is_strongly_redundant_component(bc):
                        bcs_to_remove.append(bc)
                    if suppress_redundant == 'all' and self._is_redundant_component(bc):
                        bcs_to_remove.append(bc)
                removed_nodes = []
                for bc in bcs_to_remove:
                    changes = True
                    bc = [shifted_node_number(node, removed_nodes) for node in bc]
                    self._suppress_component(bc)
                    removed_nodes.extend(bc)

    def _stable_ancestors_of_nodes(self, node_numbers: list):
        leafless_node_numbers = []
        for node_number in node_numbers:
            if self._is_leaf_node(node_number):
                parent_of_leaf = self._nodes_above_nodes({node_number}, max_height=1)
                parent_of_leaf.remove(node_number)
                node_number = parent_of_leaf.pop()
            leafless_node_numbers.append(node_number)
        partial_ordering = self._partial_ordering
        tracks_to_min_node = []
        tracks_to_max_node = []
        for track in partial_ordering:
            mi = 10 ** 10
            ma = -1
            for node in leafless_node_numbers:
                try:
                    mi = min(mi, track.index(node))
                    ma = max(ma, track.index(node))
                except ValueError:
                    pass
            if mi != 10 ** 10:
                tracks_to_min_node.append(track[:mi + 1])
                tracks_to_max_node.append(track[:ma + 1])

        stable_ancestors = (set(tracks_to_min_node[0]).intersection(*tracks_to_min_node))
        stable_ancestors_indicis_dict = {tracks_to_min_node[0].index(sa): sa for sa in stable_ancestors}
        lowest_stable_ancestor_index = max(stable_ancestors_indicis_dict.keys())
        lowest_stable_ancestor = stable_ancestors_indicis_dict[lowest_stable_ancestor_index]

        tracks_between_nodes = []
        for track in tracks_to_max_node:
            lowest_stable_ancestor_index = track.index(lowest_stable_ancestor)
            tracks_between_nodes.append(track[lowest_stable_ancestor_index:])

        nodes_between = set(tracks_between_nodes[0]).union(*tracks_between_nodes)
        nodes_between.update(leafless_node_numbers)
        nodes_between.update(node_numbers)
        return stable_ancestors_indicis_dict, lowest_stable_ancestor, nodes_between

    def stable_ancestors_of_nodes(self, node_names: list):
        node_numbers = self.get_node_numbers_of(node_names)
        stable_ancestors_indicis_dict, lowest_stable_ancestor, nodes_between = self._stable_ancestors_of_nodes(node_numbers)
        # TODO: To names

    def visualize(self, file_path: str = None, internal_node_labels=True, edge_labels=False, rankdir='TB', format='png'):
        """Visualize network."""
        dot = Digraph()
        dot.graph_attr["rankdir"] = rankdir
        dot.engine = 'dot'
        edges = self.edges
        for node_name in self.node_name_map:
            if internal_node_labels or self.is_leaf_node(node_name):
                dot.node(name=node_name, label=node_name, **{'width': str(0.5), 'height': str(0.5)})
            else:
                dot.node(name=node_name, label="", **{'width': str(0.1), 'height': str(0.1), 'fillcolor': 'black', 'style': 'filled'})
        for index, edge in enumerate(edges):
            if edge_labels:
                dot.edge(edge[0], edge[1], label=str(index))
            else:
                dot.edge(edge[0], edge[1])
        if file_path:
            dot.format = format
            if file_path[-4:] == f'.{format}':
                dot.render(filename=file_path[:-4])
            else:
                dot.render(filename=file_path)
        else:
            dot.render(view=True)
            time.sleep(0.2)

    def to_df(self, directed: bool = True) -> pd.DataFrame:
        """Retrieve dataframe representation of network."""
        self.logger.debug("Retrieve {}directed dataframe representation of network.".format("" if directed else "un-"))
        ordered_node_names = self.node_names
        mask = self.adj_matrix > 0
        data = self.adj_matrix * mask if directed else self.adj_matrix
        length = self.number_of_nodes - self.number_of_leaves if directed else self.number_of_nodes
        result = pd.DataFrame(columns=ordered_node_names, index=ordered_node_names[:length], data=data[:length])
        return result.astype(int)

    # def network_set_consistency(self, other, n=3, max_processes=1, progress_bar=False):
    #     self_network_info_list = self.network_set(n=n, max_processes=max_processes, progress_bar=progress_bar)
    #     other_network_info_list = other.network_set(n=n, max_processes=max_processes, progress_bar=progress_bar)
    #
    #     if progress_bar:
    #         count = 0
    #         for self_network_info in tqdm(self_network_info_list, desc=f"Comparing network sets"):
    #             for other_network_info in other_network_info_list:
    #                 equal, _, _ = self_network_info.network.equal_structure(other_network_info.network, optimize=True, equal_naming=True)
    #                 if equal:
    #                     count += 1
    #                     break
    #     else:
    #         count = 0
    #         for self_network_info in self_network_info_list:
    #             for other_network_info in other_network_info_list:
    #                 equal, _, _ = self_network_info.network.equal_structure(other_network_info.network, optimize=True, equal_naming=True)
    #                 if equal:
    #                     count += 1
    #                     break
    #
    #     return count, len(self_network_info_list)
    #
    # def full_network_set_consistency(self, other, max_processes=1, progress_bar=False):
    #     total = 0
    #     count = 0
    #     for n in range(2, self.number_of_leaves + 1):
    #         extra_count, extra_total = self.network_set_consistency(other, n, max_processes, progress_bar)
    #         total += extra_total
    #         count += extra_count
    #     return count, total
    #
    # def tree_set_consistency(self, other, n=3, max_processes=1, progress_bar=False):
    #     self_tree_info_list = self.tree_set(n=n, max_processes=max_processes, progress_bar=progress_bar)
    #     other_tree_info_list = other.tree_set(n=n, max_processes=max_processes, progress_bar=progress_bar)
    #
    #     if progress_bar:
    #         count = 0
    #         for self_tree_info in tqdm(self_tree_info_list, desc=f"Comparing network sets"):
    #             for other_tree_info in other_tree_info_list:
    #                 equal, _, _ = self_tree_info.network.equal_structure(other_tree_info.network, optimize=True, equal_naming=True)
    #                 if equal:
    #                     count += 1
    #                     break
    #     else:
    #         count = 0
    #         for self_tree_info in self_tree_info_list:
    #             for other_tree_info in other_tree_info_list:
    #                 equal, _, _ = self_tree_info.network.equal_structure(other_tree_info.network, optimize=True, equal_naming=True)
    #                 if equal:
    #                     count += 1
    #                     break
    #
    #     return count, len(self_tree_info_list)
    #
    # def full_tree_set_consistency(self, other, max_processes=1, progress_bar=False):
    #     total = 0
    #     count = 0
    #     for n in range(2, self.number_of_leaves + 1):
    #         extra_count, extra_total = self.tree_set_consistency(other, n, max_processes, progress_bar)
    #         total += extra_total
    #         count += extra_count
    #     return count, total

    def equal_structure(self, other, optimize=True, equal_naming=False, progress_bar=False) -> (bool, list, list):
        if self.number_of_nodes != other.number_of_nodes:
            return False, [], []
        if self.number_of_leaves != other.number_of_leaves:
            return False, [], []
        if equal_naming and (set(self.leaf_names) != set(other.leaf_names)):
            return False, [], []
        if optimize or (self.optimization_available() and other.optimization_available()):
            self_cut_arc_sets = self._cut_arc_sets
            other_cut_arc_sets = other._cut_arc_sets
            if len(self_cut_arc_sets) != len(other_cut_arc_sets):
                return False, [], []
            if collections.Counter([len(cut_arc_set) for cut_arc_set in self_cut_arc_sets]) != \
                    collections.Counter([len(cut_arc_set) for cut_arc_set in other_cut_arc_sets]):
                return False, [], []

        if progress_bar:
            pbar = tqdm(total=len(self.node_name_map), desc=f"Checking if network {self.uid} is equal to network {other.uid}")
        else:
            pbar = None
        self_root = self._root
        other_root = other._root
        translation_dict = bidict()
        translation_dict[self_root] = other_root

        translation_dicts = self._equal_structure(self_current_node=self_root, other=other, other_current_node=other_root, translation_dicts=[translation_dict],
                                                  pbar=pbar, optimize=optimize, equal_naming=equal_naming)
        if progress_bar:
            pbar.close()
        translation_dicts = [bidict({self.get_node_name_of(key): other.get_node_name_of(item) for key, item in translation_dict.items()})
                             for translation_dict in translation_dicts]
        leaf_translation_dicts = [bidict({key: item for key, item in translation_dict.items() if self.is_leaf_node(key)}) for translation_dict in
                                  translation_dicts]
        return len(translation_dicts) > 0, translation_dicts, leaf_translation_dicts

    def _equal_structure(self, self_current_node, other, other_current_node, translation_dicts, pbar, optimize, equal_naming):
        if len(translation_dicts) == 0:
            return translation_dicts  # False

        if pbar:
            pbar.n = max([len(translation_dict) for translation_dict in translation_dicts])
            pbar.refresh()
        if self._is_leaf_node(self_current_node) and other._is_leaf_node(other_current_node):
            if equal_naming and self.get_node_name_of(self_current_node) != other.get_node_name_of(other_current_node):
                return []  # False
            for translation_dict in translation_dicts:
                if len(translation_dict) == self.number_of_nodes:
                    return [translation_dict]
            return translation_dicts  # True

        self_node_children = self._nodes_below_nodes({self_current_node}, max_depth=1).difference({self_current_node})
        other_node_children = other._nodes_below_nodes({other_current_node}, max_depth=1).difference({other_current_node})
        if len(other_node_children) != len(self_node_children):
            return []  # False
        if optimize or (self.optimization_available() and other.optimization_available()):
            self_cut_arc_sets = self._cut_arc_sets
            other_cut_arc_sets = other._cut_arc_sets
            self_node_cut_arc_sets = [len(cut_arc_set) for cut_arc_set in self_cut_arc_sets if self_current_node in cut_arc_set]
            other_node_cut_arc_sets = [len(cut_arc_set) for cut_arc_set in other_cut_arc_sets if other_current_node in cut_arc_set]
            if collections.Counter(self_node_cut_arc_sets) != collections.Counter(other_node_cut_arc_sets):
                return []  # False

        bijections_iterator = itertools.permutations(range(len(self_node_children)))
        self_node_children = list(self_node_children)
        other_node_children = list(other_node_children)
        translation_dicts_to_return = []
        for bijection in bijections_iterator:
            extra_translation_dict = bidict({self_node_children[index]: other_node_children[value] for index, value in enumerate(bijection)})
            current_translation_dicts = []
            for translation_dict in translation_dicts:
                if self.bijection_contradicts_bijection(translation_dict, extra_translation_dict):
                    continue
                current_translation_dicts.append(bidict({**translation_dict, **extra_translation_dict}))
                for new_self_node, new_other_node in extra_translation_dict.items():
                    next_translation_dicts = []
                    for current_translation_dict in current_translation_dicts:
                        new_dicts = self._equal_structure(new_self_node, other, new_other_node, [current_translation_dict], pbar, optimize=optimize,
                                                          equal_naming=equal_naming)
                        for new_dict in new_dicts:
                            if len(new_dict) == self.number_of_nodes:
                                return [new_dict]
                        next_translation_dicts.extend(new_dicts)
                    current_translation_dicts = next_translation_dicts
            translation_dicts_to_return.extend(current_translation_dicts)
        return translation_dicts_to_return

    @staticmethod
    def bijection_contradicts_bijection(bijection_1, bijection_2):
        if len(bijection_1) > len(bijection_2):
            bijection_1, bijection_2 = bijection_2, bijection_1
        for key, value in bijection_1.items():
            if key in bijection_2:
                if bijection_2[key] != value:
                    return True
            if value in bijection_2.inverse:
                if bijection_2.inverse[value] != key:
                    return True
        return False

    def evolve_times(self, times: int, recombination_chance: int = 0, progress_bar=False) -> str:
        result = ""
        pbar = tqdm(total=times, desc=f"Evolving network")
        while times > 0:
            if times == 1:
                recombined, res = self.evolve(0)
            else:
                recombined, res = self.evolve(recombination_chance)
            times -= 1 + int(recombined)
            result += res
            pbar.update(1 + int(recombined))
        pbar.close()

        changes = True
        while changes:
            changes = False
            in_degrees = np.array(self.get_in_degrees())
            biconnected_components = self._biconnected_components
            for bc in biconnected_components:
                in_sum_bc = in_degrees[bc]
                bc_reticulations = {bc[node_number] for node_number in np.where(in_sum_bc == 2)[0]}
                if len(bc_reticulations) - self.level > 0:
                    changes = True
                    reticulation_to_remove = bc_reticulations.pop()
                    parent_reticulation = self._nodes_above_nodes({reticulation_to_remove}, max_height=1).difference([reticulation_to_remove]).pop()
                    child_reticulation = self._nodes_below_nodes({reticulation_to_remove}, max_depth=1).difference([reticulation_to_remove]).pop()
                    self._remove_arc(parent_reticulation, reticulation_to_remove)
                    self._remove_arc(reticulation_to_remove, child_reticulation)
                    self._add_arc(parent_reticulation, child_reticulation)
                    self._add_leaf_status_node(reticulation_to_remove)

        self.prune()
        return result

    def evolve(self, recombination_chance: int = 0) -> (bool, str):
        # Note: does not keep level intact perfectly
        assert 0 <= recombination_chance <= 1
        recombined = False
        if np.random.uniform() < recombination_chance:
            try:
                node_levels_dict = {0: set(self.leaf_numbers)}
                in_degrees = np.array(self.get_in_degrees())
                biconnected_components = self._biconnected_components
                for bc in biconnected_components:
                    in_sum_bc = in_degrees[bc]
                    bc_level = sum(in_sum_bc > 1)
                    bc = self._leaves_below_nodes(set(bc), max_depth=1)
                    try:
                        node_levels_dict[bc_level].update(bc)
                    except KeyError:
                        node_levels_dict[bc_level] = set(bc)
                    if bc_level > 0:
                        node_levels_dict[0].difference_update(bc)
                first_leaf = set().union(*node_levels_dict.values()).pop()
                first_leaf_name = self.get_node_name_of(first_leaf)
                max_level = max(node_levels_dict.keys())
                for level in set(range(max(max_level, self.level))).difference(node_levels_dict.keys()):
                    node_levels_dict[level] = {}
                first_leaf_level = np.where(np.array([int(first_leaf in node_levels_dict[level]) for level in range(max_level + 1)]) == 1)[0][0]
                max_second_leaf_level = self.level - first_leaf_level - 1
                second_leaf = set(itertools.chain.from_iterable([node_levels_dict[level] for level in range(max_second_leaf_level + 1)])).difference(
                    [first_leaf]).pop()
                second_leaf_name = self.get_node_name_of(second_leaf)
                self._recombine_leaves(first_leaf, second_leaf)
                result = f"[{first_leaf_name} , {second_leaf_name}]-"
            except KeyError:
                result = ""
                return recombined, result
            recombined = True
        else:
            leaf_number = random.choice(self.leaf_numbers)
            self._split_leaf(leaf_number)
            result = f"[{self.get_node_name_of(leaf_number)}]-"

        self.reset_optimization_variables()

        return recombined, result

    def terminate(self, termination_probability: int = 0.1):
        to_terminate = []
        for leaf_number in self.leaf_numbers:
            if np.random.uniform(low=0.0, high=1.0, size=1)[0] <= termination_probability:
                to_terminate.append(leaf_number)
        for leaf_number in sorted(to_terminate, reverse=True):
            self._terminate_node(leaf_number)
        self.prune()
        self.standardize_node_names()

    def enewick(self):
        return self._enewick(self._root)

    def _enewick(self, current_node, traversed_nodes=None):
        traversed_nodes = coalesce(traversed_nodes, [])
        if current_node in traversed_nodes:
            return self.get_node_name_of(current_node)
        traversed_nodes.append(current_node)
        current_node_children = self._nodes_below_nodes({current_node}, max_depth=1).difference([current_node])
        if len(current_node_children) == 0:
            return self.get_node_name_of(current_node)
        elif len(current_node_children) == 1:
            child = current_node_children.pop()
            return "(" + self._enewick(child, traversed_nodes) + ")" + self.get_node_name_of(current_node)
        elif len(current_node_children) == 2:
            child_1 = current_node_children.pop()
            child_2 = current_node_children.pop()
            return "(" + self._enewick(child_1, traversed_nodes) + ", " + self._enewick(child_2, traversed_nodes) + ")" + self.get_node_name_of(current_node)

    def __str__(self):
        return str(self.to_df())

    def __getstate__(self):
        self.logger = 'network.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        self.logger = logging.getLogger(self.logger)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger)
        return self.__dict__


""" ############################################################################################################################################################
                                                                      NETWORK INFO
############################################################################################################################################################ """


class NetworkInfo:
    def __init__(self, network: RootedLevelKNetwork, info: dict = None):
        self.network = network
        self.info = coalesce(info, dict())
        self.multiplicity = 1

    def collapse(self, leaf_set) -> list:
        new_networks = self.network.collapse(leaf_set)
        return new_networks

    def calculate_info(self):
        self['cut_arc_sets'] = self.network.cut_arc_sets
        self['strict_level'] = self.network.strict_level()
        self['leaf_order'] = self.network.leaf_ordering

    def add_info(self, network_info):
        for key, value in network_info.info.items():
            self.info[key] = value

    @classmethod
    def restrict(cls, network_info, leaf_names):
        network = RootedLevelKNetwork.restrict(network_info.network, leaf_names=leaf_names, suppress_parallel=True, suppress_redundant='all')
        result = cls(network)
        return result

    def __getitem__(self, item):
        return self.info[item]

    def __setitem__(self, item, value):
        self.info[item] = value

    def __copy__(self):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.network = copy.copy(self.network)
        cp.info = copy.copy(self.info)
        cp.multiplicity = copy.copy(self.multiplicity)
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.network = copy.deepcopy(self.network)
        cp.info = copy.deepcopy(self.info)
        cp.multiplicity = copy.copy(self.multiplicity)
        return cp

    def __str__(self):
        return str(self.network.leaf_names) + "\n" + pp.pformat(self.info)


""" ############################################################################################################################################################
                                                                  NETWORK INFO LIST
############################################################################################################################################################ """


class NetworkInfoList:
    def __init__(self, network_size: int = None, lst: list = None):
        # Set up logging
        self.uid = guid()
        self.logger = logging.getLogger('network_info_list.{}'.format(self.uid))
        self.logger.debug(f"Created NetworkInfoList ({self.uid})")

        self.list = coalesce(lst, [])
        self.network_size = network_size

    @classmethod
    def collapse(cls, network_info_list, leaf_set: set):
        result = cls(network_size=network_info_list.network_size)
        for network_info in network_info_list:
            for new_network in network_info.collapse(leaf_set):
                new_network_info = NetworkInfo(new_network)
                new_network_info.multiplicity = network_info.multiplicity
                result.append(new_network_info)
        return result

    @classmethod
    def induced_networks_of_network_info_list(cls, network_info_list, network_size: int, max_processes: int = 1, progress_bar: bool = False):
        """ Network list containing all exhibited networks of networks in list of certain size """
        assert network_size < network_info_list.network_size, "Can only create smaller network info lists with smaller network size."
        result = cls(network_size)
        for network_info in network_info_list:
            result.extend(cls.induced_strict_network_set(network_info.network, network_size, max_processes, progress_bar=False))
        return result

    @classmethod
    def induced_strict_network_set(cls, network: RootedLevelKNetwork, network_size: int, max_processes: int = 1, progress_bar: bool = False,
                                   method='Recursive'):
        """ Network list containing all exhibited network of network of certain size """
        if method == 'Recursive':
            return cls._induced_strict_network_set_recursive(network, network_size, max_processes, progress_bar, network.leaf_names, 0)
        elif method == 'Iterative':
            return cls.induced_strict_network_set_iterative(network, network_size, max_processes, progress_bar)

    @classmethod
    def _induced_strict_network_set_recursive(cls, network: RootedLevelKNetwork, network_size: int, max_processes: int, progress_bar: bool, leaf_names: list,
                                              index: int):
        result = cls(network_size=network_size)
        number_of_leaves = network.number_of_leaves
        if number_of_leaves == network_size:
            network.prune()
            result.append(NetworkInfo(network))
            return result
        else:
            terminate_leaf_name_iterator = [leaf_names[index] for index in range(index, len(leaf_names) - (number_of_leaves - network_size) + 1)]
            if max_processes == 1:
                if progress_bar:
                    for index_2, terminate_leaf_name in tqdm(enumerate(terminate_leaf_name_iterator), total=len(terminate_leaf_name_iterator),
                                                             desc="Computing exhibited trinets recursively"):
                        new_network = copy.deepcopy(network)
                        new_network.terminate_leaf(terminate_leaf_name)
                        result.extend(
                            cls._induced_strict_network_set_recursive(new_network, network_size, max_processes=1, progress_bar=False, leaf_names=leaf_names,
                                                                      index=index + index_2 + 1))
                else:
                    for index_2, terminate_leaf_name in enumerate(terminate_leaf_name_iterator):
                        new_network = copy.deepcopy(network)
                        new_network.terminate_leaf(terminate_leaf_name)
                        result.extend(
                            cls._induced_strict_network_set_recursive(new_network, network_size, max_processes=1, progress_bar=False, leaf_names=leaf_names,
                                                                      index=index + index_2 + 1))
            else:
                pool = multiprocessing.Pool(max_processes)
                index_2 = 0
                if progress_bar:
                    for new_network in tqdm(
                            pool.imap_unordered(cls._induced_strict_network_set_recursive_helper, zip(terminate_leaf_name_iterator, itertools.repeat(network))),
                            total=len(terminate_leaf_name_iterator), desc=f"Computing exhibited trinets recursively using {max_processes} processes"):
                        result.extend(cls._induced_strict_network_set_recursive(
                            new_network
                            , network_size
                            , max_processes=1
                            , progress_bar=False
                            , leaf_names=leaf_names
                            , index=index + index_2 + 1))
                        index_2 += 1
                else:
                    for new_network in pool.imap_unordered(
                            cls._induced_strict_network_set_recursive_helper, zip(terminate_leaf_name_iterator, itertools.repeat(network))):
                        result.extend(cls._induced_strict_network_set_recursive(
                            new_network
                            , network_size
                            , max_processes=1
                            , progress_bar=False
                            , leaf_names=leaf_names
                            , index=index + index_2 + 1))
                        index_2 += 1
                pool.close()
                pool.join()
            return result

    @classmethod
    def _induced_strict_network_set_recursive_helper(cls, leaf_network):
        terminate_leaf_name, network = leaf_network
        new_network = copy.deepcopy(network)
        new_network.terminate_leaf(terminate_leaf_name)
        return new_network

    @classmethod
    def induced_strict_network_set_iterative(cls, network: RootedLevelKNetwork, network_size: int, max_processes: int = 1, progress_bar: bool = False):
        result = cls(network_size=network_size)
        n_set_iterator = itertools.combinations(network.leaf_numbers, 3)
        pool = multiprocessing.Pool(max_processes)
        if progress_bar:
            for new_network in tqdm(
                    pool.imap_unordered(cls._induced_strict_network_set_iterative_helper, zip(n_set_iterator, itertools.repeat(network))),
                    total=ncr(network.number_of_leaves, network_size), desc=f"Computing exhibited trinets iteratively using {max_processes} processes"):
                result.append(NetworkInfo(new_network))
        else:
            for new_network in pool.imap_unordered(cls._induced_strict_network_set_iterative_helper, zip(n_set_iterator, itertools.repeat(network))):
                result.append(NetworkInfo(new_network))
        pool.close()
        pool.join()
        return result

    @classmethod
    def _induced_strict_network_set_iterative_helper(cls, n_set_network):
        n_set, network = n_set_network
        return RootedLevelKNetwork._restrict(network, n_set)

    @classmethod
    def displayed_trees(cls, network: RootedLevelKNetwork, max_processes: int = 1, progress_bar: bool = False):
        network_info_list = cls(network_size=network.number_of_leaves)
        network_info_list.append(NetworkInfo(network))
        return cls._displayed_trees_recursive(network_info_list, max_processes, progress_bar, network.reticulations, 0)

    @classmethod
    def _displayed_trees_recursive(cls, network_info_list, max_processes: int, progress_bar: bool, reticulation_names: list, index: int):
        result = cls(network_size=network_info_list.network_size)
        number_of_reticulations = network_info_list[0].network.number_of_reticulations
        if number_of_reticulations == 0:
            for network_info in network_info_list:
                network_info.network.prune()
                result.append(network_info)
            return result
        else:
            suppress_reticulation_name_iterator = [reticulation_names[index] for index in
                                                   range(index, len(reticulation_names) - (number_of_reticulations - 0) + 1)]
            if max_processes == 1:
                if progress_bar:
                    for index_2, reticulation_name in tqdm(enumerate(suppress_reticulation_name_iterator), total=len(suppress_reticulation_name_iterator),
                                                           desc="Computing displayed trees"):
                        for network_info in network_info_list:
                            new_network_info_list = cls.suppress_reticulation(network_info.network, reticulation_name)
                            result.extend(cls._displayed_trees_recursive(
                                new_network_info_list
                                , max_processes=1
                                , progress_bar=False
                                , reticulation_names=reticulation_names
                                , index=index + index_2 + 1))
                else:
                    for index_2, reticulation_name in enumerate(suppress_reticulation_name_iterator):
                        for network_info in network_info_list:
                            new_network_info_list = cls.suppress_reticulation(network_info.network, reticulation_name)
                            result.extend(cls._displayed_trees_recursive(
                                new_network_info_list
                                , max_processes=1
                                , progress_bar=False
                                , reticulation_names=reticulation_names
                                , index=index + index_2 + 1))
            else:
                pool = multiprocessing.Pool(max_processes)
                index_2 = 0
                if progress_bar:
                    for new_network_info_list in tqdm(
                            pool.imap_unordered(
                                cls._displayed_trees_recursive_helper, zip(suppress_reticulation_name_iterator, itertools.repeat(network_info_list))),
                            total=len(suppress_reticulation_name_iterator), desc=f"Computing displayed trees using {max_processes} processes"):
                        result.extend(
                            cls._displayed_trees_recursive(
                                new_network_info_list
                                , max_processes=1
                                , progress_bar=False
                                , reticulation_names=reticulation_names
                                , index=index + index_2 + 1))
                        index_2 += 1
                else:
                    for new_network_info_list in pool.imap_unordered(
                            cls._displayed_trees_recursive_helper, zip(suppress_reticulation_name_iterator, itertools.repeat(network_info_list))):
                        result.extend(cls._displayed_trees_recursive(
                            new_network_info_list
                            , max_processes=1
                            , progress_bar=False
                            , reticulation_names=reticulation_names
                            , index=index + index_2 + 1
                        ))
                        index_2 += 1
                pool.close()
                pool.join()
            return result

    @classmethod
    def _displayed_trees_recursive_helper(cls, reticulation_network_info_list):
        reticulation_name, network_info_list = reticulation_network_info_list
        result = cls(network_size=network_info_list.network_size)
        for network_info in network_info_list:
            new_network_info_list = cls.suppress_reticulation(network_info.network, reticulation_name)
            result.extend(new_network_info_list)
        return result

    @classmethod
    def suppress_reticulation(cls, network, reticulation_name):
        result = cls(network_size=network.number_of_leaves)
        reticulation_number = network.get_node_number_of(reticulation_name)
        parents_current_reticulation = network._nodes_above_nodes(children_numbers={reticulation_number}, max_height=1).difference([reticulation_number])
        for current_parent in parents_current_reticulation:
            removed_node_numbers = []
            shifted_reticulation = reticulation_number
            new_network = copy.deepcopy(network)
            for other_parent in set(parents_current_reticulation).difference([current_parent]):
                shifted_other_parent = shifted_node_number(other_parent, removed_node_numbers)
                shifted_reticulation = shifted_node_number(reticulation_number, removed_node_numbers)
                new_network._remove_arc(shifted_other_parent, shifted_node_number(shifted_reticulation, removed_node_numbers))
                in_degree, out_degree = new_network._in_degree_node(shifted_other_parent), new_network._out_degree_node(shifted_other_parent)
                if out_degree == 0:
                    removed_node_numbers.extend(new_network._terminate_node(shifted_other_parent))
                elif out_degree == 1 and in_degree == 1:
                    new_network._suppress_node(shifted_other_parent)
                    removed_node_numbers.append(shifted_other_parent)
                elif in_degree == 0 and out_degree == 1:
                    removed_node_numbers.extend(new_network._terminate_node(shifted_other_parent))
            new_network._suppress_node(shifted_node_number(shifted_reticulation, removed_node_numbers))
            result.append(NetworkInfo(new_network))
        return result

    @classmethod
    def induced_strict_tree_set(cls, network: RootedLevelKNetwork, network_size: int, max_processes: int = 1, progress_bar: bool = False):
        result = cls(network_size=network_size)
        displayed_trees = cls.displayed_trees(network, max_processes, progress_bar)
        pool = multiprocessing.Pool(max_processes)
        if progress_bar:
            for new_tree_info_list in tqdm(pool.imap_unordered(cls._induced_strict_tree_set_helper, zip(displayed_trees, itertools.repeat(network_size)))
                    , total=len(displayed_trees)
                    , desc=f"Computing exhibited trees using {max_processes} processes"):
                result.extend(new_tree_info_list)
        else:
            for new_tree_info_list in pool.imap_unordered(cls._induced_strict_tree_set_helper, zip(displayed_trees, itertools.repeat(network_size))):
                result.extend(new_tree_info_list)
        pool.close()
        pool.join()

        return result

    @classmethod
    def _induced_strict_tree_set_helper(cls, tree_network_info):
        tree_info, network_size = tree_network_info
        return cls.induced_strict_network_set(tree_info.network, network_size, 1, False)

    @classmethod
    def induced_cluster_set(cls, network: RootedLevelKNetwork, max_processes: int = 1, progress_bar: bool = False):
        result = cls(network_size=network.number_of_leaves)
        displayed_trees = cls.displayed_trees(network, max_processes, progress_bar)
        pool = multiprocessing.Pool(max_processes)
        if progress_bar:
            for new_tree_info_list in tqdm(pool.imap_unordered(cls._induced_cluster_set_recursive, displayed_trees)
                    , total=len(displayed_trees)
                    , desc=f"Computing exhibited clusters using {max_processes} processes"):
                result.extend(new_tree_info_list)
        else:
            for new_tree_info_list in pool.imap_unordered(cls._induced_cluster_set_recursive, displayed_trees):
                result.extend(new_tree_info_list)
        pool.close()
        pool.join()
        return result

    @classmethod
    def _induced_cluster_set_recursive(cls, tree_info: NetworkInfo):
        result = cls(network_size=tree_info.network.number_of_leaves, lst=[tree_info])
        sub_tree_info_list = cls._remove_tree_root(tree_info)
        for sub_tree_info in sub_tree_info_list:
            result.extend(cls._induced_cluster_set_recursive(sub_tree_info))
        return result

    @classmethod
    def _remove_tree_root(cls, tree_info: NetworkInfo):
        result = cls(network_size=tree_info.network.number_of_leaves)
        root = tree_info.network._root
        root_children = tree_info.network._nodes_below_nodes({root}, max_depth=1).difference([root])
        for root_child in root_children:
            sub_tree = RootedLevelKNetwork._get_network_below_node(tree_info.network, root_child)
            result.append(NetworkInfo(sub_tree))
        return result

    # -------------------------------- NOISE METHODS -------------------------------#
    def replace_network_distort(self, prob, network_pool):
        assert 0 <= prob <= 1, "Probability must be between 0 and 1"
        amount = np.random.binomial(len(self), prob)
        networks_to_distort = np.random.choice(range(len(self)), amount, replace=False)
        for index in networks_to_distort:
            replacement_network = copy.deepcopy(np.random.choice(network_pool))
            replacement_network.network.standardize_node_names(char_type='uid')
            for (old_name, new_name) in zip(replacement_network.network.leaf_names, self[index].network.leaf_names):
                replacement_network.network.rename_node(old_name, new_name)
            self[index] = replacement_network

    def remove_network_distort(self, prob):
        pass
        # TODO

    def tail_move_distort(self, prob):
        pass
        # TODO

    def _analyse_networks(self) -> (dict, float, float, float):
        all_leaves = self.represented_leaves()
        n_set_iterator = itertools.combinations(all_leaves, self.network_size)
        covered_networks = dict()
        for n_set in n_set_iterator:
            covered_networks[tuple(sorted(n_set))] = []

        for _, sub_network_list in self.per_leaf_set_iterator():
            total = sub_network_list.volume()
            for network_info in sub_network_list:
                taxa = tuple(sorted(network_info.network.leaf_names))
                try:
                    networks = [network_count['network'] for network_count in covered_networks[taxa]]
                except KeyError:
                    continue  # In case taxa is not in covered_networks. Happens when there are extra_taxa
                for index, network in enumerate(networks):
                    equal, _, _ = network_info.network.equal_structure(network, equal_naming=True)
                    if equal:
                        covered_networks[taxa][index]['count'] += network_info.multiplicity / total
                        break
                else:
                    covered_networks[taxa].append({'network': network_info.network, 'count': network_info.multiplicity / total})

        unique_networks = 0.
        unique_n_sets = 0.
        for taxa, network_counts in covered_networks.items():
            l = len(network_counts)
            unique_networks += l
            unique_n_sets += int(l > 0)

        max_number_of_networks = ncr(len(all_leaves), self.network_size)
        redundancy = (len(self) - unique_networks) / max_number_of_networks
        density = unique_n_sets / max_number_of_networks
        inconsistency = (unique_networks - unique_n_sets) / unique_n_sets

        return covered_networks, redundancy, density, inconsistency

    def summary(self, per_network=False) -> str:
        all_leaves = self.represented_leaves()

        covered_networks, redundancy, density, inconsistency_0 = self._analyse_networks()

        network_info_list_1 = NetworkInfoList.induced_networks_of_network_info_list(self, self.network_size - 1)
        covered_networks_1, _, _, inconsistency_1 = network_info_list_1._analyse_networks()

        summary = \
            f"NetworkInfoList ({self.uid}) the following properties for taxa \n" \
            f"        {all_leaves}: \n" \
            f" - and {len(self)} networks of size {self.network_size} \n"

        summary += f" - redundancy = {redundancy} \n" \
                   f" - density = {density} \n" \
                   f" - inconsistency-0 = {inconsistency_0} \n" \
                   f" - inconsistency-1 = {inconsistency_1} \n" \
                   f" - consistent = {bool(redundancy == 0.0 and density == 1.0 and inconsistency_0 == 0.0 and inconsistency_1 == 0.0)}"

        if per_network:
            per_network_summary = "taxa - [ network - count | ... | network - count ] \n"
            for taxa, network_counts in covered_networks.items():
                current_network_summary = [str(network_counts['network'].uid) + " - " + str(network_counts['count']) for network_counts in network_counts]
                per_network_summary += f"{taxa} : {'|'.join(current_network_summary)} \n"

            summary += f" \n " \
                       f" Per Network Summary: \n" \
                       f"{per_network_summary}"

        return summary

    def is_consistent(self) -> bool:
        _, redundancy, density, inconsistency_0 = self._analyse_networks()
        network_info_list_1 = NetworkInfoList.induced_networks_of_network_info_list(self, self.network_size - 1)
        covered_networks_1, _, _, inconsistency_1 = network_info_list_1._analyse_networks()
        return bool(redundancy == 0.0 and density == 1.0 and inconsistency_0 == 0.0 and inconsistency_1 == 0.0)

    def calculate_info(self):
        for network_info in self:
            network_info.calculate_info()

    def add_info(self, network_info_list):
        for network_info in self:
            equal_structured_networks, relation_dicts = network_info_list.find_equal_structured_network(network_info.network)
            try:
                equal_structured_network = equal_structured_networks[0]
                relation_dict = relation_dicts[0]
                network_info.add_info(equal_structured_network)
                network_info['cut_arc_sets'] = [[relation_dict[node] for node in cut_arc_set]
                                                for cut_arc_set in equal_structured_network['cut_arc_sets']]
                network_info['relation_dict'] = relation_dict
                network_info['leaf_order'] = network_info.network.leaf_ordering
            except IndexError:
                network_info.calculate_info()

    # ------------------ ADD METHODS -------------------- #
    def append(self, other: NetworkInfo):
        assert type(other) == NetworkInfo, "Can only add object of type NetworkInfo to a NetworkInfoList"
        assert len(other.network.leaf_names) <= self.network_size, \
            f"Can only add network with less or equal than {self.network_size} leaves to this NetworkInfoList, this network has {len(other.network.leaf_names)} leaves"

        for network_info in self:
            are_equal, relation_dict, leaf_translation_dicts = network_info.network.equal_structure(other.network, equal_naming=True, optimize=False)
            if are_equal:
                network_info.multiplicity += other.multiplicity
                return
        self.list.append(other)

    def extend(self, other):
        assert type(other) == NetworkInfoList, "Can only extend using a NetworkInfoList"
        for o in other:
            self.append(o)

    def __add__(self, other):
        result = copy.deepcopy(self)
        result.extend(other)
        return result

    def remove(self, other):
        self.list.remove(other)

    def pop(self, index):
        self.list.pop(index)

    # ----------------- FIND NETWORK METHODS ------------------- #
    def find_equal_structured_network(self, network):
        result = NetworkInfoList(network_size=self.network_size)
        relation_dicts = []
        for network_info in self:
            are_equal, relation_dict, leaf_translation_dicts = network_info.network.equal_structure(network, optimize=True)
            if are_equal:
                result.append(network_info)
                relation_dicts.append(leaf_translation_dicts[0])
        return result, relation_dicts

    def find_equal_network(self, network):
        for network_info in self:
            if network_info.network == network:
                return network_info
        return False

    def find_equal_leaf_names_network(self, network):
        return self.find_network_by_leaf_names(network.leaf_names)

    def find_network_by_leaf_names(self, leaf_names):
        return NetworkInfoList(self.network_size, [network_info for network_info in self.list if set(network_info.network.leaf_names) == set(leaf_names)])

    def contains_network_with_leaf_names(self, leaf_names):
        for network_info in self:
            if set(network_info.network.leaf_names) == set(leaf_names):
                return True
        return False

    def remove_network_by_leaf_names(self, leaf_names):
        network_infos = self.find_network_by_leaf_names(leaf_names)
        for network_info in network_infos:
            self.remove(network_info)
        return network_infos

    def find_networks_with_leaf_names_in(self, leaf_names):
        result = NetworkInfoList(network_size=self.network_size)
        for network_info in self:
            overlap = set(leaf_names).intersection(network_info.network.leaf_names)
            if len(overlap) >= 2:
                result.append(NetworkInfo.restrict(network_info, list(overlap)))
        return result

    def networks_where(self, category, value):
        result = NetworkInfoList(network_size=self.network_size)
        for network_info in self:
            if network_info[category] == value:
                result.append(network_info)
        return result

    def networks_with_category(self, category):
        result = NetworkInfoList(network_size=self.network_size)
        for network_info in self:
            if category in network_info.info.keys():
                result.append(network_info)
        return result

    def networks_of_size(self, size):
        result = NetworkInfoList(network_size=size)
        for network_info in self:
            if len(network_info.network.leaf_numbers) == size:
                result.append(network_info)
        return result

    def biconnected_networks(self):
        result = NetworkInfoList(network_size=self.network_size)
        for network_info in self:
            if network_info.network.biconnected:
                result.append(network_info)
        return result

    def represented_leaves(self) -> set:
        result = set()
        for network_info in self:
            result.update(network_info.network.leaf_names)
        return result

    def __getitem__(self, key):
        return self.list[key]

    def __setitem__(self, key, value):
        self.list[key] = value

    def __iter__(self):
        return iter(self.list)

    # TODO: put method in solver
    def per_leaf_set_iterator(self, method='max'):
        result = {}
        for network_info in self:
            sorted_leaf_names = tuple(sorted(network_info.network.leaf_names))
            if sorted_leaf_names in result.keys():
                result[sorted_leaf_names].append(network_info)
            else:
                result[sorted_leaf_names] = NetworkInfoList(network_size=self.network_size, lst=[network_info])
        if method == 'max':
            max_result = dict()
            for sorted_leaf_names, network_info_list in result.items():
                multiplicities = np.array([network_info.multiplicity for network_info in network_info_list])
                max_index = np.argmax(multiplicities)
                max_result[sorted_leaf_names] = NetworkInfoList(network_size=self.network_size, lst=[network_info_list[max_index]])
            return max_result.items()
        return result.items()

    def __len__(self):
        return len(self.list)

    # TODO naming
    def volume(self):
        return sum([network_info.multiplicity for network_info in self])

    def __copy__(self):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.list = copy.copy(self.list)
        cp.network_size = copy.copy(self.network_size)
        cp.uid = guid()
        cp.logger = logging.getLogger('network_info_list.{}'.format(cp.uid))
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.list = copy.deepcopy(self.list)
        cp.network_size = copy.deepcopy(self.network_size)
        cp.uid = guid()
        cp.logger = logging.getLogger('network_info_list.{}'.format(cp.uid))
        return cp

    def __getstate__(self):
        self.logger = 'network_info_list.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        self.logger = logging.getLogger(self.logger)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger)
        return self.__dict__


class Omega:
    def __init__(self, arc_weights: dict = None, node_set: set = None):
        self.node_set = coalesce(node_set, set())
        self.arc_weights = coalesce(arc_weights, dict())
        assert set(itertools.chain.from_iterable(self.arc_weights.keys())).issubset(self.node_set)
        arc_iterator = itertools.permutations(self.node_set, 2)
        for arc in arc_iterator:
            if arc not in self.arc_weights.keys():
                self.arc_weights[arc] = 0

    @classmethod
    def from_network_info_list(cls, network_info_list: NetworkInfoList):
        result = cls(node_set=network_info_list.represented_leaves())
        for _, sub_network_list in network_info_list.per_leaf_set_iterator():
            total = sub_network_list.volume()
            for network_info in sub_network_list:
                cut_arc_sets = network_info['cut_arc_sets']
                for cut_arc_set in cut_arc_sets:
                    if len(cut_arc_set) == 2:
                        x = cut_arc_set[0]
                        y = cut_arc_set[1]
                        z = [Z for Z in network_info.network.leaf_names if Z not in cut_arc_set][0]
                        for i in (x, y):
                            result.increase_weight(i, z, network_info.multiplicity / total)
        return result

    def add_node(self, node):
        assert node not in self.node_set, f"Node {node} is already in this graph"
        for other_node in self.node_set:
            self.arc_weights[(node, other_node)] = 0
            self.arc_weights[(other_node, node)] = 0
        self.node_set.add(node)

    def increase_weight(self, from_node, to_node, w: int = 1):
        assert {from_node, to_node}.issubset(self.node_set), f"Can not add weight to arc {(from_node, to_node)} as it is not in the graph"
        self.arc_weights[(from_node, to_node)] += w

    def augmented_graph(self, threshold=0):
        auxiliary_digraph = dict()
        arc_iterator = itertools.permutations(self.node_set, 2)
        for arc in arc_iterator:
            weight = self.arc_weights[arc]
            if weight <= threshold:
                if arc[0] in auxiliary_digraph:
                    auxiliary_digraph[arc[0]].add(arc[1])
                else:
                    auxiliary_digraph[arc[0]] = {arc[1]}
        return auxiliary_digraph

    def minimum_weight(self):
        arc_iterator = itertools.permutations(self.node_set, 2)
        return min([self.arc_weights[arc] for arc in arc_iterator])

    def minimal_sink_sets(self):
        threshold = self.minimum_weight()
        sub_graph = self.augmented_graph(threshold)
        strongly_connected_components = tarjan(sub_graph)
        minimal_sink_sets = self.strongly_connected_components_to_minimal_sink_sets(strongly_connected_components, sub_graph)
        if len(minimal_sink_sets) != 0:
            score = 1 - threshold / min([len(mss) for mss in minimal_sink_sets])
            return minimal_sink_sets, int(100 * score)

    @staticmethod
    def strongly_connected_components_to_minimal_sink_sets(strongly_connected_components, sub_graph):
        minimal_sink_sets = []
        for strongly_connected_component in strongly_connected_components:
            for node in strongly_connected_component:
                if node in sub_graph.keys():
                    to_leaves = sub_graph[node]
                    if not set(to_leaves).issubset(strongly_connected_component):
                        break
            else:
                minimal_sink_sets.append(strongly_connected_component)
        return [mss for mss in minimal_sink_sets if len(mss) > 1]

    def visualize(self, file_path: str = None, edge_labels=True, threshold: int = 0):
        dot = Digraph()
        dot.engine = 'dot'

        for node_name in self.node_set:
            dot.node(name=node_name, label=node_name, **{'width': str(0.5), 'height': str(0.5)})
        for arc, weight in self.arc_weights.items():
            if weight > threshold:
                continue
            if edge_labels:
                dot.edge(arc[0], arc[1], label=str(weight))
            else:
                dot.edge(arc[0], arc[1])
        if file_path:
            dot.format = 'png'
            if file_path[-4:] == '.png':
                dot.render(filename=file_path[:-4], format='png')
            else:
                dot.render(filename=file_path, format='png')
        else:
            dot.render(view=True)
            time.sleep(0.2)
