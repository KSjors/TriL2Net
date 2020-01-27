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
import typing
from config import settings
from data.trinet_format import LEVEL_1_TRINET_NAMES
import operator

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
        self.node_name_map = node_name_map  # name: number

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
    def from_rooted_level_k_generator(cls, generator, ):
        adj_matrix = copy.deepcopy(generator.adj_matrix)
        node_name_map = copy.deepcopy(generator.node_name_map)
        leaf_names = copy.deepcopy(generator.leaf_names)
        result = cls(adj_matrix, node_name_map, leaf_names=leaf_names)
        return result

    @classmethod
    def subnet_from_rooted_level_k_generator(cls, generator, leaves_per_internal_arc):
        adj_matrix = copy.deepcopy(generator.adj_matrix)
        node_name_map = copy.deepcopy(generator.node_name_map)
        leaf_names = copy.deepcopy(generator.leaf_names)
        result = cls(adj_matrix, node_name_map, leaf_names=leaf_names)
        result.standardize_node_names()
        for arc in result._internal_arcs:
            from_node, to_node = arc
            for _ in range(leaves_per_internal_arc):
                from_node, _ = result._add_leaf_to_arc((from_node, to_node))

        for leaf in result._reticulation_leaves:
            result._split_leaf(leaf)
        result.standardize_node_names()
        return result

    @classmethod
    def from_enewick(cls, string, check_valid=True):
        adjacency_dict, leaf_names = enewick(string)
        result = cls.from_connections_dict(adjacency_dict, leaf_names=list(leaf_names), check_valid=check_valid)
        return result

    @classmethod
    def from_trinet_format(cls, x):
        x = x[3:-1]
        trinet_type = x[-2:]
        x = x[:-5]
        taxa = x.split(', ')

        enewick = LEVEL_1_TRINET_NAMES[trinet_type]
        result = cls.from_enewick(enewick)
        for old_name, new_name in zip('ABC'[:len(taxa)], taxa):
            result.rename_node(old_name, new_name)

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
        network = copy.deepcopy(network)
        node_numbers_below_node = sorted(list(network._nodes_below_nodes({node_number})))
        node_names_below_node = network.get_node_names_of(node_numbers_below_node)
        node_name_map = bidict({name: number for number, name in enumerate(node_names_below_node)})
        adj_matrix = network.adj_matrix[:, node_numbers_below_node][node_numbers_below_node, :]
        leaf_names = [name for name in node_names_below_node if network.is_leaf_node(name)]
        new_network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_names=leaf_names, level=network.level, dimension=network.dimension)

        # get in-degree1 and out-degree1 nodes
        id, od = np.array(new_network.get_in_degrees()), np.array(new_network.get_out_degrees())
        nodes_to_terminate = np.where((id == 1.) & (od == 1))[0]
        new_network._terminate_nodes(set(nodes_to_terminate))

        return new_network

    @classmethod
    def random(cls, number_of_leaves: int, recombination_chance: float, level: int):
        enewick_string = '(a, b)0'
        network = cls.from_enewick(enewick_string)
        network.level = level
        network.evolve_times(number_of_leaves - 2, recombination_chance, progress_bar=False)
        network.standardize_node_names()
        return network

    @classmethod
    def from_subnets(cls, subnets, number_of_subnets, n):
        l = list(subnets.per_network_info())
        result = copy.deepcopy(random.sample(l, 1)[0]).network
        for _ in range(number_of_subnets - 1):
            subnet = random.sample(l, 1)[0].network
            leaf_to_replace = random.sample(result.leaf_numbers, 1)[0]
            result._replace_leaf_with_network(leaf_to_replace, subnet, replace_names=True)
        leaves_to_many = result.number_of_leaves - n
        assert leaves_to_many >= 0, f"Can not generate network with {n} leaves using {number_of_subnets} of these subnets"
        leaves_to_terminate = random.sample(result.leaf_numbers, leaves_to_many)
        result._terminate_nodes(set(leaves_to_terminate))
        result.standardize_node_names()
        return result

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
    def get_node_names_of(self, node_number_list: typing.Iterable[int]) -> list:
        """Retrieve names corresponding to node numbers."""
        return list(map(self.get_node_name_of, node_number_list))

    def get_node_numbers_of(self, node_name_list: typing.Iterable[str]) -> list:
        """Retrieve numbers corresponding to node names."""
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
        return self._roots[0]

    @property
    def roots(self) -> list:
        return self.get_node_names_of(self._roots)

    @property
    def _roots(self) -> list:
        return np.where(sum(self.adj_matrix > 0) == 0)[0]

    @property
    def reticulations(self) -> list:
        """Retrieve node names of all reticulations."""
        return self.get_node_names_of(self._reticulations)

    @property
    def _reticulations(self) -> list:
        return list(np.where(np.array(self.get_in_degrees()) > 1)[0])

    @property
    def reticulation_leaves(self) -> list:
        """Retrieve node names of all reticulations."""
        return self.get_node_names_of(self._reticulation_leaves)

    @property
    def _reticulation_leaves(self) -> list:
        result = []
        reticulation_nodes = self._reticulations
        for reticulation_node in self._reticulations:
            child_node = self._nodes_below_nodes({reticulation_node}, max_depth=1).difference({reticulation_node}).pop()
            if self._is_leaf_node(child_node):
                result.append(child_node)
        return result

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
        return self._o_cut_arc_sets

    @property
    def cut_arc_matrix(self) -> np.ndarray:
        """Compute indicator matrix for arcs which are cut-arcs."""
        if self._o_cut_arc_matrix is None:
            self._compute_cut_arc_matrix()
        return self._o_cut_arc_matrix

    @property
    def partial_ordering(self) -> tuple:
        return tuple(tuple(self.get_node_names_of(partial_ordering)) for partial_ordering in self._partial_ordering)

    @property
    def _partial_ordering(self) -> tuple:
        if self._o_partial_ordering is None:
            self._compute_partial_ordering()
        return self._o_partial_ordering

    @property
    def leaf_order(self) -> tuple:
        return tuple(tuple(self.get_node_names_of(leaf_ordering)) for leaf_ordering in self._leaf_order)

    @property
    def _leaf_order(self) -> tuple:
        if self._o_leaf_ordering is None:
            self._compute_leaf_ordering()
        return self._o_leaf_ordering

    @property
    def reticulationless_leaf_order(self) -> tuple:
        return tuple(tuple(self.get_node_names_of(leaf_ordering)) for leaf_ordering in self._reticulationless_leaf_order)

    @property
    def _reticulationless_leaf_order(self) -> tuple:
        reticulation_leaves = self._reticulation_leaves
        return tuple((leaf for leaf in leaf_ordering if not leaf in reticulation_leaves) for leaf_ordering in self._leaf_order)

    @property
    def number_of_arcs(self):
        mask = self.adj_matrix > 0
        return np.sum(np.multiply(self.adj_matrix, mask))

    @property
    def arcs(self) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        return [(self.get_node_name_of(from_node), self.get_node_name_of(to_node)) for from_node, to_node in self._arcs]

    @property
    def _arcs(self) -> list:
        i = 1
        from_nodes, to_nodes = [], []
        while np.any(self.adj_matrix >= i):
            extra_from_nodes, extra_to_nodes = np.where(self.adj_matrix >= i)
            to_nodes.extend(extra_to_nodes)
            from_nodes.extend(extra_from_nodes)
            i += 1
        return [(from_node, to_node) for from_node, to_node in zip(np.array(from_nodes), np.array(to_nodes))]

    # TODO: to iterators
    @property
    def internal_arcs(self) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        return [(self.get_node_name_of(from_node), self.get_node_name_of(to_node)) for from_node, to_node in self._internal_arcs]

    @property
    def _internal_arcs(self) -> list:
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

    @property
    def strict_level(self) -> int:
        in_degrees = np.array(self.get_in_degrees())

        level = 0
        biconnected_components = self._biconnected_components
        for bc in biconnected_components:
            in_sum_bc = in_degrees[bc]
            level = max(level, sum(in_sum_bc > 1))

        return level

    @property
    def movable_arcs(self):
        return [arc for arc in self.arcs if self.is_movable_arc(arc)]

    @property
    def _movable_arcs(self):
        return [arc for arc in self._arcs if self._is_movable_arc(arc)]

    def summary(self):
        result = dict()
        result['cut arc sets'] = [len(ca) for ca in self._cut_arc_sets]
        result['n'] = self.number_of_leaves
        result['biconnected components'] = {'0': {}, '1': {}, '2': {}}
        result['biconnected components leaves'] = {'0': {}, '1': {}, '2': {}}
        bcs = self._biconnected_components
        total_reticulations = 0
        total_reticulation_leaves = 0
        for bc in bcs:
            bc_reticulations = len([node for node in bc if self._is_reticulation_node(node)])
            total_reticulations += bc_reticulations
            bc_leaves = [node for node in self._nodes_below_nodes(set(bc), max_depth=1) if self._is_leaf_node(node)]
            bc_reticulation_leaves = len([node for node in bc_leaves if self._is_reticulation_leaf(node)])
            total_reticulation_leaves += bc_reticulation_leaves
            try:
                result['biconnected components'][str(bc_reticulations)][len(bc_leaves)] += 1
            except KeyError:
                result['biconnected components'][str(bc_reticulations)][len(bc_leaves)] = 1
            try:
                result['biconnected components leaves'][str(bc_reticulation_leaves)][len(bc_leaves)] += 1
            except KeyError:
                result['biconnected components'][str(bc_reticulations)][len(bc_leaves)] = 1
        result['total reticulations'] = total_reticulations
        result['total reticulation leaves'] = total_reticulation_leaves
        return result

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

    def first_unused_leaf_name(self, char_type: str = 'alph') -> str:
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
        return self._suppress_component([node_number])

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
                return True
        return False

    # --------------------------------------------------------- ALTER NETWORK METHODS --------------------------------------------------------- #
    def terminate_leaves(self, leaf_names_to_terminate: set = None, leaf_names_to_keep: set = None):
        assert bool(leaf_names_to_terminate is not None) or bool(leaf_names_to_keep is not None), "Exactly one of these parameters has to be set"
        leaf_numbers_to_terminate = set(self.get_node_numbers_of(coalesce(leaf_names_to_terminate, set())))
        leaf_numbers_to_keep = set(self.get_node_numbers_of(coalesce(leaf_names_to_keep, set())))
        self._terminate_leaves(leaf_numbers_to_terminate=leaf_numbers_to_terminate, leaf_numbers_to_keep=leaf_numbers_to_keep)
        self.reset_optimization_variables()

    def _terminate_leaves(self, leaf_numbers_to_terminate: set = None, leaf_numbers_to_keep: set = None):
        leaf_numbers_to_terminate = coalesce(leaf_numbers_to_terminate, set())
        leaf_numbers_to_keep = coalesce(leaf_numbers_to_keep, set())
        if len(leaf_numbers_to_terminate) == 0:
            leaf_numbers_to_terminate = set(self.leaf_numbers).difference(leaf_numbers_to_keep)
        self._terminate_nodes(leaf_numbers_to_terminate)

    def _terminate_leaf(self, leaf_number: int):
        """Remove leaf and prune to make a proper network again"""
        self._terminate_nodes({leaf_number})

    def terminate_leaf(self, leaf_name: str):
        """Remove leaf and prune to make a proper network again"""
        assert self.is_leaf_node(leaf_name), "Node must be a leaf"
        self._terminate_leaf(self.get_node_number_of(leaf_name))
        self.reset_optimization_variables()

    def add_leaf_to_arc(self, arc: list, leaf_name: str = None, char_type='alph') -> (str, str):
        """Add node between parent and child and attach leaf with leaf_name to it."""
        res = self.get_node_names_of(list(self._add_leaf_to_arc(self.get_node_numbers_of(arc), leaf_name, char_type)))
        self.reset_optimization_variables()
        return res

    def _add_leaf_to_arc(self, arc, leaf_name: str = None, char_type='alph') -> (int, int):
        parent = arc[0]
        child = arc[1]

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
        self.reset_optimization_variables()

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
        self.reset_optimization_variables()

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
        res = self._replace_leaf_with_network(leaf_number_to_replace, replacement_network, replace_names, char_type)
        self.reset_optimization_variables()
        return res

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

        # Replace arc from parent of leaf to leaf with arc from parent of leaf to root of replacement network
        parent_of_leaf = self._in_nodes_of_node(leaf_number_to_replace)[0]
        self._remove_arc(parent_of_leaf, leaf_number_to_replace)
        self._add_arc(parent_of_leaf, self.get_node_number_of(replacement_network_root_name + "*"))

        # Add all connections from replacement network to current network
        for internal_node_name in replacement_network_internal_node_names:
            to_node_names = replacement_network.out_nodes_of_node(internal_node_name)
            internal_node_number = self.get_node_number_of(internal_node_name + "*")
            for to_node_name in to_node_names:
                if to_node_name in replacement_network_internal_node_names or replace_names:
                    to_node_name += "*"
                self._add_arc(internal_node_number, self.get_node_number_of(to_node_name))

        # Remove leaf from network
        self._remove_node_from_network(leaf_number_to_replace)

        if replace_names:
            self.standardize_node_names(char_type)

    def tail_move(self, arc_1, arc_2):
        assert arc_1 != arc_2, "Arcs have to be unequal to each other"
        assert self.is_movable_arc(arc_1), "Can not move this arc"
        self._tail_move(self.get_node_numbers_of(arc_1), self.get_node_numbers_of(arc_2))
        self.reset_optimization_variables()

    def _tail_move(self, moving_arc, target_arc):
        assert moving_arc != target_arc, "Moving and target arc are equal"
        # Create new place for tail
        new_node_name, new_node_number = self._add_node_to_network()
        self._remove_arc(target_arc[0], target_arc[1])
        self._add_arc(target_arc[0], new_node_number)
        self._add_arc(new_node_number, target_arc[1])

        # Remove tail
        self._remove_arc(moving_arc[0], moving_arc[1])
        removed_nodes = self._terminate_nodes({moving_arc[0]})
        new_node_number = shifted_node_number(new_node_number, removed_nodes)

        # Add tail
        self._add_arc(new_node_number, shifted_node_number(moving_arc[1], removed_nodes))
        return new_node_name, new_node_number

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

        for i in self._dfs_nodes():
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
            ca_set = [node for node in self._nodes_below_nodes({to_node}) if self._is_leaf_node(node)]
            if len(ca_set) != 1:
                self._o_cut_arc_sets.append(ca_set)
        self._o_cut_arc_sets.append(list(self.leaf_numbers))

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
        self._o_partial_ordering = tuple(tuple(track) for track in self._o_partial_ordering)

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
                elif len(children) == 2:
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
                else:
                    self.visualize()
                    raise NotImplementedError
        self._o_leaf_ordering = tuple(tuple(track) for track in self._o_leaf_ordering)

    # ---------------------------------------------------------------- NODE TYPE METHODS ---------------------------------------------------------------- #
    def is_leaf_node(self, node_name: str) -> bool:
        return self._is_leaf_node(self.get_node_number_of(node_name))

    def _is_leaf_node(self, node_number: int) -> bool:
        return node_number in self.leaf_numbers

    def is_reticulation_node(self, node_name: str) -> bool:
        return self._is_reticulation_node(self.get_node_number_of(node_name))

    def _is_reticulation_node(self, node_number: int) -> bool:
        return self._in_degree_node(node_number) > 1

    def is_reticulation_leaf(self, node_name: str) -> bool:
        return self._is_reticulation_leaf(self.get_node_number_of(node_name))

    def _is_reticulation_leaf(self, node_number: int) -> bool:
        if self._is_leaf_node(node_number):
            parent_leaf = self._nodes_above_nodes({node_number}).difference({node_number}).pop()
            if self._is_reticulation_node(parent_leaf):
                return True
        return False

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

    # ---------------------------------------------------------------- ARC TYPE METHODS ---------------------------------------------------------------- #
    def is_movable_arc(self, arc):
        return self._is_movable_arc(self.get_node_numbers_of(arc))

    def _is_movable_arc(self, arc, keep_tier=False):
        from_node = arc[0]
        to_node = arc[1]
        if self._is_root_node(from_node):
            return False
        if self._is_reticulation_node(from_node):
            return False

        # Check for triangle (as defined in "Exploring the Tiers of Rooted Phylogenetic Network Space Using Tail Moves")
        if keep_tier:
            other_children_from_node = self._nodes_below_nodes({from_node}, max_depth=1).difference({from_node, to_node})
            parent_from_node = self._nodes_above_nodes({from_node}, max_height=1).difference({from_node}).pop()
            children_parent_from_node = self._nodes_below_nodes({parent_from_node}, max_depth=1).difference({parent_from_node, from_node})
            if len(children_parent_from_node) > 0 and children_parent_from_node.issubset(other_children_from_node):
                return False
        return True

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
    def _terminate_nodes(self, node_numbers_to_remove: set) -> list:
        all_removed_node_numbers = []
        while node_numbers_to_remove:
            node_number = node_numbers_to_remove.pop()
            removed_node_numbers = []

            children_numbers = self._nodes_below_nodes({node_number}, max_depth=1).difference([node_number])
            parent_numbers = self._nodes_above_nodes({node_number}, max_height=1).difference([node_number])
            if len(children_numbers) == 0:
                self._remove_node_from_network(node_number)
                removed_node_numbers.append(node_number)
                node_numbers_to_remove.update(parent_numbers)
            elif len(parent_numbers) == 0 and len(children_numbers) == 1:
                self._remove_node_from_network(node_number)
                removed_node_numbers.append(node_number)
            elif len(parent_numbers) == 1 and len(children_numbers) == 1:
                arc_added = self._suppress_node(node_number)
                removed_node_numbers.append(node_number)
                if not arc_added:
                    node_numbers_to_remove.update(children_numbers)
                    node_numbers_to_remove.update(parent_numbers)

            # Shift all nodes numbers
            node_numbers_to_remove = {shifted_node_number(node_number, removed_node_numbers) for node_number in node_numbers_to_remove}
            all_removed_node_numbers.extend(removed_node_numbers)

        return original_node_numbers(all_removed_node_numbers)

    # TODO rename
    def prune(self, suppress_redundant: str = 'all', suppress_parallel: bool = True):
        """Suppress al redundant components and parallel arcs in network."""
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
                self.reset_optimization_variables()
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
        self.reset_optimization_variables()

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

    def visualize(self, file_path: str = None, internal_node_labels=True, arc_labels=False, rankdir='TB', format='png'):
        """Visualize network."""
        dot = Digraph()
        dot.graph_attr["rankdir"] = rankdir
        dot.engine = 'dot'
        arcs = self.arcs
        for node_name in self.node_name_map:
            if internal_node_labels or self.is_leaf_node(node_name):
                dot.node(name=node_name, label=node_name, **{'width': str(0.5), 'height': str(0.5)})
            else:
                dot.node(name=node_name, label="", **{'width': str(0.1), 'height': str(0.1), 'fillcolor': 'black', 'style': 'filled'})
        for index, arc in enumerate(arcs):
            if arc_labels:
                dot.edge(arc[0], arc[1], label=str(index))
            else:
                dot.edge(arc[0], arc[1])
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

    # --------------------- COMPARE METHODS ------------------------ #
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

        # Update progressbar
        if pbar:
            pbar.n = max([len(translation_dict) for translation_dict in translation_dicts])
            pbar.refresh()

        # If at leaves no more recursion
        if self._is_leaf_node(self_current_node) and other._is_leaf_node(other_current_node):
            # In case leaf names should be equal, check if this is the case
            if equal_naming and self.get_node_name_of(self_current_node) != other.get_node_name_of(other_current_node):
                return []

            # Check if one of the homomorphisms has full size
            for translation_dict in translation_dicts:
                if len(translation_dict) == self.number_of_nodes:
                    return [translation_dict]

            # Return part homomorphisms
            return translation_dicts

        # Get children of nodes
        self_node_children = self._nodes_below_nodes({self_current_node}, max_depth=1).difference({self_current_node})
        other_node_children = other._nodes_below_nodes({other_current_node}, max_depth=1).difference({other_current_node})

        # If not same amount of children, this recursion branch fails
        if len(other_node_children) != len(self_node_children):
            return []

            # If using cut-arc sets, check if these are same size, if not, this branch fails
        if optimize or (self.optimization_available() and other.optimization_available()):
            self_cut_arc_sets = self._cut_arc_sets
            other_cut_arc_sets = other._cut_arc_sets
            self_node_cut_arc_sets = [len(cut_arc_set) for cut_arc_set in self_cut_arc_sets if self_current_node in cut_arc_set]
            other_node_cut_arc_sets = [len(cut_arc_set) for cut_arc_set in other_cut_arc_sets if other_current_node in cut_arc_set]
            if collections.Counter(self_node_cut_arc_sets) != collections.Counter(other_node_cut_arc_sets):
                return []

        # Check all bijections between the children
        bijections_iterator = itertools.permutations(range(len(self_node_children)))
        self_node_children = list(self_node_children)
        other_node_children = list(other_node_children)
        translation_dicts_to_return = []
        for bijection in bijections_iterator:
            extra_translation_dict = bidict({self_node_children[index]: other_node_children[value] for index, value in enumerate(bijection)})
            current_translation_dicts = []

            # Iterate over all homormorphisms in this branch
            for translation_dict in translation_dicts:
                # Check if new bijection is in line with homomorphism
                if self.bijection_contradicts_bijection(translation_dict, extra_translation_dict):
                    continue

                # Extend homomorphism with bijection
                current_translation_dicts.append(bidict({**translation_dict, **extra_translation_dict}))
                for new_self_node, new_other_node in extra_translation_dict.items():
                    next_translation_dicts = []

                    # Recurse for each possible branch
                    for current_translation_dict in current_translation_dicts:
                        new_dicts = self._equal_structure(new_self_node, other, new_other_node, [current_translation_dict], pbar, optimize=optimize,
                                                          equal_naming=equal_naming)

                        next_translation_dicts.extend(new_dicts)

                    # Add new branches to recurse
                    current_translation_dicts = next_translation_dicts

                # In case there is a full homomorphism, finish
                for current_translation_dict in current_translation_dicts:
                    if len(current_translation_dict) == self.number_of_nodes:
                        return [current_translation_dict]

            # Add finishes branches
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

    def level_1_trinet_name(self, level_1_trinets=None, leaf_name_sep=' '):
        if level_1_trinets is None:
            level_1_trinets = {name: RootedLevelKNetwork.from_enewick(eNewick_string) for name, eNewick_string in LEVEL_1_TRINET_NAMES.items()}
        for name, level_1_trinet in level_1_trinets.items():
            equal, _, leaf_translation_dicts = level_1_trinet.equal_structure(self)
            if equal:
                leaf_translation_dict = leaf_translation_dicts[0]
                return name, f"{leaf_translation_dict['A']}{leaf_name_sep}{leaf_translation_dict['B']}{leaf_name_sep}{leaf_translation_dict['C']}"
        return None, None

    def cut_arc_set_consistency(self, other):
        cut_arc_sets = {tuple(sorted(ca_set)) for ca_set in self.cut_arc_sets}
        other_cut_arc_sets = {tuple(sorted(ca_set)) for ca_set in other.cut_arc_sets}
        return sum([cut_arc_set in other_cut_arc_sets for cut_arc_set in cut_arc_sets]) / len(cut_arc_sets)

    # ------------------- EVOLVE METHODS ------------------ #
    def evolve_times(self, times: int, recombination_chance: float = 0, progress_bar=False):
        if progress_bar:
            for _ in tqdm(range(times), desc="Evolving network"):
                self.evolve(recombination_chance)
        else:
            for _ in range(times):
                self.evolve(recombination_chance)

        # Remove reticulations from biconnected components that are above level at random
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
                    self._terminate_nodes({reticulation_to_remove})
                    break

        self.prune()

    def evolve(self, recombination_chance: float = 0) -> (bool, str):
        # Note: does not keep level intact perfectly
        assert 0 <= recombination_chance < 1, "Recombincation chance not between 0 and 1"
        if np.random.uniform() < recombination_chance:
            leaf_1, leaf_2 = random.sample(self.leaf_numbers, 2)
            self._recombine_leaves(leaf_1, leaf_2)
        else:
            leaf_number = random.choice(self.leaf_numbers)
            self._split_leaf(leaf_number)

        self.reset_optimization_variables()

    def terminate_percentage_leaves(self, termination_percentage: int = 0.1):
        leaves_to_terminate = random.sample(self.leaf_numbers, int(self.number_of_leaves * termination_percentage))
        self._terminate_nodes(set(leaves_to_terminate))
        # self.prune()
        self.standardize_node_names()

    def enewick(self):
        reticulation_name_map = {reticulation_name: f'#H{index}' for index, reticulation_name in enumerate(self.reticulations)}
        return self._enewick(self._root, reticulation_name_map=reticulation_name_map) + 'root;'

    def _enewick(self, current_node, traversed_nodes=None, reticulation_name_map=None):
        traversed_nodes = coalesce(traversed_nodes, [])

        # add tag to reticulations
        node_name = self.get_node_name_of(current_node)
        if reticulation_name_map and node_name in reticulation_name_map:
            node_name_tag = reticulation_name_map[node_name]
        elif node_name in self.leaf_names:
            node_name_tag = node_name
        else:
            node_name_tag = ''

        # Stop if traversed
        if current_node in traversed_nodes:
            return node_name_tag

        # update traversed node list
        traversed_nodes.append(current_node)

        # Recurse children
        current_node_children = self._nodes_below_nodes({current_node}, max_depth=1).difference([current_node])
        if len(current_node_children) == 0:
            return node_name_tag
        elif len(current_node_children) == 1:
            child = current_node_children.pop()
            return "(" + self._enewick(child, traversed_nodes, reticulation_name_map=reticulation_name_map) + ")" + node_name_tag
        elif len(current_node_children) == 2:
            child_1 = current_node_children.pop()
            child_2 = current_node_children.pop()
            return "(" + self._enewick(child_1, traversed_nodes, reticulation_name_map=reticulation_name_map) + "," + self._enewick(child_2, traversed_nodes,
                                                                                                                                     reticulation_name_map=reticulation_name_map) + ")" + node_name_tag
        raise ValueError

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

    def dfs_nodes(self, current_node=None):
        for node_number in self._dfs_nodes(current_node):
            yield self.get_node_name_of(node_number)

    def _dfs_nodes(self, current_node=None):
        current_node = coalesce(current_node, self._root)
        yield current_node
        for child in self._nodes_below_nodes({current_node}, max_depth=1).difference({current_node}):
            yield from self._dfs_nodes(child)


class RootedLevelKGenerator(RootedLevelKNetwork):
    def __init__(self, name, dir_adj_matrix: np.ndarray, symmetrical_nodes: bidict, level: int = 2, dimension: int = 2, check_valid: bool = True,
                 char_type='ALPH'):
        network = RootedLevelKNetwork.from_dir_adj_matrix(dir_adj_matrix=dir_adj_matrix, level=level, dimension=dimension, check_valid=check_valid,
                                                          char_type=char_type)
        super().__init__(network.adj_matrix, network.node_name_map, leaf_numbers=network.leaf_numbers, level=network.level, dimension=network.dimension)
        self.name = name
        self.symmetrical_nodes = symmetrical_nodes

    def build_trinets(self):
        base_net = self

        reticulations = self.leaves_below_nodes(set(self.reticulations), 1)

        # --------- Trinets -----------
        # Create iterator of possible combinations of leaves to add
        internal_arcs = base_net.internal_arcs
        number_of_generator_leaves = len(base_net.leaf_numbers)
        internal_arc_iterator = itertools.combinations(internal_arcs, 3 - number_of_generator_leaves)

        # For each possible combination, create trinet and save it to trinets_gen_sides list
        trinet_info_list = NetworkSet(network_size=3)

        for arcs in internal_arc_iterator:
            extra_leaf_dict = {}
            # Add extra leaves to base net and save together with underlying generator (self) and added edges
            current_trinet = RootedLevelKNetwork.from_rooted_level_k_generator(self)
            for arc in arcs:
                _, leaf_name = current_trinet.add_leaf_to_arc(arc, char_type='ALPH')
                extra_leaf_dict[leaf_name] = arc
            current_trinet.prune()
            if current_trinet.number_of_reticulations != self.level:
                continue
            current_trinet.reset_optimization_variables()
            trinet_info = NetworkInfo(current_trinet, {'generator'      : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                       'extra_leaf_dict': copy.deepcopy(extra_leaf_dict), 'symmetrical_nodes': self.symmetrical_nodes})
            trinet_info_list.append(trinet_info)

        # In case generator has only one leaf, also add two leaves to the same edge (only need to do one node of every symmetry pair)
        if number_of_generator_leaves == 1:
            extra_leaf_dict = {}
            internal_arc_iterator = itertools.combinations(internal_arcs, 1)
            for arc in internal_arc_iterator:
                current_trinet = RootedLevelKNetwork.from_rooted_level_k_generator(self)
                new_node_name, leaf_name_1 = current_trinet.add_leaf_to_arc(arc[0], char_type='ALPH')
                extra_leaf_dict[leaf_name_1] = arc[0]
                _, leaf_name_2 = current_trinet.add_leaf_to_arc([new_node_name, arc[0][1]], char_type='ALPH')
                extra_leaf_dict[leaf_name_2] = arc[0]

                current_trinet.prune()
                if current_trinet.number_of_reticulations != self.level:
                    continue
                current_trinet.reset_optimization_variables()
                trinet_info = NetworkInfo(current_trinet, {'generator'      : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                           'extra_leaf_dict': copy.deepcopy(extra_leaf_dict), 'symmetrical_nodes': self.symmetrical_nodes})
                trinet_info_list.append(trinet_info)
        return trinet_info_list

    def build_binets(self):
        base_net = self

        reticulations = self.leaves_below_nodes(set(self.reticulations), 1)

        # --------- Binets -----------
        # Create iterator of possible combinations of leaves to add
        internal_arcs = base_net.internal_arcs
        number_of_generator_leaves = len(base_net.leaf_numbers)
        internal_arc_iterator = itertools.combinations(internal_arcs, 2 - number_of_generator_leaves)

        # For each possible combination, create binet and save it to trinets_gen_sides list
        binet_info_list = NetworkSet(network_size=2)
        for arcs in internal_arc_iterator:
            extra_leaf_dict = {}
            # Add extra leaves to base net and save together with underlying generator (self) and added edges
            current_binet = RootedLevelKNetwork.from_rooted_level_k_generator(self)
            for arc in arcs:
                _, leaf_name = current_binet.add_leaf_to_arc(arc, char_type='ALPH')
                extra_leaf_dict[leaf_name] = arc
            current_binet.prune()
            if current_binet.number_of_reticulations != self.level:
                continue
            current_binet.reset_optimization_variables()
            current_binet.calculate_optimization_variables()
            binet_info = NetworkInfo(current_binet, {'generator'      : self, 'generator_name': self.name, 'reticulations': reticulations,
                                                     'extra_leaf_dict': copy.deepcopy(extra_leaf_dict), 'symmetrical_nodes': self.symmetrical_nodes})
            binet_info_list.append(binet_info)
        return binet_info_list

    @property
    def sides(self):
        return self.arc_sides + self.reticulation_sides

    @property
    def arc_sides(self):
        return self.internal_arcs

    @property
    def reticulation_sides(self):
        reticulation_sides = []
        reticulations = self.reticulations
        for reticulation in reticulations:
            reticulation_child = self.nodes_below_nodes({reticulation}, max_depth=1).difference({reticulation}).pop()
            if self.is_leaf_node(reticulation_child):
                reticulation_sides.append(reticulation)
        return reticulation_sides

    @property
    def sets_of_symmetric_arc_sides(self):
        # TODO: make general
        symmetric_sides_sets_per_generator = {
            '1' : {
                ('0', '1'): [('0', '1'), ('0', '1')]
            },
            '2a': {
                ('0', '1'): [('0', '1')],
                ('0', '2'): [('0', '2')],
                ('1', '2'): [('1', '2')],
                ('1', '3'): [('1', '3')],
                ('2', '3'): [('2', '3')]
            },
            '2b': {
                ('0', '1'): [('0', '1')],
                ('0', '2'): [('0', '2')],
                ('1', '3'): [('1', '3')],
                ('1', '4'): [('1', '4')],
                ('3', '2'): [('3', '2')],
                ('3', '4'): [('3', '4')]
            },
            '2c': {
                ('0', '1'): [('0', '1'), ('0', '2')],
                ('1', '3'): [('1', '3'), ('2', '3')],
                ('1', '4'): [('1', '4'), ('2', '4')]
            },
            '2d': {
                ('0', '1'): [('0', '1')],
                ('0', '2'): [('0', '2')],
                ('1', '3'): [('1', '3'), ('1', '3')],
                ('3', '2'): [('3', '2')]
            }
        }
        return symmetric_sides_sets_per_generator[self.name]

    @property
    def sets_of_symmetric_reticulation_sides(self):
        return {key: [key] for key in self.reticulation_sides}

    @property
    def sets_of_symmetric_sides(self):
        return {**self.sets_of_symmetric_arc_sides, **self.sets_of_symmetric_reticulation_sides}

    @property
    def crucial_sets_of_symmetric_sides(self):
        result = []
        for key, sides in self.sets_of_symmetric_sides:
            if len(sides) > 1 and sides[0] == sides[1]:
                result.append(key)
        return result

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
        cp.name = copy.copy(self.name)
        cp.symmetrical_nodes = copy.copy(self.symmetrical_nodes)
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
        cp.name = copy.copy(self.name)
        cp.symmetrical_nodes = copy.deepcopy(self.symmetrical_nodes)
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp._o_cut_arc_matrix = copy.deepcopy(self._o_cut_arc_matrix)
        cp._o_cut_arc_sets = copy.deepcopy(self._o_cut_arc_sets)
        cp._o_biconnected_components = copy.deepcopy(self._o_biconnected_components)
        cp._o_partial_ordering = copy.deepcopy(self._o_partial_ordering)
        cp._o_leaf_ordering = copy.deepcopy(self._o_leaf_ordering)
        return cp


""" ############################################################################################################################################################
                                                                      NETWORK INFO
############################################################################################################################################################ """


class NetworkInfo:
    def __init__(self, network: RootedLevelKNetwork, info: dict = None, multiplicity: int = 1):
        assert type(network) != self.__class__, "Creating network info on wrong type"
        # TODO--> uncomment and build network from generator
        # Set up logging
        self.uid = guid()
        self.logger = logging.getLogger('network_info.{}'.format(self.uid))
        self.logger.debug(f"Created NetworkInfo ({self.uid})")

        self.network = network
        self.info = coalesce(info, dict())
        assert type(self.info) == dict, "Info dict of network info has wrong type"
        self.multiplicity = multiplicity

    def name(self):
        return tuple(sorted(self.network.leaf_names))

    def calculate_info(self):
        self['cut_arc_sets'] = self.network.cut_arc_sets
        self['strict_level'] = self.network.strict_level
        self['reticulationless_leaf_order'] = self.network.reticulationless_leaf_order
        self['leaf_order'] = self.network.leaf_order
        self['number_of_reticulations'] = self.network.number_of_reticulations
        self['biconnected'] = bool(len(self['cut_arc_sets']) == 1)

    def add_structure_info(self, network_info, relation_dict):
        self.info['relation_dict'] = relation_dict
        for key in ['generator', 'generator_name', 'reticulations', 'extra_leaf_dict', 'symmetrical_nodes']:
            value = network_info[key]
            if key == 'extra_leaf_dict':
                self.info[key] = {relation_dict[leaf_name]: arc for leaf_name, arc in value.items()}
            else:
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
        cp.uid = guid()
        cp.logger = logging.getLogger('network_info.{}'.format(cp.uid))
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.network = copy.deepcopy(self.network)
        cp.info = copy.deepcopy(self.info)
        cp.multiplicity = copy.copy(self.multiplicity)
        cp.uid = guid()
        cp.logger = logging.getLogger('network_info.{}'.format(cp.uid))
        return cp

    def __str__(self):
        return str(self.network.leaf_names) + "\n" + pp.pformat(self.info)

    def __getstate__(self):
        self.logger = 'network_info.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        self.logger = logging.getLogger(self.logger)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger)
        return self.__dict__


""" ############################################################################################################################################################
                                                                  NETWORK INFO LIST
############################################################################################################################################################ """


class NetworkSet:
    def __init__(self, network_size: int = None, network_list: list = None, network_info_list: list = None, equal_naming: bool = True):
        # Set up logging
        self.uid = guid()
        self.logger = logging.getLogger('network_set.{}'.format(self.uid))
        self.logger.debug(f"Created NetworkSet ({self.uid})")

        # Set up structure
        self.dictionary = dict()  # leaf_names: {volume: _, network_info_set: {_, _, ...} }
        self.equal_naming = equal_naming

        # Fill structure
        self.network_size = network_size
        network_list = coalesce(network_list, [])
        for network in network_list:
            self.append(NetworkInfo(network))

        network_info_list = coalesce(network_info_list, [])
        for network_info in network_info_list:
            self.append(network_info)

    @classmethod
    def networks_with_reticulations(cls, network_set, reticulation_leaves):
        result = cls(network_set.network_size)
        for network_info, weight in network_set.per_network_iterator(method=settings.WEIGHTED_SUM):
            if set(network_info.network.reticulation_leaves).issubset(reticulation_leaves):
                result.append(network_info)
        return result

    @classmethod
    def from_enewick(cls, enewick, network_size):
        return cls(network_info_list=[NetworkInfo(RootedLevelKNetwork.from_enewick(ene), multiplicity=multiplicity) for ene, multiplicity in enewick.items()],
                   network_size=network_size)

    @classmethod
    def collapse_network(cls, network: RootedLevelKNetwork, leaf_set: set, multiplicity: int = 1):
        """Returns list of shrunken networks"""
        result = cls()
        intersection = set(network.leaf_names).intersection(leaf_set)
        if len(intersection) == 0:
            result.append(NetworkInfo(copy.deepcopy(network), multiplicity=multiplicity))
        elif network.number_of_leaves - len(intersection) >= 1:
            leaf_numbers_set = set(network.get_node_numbers_of(list(intersection)))
            mss_name = mss_leaf_name(leaf_set)
            for leaf_to_keep in leaf_numbers_set:
                new_network = copy.deepcopy(network)
                new_network._rename_node(leaf_to_keep, mss_name)
                new_network._terminate_nodes(leaf_numbers_set.difference([leaf_to_keep]))
                new_network.prune()
                result.append(NetworkInfo(new_network, multiplicity=multiplicity))
        return result

    @classmethod
    def collapse_network_info_list(cls, network_info_list, leaf_set: set):
        result = cls()
        for network_info in network_info_list.per_network_info():
            result.extend(cls.collapse_network(network_info.network, leaf_set, network_info.multiplicity))
        return result

    # ---------------- 0
    @classmethod
    def induced_tree_set_of_network_set(cls, network_set, max_processes: int = 1, progress_bar: bool = False):
        result = cls()
        if progress_bar:
            for network_info in tqdm(network_set.per_network_info()):
                induced_trees = cls.induced_strict_tree_set(network_info.network, network_size=network_info.network.number_of_leaves,
                                                            max_processes=max_processes,
                                                            multiplicity=network_info.multiplicity)
                result.extend(induced_trees)
        else:
            for network_info in network_set.per_network_info():
                induced_trees = cls.induced_strict_tree_set(network_info.network, network_size=network_info.network.number_of_leaves,
                                                            max_processes=max_processes,
                                                            multiplicity=network_info.multiplicity)
                result.extend(induced_trees)
        return result

    # ---------------- 1
    @classmethod
    def induced_network_set_of_network_set(cls, network_info_list, network_size: int, progress_bar: bool = False):
        """ Network list containing all exhibited networks of networks in list of certain size """
        assert network_info_list.network_size is None or network_size <= network_info_list.network_size, "Can only create smaller network info lists with smaller network size."
        result = cls(network_size)

        if progress_bar:
            for network_info, weight in tqdm(network_info_list.per_network_iterator(method=settings.WEIGHTED_SUM)):
                result.extend(
                    cls.induced_strict_network_set(network_info.network, network_size, 1, progress_bar=False, multiplicity=weight))
        else:
            for network_info, weight in network_info_list.per_network_iterator(method=settings.WEIGHTED_SUM):
                result.extend(
                    cls.induced_strict_network_set(network_info.network, network_size, 1, progress_bar=False, multiplicity=weight))
        return result

    # ---------------- 1.1
    @classmethod
    def induced_strict_network_set(cls, network: RootedLevelKNetwork, network_size: int, max_processes: int = 1, progress_bar: bool = False,
                                   method=settings.ITERATIVE, multiplicity: int = 1):
        """ Network list containing all exhibited network of network of certain size """
        if method == settings.RECURSIVE:
            return cls._induced_strict_network_set_recursive(network, network_size, max_processes, progress_bar, network.leaf_names, 0, multiplicity)
        elif method == settings.ITERATIVE:
            return cls.induced_strict_network_set_iterative(network, network_size, max_processes, progress_bar, multiplicity)
        else:
            raise ValueError

    # ---------------- 1.1.1
    @classmethod
    def _induced_strict_network_set_recursive(cls, network: RootedLevelKNetwork, network_size: int, max_processes: int, progress_bar: bool, leaf_names: list,
                                              index: int, multiplicity: int = 1):
        raise NotImplementedError
        result = cls(network_size=network_size)
        number_of_leaves = network.number_of_leaves
        if number_of_leaves == network_size:
            # network.prune()
            result.append(NetworkInfo(network, multiplicity=multiplicity))
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
                                                                      index=index + index_2 + 1, multiplicity=multiplicity))
                else:
                    for index_2, terminate_leaf_name in enumerate(terminate_leaf_name_iterator):
                        new_network = copy.deepcopy(network)
                        new_network.terminate_leaf(terminate_leaf_name)
                        result.extend(
                            cls._induced_strict_network_set_recursive(new_network, network_size, max_processes=1, progress_bar=False, leaf_names=leaf_names,
                                                                      index=index + index_2 + 1, multiplicity=multiplicity))
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
                            , index=index + index_2 + 1
                            , multiplicity=multiplicity))
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
                            , index=index + index_2 + 1
                            , multiplicity=multiplicity))
                        index_2 += 1
                pool.close()
                pool.join()
            return result

    # ---------------- 1.1.1.1
    @classmethod
    def _induced_strict_network_set_recursive_helper(cls, leaf_network):
        terminate_leaf_name, network = leaf_network
        new_network = copy.deepcopy(network)
        new_network.terminate_leaf(terminate_leaf_name)
        return new_network

    # ---------------- 1.1.2
    @classmethod
    def induced_strict_network_set_iterative(cls, network: RootedLevelKNetwork, network_size: int, max_processes: int = 1, progress_bar: bool = False,
                                             multiplicity: int = 1):
        result = cls(network_size=network_size)
        n_set_iterator = itertools.combinations(network.leaf_numbers, network_size)
        if max_processes > 1:
            pool = multiprocessing.Pool(max_processes)
            if progress_bar:
                for new_network in tqdm(
                        pool.imap_unordered(cls._induced_strict_network_set_iterative_helper, zip(n_set_iterator, itertools.repeat(network))),
                        total=ncr(network.number_of_leaves, network_size), desc=f"Computing exhibited trinets iteratively using {max_processes} processes"):
                    result.append(NetworkInfo(new_network, multiplicity=multiplicity))
            else:
                for new_network in pool.imap_unordered(cls._induced_strict_network_set_iterative_helper, zip(n_set_iterator, itertools.repeat(network))):
                    result.append(NetworkInfo(new_network, multiplicity=multiplicity))
            pool.close()
            pool.join()
        else:
            if progress_bar:
                for n_set in tqdm(n_set_iterator, total=ncr(network.number_of_leaves, network_size),
                                  desc=f"Computing exhibited trinets iteratively using {max_processes} processes"):
                    new_network = RootedLevelKNetwork._restrict(network, n_set)
                    result.append(NetworkInfo(new_network, multiplicity=multiplicity))
            else:
                for n_set in n_set_iterator:
                    new_network = RootedLevelKNetwork._restrict(network, n_set)
                    result.append(NetworkInfo(new_network, multiplicity=multiplicity))
        return result

    # ---------------- 1.1.2.1
    @classmethod
    def _induced_strict_network_set_iterative_helper(cls, n_set_network):
        n_set, network = n_set_network
        return RootedLevelKNetwork._restrict(network, n_set)

    # ---------------- 2
    @classmethod
    def displayed_trees(cls, network: RootedLevelKNetwork, max_processes: int = 1, progress_bar: bool = False, multiplicity: int = 1):
        network_info_list = cls()
        network_info_list.append(NetworkInfo(network, multiplicity=multiplicity))
        return cls._displayed_trees_recursive(network_info_list, network.reticulations, max_processes, progress_bar, )

    # ---------------- 2.1
    @classmethod
    def _displayed_trees_recursive(cls, network_info_list, reticulation_names: list, max_processes: int = 1, progress_bar: bool = False):
        result = cls(network_size=network_info_list.network_size)
        if len(reticulation_names) == 0:
            for network_info in network_info_list.per_network_info():
                network_info.network.prune()
            return network_info_list
        else:
            if max_processes == 1:
                if progress_bar:
                    for network_info in network_info_list.per_network_info():
                        new_network_info_list = cls.suppress_reticulation(network_info.network, reticulation_names[0], network_info.multiplicity)
                        result.extend(cls._displayed_trees_recursive(new_network_info_list, reticulation_names=reticulation_names[1:]))
                else:
                    for network_info in network_info_list.per_network_info():
                        new_network_info_list = cls.suppress_reticulation(network_info.network, reticulation_names[0], network_info.multiplicity)
                        result.extend(cls._displayed_trees_recursive(new_network_info_list, reticulation_names=reticulation_names[1:]))
            else:
                pool = multiprocessing.Pool(max_processes)
                index_2 = 0
                if progress_bar:
                    for new_network_info_list in tqdm(
                            pool.imap_unordered(
                                cls._displayed_trees_recursive_helper, zip(network_info_list.per_network_info(), itertools.repeat(reticulation_names[0]))),
                            desc=f"Computing displayed trees using {max_processes} processes"):
                        result.extend(cls._displayed_trees_recursive(new_network_info_list, reticulation_names=reticulation_names[1:]))
                        index_2 += 1
                else:
                    for new_network_info_list in pool.imap_unordered(cls._displayed_trees_recursive_helper,
                                                                     zip(network_info_list.per_network_info(), itertools.repeat(reticulation_names[0]))):
                        result.extend(cls._displayed_trees_recursive(new_network_info_list, reticulation_names=reticulation_names[1:]))
                        index_2 += 1
                pool.close()
                pool.join()
        return result

    # ---------------- 2.1.1
    @classmethod
    def _displayed_trees_recursive_helper(cls, network_info_reticulation):
        network_info, reticulation_name = network_info_reticulation
        return cls.suppress_reticulation(network_info.network, reticulation_name, network_info.multiplicity)

    # ---------------- 2.1.1.1
    @classmethod
    def suppress_reticulation(cls, network: RootedLevelKNetwork, reticulation_name, multiplicity: int = 1):
        result = cls()

        # Check if reticulation has been removed
        if reticulation_name not in network.node_name_map:
            new_network = copy.deepcopy(network)
            result.append(NetworkInfo(new_network, multiplicity=multiplicity * 2))  # binary --> times 2
            return result
        reticulation_number = network.get_node_number_of(reticulation_name)

        # If node is not a reticulation, do nothing
        if not network._is_reticulation_node(reticulation_number):
            new_network = copy.deepcopy(network)
            result.append(NetworkInfo(new_network, multiplicity=multiplicity * 2))
            return result

        # Else, create a network where only the edge between the reticulation and one of its parents is kept
        parents_current_reticulation = network._nodes_above_nodes(children_numbers={reticulation_number}, max_height=1).difference([reticulation_number])
        # Pick a parent
        for current_parent in parents_current_reticulation:
            removed_node_numbers = []
            shifted_reticulation = reticulation_number
            new_network = copy.deepcopy(network)
            # Remove edges to other parents
            for other_parent in set(parents_current_reticulation).difference([current_parent]):
                # Shift node numbers
                shifted_other_parent = shifted_node_number(other_parent, removed_node_numbers)
                shifted_reticulation = shifted_node_number(reticulation_number, removed_node_numbers)
                new_network._remove_arc(shifted_other_parent, shifted_node_number(shifted_reticulation, removed_node_numbers))
                in_degree, out_degree = new_network._in_degree_node(shifted_other_parent), new_network._out_degree_node(shifted_other_parent)
                # If other parent has out_degree 0, can be suppressed, or is root with out_degree 1, it has to be terminated
                if out_degree == 0:
                    removed_node_numbers.extend(new_network._terminate_nodes({shifted_other_parent}))
                elif out_degree == 1 and in_degree == 1:
                    new_network._suppress_node(shifted_other_parent)
                    removed_node_numbers.append(shifted_other_parent)
                elif in_degree == 0 and out_degree == 1:
                    removed_node_numbers.extend(new_network._terminate_nodes({shifted_other_parent}))
            new_network._suppress_node(shifted_node_number(shifted_reticulation, removed_node_numbers))
            new_network.reset_optimization_variables()
            result.append(NetworkInfo(new_network, multiplicity=multiplicity))
        return result

    # ----------------- 3
    @classmethod
    def induced_strict_tree_set(cls, network: RootedLevelKNetwork, network_size: int, max_processes: int = 1, progress_bar: bool = False,
                                method: int = settings.ITERATIVE, multiplicity: int = 1):
        network_set = NetworkSet.induced_strict_network_set(network, network_size, max_processes, progress_bar, method=method, multiplicity=multiplicity)
        return cls.induced_strict_tree_set_of_network_set(network_set, network_size, max_processes, progress_bar)

    # ----------------- 3.1
    @classmethod
    def induced_strict_tree_set_of_network_set(cls, network_set, network_size: int, max_processes: int = 1, progress_bar: bool = False):
        network_set_iterator = network_set.per_network_info()
        result = cls(network_size=network_size)
        if max_processes > 1:
            pool = multiprocessing.Pool(max_processes)
            if progress_bar:
                for new_tree_info_list in tqdm(
                        pool.imap_unordered(cls._induced_strict_tree_set_helper, network_set_iterator)
                        , desc=f"Computing exhibited trees using {max_processes} processes"):
                    result.extend(new_tree_info_list)
            else:
                for new_tree_info_list in pool.imap_unordered(cls._induced_strict_tree_set_helper, network_set_iterator):
                    result.extend(new_tree_info_list)
            pool.close()
            pool.join()
        else:
            if progress_bar:
                for tree_info in tqdm(network_set_iterator, desc=f"Computing exhibited trees using {max_processes} processes"):
                    result.extend(cls.displayed_trees(tree_info.network, multiplicity=tree_info.multiplicity))
            else:
                for tree_info in network_set_iterator:
                    result.extend(cls.displayed_trees(tree_info.network, multiplicity=tree_info.multiplicity))
        return result

    # ----------------- 3.1.1
    @classmethod
    def _induced_strict_tree_set_helper(cls, tree_info):
        return cls.displayed_trees(tree_info.network, multiplicity=tree_info.multiplicity, progress_bar=False)

    # ----------------- 4
    @classmethod
    def induced_cluster_set(cls, network: RootedLevelKNetwork, max_processes: int = 1, progress_bar: bool = False):
        result = cls()
        displayed_trees = cls.displayed_trees(network, max_processes, progress_bar).per_network_info()
        pool = multiprocessing.Pool(max_processes)
        if progress_bar:
            for new_tree_info_list in tqdm(pool.imap_unordered(cls._induced_cluster_set_recursive, displayed_trees)
                    , desc=f"Computing exhibited clusters using {max_processes} processes"):
                result.extend(new_tree_info_list)
        else:
            for new_tree_info_list in pool.imap_unordered(cls._induced_cluster_set_recursive, displayed_trees):
                result.extend(new_tree_info_list)
        pool.close()
        pool.join()
        # result.set_multiplicities_to_one()
        return result

    # ----------------- 4.1
    @classmethod
    def _induced_cluster_set_recursive(cls, tree_info):
        result = cls(network_info_list=[tree_info])
        sub_tree_info_list = cls._remove_tree_root(tree_info)
        for sub_tree_info in sub_tree_info_list.per_network_info():
            result.extend(cls._induced_cluster_set_recursive(sub_tree_info))
        return result

    # ----------------- 4.1.1
    @classmethod
    def _remove_tree_root(cls, tree_info: NetworkInfo):
        result = cls()
        root = tree_info.network._root
        root_children = tree_info.network._nodes_below_nodes({root}, max_depth=1).difference([root])
        for root_child in root_children:
            sub_tree = RootedLevelKNetwork._get_network_below_node(tree_info.network, root_child)
            if sub_tree.number_of_leaves >= 2:
                result.append(NetworkInfo(sub_tree, multiplicity=tree_info.multiplicity))
        return result

    # ---------------------..
    @classmethod
    def subnets_from_generators(cls, generators, leaves_per_internal_arc):
        result = cls()
        for generator in generators:
            result.append(NetworkInfo(RootedLevelKNetwork.subnet_from_rooted_level_k_generator(generator, leaves_per_internal_arc)))
        return result

    # -------------------------------- NOISE METHODS -------------------------------#
    @classmethod
    def distort(cls, network_set, Q_u, Q_tm, Q_d, network_info_pool, max_replacement_level):
        assert Q_u + Q_tm + Q_d <= 1
        result = cls(network_size=network_set.network_size)
        max_replacement_level_pool = cls.networks_where(network_info_pool, 'number_of_reticulations', max_replacement_level, lambda x, y: x <= y)
        for network_info in network_set.per_network_info():
            Q = random.uniform(0, 1)
            if Q <= Q_u:
                random_network_info = np.random.choice(list(max_replacement_level_pool.per_network_info()))
                replacement_network_info = copy.deepcopy(random_network_info)
                leaf_names = network_info.network.leaf_names
                replacement_leaf_names = replacement_network_info.network.leaf_names
                for replacement_leaf_name, leaf_name in zip(replacement_leaf_names, leaf_names):
                    replacement_network_info.network.rename_node(replacement_leaf_name, leaf_name)
                result.append(replacement_network_info)
            elif Q <= Q_u + Q_tm:
                network = copy.deepcopy(network_info.network)
                while True:
                    moving_arc = random.sample(network.movable_arcs, 1)[0]
                    nodes_below_moving_arc = network.nodes_below_nodes({moving_arc[0]})
                    all_arcs = network.arcs
                    all_arcs.remove(moving_arc)
                    all_arcs = [arc for arc in all_arcs if
                                arc[1] != moving_arc[1] and
                                arc[0] not in nodes_below_moving_arc]
                    if network.is_tree_node(moving_arc[0]):
                        parent_moving_arc = network.nodes_above_nodes({moving_arc[0]}, max_height=1).difference({moving_arc[0]}).pop()
                        all_arcs.remove((parent_moving_arc, moving_arc[0]))
                    if len(all_arcs) > 0:
                        target_arc = random.sample(all_arcs, 1)[0]
                        network.tail_move(moving_arc, target_arc)
                        if network.strict_level > network.level:
                            network = copy.deepcopy(network_info.network)
                            continue
                        result.append(NetworkInfo(network, multiplicity=network_info.multiplicity))
                        break
            elif Q <= Q_u + Q_tm + Q_d:
                pass
            else:
                result.append(copy.deepcopy(network_info))
        return result

    @classmethod
    def uniform_distort(cls, network_set, prob, network_info_pool):
        assert 0 <= prob <= 1, "Probability must be between 0 and 1"
        result = cls(network_size=network_set.network_size)
        for network_info in network_set.per_network_info():
            if random.uniform(0, 1) <= prob:
                random_network_info = np.random.choice(list(network_info_pool.per_network_info()))
                replacement_network_info = copy.deepcopy(random_network_info)
                leaf_names = network_info.network.leaf_names
                replacement_leaf_names = replacement_network_info.network.leaf_names
                for replacement_leaf_name, leaf_name in zip(replacement_leaf_names, leaf_names):
                    replacement_network_info.network.rename_node(replacement_leaf_name, leaf_name)
                result.append(replacement_network_info)
            else:
                result.append(copy.deepcopy(network_info))
        return result

    @classmethod
    def deletion_distort(cls, network_set, prob):
        assert 0 <= prob <= 1, "Probability must be between 0 and 1"
        result = cls(network_size=network_set.network_size)
        for network_info in network_set.per_network_info():
            if random.uniform(0, 1) <= prob:
                pass
            else:
                result.append(copy.deepcopy(network_info))
        return result

    @classmethod
    def tail_move_distort(cls, network_set, prob):
        assert 0 <= prob <= 1, "Probability must be between 0 and 1"
        result = cls(network_size=network_set.network_size)
        for network_info in network_set.per_network_info():
            network = copy.deepcopy(network_info.network)
            if random.uniform(0, 1) <= prob:
                moving_arc = random.sample(network.movable_arcs, 1)[0]
                nodes_below_moving_arc = network.nodes_below_nodes({moving_arc[0]})
                all_arcs = network.arcs
                all_arcs.remove(moving_arc)
                all_arcs = [arc for arc in all_arcs if
                            arc[1] != moving_arc[1] and
                            arc[0] not in nodes_below_moving_arc]
                if network.is_tree_node(moving_arc[0]):
                    parent_moving_arc = network.nodes_above_nodes({moving_arc[0]}, max_height=1).difference({moving_arc[0]}).pop()
                    all_arcs.remove((parent_moving_arc, moving_arc[0]))
                if len(all_arcs) == 0:
                    raise NotImplementedError("Could not find a target arc")
                target_arc = random.sample(all_arcs, 1)[0]
                network.tail_move(moving_arc, target_arc)

                if len(network.roots) != 1:
                    raise NotImplementedError("Tail move results in more than 1 root")
            result.append(NetworkInfo(network, multiplicity=network_info.multiplicity))

        return result

    @classmethod
    def from_enewick_format(cls, file_path, network_size=3, sep=';\n'):
        result = cls(network_size)
        with open(file_path, 'r') as f:
            data = f.read()
        data = data.split(sep)
        for d in data:
            try:
                network = RootedLevelKNetwork.from_enewick(d)
            except:
                print(d)
                raise NotImplementedError
            result.append(NetworkInfo(network))
        return result

    @classmethod
    def from_named_trinet_format(cls, file_path, progress_bar=False):
        result = cls(3)
        with open(file_path, 'r') as f:
            data = f.read()
        data = data.split('\n')
        if progress_bar:
            for d in tqdm(data):
                index = d.index('Tr{')
                d = d[index:]
                network = RootedLevelKNetwork.from_trinet_format(d)
                result.append(NetworkInfo(network))
        else:
            for d in data:
                index = d.index('Tr{')
                d = d[index:]
                network = RootedLevelKNetwork.from_trinet_format(d)
                result.append(NetworkInfo(network))
        return result

    def summary(self, depth=None):
        depth = coalesce(depth, [0, 1])
        # Compute properties of this network set
        if 0 in depth:
            summary_0 = self._analyse_networks()
        else:
            summary_0 = dict()
        if 1 in depth:
            # Compute properties of network set of one size smaller
            network_info_list_1 = NetworkSet.induced_network_set_of_network_set(self, self.network_size - 1, progress_bar=False)
            summary_1 = network_info_list_1._analyse_networks()

            summary_0['inconsistency-1'] = summary_1['inconsistency-0']
            summary_0['certainty-1'] = summary_1['certainty-0']

        return summary_0

    def _analyse_networks(self, method: int = settings.WEIGHTED_AVERAGE) -> (dict, float, float, float):
        assert self.network_size is not None, "Can not create summary of network set without a given size"
        # Number of different n-sets
        represented_leaves = self.represented_leaves()
        number_of_represented_leaves = len(represented_leaves)
        number_of_n_sets = ncr(len(represented_leaves), self.network_size)

        # Number of n-sets described by self
        described_n_sets = len(self)

        # Total number of networks (including multiplicity)
        total_networks = self.volume()

        # Number of unique networks (multiplicity = 1)
        unique_networks = len([1 for _ in self.per_network_info()])

        # Properties
        duplicity = (total_networks - unique_networks) / described_n_sets
        density = described_n_sets / number_of_n_sets
        inconsistency = (unique_networks - described_n_sets) / described_n_sets
        certainty = 0
        n_set_iterator = itertools.combinations(represented_leaves, self.network_size)
        for n_set in tqdm(n_set_iterator):
            name = tuple(sorted(n_set))
            n_set_certainty = 1
            empty = False
            for _, weight in self.per_network_iterator(method=method, name=name):
                n_set_certainty *= weight
                empty = True
            n_set_certainty *= (1 - empty)
            certainty += n_set_certainty
        certainty /= number_of_n_sets

        summary = {
            'duplicity-0'      : duplicity
            , 'density-0'      : density
            , 'inconsistency-0': inconsistency
            , 'certainty-0'    : certainty
            , 'n'              : number_of_represented_leaves
        }
        return summary

    def is_consistent(self) -> bool:
        pass
        # _, redundancy, density, inconsistency_0 = self._analyse_networks()
        # network_info_list_1 = NetworkSet.induced_network_set_of_network_set(self, self.network_size - 1)
        # covered_networks_1, _, _, inconsistency_1 = network_info_list_1._analyse_networks()
        # return bool(redundancy == 0.0 and density == 1.0 and inconsistency_0 == 0.0 and inconsistency_1 == 0.0)

    def calculate_info(self):
        for network_info in self.per_network_info():
            network_info.calculate_info()

    def add_structure_info(self, network_info_list):
        for network_info in self.per_network_info():
            equal_structured_network, relation_dict = network_info_list.find_equal_structured_network(network_info)
            if equal_structured_network:
                network_info.add_structure_info(equal_structured_network, relation_dict)
            else:
                network_info.network.visualize()
                raise NotImplementedError

    # ------------------ CONSISTENCY SCORE -------------- #
    def consistency_score(self, other, method):
        intersection = 0
        network_set_len = 0
        for network_info, w in self.per_network_iterator(method=method):
            network_set_len += w
            name = network_info.name()
            p = 0
            for other_network_info, other_w in other.per_network_iterator(method=method, name=name):
                equal, _, _ = network_info.network.equal_structure(other_network_info.network, optimize=True, equal_naming=True)
                if equal:
                    p = min(p + other_w, w)
                    if p == w:
                        break
            intersection += p
        return intersection / network_set_len

    # ------------------ COUNT METHODS ------------------ #
    def per_name_iterator(self, leaf_names, method):
        assert self.network_size, "Can only iterate per name if network_size is defined"
        name_iterator = itertools.combinations(leaf_names, self.network_size)
        for name in name_iterator:
            name = tuple(sorted(name))
            yield from self.per_network_iterator(method=method, name=name, empty_names=True)

    def per_network_iterator(self, method, name=None, empty_names: bool = False):
        if method == settings.MAXIMUM_MULTIPLICITY:
            return self.maximum_multiplicity_iterator(name, empty_names)
        elif method == settings.WEIGHTED_AVERAGE:
            return self.weighted_average_iterator(name, empty_names)
        elif method == settings.WEIGHTED_SUM:
            return self.weighted_sum_iterator(name, empty_names)
        else:
            raise ValueError

    def maximum_multiplicity_iterator(self, name=None, empty_names=False):
        if name:
            if name in self.dictionary.keys():
                dictionary = self[name]
                maximum_multiplicity = -np.inf
                corresponding_network_info = None
                for network_info in dictionary['network_info_set']:
                    if network_info.multiplicity > maximum_multiplicity:
                        maximum_multiplicity = network_info.multiplicity
                        corresponding_network_info = network_info
                yield corresponding_network_info, 1
            elif empty_names:
                yield name, None
        else:
            for name, dictionary in self:
                maximum_multiplicity = -np.inf
                corresponding_network_info = None
                for network_info in dictionary['network_info_set']:
                    if network_info.multiplicity > maximum_multiplicity:
                        maximum_multiplicity = network_info.multiplicity
                        corresponding_network_info = network_info
                yield corresponding_network_info, 1

    def weighted_average_iterator(self, name=None, empty_names: bool = False):
        if name:
            if name in self.dictionary.keys():
                dictionary = self[name]
                volume = dictionary['volume']
                for network_info in dictionary['network_info_set']:
                    yield network_info, network_info.multiplicity / volume
            elif empty_names:
                yield name, None
        else:
            for name, dictionary in self:
                volume = dictionary['volume']
                for network_info in dictionary['network_info_set']:
                    yield network_info, network_info.multiplicity / volume

    def weighted_sum_iterator(self, name=None, empty_names: bool = False):
        if name:
            if name in self.dictionary.keys():
                dictionary = self[name]
                for network_info in dictionary['network_info_set']:
                    yield network_info, network_info.multiplicity
            elif empty_names:
                yield name, None
        else:
            for name, dictionary in self:
                for network_info in dictionary['network_info_set']:
                    yield network_info, network_info.multiplicity

    # ------------------------------- COUNTING METHODS ----------------------------------- #
    def count_category(self, category: str, method: int):
        result = {}
        for network_info, w in self.per_network_iterator(method=method):
            key = network_info[category]
            if type(key) == dict:
                for key_2, value in key.items():
                    if key_2 in result.keys():
                        if value in result[value].keys():
                            result[key_2][value] += w
                        else:
                            result[key_2][value] = w
                    else:
                        result[key_2] = {value: w}
            else:
                if key in result.keys():
                    result[key] += w
                else:
                    result[key] = w
        return result

    # ------------------ ADD METHODS -------------------- #
    def append(self, other: NetworkInfo):
        assert type(other) == NetworkInfo, "Can only add object of type NetworkInfo to a NetworkInfoList"
        assert self.network_size is None or other.network.number_of_leaves == self.network_size, \
            f"Can only add network equal to {self.network_size} leaves to this NetworkInfoList, this network has {len(other.network.leaf_names)} leaves"
        name = other.name()
        if name in self.dictionary.keys():
            network_info_set = self.dictionary[name]['network_info_set']
            self.dictionary[name]['volume'] += other.multiplicity
            for network_info in network_info_set:
                equal, _, _ = network_info.network.equal_structure(other.network, equal_naming=self.equal_naming, optimize=False)
                if equal:
                    network_info.multiplicity += other.multiplicity
                    break
            else:
                network_info_set.add(other)
        else:
            self.dictionary[name] = {'volume': other.multiplicity, 'network_info_set': {other}}

    def extend(self, other):
        assert type(other) == NetworkSet, "Can only extend using a NetworkSet"
        for name, dictionary in other:
            for network_info in dictionary['network_info_set']:
                self.append(network_info)

    def __add__(self, other):
        result = copy.deepcopy(self)
        result.extend(other)
        return result

    def __setitem__(self, name, value):
        self.dictionary[name] = value

    def __getitem__(self, name):
        return self.dictionary[name]

    def any(self):
        return next(self.per_network_info())

    # ----------------- MISC -------------------------#
    def set_multiplicities_to_one(self):
        for name in self.dictionary.keys():
            network_info_set = self.dictionary[name]['network_info_set']
            self.dictionary[name]['volume'] += len(network_info_set)
            for network_info in network_info_set:
                network_info.multiplicity = 1

    # ----------------- FIND NETWORK METHODS ------------------- #
    def find_equal_structured_network(self, network_info: NetworkInfo):
        for other_network_info in self.per_network_info():
            are_equal, relation_dict, leaf_translation_dicts = other_network_info.network.equal_structure(network_info.network, optimize=True)
            if are_equal:
                return other_network_info, relation_dict[0]
        return None, None

    def find_equal_network(self, network_info: NetworkInfo):
        name = network_info.name()
        if name in self:
            dictionary = self[name]
            for other_network_info in dictionary['network_info_set']:
                if other_network_info.network == network_info.network:
                    return other_network_info
        return False

    @classmethod
    def restrict(cls, network_set, leaf_names):
        result = cls(network_size=network_set.network_size)
        for name, dictionary in network_set:
            overlap = set(leaf_names).intersection(list(name))
            if len(overlap) == network_set.network_size:
                result[name] = copy.deepcopy(dictionary)
        return result

    @classmethod
    def networks_where(cls, network_set, category, value, relation=lambda x, y: x == y):
        result = cls(network_size=network_set.network_size)
        for name, dictionary in network_set:
            new_dictionary = {'volume': 0, 'network_info_set': set()}
            for network_info in dictionary['network_info_set']:
                if relation(network_info[category], value):
                    new_dictionary['network_info_set'].add(network_info)
                    new_dictionary['volume'] += network_info.multiplicity
            if new_dictionary['volume'] != 0:
                result[name] = new_dictionary
        return result

    @classmethod
    def networks_with_category(cls, network_set, category):
        result = cls(network_size=network_set.network_size)
        for name, dictionary in network_set:
            new_dictionary = {'volume': 0, 'network_info_set': set()}
            for network_info in dictionary['network_info_set']:
                if category in network_info.info.keys():
                    new_dictionary['network_info_set'].add(network_info)
                    new_dictionary['volume'] += network_info.multiplicity
            if new_dictionary['volume'] != 0:
                result[name] = new_dictionary
        return result

    @classmethod
    def networks_of_size(cls, network_set, size):
        result = cls(network_size=size)
        for name, dictionary in network_set:
            if len(name) == size:
                result[name] = copy.deepcopy(dictionary)
        return result

    def biconnected_networks(self):
        result = NetworkSet(network_size=self.network_size)
        for network_info in self:
            if network_info.network.biconnected:
                result.append(network_info)
        return result

    def represented_leaves(self) -> set:
        result = set()
        for name, _ in self:
            result.update(name)
        return result

    def to_enewick(self):
        return {network_info.network.enewick(): network_info.multiplicity for network_info in self.per_network_info()}

    def save_to_file(self, file_path, frmt=settings.FORMAT_eNewick_multiplicity):
        if frmt == settings.FORMAT_eNewick_multiplicity:
            data = self.to_enewick()
            text_data = '\n'.join([f"{enewick} {multiplicity}" for enewick, multiplicity in data.items()])
        elif frmt == settings.FORMAT_eNewick:
            data = self.to_enewick()
            text_data = '\n'.join([f"{enewick}" for enewick, _ in data.items()])
        # elif frmt == settings.FORMAT_tnet:
        #     data = []
        #     for i, network_info in enumerate(self.per_network_info()):
        #         level_1_trinet_name, leaf_names = network_info.network.level_1_trinet_name(leaf_name_sep=', ')
        #         if level_1_trinet_name is None:
        #             raise ValueError("Trinet is not level-1")
        #         data.append(f"Tr{i} = Tr{{{leaf_names} : {level_1_trinet_name}}}")
        #     text_data = '\n'.join(data)
        elif frmt == settings.FORMAT_tnets:
            data = []
            for network_info in self.per_network_info():
                level_1_trinet_name, leaf_names = network_info.network.level_1_trinet_name(leaf_name_sep=' ')
                if level_1_trinet_name is None:
                    raise ValueError("Trinet is not level-1")
                data.append(f"{leaf_names} {level_1_trinet_name}")
            text_data = '\n'.join(data)
        else:
            raise ValueError('Format not defined')
        with open(f'{file_path}.{settings.FORMAT_extension[frmt]}', 'w+') as f:
            f.write(text_data)

    def per_network_info(self):
        for _, dictionary in self:
            for network_info in dictionary['network_info_set']:
                yield network_info

    def __iter__(self):
        return iter(self.dictionary.items())

    def __len__(self):
        return len([[1 for _ in name] for name in self])

    # TODO naming
    def volume(self):
        return sum([network_info.multiplicity for network_info in self.per_network_info()])

    def __copy__(self):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.dictionary = copy.copy(self.dictionary)
        cp.network_size = copy.copy(self.network_size)
        cp.uid = guid()
        cp.logger = logging.getLogger('network_set.{}'.format(cp.uid))
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.dictionary = copy.deepcopy(self.dictionary)
        cp.network_size = copy.deepcopy(self.network_size)
        cp.uid = guid()
        cp.logger = logging.getLogger('network_set.{}'.format(cp.uid))
        return cp

    def __getstate__(self):
        self.logger = 'network_set.{}'.format(self.uid)
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
    def from_network_info_list(cls, network_info_list: NetworkSet, count_method, fill_gaps):
        result = cls(node_set=network_info_list.represented_leaves())
        # TODO: Test
        for network_info, weight in network_info_list.per_name_iterator(leaf_names=network_info_list.represented_leaves(), method=count_method):
            if weight:
                cut_arc_sets = network_info['cut_arc_sets']
                for cut_arc_set in cut_arc_sets:
                    if len(cut_arc_set) == 2:
                        x = cut_arc_set[0]
                        y = cut_arc_set[1]
                        z = [Z for Z in network_info.network.leaf_names if Z not in cut_arc_set][0]
                        for i in (x, y):
                            result.increase_weight(i, z, weight)
            elif fill_gaps:
                leaves = network_info
                two_leaf_iterator = itertools.permutations(leaves, 2)
                for two_leaves in two_leaf_iterator:
                    i, z = two_leaves
                    result.increase_weight(i, z, fill_gaps / 3)
        return result

    def add_node(self, node):
        assert node not in self.node_set, f"Node {node} is already in this graph"
        for other_node in self.node_set:
            self.arc_weights[(node, other_node)] = 0
            self.arc_weights[(other_node, node)] = 0
        self.node_set.add(node)

    def increase_weight(self, from_node, to_node, w: float = 1):
        assert {from_node, to_node}.issubset(self.node_set), f"Can not add weight to arc {(from_node, to_node)} as it is not in the graph"
        self.arc_weights[(from_node, to_node)] += w

    def augmented_graph(self, threshold=0):
        auxiliary_digraph = {node: set() for node in self.node_set}
        arc_iterator = itertools.permutations(self.node_set, 2)
        for arc in arc_iterator:
            weight = self.arc_weights[arc]
            if weight <= threshold:
                auxiliary_digraph[arc[0]].add(arc[1])
        return auxiliary_digraph

    def minimum_weight(self):
        return min(self.arc_weights.values())

    def minimal_sink_sets(self, method: int):
        # TODO check complexity and maybe do something smarter
        if method == settings.FIRST_SCC_THAT_IS_MSS:
            return self.minimal_sink_sets_default()
        elif method == settings.EXPAND_FIRST_SCC:
            return self.minimal_sink_sets_expand()
        else:
            raise ValueError

    def minimal_sink_sets_default(self):
        weights = self.arc_weights.values()
        for threshold in sorted(weights):
            sub_graph = self.augmented_graph(threshold)
            strongly_connected_components = tarjan(sub_graph)
            minimal_sink_set = self.strongly_connected_components_to_minimal_sink_set_default(strongly_connected_components, sub_graph)
            if minimal_sink_set:
                return minimal_sink_set
        self.visualize(threshold=100)
        raise NotImplementedError("No threshold gives a minimal sink-set")

    @staticmethod
    def strongly_connected_components_to_minimal_sink_set_default(strongly_connected_components, sub_graph):
        # Find minimal sink-sets
        minimal_sink_sets = []
        for strongly_connected_component in strongly_connected_components:
            for node in strongly_connected_component:
                if node in sub_graph.keys():
                    to_leaves = sub_graph[node]
                    if not set(to_leaves).issubset(strongly_connected_component):
                        break
            else:
                if len(strongly_connected_component) > 1:
                    minimal_sink_sets.append(strongly_connected_component)
        if len(minimal_sink_sets) > 0:
            smallest_mss_index = int(np.argmin([len(mss) for mss in minimal_sink_sets]))
            return minimal_sink_sets[smallest_mss_index]
        else:
            return None

    def minimal_sink_sets_expand(self):
        # Maybe maximum cut in undirected network where weight is sum of directed arc weights
        threshold = self.minimum_weight()
        sub_graph = self.augmented_graph(threshold)
        strongly_connected_components = tarjan(sub_graph)
        minimal_sink_set = self.strongly_connected_components_to_minimal_sink_set_expand(strongly_connected_components, sub_graph)
        return minimal_sink_set

    @staticmethod
    def strongly_connected_components_to_minimal_sink_set_expand(strongly_connected_components, sub_graph):
        # Create condensed graph
        condensed_graph = {}
        for i, SCC in enumerate(strongly_connected_components):
            condensed_graph[i] = set()
            for node in SCC:
                to_nodes = sub_graph[node]
                for to_node in to_nodes:
                    j = next(i for i, SCC in enumerate(strongly_connected_components) if to_node in SCC)
                    condensed_graph[i].add(j)

        # Find: SCC of size at least two with out-degree 0 | SCC with out-degree > 1
        parent_sink_sets = []
        minimal_sink_sets = []
        for i, SCC in enumerate(strongly_connected_components):
            if len(condensed_graph[i]) == 0 and len(SCC) > 1:
                minimal_sink_sets.append(SCC)
            elif len(condensed_graph[i]) >= 1:
                parent_sink_sets.append(i)

        # Return smallest minimal sink-set if one exists
        if len(minimal_sink_sets) != 0:
            smallest_mss_index = int(np.argmin([len(mss) for mss in minimal_sink_sets]))
            return minimal_sink_sets[smallest_mss_index]

        # Else extend an internal node of the condensed graph with its descendants
        else:
            for i in parent_sink_sets:
                MSS = {i}
                new_nodes = condensed_graph[i]
                while not new_nodes.issubset(MSS):
                    MSS.update(new_nodes)
                    new_nodes = set(itertools.chain.from_iterable(condensed_graph[SCC] for SCC in new_nodes))
                minimal_sink_sets.append(list(itertools.chain.from_iterable([strongly_connected_components[node] for node in MSS])))

            # Return smallest sink-set
            smallest_mss_index = int(np.argmin([len(mss) for mss in minimal_sink_sets]))
            return minimal_sink_sets[smallest_mss_index]

    def visualize(self, file_path: str = None, arc_labels=True, threshold: int = 0, format='png'):
        dot = Digraph()
        dot.engine = 'dot'

        for node_name in self.node_set:
            dot.node(name=node_name, label=node_name, **{'width': str(0.5), 'height': str(0.5)})
        for arc, weight in self.arc_weights.items():
            if weight > threshold:
                continue
            if arc_labels:
                dot.edge(arc[0], arc[1], label=str(weight))
            else:
                dot.edge(arc[0], arc[1])
        if file_path:
            dot.format = format
            if file_path[-4:] == f'.{format}':
                dot.render(filename=file_path[:-4], format=format)
            else:
                dot.render(filename=file_path, format=format)
        else:
            dot.render(view=True)
            time.sleep(0.2)
