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
        self._biconnected_components = None
        self._partial_ordering = None
        self._cut_arc_matrix = None
        self._cut_arc_sets = None

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
        cp._biconnected_components = copy.copy(self._biconnected_components)
        cp._partial_ordering = copy.copy(self._partial_ordering)
        cp._cut_arc_matrix = copy.copy(self._cut_arc_matrix)
        cp._cut_arc_sets = copy.copy(self._cut_arc_sets)
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
        cp._biconnected_components = copy.deepcopy(self._biconnected_components)
        cp._partial_ordering = copy.deepcopy(self._partial_ordering)
        cp._cut_arc_matrix = copy.copy(self._cut_arc_matrix)
        cp._cut_arc_sets = copy.deepcopy(self._cut_arc_sets)
        return cp

    @classmethod
    def from_enewick(cls, string, check_valid=True):
        adjacency_dict = enewick(string)
        result = cls.from_connections_dict(adjacency_dict, check_valid=check_valid)
        return result

    @classmethod
    def trinet_from_network(cls, network, node_names: list = None, node_numbers: list = None):
        """Create trinet from network."""
        assert (node_numbers is not None) or (node_names is not None), "At least one of node_names and node_numbers has to be given"
        if node_numbers is not None:
            taxa = node_numbers
        else:
            taxa = network.get_node_numbers_of(node_names)
        assert len(taxa) <= 3, "Can not create trinet from network {} using as this {} are more than 3 leaves.".format(network.uid, taxa)
        trinet = cls.from_network(network, node_numbers=taxa, suppress_redundant='all', suppress_parallel=True)
        return trinet

    @classmethod
    def from_network(cls, network, node_numbers: list = None, node_names: list = None, suppress_redundant='none', suppress_parallel=False):
        """Create sub-network from network."""
        assert (node_numbers is not None) or (node_names is not None), "At least one of node_names and node_numbers has to be given"
        if node_numbers is not None:
            taxa = node_numbers
        else:
            taxa = network.get_node_numbers_of(node_names)
        logging.debug("Creating sub=network from {} using {}.".format(network.uid, taxa))
        network.logger.debug("Creating sub-network using {}.".format(taxa))
        assert set(taxa).issubset(set(network.leaf_numbers)), "Can not create sub-network of network {} using {} as they are not leaves of network.".format(
            network.uid,
            taxa)
        adj_matrix = copy.deepcopy(network.adj_matrix)
        node_name_map = copy.deepcopy(network.node_name_map)
        level = copy.copy(network.level)
        dimension = copy.copy(network.dimension)
        taxa = copy.copy(taxa)
        new_network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_numbers=taxa, level=level, dimension=dimension)
        new_network.prune(suppress_redundant=suppress_redundant, suppress_parallel=suppress_parallel)
        new_network.logger.debug("Sub-network of network {} with taxa {}.".format(network.uid, taxa))
        return new_network

    @classmethod
    def from_dir_adj_matrix(cls, dir_adj_matrix: np.ndarray, level=2, dimension=2, check_valid=True, char_type='alph'):
        """Create network from directed adjacency matrix. Assumes only internal nodes have rows. Other nodes are leaves."""
        logging.debug("Creating network from directed adjacency matrix.")
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
        network.logger.debug("Created from directed adjacency matrix.")
        return network

    @classmethod
    def from_connections_dict(cls, connections_dict: dict, level=2, dimension=2, check_valid=True):
        """Create network from connection dictionary. Assumes nodes without outgoing arcs are leaves."""
        node_name_map = bidict()
        leaf_names = list()

        # Empty adjacency matrix
        adj_matrix = np.zeros((0, 0)).astype(int)
        network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_names=leaf_names, level=level, dimension=dimension)

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
                network._add_connection(from_node_number, to_node_number)
        assert (not check_valid) or network.is_valid(), "Connection dictionary results in an invalid network."
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
        biconnected_components = self._get_biconnected_components(leafless=True)
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

    # --------------------------------------------------------- ? --------------------------------------------------------- #
    @property
    def node_names(self) -> list:
        """Retrieve list of node names ordered by their number"""
        return self.get_node_names_of(sorted(list(self.node_name_map.values())))

    @property
    def leaf_names(self) -> list:
        """Retrieve list of leaf names ordered by their number"""
        return self.get_node_names_of(sorted(self.leaf_numbers))

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

    def is_connected_node(self, node_name: str) -> bool:
        """Check whether node is connected to the graph."""
        return self._is_connected_node(self.get_node_number_of(node_name))

    def _is_connected_node(self, node_number: int) -> bool:
        """Check whether node is connected to the graph."""
        return sum(self.adj_matrix[node_number, :] != 0) > 0

    # --------------------------------------------------------- ALTER NETWORK METHODS --------------------------------------------------------- #
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

    def rename_node(self, old_name: str, new_name: str) -> None:
        """Replace name of node with name old_name by new_name."""
        assert new_name not in self.node_names, f"There already exists a node with node name {new_name}."

        # Rename node in node_names
        node_number = self.node_name_map.pop(old_name)
        self.node_name_map.put(new_name, node_number)

        # Rename node in leaf_names if node is a leaf
        try:
            self.leaf_numbers.remove(old_name)
            self.leaf_numbers.append(new_name)
        except ValueError:
            pass

        self._biconnected_components = None
        self._partial_ordering = None
        self._cut_arc_sets = None

    def _rename_node(self, node_number: int, new_name: str) -> None:
        assert new_name not in self.node_names, f"There already exists a node with node name {new_name}."

        # Rename node in node_names
        old_name = self.node_name_map.inverse[node_number]
        self.node_name_map.pop(old_name)
        self.node_name_map.put(new_name, node_number)

        # Rename node in leaf_names if node is a leaf
        try:
            self.leaf_numbers.remove(old_name)
            self.leaf_numbers.append(new_name)
        except ValueError:
            pass

        self._biconnected_components = None
        self._partial_ordering = None
        self._cut_arc_sets = None

    def shrink(self, leaf_names_set):
        """Shrink set of leaves"""
        intersection = set(self.leaf_names).intersection(leaf_names_set)
        if len(intersection) == len(self.leaf_names):
            return False
        if len(intersection) == 0:
            return True
        leaf_numbers_set = set(self.get_node_numbers_of(list(intersection)))
        mss_name = mss_leaf_name(leaf_names_set)
        return self._shrink(leaf_numbers_set, mss_name)

    def _shrink(self, leaf_numbers_set, mss_name):
        leaf_to_keep = leaf_numbers_set.pop()
        for leaf in leaf_numbers_set:
            self._remove_leaf_status_node(leaf)
        self._rename_node(leaf_to_keep, mss_name)
        self.prune()
        return mss_name

    def _add_connection(self, from_node_number: int, to_node_number: int):
        """Add connection between from_node_name and to_node_name."""
        self.adj_matrix[from_node_number][to_node_number] += 1
        self.adj_matrix[to_node_number][from_node_number] -= 1

        self._biconnected_components = None
        self._partial_ordering = None
        self._cut_arc_matrix = None
        self._cut_arc_sets = None

    def _remove_connection(self, from_node_number: int, to_node_number: int):
        """Remove connection between from_node_name and to_node_name."""
        assert self.adj_matrix[from_node_number][to_node_number] > 0 > self.adj_matrix[to_node_number][
            from_node_number], "Cannot remove non-existent connection"

        self.adj_matrix[from_node_number][to_node_number] -= 1
        self.adj_matrix[to_node_number][from_node_number] += 1

        self._biconnected_components = None
        self._partial_ordering = None
        self._cut_arc_matrix = None
        self._cut_arc_sets = None

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

    def _remove_node_from_network(self, node_number: int, force: bool = False):
        """Remove node fully from network."""
        if self._is_leaf_node(node_number):
            self._remove_leaf_status_node(node_number)

        self._remove_node_from_adj_matrix(node_number)
        self._remove_node_from_dict(node_number)
        self.number_of_nodes -= 1

        self.reset_optimization_variables()

    def _remove_component(self, component: list) -> None:
        """Remove all nodes in component."""
        for node_number in sorted(component, reverse=True):
            self._remove_node_from_network(node_number=node_number, force=True)

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
        return node_name, node_number

    def _suppress_node(self, node_number: int):
        """Remove node and connect its parent to its child if they both exist."""
        self._suppress_component([node_number])

    def _suppress_component(self, component: list):
        """Remove all nodes in component and connect its parent to its child if they both exist."""
        in_nodes = self._get_in_nodes_component(component)
        out_nodes = self._get_out_nodes_component(component)
        self._remove_component(component)

        if len(in_nodes) == 1 and len(out_nodes) == 1:
            out_node = out_nodes[0]
            in_node = in_nodes[0]
            out_node -= sum(np.array(component) < out_node)
            in_node -= sum(np.array(component) < in_node)
            self._add_connection(in_node, out_node)

    def remove_leaf(self, leaf_name: str):
        """Remove leaf and prune to make a proper network again"""
        leaf_number = self.get_node_number_of(leaf_name)
        self._remove_leaf_status_node(leaf_number)
        self.prune()

    def _remove_leaf_status_node(self, node_number: int):
        """Remove leaf with name node_name"""
        self.leaf_numbers.remove(node_number)

    def _add_leaf_status_node(self, node_number: int):
        """Remove leaf with name node_name"""
        self.leaf_numbers.append(node_number)

    def add_leaf_to_edge(self, edge, leaf_name: str = None, char_type='alph') -> (str, str):
        """Add node between parent and child and attach leaf with leaf_name to it."""
        return self.get_node_names_of(list(self._add_leaf_to_edge(self.get_node_numbers_of(edge), leaf_name, char_type)))

    def _add_leaf_to_edge(self, edge, leaf_name: str = None, char_type='alph') -> (int, int):
        parent = edge[0]
        child = edge[1]

        internal_name, internal_number = self._add_node_to_network()
        leaf_name, leaf_number = self._add_node_to_network(leaf_name, leaf=True, char_type=char_type)

        self._remove_connection(parent, child)

        self._add_connection(parent, internal_number)
        self._add_connection(internal_number, child)
        self._add_connection(internal_number, leaf_number)
        return internal_number, leaf_number

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

        self._add_connection(leaf_number_1, new_leaf_number_1)
        self._add_connection(leaf_number_2, new_leaf_number_2)
        self._add_connection(leaf_number_1, internal_leaf_number)
        self._add_connection(leaf_number_2, internal_leaf_number)
        self._add_connection(internal_leaf_number, recombined_leaf_number)

    def split_leaf(self, leaf_name_to_split: str):
        """ Split leaf into two leaves"""
        assert self.is_leaf_node(leaf_name_to_split), "Can not split non-leaf-node"
        self._split_leaf(self.get_node_number_of(leaf_name_to_split))

    def _split_leaf(self, leaf_number_to_split: int):
        self._remove_leaf_status_node(leaf_number_to_split)
        leaf_name_to_split = self.get_node_name_of(leaf_number_to_split)
        _, left_leaf_number = self._add_node_to_network(node_name='(' + leaf_name_to_split + "-I)", leaf=True)
        _, right_leaf_number = self._add_node_to_network(node_name='(' + leaf_name_to_split + "-II)", leaf=True)
        self._add_connection(leaf_number_to_split, left_leaf_number)
        self._add_connection(leaf_number_to_split, right_leaf_number)

    def replace_leaf_with_network(self, leaf_name_to_replace: str, replacement_network, replace_names: bool = False, char_type='alph'):
        assert self.is_leaf_node(leaf_name_to_replace), f"Cannot replace leaf {leaf_name_to_replace} with network as it is not a leaf."
        assert replace_names or set(self.leaf_names).isdisjoint(
            replacement_network.leaf_names), f"Cannot replace leaf {leaf_name_to_replace} with network \n \n {replacement_network}  \n \n as it has some leafs same as \n \n {self}"
        leaf_number_to_replace = self.get_node_number_of(leaf_name_to_replace)
        return self._replace_leaf_with_network(leaf_number_to_replace, replacement_network, replace_names, char_type)

    def _replace_leaf_with_network(self, leaf_number_to_replace: int, replacement_network, replace_names: bool = False, char_type='alph'):
        replacement_network = copy.deepcopy(replacement_network)
        replacement_network_internal_node_names = set(replacement_network.internal_node_names)
        replacement_network_root_name = replacement_network.get_root()

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
        parent_of_leaf = self._get_in_nodes_node(leaf_number_to_replace)[0]
        self._remove_node_from_network(leaf_number_to_replace)
        self._add_connection(parent_of_leaf, self.get_node_number_of(replacement_network_root_name + "*"))

        # Add all connections from replacement network to current network
        for internal_node_name in replacement_network_internal_node_names:
            to_node_names = replacement_network.get_out_nodes_node(internal_node_name)
            internal_node_number = self.get_node_number_of(internal_node_name + "*")
            for to_node_name in to_node_names:
                if to_node_name in replacement_network_internal_node_names or replace_names:
                    to_node_name += "*"
                self._add_connection(internal_node_number, self.get_node_number_of(to_node_name))

        if replace_names:
            self.standardize_node_names(char_type)

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

    def _get_root(self) -> int:
        return np.where(sum(self.adj_matrix > 0) == 0)[0][0]

    def get_root(self) -> str:
        """Retrieve the name of the root."""
        return self.get_node_name_of(self._get_root())

    def get_reticulations(self) -> list:
        """Retrieve node names of all reticulations."""
        return self.get_node_names_of(self._get_reticulations())

    def _get_reticulations(self) -> list:
        return list(np.where(np.array(self.get_in_degrees()) > 1)[0])

    # ---------------------------------------------------------------- DEGREE METHODS ---------------------------------------------------------------- #
    def get_in_degree_node(self, node_name: str, leafless: bool = False) -> int:
        """Retrieve number of arcs entering node_name."""
        return self._get_in_degree_node(self.get_node_number_of(node_name))

    def _get_in_degree_node(self, node_number: int, leafless: bool = False) -> int:
        """Retrieve number of arcs entering node_name."""
        mask = self.leaf_mask if leafless else np.ones(self.number_of_nodes)
        return sum(np.multiply(np.multiply(self.adj_matrix[node_number], self.adj_matrix[node_number] > 0), mask))

    def get_out_degree_node(self, node_name: str, leafless: bool = False) -> int:
        """Retrieve number of arcs exiting node_name."""
        return self._get_out_degree_node(self.get_node_number_of(node_name))

    def _get_out_degree_node(self, node_number: int, leafless: bool = False) -> int:
        mask = self.leaf_mask if leafless else np.ones(self.number_of_nodes)
        return sum(np.multiply(np.multiply(self.adj_matrix.T[node_number], self.adj_matrix.T[node_number] < 0), mask))

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

    def get_in_out_degree_node(self, node_name: str, leafless: bool = False):
        return self._get_in_out_degree_node(self.get_node_number_of(node_name), leafless)

    def _get_in_out_degree_node(self, node_number: int, leafless: bool = False):
        return self._get_in_degree_node(node_number, leafless), self._get_out_degree_node(node_number, leafless)

    # ---------------------------------------------------------------- LEVEL METHODS ---------------------------------------------------------------- #
    def strict_level(self) -> int:
        in_degrees = np.array(self.get_in_degrees())

        level = 0
        biconnected_components = self._get_biconnected_components(leafless=True)
        for bc in biconnected_components:
            in_sum_bc = in_degrees[bc]
            level = max(level, sum(in_sum_bc > 1))

        return level

    # ---------------------------------------------------------------- CONNECTION METHODS ---------------------------------------------------------------- #
    def get_in_nodes_node(self, node_name: str) -> list:
        """Retrieve all nodes which enter node node_name."""
        return self.get_node_names_of(self._get_in_nodes_node(self.get_node_number_of(node_name)))

    def _get_in_nodes_node(self, node_number: int) -> list:
        return list(np.where(self.adj_matrix[node_number, :] == -1)[0]) + list(np.where(self.adj_matrix[node_number, :] == -2)[0]) * 2

    def get_out_nodes_node(self, node_name: str) -> list:
        """Retrieve all nodes which exit node node_name."""
        return self.get_node_names_of(self._get_out_nodes_node(self.get_node_number_of(node_name)))

    def _get_out_nodes_node(self, node_number: int) -> list:
        return list(np.where(self.adj_matrix[node_number, :] == 1)[0]) + list(np.where(self.adj_matrix[node_number, :] == 2)[0]) * 2

    def get_in_nodes_component(self, component: list) -> list:
        """Retrieve all nodes which enter component."""
        return self.get_node_names_of(self._get_in_nodes_component(self.get_node_numbers_of(component)))

    def _get_in_nodes_component(self, component: list) -> list:
        in_nodes = np.array([])
        for node_number in component:
            in_nodes = np.array(list(itertools.chain(in_nodes, self._get_in_nodes_node(node_number))))
        return [node for node in in_nodes if node not in component]

    def get_out_nodes_component(self, component: list) -> list:
        """Retrieve all nodes which exit component."""
        return self.get_node_names_of(self._get_out_nodes_component(self.get_node_numbers_of(component)))

    def _get_out_nodes_component(self, component: list) -> list:
        in_nodes = np.array([])
        for node_number in component:
            in_nodes = np.array(list(itertools.chain(in_nodes, self._get_out_nodes_node(node_number))))
        return [node for node in in_nodes if node not in component]

    def get_connections_node(self, node_name: str) -> list:
        """Retrieve all nodes which enter and exit node."""
        return self.get_node_names_of(self._get_connections_node(self.get_node_number_of(node_name)))

    def _get_connections_node(self, node_number: int) -> list:
        return self._get_in_nodes_node(node_number) + self._get_out_nodes_node(node_number)

    def get_leaf_children(self, parent_names: set, max_depth: int = -1) -> set:
        """Retrieve all leaves max_depth below parents."""
        parent_numbers = set(self.get_node_numbers_of(list(parent_names)))
        leaf_numbers = list(self._get_leaf_children(parent_numbers, max_depth))
        return set(self.get_node_names_of(leaf_numbers))

    def _get_leaf_children(self, parent_names: set, max_depth: int = -1) -> set:
        """Retrieve all leaves max_depth below parents."""
        return set([child for child in self._get_children(parent_names, max_depth) if self._is_leaf_node(child)])

    def get_children(self, parent_names: set, max_depth: int = -1) -> set:
        """Retrieve all nodes max_depth below parents."""
        parent_numbers = set(self.get_node_numbers_of(list(parent_names)))
        children_numbers = self._get_children(parent_numbers, max_depth)
        return set(self.get_node_names_of(list(children_numbers)))

    def _get_children(self, parent_numbers: set, max_depth: int = -1) -> set:
        """Retrieve all nodes max_depth below parents."""
        if max_depth == 0:
            return parent_numbers
        next_generation = set()
        for parent_number in parent_numbers:
            next_generation.update(set(self._get_out_nodes_node(parent_number)))
        if next_generation:
            parent_numbers.update(self._get_children(next_generation, max_depth - 1))
            return parent_numbers
        return parent_numbers

    def get_parents(self, children_names: set, max_height: int = -1) -> set:
        """Retrieve all nodes max_depth above children."""
        children_numbers = set(self.get_node_numbers_of(list(children_names)))
        parent_numbers = self._get_parents(children_numbers, max_height)
        return set(self.get_node_names_of(list(parent_numbers)))

    def _get_parents(self, children_numbers: set, max_height: int = -1) -> set:
        if max_height == 0:
            return children_numbers
        prev_generation = set()
        for child_number in children_numbers:
            prev_generation.update(set(self._get_in_nodes_node(child_number)))
        if prev_generation:
            children_numbers.update(self._get_parents(prev_generation, max_height - 1))
            return children_numbers
        return children_numbers

    def get_connected_components(self, leafless: bool = False) -> list:
        """Retrieve the connected components (excluding leaves) of the network."""
        connected_components = self._get_connected_components(leafless)
        return [self.get_node_names_of(connected_component) for connected_component in connected_components]

    def _get_connected_components(self, leafless: bool = False) -> list:
        unchecked_nodes = set(self.node_name_map.inverse)
        components = []
        while unchecked_nodes:
            node_number = unchecked_nodes.pop()
            component = [node_number]
            unchecked_connections = set(self._get_connections_node(node_number)).intersection(unchecked_nodes)
            while unchecked_connections:
                connection_number = unchecked_connections.pop()
                unchecked_nodes.remove(connection_number)
                new_connections = set(self._get_connections_node(connection_number))
                unchecked_connections = (unchecked_connections.union(new_connections)).intersection(unchecked_nodes)
                component.append(connection_number)
            if not (leafless and len(component) == 1 and self._is_leaf_node(component[0])):
                components.append(component)
        return components

    def optimization_available(self) -> bool:
        return self._biconnected_components is not None and self._partial_ordering is not None and self._cut_arc_matrix is not None and self._cut_arc_sets is not None

    def calculate_optimzation_variables(self) -> None:
        self.cut_arc_matrix
        self._get_cut_arc_sets()
        self._get_biconnected_components()
        self._get_partial_ordering()

    def reset_optimization_variables(self) -> None:
        self._biconnected_components = None
        self._partial_ordering = None
        self._cut_arc_matrix = None
        self._cut_arc_sets = None

    @property
    def biconnected_components(self) -> list:
        """Compute biconnected components."""
        return [self.get_node_names_of(biconnected_component) for biconnected_component in self._get_biconnected_components(leafless=False)]

    def _get_biconnected_components(self, leafless: bool = False) -> list:
        """Compute biconnected components."""
        if self._biconnected_components is None:
            cut_arc_matrix = self.cut_arc_matrix
            self.adj_matrix -= cut_arc_matrix - cut_arc_matrix.T
            self._biconnected_components = self._get_connected_components(leafless=False)
            self.adj_matrix += cut_arc_matrix - cut_arc_matrix.T

        if leafless:
            return [component for component in self._biconnected_components if not (len(component) == 1 and self._is_leaf_node(component[0]))]
        else:
            return self._biconnected_components

    def is_biconnected(self, leafless: bool = False) -> bool:
        """Check if network is biconnected"""
        return len(self._get_biconnected_components(leafless=leafless)) == 1

    @property
    def cut_arc_sets(self) -> list:
        """Retrieve cut-arc sets of network."""
        return [self.get_node_names_of(cut_arc_set) for cut_arc_set in self._get_cut_arc_sets()]

    def _get_cut_arc_sets(self, leaves_only=True) -> list:
        """Retrieve cut-arc sets of network."""
        if self._cut_arc_sets is None:
            cut_arc_matrix = self.cut_arc_matrix
            _, to_nodes = np.where(cut_arc_matrix == 1)
            self._cut_arc_sets = []
            for to_node in to_nodes:
                ca_set = list(self._get_children({to_node}))
                self._cut_arc_sets.append(ca_set)
        if leaves_only:
            return [[node_number for node_number in cut_arc_set if self._is_leaf_node(node_number)] for cut_arc_set in self._cut_arc_sets]
        return self._cut_arc_sets

    @property
    def cut_arc_matrix_old(self) -> np.ndarray:
        if self._cut_arc_matrix is None:
            self._cut_arc_matrix = np.zeros((self.number_of_nodes, self.number_of_nodes))
            edges = self._get_edges()
            for edge in edges:
                from_node, to_node = edge
                self.adj_matrix[from_node][to_node] -= 1
                self.adj_matrix[to_node][from_node] += 1
                components = self.get_connected_components()
                if len(components) > 1:
                    self._cut_arc_matrix[from_node][to_node] = 1
                    self._cut_arc_matrix[to_node][from_node] = 1
                self.adj_matrix[from_node][to_node] += 1
                self.adj_matrix[to_node][from_node] -= 1
        return self._cut_arc_matrix

    # TODO fix this
    @property
    def cut_arc_matrix(self) -> np.ndarray:
        """Compute indicator matrix for arcs which are cut-arcs."""
        if self._cut_arc_matrix is None:
            visited = [False] * self.number_of_nodes
            disc = [-1] * self.number_of_nodes
            low = [-1] * self.number_of_nodes
            parent = [None] * self.number_of_nodes
            self._cut_arc_matrix = np.zeros((self.number_of_nodes, self.number_of_nodes))

            for i in range(self.number_of_nodes):
                if not visited[i]:
                    self.cut_arc_helper(i, visited, disc, low, parent, self._cut_arc_matrix)

        return self._cut_arc_matrix

    def cut_arc_helper(self, u_number, visited, disc, low, parent, cut_arc_matrix, t=0):
        # TODO source: https://www.geeksforgeeks.org/bridge-in-a-graph/
        visited[u_number] = True
        disc[u_number] = t
        low[u_number] = t
        t += 1

        for v_number in self._get_connections_node(u_number):
            if not visited[v_number]:
                parent[v_number] = u_number
                self.cut_arc_helper(v_number, visited, disc, low, parent, cut_arc_matrix, t)

                low[u_number] = min(low[u_number], low[v_number])

                if low[v_number] > disc[u_number]:
                    cut_arc_matrix[u_number][v_number] = 1
            elif v_number != parent[u_number]:
                low[u_number] = min(low[u_number], disc[v_number])

    # ---------------------------------------------------------------- NODE TYPE METHODS ---------------------------------------------------------------- #
    def is_leaf_node(self, node_name: str) -> bool:
        """Check if node_name is leaf."""
        return self._is_leaf_node(self.get_node_number_of(node_name))

    def _is_leaf_node(self, node_number: int) -> bool:
        return node_number in self.leaf_numbers

    def is_reticulation_node(self, node_name: str) -> bool:
        """Check if node_name is reticulation."""
        return self._is_reticulation_node(self.get_node_number_of(node_name))

    def _is_reticulation_node(self, node_number: int) -> bool:
        """Check if node_name is reticulation."""
        return self._get_in_degree_node(node_number) > 1

    def is_root_node(self, node_name: str) -> bool:
        """Check if node_name is root."""
        return self._is_root_node(self.get_node_number_of(node_name))

    def _is_root_node(self, node_number: int) -> bool:
        return self._get_in_degree_node(node_number) == 0

    def is_structure_node(self, node_name: str) -> bool:
        """Check if node_name is necessary for the structure: not directly connected to any leaves."""
        return self._is_structure_node(self.get_node_number_of(node_name))

    def _is_structure_node(self, node_number: int) -> bool:
        """Check if node_name is necessary for the structure: not directly connected to any leaves."""
        return (not self._is_leaf_node(node_number)) and (
                self._get_out_degree_node(node_number, leafless=False) - self._get_out_degree_node(node_number, leafless=True) == 0)

    def is_connection_node(self, node_name: str) -> bool:
        """Check if node node_name is not part of underlying generator or leaf"""
        return self._is_connection_node(self.get_node_number_of(node_name))

    def _is_connection_node(self, node_number: int) -> bool:
        """Check if node node_name is not part of underlying generator or leaf"""
        return self._get_node_type(node_number) == "connection"

    def is_generator_node(self, node_name: str) -> bool:
        """Check if node node_name is not part of underlying generator or leaf"""
        return self._is_generator_node(self.get_node_number_of(node_name))

    def _is_generator_node(self, node_number: int) -> bool:
        """Check if node_name is part of the underlying generator."""
        self.logger.debug("Checking if {} is a generator node.".format(node_number))
        return self._get_node_type(node_number) in ("structure", "reticulation", "root")

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
        if self._is_structure_node(node_number):
            return "structure"
        return "connection"

    def get_generator_nodes(self) -> list:
        """Retrieve all generator nodes"""
        return self.get_node_names_of(self._get_generator_nodes())

    def _get_generator_nodes(self) -> list:
        result = []
        for node in self.node_name_map.values():
            if self._is_generator_node(node):
                result.append(node)

    def get_edges(self, leafless: bool = False) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        return [(self.get_node_name_of(from_node), self.get_node_name_of(to_node)) for from_node, to_node in self._get_edges(leafless)]

    def _get_edges(self, leafless: bool = False) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        mask = self.leaf_mask if leafless else np.ones(self.number_of_nodes)
        i = 1
        from_nodes, to_nodes = [], []
        while np.any(np.multiply(self.adj_matrix, mask) >= i):
            extra_from_nodes, extra_to_nodes = np.where(np.multiply(self.adj_matrix, mask) >= i)
            to_nodes.extend(extra_to_nodes)
            from_nodes.extend(extra_from_nodes)
            i += 1
        return [(from_node, to_node) for from_node, to_node in zip(np.array(from_nodes), np.array(to_nodes))]

    def prune(self, suppress_redundant: str = 'all', suppress_parallel: bool = True):
        """Suppress al unnecessary/redundant/parallel nodes/edges in network."""
        stable_ancestors, lowest_stable_ancestor, nodes_between = self._stable_ancestors(self.leaf_numbers)

        for node_number in sorted(list(self.node_name_map.values()), reverse=True):
            if node_number not in nodes_between:
                self._remove_node_from_network(node_number, force=True)

        changes = True
        while changes:
            changes = False
            leaf_numbers = set(self.leaf_numbers)
            in_degrees, out_degrees = np.array(self.get_in_degrees()), np.array(self.get_out_degrees())

            # Get non-leaf nodes without children
            childless_nodes = set(np.where(out_degrees == 0)[0])
            childless_nodes.difference_update(leaf_numbers)

            # Get nodes with no entering arcs and only one exiting arc
            parentless_nodes = set(np.where(in_degrees == 0)[0])
            one_exiting_arc_nodes = set(np.where(out_degrees == 1)[0])
            parentless_one_exiting_arc_nodes = parentless_nodes.intersection(one_exiting_arc_nodes)

            # Get nodes with one in and one out node
            one_entering_arc_nodes = set(np.where(in_degrees == 1)[0])
            one_entering_and_exiting_arc_nodes = one_entering_arc_nodes.intersection(one_exiting_arc_nodes)

            nodes_to_remove = list(set(list(childless_nodes) + list(parentless_one_exiting_arc_nodes) + list(one_entering_and_exiting_arc_nodes)))

            for node in sorted(nodes_to_remove, reverse=True):
                changes = True
                self._suppress_node(node)

            # Check for parallel arcs
            if suppress_parallel:
                parallel_arcs = np.where(self.adj_matrix >= 2)
                parallel_arcs = zip(parallel_arcs[0], parallel_arcs[1])
                for from_node, to_node in parallel_arcs:
                    changes = True
                    self._remove_connection(from_node, to_node)

            # Check for redundant components:
            if suppress_redundant != 'none':
                bcs_to_remove = []
                bcs = self._get_biconnected_components(leafless=True)
                for bc in bcs:
                    if suppress_redundant == 'strongly' and self._is_strongly_redundant_component(bc):
                        bcs_to_remove.append(bc)
                    if suppress_redundant == 'all' and self._is_redundant_component(bc):
                        bcs_to_remove.append(bc)
                removed_nodes = []
                for bc in bcs_to_remove:
                    changes = True
                    bc = [node - sum(np.array(removed_nodes) < node) for node in bc]
                    self._suppress_component(bc)
                    removed_nodes.extend(bc)

    def _is_redundant_component(self, component: list) -> bool:
        """Check if component is redundant."""
        self.logger.debug("Checking if component {} is redundant.".format(component))
        out_nodes = self._get_out_nodes_component(component)
        if len(out_nodes) != 1:
            return False
        return True

    def _is_strongly_redundant_component(self, component: list) -> bool:
        """Check if component is strongly redundant."""
        self.logger.debug("Checking if component {} is strongly redundant.".format(component))
        out_nodes = self.get_out_nodes_component(component)
        if len(out_nodes) != 1:
            return False
        temp_network = copy.deepcopy(self)
        temp_network._remove_component(component)
        connected_comp = temp_network.get_connected_components()
        if len(connected_comp) > 1:
            return False
        return True

    def contains_leaf(self, component: list) -> bool:
        """Check if component contains a leaf."""
        self.logger.debug("Checking if component {} contains a leaf.".format(component))
        for node in component:
            if self.is_leaf_node(node):
                return True
        return False

    def number_of_internals_leaves_reticulations(self):
        """Retrieve number of internal, leaves, reticulation nodes."""
        self.logger.debug("Retrieving number of internal, leaves, reticulation nodes.")
        col_sum = sum(self.adj_matrix > 0)
        if type(col_sum) == int:
            col_sum = np.array([col_sum])
        number_of_reticulations = sum(col_sum > 1)
        return number_of_reticulations - len(self.leaf_numbers), len(self.leaf_numbers), number_of_reticulations

    def get_exhibited_trinets(self, max_processes=1, progress_bar=False):
        """Retrieve all trinets exhibited by network."""
        self.logger.debug("Retrieving all exhibited trinets.")

        leaves = self.leaf_numbers
        trinet_info_list = NetworkInfoList(network_size=3)
        total = int(scipy.special.comb(len(leaves), 3))
        triplets = list(itertools.combinations(leaves, 3))

        pool = multiprocessing.Pool(max_processes)
        trinet_info_list = NetworkInfoList(network_size=3)
        if progress_bar:
            for x in tqdm(pool.imap_unordered(self._exhibited_trinets_helper, triplets), total=total):
                trinet_info_list.append(x)
        else:
            for x in pool.imap_unordered(self._exhibited_trinets_helper, triplets):
                trinet_info_list.append(x)
        pool.close()
        pool.join()
        return trinet_info_list

    def _exhibited_trinets_helper(self, triplet):
        result = NetworkInfo(self.trinet_from_network(self, node_numbers=list(triplet)))
        return result

    def get_exhibited_trees(self):
        reticulation_numbers = sorted(self._get_reticulations())
        new_networks = [copy.deepcopy(self)]
        while reticulation_numbers:
            old_networks = new_networks
            new_networks = []
            current_reticulation_number = reticulation_numbers.pop(-1)
            parents_current_reticulation = self._get_parents(children_numbers=set(current_reticulation_number), max_height=1).difference(
                current_reticulation_number)
            for current_parent in parents_current_reticulation:
                for old_network in old_networks:
                    new_network = copy.deepcopy(old_network)
                    new_network._remove_connection(current_parent, current_reticulation_number)
                    new_networks.append(new_network)

        network_info_list = NetworkInfoList(network_size=len(self.leaf_numbers))
        for network in new_networks:
            network.prune()
            network_info_list.append(network)
        network_info_list.uniquify()
        return network_info_list

    @property
    def partial_ordering(self) -> list:
        return [self.get_node_names_of(partial_ordering) for partial_ordering in self._get_partial_ordering()]

    def _get_partial_ordering(self) -> list:
        if self._partial_ordering is None:
            root = self._get_root()
            if self._is_leaf_node(root):
                self._partial_ordering = [[]]
                return self._partial_ordering
            new_result = [[root]]
            self._partial_ordering = []
            changes = True
            while changes:
                changes = False
                self._partial_ordering = new_result
                new_result = []
                for track in self._partial_ordering:
                    new_tracks = []
                    children = self._get_out_nodes_node(track[-1])
                    for child in children:
                        if not self._is_leaf_node(child):
                            changes = True
                            new_tracks.append(track + [child])
                    for new_track in new_tracks:
                        new_result.append(new_track)
                    if len(new_tracks) == 0:
                        new_result.append(track)
        return self._partial_ordering

    def get_leaf_ordering(self) -> list:
        return [self.get_node_names_of(leaf_ordering) for leaf_ordering in self._get_leaf_ordering()]

    def _get_leaf_ordering(self):
        partial_node_ordering = self._get_partial_ordering()
        result = []
        for track in partial_node_ordering:
            leaf_track = []
            for node in track:
                children = self._get_out_nodes_node(node)
                if len(children) == 1:
                    if self._is_leaf_node(children[0]):
                        leaf_track.append(children[0])
                        result.append(leaf_track)
                else:
                    count = 0
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
                        result.append(leaf_track_1)
                        result.append(leaf_track_2)
        return [list(x) for x in set(tuple(x) for x in result)]

    def _stable_ancestors(self, node_numbers: list):
        leafless_node_numbers = []
        for node_number in node_numbers:
            if self._is_leaf_node(node_number):
                parent_of_leaf = self._get_parents({node_number}, max_height=1)
                parent_of_leaf.remove(node_number)
                node_number = parent_of_leaf.pop()
            leafless_node_numbers.append(node_number)
        partial_ordering = self._get_partial_ordering()
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

    def stable_ancestors(self, node_names: list):
        node_numbers = self.get_node_numbers_of(node_names)
        stable_ancestors_indicis_dict, lowest_stable_ancestor, nodes_between = self._stable_ancestors(node_numbers)
        # TODO: To names

    def visualize(self, file_path: str = None, internal_node_labels=True, edge_labels=False):
        """Visualize network."""
        self.logger.debug("Visualizing network.")
        dot = Digraph()
        dot.engine = 'dot'
        edges = self.get_edges()
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
            dot.format = 'png'
            if file_path[-4:] == '.png':
                dot.render(filename=file_path[:-4], format='png')
            else:
                dot.render(filename=file_path, format='png')
        else:
            dot.render(view=True)
            time.sleep(0.2)

    def to_df(self, directed: bool = True) -> pd.DataFrame:
        """Retrieve dataframe representation of network."""
        self.logger.debug("Retrieve {}directed dataframe representation of network.".format("" if directed else "un-"))
        ordered_node_names = self.node_names
        mask = self.adj_matrix > 0
        data = self.adj_matrix * mask if directed else self.adj_matrix
        length = self.number_of_nodes - len(self.leaf_numbers) if directed else self.number_of_nodes
        result = pd.DataFrame(columns=ordered_node_names, index=ordered_node_names[:length], data=data[:length])
        return result.astype(int)

    def _get_component_list_names(self, component_list: list) -> list:
        """Retrieve names of list of lists of node numbers."""
        self.logger.debug("Retrieve node names of {}.".format(component_list))
        return [self.get_node_name_of(comp) for comp in component_list]

    def equal_structure_old(self, other):
        self_partial_ordering = self._get_partial_ordering()
        other_partial_ordering = other._get_partial_ordering()

        if len(self_partial_ordering) != len(other_partial_ordering):
            return False, bidict()

        s1 = sum([len(row) for row in self_partial_ordering])
        s2 = sum([len(row) for row in other_partial_ordering])
        if s1 != s2:
            return False, bidict()

        self_numbers = set().union(*self_partial_ordering)
        other_numbers = set().union(*other_partial_ordering)
        if len(self_numbers) != len(other_numbers):
            return False, bidict()

        other_number_iterator = itertools.permutations(other_numbers)
        for other_name_iter in other_number_iterator:
            translation_dict = {other_number: self_number for other_number, self_number in zip(other_name_iter, self_numbers)}
            temp_other_partial_ordering = [[translation_dict[other_number] for other_number in row] for row in other_partial_ordering]
            if self._equal_structure_old(self_partial_ordering, temp_other_partial_ordering):
                self_leaf_numbers = self.leaf_numbers
                translation_dict = bidict(translation_dict).inverse

                relation_dict = bidict()
                used_leaves = []
                for self_leaf_number in self_leaf_numbers:
                    parents_self_leaf = self._get_parents({self_leaf_number}, max_height=1)
                    parents_self_leaf.remove(self_leaf_number)
                    parent_self_leaf = parents_self_leaf.pop()

                    parent_other_leaf = translation_dict[parent_self_leaf]
                    other_leaves = other._get_children({parent_other_leaf}, max_depth=1)
                    other_leaves.remove(parent_other_leaf)

                    chosen_leaf = other_leaves.pop()

                    if chosen_leaf in used_leaves or not other._is_leaf_node(chosen_leaf):
                        chosen_leaf = other_leaves.pop()
                    used_leaves.append(chosen_leaf)
                    relation_dict.put(self_leaf_number, chosen_leaf)

                return True, bidict(
                    {self.get_node_name_of(self_leaf_number): other.get_node_name_of(other_leaf_number) for self_leaf_number, other_leaf_number in
                     relation_dict.items()})
        return False, bidict()

    @staticmethod
    def _equal_structure_old(partial_ordering, other_partial_ordering):
        n = len(partial_ordering)
        it = itertools.permutations(range(n))
        for permutation in it:
            for i, j in enumerate(permutation):
                if partial_ordering[i] != other_partial_ordering[j]:
                    break
            else:
                return True
        return False

    def equal_structure(self, other, optimize=True, equal_naming=False, progress_bar=False) -> (bool, list, list):
        if self.number_of_nodes != other.number_of_nodes:
            return False, [], []
        if len(self.leaf_numbers) != len(other.leaf_numbers):
            return False, [], []
        if optimize or (self.optimization_available() and other.optimization_available()):
            self_cut_arc_sets = self._get_cut_arc_sets(leaves_only=False)
            other_cut_arc_sets = other._get_cut_arc_sets(leaves_only=False)
            if len(self_cut_arc_sets) != len(other_cut_arc_sets):
                return False, [], []
            if collections.Counter([len(cut_arc_set) for cut_arc_set in self_cut_arc_sets]) != \
                    collections.Counter([len(cut_arc_set) for cut_arc_set in other_cut_arc_sets]):
                return False, [], []

        if progress_bar:
            pbar = tqdm(total=len(self.node_name_map))
        else:
            pbar = None
        queue = multiprocessing.Queue()
        self_root = self._get_root()
        other_root = other._get_root()
        translation_dict = bidict()
        translation_dict[self_root] = other_root
        # process = multiprocessing.Process(self._equal_structure, kwargs={'self_current_node': self_root, 'other': other, 'other_current_node': other_root,
        #                                                                  'translation_dicts': [translation_dict], 'pbar': pbar, 'optimize': optimize,
        #                                                                  'equal_naming'     : equal_naming, 'queue': queue})

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

        self_node_children = self._get_children({self_current_node}, max_depth=1).difference({self_current_node})
        other_node_children = other._get_children({other_current_node}, max_depth=1).difference({other_current_node})
        if len(other_node_children) != len(self_node_children):
            return []  # False
        if optimize or (self.optimization_available() and other.optimization_available()):
            self_cut_arc_sets = self._get_cut_arc_sets(leaves_only=False)
            other_cut_arc_sets = other._get_cut_arc_sets(leaves_only=False)
            self_node_cut_arc_sets = [len(cut_arc_set) for cut_arc_set in self_cut_arc_sets if self_current_node in cut_arc_set]
            other_node_cut_arc_sets = [len(cut_arc_set) for cut_arc_set in other_cut_arc_sets if other_current_node in cut_arc_set]
            if collections.Counter(self_node_cut_arc_sets) != collections.Counter(other_node_cut_arc_sets):
                return []  # False
        if len(self_node_children) == 1:
            self_child = self_node_children.pop()
            other_child = other_node_children.pop()
            new_translation_dicts = []
            for translation_dict in translation_dicts:
                if not check_bidict(translation_dict, self_child, other_child):
                    continue
                new_translation_dict = copy.copy(translation_dict)
                new_translation_dict[self_child] = other_child
                new_dicts = self._equal_structure(self_child, other, other_child, [new_translation_dict], pbar, optimize=optimize, equal_naming=equal_naming)
                for new_dict in new_dicts:
                    if len(new_dict) == self.number_of_nodes:
                        return [new_dict]
                new_translation_dicts.extend(new_dicts)
            return new_translation_dicts
        if len(self_node_children) == 2:
            translation_dicts_to_return = []
            self_child_1 = self_node_children.pop()
            self_child_2 = self_node_children.pop()
            other_child_1 = other_node_children.pop()
            other_child_2 = other_node_children.pop()

            # Option 1
            for (other_child_option_1, other_child_option_2) in [(other_child_1, other_child_2), (other_child_2, other_child_1)]:
                for translation_dict in translation_dicts:
                    if not (check_bidict(translation_dict, self_child_1, other_child_option_1) and check_bidict(translation_dict, self_child_2,
                                                                                                                other_child_option_2)):
                        continue
                    new_translation_dict = copy.copy(translation_dict)
                    new_translation_dict[self_child_1] = other_child_option_1
                    new_translation_dict[self_child_2] = other_child_option_2
                    returned_translation_dicts_1 = self._equal_structure(self_child_1, other, other_child_option_1, [new_translation_dict], pbar,
                                                                         optimize=optimize, equal_naming=equal_naming)
                    for returned_translation_dict_1 in returned_translation_dicts_1:
                        new_dicts = self._equal_structure(self_child_2, other, other_child_option_2, [returned_translation_dict_1], pbar, optimize=optimize,
                                                          equal_naming=equal_naming)
                        for new_dict in new_dicts:
                            if len(new_dict) == self.number_of_nodes:
                                return [new_dict]
                        translation_dicts_to_return.extend(new_dicts)
            return translation_dicts_to_return

    def evolve_times(self, times: int, recombination_chance: int = 0) -> str:
        result = ""
        while times > 0:
            recombined = False
            if times == 1:
                recombined, res = self.evolve(0)
            else:
                recombined, res = self.evolve(recombination_chance)
            times -= 1 + int(recombined)
            result += res

        changes = True
        while changes:
            changes = False
            in_degrees = np.array(self.get_in_degrees())
            biconnected_components = self._get_biconnected_components(leafless=True)
            for bc in biconnected_components:
                in_sum_bc = in_degrees[bc]
                bc_reticulations = {bc[node_number] for node_number in np.where(in_sum_bc == 2)[0]}
                if len(bc_reticulations) - self.level > 0:
                    changes = True
                    reticulation_to_remove = bc_reticulations.pop()
                    parent_reticulation = self._get_parents({reticulation_to_remove}, max_height=1).difference([reticulation_to_remove]).pop()
                    child_reticulation = self._get_children({reticulation_to_remove}, max_depth=1).difference([reticulation_to_remove]).pop()
                    self._remove_connection(parent_reticulation, reticulation_to_remove)
                    self._remove_connection(reticulation_to_remove, child_reticulation)
                    self._add_connection(parent_reticulation, child_reticulation)
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
                biconnected_components = self._get_biconnected_components(leafless=True)
                for bc in biconnected_components:
                    in_sum_bc = in_degrees[bc]
                    bc_level = sum(in_sum_bc > 1)
                    bc = self._get_leaf_children(set(bc), max_depth=1)
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
                pass
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
            self._remove_leaf_status_node(leaf_number)
        self.prune()
        self.standardize_node_names()

    def enewick(self):
        return self._enewick(self._get_root())

    def _enewick(self, current_node, traversed_nodes=None):
        traversed_nodes = coalesce(traversed_nodes, [])
        if current_node in traversed_nodes:
            return self.get_node_name_of(current_node)
        traversed_nodes.append(current_node)
        current_node_children = self._get_children({current_node}, max_depth=1).difference([current_node])
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


class NetworkInfo:
    def __init__(self, network, info: dict = None):
        self.network = network
        self.info = coalesce(info, dict())

    def shrink(self, leaf_set):
        can_shrink = self.network.shrink(leaf_set)
        return bool(can_shrink)

    def calculate_info(self):
        self['cut_arc_sets'] = self.network.cut_arc_sets
        self['strict_level'] = self.network.strict_level()
        self['leaf_order'] = self.network.get_leaf_ordering()

    def add_info(self, network_info):
        for key, value in network_info.info.items():
            self.info[key] = value

    @classmethod
    def limit_to(cls, network_info, leaf_names):
        network = RootedLevelKNetwork.from_network(network_info.network, node_names=leaf_names, suppress_parallel=True, suppress_redundant='all')
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
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.network = copy.deepcopy(self.network)
        cp.info = copy.deepcopy(self.info)
        return cp

    def __str__(self):
        return str(self.network.leaf_names) + "\n" + pp.pformat(self.info)


class NetworkInfoList:
    def __init__(self, network_size: int, lst: list = None):
        self.uid = guid()
        self.logger = logging.getLogger('network_info_list.{}'.format(self.uid))
        self.logger.debug(f"Created NetworkInfoList ({self.uid})")
        self.list = coalesce(lst, [])
        self.network_size = network_size

    @classmethod
    def induced_n_network_info_list(cls, network_info_list, network_size):
        assert network_size < network_info_list.network_size, "Can only create smaller network info lists with smaller network size."
        result = cls(network_size)
        for network_info in network_info_list:
            all_leaf_numbers = network_info.network.leaf_numbers
            leaf_numbers_iterator = itertools.combinations(all_leaf_numbers, network_size)
            for leaf_numbers in leaf_numbers_iterator:
                network = RootedLevelKNetwork.from_network(network=network_info.network, node_numbers=list(leaf_numbers), suppress_parallel=True,
                                                           suppress_redundant='all')
                result.append(NetworkInfo(network))
        return result

    def uniquify(self, equal_naming=False):
        indici_to_remove = set()
        for index_1, network_info_1 in enumerate(self):
            if index_1 in indici_to_remove:
                continue
            for index_2, network_info_2 in enumerate(self[index_1 + 1:]):
                equal_structure, _, _ = network_info_1.network.equal_structure(network_info_2.network, equal_naming=equal_naming)
                if equal_structure:
                    indici_to_remove.add(index_1 + index_2 + 1)

        for index in sorted(list(indici_to_remove), reverse=True):
            self.pop(index)

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

    def move_leaf_distort(self, prob):
        assert 0 <= prob <= 1, "Probability must be between 0 and 1"
        amount = np.random.binomial(len(self), prob)
        networks_to_distort = np.random.choice(range(len(self)), amount, replace=False)
        for index in networks_to_distort:
            leaf_to_move = np.random.choice(self[index].network.leaf_names)
            self[index].network.remove_leaf(leaf_to_move)
            edges = self[index].network.get_edges()
            edge_to_add_leaf_to_index = np.random.choice(range(len(edges)))
            edge_to_add_leaf_to = edges[edge_to_add_leaf_to_index]
            self[index].network.add_leaf_to_edge(edge_to_add_leaf_to, leaf_to_move)

    def switch_leaves_distort(self, prob):
        assert 0 <= prob <= 1, "Probability must be between 0 and 1"
        amount = np.random.binomial(len(self), prob)
        networks_to_distort = np.random.choice(range(len(self)), amount, replace=False)
        for index in networks_to_distort:
            leaves_to_switch = np.random.choice(self[index].network.leaf_names, size=2, replace=False)
            self[index].network.rename_node(leaves_to_switch[0], leaves_to_switch[0] + "*")
            self[index].network.rename_node(leaves_to_switch[1], leaves_to_switch[0])
            self[index].network.rename_node(leaves_to_switch[0] + "*", leaves_to_switch[1])

    def remove_network_distort(self, prob):
        assert 0 <= prob <= 1, "Probability must be between 0 and 1"
        amount = np.random.binomial(len(self), prob)
        networks_to_distort = np.random.choice(range(len(self)), amount, replace=False)
        for index in sorted(networks_to_distort, reverse=True):
            self.pop(index)

    def _analyse_networks(self, all_taxa) -> (dict, float, float, float):
        n_set_iterator = itertools.combinations(all_taxa, self.network_size)
        covered_networks = dict()
        for n_set in n_set_iterator:
            covered_networks[tuple(sorted(n_set))] = []

        for network_info in self:
            taxa = tuple(sorted(network_info.network.leaf_names))
            try:
                networks = [network_count['network'] for network_count in covered_networks[taxa]]
            except KeyError:
                continue  # In case taxa is not in covered_networks. Happens when there are extra_taxa
            for index, network in enumerate(networks):
                equal, _, _ = network_info.network.equal_structure(network, equal_naming=True)
                if equal:
                    covered_networks[taxa][index]['count'] += 1
                    break
            else:
                covered_networks[taxa].append({'network': network_info.network, 'count': 1})

        unique_networks = 0.
        unique_n_sets = 0.
        for taxa, network_counts in covered_networks.items():
            l = len(network_counts)
            unique_networks += l
            unique_n_sets += int(l > 0)

        max_number_of_networks = ncr(len(all_taxa), self.network_size)
        redundancy = (len(self) - unique_networks) / max_number_of_networks
        density = unique_n_sets / max_number_of_networks
        inconsistency = (unique_networks - unique_n_sets) / unique_n_sets

        return covered_networks, redundancy, density, inconsistency

    def summary(self, all_taxa: list = None, per_network=False) -> str:
        represented_taxa = self.represented_leaves()
        all_taxa = coalesce(all_taxa, list(represented_taxa))
        extra_taxa = represented_taxa.difference(all_taxa)
        missing_taxa = set(all_taxa).difference(represented_taxa)

        covered_networks, redundancy, density, inconsistency_0 = self._analyse_networks(all_taxa)

        network_info_list_1 = NetworkInfoList.induced_n_network_info_list(self, self.network_size - 1)
        covered_networks_1, _, _, inconsistency_1 = network_info_list_1._analyse_networks(all_taxa)

        summary = \
            f"NetworkInfoList ({self.uid}) the following properties for taxa \n" \
            f"        {all_taxa}: \n" \
            f" - represents extra taxa: \n" \
            f"        {extra_taxa} \n" \
            f" - completely misses taxa: \n" \
            f"        {missing_taxa} \n" \
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

    def is_consistent(self, all_taxa: list = None) -> bool:
        represented_taxa = self.represented_leaves()
        all_taxa = coalesce(all_taxa, list(represented_taxa))
        _, redundancy, density, inconsistency_0 = self._analyse_networks(all_taxa)
        network_info_list_1 = NetworkInfoList.induced_n_network_info_list(self, self.network_size - 1)
        covered_networks_1, _, _, inconsistency_1 = network_info_list_1._analyse_networks(all_taxa)
        return bool(redundancy == 0.0 and density == 1.0 and inconsistency_0 == 0.0 and inconsistency_1 == 0.0)

    def calculate_info(self):
        for network_info in self:
            network_info.calculate_info()

    def add_info(self, network_info_list):
        for network_info in self:
            equal_structured_networks, relation_dicts = network_info_list.find_equal_structured_network(network_info.network)
            # TODO check if two?
            try:
                equal_structured_network = equal_structured_networks[0]
                network_info.add_info(equal_structured_network)
                network_info['relation_dict'] = relation_dicts[0]
                network_info['leaf_order'] = network_info.network.get_leaf_ordering()
            except IndexError:
                network_info.calculate_info()

    def shrink(self, leaf_set):
        to_remove = []
        for network_info in self:
            if not network_info.shrink(leaf_set):
                to_remove.append(network_info)
        for network_info in to_remove:
            self.remove(network_info)

    def append(self, other):
        assert type(other) == NetworkInfo, "Can only add object of type NetworkInfo to a NetworkInfoList"
        assert len(
            other.network.leaf_names) <= self.network_size, f"Can only add network with less or equal than {self.network_size} leaves to this NetworkInfoList, this network has {len(other.network.leaf_names)} leaves"
        self.list.append(other)

    def extend(self, other):
        for o in other:
            self.append(o)

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
        # if len(leaf_names) > 2:
        #     for trinet_info in self:
        #         if set(trinet_info.trinet.leaf_names).issubset(leaf_names):
        #             result.append(trinet_info)
        for network_info in self:
            overlap = set(leaf_names).intersection(network_info.network.leaf_names)
            if len(overlap) >= 2:
                result.append(NetworkInfo.limit_to(network_info, list(overlap)))
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

    def represented_leaves(self) -> set:
        result = set()
        for network_info in self:
            result.update(network_info.network.leaf_names)
        return result

    def weighted_auxiliary_digraph(self):
        weighted_auxiliary_digraph = dict()
        for network_info in self:
            cut_arc_sets = network_info['cut_arc_sets']
            for cut_arc_set in cut_arc_sets:
                if len(cut_arc_set) == 2:
                    x = cut_arc_set[0]
                    y = cut_arc_set[1]
                    z = [Z for Z in network_info.network.leaf_names if Z not in cut_arc_set][0]
                    for i in (x, y):
                        try:
                            weighted_auxiliary_digraph[(i, z)] += 1
                        except KeyError:
                            weighted_auxiliary_digraph[(i, z)] = 1
        return weighted_auxiliary_digraph

    def auxiliary_digraph(self, weighted_auxiliary_digraph, level):
        leaf_set = self.represented_leaves()
        auxiliary_digraph = dict()
        arc_iterator = itertools.permutations(leaf_set, 2)
        for arc in arc_iterator:
            weight = 0
            try:
                weight = weighted_auxiliary_digraph[arc]
            except KeyError:
                pass
            if weight <= level:
                try:
                    auxiliary_digraph[arc[0]].add(arc[1])
                except KeyError:
                    auxiliary_digraph[arc[0]] = {arc[1]}
        return auxiliary_digraph

    def minimal_sink_set(self, strongly_connected_components, auxiliary_digraph):
        msss = []
        for strongly_connected_component in strongly_connected_components:
            for node in strongly_connected_component:
                try:
                    to_leaves = auxiliary_digraph[node]
                    if not set(to_leaves).issubset(strongly_connected_component):
                        break
                except KeyError:
                    pass
            else:
                msss.append(strongly_connected_component)
        return [mss for mss in msss if len(mss) > 1]

    def get_minimal_sink_sets(self, level=0):
        n = len(self.represented_leaves())
        weighted_auxiliary_digraph = self.weighted_auxiliary_digraph()
        auxiliary_digraph = self.auxiliary_digraph(weighted_auxiliary_digraph, level)
        strongly_connected_components = tarjan(auxiliary_digraph)
        msss = self.minimal_sink_set(strongly_connected_components, auxiliary_digraph)
        while len(msss) == 0 and level <= n:
            level += 1
            auxiliary_digraph = self.auxiliary_digraph(weighted_auxiliary_digraph, level)
            strongly_connected_components = tarjan(auxiliary_digraph)
            msss = self.minimal_sink_set(strongly_connected_components, auxiliary_digraph)
        score = 1 - level / min([len(mss) for mss in msss])
        return msss, score

    def max_level(self):
        return max([0] + [network_info['strict_level'] for network_info in self])

    def best_level(self, max_level: int = 2):
        if max_level > 2:
            self.logger.warning("This code is not intended to work for max_level > 2")
        level_count = {level: 0 for level in range(max_level + 1)}
        for network_info in self:
            try:
                level_count[network_info['strict_level']] += 1
            except KeyError:  # when level is higher than max_level
                pass

        number_of_leaves = len(self.represented_leaves())
        number_of_possible_networks = ncr(number_of_leaves, 3)
        number_of_networks = len(self)
        level_percentage = {level: count / number_of_networks for level, count in level_count.items()}
        # TODO better naming and better boundaries

        min_percentage_matrix = np.zeros((3, 3))
        max_percentage_matrix = np.zeros((3, 3))
        # Level 0
        min_percentage_matrix[0][0] = 1.
        max_percentage_matrix[0][0] = 1.

        # Level 1
        temp = ncr(number_of_leaves - 1, 2) / number_of_possible_networks
        min_percentage_matrix[1][0] = 1 - temp
        max_percentage_matrix[1][0] = 1 - temp
        min_percentage_matrix[1][1] = temp
        max_percentage_matrix[1][1] = temp

        # Level 2
        min_percentage_matrix[2][0] = min(
            (ncr(int(math.ceil((number_of_leaves - 1) / 2.)), 3) + ncr(int(math.floor((number_of_leaves - 1) / 2.)), 3)),
            ncr(number_of_leaves - 2, 3)
        ) / number_of_possible_networks
        max_percentage_matrix[2][0] = ncr(number_of_leaves - 1, 3) / number_of_possible_networks
        min_percentage_matrix[2][2] = ncr(number_of_leaves - 2, 1) / number_of_possible_networks
        max_percentage_matrix[2][2] = ncr(number_of_leaves - 1, 2) / number_of_possible_networks
        min_percentage_matrix[2][1] = 1 - max_percentage_matrix[2][0] - max_percentage_matrix[2][2]
        max_percentage_matrix[2][1] = 1 - min_percentage_matrix[2][0] - min_percentage_matrix[2][2]

        shortage_matrix = np.zeros((3, 3))
        surplus_matrix = np.zeros((3, 3))
        for network_level in range(max_level + 1):
            for trinet_level in range(3):
                shortage_matrix[network_level][trinet_level] = max(0, min_percentage_matrix[network_level][trinet_level] - level_percentage[trinet_level])
                surplus_matrix[network_level][trinet_level] = max(0, level_percentage[trinet_level] - max_percentage_matrix[network_level][trinet_level])

        alginment_matrix = shortage_matrix + surplus_matrix
        score_array = list(np.sum(alginment_matrix, axis=1))
        best_level = score_array.index(min(score_array))
        best_score = 1-score_array[best_level]

        return best_level, best_score

    def best_generator(self):
        # TODO: in case of draws/near draws, use same heuristics as in best_level
        generators = []
        generator_count = []
        for network_info in self:
            try:
                generator_index = generators.index(network_info['generator'])
                generator_count[generator_index] += 1
            except ValueError:
                generators.append(network_info['generator'])
                generator_count.append(1)

        max_index = generator_count.index(max(generator_count))
        generator = generators[max_index]
        return copy.deepcopy(generator)

    def leaf_order_info(self):
        leaf_set = bidict({leaf: index for leaf, index in enumerate(self.represented_leaves())})
        leaf_order_matrix = np.zeros((len(leaf_set), len(leaf_set)))
        for network_info in self:

            leaf_order = network_info['leaf_order']
            for leaf_ord in leaf_order:
                for low_leaf_index in range(len(leaf_ord)):
                    for high_leaf_index in range(low_leaf_index + 1, len(leaf_ord)):
                        low_leaf = leaf_set.inverse[leaf_ord[low_leaf_index]]
                        high_leaf = leaf_set.inverse[leaf_ord[high_leaf_index]]
                        leaf_order_matrix[low_leaf, high_leaf] += 1
        return leaf_set, leaf_order_matrix

    def __add__(self, other):
        result = copy.deepcopy(self)
        for oth in other:
            result.append(oth)
        return result

    def remove(self, other):
        self.list.remove(other)

    def pop(self, index):
        self.list.pop(index)

    def __getitem__(self, key):
        return self.list[key]

    def __setitem__(self, key, value):
        self.list[key] = value

    def __iter__(self):
        return iter(self.list)

    def __len__(self):
        return len(self.list)

    def __copy__(self):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.list = copy.copy(self.list)
        cp.network_size = copy.copy(self.network_size)
        cp.uid = guid()
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.list = copy.deepcopy(self.list)
        cp.network_size = copy.deepcopy(self.network_size)
        cp.uid = guid()
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

    def get_leaf_locations(self):
        assert self.network_size <= 3, "This method only works for trinet and binet sets"
        # TODO: assert all networks have the same generator (? or use all trinets etc)
        # TODO: For generators with symmetric reticulations leaves score will not be 100%
        all_leaves = list(self.represented_leaves())
        number_of_leaves = len(all_leaves)

        # Leaf --> Edge --> (how many times leaf is on this edge)
        all_edges = list(set(self[0]['generator'].get_edges(leafless=True)))
        number_of_edges = len(all_edges)
        leaf_on_edge_count_dict = {leaf: {edge: 0 for edge in all_edges} for leaf in all_leaves}
        for network_info in self:
            relation_dict = network_info['relation_dict']
            extra_leaf_dict = network_info['extra_leaf_dict']
            for leaf, edge in extra_leaf_dict.items():
                translated_leaf = relation_dict[leaf]
                leaf_on_edge_count_dict[translated_leaf][edge] += 1
        # Put in matrix
        leaf_on_edge_count_matrix = np.zeros((number_of_edges, number_of_leaves))
        for leaf_index, leaf in enumerate(all_leaves):
            for edge_index, edge in enumerate(all_edges):
                leaf_on_edge_count_matrix[edge_index][leaf_index] = leaf_on_edge_count_dict[leaf][edge]

        # Generator reticulation --> leaf --> (how many times leaf is this reticulation leaf}
        generator_reticulations = list(self[0]['reticulations'])
        number_of_reticulations = len(generator_reticulations)
        leaf_is_reticulation_dict = {leaf: {ret: 0 for ret in generator_reticulations} for leaf in all_leaves}
        for network_info in self:
            relation_dict = network_info['relation_dict']
            for generator_leaf, leaf in relation_dict.items():
                if generator_leaf in generator_reticulations:
                    leaf_is_reticulation_dict[leaf][generator_leaf] += 1
        # Put in matrix
        leaf_is_reticulation_matrix = np.zeros((number_of_reticulations, number_of_leaves))
        for leaf_index, leaf in enumerate(all_leaves):
            for reticulation_index, generator_reticulation in enumerate(generator_reticulations):
                leaf_is_reticulation_matrix[reticulation_index][leaf_index] = leaf_is_reticulation_dict[leaf][generator_reticulation]

        # Create simplex problem
        cost_matrix = np.vstack((leaf_on_edge_count_matrix, leaf_is_reticulation_matrix))
        for j in range(number_of_leaves):
            norm = sum(cost_matrix[:, j])
            if norm != 0:
                cost_matrix[:, j] /= norm

        cost_array = cost_matrix.flatten()
        eye_list = [np.eye(number_of_leaves) for _ in range(number_of_reticulations + number_of_edges)]
        constraint_matrix_p1 = np.hstack(eye_list)
        constraint_matrix_p2 = np.zeros((number_of_reticulations, number_of_leaves * (number_of_reticulations + number_of_edges)))
        for i in range(number_of_edges, number_of_edges + number_of_reticulations):
            for j in range(number_of_leaves * i, number_of_leaves * (i + 1)):
                constraint_matrix_p2[i - number_of_edges][j] = 1
        constraint_matrix = np.vstack((constraint_matrix_p1, constraint_matrix_p2))
        b_array = [1] * (number_of_leaves + number_of_reticulations)

        solution, score = simplex_to_ILP(c=cost_array, A_eq=constraint_matrix, b_eq=b_array)

        solution = np.where(np.array(solution) == 1.0)[0]
        assert len(solution) == number_of_leaves
        edge_leaf_dict = {edge: [] for edge in all_edges}
        chosen_reticulation_relation_dict = {}
        for index in solution:
            part = int(index / number_of_leaves)
            if part < number_of_edges:
                edge_leaf_dict[all_edges[part]].append(all_leaves[index % number_of_leaves])
            else:
                chosen_reticulation_relation_dict[all_leaves[index % number_of_leaves]] = generator_reticulations[part - number_of_edges]

        return edge_leaf_dict, chosen_reticulation_relation_dict, score / number_of_leaves
