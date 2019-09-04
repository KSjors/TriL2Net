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
import time
import os
import pprint
import random
import scipy
import collections
import multiprocessing

pp = pprint.PrettyPrinter(indent=4)

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'


# GLOBALLOCK = multiprocessing.Lock()


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
        cp.adj_matrix = copy.deepcopy(self.adj_matrix)
        cp.node_name_map = copy.deepcopy(self.node_name_map)
        cp.leaf_numbers = copy.deepcopy(self.leaf_numbers)
        cp.number_of_nodes = copy.deepcopy(self.number_of_nodes)
        cp.level = copy.deepcopy(self.level)
        cp.dimension = copy.deepcopy(self.dimension)
        cp.uid = guid()
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp._biconnected_components = copy.deepcopy(self._biconnected_components)
        cp._partial_ordering = copy.deepcopy(self._partial_ordering)
        cp._cut_arc_matrix = copy.deepcopy(self._cut_arc_matrix)
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
        new_network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_numbers=taxa, level=level, dimension=dimension)
        new_network.prune(suppress_redundant=suppress_redundant, suppress_parallel=suppress_parallel)
        new_network.logger.debug("Sub-network of network {} with taxa {}.".format(network.uid, taxa))
        return new_network

    @classmethod
    def from_dir_adj_matrix(cls, dir_adj_matrix: np.ndarray, level=2, dimension=2, check_valid=True):
        """Create network from directed adjacency matrix. Assumes only internal nodes have rows. Other nodes are leaves."""
        logging.debug("Creating network from directed adjacency matrix.")
        shape = dir_adj_matrix.shape
        adj_matrix = np.zeros((shape[1], shape[1])).astype(float)
        adj_matrix[:shape[0], :shape[1]] += dir_adj_matrix
        adj_matrix[:shape[1], :shape[0]] -= dir_adj_matrix.T
        node_name_map = bidict()
        leaf_names = list()
        ln_iterator = leaf_name_iterator(1, 100)
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

    def standardize_leaf_node_names(self) -> None:
        for leaf_number in self.leaf_numbers:
            self._rename_node(leaf_number, self.get_node_name_of(leaf_number) + "^")

        for leaf_node in self.leaf_numbers:
            self._rename_node(leaf_node, self.first_unused_leaf_name())

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

        self._biconnected_components = None
        self._partial_ordering = None
        self._cut_arc_matrix = None
        self._cut_arc_sets = None

    def _remove_component(self, component: list) -> None:
        """Remove all nodes in component."""
        for node_number in sorted(component, reverse=True):
            self._remove_node_from_network(node_number=node_number, force=True)

    def _add_node_to_dict(self, node_name: str = None, leaf=False) -> (str, int):
        """Add internal node to node_names dictionary. Returns its number."""

        # Get next number, and get a node_name if no name is given
        node_number = self.number_of_nodes
        if leaf:
            node_name = self.first_unused_leaf_name() if node_name is None else node_name
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

    def _add_node_to_network(self, node_name: str = None, leaf: bool = False) -> (str, int):
        """Add internal node to network."""
        self.logger.debug("Adding internal node to network.")

        node_name, node_number = self._add_node_to_dict(node_name, leaf=leaf)
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

    def add_leaf_to_edge(self, edge, leaf_name: str = None) -> (str, str):
        """Add node between parent and child and attach leaf with leaf_name to it."""
        return self.get_node_names_of(list(self._add_leaf_to_edge(self.get_node_numbers_of(edge), leaf_name)))

    def _add_leaf_to_edge(self, edge, leaf_name: str = None) -> (int, int):
        parent = edge[0]
        child = edge[1]

        internal_name, internal_number = self._add_node_to_network()
        leaf_name, leaf_number = self._add_node_to_network(leaf_name, leaf=True)

        self._remove_connection(parent, child)

        self._add_connection(parent, internal_number)
        self._add_connection(internal_number, child)
        self._add_connection(internal_number, leaf_number)
        return internal_number, leaf_number

    def recombine_leaves(self, leaf_name_1: str, leaf_name_2: str) -> None:
        assert self.is_leaf_node(leaf_name_1) and self.is_leaf_node(leaf_name_2), "Can not recombine non-leaf nodes"
        self.get_node_name_of(self._recombine_leaves(self.get_node_number_of(leaf_name_1), self.get_node_number_of(leaf_name_2)))

    def _recombine_leaves(self, leaf_number_1: int, leaf_number_2: int) -> None:
        assert self._is_leaf_node(leaf_number_1) and self._is_leaf_node(leaf_number_2), "Can not recombine non-leaf nodes"
        self._remove_leaf_status_node(leaf_number_1)
        self._remove_leaf_status_node(leaf_number_2)

        leaf_name_1 = self.get_node_name_of(leaf_number_1)
        leaf_name_2 = self.get_node_name_of(leaf_number_2)

        new_leaf_name_1 = "(" + leaf_name_1 + "-R)"
        new_leaf_name_2 = "(" + leaf_name_2 + "-R)"
        recombined_leaf_name = "(" + leaf_name_1 + "x" + leaf_name_2 + ")"

        new_leaf_name_1, new_leaf_number_1 = self._add_node_to_network(new_leaf_name_1, leaf=True)
        new_leaf_name_2, new_leaf_number_2 = self._add_node_to_network(new_leaf_name_2, leaf=True)
        internal_leaf_name, internal_leaf_number = self._add_node_to_network(leaf=False)
        recombined_leaf_name, recombined_leaf_number = self._add_node_to_network(recombined_leaf_name, leaf=True)

        self._add_connection(leaf_number_1, new_leaf_number_1)
        self._add_connection(leaf_number_2, new_leaf_number_2)
        self._add_connection(leaf_number_1, internal_leaf_number)
        self._add_connection(leaf_number_2, internal_leaf_number)
        self._add_connection(internal_leaf_number, recombined_leaf_number)

        # parent_number_1 = self._get_parents({leaf_number_1}, max_height=1)
        # parent_number_1.remove(leaf_number_1)
        # parent_number_1 = parent_number_1.pop()
        # parent_number_2 = self._get_parents({leaf_number_2}, max_height=1)
        # parent_number_2.remove(leaf_number_2)
        # parent_number_2 = parent_number_2.pop()
        # self._remove_node_from_network(leaf_number_1, force=True)
        # self._remove_node_from_network(leaf_number_2 - int(leaf_number_1 < leaf_number_2), force=True)
        #
        # internal_node_name, internal_node_number = self._add_node_to_network(leaf=False)
        # leaf_node_name, leaf_node_number = self._add_node_to_network(recombined_leaf_name, leaf=True)
        # self._add_connection(parent_number_1 - int(leaf_number_1 < parent_number_1) - int(leaf_number_2 < parent_number_1), internal_node_number)
        # self._add_connection(parent_number_2 - int(leaf_number_1 < parent_number_2) - int(leaf_number_2 < parent_number_2), internal_node_number)
        # self._add_connection(internal_node_number, leaf_node_number)
        #
        # return leaf_node_number - int(leaf_number_1 < leaf_node_number) - int(leaf_number_2 < leaf_node_number)

    def split_leaf(self, leaf_name_to_split: str):
        """ Split leaf into two leaves"""
        assert self.is_leaf_node(leaf_name_to_split), "Can not split non-leaf-node"
        return self.get_node_names_of(list(self._split_leaf(self.get_node_number_of(leaf_name_to_split))))

    def _split_leaf(self, leaf_number_to_split: int) -> (int, int):
        self._remove_leaf_status_node(leaf_number_to_split)
        leaf_name_to_split = self.get_node_name_of(leaf_number_to_split)
        _, left_leaf_number = self._add_node_to_network(node_name='(' + leaf_name_to_split + "-I)", leaf=True)
        _, right_leaf_number = self._add_node_to_network(node_name='(' + leaf_name_to_split + "-II)", leaf=True)
        self._add_connection(leaf_number_to_split, left_leaf_number)
        self._add_connection(leaf_number_to_split, right_leaf_number)
        return left_leaf_number, right_leaf_number

    def replace_leaf_with_network(self, leaf_name_to_replace: str, replacement_network):
        assert self.is_leaf_node(leaf_name_to_replace), f"Cannot replace leaf {leaf_name_to_replace} with network as it is not a leaf."
        assert set(self.leaf_names).isdisjoint(
            replacement_network.leaf_names), f"Cannot replace leaf {leaf_name_to_replace} with network {replacement_network} as {replacement_network} has some leafs same as {self}"
        leaf_number_to_replace = self.get_node_number_of(leaf_name_to_replace)
        return self._replace_leaf_with_network(leaf_number_to_replace, replacement_network)

    def _replace_leaf_with_network(self, leaf_number_to_replace: int, replacement_network):
        assert self._is_leaf_node(leaf_number_to_replace), \
            f"Cannot replace leaf {leaf_number_to_replace} with network as it is not a leaf."
        assert set(self.leaf_names).isdisjoint(replacement_network.leaf_names), \
            f"Cannot replace leaf {leaf_number_to_replace} with network {replacement_network} as it has some leafs same as {self}"

        replacement_network = copy.deepcopy(replacement_network)
        replacement_network_internal_node_names = set(replacement_network.internal_node_names)
        replacement_network_root_name = replacement_network.get_root_name()

        # Add internal nodes from replacement network to current network with primed names
        for node_name in replacement_network.internal_node_names:
            self._add_node_to_network(node_name + "*", leaf=False)

        # Add leaves from replacement network to current network
        for leaf_name in replacement_network.leaf_names:
            self._add_node_to_network(leaf_name, leaf=True)

        # Replace leaf_name with root of replacement network
        parent_of_leaf = self._get_in_nodes_node(leaf_number_to_replace)[0]
        self._remove_node_from_network(leaf_number_to_replace)
        self._add_connection(parent_of_leaf, self.get_node_number_of(replacement_network_root_name + "*"))

        # Add all connections from replacement network to current network
        for internal_node_name in replacement_network_internal_node_names:
            to_node_names = replacement_network.get_out_nodes_node(internal_node_name)
            internal_node_number = self.get_node_number_of(internal_node_name + "*")
            for to_node_name in to_node_names:
                if to_node_name in replacement_network_internal_node_names:
                    to_node_name += "*"
                self._add_connection(internal_node_number, self.get_node_number_of(to_node_name))

    def first_unused_internal_name(self) -> str:
        """Find first unused leaf name in numerical order."""
        i = 0
        while str(i) in self.node_name_map:
            i += 1
        return str(i)

    def first_unused_leaf_name(self) -> str:
        """Find first unused leaf name in alphabetical order."""
        leaf_name_length = 1
        found = False
        leaf_name = None
        while not found:
            leaf_names = leaf_name_iterator(leaf_name_length, leaf_name_length)
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

    def get_root_number(self) -> int:
        """Retrieve the number of the root"""
        col_sum = sum(self.adj_matrix > 0)
        roots = np.where(col_sum == 0)[0]
        return roots[0]

    def get_root_name(self) -> str:
        """Retrieve the name of the root."""
        self.logger.debug("Retrieving name of root node.")
        col_sum = sum(self.adj_matrix > 0)
        roots = np.where(col_sum == 0)[0]
        return self.get_node_name_of(roots[0])

    def get_reticulations(self) -> list:
        """Retrieve node names of all reticulations."""
        self.logger.debug("Retrieving names of all reticulations.")
        in_degrees = np.array(self.get_in_degrees())
        reticulations = list(np.where(in_degrees > 1)[0])
        return self.get_node_names_of(reticulations)

    def get_leaves(self) -> list:
        """Retrieve node names of all leaves."""
        self.logger.debug("Retrieving names of all leaves.")
        return self.leaf_numbers

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
        connected_components = self._get_connections_node(leafless)
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

    def _get_cut_arc_sets(self) -> list:
        """Retrieve cut-arc sets of network."""
        if self._cut_arc_sets is None:
            cut_arc_matrix = self.cut_arc_matrix
            _, to_nodes = np.where(cut_arc_matrix == 1)
            self._cut_arc_sets = []
            for to_node in to_nodes:
                ca_set = list(self._get_leaf_children({to_node}))
                ca_set.sort()
                self._cut_arc_sets.append(ca_set)
        return self._cut_arc_sets

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
        return result

    def get_edges(self, leafless: bool = False) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        return [(self.get_node_name_of(from_node), self.get_node_name_of(to_node)) for from_node, to_node in self._get_edges()]

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

    def get_exhibited_trinets(self, max_processes=1):
        """Retrieve all trinets exhibited by network."""
        self.logger.debug("Retrieving all exhibited trinets.")

        leaves = self.leaf_numbers
        trinet_info_list = TrinetInfoList()
        total = int(scipy.special.comb(len(leaves), 3))
        triplets = list(itertools.combinations(leaves, 3))

        pool = multiprocessing.Pool(max_processes)
        trinet_info_list = TrinetInfoList()
        for x in tqdm(pool.imap_unordered(self._exhibited_trinets_helper, triplets), total=total):
            trinet_info_list.append(x)
        pool.close()
        return trinet_info_list

    def _exhibited_trinets_helper(self, triplet):
        result = TrinetInfo(self.trinet_from_network(self, node_numbers=list(triplet)))
        return result

    @property
    def partial_ordering(self) -> list:
        return [self.get_node_names_of(partial_ordering) for partial_ordering in self._get_partial_ordering()]

    def _get_partial_ordering(self) -> list:
        if self._partial_ordering is None:
            root = self.get_root_number()
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

    def visualize(self, file_path: str = None):
        """Visualize network."""
        self.logger.debug("Visualizing network.")
        dot = Digraph()
        dot.engine = 'dot'
        edges = self.get_edges()
        for node_name in self.node_name_map:
            dot.node(node_name)
        for edge in edges:
            dot.edge(edge[0], edge[1])
        if file_path:
            if file_path[-4:] == '.png':
                dot.save(file_path)
            else:
                dot.save(file_path + '.png')
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

    def equal_structure(self, other, optimize=True, equal_naming=False) -> (bool, list, list):
        if self.number_of_nodes != other.number_of_nodes:
            return False, [], []
        if len(self.leaf_numbers) != len(other.leaf_numbers):
            return False, [], []
        if optimize or (self.optimization_available() and other.optimization_available()):
            if len(self._get_cut_arc_sets()) != len(other._get_cut_arc_sets()):
                return False, [], []
            if collections.Counter([len(cut_arc_set) for cut_arc_set in self._get_cut_arc_sets()]) != \
                    collections.Counter([len(cut_arc_set) for cut_arc_set in other._get_cut_arc_sets()]):
                return False, [], []

        self_root = self.get_root_number()
        other_root = other.get_root_number()
        translation_dict = bidict()
        translation_dict[self_root] = other_root
        translation_dicts = self._equal_structure(self_current_node=self_root, other=other, other_current_node=other_root,
                                                  translation_dicts=[translation_dict], optimize=optimize, equal_naming=equal_naming)
        translation_dicts = [bidict({self.get_node_name_of(key): other.get_node_name_of(item) for key, item in translation_dict.items()})
                             for translation_dict in translation_dicts]
        leaf_translation_dicts = [bidict({key: item for key, item in translation_dict.items() if self.is_leaf_node(key)}) for translation_dict in
                                  translation_dicts]
        return len(translation_dicts) > 0, translation_dicts, leaf_translation_dicts

    def _equal_structure(self, self_current_node, other, other_current_node, translation_dicts, optimize, equal_naming):
        if len(translation_dicts) == 0:
            return translation_dicts  # False
        if self._is_leaf_node(self_current_node) and other._is_leaf_node(other_current_node):
            if (not equal_naming) or self.get_node_name_of(self_current_node) == other.get_node_name_of(other_current_node):
                for translation_dict in translation_dicts:
                    if len(translation_dict) == self.number_of_nodes:
                        return [translation_dict]
                return translation_dicts  # True

        self_node_children = self._get_children({self_current_node}, max_depth=1).difference({self_current_node})
        other_node_children = other._get_children({other_current_node}, max_depth=1).difference({other_current_node})
        if len(other_node_children) != len(self_node_children):
            return []  # False
        if optimize or (self.optimization_available() and other.optimization_available()):
            self_node_cut_arc_sets = [len(cut_arc_set) for cut_arc_set in self._get_cut_arc_sets() if self_current_node in cut_arc_set]
            other_node_cut_arc_sets = [len(cut_arc_set) for cut_arc_set in other._get_cut_arc_sets() if other_current_node in cut_arc_set]
            if collections.Counter(self_node_cut_arc_sets) != collections.Counter(other_node_cut_arc_sets):
                return []  # False
        if len(self_node_children) == 1:
            self_child = self_node_children.pop()
            other_child = other_node_children.pop()
            new_translation_dicts = []
            for translation_dict in translation_dicts:
                if not check_bidict(translation_dict, self_child, other_child):
                    continue
                new_translation_dict = copy.deepcopy(translation_dict)
                new_translation_dict[self_child] = other_child
                new_dicts = self._equal_structure(self_child, other, other_child, [new_translation_dict], optimize=optimize, equal_naming=equal_naming)
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
            for translation_dict in translation_dicts:
                if not (check_bidict(translation_dict, self_child_1, other_child_1) and check_bidict(translation_dict, self_child_2, other_child_2)):
                    continue
                new_translation_dict = copy.deepcopy(translation_dict)
                new_translation_dict[self_child_1] = other_child_1
                new_translation_dict[self_child_2] = other_child_2
                returned_translation_dicts_1 = self._equal_structure(self_child_1, other, other_child_1, [new_translation_dict], optimize=optimize, equal_naming=equal_naming)
                for returned_translation_dict_1 in returned_translation_dicts_1:
                    new_dicts = self._equal_structure(self_child_2, other, other_child_2, [returned_translation_dict_1], optimize=optimize, equal_naming=equal_naming)
                    for new_dict in new_dicts:
                        if len(new_dict) == self.number_of_nodes:
                            return [new_dict]
                    translation_dicts_to_return.extend(new_dicts)

            # Option 2
            for translation_dict in translation_dicts:
                if not (check_bidict(translation_dict, self_child_1, other_child_2) and check_bidict(translation_dict, self_child_2, other_child_1)):
                    continue
                new_translation_dict = copy.deepcopy(translation_dict)
                new_translation_dict[self_child_1] = other_child_2
                new_translation_dict[self_child_2] = other_child_1
                returned_translation_dicts_1 = self._equal_structure(self_child_1, other, other_child_2, [new_translation_dict], optimize=optimize, equal_naming=equal_naming)
                for returned_translation_dict_1 in returned_translation_dicts_1:
                    new_dicts = self._equal_structure(self_child_2, other, other_child_1, [returned_translation_dict_1], optimize=optimize, equal_naming=equal_naming)
                    for new_dict in new_dicts:
                        if len(new_dict) == self.number_of_nodes:
                            return [new_dict]
                    translation_dicts_to_return.extend(new_dicts)
            return translation_dicts_to_return

    def evolve_times(self, times: int, recombination_chance: int = 0):
        for time in range(times):
            self.evolve(recombination_chance)

        changes = True
        while changes:
            changes = False
            in_degrees = np.array(self.get_in_degrees())
            biconnected_components = self._get_biconnected_components(leafless=True)
            for bc in biconnected_components:
                in_sum_bc = in_degrees[bc]
                bc_reticulations = {bc[node_number] for node_number in np.where(in_sum_bc == 2)[0]}
                if len(bc_reticulations) - self.level > 0:
                    changes=True
                    reticulation_to_remove = bc_reticulations.pop()
                    children_of_reticulation = sorted(list(self._get_children({reticulation_to_remove})), reverse=True)
                    for node_number in children_of_reticulation:
                        self._remove_node_from_network(node_number)
                    break
        self.standardize_leaf_node_names()
        self.standardize_internal_node_names()
        self.prune()

    def evolve(self, recombination_chance: int = 0):
        # Note: does not keep level intact perfectly
        assert 0 <= recombination_chance <= 1
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
                max_level = max(node_levels_dict.keys())
                for level in set(range(max(max_level, self.level))).difference(node_levels_dict.keys()):
                    node_levels_dict[level] = {}
                first_leaf_level = np.where(np.array([int(first_leaf in node_levels_dict[level]) for level in range(max_level + 1)]) == 1)[0][0]
                max_second_leaf_level = self.level - first_leaf_level - 1
                second_leaf = set(itertools.chain.from_iterable([node_levels_dict[level] for level in range(max_second_leaf_level + 1)])).difference(
                    [first_leaf]).pop()
                self._recombine_leaves(first_leaf, second_leaf)
            except KeyError:
                pass
        else:
            leaf_number = random.choice(self.leaf_numbers)
            left_leaf_number, right_leaf_number = self._split_leaf(leaf_number)

        self._biconnected_components = None
        self._partial_ordering = None
        self._cut_arc_matrix = None
        self._cut_arc_sets = None

    def enewick(self):
        return self._enewick(self.get_root_number())

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


class TrinetInfo:
    def __init__(self, trinet, info: dict = None):
        self.trinet = trinet
        self.info = coalesce(info, dict())

    def shrink(self, leaf_set):
        can_shrink = self.trinet.shrink(leaf_set)
        return bool(can_shrink)

    def calculate_info(self):
        self.info['cut_arc_sets'] = self.trinet.cut_arc_sets

    def add_info(self, trinet_info):
        for key, value in trinet_info.info.items():
            self.info[key] = value

    @classmethod
    def limit_to(cls, trinet_info, leaf_names):
        trinet = RootedLevelKNetwork.from_network(trinet_info.trinet, node_names=leaf_names, suppress_parallel=True, suppress_redundant='all')
        result = cls(trinet)
        return result

    def __getitem__(self, item):
        return self.info[item]

    def __setitem__(self, item, value):
        self.info[item] = value

    def __copy__(self):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.trinet = copy.copy(self.trinet)
        cp.info = copy.copy(self.info)
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.trinet = copy.deepcopy(self.trinet)
        cp.info = copy.deepcopy(self.info)
        return cp

    def __str__(self):
        return str(self.trinet.leaf_names) + "\n" + pp.pformat(self.info)


class TrinetInfoList:
    def __init__(self, lst: list = None):
        self.list = coalesce(lst, [])

    def calculate_info(self):
        for trinet_info in self:
            trinet_info.calculate_info()

    def add_info(self, trinet_info_list):
        for trinet_info in self:
            equal_structured_trinets, relation_dicts = trinet_info_list.find_equal_structured_trinet(trinet_info.trinet)
            # TODO check if two?
            try:
                equal_structured_trinet = equal_structured_trinets[0]
                trinet_info.add_info(equal_structured_trinet)
                trinet_info['relation_dict'] = relation_dicts[0]
                trinet_info['leaf_order'] = trinet_info.trinet.get_leaf_ordering()
            except IndexError:
                pass

    def shrink(self, leaf_set):
        to_remove = []
        for trinet_info in self:
            if not trinet_info.shrink(leaf_set):
                to_remove.append(trinet_info)
        for trinet_info in to_remove:
            self.remove(trinet_info)
        # TODO, remove duplicates?

        # to_keep = TrinetInfoList()
        # leaf_iterator = itertools.combinations(self.represented_leaves(), 3)
        # for triplet in leaf_iterator:
        #     trinet_info_list = self.find_trinet_by_leaf_names(triplet)
        #     for trinet_info in trinet_info_list:
        #         for trinet_info_2 in to_keep:
        #             if set(trinet_info.trinet.leaf_names) == set(trinet_info_2.trinet.leaf_names) and trinet_info.trinet.equal_structure(trinet_info_2.trinet):
        #                 break
        #         else:
        #             to_keep.append(trinet_info)
        # self.list = to_keep.list

    def append(self, other):
        assert type(other) == TrinetInfo, "Can only add object of type TrinetInfo to a TrinetInfoList"
        self.list.append(other)

    def extend(self, other):
        self.list.extend(other)

    def find_equal_structured_trinet(self, trinet):
        result = TrinetInfoList()
        relation_dicts = []
        for trinet_info in self:
            are_equal, relation_dict, leaf_translation_dicts = trinet_info.trinet.equal_structure(trinet, optimize=True)
            if are_equal:
                result.append(trinet_info)
                relation_dicts.append(leaf_translation_dicts[0])
        return result, relation_dicts

    def find_equal_trinet(self, trinet):
        for trinet_info in self:
            if trinet_info.trinet == trinet:
                return trinet_info
        return False

    def find_equal_leaf_names_trinet(self, trinet):
        return self.find_trinet_by_leaf_names(trinet.leaf_names)

    def find_trinet_by_leaf_names(self, leaf_names):
        return TrinetInfoList([trinet_info for trinet_info in self.list if set(trinet_info.trinet.leaf_names) == set(leaf_names)])

    def contains_trinet_with_leaf_names(self, leaf_names):
        for trinet_info in self:
            if set(trinet_info.trinet.leaf_names) == set(leaf_names):
                return True
        return False

    def remove_trinet_by_leaf_names(self, leaf_names):
        trinet_infos = self.find_trinet_by_leaf_names(leaf_names)
        for trinet_info in trinet_infos:
            self.remove(trinet_info)
        return trinet_infos

    def find_trinets_with_leaf_names_in(self, leaf_names):
        # assert len(leaf_names) >= 2, "Leaf_names should always include at least two leaves"
        result = TrinetInfoList()
        if len(leaf_names) > 2:
            for trinet_info in self:
                if set(trinet_info.trinet.leaf_names).issubset(leaf_names):
                    result.append(trinet_info)
        else:
            for trinet_info in self:
                if set(leaf_names).issubset(trinet_info.trinet.leaf_names):
                    result.append(TrinetInfo.limit_to(trinet_info, leaf_names))
        return result

    def trinets_where(self, category, value):
        result = TrinetInfoList()
        for trinet_info in self:
            if trinet_info[category] == value:
                result.append(trinet_info)
        return result

    def trinets_with_category(self, category):
        result = TrinetInfoList()
        for trinet_info in self:
            if category in trinet_info.info.keys():
                result.append(trinet_info)
        return result

    def represented_leaves(self):
        result = set()
        for trinet_info in self:
            result.update(trinet_info.trinet.leaf_names)
        return result

    def get_minimal_sink_sets(self, level=0):
        leaf_set = self.represented_leaves()
        arc_count_dict = dict()
        for trinet_info in self:
            cut_arc_sets = trinet_info['cut_arc_sets']
            for cut_arc_set in cut_arc_sets:
                if len(cut_arc_set) == 2:
                    x = cut_arc_set[0]
                    y = cut_arc_set[1]
                    z = [Z for Z in trinet_info.trinet.leaf_names if Z not in cut_arc_set][0]
                    for i in (x, y):
                        try:
                            arc_count_dict[(i, z)] += 1
                        except KeyError:
                            arc_count_dict[(i, z)] = 1

        adjacency_dict = dict()
        arc_iterator = itertools.permutations(leaf_set, 2)
        for arc in arc_iterator:
            try:
                count = arc_count_dict[arc]
            except KeyError:
                count = 0
            if count <= level:
                try:
                    adjacency_dict[arc[0]].add(arc[1])
                except KeyError:
                    adjacency_dict[arc[0]] = {arc[1]}

        strongly_connected_components = tarjan(adjacency_dict)
        msss = []
        for strongly_connected_component in strongly_connected_components:
            for node in strongly_connected_component:
                to_leaves = adjacency_dict[node]
                if not set(to_leaves).issubset(strongly_connected_component):
                    break
            else:
                msss.append(strongly_connected_component)

        return [mss for mss in msss if len(mss) > 1]  # and len(mss) != len(leaf_set)

    def max_level(self):
        level = 0
        for trinet_info in self:
            level = max(level, trinet_info['level'])
        return level

    def best_generator(self):
        generators = []
        generator_count = []
        for trinet_info in self:
            try:
                generator_index = generators.index(trinet_info['generator'])
                generator_count[generator_index] += 1
            except ValueError:
                generators.append(trinet_info['generator'])
                generator_count.append(1)

        max_index = generator_count.index(max(generator_count))
        generator = generators[max_index]
        return copy.deepcopy(generator)

    def get_leaf_locations(self):
        leaf_set = copy.deepcopy(self.represented_leaves())
        # Create: Leaf --> Edge --> (how many times leaf is on this edge)
        leaf_on_edge_count_dict = {leaf: {} for leaf in leaf_set}
        relation_dict_count = {leaf: {'A': 0, 'B': 0, 'C': 0} for leaf in leaf_set}
        for trinet_info in self:
            relation_dict = trinet_info['relation_dict']
            for key, value in relation_dict.items():
                try:
                    relation_dict_count[value][key] += 1
                except KeyError:
                    trinet_info.trinet.visualize()
            for leaf, edge in trinet_info['extra_leaf_dict'].items():
                translated_leaf = relation_dict[leaf]
                try:
                    leaf_on_edge_count_dict[translated_leaf][edge] += 1
                except KeyError:
                    leaf_on_edge_count_dict[translated_leaf][edge] = 1

        # Create: Edge --> (leaves which are on this edge)
        edge_leaf_dict = dict()
        for leaf, edge_count_dict in leaf_on_edge_count_dict.items():
            try:
                max_count = max(edge_count_dict.values())
            except ValueError:
                continue
            best_edges = [edge for edge, count in edge_count_dict.items() if count == max_count]
            best_edge = best_edges[0]

            try:
                edge_leaf_dict[best_edge].add(leaf)
            except KeyError:
                edge_leaf_dict[best_edge] = {leaf}
            leaf_set.remove(leaf)
            relation_dict_count.pop(leaf)

        # TODO: Need to keep n (=number of reticulations in generator) in relation_dict
        # Choose the ones which are the least number of times on a side

        # Find which reticulation is which
        reticulation_relation_dict = {reticulation: set() for reticulation in relation_dict_count}
        for reticulation, counts in relation_dict_count.items():
            m = max(counts.values())
            for option, count in counts.items():
                if count == m:
                    reticulation_relation_dict[reticulation].add(option)

        # TODO: Make sure that not the same option is chosen for different reticulations
        reticulation_relation_dict = {reticulation: options.pop() for reticulation, options in reticulation_relation_dict.items()}

        return edge_leaf_dict, reticulation_relation_dict

    def leaf_order_info(self):
        leaf_set = bidict({leaf: index for leaf, index in enumerate(self.represented_leaves())})
        leaf_order_matrix = np.zeros((len(leaf_set), len(leaf_set)))
        for trinet_info in self:
            relation_dict = trinet_info['relation_dict']

            leaf_order = trinet_info['leaf_order']
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

    def __getitem__(self, item):
        return self.list[item]

    def __iter__(self):
        return iter(self.list)

    def __len__(self):
        return len(self.list)
