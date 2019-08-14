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

pp = pprint.PrettyPrinter(indent=4)

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'


class RootedLevelKNetwork:
    def __init__(self, adj_matrix: np.ndarray, node_name_map: bidict, leaf_names: list, level: int, dimension: int):
        self.uid = guid()
        self.logger = logging.getLogger('network.{}'.format(self.uid))
        self.logger.debug("Creating new rooted level {} network of dimension {}.".format(level, dimension))
        shape = adj_matrix.shape
        assert len(node_name_map) == shape[0], "Number of names does not match number of nodes."
        assert shape[0] == shape[1], "Adjacency matrix is not square."

        self.adj_matrix = adj_matrix
        self.node_name_map = node_name_map
        self._leaf_names = leaf_names

        # Set shape numbers
        self.number_of_leaves = len(self._leaf_names)
        self.number_of_nodes = shape[1]
        self.number_of_internals = self.number_of_nodes - self.number_of_leaves

        # Set network properties
        self.level = level
        self.dimension = dimension

        # optimization variables
        self.biconnected_components = None
        self.partial_ordering = None
        self.cut_arc_matrix = None
        self.cut_arc_sets = None

    # ---------------------------------------------------------------- CLASS METHODS ---------------------------------------------------------------- #
    def __copy__(self):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.adj_matrix = copy.copy(self.adj_matrix)
        cp.node_name_map = copy.copy(self.node_name_map)
        cp.leaf_names = copy.copy(self.leaf_names)
        cp.number_of_leaves = copy.copy(self.number_of_leaves)
        cp.number_of_nodes = copy.copy(self.number_of_nodes)
        cp.number_of_internals = copy.copy(self.number_of_internals)
        cp.level = copy.copy(self.level)
        cp.dimension = copy.copy(self.dimension)
        cp.uid = guid()
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp.biconnected_components = copy.copy(self.biconnected_components)
        cp.partial_ordering = copy.copy(self.partial_ordering)
        cp.cut_arc_matrix = copy.copy(self.cut_arc_matrix)
        cp.cut_arc_sets = copy.copy(self.cut_arc_sets)
        return cp

    def __deepcopy__(self, memo):
        cls = self.__class__
        cp = cls.__new__(cls)
        cp.adj_matrix = copy.deepcopy(self.adj_matrix)
        cp.node_name_map = copy.deepcopy(self.node_name_map)
        cp.leaf_names = copy.deepcopy(self.leaf_names)
        cp.number_of_leaves = copy.deepcopy(self.number_of_leaves)
        cp.number_of_nodes = copy.deepcopy(self.number_of_nodes)
        cp.number_of_internals = copy.deepcopy(self.number_of_internals)
        cp.level = copy.deepcopy(self.level)
        cp.dimension = copy.deepcopy(self.dimension)
        cp.uid = guid()
        cp.logger = logging.getLogger('network.{}'.format(self.uid))
        cp.biconnected_components = copy.deepcopy(self.biconnected_components)
        cp.partial_ordering = copy.deepcopy(self.partial_ordering)
        cp.cut_arc_matrix = copy.deepcopy(self.cut_arc_matrix)
        cp.cut_arc_sets = copy.deepcopy(self.cut_arc_sets)
        return cp

    @classmethod
    def trinet_from_network(cls, network, taxa: list):
        """Create trinet from network."""
        logging.debug("Creating trinet from network {} using {}.".format(network.uid, taxa))
        network.logger.debug("Creating trinet using {}.".format(taxa))
        assert len(taxa) <= 3, "Can not create trinet from network {} using as this {} are more than 3 leaves.".format(network.uid, taxa)
        trinet = cls.from_network(network, taxa, suppress_redundant='all', suppress_parallel=True)
        trinet.logger.debug("Trinet of network {} with taxa {}.".format(network.uid, taxa))
        return trinet

    @classmethod
    def from_network(cls, network, taxa: list, suppress_redundant='none', suppress_parallel=False):
        """Create sub-network from network."""
        logging.debug("Creating sub=network from {} using {}.".format(network.uid, taxa))
        network.logger.debug("Creating sub-network using {}.".format(taxa))
        assert set(taxa).issubset(set(network.leaf_names)), "Can not create sub-network of network {} using {} as they are not leaves of network.".format(
            network.uid,
            taxa)
        adj_matrix = copy.deepcopy(network.adj_matrix)
        node_name_map = copy.deepcopy(network.node_name_map)
        level = copy.copy(network.level)
        dimension = copy.copy(network.dimension)
        new_network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_names=taxa, level=level, dimension=dimension)
        new_network.prune(suppress_redundant=suppress_redundant, suppress_parallel=suppress_parallel)
        new_network.logger.debug("Sub-network of network {} with taxa {}.".format(network.uid, taxa))
        return new_network

    @classmethod
    def from_dir_adj_matrix(cls, dir_adj_matrix: np.ndarray, level=2, dimension=2):
        """Create network from directed adjacency matrix. Assumes only internal nodes have rows. Other nodes are leaves."""
        logging.debug("Creating network from directed adjacency matrix.")
        shape = dir_adj_matrix.shape
        adj_matrix = np.zeros((shape[1], shape[1])).astype(int)
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
        assert network.is_valid(), "Connection dictionary results in a invalid network."
        network.logger.debug("Created from directed adjacency matrix.")
        return network

    @classmethod
    def from_connections_dict(cls, connections_dict: dict, level=2, dimension=2):
        """Create network from connection dictionary. Assumes nodes without outgoing arcs are leaves."""
        logging.debug("Creating network from connections dictionary.")
        node_name_map = bidict()
        leaf_names = list()

        # Empty adjacency matrix
        adj_matrix = np.zeros((0, 0)).astype(int)
        network = cls(adj_matrix=adj_matrix, node_name_map=node_name_map, leaf_names=leaf_names, level=level, dimension=dimension)

        network.logger.debug("Adding internal nodes to network.")
        keys = sorted(connections_dict.keys())
        keys.sort()
        for key in keys:
            network._add_internal_node_to_network(str(key))

        network.logger.debug("Adding leaf nodes to network.")
        for from_node in connections_dict.keys():
            for to_node in connections_dict[from_node]:
                if str(to_node) not in network.node_name_map:
                    network._add_leaf_to_network(to_node)

        network.logger.debug("Adding connections between nodes.")
        for from_node in connections_dict.keys():
            for to_node in connections_dict[from_node]:
                network._add_connection(str(from_node), str(to_node))
        assert network.is_valid(), "Connection dictionary results in an invalid network."
        network.logger.debug("Created from connections dictionary.")
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
        biconnected_components = self.get_biconnected_components(leafless=True)
        for bc in biconnected_components:
            node_numbers_bc = self.get_node_numbers_of(bc)
            in_sum_bc = in_degrees[node_numbers_bc]
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
    def leaf_numbers(self):
        """Retrieve numbers corresponding to leaves."""
        return self.get_node_numbers_of(self._leaf_names)

    @property
    def node_names(self) -> list:
        """Retrieve list of node names ordered by their number"""
        return self.get_node_names_of(sorted(list(self.node_name_map.values())))

    @property
    def leaf_names(self) -> list:
        """Retrieve list of leaf names ordered by their number"""
        return self.get_node_names_of(sorted(self.leaf_numbers))

    @leaf_names.setter
    def leaf_names(self, leaf_names: list):
        assert set(leaf_names).issubset(self.node_name_map), "Not all leaf names are in the network."
        self._leaf_names = list(leaf_names)

    def is_connected_node(self, node_name: str) -> bool:
        """Check whether node is connected to the graph."""
        return sum(self.adj_matrix[self.get_node_number_of(node_name), :] != 0) > 0

    # --------------------------------------------------------- ALTER NETWORK METHODS --------------------------------------------------------- #
    def standardize_internal_node_names(self) -> None:
        self.logger.debug(f"Standardizing internal node names")
        # Temporary names
        for internal_node in self.get_node_names_of(list(range(self.number_of_internals))):
            self.rename_node(internal_node, internal_node + "^")

        for index, internal_node in enumerate(self.get_node_names_of(list(range(self.number_of_internals)))):
            self.rename_node(internal_node, str(index))

    def rename_node(self, old_name: str, new_name: str) -> None:
        """Replace name of node with name old_name by new_name."""
        assert new_name not in self.node_names, f"There already exists a node with node name {new_name}."

        translation_dict = {node_name: node_name for node_name in self.node_name_map.keys()}
        translation_dict[old_name] = new_name

        # Rename node in node_names
        node_number = self.node_name_map.pop(old_name)
        self.node_name_map.put(new_name, node_number)

        # Rename node in leaf_names if node is a leaf
        try:
            self._leaf_names.remove(old_name)
            self._leaf_names.append(new_name)
        except ValueError:
            pass

        self.biconnected_components = None
        self.partial_ordering = None
        self.cut_arc_sets = None

    def shrink(self, leaf_set):
        """Shrink set of leaves"""
        intersection = set(self.leaf_names).intersection(leaf_set)
        if len(intersection) == len(self.leaf_names):
            return False
        if len(intersection) == 0:
            return True
        leaf_to_keep = intersection.pop()
        for leaf in intersection:
            self._remove_leaf_status_node(leaf)
        self.prune()
        mss_name = mss_leaf_name(leaf_set)
        self.rename_node(leaf_to_keep, mss_name)

        return mss_name

    def _add_connection(self, from_node_name: str, to_node_name: str):
        """Add connection between from_node_name and to_node_name."""
        self.logger.debug(f"Adding connection between {from_node_name} and {to_node_name}.")

        from_node_number = self.get_node_number_of(from_node_name)
        to_node_number = self.get_node_number_of(to_node_name)

        self.adj_matrix[from_node_number][to_node_number] += 1
        self.adj_matrix[to_node_number][from_node_number] -= 1

        self.biconnected_components = None
        self.partial_ordering = None
        self.cut_arc_matrix = None
        self.cut_arc_sets = None

    def _remove_connection(self, from_node_name: str, to_node_name: str) -> bool:
        """Remove connection between from_node_name and to_node_name."""
        self.logger.debug(f"Removing connection between {from_node_name} and {to_node_name}.")

        from_node_number = self.get_node_number_of(from_node_name)
        to_node_number = self.get_node_number_of(to_node_name)

        assert self.adj_matrix[from_node_number][to_node_number] > 0 > self.adj_matrix[to_node_number][
            from_node_number], "Cannot remove non-existent connection"

        self.adj_matrix[from_node_number][to_node_number] -= 1
        self.adj_matrix[to_node_number][from_node_number] += 1

        self.biconnected_components = None
        self.partial_ordering = None
        self.cut_arc_matrix = None
        self.cut_arc_sets = None

    def _remove_node_name_from_dict(self, node_name: str) -> int:
        """Remove node with name node_name from node_names dictionary."""
        # Remove node name from node_names
        node_number = self.node_name_map.pop(node_name)

        # Decrease node number of nodes with number higher than node_name's
        for y in range(node_number + 1, self.number_of_nodes):
            self.node_name_map[self.node_name_map.inverse[y]] -= 1

        # Remove node name from leaf_names if node is a leaf
        try:
            self._leaf_names.remove(node_name)
        except ValueError:
            pass

        return int(node_number)

    def _remove_node_name_from_adj_matrix(self, node_name: str):
        """Remove node from adj_matrix."""
        node_number = self.get_node_number_of(node_name)
        mask = np.ones(self.number_of_nodes, dtype=bool)
        mask[node_number] = False
        self.adj_matrix = self.adj_matrix[mask, :][:, mask]

    def _remove_node_name_from_network(self, node_name: str, force: bool = False):
        """Remove node fully from network."""
        self.logger.debug(f"Removing node {node_name} from network.")

        is_connected = force or self.is_connected_node(node_name)
        assert is_connected, f"Can not remove node {node_name} as it is still connected"
        node_is_leaf = self.is_leaf_node(node_name)
        self._remove_node_name_from_adj_matrix(node_name)
        self._remove_node_name_from_dict(node_name)
        self.number_of_nodes -= 1
        if node_is_leaf:
            self.number_of_leaves -= 1
        else:
            self.number_of_internals -= 1

        self.biconnected_components = None
        self.partial_ordering = None
        self.cut_arc_matrix = None
        self.cut_arc_sets = None

    def _remove_component(self, component: list):
        """Remove all nodes in component."""
        self.logger.debug(f"Removing nodes in {component} from network.")

        for node_name in component:
            self._remove_node_name_from_network(node_name=node_name, force=True)

    def _add_internal_node_to_dict(self, node_name: str = None) -> str:
        """Add internal node to node_names dictionary. Returns its number."""
        self.logger.debug(f"Adding internal node {node_name} to dictionary.")

        # Get next number, and get a node_name if no name is given
        new_number = self.number_of_internals
        node_name = str(self.number_of_internals) if node_name is None else node_name

        # Increase node number of leaves by 1, as internal node comes before leaves
        for y in range(self.number_of_nodes - 1, new_number - 1, -1):
            self.node_name_map[self.node_name_map.inverse[y]] += 1
        self.node_name_map.put(node_name, new_number)
        return node_name

    def _add_internal_node_to_adj_matrix(self):
        """Add internal node to adj matrix """
        self.logger.debug("Adding internal node to adjacency matrix.")

        v1 = self.adj_matrix[:self.number_of_internals]
        v2 = np.zeros((1, self.number_of_nodes))
        v3 = self.adj_matrix[self.number_of_internals:]
        self.adj_matrix = np.concatenate((v1, v2, v3), axis=0)

        h1 = self.adj_matrix[:, :self.number_of_internals]
        h2 = np.zeros((self.number_of_nodes + 1, 1))
        h3 = self.adj_matrix[:, self.number_of_internals:]
        self.adj_matrix = np.concatenate((h1, h2, h3), axis=1)

    def _add_internal_node_to_network(self, node_name: str = None) -> str:
        """Add internal node to network."""
        self.logger.debug("Adding internal node to network.")

        node_name = self._add_internal_node_to_dict(node_name)
        self.logger.debug(f"Internal node name is {node_name}.")
        self._add_internal_node_to_adj_matrix()

        # Update numbers
        self.number_of_internals += 1
        self.number_of_nodes += 1
        return node_name

    def _suppress_node(self, node_name: str):
        """Remove node and connect its parent to its child if they both exist."""
        self._suppress_component([node_name])

    def _suppress_component(self, component: list):
        """Remove all nodes in component and connect its parent to its child if they both exist."""
        self.logger.debug(f"Suppressing component {component}.")

        in_nodes = self.get_in_nodes_component(component)
        out_nodes = self.get_out_nodes_component(component)
        assert len(out_nodes) == 0 or len(out_nodes) <= 1 and len(
            in_nodes) <= 1, f"Can not suppress node {component} as it does not have at most one in and at most one out node"
        self._remove_component(component)

        if len(in_nodes) == 1 and len(out_nodes) == 1:
            out_node = out_nodes[0]
            in_node = in_nodes[0]
            self._add_connection(in_node, out_node)

    def _remove_leaf_status_node(self, node_name: str):
        """Remove leaf with name node_name"""
        self.logger.debug(f"Changing node {node_name} from leaf status to non-leaf status.")

        assert self.is_leaf_node(node_name), f"Can not remove leaf status of node {node_name} as it is not a leaf."
        self._leaf_names.remove(node_name)
        self.number_of_leaves -= 1
        self.number_of_internals += 1

    def _add_leaf_to_network(self, leaf_name: str = None) -> str:
        """Add leaf name with name leaf_name or the first available name to network."""
        self.logger.debug(f"Adding node {leaf_name} as leaf to network.")

        leaf_name = self._add_leaf_name_to_dict(leaf_name)
        self._add_leaf_to_adj_matrix()

        # Update numbers
        self.number_of_leaves += 1
        self.number_of_nodes += 1
        return leaf_name

    def _add_leaf_to_adj_matrix(self):
        """Increase the size of the adj_matrix by one at the end."""
        self.logger.debug("Increasing size of adjacency matrix by one for leaf.")

        v1 = self.adj_matrix
        v2 = np.zeros((1, self.number_of_nodes))
        self.adj_matrix = np.concatenate((v1, v2), axis=0)

        h1 = self.adj_matrix
        h2 = np.zeros((self.number_of_nodes + 1, 1))
        self.adj_matrix = np.concatenate((h1, h2), axis=1)

    def _add_leaf_name_to_dict(self, leaf_name: str) -> str:
        """Add leaf_name or first available name to node_names dict."""
        self.logger.debug(f"Adding leaf {leaf_name} to dictionary.")

        leaf_name = leaf_name if leaf_name else self.first_unused_leaf_name()
        new_number = self.number_of_nodes
        self.node_name_map.put(leaf_name, new_number)
        self._leaf_names.append(leaf_name)
        return leaf_name

    def add_leaf_to_edge(self, edge, leaf_name: str = None) -> (str, str):
        """Add node between parent and child and attach leaf with leaf_name to it."""
        parent = edge[0]
        child = edge[1]

        internal_name = self._add_internal_node_to_network()
        leaf_name = self._add_leaf_to_network(leaf_name)

        self._remove_connection(parent, child)

        self._add_connection(parent, internal_name)
        self._add_connection(internal_name, child)
        self._add_connection(internal_name, leaf_name)
        return internal_name, leaf_name

    def add_recombination_to_edges(self, edge):
        # TODO
        pass

    def split_leaf(self, leaf_to_split: str):
        self.logger.debug(f"Replacing leaf {leaf_to_split} with cherry.")

        assert self.is_leaf_node(leaf_to_split), f"Cannot replace leaf {leaf_to_split} with cherry as it is not a leaf."
        self._remove_leaf_status_node(leaf_to_split)
        self._add_leaf_to_network(leaf_name=leaf_to_split + "-I")
        self._add_leaf_to_network(leaf_name=leaf_to_split + "-II")
        self._add_connection(leaf_to_split, leaf_to_split + "-I")
        self._add_connection(leaf_to_split, leaf_to_split + "-II")

    def replace_leaf_with_network(self, leaf_name_to_replace: str, replacement_network):
        self.logger.debug(f"Replacing leaf {leaf_name_to_replace} with network {replacement_network}.")

        assert self.is_leaf_node(leaf_name_to_replace), f"Cannot replace leaf {leaf_name_to_replace} with network as it is not a leaf."
        assert set(self.leaf_names).isdisjoint(
            replacement_network.leaf_names), f"Cannot replace leaf {leaf_name_to_replace} with network {replacement_network} as {replacement_network} has some leafs same as {self}"
        replacement_network = copy.deepcopy(replacement_network)
        replacement_network_internal_nodes = set(replacement_network.node_names).difference(replacement_network.leaf_names)
        replacement_network_root = replacement_network.get_root_name() + "*"

        # Add internal nodes from replacement network to current network with primed names
        for node_name in replacement_network_internal_nodes:
            self._add_internal_node_to_network(node_name + "*")

        # Add leaves from replacement network to current network
        for leaf_name in replacement_network.leaf_names:
            self._add_leaf_to_network(leaf_name)

        # Replace leaf_name with root of replacement network
        parent_of_leaf = self.get_in_nodes_node(leaf_name_to_replace)[0]
        self._remove_node_name_from_network(leaf_name_to_replace)
        self._add_connection(parent_of_leaf, replacement_network_root)

        # Add all connections from replacement network to current network
        for internal_node in replacement_network_internal_nodes:
            to_nodes = replacement_network.get_out_nodes_node(internal_node)
            for to_node in to_nodes:
                if to_node in replacement_network_internal_nodes:
                    self._add_connection(internal_node + "*", to_node + "*")
                else:
                    self._add_connection(internal_node + "*", to_node)

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
        return self._leaf_names

    def get_in_degree_node(self, node_name: str, leafless: bool = False) -> int:
        """Retrieve number of arcs entering node_name."""
        self.logger.debug("Getting {}in degree of node {}.".format("leafless " if leafless else "", node_name))
        node_number = self.get_node_number_of(node_name)
        leafless_indicator = self.number_of_internals if leafless else self.number_of_nodes
        mask = self.adj_matrix[:leafless_indicator, node_number] > 0
        in_degree = sum(self.adj_matrix[:leafless_indicator, node_number][mask])
        return int(in_degree)

    def get_out_degree_node(self, node_name: str, leafless: bool = False) -> int:
        """Retrieve number of arcs exiting node_name."""
        self.logger.debug("Getting {}out degree of node {}.".format("leafless " if leafless else "", node_name))
        node_number = self.get_node_number_of(node_name)
        leafless_indicator = self.number_of_internals if leafless else self.number_of_nodes
        mask = self.adj_matrix[:leafless_indicator, node_number] < 0
        out_degree = abs(sum(self.adj_matrix[:leafless_indicator, node_number][mask]))
        return int(out_degree)

    def get_in_degrees(self, leafless: bool = False) -> list:
        """Retrieve number of arcs entering each node."""
        self.logger.debug("Getting {}in degree of all node.".format("leafless " if leafless else ""))
        leafless_indicator = self.number_of_internals if leafless else self.number_of_nodes
        mask = self.adj_matrix[:leafless_indicator, :leafless_indicator] > 0
        in_degree = sum(self.adj_matrix[:leafless_indicator, :leafless_indicator] * mask)
        try:
            return list(in_degree)
        except TypeError:
            return [in_degree]

    def get_out_degrees(self, leafless: bool = False) -> list:
        """Retrieve number of arcs entering each node."""
        self.logger.debug("Getting {}out degree of all nodes.".format("leafless " if leafless else ""))
        leafless_indicator = self.number_of_internals if leafless else self.number_of_nodes
        mask = self.adj_matrix[:leafless_indicator, :leafless_indicator] < 0
        out_degree = -sum(self.adj_matrix[:leafless_indicator, :leafless_indicator] * mask)
        try:
            return list(out_degree)
        except TypeError:
            return [out_degree]

    def get_in_nodes_node(self, node_name: str) -> list:
        """Retrieve all nodes which enter node node_name."""
        self.logger.debug("Computing nodes entering node {}.".format(node_name))

        node_number = self.get_node_number_of(node_name)
        in_nodes = []
        for i in (1, 2):
            in_nodes += list(np.where(self.adj_matrix[node_number, :] == -i)[0]) * (i)
        return self.get_node_names_of(in_nodes)

    def get_out_nodes_node(self, node_name: str) -> list:
        """Retrieve all nodes which exit node node_name."""
        self.logger.debug("Computing nodes exiting node {}.".format(node_name))
        node_number = self.get_node_number_of(node_name)
        out_nodes = list(np.where(self.adj_matrix[node_number, :] == 1)[0]) + list(np.where(self.adj_matrix[node_number, :] == 2)[0])*2
        out_nodes_names = self.get_node_names_of(out_nodes)
        return out_nodes_names

    def get_in_nodes_component(self, component: list) -> list:
        """Retrieve all nodes which enter component."""
        self.logger.debug("Computing nodes entering component {}.".format(component))
        in_nodes = np.array([])
        for node_name in component:
            in_nodes = np.array(list(itertools.chain(in_nodes, self.get_in_nodes_node(node_name))))
        return [node for node in in_nodes if node not in component]

    def get_out_nodes_component(self, component: list) -> list:
        """Retrieve all nodes which exit component."""
        self.logger.debug("Computing nodes exiting component {}.".format(component))
        out_nodes = []
        for node_name in component:
            out_nodes = list(itertools.chain(out_nodes, self.get_out_nodes_node(node_name)))
        return [node for node in out_nodes if node not in component]

    def get_in_out_degree_node(self, node_name: str, leafless: bool = False):
        self.logger.debug("Getting {}in and out degree of node {}.".format("leafless " if leafless else "", node_name))
        return self.get_in_degree_node(node_name, leafless), self.get_out_degree_node(node_name, leafless)

    def get_connections_node(self, node_name: str) -> list:
        """Retrieve all nodes which enter and exit node."""
        self.logger.debug("Computing nodes entering and exiting {}.".format(node_name))
        return self.get_in_nodes_node(node_name) + self.get_out_nodes_node(node_name)

    def get_leaf_children(self, parents: set, max_depth: int = -1) -> set:
        """Retrieve all leaves max_depth below parents."""
        self.logger.debug("Getting leaf children of {}{}.".format(parents, " up to depth {}".format(max_depth) if max_depth >= 0 else ""))
        children = self.get_children(parents, max_depth)
        leaves = set([child for child in children if self.is_leaf_node(child)])
        return leaves

    def get_children(self, parents: set, max_depth: int = -1) -> set:
        """Retrieve all nodes max_depth below parents."""
        self.logger.debug("Getting children of {}{}.".format(parents, " up to depth {}".format(max_depth) if max_depth >= 0 else ""))
        if max_depth == 0:
            return parents
        next_generation = set()
        for parent in parents:
            next_generation.update(set(self.get_out_nodes_node(parent)))
        if next_generation:
            parents.update(self.get_children(next_generation, max_depth - 1))
            return parents
        return parents

    def get_parents(self, children: set, max_height: int = -1) -> set:
        self.logger.debug("Getting parents of {}{}.".format(children, " up to depth {}".format(max_height) if max_height >= 0 else ""))
        """Retrieve all nodes max_depth above children."""
        if max_height == 0:
            return children
        prev_generation = set()
        for child in children:
            prev_generation.update(set(self.get_in_nodes_node(child)))
        if prev_generation:
            children.update(self.get_parents(prev_generation, max_height - 1))
            return children
        return children

    def get_connected_components(self, leafless: bool = False) -> list:
        """Retrieve the connected components (excluding leaves) of the network."""
        self.logger.debug("Computing connected components.")
        # Finds set of connected components
        # Start with a node
        # Iteratively checks which nodes it is connected to through 'unchecked_connections'
        unchecked_nodes = set(self.node_name_map)
        components = []
        while unchecked_nodes:
            node_name = unchecked_nodes.pop()
            component = [node_name]
            unchecked_connections = set(self.get_connections_node(node_name)).intersection(unchecked_nodes)
            while unchecked_connections:
                connection_name = unchecked_connections.pop()
                unchecked_nodes.remove(connection_name)
                new_connections = set(self.get_connections_node(connection_name))
                unchecked_connections = (unchecked_connections.union(new_connections)).intersection(unchecked_nodes)
                component.append(connection_name)
            if not (leafless and len(component) == 1 and self.is_leaf_node(component[0])):
                components.append(component)
        self.logger.debug("Connected components are {}.".format(components))
        return components

    def get_biconnected_components(self, leafless: bool = False) -> list:
        """Compute biconnected components."""
        if self.biconnected_components is None:
            cut_arc_matrix = self.get_cut_arc_matrix()
            self.adj_matrix -= cut_arc_matrix - cut_arc_matrix.T
            self.biconnected_components = self.get_connected_components(leafless=False)
            self.adj_matrix += cut_arc_matrix - cut_arc_matrix.T

        if leafless:
            return [component for component in self.biconnected_components if not (len(component) == 1 and self.is_leaf_node(component[0]))]
        else:
            return self.biconnected_components

    def is_biconnected(self, leafless: bool = False) -> bool:
        """Check if network is biconnected"""
        return len(self.get_biconnected_components(leafless=leafless)) == 1

    def get_cut_arc_sets(self) -> list:
        """Retrieve cut-arc sets of network."""
        if self.cut_arc_sets is None:
            cut_arc_matrix = self.get_cut_arc_matrix()
            _, to_nodes = np.where(cut_arc_matrix == 1)
            self.cut_arc_sets = []
            for to_node in to_nodes:
                ca_set = list(self.get_leaf_children({self.get_node_name_of(to_node)}))
                ca_set.sort()
                self.cut_arc_sets.append(ca_set)
        return self.cut_arc_sets

    def get_cut_arc_matrix(self) -> np.ndarray:
        """Compute indicator matrix for arcs which are cut-arcs."""
        visited = [False] * self.number_of_nodes
        disc = [-1] * self.number_of_nodes
        low = [-1] * self.number_of_nodes
        parent = [None] * self.number_of_nodes
        cut_arc_matrix = np.zeros((self.number_of_nodes, self.number_of_nodes))

        for i in range(self.number_of_nodes):
            if not visited[i]:
                self.cut_arc_helper(i, visited, disc, low, parent, cut_arc_matrix)

        return cut_arc_matrix

    def cut_arc_helper(self, u_number, visited, disc, low, parent, cut_arc_matrix, t=0):
        # TODO source: https://www.geeksforgeeks.org/bridge-in-a-graph/
        visited[u_number] = True
        disc[u_number] = t
        low[u_number] = t
        t += 1
        u_name = self.get_node_name_of(u_number)

        for v_name in self.get_connections_node(u_name):
            v_number = self.get_node_number_of(v_name)
            if not visited[v_number]:
                parent[v_number] = u_number
                self.cut_arc_helper(v_number, visited, disc, low, parent, cut_arc_matrix, t)

                low[u_number] = min(low[u_number], low[v_number])

                if low[v_number] > disc[u_number]:
                    cut_arc_matrix[u_number][v_number] = 1
            elif v_number != parent[u_number]:
                low[u_number] = min(low[u_number], disc[v_number])

    def is_leaf_node(self, node_name: str) -> bool:
        """Check if node_name is leaf."""
        self.logger.debug("Checking if {} is a leaf.".format(node_name))
        return node_name in self._leaf_names

    def is_reticulation_node(self, node_name: str) -> bool:
        """Check if node_name is reticulation."""
        self.logger.debug("Checking if {} is a reticulation.".format(node_name))
        return self.get_in_degree_node(node_name) > 1

    def is_root_node(self, node_name: str) -> bool:
        """Check if node_name is root."""
        self.logger.debug("Checking if {} is a root.".format(node_name))
        return self.get_in_degree_node(node_name) == 0

    def is_structure_node(self, node_name: str) -> bool:
        """Check if node_name is necessary for the structure: not directly connected to any leaves."""
        self.logger.debug("Checking if {} is a structure node.".format(node_name))
        return (not self.is_leaf_node(node_name)) and (
                self.get_out_degree_node(node_name, leafless=False) - self.get_out_degree_node(node_name, leafless=True) == 0)

    def is_connection_node(self, node_name: str) -> bool:
        """Check if node node_name is not part of underlying generator or leaf"""
        self.logger.debug("Checking if {} is a connection node.".format(node_name))
        return self.get_node_type(node_name) == "connection"

    def get_node_type(self, node_name: str) -> str:
        """Find out what type of node node_name is."""
        self.logger.debug("Compute what type of node node {} is.".format(node_name))
        if self.is_leaf_node(node_name):
            return "leaf"
        if self.is_reticulation_node(node_name):
            return "reticulation"
        if self.is_root_node(node_name):
            return "root"
        if self.is_structure_node(node_name):
            return "structure"
        return "connection"

    def is_generator_node(self, node_name: str) -> bool:
        """Check if node_name is part of the underlying generator."""
        self.logger.debug("Checking if {} is a generator node.".format(node_name))
        return self.get_node_type(node_name) in ("structure", "reticulation", "root")

    def get_generator_nodes(self) -> list:
        """Retrieve all generator nodes"""
        self.logger.debug("Retrieving all generator nodes.")
        result = []
        for node in self.node_name_map.keys():
            if self.is_generator_node(node):
                result.append(node)
        return result

    def is_generator(self) -> bool:
        """Check if the network is isomorphic to a generator."""
        self.logger.debug("Checking if network is isomorphic to a generator.")
        for node_name in self.node_name_map.keys():
            if not self.is_connection_node(node_name):
                return False
        return True

    def get_generator_sub_network(self):
        """Retrieve the underlying generator of the network."""
        self.logger.debug("Retrieving the underlying generator of the network.")
        generator_nodes = self.get_generator_nodes()
        generator_nodes_leaves = self.get_leaf_children(set(generator_nodes), 1)
        return RootedLevelKNetwork.from_network(self, list(generator_nodes_leaves), suppress_redundant='none', suppress_parallel=False)

    def get_edges(self, leafless: bool = False) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        self.logger.debug("Retrieving all the {}edges in the network.".format("leafless " if leafless else ""))
        check_size = self.number_of_internals if leafless else self.number_of_nodes
        i = 1
        from_nodes, to_nodes = [], []
        while np.any(self.adj_matrix[:check_size, :check_size] >= i):
            extra_from_nodes, extra_to_nodes = np.where(self.adj_matrix[:check_size, :check_size] >= i)
            from_nodes += list(self.get_node_names_of(extra_from_nodes))
            to_nodes += list(self.get_node_names_of(extra_to_nodes))
            i += 1
        return [(from_node, to_node) for from_node, to_node in zip(np.array(from_nodes), np.array(to_nodes))]

    def prune(self, suppress_redundant: str = 'all', suppress_parallel: bool = True):
        """Suppress al unnecessary/redundant/parallel nodes/edges in network."""
        self.logger.debug("Pruning: Suppressing {} redundant components{}.".format(suppress_redundant, " and parallel edges" if suppress_parallel else ""))

        assert suppress_redundant in ('none', 'strongly', 'all'), "suppress_redundant parameter must be one of none, strongly, all."

        stable_ancestors, lowest_stable_ancestor, nodes_between = self.stable_ancestors(self.leaf_names)

        for node_name in self.node_names:
            if node_name not in nodes_between:
                self._remove_node_name_from_network(node_name, force=True)

        changes = True
        while changes:
            changes = False
            leaf_numbers = set(self.leaf_numbers)
            in_degrees, out_degrees = np.array(self.get_in_degrees()), np.array(self.get_out_degrees())

            # Get non-leaf nodes without children
            childless_nodes = set(np.where(out_degrees == 0)[0])
            childless_nodes.difference_update(leaf_numbers)
            childless_nodes = self.get_node_names_of(list(childless_nodes))

            # Get nodes with no entering arcs and only one exiting arc
            parentless_nodes = set(np.where(in_degrees == 0)[0])
            one_exiting_arc_nodes = set(np.where(out_degrees == 1)[0])
            parentless_one_exiting_arc_nodes = parentless_nodes.intersection(one_exiting_arc_nodes)
            parentless_one_exiting_arc_nodes = self.get_node_names_of(list(parentless_one_exiting_arc_nodes))

            # Get nodes with one in and one out node
            one_entering_arc_nodes = set(np.where(in_degrees == 1)[0])
            one_entering_and_exiting_arc_nodes = one_entering_arc_nodes.intersection(one_exiting_arc_nodes)
            one_entering_and_exiting_arc_nodes = self.get_node_names_of(list(one_entering_and_exiting_arc_nodes))

            nodes_to_remove = list(set(childless_nodes + parentless_one_exiting_arc_nodes + one_entering_and_exiting_arc_nodes))
            for node in nodes_to_remove:
                changes = True
                self._suppress_node(node)

            # Check for parallel arcs
            if suppress_parallel:
                parallel_arcs = np.where(self.adj_matrix >= 2)
                parallel_arcs = zip(self.get_node_names_of(parallel_arcs[0]), self.get_node_names_of(parallel_arcs[1]))
                for from_node, to_node in parallel_arcs:
                    changes = True
                    self._remove_connection(from_node, to_node)

            # Check for redundant components:
            if suppress_redundant != 'none':
                bcs_to_remove = []
                bcs = self.get_biconnected_components(leafless=True)
                for bc in bcs:
                    if suppress_redundant == 'strongly' and self._is_strongly_redundant_component(bc):
                        bcs_to_remove.append(bc)
                    if suppress_redundant == 'all' and self._is_redundant_component(bc):
                        bcs_to_remove.append(bc)
                for bc in bcs_to_remove:
                    changes = True
                    self._suppress_component(bc)

    def _is_redundant_component(self, component: list) -> bool:
        """Check if component is redundant."""
        self.logger.debug("Checking if component {} is redundant.".format(component))
        out_nodes = self.get_out_nodes_component(component)
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
        number_of_reticulations = sum(col_sum > 1)
        return self.number_of_internals, self.number_of_leaves, number_of_reticulations

    def get_exhibited_trinets(self):
        """Retrieve all trinets exhibited by network."""
        self.logger.debug("Retrieving all exhibited trinets.")
        leaves = self.get_node_names_of(list(range(self.number_of_internals, self.number_of_nodes)))
        triplets = itertools.combinations(leaves, 3)
        trinet_info_list = TrinetInfoList()
        for triplet in tqdm(triplets):
            current_trinet = RootedLevelKNetwork.trinet_from_network(self, list(triplet))
            trinet_info = TrinetInfo(current_trinet)
            trinet_info_list.append(trinet_info)
        return trinet_info_list

    def get_partial_ordering(self) -> list:

        root = self.get_root_name()
        new_result = [[root]]
        result = []
        changes = True
        while changes:
            changes = False

            result = copy.deepcopy(new_result)
            new_result = []
            for track in result:
                new_tracks = []
                children = self.get_out_nodes_node(track[-1])
                for child in children:
                    if not self.is_leaf_node(child):
                        changes = True
                        new_tracks.append(track + [child])
                for new_track in new_tracks:
                    new_result.append(new_track)
                if len(new_tracks) == 0:
                    new_result.append(track)
        return result

    def get_leaf_ordering(self):
        partial_node_ordering = self.get_partial_ordering()
        result = []
        for track in partial_node_ordering:
            leaf_track = []
            for node in track:
                children = self.get_out_nodes_node(node)
                if len(children) == 1:
                    if self.is_leaf_node(children[0]):
                        leaf_track.append(children[0])
                        result.append(leaf_track)
                else:
                    count = 0
                    if self.is_leaf_node(children[0]):
                        leaf_track_1 = copy.copy(leaf_track)
                        leaf_track_1.append(children[0])
                        count += 1
                    if self.is_leaf_node(children[1]):
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

    def stable_ancestors(self, node_names):
        leafless_node_names = []
        for node_name in node_names:
            if self.is_leaf_node(node_name):
                parent_of_leaf = self.get_parents({node_name}, max_height=1)
                parent_of_leaf.remove(node_name)
                node_name = parent_of_leaf.pop()
            leafless_node_names.append(node_name)
        partial_ordering = self.get_partial_ordering()
        tracks_to_min_node = []
        tracks_to_max_node = []
        for track in partial_ordering:
            mi = 10 ** 10
            ma = -1
            for node in leafless_node_names:
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
        nodes_between.update(leafless_node_names)
        nodes_between.update(node_names)
        return stable_ancestors_indicis_dict, lowest_stable_ancestor, nodes_between

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
        length = self.number_of_internals if directed else self.number_of_nodes
        result = pd.DataFrame(columns=ordered_node_names, index=ordered_node_names[:length], data=data[:length])
        return result.astype(int)

    def _get_component_list_names(self, component_list: list) -> list:
        """Retrieve names of list of lists of node numbers."""
        self.logger.debug("Retrieve node names of {}.".format(component_list))
        return [self.get_node_name_of(comp) for comp in component_list]

    def equal_structure(self, other):
        self_partial_ordering = self.get_partial_ordering()
        other_partial_ordering = other.get_partial_ordering()

        if len(self_partial_ordering) != len(other_partial_ordering):
            return False, bidict()

        s1 = sum([len(row) for row in self_partial_ordering])
        s2 = sum([len(row) for row in other_partial_ordering])
        if s1 != s2:
            return False, bidict()

        self_names = set().union(*self_partial_ordering)
        other_names = set().union(*other_partial_ordering)
        if len(self_names) != len(other_names):
            return False, bidict()

        other_name_iterator = itertools.permutations(other_names)
        for other_name_iter in other_name_iterator:
            translation_dict = {other_name: self_name for other_name, self_name in zip(other_name_iter, self_names)}
            temp_other_partial_ordering = [[translation_dict[other_name] for other_name in row] for row in other_partial_ordering]
            if self._equal_structure(self_partial_ordering, temp_other_partial_ordering):
                self_leaves = self.leaf_names
                translation_dict = bidict(translation_dict).inverse

                relation_dict = bidict()
                used_leaves = []
                for self_leaf in self_leaves:
                    parents_self_leaf = self.get_parents({self_leaf}, max_height=1)
                    parents_self_leaf.remove(self_leaf)
                    parent_self_leaf = parents_self_leaf.pop()
                    parent_other_leaf = translation_dict[parent_self_leaf]
                    other_leaves = other.get_children({parent_other_leaf}, max_depth=1)
                    other_leaves.remove(parent_other_leaf)
                    chosen_leaf = other_leaves.pop()
                    if chosen_leaf in used_leaves or not other.is_leaf_node(chosen_leaf):
                        chosen_leaf = other_leaves.pop()
                    used_leaves.append(chosen_leaf)
                    relation_dict.put(self_leaf, chosen_leaf)

                return True, relation_dict
        return False, bidict()

    @staticmethod
    def _equal_structure(partial_ordering, other_partial_ordering):
        n = len(partial_ordering)
        it = itertools.permutations(range(n))
        for permutation in it:
            for i, j in enumerate(permutation):
                if partial_ordering[i] != other_partial_ordering[j]:
                    break
            else:
                return True
        return False

    def equal_structure2(self, other):
        if not self.number_of_internals == other.number_of_internals:
            return False, bidict()
        if not self.adj_matrix.shape == other.adj_matrix.shape:
            return False, bidict()

        self_root = self.get_root_name()
        other_root = other.get_root_name()

        self_root_number = self.get_node_number_of(self_root)
        other_root_number = other.get_node_number_of(other_root)

        self_matrix = (self.adj_matrix > 0).astype(int)
        other_matrix = (other.adj_matrix > 0).astype(int)

        equal, relation_dict = self._equal_structure(self_matrix, other_matrix, self_root_number, other_root_number, dict())
        if equal:
            return True, bidict({self.get_node_name_of(key): other.get_node_name_of(value) for key, value in
                                 relation_dict.items() if self.is_leaf_node(self.get_node_name_of(key))})
        else:
            return False, bidict()

    def _equal_structure2(self, matrix_1, matrix_2, root_1, root_2, relation_dict, level=0):
        if sum(matrix_1[root_1]) != sum(matrix_2[root_2]):
            return False, dict()
        if sum(matrix_1[root_1]) == 0:
            relation_dict[root_1] = root_2
            return True, relation_dict

        to_nodes_1 = [index for index, value in enumerate(matrix_1[root_1]) if value == 1]
        to_nodes_2 = [index for index, value in enumerate(matrix_2[root_2]) if value == 1]
        equal = False
        if len(to_nodes_1) == 1:
            to_node_1 = to_nodes_1[0]
            to_node_2 = to_nodes_2[0]
            print(f"Comparing {to_node_1} to {to_node_2}")
            if to_node_1 in relation_dict.keys():
                if relation_dict[to_node_1] != to_node_2:
                    print(f"{to_node_1} is already in relat")
                    return False, dict()
            new_relation_dict = copy.deepcopy(relation_dict)
            new_relation_dict[to_node_1] = to_node_2
            are_equal, returned_relation_dict = self._equal_structure(matrix_1, matrix_2, to_node_1, to_node_2, new_relation_dict, level + 1)
            if are_equal:
                relation_dict = returned_relation_dict
                equal = True
        elif len(to_nodes_1) == 2:
            # to_nodes_1 = (a,b), to_nodes_2 = (x,y)
            # Try: a == x ^ b == y
            # a == x
            to_node_1 = to_nodes_1[0]
            to_node_2 = to_nodes_2[0]
            print(f"Comparing {to_node_1} to {to_node_2}")
            are_equal = False
            if to_node_1 in relation_dict.keys() and relation_dict[to_node_1] == to_node_2:
                new_relation_dict = copy.deepcopy(relation_dict)
                new_relation_dict[to_node_1] = to_node_2
                are_equal, returned_relation_dict = self._equal_structure(matrix_1, matrix_2, to_node_1, to_node_2, new_relation_dict, level + 1)
            if are_equal:
                # b == y
                to_node_1 = to_nodes_1[1]
                to_node_2 = to_nodes_2[1]

                are_equal = False
                if to_node_1 in returned_relation_dict.keys() and returned_relation_dict[to_node_1] == to_node_2:
                    new_relation_dict = copy.deepcopy(returned_relation_dict)
                    new_relation_dict[to_node_1] = to_node_2
                    are_equal, returned_relation_dict = self._equal_structure(matrix_1, matrix_2, to_node_1, to_node_2, new_relation_dict, level + 1)
                if are_equal:
                    relation_dict = returned_relation_dict
                    equal = True
            else:
                # Try: a == y ^ b == x
                # a == y
                to_node_1 = to_nodes_1[0]
                to_node_2 = to_nodes_2[1]
                print(f"Comparing {to_node_1} to {to_node_2}")
                are_equal = False
                if to_node_1 in relation_dict.keys() and relation_dict[to_node_1] == to_node_2:
                    new_relation_dict = copy.deepcopy(relation_dict)
                    new_relation_dict[to_node_1] = to_node_2
                    are_equal, returned_relation_dict = self._equal_structure(matrix_1, matrix_2, to_node_1, to_node_2, new_relation_dict, level + 1)
                if are_equal:
                    # b == x
                    to_node_1 = to_nodes_1[1]
                    to_node_2 = to_nodes_2[0]
                    print(f"Comparing {to_node_1} to {to_node_2}")
                    are_equal = False
                    if to_node_1 in returned_relation_dict.keys() and returned_relation_dict[to_node_1] == to_node_2:
                        new_relation_dict = copy.deepcopy(returned_relation_dict)
                        new_relation_dict[to_node_1] = to_node_2
                        are_equal, returned_relation_dict = self._equal_structure(matrix_1, matrix_2, to_node_1, to_node_2, new_relation_dict, level + 1)
                    if are_equal:
                        relation_dict = returned_relation_dict
                        equal = True

        if equal:
            return True, relation_dict
        else:
            return False, dict()

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
        self.info['cut_arc_sets'] = self.trinet.get_cut_arc_sets()

    def add_info(self, trinet_info):
        for key, value in trinet_info.info.items():
            self.info[key] = value

    @classmethod
    def limit_to(cls, trinet_info, leaf_names):
        trinet = RootedLevelKNetwork.from_network(trinet_info.trinet, leaf_names, suppress_parallel=True, suppress_redundant='all')
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


class TrinetInfoList(list):
    def __init__(self):
        super().__init__()

    def calculate_info(self):
        for trinet_info in self:
            trinet_info.calculate_info()

    def add_info(self, trinet_info_list):
        for trinet_info in self:
            equal_structured_trinets, relation_dicts = trinet_info_list.find_equal_structured_trinet(trinet_info.trinet)
            try:
                equal_structured_trinet = equal_structured_trinets[0]
                trinet_info.add_info(equal_structured_trinet)
                trinet_info['relation_dict'] = relation_dicts[0]
                trinet_info['leaf_order'] = trinet_info.trinet.get_leaf_ordering()
            except IndexError:
                pass

    def shrink(self, leaf_set):
        to_remove = []
        for trinet_info in iter(self):
            if not trinet_info.shrink(leaf_set):
                to_remove.append(trinet_info)
        for trinet_info in to_remove:
            self.remove(trinet_info)
        # TODO, remove duplicates?

    def append(self, other):
        assert type(other) == TrinetInfo, "Can only add object of type TrinetInfo to a TrinetInfoList"
        super().append(other)

    def find_equal_structured_trinet(self, trinet):
        result = TrinetInfoList()
        relation_dicts = []
        for trinet_info in iter(self):
            are_equal, relation_dict = trinet_info.trinet.equal_structure(trinet)
            if are_equal:
                result.append(trinet_info)
                relation_dicts.append(relation_dict)
        return result, relation_dicts

    def find_equal_trinet(self, trinet):
        for trinet_info in iter(self):
            if trinet_info.trinet == trinet:
                return trinet_info
        return False

    def find_equal_leaf_names_trinet(self, trinet):
        return self.find_trinet_by_leaf_names(trinet.leaf_names)

    def find_trinet_by_leaf_names(self, leaf_names):
        result = TrinetInfoList()
        for trinet_info in iter(self):
            if set(trinet_info.trinet.leaf_names) == set(leaf_names):
                result.append(trinet_info)
        return result

    def contains_trinet_with_leaf_names(self, leaf_names):
        for trinet_info in iter(self):
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
            for trinet_info in iter(self):
                if set(trinet_info.trinet.leaf_names).issubset(leaf_names):
                    result.append(trinet_info)
        else:
            for trinet_info in iter(self):
                if set(leaf_names).issubset(trinet_info.trinet.leaf_names):
                    result.append(TrinetInfo.limit_to(trinet_info, leaf_names))
        return result

    def trinets_where(self, category, value):
        result = TrinetInfoList()
        for trinet_info in iter(self):
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
                relation_dict_count[value][key] += 1
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

    def order_leaves(self, edge_leaf_dict, generator):
        symmetrical_nodes = generator.symmetrical_nodes

    def __add__(self, other):
        result = copy.deepcopy(self)
        for oth in other:
            result.append(oth)
        return result
