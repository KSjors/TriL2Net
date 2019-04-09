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
from tqdm import tqdm
from bidict import bidict
from graphviz import Digraph
from utils.help_functions import leaf_name_iterator, guid
import time
import os

logging.basicConfig(level=logging.DEBUG)

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'


class RootedLevelKNetwork:
    def __init__(self, adj_matrix: np.ndarray, node_names: bidict, leaf_names: set, level: int, dimension: int):
        logging.debug("Creating new rooted level {} network of dimension {}.".format(level, dimension))
        shape = adj_matrix.shape
        assert len(node_names) == shape[0], "Number of names does not match number of nodes."
        assert shape[0] == shape[1], "Adjacency matrix is not square."
        assert leaf_names.issubset(node_names)

        self.uid = guid()
        logging.debug("Network has uid {}.".format(self.uid))

        self.adj_matrix = adj_matrix
        self.node_names = node_names
        self.leaf_names = leaf_names

        # Set shape numbers
        self.number_of_leaves = len(self.leaf_names)
        self.number_of_nodes = shape[1]
        self.number_of_internals = self.number_of_nodes - self.number_of_leaves

        # Set network properties
        self.level = level
        self.dimension = dimension

    @classmethod
    def copy_network(cls, network):
        """Copy network."""
        logging.debug("Copying network {}.".format(network.uid))
        adj_matrix = copy.deepcopy(network.adj_matrix)
        node_names = copy.deepcopy(network.node_names)
        leaf_names = copy.deepcopy(network.leaf_names)
        level = copy.copy(network.level)
        dimension = copy.copy(network.dimension)
        network_copy = cls(adj_matrix=adj_matrix, node_names=node_names, leaf_names=leaf_names, level=level, dimension=dimension)
        return network_copy

    @classmethod
    def trinet_from_network(cls, network, taxa: set):
        """Create trinet from network."""
        logging.debug("Creating trinet from network {} using {}.".format(network.uid, taxa))
        assert len(taxa), "Can not create trinet from network {} using as this {} are more than 3 leaves.".format(network.uid, taxa)
        trinet = RootedLevelKNetwork.from_network(network, taxa, suppress_redundant='strongly', suppress_parallel=True)
        # trinet.to_standard_form()
        return trinet

    @classmethod
    def from_network(cls, network, taxa: set, suppress_redundant='none', suppress_parallel=False):
        """Create sub-network from network."""
        logging.debug("Creating sub-network from network {} using {}.".format(network.uid, taxa))
        assert taxa.issubset(network.get_leaves()), "Can not create sub-network of network {} using {} as they are not leaves of network.".format(network.uid,
                                                                                                                                                  taxa)
        adj_matrix = copy.deepcopy(network.adj_matrix)
        node_names = copy.deepcopy(network.node_names)
        level = copy.copy(network.level)
        dimension = copy.copy(network.dimension)
        new_network = cls(adj_matrix=adj_matrix, node_names=node_names, leaf_names=taxa, level=level, dimension=dimension)
        new_network.prune(suppress_redundant=suppress_redundant, suppress_parallel=suppress_parallel)
        return new_network

    @classmethod
    def trinet_from_adj_matrix(cls, dir_adj_matrix: np.ndarray, level=2, dimension=2):
        """Create network from directed adjacency matrix. Assumes only internal nodes have rows. Other nodes are leaves."""
        logging.debug("Creating trinet from directed adjacency matrix.")

    @classmethod
    def from_dir_adj_matrix(cls, dir_adj_matrix: np.ndarray, level=2, dimension=2):
        """Create network from directed adjacency matrix. Assumes only internal nodes have rows. Other nodes are leaves."""
        logging.debug("Creating network from directed adjacency matrix.")
        shape = dir_adj_matrix.shape
        adj_matrix = np.zeros((shape[1], shape[1])).astype(int)
        adj_matrix[:shape[0], :shape[1]] += dir_adj_matrix
        adj_matrix[:shape[1], :shape[0]] -= dir_adj_matrix.T
        node_names = bidict()
        leaf_names = set()
        ln_iterator = leaf_name_iterator(1, 100)
        for i in range(shape[1]):
            if i < shape[0]:
                node_names.put(str(i), i)
            else:
                leaf_name = "".join(next(ln_iterator))
                node_names.put(leaf_name, i)
                leaf_names.add(leaf_name)
        network = cls(adj_matrix=adj_matrix, node_names=node_names, leaf_names=leaf_names, level=level, dimension=dimension)
        assert network.is_valid(), "Connection dictionary results in a invalid network."
        return network

    @classmethod
    def from_connections_dict(cls, connections_dict: dict, level=2, dimension=2):
        """Create network from connection dictionary. Assumes nodes without outgoing arcs are leaves."""
        node_names = bidict()
        leaf_names = set()

        # Fill in internal node names
        count = 0
        keys = sorted(connections_dict.keys())
        keys.sort()
        for key in keys:
            node_names.put(str(key), count)
            count += 1

        # Fill in leaf node names
        for from_node in connections_dict.keys():
            for to_node in connections_dict[from_node]:
                if str(to_node) not in node_names:
                    node_names.put(str(to_node), int(count))
                    leaf_names.add(to_node)
                    count += 1
        number_of_nodes = count

        # Empty adjacency matrix
        adj_matrix = np.zeros((number_of_nodes, number_of_nodes)).astype(int)
        network = cls(adj_matrix=adj_matrix, node_names=node_names, leaf_names=leaf_names, level=level, dimension=dimension)

        # Fill in adjacency matrix
        for from_node in connections_dict.keys():
            for to_node in connections_dict[from_node]:
                network._add_connection(str(from_node), str(to_node))
        assert network.is_valid(), "Connection dictionary results in a invalid network."
        return network

    def is_valid(self) -> bool:
        """Check if network is valid: has right degrees and number of roots and reticulations"""
        logging.debug("Checking if network {} is valid.".format(self.uid))
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

        total_number_of_reticulations = sum(in_degrees > 1)
        if total_number_of_reticulations <= self.level:
            return True
        biconnected_components = self.get_biconnected_components()
        for bc in biconnected_components:
            node_numbers_bc = self._get_node_numbers(bc)
            in_sum_bc = in_degrees[node_numbers_bc]
            number_of_reticulations = sum(in_sum_bc > 1)
            assert number_of_reticulations > self.level, "Network {} is not valid, biconnected component {} has to many reticulations ({} > {}).".format(
                self.uid, bc, number_of_reticulations, self.level)
        return True

    def _to_block_form(self) -> int:
        """Put adjacency matrix in block form. Returns number of generator nodes."""
        logging.debug("Putting adjacency matrix of trinet {} in block form".format(self.uid))
        generator_nodes = self.get_generator_nodes()
        count = 0
        for node in generator_nodes:
            self.set_node_as_number(node, count)
            count += 1
        return count

    def to_standard_form(self):
        logging.debug("Putting adjacency matrix of trinet {} in standard form.".format(self.uid))
        ints, leaves, rets = self.number_of_internals_leaves_reticulations()
        assert 0 <= rets <= 2, "Trinet {} has to many reticulations ({}).".format(self.uid, rets)
        assert leaves == 3, "Trinet does not have exactly three leaves"
        if rets == 0:
            self._to_standard_form_0()
        elif rets == 1:
            self._to_standard_form_1()
        else:
            self._to_standard_form_2()

    def _to_standard_form_0(self) -> list:
        """Put adjacency matrix of level-0 trinet in standard form."""
        logging.debug("Putting adjacency matrix of level-0 trinet {} in standard form".format(self.uid))
        root = self.get_root_name()
        self.set_node_as_number(root, 0)
        root_children = self.get_children({root}, 1).difference({root})
        child_a = root_children.pop()
        child_b = root_children.pop()
        transformations = []
        if not self._is_leaf_node(child_a):
            transformations.append(self.set_node_as_number(child_a, 1))
            transformations.append(self.set_node_as_number(child_b, 2))
        else:
            transformations.append(self.set_node_as_number(child_b, 1))
            transformations.append(self.set_node_as_number(child_a, 2))
        return transformations

    def _to_standard_form_1(self) -> list:
        """Put adjacency matrix of level-1 trinet in standard form."""
        logging.debug("Putting adjacency matrix of level-1 trinet {} in standard form".format(self.uid))

        if not self.is_biconnected():
            return []

        number_of_generator_nodes = self._to_block_form()

        generator = self.get_generator_sub_network()
        transformations = generator.to_standard_form_gen_1()
        for transformation in transformations:
            if len(transformation) == 2:
                self.switch_nodes_in_adj_matrix(transformation[0], transformation[1])
        transformations += self._sort_extra_nodes(number_of_generator_nodes)
        return transformations

    def to_standard_form_gen_1(self) -> list:
        """Put level-1 generator adjacency matrix in standard form."""
        logging.debug("Putting adjacency matrix of level-1 generator {} in standard form".format(self.uid))
        transformations = []
        root = self.get_root_name()
        transformations.append(self.set_node_as_number(root, 0))
        root_children = set(self.get_children({root}, 1).difference({root}))
        assert len(root_children) == 1, "Level-1 generator root has more than one child."

        child_a = root_children.pop()
        transformations.append(self.set_node_as_number(child_a, 1))
        return transformations

    def _to_standard_form_2(self) -> list:
        """Put adjacency matrix of level-2 trinet in standard form."""
        logging.debug("Putting adjacency matrix of level-2 trinet {} in standard form".format(self.uid))
        if not self.is_biconnected():
            return []

        number_of_generator_nodes = self._to_block_form()

        # TODO: make a class function
        generator = self.get_generator_sub_network()
        transformations = generator.to_standard_form_gen_2()

        for transformation in transformations:
            if len(transformation) == 2:
                self.switch_nodes_in_adj_matrix(transformation[0], transformation[1])
        transformations += self._sort_extra_nodes(number_of_generator_nodes)
        return transformations

    def to_standard_form_gen_2(self) -> list:
        """Put level-2 generator adjacency matrix in standard form."""
        logging.debug("Putting adjacency matrix of level-2 generator {} in standard form".format(self.uid))
        transformations = []

        root = self.get_root_name()
        transformations.append(self.set_node_as_number(root, 0))

        root_children = self.get_children({root}, 1).difference({root})
        assert len(root_children) == 2, "Level-2 generator root does not have two children"
        child_a = root_children.pop()
        child_b = root_children.pop()
        transformations.append(self.set_node_as_number(child_a, 1))

        # Find out whose degree is (1, 2)
        if self.get_in_out_degree_node(child_a) == (1, 2):
            child_1 = child_a
            child_2 = child_b
        else:
            child_1 = child_b
            child_2 = child_a

        # Set their numbers
        transformations.append(self.set_node_as_number(child_1, 1))
        transformations.append(self.set_node_as_number(child_2, 2))

        # Set children of 1 other than 2 as 3 (and 4). In case of 3 and 4, child with degree (1, 2) as 3
        children = self.get_children({child_1}, 1).difference({child_1, child_2})
        child_c = children.pop()
        try:
            child_d = root_children.pop()
            # Find out whose degree is (1, 2)
            if self.get_in_out_degree_node(child_c) == (1, 2):
                child_3 = child_c
                child_4 = child_d
            else:
                child_3 = child_d
                child_4 = child_c
            transformations.append(self.set_node_as_number(child_3, 3))
            transformations.append(self.set_node_as_number(child_4, 4))
        except KeyError:
            transformations.append(self.set_node_as_number(child_c, 3))
        return transformations

    def _sort_extra_nodes(self, number_of_generator_nodes: int) -> list:
        """Put extra nodes of trinet {} in order."""
        logging.debug("Putting extra nodes of trinet {} in order.".format(self.uid))
        # TODO: optimize (non-recursive)
        # Ordering extra block based on depth of parent (first) and child (second)
        extra_nodes = self._get_node_names(list(range(number_of_generator_nodes, self.number_of_internals)))
        extra_nodes_copy = copy.deepcopy(extra_nodes)
        parent_number = [np.where(self.adj_matrix[:, i] >= 1)[0][0] for i in
                         range(number_of_generator_nodes, self.number_of_internals)]
        child_number = [np.where(self.adj_matrix[:, i] <= -1)[0][0] for i in
                        range(number_of_generator_nodes, self.number_of_internals)]
        combined = zip(parent_number, child_number, extra_nodes)
        minor_order = sorted(combined, key=lambda x: x[1])
        major_order = sorted(minor_order, key=lambda x: x[0])
        _, _, ordering = zip(*major_order)
        transformations = []
        if extra_nodes_copy != list(ordering):
            count = 0
            for extra_node in ordering:
                transformations.append(self.set_node_as_number(extra_node, number_of_generator_nodes + count))
                count += 1
            transformations += self._sort_extra_nodes(number_of_generator_nodes)
        return transformations

    def _rename_node(self, old_name: str, new_name: str):
        """Replace name of node with name old_name by new_name."""
        logging.debug("Renaming node {} with name {}.".format(old_name, new_name))
        node_number = self.node_names.pop(old_name)
        self.node_names.put(new_name, node_number)
        if old_name in self.leaf_names:
            self.leaf_names.remove(old_name)
            self.leaf_names.add(new_name)

    def _get_node_names(self, node_number_list: list) -> list:
        """Retrieve names corresponding to node numbers."""
        logging.debug("Retrieving node names of {}.".format(node_number_list))
        result = []
        for node_number in node_number_list:
            result.append(self._get_node_name(node_number))
        return result

    def _get_node_name(self, node_number: int) -> str:
        """Retrieve name corresponding to node number"""
        logging.debug("Retrieving node name of {}.".format(node_number))
        return self.node_names.inverse[node_number]

    def _get_node_numbers(self, node_name_list: list) -> list:
        """Retrieve numbers corresponding to node names."""
        logging.debug("Retrieving node number of {}.".format(node_name_list))
        result = []
        for node_name in node_name_list:
            result.append(self._get_node_number(node_name))
        return result

    def _get_node_number(self, node_name: str) -> int:
        """Retrieve number corresponding to node name."""
        logging.debug("Retrieving node number of {}.".format(node_name))
        return int(self.node_names[node_name])

    def _get_leaf_numbers(self):
        leaves = list(self.get_leaves())
        leaf_numbers = self._get_node_numbers(leaves)
        return leaf_numbers

    def _get_leaf_mask(self) -> list:
        """Return indicator array (node == leaf)."""
        logging.debug("Retrieving leaf mask.")
        leaf_numbers = self._get_leaf_numbers()
        mask = [(i in leaf_numbers) for i in range(self.number_of_nodes)]
        return mask

    def _remove_node_name_from_dict(self, node_name: str) -> int:
        """Remove node with name node_name from node_names dictionary."""
        logging.debug("Removing node {} from dictionary.".format(node_name))
        node_number = self.node_names.pop(node_name)
        for y in range(node_number + 1, self.number_of_nodes):
            self.node_names[self.node_names.inverse[y]] -= 1
        if self._is_leaf_node(node_name):
            self.leaf_names.remove(node_name)
        return node_number

    def _add_internal_node_to_dict(self) -> str:
        """Add internal node to node_names dictionary. Returns its number."""
        logging.debug("Adding internal node to dictionary.")
        new_number = self.number_of_internals
        node_name = str(self.number_of_internals)
        for y in range(self.number_of_nodes - 1, new_number - 1, -1):
            self.node_names[self.node_names.inverse[y]] += 1
        self.node_names.put(node_name, new_number)
        return node_name

    def _add_internal_node_to_adj_matrix(self):
        """Add internal node to adj matrix """
        logging.debug("Adding internal node to adjacency matrix.")
        v1 = self.adj_matrix[:self.number_of_internals]
        v2 = np.zeros((1, self.number_of_nodes))
        v3 = self.adj_matrix[self.number_of_internals:]
        self.adj_matrix = np.concatenate((v1, v2, v3), axis=0)

        h1 = self.adj_matrix[:, :self.number_of_internals]
        h2 = np.zeros((self.number_of_nodes + 1, 1))
        h3 = self.adj_matrix[:, self.number_of_internals:]
        self.adj_matrix = np.concatenate((h1, h2, h3), axis=1)

    def _get_ordered_node_names(self) -> list:
        """Retrieve list of node names ordered by their number"""
        logging.debug("Retrieving order node name list.")
        node_numbers = list(self.node_names.inverse)
        node_numbers.sort()
        return self._get_node_names(node_numbers)

    def _add_connection(self, from_node_name: str, to_node_name: str):
        """Add connection between from_node_name and to_node_name."""
        logging.debug("Adding connection between {} and {}.".format(from_node_name, to_node_name))
        from_node_number = self._get_node_number(from_node_name)
        to_node_number = self._get_node_number(to_node_name)

        self.adj_matrix[from_node_number][to_node_number] += 1
        self.adj_matrix[to_node_number][from_node_number] -= 1

    def _remove_connection(self, from_node_name: str, to_node_name: str):
        """Remove connection between from_node_name and to_node_name."""
        logging.debug("Removing connection between {} and {}.".format(from_node_name, to_node_name))
        from_node_number = self._get_node_number(from_node_name)
        to_node_number = self._get_node_number(to_node_name)

        self.adj_matrix[from_node_number][to_node_number] -= 1
        self.adj_matrix[to_node_number][from_node_number] += 1

    def _is_connected_node(self, node_name: str) -> bool:
        """Check whether node is connected to the graph."""
        logging.debug("Checking whether {} is connected to the network.".format(node_name))
        node_number = self._get_node_number(node_name)
        return sum(self.adj_matrix[node_number, :] != 0) > 0

    def _remove_node_name_from_adj_matrix(self, node_name: str):
        """Remove node from adj_matrix."""
        logging.debug("Removing node {}  from adjacency matrix.".format(node_name))
        node_number = self._get_node_number(node_name)
        mask = np.ones(self.number_of_nodes, dtype=bool)
        mask[node_number] = False
        self.adj_matrix = self.adj_matrix[mask, :][:, mask]

    def _remove_node_name_from_network(self, node_name: str, force: bool = False):
        """Remove node fully from network."""
        logging.debug("Removing node {} from network.".format(node_name))
        is_connected = force or self._is_connected_node(node_name)
        assert is_connected, "Can not remove node {} as it is still connected".format(node_name)
        node_is_leaf = self._is_leaf_node(node_name)
        self._remove_node_name_from_adj_matrix(node_name)
        self._remove_node_name_from_dict(node_name)
        self.number_of_nodes -= 1
        if node_is_leaf:
            self.number_of_leaves -= 1
        else:
            self.number_of_internals -= 1

    def _remove_component(self, component: list):
        """Remove all nodes in component."""
        logging.debug("Removing nodes in {} from network.".format(component))
        for node_name in component:
            self._remove_node_name_from_network(node_name, True)

    def _add_internal_node_to_network(self) -> str:
        """Add internal node to network."""
        logging.debug("Adding internal node to network.")
        node_name = self._add_internal_node_to_dict()
        logging.debug("Internal node name is {}.".format(node_name))
        self._add_internal_node_to_adj_matrix()

        # Update numbers
        self.number_of_internals += 1
        self.number_of_nodes += 1
        return node_name

    def _suppress_node(self, node_name: str):
        """Remove node and connect its parent to its child if they both exist."""
        self._suppress_component([node_name])

    def _suppress_component(self, component: list):
        """Remove all nodes in component and connect its parent to its child if the both exist."""
        logging.debug("Suppressing component {}.".format(component))
        in_nodes = self.get_in_nodes_component(component)
        out_nodes = self.get_out_nodes_component(component)
        if component == ['1']:
            print(in_nodes)
            print(out_nodes)
        assert len(out_nodes) == 0 or len(out_nodes) <= 1 and len(
            in_nodes) <= 1, "Can not suppress node {} as it does not have at most one in and at most one out node".format(
            component)

        self._remove_component(component)

        if len(in_nodes) == 1 and len(out_nodes) == 1:
            out_node = out_nodes[0]
            in_node = in_nodes[0]
            self._add_connection(in_node, out_node)

    def _remove_leaf_status_node(self, node_name: str):
        """Remove leaf with name node_name"""
        logging.debug("Changing node {} from leaf status to non-leaf status.".format(node_name))
        assert self._is_leaf_node(node_name), "Can not remove leaf status of node {} as it is not a leaf.".format(node_name)
        self.leaf_names.remove(node_name)

    def _add_leaf_to_network(self, leaf_name: str = None) -> str:
        """Add leaf name with name leaf_name or the first available name to network."""
        logging.debug("Adding node {} as leaf to network.".format(leaf_name if leaf_name else "?"))
        leaf_name = self._add_leaf_name_to_dict(leaf_name)
        self._add_leaf_to_adj_matrix()

        # Update numbers
        self.number_of_leaves += 1
        self.number_of_nodes += 1
        return leaf_name

    def _add_leaf_to_adj_matrix(self):
        """Increase the size of the adj_matrix by one at the end."""
        logging.debug("Increasing size of adjacency matrix by one for leaf.")
        v1 = self.adj_matrix
        v2 = np.zeros((1, self.number_of_nodes))
        self.adj_matrix = np.concatenate((v1, v2), axis=0)

        h1 = self.adj_matrix
        h2 = np.zeros((self.number_of_nodes + 1, 1))
        self.adj_matrix = np.concatenate((h1, h2), axis=1)

    def _add_leaf_name_to_dict(self, leaf_name: str) -> str:
        """Add leaf_name or first available name to node_names dict."""
        logging.debug("Adding leaf {} to dictionary.".format(leaf_name if leaf_name else "?"))
        leaf_name = leaf_name if leaf_name else self._first_unused_leaf_name()
        new_number = self.number_of_nodes
        self.node_names.put(leaf_name, new_number)
        self.leaf_names.add(leaf_name)
        return leaf_name

    def add_leaf_to_edge(self, parent: str, child: str, leaf_name: str = None) -> (str, str):
        """Add node between parent and child and attach leaf with leaf_name to it."""
        logging.debug("Adding leaf {} between node {} and {}.".format(leaf_name if leaf_name else "?", parent, child))
        internal_name = self._add_internal_node_to_network()
        leaf_name = self._add_leaf_to_network(leaf_name)

        self._remove_connection(parent, child)
        self._add_connection(parent, internal_name)
        self._add_connection(internal_name, child)
        self._add_connection(internal_name, leaf_name)
        return internal_name, leaf_name

    def _first_unused_leaf_name(self) -> str:
        """Find first unused leaf name in alphabetical order."""
        logging.debug("Retrieving first unused leaf name.")
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
                    if not self._is_leaf_node(leaf_name):
                        found = True
                        break
                except KeyError:
                    found = True
                    break
        assert leaf_name is not None, "Could not find a leaf name."
        logging.debug("First unused leaf name is {}.".format(leaf_name))
        return leaf_name

    def get_root_name(self) -> str:
        """Retrieve the name of the root."""
        logging.debug("Retrieving name of root node.")
        col_sum = sum(self.adj_matrix > 0)
        roots = np.where(col_sum == 0)[0]
        return self._get_node_name(roots[0])

    def get_reticulations(self) -> list:
        """Retrieve node names of all reticulations."""
        logging.debug("Retrieving names of all reticulations.")
        in_degrees = np.array(self.get_in_degrees())
        reticulations = list(np.where(in_degrees > 1)[0])
        return self._get_node_names(reticulations)

    def get_leaves(self) -> set:
        """Retrieve node names of all leaves."""
        logging.debug("Retrieving names of all leaves.")
        return self.leaf_names

    def get_in_degree_node(self, node_name: str, leafless: bool = False) -> int:
        """Retrieve number of arcs entering node_name."""
        logging.debug("Getting {}in degree of node {}.".format("leafless " if leafless else "", node_name))
        node_number = self._get_node_number(node_name)
        leafless_indicator = self.number_of_internals if leafless else self.number_of_nodes
        mask = self.adj_matrix[:leafless_indicator, node_number] > 0
        in_degree = sum(self.adj_matrix[:leafless_indicator, node_number][mask])
        return in_degree

    def get_out_degree_node(self, node_name: str, leafless: bool = False) -> int:
        """Retrieve number of arcs exiting node_name."""
        logging.debug("Getting {}out degree of node {}.".format("leafless " if leafless else "", node_name))
        node_number = self._get_node_number(node_name)
        leafless_indicator = self.number_of_internals if leafless else self.number_of_nodes
        mask = self.adj_matrix[:leafless_indicator, node_number] < 0
        out_degree = abs(sum(self.adj_matrix[:leafless_indicator, node_number][mask]))
        return out_degree

    def get_in_degrees(self, leafless: bool = False) -> list:
        """Retrieve number of arcs entering each node."""
        logging.debug("Getting {}in degree of all node.".format("leafless " if leafless else ""))
        leafless_indicator = self.number_of_internals if leafless else self.number_of_nodes
        mask = self.adj_matrix[:leafless_indicator, :leafless_indicator] > 0
        in_degree = sum(self.adj_matrix[:leafless_indicator, :leafless_indicator] * mask)
        return list(in_degree)

    def get_out_degrees(self, leafless: bool = False) -> list:
        """Retrieve number of arcs entering each node."""
        logging.debug("Getting {}out degree of all nodes.".format("leafless " if leafless else ""))
        leafless_indicator = self.number_of_internals if leafless else self.number_of_nodes
        mask = self.adj_matrix[:leafless_indicator, :leafless_indicator] < 0
        out_degree = -sum(self.adj_matrix[:leafless_indicator, :leafless_indicator] * mask)
        return list(out_degree)

    def get_in_nodes_node(self, node_name: str) -> list:
        """Retrieve all nodes which enter node node_name."""
        logging.debug("Computing nodes entering node {}.".format(node_name))
        node_number = self._get_node_number(node_name)
        in_nodes = []
        i = -1
        while np.any(self.adj_matrix[node_number, :] <= i):
            in_nodes += list(np.where(self.adj_matrix[node_number, :] <= i)[0])
            i -= 1
        return self._get_node_names(in_nodes)

    def get_out_nodes_node(self, node_name: str) -> list:
        """Retrieve all nodes which exit node node_name."""
        logging.debug("Computing nodes exiting node {}.".format(node_name))
        out_nodes = []
        node_number = self._get_node_number(node_name)
        i = 1
        while np.any(self.adj_matrix[node_number, :] >= i):
            out_nodes += list(np.where(self.adj_matrix[node_number, :] >= i)[0])
            i += 1
        return self._get_node_names(out_nodes)

    def get_in_nodes_component(self, component: list) -> list:
        """Retrieve all nodes which enter component."""
        logging.debug("Computing nodes entering component {}.".format(component))
        in_nodes = np.array([])
        for node_name in component:
            in_nodes = np.array(list(itertools.chain(in_nodes, self.get_in_nodes_node(node_name))))
        return [node for node in in_nodes if node not in component]

    def get_out_nodes_component(self, component: list) -> list:
        """Retrieve all nodes which exit component."""
        logging.debug("Computing nodes exiting component {}.".format(component))
        out_nodes = []
        for node_name in component:
            out_nodes = list(itertools.chain(out_nodes, self.get_out_nodes_node(node_name)))
        return [node for node in out_nodes if node not in component]

    def get_in_out_degree_node(self, node_name: str, leafless: bool = False):
        logging.debug("Getting {}in and out degree of node {}.".format("leafless " if leafless else "", node_name))
        return self.get_in_degree_node(node_name, leafless), self.get_out_degree_node(node_name, leafless)

    def get_connections_node(self, node_name: str) -> list:
        """Retrieve all nodes which enter and exit node."""
        logging.debug("Computing nodes entering and exiting {}.".format(node_name))
        return self.get_in_nodes_node(node_name) + self.get_out_nodes_node(node_name)

    def get_leaf_children(self, parents: set, max_depth: int = -1) -> set:
        """Retrieve all leaves max_depth below parents."""
        logging.debug("Getting leaf children of {}{}.".format(parents, " up to depth {}".format(max_depth) if max_depth >= 0 else ""))
        children = self.get_children(parents, max_depth)
        leaves = set([child for child in children if self._is_leaf_node(child)])
        return leaves

    def get_children(self, parents: set, max_depth: int = -1) -> set:
        """Retrieve all nodes max_depth below parents."""
        logging.debug("Getting children of {}{}.".format(parents, " up to depth {}".format(max_depth) if max_depth >= 0 else ""))
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
        logging.debug("Getting parents of {}{}.".format(children, " up to depth {}".format(max_height) if max_height >= 0 else ""))
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

    def _get_connected_components(self) -> list:
        """Retrieve the connected components (excluding leaves) of the network."""
        logging.debug("Computing connected components.")
        # Finds set of connected components
        # Start with a node
        # Iteratively checks which nodes it is connected to through 'unchecked_connections'
        unchecked_nodes = set(self.node_names)
        components = []
        while unchecked_nodes:
            node_name = unchecked_nodes.pop()
            component = [node_name]
            unchecked_connections = set(self.get_connections_node(node_name)).intersection(unchecked_nodes)
            while unchecked_connections:
                connection_name = unchecked_connections.pop()
                unchecked_nodes.remove(connection_name)
                new_connections = set(self.get_connections_node(connection_name))
                unchecked_connections = unchecked_connections.union(new_connections).intersection(unchecked_nodes)
                component.append(connection_name)
                components.append(component)
        return components

    def get_biconnected_components(self) -> list:
        """Compute biconnected components."""
        logging.debug("Computing biconnected components")
        cut_arc_matrix = self._get_cut_arc_matrix()
        cut_adj_matrix = copy.deepcopy(self.adj_matrix) - cut_arc_matrix
        cut_network = RootedLevelKNetwork(adj_matrix=cut_adj_matrix, node_names=self.node_names, leaf_names=self.leaf_names, level=self.level,
                                          dimension=self.dimension)
        return cut_network._get_connected_components()

    def is_biconnected(self):
        """Check if network is biconnected"""
        logging.debug("Checking if network is biconnected.")
        biconnected_components = self.get_biconnected_components()
        return len(biconnected_components) == 1

    def cut_arc_sets(self) -> list:
        """Retrieve cut-arc sets of network."""
        logging.debug("Computing cut-arc sets of network.")
        cut_arc_matrix = self._get_cut_arc_matrix()
        _, to_nodes = np.where(cut_arc_matrix == 1)
        result = []
        for to_node in to_nodes:
            ca_set = np.array(self.get_leaf_children({self._get_node_name(to_node)}))
            ca_set.sort()
            result.append(ca_set)
        return result

    def _get_cut_arc_matrix(self) -> np.ndarray:
        """Compute indicator matrix for arcs which are cut-arcs."""
        logging.debug("Computing cut-arc matrix.")
        # Finds set of cut_arcs components
        # Creates similar network from current network without one edge
        # Checks their connected components

        # Find edges in network
        cut_arc_matrix = np.zeros_like(self.adj_matrix)
        edges = self.get_edges()

        # Loop through all edges
        for from_node_name, to_node_name in edges:
            if self._is_cut_arc(from_node_name, to_node_name):
                from_node_number, to_node_number = self._get_node_number(from_node_name), self._get_node_number(to_node_name)
                cut_arc_matrix[from_node_number][to_node_number] = 1
                cut_arc_matrix[to_node_number][from_node_number] = -1
        return cut_arc_matrix.astype(int)

    def _is_cut_arc(self, from_node_name: str, to_node_name: str) -> bool:
        """Check if arc (from_node, to_node) is a cut-arc."""
        logging.debug("Checking if arc from {} to {} is a cut-arc.".format(from_node_name, to_node_name))
        self._remove_connection(from_node_name, to_node_name)
        new_components = self._get_connected_components()
        self._add_connection(from_node_name, to_node_name)
        if len(new_components) > 1:
            return True
        return False

    def _is_leaf_node(self, node_name: str) -> bool:
        """Check if node_name is leaf."""
        logging.debug("Checking if {} is a leaf.".format(node_name))
        return node_name in self.leaf_names

    def _is_reticulation_node(self, node_name: str) -> bool:
        """Check if node_name is reticulation."""
        logging.debug("Checking if {} is a reticulation.".format(node_name))
        return self.get_in_degree_node(node_name) > 1

    def _is_root_node(self, node_name: str) -> bool:
        """Check if node_name is root."""
        logging.debug("Checking if {} is a root.".format(node_name))
        return self.get_in_degree_node(node_name) == 0
        return self.get_in_degree_node(node_name) == 0

    def _is_structure_node(self, node_name: str) -> bool:
        """Check if node_name is necessary for the structure: not directly connected to any leaves."""
        logging.debug("Checking if {} is a structure node.".format(node_name))
        return (not self._is_leaf_node(node_name)) and self.get_out_degree_node(node_name, leafless=True) == 0

    def _is_connection_node(self, node_name: str) -> bool:
        """Check if node node_name is not part of underlying generator or leaf"""
        logging.debug("Checking if {} is a connection node.".format(node_name))
        return self.get_node_type(node_name) == "connection"

    def get_node_type(self, node_name: str) -> str:
        """Find out what type of node node_name is."""
        logging.debug("Compute what type of node node {} is.".format(node_name))
        if self._is_leaf_node(node_name):
            return "leaf"
        if self._is_reticulation_node(node_name):
            return "reticulation"
        if self._is_root_node(node_name):
            return "root"
        if self._is_structure_node(node_name):
            return "structure"
        return "connection"

    def _is_generator_node(self, node_name: str) -> bool:
        """Check if node_name is part of the underlying generator."""
        logging.debug("Checking if {} is a generator node.".format(node_name))
        return self.get_node_type(node_name) in ("structure", "reticulation", "root")

    def get_generator_nodes(self) -> list:
        """Retrieve all generator nodes"""
        logging.debug("Retrieving all generator nodes.")
        result = []
        for node in self.node_names.keys():
            if self._is_generator_node(node):
                result.append(node)
        return result

    def is_generator(self) -> bool:
        """Check if the network is isomorphic to a generator."""
        logging.debug("Checking if network is isomorphic to a generator.")
        for node_name in self.node_names.keys():
            if not self._is_connection_node(node_name):
                return False
        return True

    def get_generator_sub_network(self):
        """Retrieve the underlying generator of the network."""
        logging.debug("Retrieving the underlying generator of the network.")
        generator_nodes = self.get_generator_nodes()
        generator_nodes_leaves = self.get_leaf_children(set(generator_nodes), 1)
        return RootedLevelKNetwork.from_network(self, generator_nodes_leaves)

    def get_edges(self, leafless: bool = False) -> list:
        """Retrieve all the edges (from_node, to_node) in the network."""
        logging.debug("Retrieving all the {}edges in the network.".format("leafless " if leafless else ""))
        check_size = self.number_of_internals if leafless else self.number_of_nodes
        i = 1
        from_nodes, to_nodes = [], []
        while np.any(self.adj_matrix[:check_size, :check_size] >= i):
            extra_from_nodes, extra_to_nodes = np.where(self.adj_matrix[:check_size, :check_size] >= i)
            from_nodes += list(self._get_node_names(extra_from_nodes))
            to_nodes += list(self._get_node_names(extra_to_nodes))
            i += 1
        return [[from_node, to_node] for from_node, to_node in zip(np.array(from_nodes), np.array(to_nodes))]

    def prune(self, suppress_redundant: str = 'all', suppress_parallel: bool = True):
        """Suppress al unnecessary/redundant/parallel nodes/edges in network."""
        logging.debug("Pruning: Suppressing {} redundant components{}.".format(suppress_redundant, " and parallel edges" if suppress_parallel else ""))
        assert suppress_redundant in ('none', 'strongly', 'all'), "suppress_redundant parameter must be one of none, strongly, all."

        changes = True

        while changes:
            changes = False
            leaf_numbers = set(self._get_leaf_numbers())
            in_degrees, out_degrees = np.array(self.get_in_degrees()), np.array(self.get_out_degrees())

            # Get non-leaf nodes without children
            childless_nodes = set(np.where(out_degrees == 0)[0])
            childless_nodes.difference_update(leaf_numbers)
            childless_nodes = self._get_node_names(list(childless_nodes))

            # Get nodes with no entering arcs and only one exiting arc
            parentless_nodes = set(np.where(in_degrees == 0)[0])
            one_exiting_arc_nodes = set(np.where(out_degrees == 1)[0])
            parentless_one_exiting_arc_nodes = parentless_nodes.intersection(one_exiting_arc_nodes)
            parentless_one_exiting_arc_nodes = self._get_node_names(list(parentless_one_exiting_arc_nodes))

            # Get nodes with one in and one out node
            one_entering_arc_nodes = set(np.where(in_degrees == 1)[0])
            one_entering_and_exiting_arc_nodes = one_entering_arc_nodes.intersection(one_exiting_arc_nodes)
            one_entering_and_exiting_arc_nodes = self._get_node_names(list(one_entering_and_exiting_arc_nodes))
            #
            # print("childless", childless_nodes)
            # print("parentless one exiting", parentless_one_exiting_arc_nodes)
            # print("one enter one exit", one_entering_and_exiting_arc_nodes)

            nodes_to_remove = list(set(childless_nodes + parentless_one_exiting_arc_nodes + one_entering_and_exiting_arc_nodes))
            for node in nodes_to_remove:
                changes = True
                self._suppress_node(node)

            # Check for parallel arcs
            if suppress_parallel:
                parallel_arcs = np.where(self.adj_matrix >= 2)
                parallel_arcs = zip(self._get_node_names(parallel_arcs[0]), self._get_node_names(parallel_arcs[1]))
                # print("Parallel", parallel_arcs)
                for from_node, to_node in parallel_arcs:
                    changes = True
                    self._remove_connection(from_node, to_node)

            # Check for redundant components:
            if suppress_redundant != 'none':
                bcs_to_remove = []
                bcs = self.get_biconnected_components()
                print("---------", bcs)
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
        logging.debug("Checking if component {} is redundant.".format(component))
        out_nodes = self.get_out_nodes_component(component)
        if len(out_nodes) != 1:
            return False
        return True

    def _is_strongly_redundant_component(self, component: list) -> bool:
        """Check if component is strongly redundant."""
        logging.debug("Checking if component {} is strongly redundant.".format(component))
        out_nodes = self.get_out_nodes_component(component)
        if len(out_nodes) != 1:
            return False
        temp_network = copy.deepcopy(self)
        temp_network._remove_component(component)
        connected_comp = temp_network._get_connected_components()
        if len(connected_comp) > 1:
            return False
        return True

    def contains_leaf(self, component: list) -> bool:
        """Check if component contains a leaf."""
        logging.debug("Checking if component {} contains a leaf.".format(component))
        for node in component:
            if self._is_leaf_node(node):
                return True
        return False

    def number_of_internals_leaves_reticulations(self):
        """Retrieve number of internal, leaves, reticulation nodes."""
        logging.debug("Retrieving number of internal, leaves, reticulation nodes.")
        col_sum = sum(self.adj_matrix > 0)
        number_of_reticulations = sum(col_sum > 1)
        return self.number_of_internals, self.number_of_leaves, number_of_reticulations

    def switch_nodes_in_adj_matrix(self, node_name_1: str, node_name_2: str):
        """Switch rows and columns corresponding to nodes from place in adjacency matrix"""
        logging.debug("Switch rows and columns corresponding to nodes {} and {} from place in adjacency matrix.".format(node_name_1, node_name_2))
        # Switches node place in matrix
        assert not (self._is_leaf_node(node_name_1) ^ self._is_leaf_node(node_name_2)), "WARNING: can not switch nodes {} and {} as exactly one is a leaf".format(node_name_1, node_name_2)

        x = self._get_node_number(node_name_1)
        y = self._get_node_number(node_name_2)

        # Adjacency matrix: Switching rows
        row_x = copy.deepcopy(self.adj_matrix[x, :])
        row_y = copy.deepcopy(self.adj_matrix[y, :])
        self.adj_matrix[x, :] = row_y
        self.adj_matrix[y, :] = row_x

        # Adjacency matrix: Switching cols
        col_x = copy.deepcopy(self.adj_matrix[:, x])
        col_y = copy.deepcopy(self.adj_matrix[:, y])
        self.adj_matrix[:, x] = col_y
        self.adj_matrix[:, y] = col_x

        # Switch node numbers in dictionary
        x = self.node_names.pop(node_name_1)
        y = self.node_names.pop(node_name_2)
        self.node_names.put(node_name_1, y)
        self.node_names.put(node_name_2, x)

    def set_node_as_number(self, node_name: str, new_number: int) -> list:
        """Set row and column of node as the new_number-th in adjacency matrix."""
        logging.debug("Setting row and column of node {} as the {}-th in adjacency matrix.".format(node_name, new_number))
        current_node_number = self._get_node_number(node_name)
        if current_node_number != new_number:
            # Change position to new_number --> new_number currently occupied by Y, have to switch X and Y
            node_name_of_new_number = self._get_node_name(new_number)
            self.switch_nodes_in_adj_matrix(node_name, node_name_of_new_number)
            return [node_name, node_name_of_new_number]
        return [node_name]

    def get_exhibited_trinets(self) -> list:
        """Retrieve all trinets exhibited by network."""
        logging.debug("Retrieving all trinets exhibited by network {}.".format(self.uid))
        leaves = self._get_node_names(list(range(self.number_of_internals, self.number_of_nodes)))
        triplets = itertools.combinations(leaves, 3)
        trinets = []
        for triplet in tqdm(triplets):
            current_trinet = RootedLevelKNetwork.trinet_from_network(self, set(triplet))
            current_trinet.to_standard_form()
            trinets.append((triplet, current_trinet))
        return trinets

    def visualize(self):
        """Visualize network."""
        logging.debug("Visualizing network.")
        dot = Digraph()
        dot.engine = 'dot'
        edges = self.get_edges()
        for edge in edges:
            dot.edge(edge[0], edge[1])

        dot.render(view=True)

    def to_df(self, directed: bool = True):
        """Retrieve dataframe representation of network."""
        logging.debug("Retrieve {}directed dataframe representation of network.".format("" if directed else "un-"))
        ordered_node_names = self._get_ordered_node_names()
        mask = self.adj_matrix > 0
        data = self.adj_matrix * mask if directed else self.adj_matrix
        length = self.number_of_internals if directed else self.number_of_nodes
        result = pd.DataFrame(columns=ordered_node_names, index=ordered_node_names[:length], data=data[:length])
        return result

    def _get_component_list_names(self, component_list: list) -> list:
        """Retrieve names of list of lists of node numbers."""
        logging.debug("Retrieve node names of {}.".format(component_list))
        return [self._get_node_name(comp) for comp in component_list]

    def __eq__(self, other):
        if not self.__class__ == other.__class__:
            return False
        if not self.number_of_internals == other.number_of_internals:
            return False
        if not np.array_equal(self.adj_matrix[:, :self.number_of_internals],
                              other.adj_matrix[:, :other.number_of_internals]):
            return False
        return True

    def __str__(self):
        return str(self.node_names) + str("\n") + str(self.adj_matrix)
