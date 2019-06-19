import logging
from datastructures.omega import *
from datastructures.rooted_level_k_network import *
from data.all_trinets import *
from datastructures.trinet_set import TrinetSet, MssTrinetSet
from utils.help_functions import guid, leaf_names_to_identifier, mss_leaf_name
import copy
import time

# regenerate_trinets()
all_generators, trinet_lookup_dict = get_trinets()


class Solver:
    def __init__(self, trinet_set: TrinetSet):
        self.uid = guid()
        self.logger_name = f"solver.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)
        self.logger.debug(f"Creating new solver object from trinet set {trinet_set.uid}.")
        self.trinet_sets = [trinet_set]
        self.transformations = {}
        self.components = {}

    def solve(self):
        while self.next_transform():
            continue
        next(iter(self.components.keys())).visualize()

    def next_transform(self):
        self.logger.info("Performing the next transformation.")
        if len(self.trinet_sets[-1].taxa_names) == 2:
            self.last_transform()
            return False
        elif  len(self.trinet_sets[-1].taxa_names) == 1:
            return False
        mss = self.get_next_mss()
        biconnected_component = self.transform(mss)
        # biconnected_component.visualize()
        component, taxa = self.expand_mss(biconnected_component)
        self.components[component] = taxa
        self.logger.info(f"Created component {biconnected_component.uid} for minimal sink-set {mss}.")
        return True

    def last_transform(self):
        self.logger.info("Performing last transformation.")
        # If there are only two taxa [ie.:  '(AB)C' and '(DE)'] left, no trinet can be found. Therefore take one leaf from each and a random third [ie.: A,D and C]
        taxa_sets = list(self.components.values())
        components = list(self.components.keys())

        # One leaf from each
        important_taxa = {}
        for taxa_set, component in zip(taxa_sets, components):
            for taxus in taxa_set:
                if '(' not in taxus:
                    important_taxa[component] = taxus
                    break

        # Random leaf
        original_taxa = list(self.trinet_sets[0].taxa_names.keys())
        random_leaf = list(set(original_taxa).difference(set(important_taxa.values())).pop())

        # Triplet
        triplet = list(important_taxa.values()) + random_leaf
        trinet_identifier = leaf_names_to_identifier(triplet)
        trinet = self.trinet_sets[0].trinet_dict[trinet_identifier]
        final_network = RootedLevelKNetwork.from_network(trinet, list(important_taxa.values()), suppress_redundant='all', suppress_parallel=True)
        for taxus in important_taxa.values():
            final_network.rename_node(taxus, taxus + '_temp')

        # final_network.visualize()

        for component, taxus in important_taxa.items():
            final_network.replace_leaf_with_network(taxus + "_temp", component)
            final_network.standardize_internal_node_names()

        final_network.visualize()

    def expand_mss(self, biconnected_component):
        self.logger.info(f"Looking for leaves in component {biconnected_component.uid} to expand.")
        replace_dict = {}
        all_taxa = biconnected_component.leaf_names
        self.logger.info(f"Component has taxa {all_taxa}.")
        for leaf_name in biconnected_component.leaf_names:
            if '(' in leaf_name:
                self.logger.info(f"Need to replace {leaf_name}")
                for component, taxa in self.components.items():
                    if self.transformations[leaf_name] == taxa:
                        self.logger.info(f"Found {component.uid} to replace {leaf_name}.")
                        replace_dict[leaf_name] = component

        if replace_dict:
            for leaf_name, component in replace_dict.items():
                self.logger.info(f"Expanding {leaf_name} using {component.uid}")
                biconnected_component.replace_leaf_with_network(leaf_name, component)
                self.components.pop(component)
        else:
            self.logger.info("No leaves found to expand.")

        biconnected_component.standardize_internal_node_names()

        return biconnected_component, sorted(all_taxa)

    def get_next_mss(self):
        # TODO all found mss in one go
        current_omega = Omega(self.trinet_sets[-1])
        current_mss = current_omega.minimal_sink_sets()
        current_mss.sort(key=len)
        for mss in current_mss:
            if len(mss) >= 2:
                self.logger.info(f"Next minimal sink-set is {mss}")
                return mss

    def transform(self, mss: list):
        self.logger.info(f"Transforming {mss}.")
        if len(mss) == 2:
            # TODO: discard all faulty trinets first
            # Find trinet of which mss is part
            current_taxa = list(self.trinet_sets[-1].taxa_names)
            third_leaf = current_taxa.pop(-1)
            while third_leaf in mss:
                third_leaf = current_taxa.pop(-1)
            trinet_identifier = leaf_names_to_identifier(mss + [third_leaf])
            trinet_of_mss = self.trinet_sets[-1].trinet_dict[trinet_identifier]
            biconnected_component = RootedLevelKNetwork.from_network(trinet_of_mss, mss)
        elif len(mss) == 3:
            # TODO what if this is the only trinet that disagrees
            trinet_identifier = leaf_names_to_identifier(mss)
            biconnected_component = self.trinet_sets[-1].trinet_dict[trinet_identifier]
            # mss_network.visualize()
        else:
            self.logger.info(f"MSS is represented by trinets {self.trinet_sets[-1].trinet_dict.keys()}.")
            mss_trinets = MssTrinetSet.from_trinet_set(trinet_set=self.trinet_sets[-1], mss=mss)
            # mss_trinets.set_underlying_generator()
            biconnected_component = mss_trinets.get_underlying_network()
            # biconnected_component = self.get_network(mss_trinets.trinet_dict, mss)

        # TODO: suppress in @classmethod?
        # TODO: save network of mss
        new_trinet_set = TrinetSet.copy_trinet_set(self.trinet_sets[-1])
        new_name = new_trinet_set.suppress_minimal_sink_set(mss)
        self.transformations[new_name] = sorted(mss)
        self.trinet_sets.append(new_trinet_set)
        return biconnected_component

    def state(self, number=-1):
        assert -1 <= number < len(self.trinet_sets), "number not in range"
        current_omega = Omega(self.trinet_sets[number])
        current_omega.visualize()

    def __str__(self):
        return str(self.trinet_sets[-1]) + " " + str(self.transformations[-1])

    def __getstate__(self):
        self.logger_name = 'solver.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger_name)
        return self.__dict__

    @staticmethod
    def get_network(trinets, mss):
        network_info = {}
        for level, generator_list in all_generators.items():
            network_info[level] = {}
            for gen in generator_list:
                network_info[level][gen] = {}

        for trinet_taxa, trinet in trinets.items():
            # Find out what the trinet structure of the trinet is
            trinet_info = trinet_lookup_dict[trinet]
            network_info[trinet_info['generator'].level][trinet_info['generator']][trinet] = trinet_info['on_edges']

        # Group the trinet structures
        generator_level_count = {}
        for level, generator_list in all_generators.items():
            generator_level_count[level] = 0
            for gen in generator_list:
                generator_level_count[level] += len(network_info[level][gen])

        # Take maximum found generator level
        found_levels = [level for level, count in generator_level_count.items() if level != 0]
        underlying_level = max(found_levels)

        # Count occurrences of generators
        trinets_per_generator = [len(trinets) for _, trinets in network_info[underlying_level].items()]
        # Take the generators with the max occurrences
        max_count_generators = [gen for gen, trinets in network_info[underlying_level].items() if len(trinets) == max(trinets_per_generator)]
        # Take the first of these (what else?)
        underlying_generator = max_count_generators[0]

        # now we know the underlying generator, we need to find out on which edges on the generator the leaf have been added
        # save on which edge the leaf has been added in each trinet
        number_of_added_leaves = len(next(iter(network_info[underlying_level][underlying_generator].values())))

        leaf_on_edge_dict = {}
        for leaf in mss:
            leaf_on_edge_dict[leaf] = {}

        for trinet, on_edges in network_info[underlying_level][underlying_generator].items():
            for leaf, edge in zip(trinet.leaf_names[3 - number_of_added_leaves:], on_edges):
                if edge in leaf_on_edge_dict[leaf].keys():
                    leaf_on_edge_dict[leaf][edge] += 1
                else:
                    leaf_on_edge_dict[leaf][edge] = 1

        # For each leaf take the edge which occurs the most often. If there are no edges it occurs onm then it is a reticulation
        reticulations = []
        actual_leaf_on_edge_dict = {}
        for leaf, on_edges in leaf_on_edge_dict.items():
            if len(on_edges.values()) == 0:
                reticulations.append(leaf)
                continue
            max_occurrence = max(on_edges.values())
            max_occurring_edges = [edge for edge, occurrence in list(on_edges.items()) if occurrence == max_occurrence]
            actual_leaf_on_edge_dict[leaf] = max_occurring_edges[0]

        # From leaf: edge dict to edge: leaf dict
        on_edge_leaf_dict = {}
        for leaf, on_edge in actual_leaf_on_edge_dict.items():
            on_edge_leaf_dict[on_edge] = on_edge_leaf_dict.get(on_edge, [])
            on_edge_leaf_dict[on_edge].append(leaf)

        # Look at trinet with all reticulations present to find the ordering of the reticualtions TODO look at all trinets with all reticulations and count etc.
        reticulation_order = next(iter(network_info[underlying_level][underlying_generator].keys())).leaf_names[0:len(reticulations)]

        biconnected_component = RootedLevelKNetwork.copy_network(underlying_generator)

        # Rename reticulations
        test_network_leaf_names = biconnected_component.leaf_names
        for old_leaf_name, new_leaf_name in zip(test_network_leaf_names, reticulation_order):
            biconnected_component.rename_node(old_leaf_name, new_leaf_name)

        # Add other leaves to edges
        for on_edge, leaves in on_edge_leaf_dict.items():
            parent = on_edge[0]
            for leaf in leaves:
                extra_internal_node, _ = biconnected_component.add_leaf_to_edge(parent, on_edge[1], leaf)
                parent = extra_internal_node
        return biconnected_component
