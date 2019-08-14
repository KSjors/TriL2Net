import numpy as np
import copy
from utils.help_functions import coalesce
from tarjan import tarjan
import itertools
from datastructures.rooted_level_k_network import RootedLevelKNetwork
import pprint

pp = pprint.PrettyPrinter(indent=4)


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
        trinet = RootedLevelKNetwork.from_network(trinet_info, leaf_names, suppress_parallel=True, suppress_redundant='all')
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
            except IndexError:
                # if len(trinet_info.trinet.leaf_names) == 2:
                    # trinet_info.trinet.visualize()
                pass

    def shrink(self, leaf_set):
        for trinet_info in iter(self):
            if not trinet_info.shrink(leaf_set):
                self.remove(trinet_info)

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
        assert len(leaf_names) >= 2, "Leaf_names should always include at least two leaves"
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
                    adjacency_dict[arc[0]] = set(arc[1])

        strongly_connected_components = tarjan(adjacency_dict)
        msss = []
        for strongly_connected_component in strongly_connected_components:
            for node in strongly_connected_component:
                to_leaves = adjacency_dict[node]
                if not set(to_leaves).issubset(strongly_connected_component):
                    break
            else:
                msss.append(strongly_connected_component)

        return [mss for mss in msss if len(mss) > 1 and len(mss) != len(leaf_set)]

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
        return generator

    def leaf_locations(self):
        leaf_set = copy.deepcopy(self.represented_leaves())
        leaf_on_edge_count_dict = {leaf: {} for leaf in leaf_set}
        for trinet_info in self:
            relation_dict = trinet_info['relation_dict']
            for leaf, edge in trinet_info['extra_leaf_dict'].items():
                translated_leaf = relation_dict[leaf]
                try:
                    leaf_on_edge_count_dict[translated_leaf][edge] += 1
                except KeyError:
                    leaf_on_edge_count_dict[translated_leaf][edge] = 1

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
                edge_leaf_dict[best_edge] = set(leaf)
            leaf_set.remove(leaf)

        return edge_leaf_dict, {leaf: relation_dict.inverse[leaf] for leaf in leaf_set}

    # def order_leaves(self, ):
