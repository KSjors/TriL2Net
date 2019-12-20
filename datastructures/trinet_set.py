# # -*- coding: utf-8 -*-
# """
# Created on Tue Feb 26 14:11:45 2019
#
# @author: Sjors
# """
#
# import logging
# from utils.help_functions import guid, leaf_names_to_identifier, mss_leaf_name
# from datastructures.rooted_level_k_network import  RootedLevelKNetwork, TrinetInfoList, TrinetInfo, RootedLevelKGenerator
# from data.all_trinets import get_trinets
# from graphviz import Digraph
# import pprint
# import itertools
#
# all_generators, standard_trinet_info_list = get_trinets()
#
# pp = pprint.PrettyPrinter(indent=4)
#
#
# class TrinetSet:
#     def __init__(self, trinet_info_list: TrinetInfoList, taxa_names: bidict):
#         self.uid = guid()
#         self.logger = logging.getLogger("Created trinet set {}".format(self.uid))
#         self.trinet_info_list = trinet_info_list
#         self.taxa_names = taxa_names
#         self.number_of_taxa = len(taxa_names)
#
#     @classmethod
#     def copy_trinet_set(cls, trinet_set):
#         logging.debug("Copying trinet_set {}.".format(trinet_set.uid))
#         trinet_set.logger.debug("Being copied.")
#         trinet_dict = copy.deepcopy(trinet_set.trinet_dict)
#         taxa_names = copy.deepcopy(trinet_set.taxa_names)
#         trinet_set_copy = cls(trinet_dict, taxa_names)
#         trinet_set_copy.logger.debug("Copy from {}.".format(trinet_set.uid))
#         return trinet_set_copy
#
#     @classmethod
#     def from_trinet_info_list(cls, trinet_info_list: TrinetInfoList):
#         logging.debug("Creating trinet set from triplet trinet list.")
#         trinet_info_list = TrinetInfoList()
#         taxa_names = bidict()
#         trinet_set = cls(trinet_info_list, taxa_names)
#         for trinet_info in trinet_info_list:
#             trinet_set.add_trinet(trinet_info)
#
#         trinet_set.logger.debug("Created from triplet trinet list.")
#         return trinet_set
#
#     def add_trinet(self, trinet_info: TrinetInfo):
#         """Add trinet to dictionary."""
#         trinet = trinet_info.trinet
#         triplet = trinet.leaf_names
#         self.logger.debug("Adding trinet {} with leaves {}.".format(trinet.uid, triplet))
#         assert len(triplet) == 3, "Can not add trinet {} to trinet_set {} as it does not have exactly three leaves.".format(trinet.uid, self.uid)
#         assert self.trinet_info_list.contains_trinet_with_leaf_names(tuple), "Can not add trinet {} to trinet_set {} as it already exists.".format(trinet.uid, self.uid)
#         self.trinet_info_list.append(trinet_info)
#         for X in triplet:
#             if X not in self.taxa_names.keys():
#                 self.taxa_names[X] = self.number_of_taxa
#                 self.number_of_taxa += 1
#         return trinet_info
#
#     def remove_trinet(self, triplet: list):
#         """Remove trinet with taxa = triplet from dictionary."""
#         self.logger.debug("Removing triplet {}.".format(triplet))
#         assert len(triplet) == 3, "Can not remove {} as it is not a triplet.".format(triplet)
#         trinet_info = self.trinet_info_list.remove_trinet_by_leaf_names(triplet)
#         trinet = trinet_info.trinet
#         self.logger.debug("Triplet {} corresponds to trinet {}.".format(triplet, trinet.uid))
#         return trinet_info
#
#     def cut_arc_sets_per_triplet(self) -> dict:
#         """Calculate cut-arc sets of each trinet."""
#         self.logger.debug("Calculating cut-arc sets per trinet.")
#         cut_arc_sets = {}
#         for trinet in self.trinet_info_list:
#             cut_arc_sets[trinet.leaf_names] = trinet.cut_arc_sets
#         return cut_arc_sets
#
#     def suppress_minimal_sink_set(self, mss: list) -> str:
#         """Suppress minimal sink set."""
#         self.logger.debug("Suppressing minimal sink set {}.".format(mss))
#         current_taxa_names = copy.deepcopy(list(self.taxa_names))
#
#         # Remove all binets and trinets that are a subset of minimal sink set mss
#         if len(mss) > 2:
#             triplet_iterator = itertools.combinations(mss, 3)
#             for triplet in triplet_iterator:
#                 try:
#                     self.trinet_info_list.remove_trinet_by_leaf_names(triplet)
#                 except KeyError:
#                     pass
#
#         binet_iterator = itertools.combinations(mss, 2)
#         for binet in binet_iterator:
#             third_leaf_iterator = itertools.combinations(set(current_taxa_names) - set(mss), 1)
#             for third_leaf in third_leaf_iterator:
#                 triplet = list(binet) + list(third_leaf)
#                 self.trinet_info_list.remove_trinet_by_leaf_names(triplet)
#
#         # Replace name of each leaf in mss with combined name
#         new_name = mss_leaf_name(mss)
#         for leaf in mss:
#             binet_iterator = itertools.combinations(set(current_taxa_names) - set(mss), 2)
#             for binet in binet_iterator:
#                 leaves = list(binet) + list([leaf])
#                 trinet_network = self.trinet_info_list.remove_trinet_by_leaf_names(leaves)
#                 trinet_network.rename_node(leaf, new_name)
#                 self.trinet_info_list.append(trinet_network)
#
#         for leaf in mss:
#             self.remove_leaf_name(leaf)
#
#         self.add_leaf_name(new_name)
#         return new_name
#
#     def remove_leaf_name(self, node_name: str) -> int:
#         """Remove leaf name from taxa_names and lower all numbers of names above by one"""
#         self.logger.debug("Removing leaf name {}.".format(node_name))
#         node_number = self.taxa_names.pop(node_name)
#         for y in range(node_number + 1, self.number_of_taxa):
#             self.taxa_names[self.taxa_names.inverse[y]] -= 1
#         self.number_of_taxa -= 1
#         return node_number
#
#     def add_leaf_name(self, node_name: str) -> str:
#         """Add leaf name."""
#         self.logger.debug("Adding leaf name {}.".format(node_name))
#         new_number = self.number_of_taxa
#         self.taxa_names.put(node_name, new_number)
#         self.number_of_taxa += 1
#         return node_name
#
#     def __str__(self) -> str:
#         return str(self.taxa_names)
#
#     def __getstate__(self):
#         self.logger = 'trinet_set.{}'.format(self.uid)
#         result = copy.deepcopy(self.__dict__)
#         self.logger = logging.getLogger(self.logger)
#         return result
#
#     def __setstate__(self, d):
#         self.__dict__ = d
#         self.logger = logging.getLogger(self.logger)
#         return self.__dict__
#
#
# class MssTrinetSet(TrinetSet):
#     def __init__(self, trinet_info_list: TrinetInfoList, taxa_names: bidict):
#         super().__init__(trinet_info_list, taxa_names)
#
#         for trinet_info in self.trinet_info_list:
#             # Find out what the trinet structure of the trinet is
#             equal_structured_trinet = standard_trinet_info_list.find_equal_structured_trinet(trinet_info)
#             trinet_info.add_info(equal_structured_trinet)
#
#         self.inconsistent_trinet_dict = dict()
#         self.missing_trinets = []
#
#         self.underlying_generator = None
#         self.underlying_generator_level = None
#
#         self.reticulations = set()
#         self.on_edge_leaf_dict = dict()
#
#     @classmethod
#     def from_trinet_set(cls, trinet_set: TrinetSet, mss: list):
#         result = cls(TrinetInfoList(), bidict())
#
#         result.logger.info(f"MSS is {mss}")
#         triplet_iterator = itertools.combinations(mss, 3)
#         for triplet in triplet_iterator:
#             try:
#                 trinet_info = trinet_set.trinet_info_list.find_trinet_by_leaf_names(triplet)
#                 result.add_trinet(trinet_info)
#             except KeyError:
#                 result.missing_trinets.append(triplet)
#
#         return result
#
#     def add_trinet(self, trinet: RootedLevelKNetwork):
#         trinet_info = super().add_trinet(trinet)
#         equal_structured_trinet = standard_trinet_info_list.find_equal_structured_trinet(trinet)
#         trinet_info.add_info(equal_structured_trinet)
#         self.logger.info(f"Added trinet with leaves {trinet.leaf_names}.")
#
#     def set_underlying_generator_level(self, method=max):
#         self.logger.info("Computing underlying generator level")
#         # Group the trinet structures
#         generator_level_count = {}
#         for trinet_info in self.trinet_info_list:
#             try:
#                 generator_level_count[trinet_info['level']] += 1
#             except KeyError:
#                 generator_level_count[trinet_info['level']] = 1
#
#         # Take maximum found generator level
#         found_levels = [level for level, count in generator_level_count.items() if count > 0]
#         underlying_level = method(found_levels)
#         self.underlying_generator_level = underlying_level
#         self.logger.info(f"Set underlying generator level to {self.underlying_generator_level}.")
#
#     def set_underlying_generator(self, method=lambda a: a[0]):
#         self.logger.info("Computing underlying generator")
#         if self.underlying_generator_level is None:
#             self.set_underlying_generator_level()
#         # Count occurrences of generators
#         trinets_per_generator = [len(trinet_identifier_dict) for _, trinet_identifier_dict in self.mss_info[self.underlying_generator_level].items()]
#         # Take the generators with the max occurrences
#         max_count_generators = [gen for gen, trinet_identifier_dict in self.mss_info[self.underlying_generator_level].items() if
#                                 len(trinet_identifier_dict) == max(trinets_per_generator)]
#         # Take the first of these (what else?)
#         underlying_generator = method(max_count_generators)
#         self.underlying_generator = underlying_generator
#
#         # Remove all trinets with same or higher underlying level but different underlying generator
#         # TODO: remove inconsistent trinets at all?
#         temp_dict = copy.deepcopy(self.mss_info)
#         for level in range(self.underlying_generator_level, 2 + 1):
#             for generator, leaf_name_on_edge_dict in temp_dict[level].items():
#                 # TODO check if this properly removes 'false' trinets
#                 if generator != underlying_generator:
#                     for leaf_names in leaf_name_on_edge_dict:
#                         trinet_identifier = leaf_names_to_identifier(list(leaf_names))
#                         self.inconsistent_trinet_dict[trinet_identifier] = self.trinet_dict.pop(trinet_identifier)
#                     self.mss_info[level].pop(generator)
#         # TODO: Remove trinets with lower level which are also inconsistent. e.g. If underlying generator is
#         #  2a with reticulation A, there can not be a trinet with reticulation A of level-1.
#         #  However, 2d with reticulation A does allow for a trinet with reticulation A of level-1
#
#     def set_on_edge_leaf_dict_and_reticulation_set(self):
#         self.logger.info("Computing edge-leaf dict and reticulation list")
#         if self.underlying_generator is None:
#             self.set_underlying_generator()
#         number_of_added_leaves = len(next(iter(self.mss_info[self.underlying_generator_level][self.underlying_generator].values())))
#
#         # Dictionary stating on which edges a leaf is found
#         leaf_on_edge_dict = {}
#         for leaf in self.taxa_names.keys():
#             leaf_on_edge_dict[leaf] = {}
#
#         # Fill dictionary
#         for leaf_names, on_edges in self.mss_info[self.underlying_generator_level][self.underlying_generator].items():
#             for leaf, edge in zip(leaf_names[3 - number_of_added_leaves:], on_edges):
#                 try:
#                     leaf_on_edge_dict[leaf][edge] += 1
#                 except KeyError:
#                     leaf_on_edge_dict[leaf][edge] = 1
#
#         # For each leaf take the edge which occurs the most often. If there are no edges it occurs on then it is a reticulation
#         actual_leaf_on_edge_dict = {}
#         for leaf, on_edges in leaf_on_edge_dict.items():
#             if len(on_edges.values()) == 0:
#                 self.reticulations.add(leaf)
#                 continue
#             max_occurrence = max(on_edges.values())
#             max_occurring_edges = [edge for edge, occurrence in list(on_edges.items()) if occurrence == max_occurrence]
#             actual_leaf_on_edge_dict[leaf] = max_occurring_edges[0]
#         # TODO: remove trinets which have this leaf on a different side?
#
#         # Invert dict. From leaf: edge dict to edge: leaf dict
#         self.on_edge_leaf_dict = {}
#         for leaf, on_edge in actual_leaf_on_edge_dict.items():
#             self.on_edge_leaf_dict[on_edge] = self.on_edge_leaf_dict.get(on_edge, set())
#             self.on_edge_leaf_dict[on_edge].add(leaf)
#
#         for edge, leaves in self.on_edge_leaf_dict.items():
#             self.on_edge_leaf_dict[edge] = self.order_leaves_on_edge(edge, leaves)
#
#     def order_leaves_on_edge(self, edge, leaves):
#         leaf_name_map = bidict()
#         i = 0
#         for leaf in leaves:
#             leaf_name_map.put(leaf, i)
#             i += 1
#         score_matrix = np.zeros((i, i))
#
#         leaves = set(leaves)
#
#         for level, generator_dict in self.mss_info.items():
#             for generator, trinet_info_dict in generator_dict.items():
#                 for leaf_names, on_edges in trinet_info_dict.items():
#                     if len(set(leaf_names).intersection(leaves)) >= 2:
#                         current_leaves = set(leaf_names).intersection(leaves)
#
#                         trinet_identifier = leaf_names_to_identifier(leaf_names)
#                         trinet = self.trinet_dict[trinet_identifier]
#
#                         # leaf_orderings = trinet.get_leaf_ordering()
#                         # for leaf_ordering in leaf_orderings:
#                         #     try:
#                         #         first_index = leaf_ordering.index(current_leaves[0])
#                         #         second_index = leaf_ordering.index(current_leaves[1])
#                         #         if first_index > second_index:
#                         #
#                         #     except:
#                         #
#                         # else:
#                         #     raise AssertionError("Trinet has more than 1 cut-arc set.")
#
#         for reticulation in self.reticulations:
#             # TODO use all trinets?
#             two_leaf_iterator = itertools.combinations(leaves, 2)
#             for two_leaves in two_leaf_iterator:
#                 triplet = list(two_leaves) + [reticulation]
#                 trinet_identifier = leaf_names_to_identifier(triplet)
#                 trinet = self.trinet_dict[trinet_identifier]
#                 level = len(trinet.get_reticulations())
#
#                 cut_arc_set = [cut_arc_set for cut_arc_set in trinet.cut_arc_sets if len(cut_arc_set) == 2]
#                 leaf_names = trinet.leaf_names
#                 leaf_names.remove(reticulation)
#                 if len(cut_arc_set) == 1:
#                     if reticulation not in cut_arc_set[0]:
#                         # No information for order
#                         continue
#                     else:
#                         undeep_leaf_number = leaf_name_map[leaf_names[0]]
#                         deep_leaf_number = leaf_name_map[leaf_names[1]]
#                 else:
#                     undeep_leaf_number = leaf_name_map[leaf_names[0]]
#                     deep_leaf_number = leaf_name_map[leaf_names[1]]
#                 score_matrix[undeep_leaf_number, deep_leaf_number] += 1
#         self.visualize_tournament(score_matrix, leaf_name_map)
#         ranking = self.determine_ranking(score_matrix=score_matrix)
#         ranked_leaves = [leaf_name_map.inverse[number] for number in ranking]
#         return ranked_leaves
#
#     @staticmethod
#     def visualize_tournament(tournament: np.ndarray, leaf_name_map: bidict):
#         dot = Digraph()
#         dot.engine = 'dot'
#
#         for leaf_name, number in leaf_name_map.items():
#             dot.node(str(number), leaf_name)
#
#         for i in range(tournament.shape[0]):
#             for j in range(tournament.shape[1]):
#                 if tournament[i, j] > 0:
#                     dot.edge(str(i), str(j))
#
#         dot.render(view=True)
#         time.sleep(0.2)
#
#     def determine_ranking(self, score_matrix, ranking=None):
#         score_matrix = copy.copy(score_matrix)
#         if ranking is None:
#             ranking = []
#         if sum(sum(score_matrix)) == 0:
#             all_players = list(range(len(score_matrix)))
#             for player in ranking:
#                 all_players.remove(player)
#             return ranking + all_players
#         player_wins = sum(score_matrix.T)
#         best_players = np.where(player_wins == player_wins.max())[0]
#         if len(best_players) == 1 or len(best_players) == len(score_matrix) - len(ranking):
#             # TODO all same ranking in one go?
#             best_player = best_players[0]
#         else:
#             intermediate_ranking = self.determine_ranking(score_matrix[best_players][:, best_players])
#             best_player = best_players[intermediate_ranking[0]]
#         score_matrix[best_player, :] = 0
#         score_matrix[:, best_player] = 0
#         ranking.append(best_player)
#         return self.determine_ranking(score_matrix, ranking)
#
#     def get_underlying_network(self):
#         self.logger.info("Computing underlying network.")
#         if self.underlying_generator is None:
#             self.set_underlying_generator()
#         network = copy.deepcopy(self.underlying_generator)
#         if len(self.on_edge_leaf_dict) + len(self.reticulations) == 0:
#             self.set_on_edge_leaf_dict_and_reticulation_set()
#
#         # Rename reticulations %TODO look at all trinets with all reticulations and count etc.
#         reticulation_order = next(iter(self.mss_info[self.underlying_generator_level][self.underlying_generator].keys()))[0:len(self.reticulations)]
#         network_leaf_names = network.leaf_names
#         for old_leaf_name, new_leaf_name in zip(network_leaf_names, reticulation_order):
#             network.rename_node(old_name=old_leaf_name, new_name=new_leaf_name)
#
#         self.logger.info(f"Leaves to add are {list(self.on_edge_leaf_dict.values())}, reticulations are {self.reticulations}")
#         # Add other leaves to edges
#         for edge, leaves in self.on_edge_leaf_dict.items():
#             parent = edge[0]
#             for leaf in leaves:
#                 extra_internal_node, _ = network.add_leaf_to_edge([parent, edge[1]], leaf)
#                 parent = extra_internal_node
#         return network
