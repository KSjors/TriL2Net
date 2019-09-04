import logging
# from datastructures.omega import *
from datastructures.rooted_level_k_network import *
from data.all_trinets import *
from utils.help_functions import *
import copy
import multiprocessing
import time


class Solver:
    def __init__(self, standard_trinet_info_list, trinet_info_list: TrinetInfoList):
        self.uid = guid()
        self.logger_name = f"solver.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)
        self.standard_trinet_info_list = standard_trinet_info_list
        self.trinet_info_list = trinet_info_list
        self.trinet_info_list.calculate_info()
        self.last_component = None
        self.transformations = {}
        self.components = {}

    def solve(self):
        self.logger.info('Started solving')
        self.logger.info('Shrinking ...')
        while self.next_shrink():
            continue
        self.logger.info('Expanding ...')
        while self.expand_mss(self.last_component):
            self.last_component.visualize()
            continue
        self.logger.info('Finished solving')
        self.last_component.visualize()
        return self.last_component

    def next_shrink(self, max_processes=1):
        self.logger.info("Performing the next shrink.")
        leaves_leftover = self.trinet_info_list.represented_leaves()
        self.logger.info(f"Leaves leftover are {leaves_leftover}")
        assert len({'A', 'B', 'C'}.intersection(leaves_leftover)) == 0, "Network can not contain leaves with names 'A', 'B', or 'C'"
        if len(leaves_leftover) <= 1:
            self.logger.info("Finished shrinking, as there is only one leaf leftover")
            return False
        elif len(leaves_leftover) == 2:
            minimal_sink_sets = [leaves_leftover]
        else:
            minimal_sink_sets = self.trinet_info_list.get_minimal_sink_sets(level=0)
        self.logger.info(f"Minimal sink sets are {minimal_sink_sets}")
        for minimal_sink_set in sorted(minimal_sink_sets, key=lambda x: len(x)):
            self.shrink_mss(minimal_sink_set)
        return True

    def shrink_mss(self, minimal_sink_set):
        self.logger.info(f"Shrinking minimal sink-set {minimal_sink_set}")

        # Get underlying structure and level
        # Get trinets which describe sink set
        trinet_info_list_mss = self.trinet_info_list.find_trinets_with_leaf_names_in(minimal_sink_set)
        self.logger.info(f"Number of trinets for this mss is {len(trinet_info_list_mss)}")
        trinet_info_list_mss.add_info(self.standard_trinet_info_list)
        if len(minimal_sink_set) == 2:
            trinet_info_list_mss.calculate_info()

        # Get trinets with information
        trinet_info_list_mss = trinet_info_list_mss.trinets_with_category('relation_dict')

        # Find level and generator
        level = trinet_info_list_mss.max_level()
        trinet_info_list_mss_gen_level = trinet_info_list_mss.trinets_where('level', level)

        self.logger.info(f"Level of underlying generator of minimal sink-set is {level}")

        if level == 0:
            generator = copy.deepcopy(trinet_info_list_mss_gen_level[0].trinet)
            self.logger.info(f"Generator of underlying generator of minimal sink-set is 0")
        elif level >= 1:
            generator = trinet_info_list_mss_gen_level.best_generator()
            self.logger.info(f"Generator of underlying generator of minimal sink-set is {generator.name}")

            # Get information about placement of leaves on sides
            # Get trinets with this generator and find on which side each leaf is
            trinet_info_list_mss_gen_level_gen = trinet_info_list_mss_gen_level.trinets_where('generator_name', generator.name)
            edge_leaf_dict, reticulation_names = trinet_info_list_mss_gen_level_gen.get_leaf_locations()

            self.logger.info(f"Found that the leaves are on the edges as follows {pp.pformat(edge_leaf_dict)}")
            self.logger.info(f"Found that the reticulations in the generator should be renamed as follows {pp.pformat(reticulation_names)}")

            # Get trinets of level-1 for ordering of leaf on sides
            trinet_info_list_mss_level_1_2 = trinet_info_list_mss.trinets_where('level', 1) + trinet_info_list_mss.trinets_where('level', 2)
            leaf_set, leaf_order_matrix_1 = trinet_info_list_mss_level_1_2.leaf_order_info()

            # Remove symmetries
            sides_per_edge = {}
            rets = list(reticulation_names.keys())
            if generator.name == '1':
                sides_per_edge = {
                    ('0', '1'): [('0', '1'), ('0', '1')]
                }
            elif generator.name == '2a':
                sides_per_edge = {
                    ('0', '1'): [('0', '1')],
                    ('0', '2'): [('0', '2')],
                    ('1', '2'): [('1', '2')],
                    ('1', '3'): [('1', '3')],
                    ('2', '3'): [('2', '3')]
                }
            elif generator.name == '2b':
                sides_per_edge = {
                    ('0', '1'): [('0', '1')],
                    ('0', '2'): [('0', '2')],
                    ('1', '3'): [('1', '3')],
                    ('1', '4'): [('1', '4')],
                    ('3', '2'): [('3', '2')],
                    ('3', '4'): [('3', '4')]
                }
            elif generator.name == '2c':
                # Nodes (1,2), (3,4) symmetrical. Need to split edge (1,3) into two sets, one for each reticulation
                edge_1_3_leaves = edge_leaf_dict[('1', '3')]
                edge_1_3_leaf_indicis = sorted([leaf_set.inverse[leaf_name] for leaf_name in edge_1_3_leaves])
                leaf_set_2, leaf_order_matrix_2 = trinet_info_list_mss_gen_level.leaf_order_info()
                x = leaf_set_2.inverse[rets[0]]
                y = leaf_set_2.inverse[rets[1]]
                leaf_per_reticulation = (np.array([np.argmax(row) for row in leaf_order_matrix_2[edge_1_3_leaf_indicis][:, [x, y]]]))
                leaves_on_ret_0 = set([leaf_set[edge_1_3_leaf_indicis[leaf_index]] for leaf_index in np.where(leaf_per_reticulation == 0.)[0]])
                leaves_on_ret_1 = set([leaf_set[edge_1_3_leaf_indicis[leaf_index]] for leaf_index in np.where(leaf_per_reticulation == 1.)[0]])
                edge_leaf_dict[('1', '3')] = leaves_on_ret_0
                edge_leaf_dict[('1', '4')] = leaves_on_ret_1
                sides_per_edge = {
                    ('0', '1'): [('0', '1'), ('0', '2')],
                    ('1', '3'): [('1', '3'), ('2', '3')],
                    ('1', '4'): [('1', '4'), ('2', '4')]
                }
            elif generator.name == '2d':
                sides_per_edge = {
                    ('0', '1'): [('0', '1')],
                    ('0', '2'): [('0', '2')],
                    ('1', '3'): [('1', '3'), ('1', '3')],
                    ('3', '2'): [('3', '2')]
                }

            for new_name in rets:
                old_name = reticulation_names[new_name]
                generator.rename_node(new_name=new_name, old_name=old_name)

            if generator.name in ('1', '2a', '2b', '2d'):
                for edge, leaves in edge_leaf_dict.items():
                    self.logger.info(f"Ordering leaves {leaves} on side {edge}")
                    sides = len(sides_per_edge[edge])
                    sub_leaf_order_matrix, sub_leaf_set = self.sub_leaf_order_matrix(leaf_set, leaf_order_matrix_1, leaves)
                    ordered_leaves = self.order_leaves(sub_leaf_set, sub_leaf_order_matrix, sides=sides)
                    self.logger.info(f"Leaf order is {ordered_leaves}")
                    for side_number, side_edge in enumerate(sides_per_edge[edge]):
                        from_node = edge[0]
                        for leaf in ordered_leaves[side_number]:
                            internal_name, leaf_name = generator.add_leaf_to_edge([from_node, edge[1]], leaf)
                            from_node = internal_name
            if generator.name == '2c':
                ret_0_parent = generator.get_parents({rets[0]}, max_height=1)
                ret_0_parent.remove(rets[0])
                ret_0_parent = ret_0_parent.pop()
                ret_1_parent = generator.get_parents({rets[1]}, max_height=1)
                ret_1_parent.remove(rets[1])
                ret_1_parent = ret_1_parent.pop()
                edges = [[('0', '1'), ('1', ret_0_parent), ('1', ret_1_parent)], [('0', '2'), ('2', ret_0_parent), ('2', ret_1_parent)]]
                leaves_on_edges = [[[], [], []], [[], [], []]]

                side_orders = []
                edge_side_leaf_dict = dict()
                for index, edge in enumerate(edges[0]):
                    leaves = edge_leaf_dict[edge]
                    sides = len(sides_per_edge[edge])
                    sub_leaf_order_matrix, sub_leaf_set = self.sub_leaf_order_matrix(leaf_set, leaf_order_matrix_1, leaves)
                    ordered_leaves = self.order_leaves(sub_leaf_set, sub_leaf_order_matrix, sides=sides)
                    side_orders.append(ordered_leaves)
                    edge_side_leaf_dict[edge] = ordered_leaves
                    for row_index, row in enumerate(ordered_leaves):
                        leaves_on_edges[row_index][index] = row

                side_separation = self.align_sides(side_orders, leaf_set, leaf_order_matrix_1)

                leaves_on_edges_ordered = copy.deepcopy(leaves_on_edges)
                for leg_index, leg in enumerate(leaves_on_edges[0]):
                    if set(leg).issubset(side_separation[1]):
                        leaves_on_edges_ordered[0][leg_index], leaves_on_edges_ordered[1][leg_index] = leaves_on_edges_ordered[1][leg_index], \
                                                                                                       leaves_on_edges_ordered[0][leg_index]

                ret_0_left_alignment = len(set(leaves_on_ret_0).intersection(leaves_on_edges_ordered[0][1] + leaves_on_edges_ordered[1][1]))
                ret_0_right_alignment = len(set(leaves_on_ret_0).intersection(leaves_on_edges_ordered[0][2] + leaves_on_edges_ordered[1][2]))
                ret_1_left_alignment = len(set(leaves_on_ret_1).intersection(leaves_on_edges_ordered[0][1] + leaves_on_edges_ordered[1][1]))
                ret_1_right_alignment = len(set(leaves_on_ret_1).intersection(leaves_on_edges_ordered[0][2] + leaves_on_edges_ordered[1][2]))
                if ret_0_right_alignment + ret_1_left_alignment > ret_0_left_alignment + ret_1_right_alignment:
                    leaves_on_edges_ordered[0][1], leaves_on_edges_ordered[0][2] = leaves_on_edges_ordered[0][2], leaves_on_edges_ordered[0][1]
                    leaves_on_edges_ordered[1][1], leaves_on_edges_ordered[1][2] = leaves_on_edges_ordered[1][2], leaves_on_edges_ordered[1][1]

                for i in range(len(leaves_on_edges_ordered)):
                    for j in range(len(leaves_on_edges_ordered[i])):
                        from_node = edges[i][j][0]
                        for leaf in leaves_on_edges_ordered[i][j]:
                            internal_name, leaf_name = generator.add_leaf_to_edge((from_node, edges[i][j][1]), leaf)
                            from_node = internal_name

        mss_name = mss_leaf_name(minimal_sink_set)
        generator.visualize()
        self.last_component = generator
        self.components[mss_name] = generator
        self.transformations[mss_name] = minimal_sink_set
        self.trinet_info_list.shrink(minimal_sink_set)
        self.trinet_info_list.calculate_info()

    def next_expansion(self):
        self.logger.info("Performing the next expansion.")
        return self.expand_mss(self.last_component)

    @staticmethod
    def add_leaves(generator, edge_side_leaf_dict, reticulation_names, side_orders, side_order_alignment):
        for edge, sides in edge_side_leaf_dict.items():
            for side in sides:
                for leaf in side:
                    internal_name, leaf_name = generator.add_leaf_to_edge(edge, leaf)
                    from_node = internal_name

    @staticmethod
    def sub_leaf_order_matrix(leaf_set, leaf_order_matrix, leaves):
        leaf_indicis = sorted([leaf_set.inverse[leaf] for leaf in leaves])

        sub_matrix = copy.deepcopy(leaf_order_matrix)
        sub_matrix = sub_matrix[leaf_indicis, :][:, leaf_indicis]

        sub_leaf_set = bidict()
        for new_index, leaf_index in enumerate(leaf_indicis):
            sub_leaf_set[new_index] = leaf_set[leaf_index]
        return sub_matrix, sub_leaf_set

    @staticmethod
    def order_leaves(leaf_set, leaf_order_matrix, sides=1):
        assert sides in (1, 2), "Only works for 1 or 2 sides"
        side_orders = ([[] for side in range(sides)])
        placed_leaves = []
        # place all leaves which have children or parents on the same side
        for _ in [0, 1]:
            # 0 --> do all leaves with children
            # 1 --> do leftover leaves with parents
            row_sum = leaf_order_matrix.sum(axis=1)
            row_sum = [s if index not in placed_leaves else 0 for index, s in enumerate(row_sum)]
            while sum(row_sum) != 0:
                # leaf with maximum out-degree
                current_leaf = np.argmax(row_sum)
                current_leaf_children = set(list(np.where(leaf_order_matrix[current_leaf] >= 1)[0]) + [current_leaf])
                L = 0.5 * len(current_leaf_children)

                # Find which side it is on
                side_alignment = np.zeros(sides)
                for side, side_order in enumerate(side_orders):
                    for leaf in side_order:
                        leaf_children = set(list(np.where(leaf_order_matrix[leaf] >= 1)[0]) + [leaf])
                        same_children = leaf_children.intersection(current_leaf_children)
                        if len(same_children) >= L:
                            side_alignment[side] += 1
                side_alignment_1 = np.zeros(sides)
                divide_boolean = True
                for side in range(sides):
                    side_alignment_1[side] += (len(placed_leaves) - len(side_orders[side])) \
                                              - (sum(side_alignment) - side_alignment[side])
                    if len(side_orders[side]) == 0:
                        divide_boolean = False

                if divide_boolean:
                    for side in range(sides):
                        side_alignment_1 /= len(side_orders[side])

                best_side = np.argmax(side_alignment_1)

                side_orders[best_side].append(current_leaf)
                placed_leaves.append(current_leaf)

                row_sum = leaf_order_matrix.sum(axis=1)
                row_sum = [s if index not in placed_leaves else 0 for index, s in enumerate(row_sum)]
            leaf_order_matrix = leaf_order_matrix.T

        if len(placed_leaves) < len(leaf_set.keys()):
            for leaf in leaf_set.keys():
                for side_order in side_orders:
                    if len(side_order) == 0:
                        side_order.append(leaf)
                        break
        return [[leaf_set[leaf] for leaf in side_order] for side_order in side_orders]

    def align_sides(self, side_orders, leaf_set, leaf_order_matrix):
        side_set = bidict()
        index = 0
        active_sides = []
        for side_order in side_orders:
            for side in side_order:
                if len(side) != 0:
                    side_set.put(index, index)
                    active_sides.append(side)
                    index += 1

        side_order_matrix = np.array([])
        for side in active_sides:
            leaf_indicis_of_side = [leaf_set.inverse[leaf] for leaf in side]
            new_row = np.sum((leaf_order_matrix[leaf_indicis_of_side]), axis=0)
            try:
                side_order_matrix = np.vstack([side_order_matrix, new_row])
            except ValueError:
                side_order_matrix = new_row

        side_order_matrix_2 = np.array([])
        for side in active_sides:
            leaf_indicis_of_side = [leaf_set.inverse[leaf] for leaf in side]
            new_col = np.sum((side_order_matrix[:, leaf_indicis_of_side]), axis=1)
            try:
                side_order_matrix_2 = np.vstack([side_order_matrix_2, new_col])
            except ValueError:
                side_order_matrix_2 = new_col

        for i in range(index):
            side_order_matrix_2[i, i] = 0

        order = self.order_leaves(side_set, side_order_matrix_2.T, sides=2)

        result = []
        for i in range(2):
            result.append([])
            for o in order[i]:
                result[-1] += active_sides[o]

        return result

    def expand_mss(self, biconnected_component):
        self.logger.info(f"Looking for leaves in component {biconnected_component.uid} to expand.")
        replace_dict = {}
        all_taxa = biconnected_component.leaf_names
        self.logger.info(f"Component has taxa {all_taxa}.")
        for leaf_name in biconnected_component.leaf_names:
            if '(' in leaf_name:
                self.logger.info(f"Need to replace {leaf_name}")
                replace_dict[leaf_name] = self.components[leaf_name]
        if replace_dict:
            for leaf_name, component in replace_dict.items():
                self.logger.info(f"Expanding {leaf_name} using {component.uid}")
                biconnected_component.replace_leaf_with_network(leaf_name, component)
                biconnected_component.standardize_internal_node_names()
                self.components.pop(leaf_name)
            return True
        else:
            self.logger.info("No leaves found to expand.")
            return False

    def __getstate__(self):
        self.logger_name = 'solver.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger_name)
        return self.__dict__
