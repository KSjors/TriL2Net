import logging
import operator
from datastructures.rooted_level_k_generator import *
from datastructures.rooted_level_k_network import *
from data.all_trinets import *
from data.generators import symmetric_sides_sets_per_generator, generator_dict
from utils.help_functions import *
import copy
from config import settings
import multiprocessing
import time


class Solver:
    def __init__(self
                 , standard_trinet_info_list: NetworkSet
                 , trinet_set: NetworkSet
                 , cut_arc_set_count_method: int = settings.MAXIMUM_MULTIPLICITY
                 , minimal_sink_set_method: int = settings.DEFAULT_SINK_SET
                 , leaf_locator_method: int = settings.GREEDY
                 , level_threshold_method: int = settings.DEFAULT_THRESHOLD
                 , level_count_method: int = settings.WEIGHTED_AVERAGE
                 , generator_count_method: int = settings.WEIGHTED_AVERAGE
                 , symmetric_sides_set_count_method: int = settings.WEIGHTED_AVERAGE
                 , leaf_order_count_method: int = settings.WEIGHTED_AVERAGE
                 , leaf_order_method: int = settings.DEFAULT_ORDER
                 , is_main=True
                 ):
        # Set up logging
        self.uid = guid()
        self.logger_name = f"solver.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)

        # Save config
        self.is_main = is_main
        self.cut_arc_set_count_method = cut_arc_set_count_method  # TODO
        self.minimal_sink_set_method = minimal_sink_set_method
        self.leaf_locator_method = leaf_locator_method
        self.level_threshold_method = level_threshold_method
        self.level_count_method = level_count_method
        self.generator_count_method = generator_count_method  # TODO
        self.symmetric_side_set_count_method = symmetric_sides_set_count_method
        self.leaf_order_count_method = leaf_order_count_method
        self.leaf_order_method = leaf_order_method

        # Set up data
        self.standard_trinet_info_list = standard_trinet_info_list
        self.trinet_set = trinet_set
        self.taxa = trinet_set.represented_leaves()
        self.scores = {}

        # Set up step output array
        self.steps = []

        # Logging
        self.logger.info(f"Created solver for list of {len(self.trinet_set)} trinets representing taxa {self.taxa}")

    def solve(self) -> (RootedLevelKNetwork, float):
        """ Solve """
        self.logger.info('Started solving')
        self.logger.info('Shrinking ...')
        self.trinet_set.add_info(self.standard_trinet_info_list)
        while self.next_shrink():
            continue
        self.logger.info('Expanding ...')
        while self.next_expansion():
            continue
        self.logger.info('Finished solving')
        return self.steps[-1]['network'], self.scores

    def next_shrink(self):
        self.logger.info("Performing the next shrink.")
        leaves_leftover = self.trinet_set.represented_leaves()
        self.logger.info(f"Leaves leftover are {leaves_leftover}")
        assert len({'A', 'B', 'C'}.intersection(leaves_leftover)) == 0, "Network can not contain leaves with names 'A', 'B', or 'C'"
        if len(leaves_leftover) <= 1:
            self.logger.info("Finished shrinking, as there is only one leaf leftover")
            return False
        elif len(leaves_leftover) == 2:
            minimal_sink_sets = [leaves_leftover]
            self.logger.info(f"Inconsistency: Score of minimal sink-set is 100%")
        else:
            auxiliary_graph = Omega.from_network_info_list(self.trinet_set.networks_of_size(3))
            minimal_sink_sets, mss_score = auxiliary_graph.minimal_sink_sets(self.minimal_sink_set_method)
            self.logger.info(f"Inconsistency: Score of minimal sink-set is {mss_score}%")
        self.logger.info(f"Minimal sink-sets are \n {pp.pformat(minimal_sink_sets)}")
        for minimal_sink_set in sorted(minimal_sink_sets, key=lambda x: len(x)):
            self.shrink_mss(minimal_sink_set)
        return True

    def next_expansion(self):
        self.logger.info("Performing the next expansion.")
        steps = self.expand_mss(copy.deepcopy(self.steps[-1]['network']))
        for step in steps:
            self.steps.append(step)
            # if self.is_main:
            #     step['network'].visualize()
        return len(steps) > 0

    def shrink_mss(self, minimal_sink_set):
        self.logger.info(f"Shrinking minimal sink-set {minimal_sink_set}")
        steps = self.compute_network_of_mss(minimal_sink_set)
        for step in steps:
            self.steps.append(step)
            if step['type'] == 'shrink':
                self.trinet_set = NetworkSet.collapse(self.trinet_set, step['leaf_set'])
        self.trinet_set.add_info(self.standard_trinet_info_list)

    def compute_reticulation_thresholds(self, number_of_leaves):
        if self.level_threshold_method == settings.DEFAULT_THRESHOLD:
            return self.compute_reticulation_thresholds_default(number_of_leaves)
        else:
            raise ValueError

    @staticmethod
    def compute_reticulation_thresholds_default(number_of_leaves):
        number_of_possible_trinets = ncr(number_of_leaves, 3)
        trinet_level_2_threshold = 0.5 * min(ncr(number_of_leaves - 1, 2), ncr(number_of_leaves - 2, 1)) / number_of_possible_trinets
        trinet_level_1_threshold = 0.5 * ncr(number_of_leaves - 1, 2) / number_of_possible_trinets
        trinet_level_thresholds = [0, trinet_level_1_threshold, trinet_level_2_threshold]
        return trinet_level_thresholds

    def compute_level_of_network_list(self, network_set: NetworkSet, max_level: int = 2):
        self.logger.info(f"Computing level ...")
        if max_level > 2:
            self.logger.warning("This code is not intended to work for max_level > 2")

        # Compute number of trinets that have a k reticulations
        reticulation_counts = self.count_category(network_set, 'strict_level', self.level_count_method)
        reticulation_percentages = {level: reticulations / len(network_set) for level, reticulations in reticulation_counts.items()}

        # Compute percentage thresholds
        reticulation_thresholds = self.compute_reticulation_thresholds(len(network_set.represented_leaves()))

        # Compute breaches and pick highest breach
        threshold_passed_index = np.array(
            [int(percentage >= threshold) for percentage, threshold in zip(reticulation_percentages.values(), reticulation_thresholds)])
        best_level = max(np.where(threshold_passed_index[:max_level + 1] == 1)[0])

        # Compute scores
        # TODO

        # Logging
        self.logger.info(f"Best level is {best_level}")
        # self.logger.info(f"Inconsistency: score of chosen level is {level_percentages}%")
        # self.logger.info(f"             : score of other levels is \n {pp.pformat(score_dict)}")

        return best_level

    def compute_generator_of_network_list(self, network_set: NetworkSet, level: int) -> RootedLevelKGenerator:
        self.logger.info("Computing generator ...")
        assert level >= 1, "This code only works for level >= 1"

        # Retrieve networks of interest
        network_set = NetworkSet.networks_where(network_set, 'strict_level', level)
        network_set = NetworkSet.networks_with_category(network_set, 'generator_name')

        # Count generators
        generator_name_count = self.weighted_average(network_set, 'generator_name')

        # Compute best generator
        best_generator_name = max(generator_name_count.items(), key=operator.itemgetter(1))[0]
        generator = copy.deepcopy(generator_dict[best_generator_name])

        # Compute scores
        # TODO

        # Logging
        self.logger.info(f"Generator of underlying generator of minimal sink-set is {generator.name}")
        return generator

    def count_leaf_symmetric_side_set(self, network_set: NetworkSet, leaf_set: bidict, generator: RootedLevelKGenerator):
        # Get trinets with information of leaf locations
        network_set = NetworkSet.networks_with_category(network_set, 'generator_name')
        network_set = NetworkSet.networks_where(network_set, 'generator_name', generator.name)

        # Count times leaf is on an edge side
        symmetric_edge_side_sets = list(generator.sets_of_symmetric_edge_sides.keys())
        # leaf_on_symmetric_side_set_count = {leaf: {symmetric_edge_side_set: 0 for symmetric_edge_side_set in symmetric_edge_side_sets} for leaf in
        #                                     leaf_set.values()}
        leaf_on_symmetric_side_set_count = self.count_category(network_set, 'extra_leaf_dict', self.symmetric_side_set_count_method)
        print(leaf_on_symmetric_side_set_count)
        # for name, dictionary in network_set:
        #     volume = dictionary['volume']
        #     for network_info in dictionary['network_info_set']:
        #         extra_leaf_dict = network_info['extra_leaf_dict']
        #         for leaf, symmetric_edge_side_set in extra_leaf_dict.items():
        #             leaf_on_symmetric_side_set_count[leaf][symmetric_edge_side_set] += network_info.multiplicity / volume

        # Put this count in matrix form
        number_of_leaves = len(leaf_set)
        number_of_symmetric_edge_side_sets = len(symmetric_edge_side_sets)
        leaf_on_symmetric_edge_side_set_count_matrix = np.zeros((number_of_symmetric_edge_side_sets, number_of_leaves))
        for leaf_index, leaf in enumerate(leaf_set.values()):
            for symmetric_edge_side_set_index, symmetric_edge_side_set in enumerate(symmetric_edge_side_sets):
                leaf_on_symmetric_edge_side_set_count_matrix[symmetric_edge_side_set_index][leaf_index] = leaf_on_symmetric_side_set_count[leaf][
                    symmetric_edge_side_set]

        # Count how many times leaf is on a reticulation side
        # TODO --> Get from generator by reticulation sides or smthng
        n1 = network_set.any()
        reticulation_sides = list(n1['reticulations'])
        number_of_reticulation_sides = len(reticulation_sides)

        leaf_is_reticulation_side_dict = {leaf: {ret_side: 0 for ret_side in reticulation_sides} for leaf in leaf_set.values()}
        for name, dictionary in network_set:
            volume = dictionary['volume']
            for network_info in dictionary['network_info_set']:
                relation_dict = network_info['relation_dict']
                for generator_leaf, leaf in relation_dict.items():
                    if generator_leaf in reticulation_sides:
                        leaf_is_reticulation_side_dict[leaf][generator_leaf] += network_info.multiplicity / volume

        # Put this count in a matrix
        leaf_on_reticulation_side_count_matrix = np.zeros((number_of_reticulation_sides, number_of_leaves))
        for leaf_index, leaf in enumerate(leaf_set.values()):
            for reticulation_side_index, reticulation_side in enumerate(reticulation_sides):
                leaf_on_reticulation_side_count_matrix[reticulation_side_index][leaf_index] = leaf_is_reticulation_side_dict[leaf][reticulation_side]

        return leaf_on_symmetric_edge_side_set_count_matrix, leaf_on_reticulation_side_count_matrix, symmetric_edge_side_sets, reticulation_sides

    def locate_leaves_ILP(self, leaf_set, leaf_on_edge_side_count_matrix, leaf_on_reticulation_side_count_matrix, edge_sides, reticulation_sides):
        self.logger.info(f"Locating using using ILP")
        number_of_edges, number_of_leaves = leaf_on_edge_side_count_matrix.shape
        number_of_reticulations = leaf_on_reticulation_side_count_matrix.shape[0]

        # Create simplex problem
        # cost array
        cost_matrix = np.vstack((leaf_on_edge_side_count_matrix, leaf_on_reticulation_side_count_matrix))
        for j in range(number_of_leaves):
            norm = sum(cost_matrix[:, j])
            if norm != 0:
                cost_matrix[:, j] /= norm
        cost_array = cost_matrix.flatten()

        # constraint coefficient array
        eye_list = [np.eye(number_of_leaves) for _ in range(number_of_reticulations + number_of_edges)]
        constraint_matrix_p1 = np.hstack(eye_list)
        constraint_matrix_p2 = np.zeros((number_of_reticulations, number_of_leaves * (number_of_reticulations + number_of_edges)))
        for i in range(number_of_edges, number_of_edges + number_of_reticulations):
            for j in range(number_of_leaves * i, number_of_leaves * (i + 1)):
                constraint_matrix_p2[i - number_of_edges][j] = 1
        constraint_matrix = np.vstack((constraint_matrix_p1, constraint_matrix_p2))

        # constraint array
        b_array = [1] * (number_of_leaves + number_of_reticulations)

        # ILP solution
        solution, score, status = simplex_to_ILP(c=cost_array, A_eq=constraint_matrix, b_eq=b_array)
        solution = np.where(np.array(solution) == 1.0)[0]

        # solution back to leaf names
        edge_side_leaf_dict = {edge_side: set() for edge_side in edge_sides}
        leaf_reticulation_side_bidict = bidict()
        for index in solution:
            part = int(index / number_of_leaves)
            if part < number_of_edges:
                edge_side_leaf_dict[edge_sides[part]].add(leaf_set[index % number_of_leaves])
            else:
                leaf_reticulation_side_bidict[leaf_set[index % number_of_leaves]] = reticulation_sides[part - number_of_edges]

        # Score adjustment
        score = int(100 * score / number_of_leaves)
        score_dict = {leaf_set[index % number_of_leaves]: f"{int(100 * cost_array[index])}%" for index in solution}

        # Logging
        self.logger.info(f"Placing leaves on side edges as follows: \n {pp.pformat(edge_side_leaf_dict)}")
        self.logger.info(f"Placing leaves on reticulation sides as follows: \n {pp.pformat(leaf_reticulation_side_bidict)}")
        self.logger.info(f"Inconsistency: ILP score is {score} %")
        self.logger.info(f"             : ILP score per leaf is \n {score_dict} %")
        return edge_side_leaf_dict, leaf_reticulation_side_bidict

    # TODO: naming
    def locate_leaves_greedy(self, leaf_set, leaf_on_edge_side_count_matrix, leaf_on_reticulation_side_count_matrix, edge_sides, reticulation_sides):
        self.logger.info(f"Locating using using normal")

        number_of_reticulations = leaf_on_reticulation_side_count_matrix.shape[0]

        # Pick best leaves to place on reticulation sides
        reticulation_side_permutation_iterator = itertools.permutations(leaf_set.keys(), number_of_reticulations)
        best_score = 0
        best_permutation = None
        for reticulation_sides_permutation in reticulation_side_permutation_iterator:
            score = sum(
                [leaf_on_reticulation_side_count_matrix[reticulation_side_index, leaf_index] for reticulation_side_index, leaf_index in
                 enumerate(reticulation_sides_permutation)])
            if score >= best_score:
                best_permutation = reticulation_sides_permutation
                best_score = score

        # Put in bidict
        score_dict = dict()
        leaf_reticulation_side_bidict = bidict()
        for reticulation_side_index, leaf_index in enumerate(best_permutation):
            leaf_reticulation_side_bidict[leaf_set[leaf_index]] = reticulation_sides[reticulation_side_index]
            score_dict[leaf_set[leaf_index]] = leaf_on_reticulation_side_count_matrix[reticulation_side_index, leaf_index]

        # Place leftover leaves on best edge side
        edge_side_leaf_dict = {edge_side: set() for edge_side in edge_sides}
        for leaf_index in set(leaf_set.keys()).difference(best_permutation):
            edge_side_index = np.argmax(leaf_on_edge_side_count_matrix[:, leaf_index])
            edge_side_leaf_dict[edge_sides[edge_side_index]].add(leaf_set[leaf_index])
            score_dict[leaf_set[leaf_index]] = leaf_on_edge_side_count_matrix[edge_side_index, leaf_index]

        # Normalize score
        score = sum(score_dict.values()) / len(leaf_set)

        # Logging
        self.logger.info(f"Found that leaves are on symmetric edge side sets as follows: \n {pp.pformat(edge_side_leaf_dict)}")
        self.logger.info(f"Found that leaves are on reticulation sides as follows: \n {pp.pformat(leaf_reticulation_side_bidict)}")
        # self.logger.info(f"Inconsistency: score is {score} %")
        # self.logger.info(f"             : score per leaf is \n {score_dict} %")

        return edge_side_leaf_dict, leaf_reticulation_side_bidict

    def fill_symmetric_side_sets(self, network_set: NetworkSet, leaf_set: bidict, generator: RootedLevelKGenerator):
        self.logger.info("Filling symmetric side sets  ...")
        # TODO: For generators with symmetric reticulations leaves score will not be 100%

        # Count occurence of leaves on symmetric side sets
        leaf_on_symmetric_edge_side_set_count_matrix, leaf_on_reticulation_side_count_matrix, symmetric_edge_side_sets, reticulation_sides \
            = self.count_leaf_symmetric_side_set(network_set, leaf_set, generator)

        # Get leaf placement
        if self.leaf_locator_method == 'ILP':
            edge_symmetric_side_set_leaf_dict, leaf_reticulation_side_bidict \
                = self.locate_leaves_ILP(leaf_set, leaf_on_symmetric_edge_side_set_count_matrix, leaf_on_reticulation_side_count_matrix,
                                         symmetric_edge_side_sets, reticulation_sides)
        elif self.leaf_locator_method == 'greedy':
            edge_symmetric_side_set_leaf_dict, leaf_reticulation_side_bidict \
                = self.locate_leaves_greedy(leaf_set, leaf_on_symmetric_edge_side_set_count_matrix, leaf_on_reticulation_side_count_matrix,
                                            symmetric_edge_side_sets, reticulation_sides)
        else:
            raise ValueError('Leaf locator method unknown')

        # Compute score
        # TODO

        # Logging
        self.logger.info(f"Found that the leaves are on the edges as follows \n {pp.pformat(edge_symmetric_side_set_leaf_dict)}")
        self.logger.info(f"Found that the reticulations in the generator should be renamed as follows \n {pp.pformat(leaf_reticulation_side_bidict)}")
        return edge_symmetric_side_set_leaf_dict, leaf_reticulation_side_bidict, reticulation_sides

    def split_symmetric_side_sets(self, network_set, generator, symmetric_edge_side_set_leaf_dict):
        self.logger.info("Splitting symmetric side sets  ...")

        # Get trinets with information on same side information
        network_set = NetworkSet.networks_with_category(network_set, 'strict_level')
        network_set = NetworkSet.networks_where(network_set, 'strict_level', 1) + NetworkSet.networks_where(network_set, 'strict_level', 2)

        # Get leaf order matrix
        leaf_order_matrix, leaf_set = self.compute_leaf_order_matrix(network_set)

        # Create leaf alignment matrix
        leaf_side_alignment_matrix = leaf_order_matrix + leaf_order_matrix.T

        # Retrieve symmetric side set info
        symmetric_sides_sets = symmetric_sides_sets_per_generator[generator.name]

        # Split each group of leaves assigned to a set of symmetric sides into groups of leaves per side
        leaves_per_side_per_symmetric_side = {symmetric_edge_side_set: None for symmetric_edge_side_set in symmetric_edge_side_set_leaf_dict}
        for symmetric_edge_side_set, leaves in symmetric_edge_side_set_leaf_dict.items():
            # Create empty list of leaves per side in set of symmetric sides
            leaves_per_side = {side_index: set() for side_index, _ in enumerate(symmetric_sides_sets[symmetric_edge_side_set])}

            # If only one side add all leaves to it
            if len(leaves_per_side.keys()) == 1:
                leaves_per_side[0] = [leaf_set.inverse[leaf_name] for leaf_name in leaves]

            # Else compute alignment scores and add leaf to best side
            else:
                for leaf_name in leaves:
                    leaf_index = leaf_set.inverse[leaf_name]

                    # Compute side alignment score
                    side_alignment_score = {side_index: sum(leaf_side_alignment_matrix[leaf_index, list(placed_leaf_indicis)]) for
                                            side_index, placed_leaf_indicis in
                                            leaves_per_side.items()}

                    # Normalize this score
                    normalized_side_alignment_score = {}
                    for side_index_0, score_0 in side_alignment_score.items():
                        temp = 2 * score_0 - len(leaves_per_side[side_index_0])
                        for side_index_1, score_1 in side_alignment_score.items():
                            temp -= score_1 - len(leaves_per_side[side_index_1])
                        normalized_side_alignment_score[side_index_0] = temp

                    # Find best score
                    best_side = max(normalized_side_alignment_score.items(), key=operator.itemgetter(1))[0]
                    leaves_per_side[best_side].add(leaf_index)

            # Save split groups
            leaves_per_side_per_symmetric_side[symmetric_edge_side_set] = {side: {leaf_set[leaf_index] for leaf_index in leaves} for side, leaves in
                                                                           leaves_per_side.items()}

        self.logger.info(f"Found that the leaves are on the edge sides as follows \n {pp.pformat(leaves_per_side_per_symmetric_side)}")

        return leaves_per_side_per_symmetric_side, leaf_order_matrix, leaf_set

    def compute_network_of_mss(self, minimal_sink_set) -> list:
        self.logger.info(f"Computing network of minimal sink-set {minimal_sink_set}")

        if len(minimal_sink_set) == 2:
            binet_set = NetworkSet.induced_network_set_of_network_set(self.trinet_set, 2, max_processes=1, progress_bar=False)
            restricted_binet_set = NetworkSet.restrict(binet_set, minimal_sink_set)
            highest_multiplicity = -np.inf
            best_binet = None
            for network_info in restricted_binet_set.per_network_info():
                multiplicity = network_info.multiplicity
                if multiplicity > highest_multiplicity:
                    highest_multiplicity = multiplicity
                    best_binet = copy.deepcopy(network_info.network)
            step = [{'network': best_binet, 'leaf_set': set(minimal_sink_set), 'type': 'shrink', 'name': mss_leaf_name(minimal_sink_set)}]
            return step

        restricted_trinet_set = NetworkSet.restrict(self.trinet_set, minimal_sink_set)
        restricted_trinet_set.add_info(self.standard_trinet_info_list)
        self.logger.info(f"Number of networks that describe this mss is {len(restricted_trinet_set)}")

        leaf_set = bidict({index: leaf for index, leaf in enumerate(minimal_sink_set)})

        # Compute level
        level = self.compute_level_of_network_list(restricted_trinet_set)

        # In case level is zero pick most found binet/trinet
        if level == 0:
            best_network = None
            if len(minimal_sink_set) == 3:
                restricted_trinet_set_t = NetworkSet.induced_network_set_of_network_set(restricted_trinet_set, len(minimal_sink_set))
                for network_info in restricted_trinet_set_t.maximum_multiplicity_iterator():
                    best_network = copy.deepcopy(network_info.network)
            else:
                if len(minimal_sink_set) == len(self.taxa):
                    # TODO: Check if this can happen (in case of ties between trinets describing he same leaf sets) and do something
                    raise NotImplementedError
                sub_solver = Solver(self.standard_trinet_info_list, restricted_trinet_set, self.leaf_locator_method, is_main=False)
                sub_solver.solve()
                return sub_solver.steps
            steps = [{'network': best_network, 'leaf_set': set(minimal_sink_set), 'type': 'shrink', 'name': mss_leaf_name(minimal_sink_set)}]
            return steps

        # Compute generator
        generator = self.compute_generator_of_network_list(restricted_trinet_set, level)
        network = RootedLevelKNetwork.from_generator(generator)

        # Leaves per symmetric side set
        edge_symmetric_side_set_leaf_dict, leaf_reticulation_side_bidict, reticulation_sides \
            = self.fill_symmetric_side_sets(restricted_trinet_set, leaf_set, generator)

        # Choose parent name for each reticulation leaf and rename reticulations in network
        leaf_below_reticulation_side_bidict = {}
        for reticulation_leaf, old_reticulation_leaf in leaf_reticulation_side_bidict.items():
            leaf_below_reticulation_side_bidict[reticulation_leaf] = network.nodes_above_nodes({old_reticulation_leaf}, max_height=1).difference(
                {old_reticulation_leaf}).pop()
            network.rename_node(new_name=reticulation_leaf, old_name=old_reticulation_leaf)

        # Split symmetric side sets into set of leaves per side
        leaves_per_side_per_symmetric_side, leaf_order_matrix, leaf_set = self.split_symmetric_side_sets(restricted_trinet_set, generator,
                                                                                                         edge_symmetric_side_set_leaf_dict)

        # Order leaves per side
        ordered_leaves_per_side_per_symmetric_side = self.compute_leaf_order(leaves_per_side_per_symmetric_side, leaf_order_matrix, leaf_set)

        # In case of generator 2c, need to align the sides in different symmetric side sets
        if generator.name == '2c':
            ordered_leaves_per_side_per_symmetric_side = self.align_sides(ordered_leaves_per_side_per_symmetric_side, leaf_order_matrix)
            ordered_leaves_per_side_per_symmetric_side = self.align_reticulations(ordered_leaves_per_side_per_symmetric_side,
                                                                                  leaf_below_reticulation_side_bidict, leaf_order_matrix, leaf_set)

        # Split symmetric sides in to sides
        ordered_leaves_per_side = self.name_sides_in_symmetric_side_sets(generator, ordered_leaves_per_side_per_symmetric_side)

        # Add leaves to edges
        self.add_leaves_to_edges(network, ordered_leaves_per_side)
        step = [{'network': network, 'leaf_set': set(minimal_sink_set), 'type': 'shrink', 'name': mss_leaf_name(minimal_sink_set)}]
        return step

    def name_sides_in_symmetric_side_sets(self, generator, ordered_leaves_per_side_per_symmetric_side):
        self.logger.info("Naming sides in symmetric side sets ...")
        symmetric_side_sets = generator.sets_of_symmetric_sides
        ordered_leaves_per_side = {}
        for symmetric_side_set, leaves_per_side in ordered_leaves_per_side_per_symmetric_side.items():
            sides = symmetric_side_sets[symmetric_side_set]
            for side_index, leaves in leaves_per_side.items():
                ordered_leaves_per_side[sides[side_index]] = leaves
        self.logger.info(f"Found that leaves should be on the sides as follows \n {pp.pformat(ordered_leaves_per_side)}")
        return ordered_leaves_per_side

    def order_leaves(self, sub_leaf_order_matrix, method=settings.DEFAULT_ORDER):
        if method == settings.DEFAULT_ORDER:
            return self.order_leaves_default(sub_leaf_order_matrix)
        else:
            raise ValueError

    @staticmethod
    def order_leaves_default(sub_leaf_order_matrix):
        leaf_align_matrix = sub_leaf_order_matrix + sub_leaf_order_matrix.T
        leaf_align_matrix[leaf_align_matrix == 0] = 1
        normalized_leaf_order_matrix = np.divide(sub_leaf_order_matrix - sub_leaf_order_matrix.T, leaf_align_matrix)

        n = normalized_leaf_order_matrix.shape[0]
        ordered_leaves = []
        for _ in range(n):
            out_degree = np.sum(normalized_leaf_order_matrix, axis=1)
            out_degree[ordered_leaves] = -np.inf
            best_leaf = int(np.argmax(out_degree))
            ordered_leaves.append(best_leaf)
        return ordered_leaves

    def compute_leaf_order(self, symmetric_edge_side_set_leaf_dict_per_side: dict, leaf_order_matrix, leaf_set: bidict):
        self.logger.info(f"Computing leaf order for side sets ...")

        symmetric_edge_side_set_ordered_leaf_dict_per_side = {}
        for symmetric_edge_side_set, leaves_per_side in symmetric_edge_side_set_leaf_dict_per_side.items():
            symmetric_edge_side_set_ordered_leaf_dict_per_side[symmetric_edge_side_set] = []
            for _, leaf_names in leaves_per_side.items():
                # Only look at leaves on this side
                leaf_indicis = sorted([leaf_set.inverse[leaf_name] for leaf_name in leaf_names])
                sub_leaf_order_matrix = leaf_order_matrix[leaf_indicis, :][:, leaf_indicis]

                # Order leaves
                ordered_leaves = self.order_leaves(sub_leaf_order_matrix, method=self.leaf_order_method)

                # Reindex leaves
                ordered_leaves = [leaf_set[leaf_indicis[leaf_index]] for leaf_index in ordered_leaves]

                # Save
                symmetric_edge_side_set_ordered_leaf_dict_per_side[symmetric_edge_side_set].append(ordered_leaves)

        self.logger.info(f"Found that the leaves are on the edges as follows \n {pp.pformat(symmetric_edge_side_set_ordered_leaf_dict_per_side)}")
        return symmetric_edge_side_set_ordered_leaf_dict_per_side

    def compute_leaf_order_matrix(self, network_set):
        # Count leaf orders
        leaf_orderings_count = self.count_category(network_set, 'leaf_order', self.leaf_order_count_method)

        # Count to matrix
        leaf_set = bidict({index: leaf_name for index, leaf_name in enumerate(self.trinet_set.represented_leaves())})
        leaf_order_matrix = np.zeros((len(leaf_set), len(leaf_set)))
        for leaf_orderings, weight in leaf_orderings_count.items():
            for leaf_order in leaf_orderings:
                for high_leaf_index in range(len(leaf_order)):
                    for low_leaf_index in range(high_leaf_index + 1, len(leaf_order)):
                        low_leaf = leaf_set.inverse[leaf_order[low_leaf_index]]
                        high_leaf = leaf_set.inverse[leaf_order[high_leaf_index]]
                        leaf_order_matrix[high_leaf, low_leaf] += weight
        return leaf_order_matrix, leaf_set

    @staticmethod
    def add_leaves_to_edges(network, ordered_leaves_per_side):
        for side, ordered_leaves in ordered_leaves_per_side.items():
            from_node = side[0]
            to_node = side[1]
            for leaf in ordered_leaves:
                internal_name, _ = network.add_leaf_to_edge(edge=[from_node, to_node], leaf_name=leaf)
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

    def align_reticulations(self, ordered_leaves_per_side_per_symmetric_side, leaf_below_reticulation_side_bidict, leaf_order_matrix, leaf_set):
        self.logger.info("Aligning reticulation sides to symmetric side sets")

        # TODO: put in thesis

        # Compute how often leaves of symmetric side set are above a reticulation
        reticulation_indicis = {leaf_set.inverse[reticulation_leaf]: reticulation for reticulation_leaf, reticulation in
                                leaf_below_reticulation_side_bidict.items()}
        symmetric_side_set_above_reticulation_score = {}
        for symmetric_sides_set, ordered_leaves_per_side in ordered_leaves_per_side_per_symmetric_side.items():
            leaf_lists = list(ordered_leaves_per_side.values())
            leaf_indicis = [leaf_set.inverse[leaf_name] for leaf_name in itertools.chain.from_iterable(leaf_lists)]
            symmetric_side_set_above_reticulation_score[symmetric_sides_set] = {reticulation: sum(leaf_order_matrix[leaf_indicis, reticulation_index]) for
                                                                                reticulation_index, reticulation in reticulation_indicis.items()}

        # Check if symmetric side set (1,3) and (1,4) need to be swapped
        left_on_first_ret_right_on_second_ret = symmetric_side_set_above_reticulation_score[('1', '3')]['3'] + \
                                                symmetric_side_set_above_reticulation_score[('1', '4')]['4']
        left_on_second_ret_right_on_first_ret = symmetric_side_set_above_reticulation_score[('1', '3')]['3'] + \
                                                symmetric_side_set_above_reticulation_score[('1', '4')]['4']
        if left_on_first_ret_right_on_second_ret < left_on_second_ret_right_on_first_ret:
            ordered_leaves_per_side_per_symmetric_side[('1', '3')], ordered_leaves_per_side_per_symmetric_side[('1', '4')] = \
                ordered_leaves_per_side_per_symmetric_side[('1', '4')], ordered_leaves_per_side_per_symmetric_side[('1', '3')]

        self.logger.info(f"Found that the leaves are on the edges as follows \n {pp.pformat(ordered_leaves_per_side_per_symmetric_side)}")
        return ordered_leaves_per_side_per_symmetric_side

    def align_sides(self, ordered_leaves_per_side_per_symmetric_side, leaf_order_matrix):
        leaf_alignment_matrix = leaf_order_matrix + leaf_order_matrix.T

        # Set up alignment score structure
        alignment_scores = np.zeros((6, 6))

        # Number symmetric side sets
        symmetric_side_set_ordering = list(ordered_leaves_per_side_per_symmetric_side.keys())

        # Compute alignment scores
        symmetric_side_set_combination_iterator = itertools.combinations([0, 1, 2], 2)

        for symmetric_side_set_combination in symmetric_side_set_combination_iterator:
            symmetric_side_set_index_0 = symmetric_side_set_ordering[symmetric_side_set_combination[0]]
            symmetric_side_set_index_1 = symmetric_side_set_ordering[symmetric_side_set_combination[1]]
            set_of_symmetric_sides_0 = ordered_leaves_per_side_per_symmetric_side[symmetric_side_set_index_0]
            set_of_symmetric_sides_1 = ordered_leaves_per_side_per_symmetric_side[symmetric_side_set_index_1]
            for side_0_index, side_0 in enumerate(set_of_symmetric_sides_0):
                for side_1_index, side_1 in enumerate(set_of_symmetric_sides_1):
                    u = sum([leaf_alignment_matrix[leaf_index_0][leaf_index_1] for leaf_index_0 in side_0 for leaf_index_1 in side_1]) - len(side_0) * len(
                        side_1)
                    matrix_index_0 = symmetric_side_set_combination[0] * 2 + side_0_index
                    matrix_index_1 = symmetric_side_set_combination[1] * 2 + side_1_index
                    alignment_scores[matrix_index_0][matrix_index_1] = u

        # Iterate through different configuration and pick best
        configurations = [[0, 1]] * 3
        configuration_iterator = itertools.product(*configurations)
        best_configuration = None
        best_score = - np.inf
        for configuration in configuration_iterator:
            two_symmetric_side_sets_iterator = itertools.combinations([0, 1, 2], 2)
            score = 0
            for two_symmetric_side_sets in two_symmetric_side_sets_iterator:
                matrix_index_0 = two_symmetric_side_sets[0] * 2 + configuration[0]
                matrix_index_1 = two_symmetric_side_sets[1] * 2 + configuration[1]
                score += alignment_scores[matrix_index_0][matrix_index_1]
            if score > best_score:
                best_score = score
                best_configuration = configuration

        # Apply best configuration
        aligned_ordered_leaves_per_side_per_symmetric_side = {}
        for symmetric_side_set_index, symmetric_side_set in enumerate(symmetric_side_set_ordering):
            left_side = best_configuration[symmetric_side_set_index]
            right_side = (best_configuration[symmetric_side_set_index] + 1) % 2
            symmetric_side_set_leaves_per_side = ordered_leaves_per_side_per_symmetric_side[symmetric_side_set]
            aligned_ordered_leaves_per_side_per_symmetric_side[symmetric_side_set] = {0: symmetric_side_set_leaves_per_side[left_side],
                                                                                      1: symmetric_side_set_leaves_per_side[right_side]}

        self.logger.info(f"Found that the leaves are on the edges as follows \n {pp.pformat(aligned_ordered_leaves_per_side_per_symmetric_side)}")
        return aligned_ordered_leaves_per_side_per_symmetric_side

    def expand_mss(self, biconnected_component):
        self.logger.info(f"Looking for leaves in component {biconnected_component.uid} to expand.")
        replace_dict = {}
        all_taxa = biconnected_component.leaf_names
        self.logger.info(f"Component has taxa {all_taxa}.")
        for leaf_name in set(biconnected_component.leaf_names).difference(self.taxa):
            self.logger.info(f"Need to replace {leaf_name}")
            for step in self.steps:
                if step['name'] == leaf_name:
                    replace_dict[leaf_name] = step['network']
                    break
        if replace_dict:
            for leaf_name, component in replace_dict.items():
                self.logger.info(f"Expanding {leaf_name} using {component.uid}")
                biconnected_component.replace_leaf_with_network(leaf_name, component)
                biconnected_component.standardize_internal_node_names()
            steps = [{'network': biconnected_component, 'leaf_set': set(all_taxa), 'type': 'expand', 'name': mss_leaf_name(all_taxa)}]
            return steps
        else:
            self.logger.info("No leaves found to expand.")
            return []

    # ------------------------------- COUNTING METHODS ----------------------------------- #
    @staticmethod
    def count_category(network_set: NetworkSet, category: str, method: int = settings.MAXIMUM_MULTIPLICITY):
        result = {}
        for network_info, w in network_set.per_network_iterator(method=method):
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

    def __getstate__(self):
        self.logger_name = 'solver.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger_name)
        return self.__dict__
