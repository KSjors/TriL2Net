import logging
from datastructures.rooted_level_k_generator import *
from datastructures.rooted_level_k_network import *
from data.all_trinets import *
from data.generators import sides_per_edge_per_generator
from utils.help_functions import *
import copy
import multiprocessing
import time


class Solver:
    def __init__(self, standard_trinet_info_list, trinet_info_list: NetworkInfoList, leaf_locator: str = 'normal', is_main=True):
        assert leaf_locator in ('normal', 'ILP')
        self.uid = guid()
        self.logger_name = f"solver.{self.uid}"
        self.is_main = is_main
        self.logger = logging.getLogger(self.logger_name)
        self.leaf_locator = leaf_locator
        self.standard_trinet_info_list = standard_trinet_info_list
        self.trinet_info_list = trinet_info_list

        self.taxa = trinet_info_list.represented_leaves()
        self.steps = []

        self.logger.info(f"Created solver for list of {len(self.trinet_info_list)} trinets representing taxa {self.taxa}")

    def solve(self) -> RootedLevelKNetwork:
        self.logger.info('Started solving')
        self.logger.info('Shrinking ...')
        self.trinet_info_list.add_info(self.standard_trinet_info_list)
        while self.next_shrink():
            continue
        self.logger.info('Expanding ...')
        while self.next_expansion():
            continue
        self.logger.info('Finished solving')
        return self.steps[-1]['network']

    def next_shrink(self):
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
            auxiliary_graph = Omega.from_network_info_list(self.trinet_info_list.networks_of_size(3))
            minimal_sink_sets, mss_score = auxiliary_graph.minimal_sink_sets()
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
            # if self.is_main:
            #     step['network'].visualize()
            if step['type'] == 'shrink':
                self.trinet_info_list = self.trinet_info_list.shrink(step['leaf_set'])
        self.trinet_info_list.add_info(self.standard_trinet_info_list)

    def compute_level_of_network_list(self, network_list: NetworkInfoList, max_level: int = 2):
        self.logger.info(f"Computing level ...")
        if max_level > 2:
            self.logger.warning("This code is not intended to work for max_level > 2")

        number_of_leaves = len(network_list.represented_leaves())

        # Trinets
        trinet_list = network_list.networks_of_size(3)
        number_of_trinets = max(len(trinet_list), 1)
        number_of_possible_trinets = ncr(number_of_leaves, 3)
        trinet_level_counts = {level: 0 for level in range(max_level + 1)}
        for _, sub_trinet_list in trinet_list.per_leaf_set_iterator():
            total = sub_trinet_list.volume()
            for network_info in sub_trinet_list:
                try:
                    trinet_level_counts[network_info['strict_level']] += network_info.multiplicity / total
                except KeyError:  # when level is higher than max_level
                    pass

        trinet_level_percentages = {level: count / number_of_trinets for level, count in trinet_level_counts.items()}
        trinet_level_2_threshold = 0.5 * min(ncr(number_of_leaves - 1, 2), ncr(number_of_leaves - 2, 1)) / number_of_possible_trinets
        trinet_level_1_threshold = 0.5 * ncr(number_of_leaves - 1, 2) / number_of_possible_trinets
        trinet_level_thresholds = [0, trinet_level_1_threshold, trinet_level_2_threshold]
        trinet_threshold_passed_index = np.array(
            [int(percentage >= threshold) for percentage, threshold in zip(trinet_level_percentages.values(), trinet_level_thresholds)])
        trinet_best_level = max(np.where(trinet_threshold_passed_index == 1)[0])

        # Binets
        binet_list = network_list.networks_of_size(2)
        number_of_binets = max(len(binet_list), 1)
        number_of_possible_binets = ncr(number_of_binets, 3)
        binet_level_counts = {level: 0 for level in range(max_level + 1)}
        for _, sub_binet_list in binet_list.per_leaf_set_iterator():
            total = sub_binet_list.volume()
            for network_info in sub_binet_list:
                try:
                    binet_level_counts[network_info['strict_level']] += network_info.multiplicity / total
                except KeyError:  # when level is higher than max_level
                    pass

        binet_level_percentages = {level: count / number_of_binets for level, count in binet_level_counts.items()}
        binet_level_2_threshold = 1 / number_of_possible_binets
        binet_level_1_threshold = 0.5 * ncr(number_of_leaves - 1, 1) / number_of_possible_binets
        binet_level_thresholds = [0, binet_level_1_threshold, binet_level_2_threshold]
        binet_threshold_passed_index = np.array(
            [int(percentage >= threshold) for percentage, threshold in zip(binet_level_percentages.values(), binet_level_thresholds)])
        binet_best_level = max(np.where(binet_threshold_passed_index == 1)[0])

        trinet_binet_treshold_passed_index = (trinet_threshold_passed_index + binet_threshold_passed_index) >= 1
        best_level = max(np.where(trinet_binet_treshold_passed_index == 1)[0])

        # # TODO better naming and better boundaries
        # min_percentage_matrix = np.zeros((3, 3))
        # max_percentage_matrix = np.zeros((3, 3))
        # # Level 0
        # min_percentage_matrix[0][0] = 1.
        # max_percentage_matrix[0][0] = 1.
        #
        # # Level 1
        # temp = ncr(number_of_leaves - 1, 2) / number_of_possible_networks
        # min_percentage_matrix[1][0] = 1 - temp
        # max_percentage_matrix[1][0] = 1 - temp
        # min_percentage_matrix[1][1] = temp
        # max_percentage_matrix[1][1] = temp
        #
        # # Level 2
        # min_percentage_matrix[2][0] = min(
        #     (ncr(int(math.ceil((number_of_leaves - 1) / 2.)), 3) + ncr(int(math.floor((number_of_leaves - 1) / 2.)), 3)),
        #     ncr(number_of_leaves - 2, 3)
        # ) / number_of_possible_networks
        # max_percentage_matrix[2][0] = ncr(number_of_leaves - 1, 3) / number_of_possible_networks
        # min_percentage_matrix[2][2] = ncr(number_of_leaves - 2, 1) / number_of_possible_networks
        # max_percentage_matrix[2][2] = ncr(number_of_leaves - 1, 2) / number_of_possible_networks
        # min_percentage_matrix[2][1] = 1 - max_percentage_matrix[2][0] - max_percentage_matrix[2][2]
        # max_percentage_matrix[2][1] = 1 - min_percentage_matrix[2][0] - min_percentage_matrix[2][2]
        #
        # shortage_matrix = np.zeros((3, 3))
        # surplus_matrix = np.zeros((3, 3))
        # for network_level in range(max_level + 1):
        #     for trinet_level in range(3):
        #         shortage_matrix[network_level][trinet_level] = max(0, min_percentage_matrix[network_level][trinet_level] - level_percentage[trinet_level])
        #         surplus_matrix[network_level][trinet_level] = max(0, level_percentage[trinet_level] - max_percentage_matrix[network_level][trinet_level])
        #
        # alignment_matrix = shortage_matrix + surplus_matrix
        # score_array = list(np.sum(alignment_matrix, axis=1))
        # best_level = score_array.index(min(score_array))
        # best_score = int(100 * (1 - score_array[best_level]))
        # score_dict = {level: int(100 * (1 - score)) for level, score in enumerate(score_array)}
        # time.sleep(1)
        # print(min_percentage_matrix)
        # print(max_percentage_matrix)
        # print(level_percentage)
        # print(alignment_matrix)
        # print(score_array)
        # time.sleep(1)

        self.logger.info(f"Best level is {best_level}")
        # self.logger.info(f"Inconsistency: score of chosen level is {level_percentages}%")
        # self.logger.info(f"             : score of other levels is \n {pp.pformat(score_dict)}")
        return best_level

    def compute_generator_of_network_list(self, network_list: NetworkInfoList, level: int):
        self.logger.info("Computing generator ...")
        assert level >= 1, "This code only works for level >= 1"
        # Get trinets for finding out generator
        network_list = network_list.networks_where('strict_level', level)
        network_list = network_list.networks_with_category('generator_name')
        # TODO: ERROR --> IT CAN HAPPEN THAT THERE ARE NO TRINETS LEFT DUE TO MSS INCONSISTENCY
        #  e.g. : data/error_networks/mss_inconsistency.pickle
        #  --> Solve problem for sub network list ?

        # Find Generator
        # TODO: in case of draws/near draws, use same heuristics as in best_level
        all_generators = []
        generator_count = []
        generator_names = []
        for _, sub_network_list in network_list.per_leaf_set_iterator():
            total = sub_network_list.volume()
            for network_info in sub_network_list:
                try:
                    generator_index = generator_names.index(network_info['generator'].name)
                    generator_count[generator_index] += network_info.multiplicity / total
                except ValueError:
                    all_generators.append(network_info['generator'])
                    generator_names.append(network_info['generator'].name)
                    generator_count.append(network_info.multiplicity / total)

        max_count = max(generator_count)
        number_of_generators = sum(generator_count)
        max_index = generator_count.index(max_count)
        generator = all_generators[max_index]
        score = int(100 * max_count / number_of_generators)
        score_array = {generator.name: int(100 * count / number_of_generators) for generator, count in zip(all_generators, generator_count)}
        self.logger.info(f"Generator of underlying generator of minimal sink-set is {generator.name}")
        self.logger.info(f"Inconsistency: score of generator is {score}%")
        self.logger.info(f"             : scores of other generators are {score_array}")
        return copy.deepcopy(generator)

    @staticmethod
    def compute_leaf_location_info_of_network_list(network_list: NetworkInfoList, leaf_set: bidict, generator: RootedLevelKGenerator):
        # Get trinets with information of leaf locations
        network_list = network_list.networks_with_category('generator_name').networks_where('generator_name', generator.name)
        network_list.uniquify(equal_naming=True)

        # TODO: assert all networks have the same generator (? or use all trinets etc)

        # Leaf --> Edge --> (how many times leaf is on this edge side)
        number_of_leaves = len(leaf_set)
        edge_sides = list(set(generator.internal_edges))
        number_of_edge_sides = len(edge_sides)
        leaf_on_edge_side_count_dict = {leaf: {side: 0 for side in edge_sides} for leaf in leaf_set.values()}
        for _, sub_network_list in network_list.per_leaf_set_iterator():
            total = sub_network_list.volume()
            for network_info in sub_network_list:
                relation_dict = network_info['relation_dict']
                extra_leaf_dict = network_info['extra_leaf_dict']
                for leaf, edge_side in extra_leaf_dict.items():
                    translated_leaf = relation_dict[leaf]
                    leaf_on_edge_side_count_dict[translated_leaf][edge_side] += network_info.multiplicity / total
        # Put in matrix
        leaf_on_edge_side_count_matrix = np.zeros((number_of_edge_sides, number_of_leaves))
        for leaf_index, leaf in enumerate(leaf_set.values()):
            for edge_side_index, edge_side in enumerate(edge_sides):
                leaf_on_edge_side_count_matrix[edge_side_index][leaf_index] = leaf_on_edge_side_count_dict[leaf][edge_side]

        # Generator reticulation --> leaf --> (how many times leaf is this reticulation leaf}
        reticulation_sides = list(network_list[0]['reticulations'])
        number_of_reticulation_sides = len(reticulation_sides)
        leaf_is_reticulation_side_dict = {leaf: {ret_side: 0 for ret_side in reticulation_sides} for leaf in leaf_set.values()}
        for _, sub_network_list in network_list.per_leaf_set_iterator():
            total = sub_network_list.volume()
            for network_info in sub_network_list:
                relation_dict = network_info['relation_dict']
                for generator_leaf, leaf in relation_dict.items():
                    if generator_leaf in reticulation_sides:
                        leaf_is_reticulation_side_dict[leaf][generator_leaf] += network_info.multiplicity / total
        # Put in matrix
        leaf_on_reticulation_side_count_matrix = np.zeros((number_of_reticulation_sides, number_of_leaves))
        for leaf_index, leaf in enumerate(leaf_set.values()):
            for reticulation_side_index, reticulation_side in enumerate(reticulation_sides):
                leaf_on_reticulation_side_count_matrix[reticulation_side_index][leaf_index] = leaf_is_reticulation_side_dict[leaf][reticulation_side]
        return leaf_on_edge_side_count_matrix, leaf_on_reticulation_side_count_matrix, edge_sides, reticulation_sides

    def locate_leaves_using_ILP(self, leaf_set, leaf_on_edge_side_count_matrix, leaf_on_reticulation_side_count_matrix, edge_sides, reticulation_sides):
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
        self.logger.info(f"Placing leaves on side edges as follows: \n {pp.pformat(edge_side_leaf_dict)}")
        self.logger.info(f"Placing leaves on reticulation sides as follows: \n {pp.pformat(leaf_reticulation_side_bidict)}")
        self.logger.info(f"Inconsistency: ILP score is {score} %")
        self.logger.info(f"             : ILP score per leaf is \n {score_dict} %")
        return edge_side_leaf_dict, leaf_reticulation_side_bidict

    # TODO: naming
    def locate_leaves_using_normal(self, leaf_set, leaf_on_edge_side_count_matrix, leaf_on_reticulation_side_count_matrix, edge_sides, reticulation_sides):
        self.logger.info(f"Locating using using normal")
        score_dict = dict()
        number_of_reticulations = leaf_on_reticulation_side_count_matrix.shape[0]
        time.sleep(1)
        print(leaf_on_edge_side_count_matrix)
        print(leaf_set)
        print(edge_sides)
        time.sleep(1)

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

        score = sum(score_dict.values()) / len(leaf_set)
        self.logger.info(f"Placing leaves on side edges as follows: \n {pp.pformat(edge_side_leaf_dict)}")
        self.logger.info(f"Placing leaves on reticulation sides as follows: \n {pp.pformat(leaf_reticulation_side_bidict)}")
        self.logger.info(f"Inconsistency: score is {score} %")
        self.logger.info(f"             : score per leaf is \n {score_dict} %")
        return edge_side_leaf_dict, leaf_reticulation_side_bidict

    def compute_leaf_locations_of_network_list(self, network_list: NetworkInfoList, leaf_set: bidict, generator: RootedLevelKGenerator):
        self.logger.info("Computing leaf locations ...")

        assert network_list.network_size <= 3, "This method only works for trinet and binet sets"
        # TODO: For generators with symmetric reticulations leaves score will not be 100%

        # Get information on leaf locations
        leaf_on_edge_count_matrix, leaf_on_reticulation_side_count_matrix, edge_sides, reticulation_sides = self.compute_leaf_location_info_of_network_list(
            network_list,
            leaf_set,
            generator)

        # Get leaf placement
        if self.leaf_locator == 'ILP':
            edge_side_leaf_dict, leaf_reticulation_side_bidict \
                = self.locate_leaves_using_ILP(leaf_set, leaf_on_edge_count_matrix, leaf_on_reticulation_side_count_matrix, edge_sides, reticulation_sides)
        else:
            edge_side_leaf_dict, leaf_reticulation_side_bidict \
                = self.locate_leaves_using_normal(leaf_set, leaf_on_edge_count_matrix, leaf_on_reticulation_side_count_matrix, edge_sides, reticulation_sides)
        reticulation_side_order = [leaf_reticulation_side_bidict.inverse[reticulation] for reticulation in reticulation_sides]

        # Remove symmetry from generator 2c
        if generator.name == '2c':
            # Nodes (1,2), (3,4) symmetrical. Need to split edge (1,3) into two sets, one for each reticulation
            edge_1_3_leaves = edge_side_leaf_dict[('1', '3')]
            edge_1_3_leaf_indicis = sorted([leaf_set.inverse[leaf_name] for leaf_name in edge_1_3_leaves])
            leaf_order_matrix_2 = self.leaf_order_info_of_network_set(network_list, leaf_set)
            time.sleep(1)
            print(leaf_order_matrix_2)
            print(leaf_set)
            time.sleep(1)
            x = leaf_set.inverse[reticulation_side_order[0]]
            y = leaf_set.inverse[reticulation_side_order[1]]
            leaf_per_reticulation = (np.array([np.argmax(row) for row in leaf_order_matrix_2[edge_1_3_leaf_indicis][:, [x, y]]]))
            # TODO: currently happens that values can be equal for each reticulation and a random one is chosen :(
            edge_side_leaf_dict[('1', '3')] = set([leaf_set[edge_1_3_leaf_indicis[leaf_index]] for leaf_index in np.where(leaf_per_reticulation == 0.)[0]])
            edge_side_leaf_dict[('1', '4')] = set([leaf_set[edge_1_3_leaf_indicis[leaf_index]] for leaf_index in np.where(leaf_per_reticulation == 1.)[0]])

        # return edge_leaf_dict, chosen_reticulation_relation_dict, score / number_of_leaves
        self.logger.info(f"Found that the leaves are on the edges as follows \n {pp.pformat(edge_side_leaf_dict)}")
        self.logger.info(f"Found that the reticulations in the generator should be renamed as follows \n {pp.pformat(leaf_reticulation_side_bidict)}")
        # self.logger.info(f"Inconsistency: total score of leaf locations is {score}%")
        # self.logger.info(f"             : score per leaf is \n {pp.pformat(score_dict)}, ")
        return edge_side_leaf_dict, leaf_reticulation_side_bidict, reticulation_sides

    @staticmethod
    def leaf_order_info_of_network_set(network_list: NetworkInfoList, leaf_set: bidict):
        leaf_order_matrix = np.zeros((len(leaf_set), len(leaf_set)))
        for _, sub_network_list in network_list.per_leaf_set_iterator():
            total = sub_network_list.volume()
            for network_info in sub_network_list:
                leaf_order = network_info['leaf_order']
                for leaf_ord in leaf_order:
                    for low_leaf_index in range(len(leaf_ord)):
                        for high_leaf_index in range(low_leaf_index + 1, len(leaf_ord)):
                            low_leaf = leaf_set.inverse[leaf_ord[low_leaf_index]]
                            high_leaf = leaf_set.inverse[leaf_ord[high_leaf_index]]
                            leaf_order_matrix[low_leaf, high_leaf] += network_info.multiplicity / total
        return leaf_order_matrix

    def compute_network_of_mss(self, minimal_sink_set) -> list:
        self.logger.info(f"Computing network of minimal sink-set {minimal_sink_set}")
        trinet_info_list_mss = self.trinet_info_list.find_networks_with_leaf_names_in(minimal_sink_set)

        trinet_info_list_mss.add_info(self.standard_trinet_info_list)
        self.logger.info(f"Number of networks that describe this mss is {len(trinet_info_list_mss)}")


        leaf_set = bidict({index: leaf for index, leaf in enumerate(minimal_sink_set)})

        # Level
        level = self.compute_level_of_network_list(trinet_info_list_mss)
        if level == 0:
            if len(minimal_sink_set) > 2:
                if len(minimal_sink_set) == len(self.taxa):
                    #TODO: Check if this can happen (in case of ties between trinets describing he same leaf sets) and do something
                    for _, sub_network_info_list in self.trinet_info_list.per_leaf_set_iterator():
                        for network_info in sub_network_info_list:
                            network_info.network.visualize()
                sub_solver = Solver(self.standard_trinet_info_list, trinet_info_list_mss, self.leaf_locator, is_main=False)
                sub_solver.solve()
                return sub_solver.steps
            else:
                best_multiplicity = 0
                best_network = None
                for network_info in trinet_info_list_mss:
                    if network_info.multiplicity > best_multiplicity:
                        best_multiplicity = network_info.multiplicity
                        best_network = copy.deepcopy(network_info.network)
            steps = [{'network': best_network, 'leaf_set': set(minimal_sink_set), 'type': 'shrink', 'name': mss_leaf_name(minimal_sink_set)}]
            return steps

        # Generator
        generator = self.compute_generator_of_network_list(trinet_info_list_mss, level)
        sides_per_edge = sides_per_edge_per_generator[generator.name]

        # Leaves per side
        edge_side_leaf_dict, leaf_reticulation_side_bidict, reticulation_sides = self.compute_leaf_locations_of_network_list(trinet_info_list_mss, leaf_set,
                                                                                                                             generator)
        # Order leaves per side
        ordered_edge_leaf_dict = self.compute_leaf_order(trinet_info_list_mss, leaf_set, generator, edge_side_leaf_dict, sides_per_edge, reticulation_sides)

        # Rename reticulations
        for reticulation in reticulation_sides:
            leaf = leaf_reticulation_side_bidict.inverse[reticulation]
            generator.rename_node(new_name=leaf, old_name=reticulation)

        # Add leaves to edges
        self.add_leaves_to_edges(generator, ordered_edge_leaf_dict)
        step = [{'network': generator, 'leaf_set': set(minimal_sink_set), 'type': 'shrink', 'name': mss_leaf_name(minimal_sink_set)}]
        return step

    def compute_leaf_order(self, network_list: NetworkInfoList, leaf_set: bidict, generator: RootedLevelKGenerator, edge_side_leaf_dict: dict,
                           sides_per_edge_side: dict,
                           reticulations: list):
        self.logger.info(f"Computing leaf order for leaves {list(leaf_set.values())}")
        # Get trinets of level-1 for ordering of leaf on sides
        network_list = network_list.networks_with_category('relation_dict')
        network_list = network_list.networks_where('strict_level', 1) + network_list.networks_where('strict_level', 2)
        network_list.uniquify(equal_naming=True, count=True)
        leaf_order_matrix_1 = self.leaf_order_info_of_network_set(network_list, leaf_set)

        ordered_edge_leaf_dict = dict()
        # Order leaves per side
        if generator.name in ('1', '2a', '2b', '2d'):
            for edge, leaves in edge_side_leaf_dict.items():
                if len(leaves) == 0:
                    ordered_edge_leaf_dict[edge] = []
                    continue
                sides = len(sides_per_edge_side[edge])
                sub_leaf_order_matrix, sub_leaf_set = self.sub_leaf_order_matrix(leaf_set, leaf_order_matrix_1, leaves)
                ordered_leaves = self.order_leaves_on_sides(sub_leaf_set, sub_leaf_order_matrix, sides=sides)
                ordered_edge_leaf_dict[edge] = ordered_leaves
        if generator.name == '2c':
            ret_0_parent = generator.nodes_above_nodes({reticulations[0]}, max_height=1)
            ret_0_parent.remove(reticulations[0])
            ret_0_parent = ret_0_parent.pop()
            ret_1_parent = generator.nodes_above_nodes({reticulations[1]}, max_height=1)
            ret_1_parent.remove(reticulations[1])
            ret_1_parent = ret_1_parent.pop()
            edges = [[('0', '1'), ('1', ret_0_parent), ('1', ret_1_parent)], [('0', '2'), ('2', ret_0_parent), ('2', ret_1_parent)]]
            leaves_on_edges = [[[], [], []], [[], [], []]]

            side_orders = []
            for index, edge in enumerate(edges[0]):
                if edge not in edge_side_leaf_dict.keys():
                    continue
                leaves = edge_side_leaf_dict[edge]
                if len(leaves) == 0:
                    continue
                sides = len(sides_per_edge_side[edge])
                sub_leaf_order_matrix, sub_leaf_set = self.sub_leaf_order_matrix(leaf_set, leaf_order_matrix_1, leaves)
                ordered_leaves = self.order_leaves_on_sides(sub_leaf_set, sub_leaf_order_matrix, sides=sides)
                side_orders.append(ordered_leaves)
                for row_index, row in enumerate(ordered_leaves):
                    leaves_on_edges[row_index][index] = row
            side_separation = self.align_sides(side_orders, leaf_set, leaf_order_matrix_1)

            leaves_on_edges_ordered = copy.deepcopy(leaves_on_edges)
            for leg_index, leg in enumerate(leaves_on_edges[0]):
                if set(leg).issubset(side_separation[1]):
                    leaves_on_edges_ordered[0][leg_index], leaves_on_edges_ordered[1][leg_index] = leaves_on_edges_ordered[1][leg_index], \
                                                                                                   leaves_on_edges_ordered[0][leg_index]
            ret_0_left_alignment = len(set(edge_side_leaf_dict[('1', '3')]).intersection(leaves_on_edges_ordered[0][1] + leaves_on_edges_ordered[1][1]))
            ret_0_right_alignment = len(set(edge_side_leaf_dict[('1', '3')]).intersection(leaves_on_edges_ordered[0][2] + leaves_on_edges_ordered[1][2]))
            ret_1_left_alignment = len(set(edge_side_leaf_dict[('1', '4')]).intersection(leaves_on_edges_ordered[0][1] + leaves_on_edges_ordered[1][1]))
            ret_1_right_alignment = len(set(edge_side_leaf_dict[('1', '4')]).intersection(leaves_on_edges_ordered[0][2] + leaves_on_edges_ordered[1][2]))
            if ret_0_right_alignment + ret_1_left_alignment > ret_0_left_alignment + ret_1_right_alignment:
                leaves_on_edges_ordered[0][1], leaves_on_edges_ordered[0][2] = leaves_on_edges_ordered[0][2], leaves_on_edges_ordered[0][1]
                leaves_on_edges_ordered[1][1], leaves_on_edges_ordered[1][2] = leaves_on_edges_ordered[1][2], leaves_on_edges_ordered[1][1]

            for i in range(len(leaves_on_edges_ordered)):
                for j in range(len(leaves_on_edges_ordered[i])):
                    ordered_edge_leaf_dict[edges[i][j]] = [leaves_on_edges_ordered[i][j]]

        self.logger.info(f"Found that leaves should be ordered as \n {ordered_edge_leaf_dict}")
        return ordered_edge_leaf_dict

    @staticmethod
    def add_leaves_to_edges(generator, edge_side_leaf_dict):
        for edge, sides in edge_side_leaf_dict.items():
            for side in sides:
                from_node = edge[0]
                to_node = edge[1]
                for leaf in side:
                    internal_name, _ = generator.add_leaf_to_edge(edge=[from_node, to_node], leaf_name=leaf)
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

    def order_leaves_on_sides(self, leaf_set, leaf_order_matrix, sides=1):
        assert sides in (1, 2), "Only works for 1 or 2 sides"
        self.logger.info(f"Ordering {list(leaf_set.values())} into {sides} sides")
        side_orders = ([[] for _ in range(sides)])
        placed_leaves = []
        leaf_order_scores = {side: dict() for side in range(sides)}

        # place all leaves with similar lineage on the same side
        for _ in [0, 1]:
            # 0 --> check lineage below
            # 1 --> check lineage above for leftover leaves (transposed leaf_order_matrix)
            while True:
                row_sum = leaf_order_matrix.sum(axis=1)
                row_sum = [s if index not in placed_leaves else 0 for index, s in enumerate(row_sum)]
                if sum(row_sum) == 0:
                    break
                else:
                    # leaf with maximum out-degree
                    leaf = np.argmax(row_sum)
                    best_side, score = self.align_leaf_on_sides(leaf, side_orders, leaf_order_matrix)
                    side_orders[best_side].append(leaf)
                    leaf_order_scores[best_side][leaf] = int(100 * score)
                    placed_leaves.append(leaf)
            leaf_order_matrix = leaf_order_matrix.T

        # Place leaf without any lineage on an empty side
        if len(placed_leaves) < len(leaf_set):
            for leaf_index in leaf_set.keys():
                if leaf_index in placed_leaves:
                    continue
                best_side_index = None
                best_side_order = None
                number_of_leaves_on_side = max([len(side_order) for side_order in side_orders])
                for side_index, side_order in enumerate(side_orders):
                    if len(side_order) <= number_of_leaves_on_side:
                        number_of_leaves_on_side = len(side_order)
                        best_side_index = side_index
                        best_side_order = side_order
                best_side_order.append(leaf_index)
                leaf_order_scores[best_side_index][leaf_index] = 0 # TODO

        leaf_order_scores = {side: {leaf_set[leaf_index]: score for leaf_index, score in leaf_score.items()} for side, leaf_score in leaf_order_scores.items()}
        leaf_order = [[leaf_set[leaf_index] for leaf_index in side_order] for side_order in side_orders]
        self.logger.info(f"Found that leaves should be ordered as \n {pp.pformat(leaf_order)}")
        self.logger.info(f"Inconsistency: score of leaf ordering is \n {pp.pformat(leaf_order_scores)}")
        return leaf_order

    @staticmethod
    def align_leaf_on_sides(leaf, side_orders, leaf_order_matrix):
        number_of_sides = len(side_orders)
        current_leaf_children = set(list(np.where(leaf_order_matrix[leaf] >= 1)[0]) + [leaf])
        L = 0.5 * len(current_leaf_children)

        # Count alignment of lineage per side
        side_alignment = np.zeros(number_of_sides)
        for side, side_order in enumerate(side_orders):
            for other_leaf in side_order:
                leaf_children = set(list(np.where(leaf_order_matrix[other_leaf] >= 1)[0]) + [other_leaf])
                same_children = leaf_children.intersection(current_leaf_children)
                if len(same_children) >= L:
                    side_alignment[side] += 1

        # Normalize alignment
        normalized_side_alignment = np.zeros(number_of_sides)
        divide_boolean = True
        for side in range(number_of_sides):
            other_side = (side + 1) % number_of_sides
            normalized_side_alignment[side] = side_alignment[side] + len(side_orders[other_side]) - side_alignment[other_side]
            if len(side_orders[side]) == 0:
                divide_boolean = False
        if divide_boolean:
            for side in range(number_of_sides):
                normalized_side_alignment /= len(side_orders[side])

        # Pick best side
        best_side = int(np.argmax(normalized_side_alignment))
        score = normalized_side_alignment[best_side] / sum(normalized_side_alignment) if sum(normalized_side_alignment) != 0 \
            else normalized_side_alignment[best_side]
        return best_side, score

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
            if side_order_matrix.ndim > 1:
                new_col = np.sum((side_order_matrix[:, leaf_indicis_of_side]), axis=1)
            else:
                new_col = np.array([[np.sum((side_order_matrix[leaf_indicis_of_side]))]])
            try:
                side_order_matrix_2 = np.vstack([side_order_matrix_2, new_col])
            except ValueError:
                side_order_matrix_2 = new_col

        for i in range(index):
            side_order_matrix_2[i, i] = 0

        order = self.order_leaves_on_sides(side_set, side_order_matrix_2.T, sides=2)

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

    def __getstate__(self):
        self.logger_name = 'solver.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger_name)
        return self.__dict__
