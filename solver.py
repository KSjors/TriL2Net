from datastructures.omega import *
from datastructures.rooted_level_k_network import *
from data.all_trinets import *
import copy
import time

# all_generators, all_trinets, all_trinets_gen_sides = get_trinets()


class Solver:
    def __init__(self, trinet_set):
        self.trinet_sets = [trinet_set]
        self.transformations = {}

    def next_transform(self):
        if len(self.trinet_sets[-1].taxa_names) == 2:
            print("WARNING: No more transformations to do")
            return None
        current_omega = Omega(self.trinet_sets[-1])
        current_mss = current_omega.minimal_sink_sets()
        current_mss.sort(key=len)
        for mss in current_mss:
            if len(mss) >= 2:
                return self.transform(mss)

    def transform(self, mss):
        mss_trinets = dict()
        if len(mss) == 2:
            # Find trinet of which mss is part
            current_taxa = list(self.trinet_sets[-1].taxa_names)
            third_leaf = current_taxa.pop(-1)
            while third_leaf in mss:
                third_leaf = current_taxa.pop(-1)
            trinet = sorted(mss + [third_leaf])
            trinet_of_mss = self.trinet_sets[-1].trinets[str(trinet)]
            mss_network = sub_network_of(trinet_of_mss, mss)
            mss_network.visualize()
        elif len(mss) == 3:
            trinet = sorted(mss)
            mss_network = self.trinet_sets[-1].trinets[str(trinet)]
            mss_network.visualize()
        else:
            trinet_iterator = itertools.combinations(mss, 3)
            for trinet in trinet_iterator:
                trinet = sorted(trinet)
                mss_trinets[str(trinet)] = self.trinet_sets[-1].trinets[str(trinet)]
            get_network(mss, mss_trinets)

        new_trinet_set = copy.deepcopy(self.trinet_sets[-1])
        new_name = new_trinet_set.suppress_minimal_sink_set(mss)
        self.transformations[new_name] = mss
        self.trinet_sets.append(new_trinet_set)

    def state(self, number=-1):
        if not -1 <= number < len(self.trinet_sets):
            print("number not in range")
        current_omega = Omega(self.trinet_sets[number])
        current_omega.visualize()

    def __str__(self):
        return str(self.trinet_sets[-1]) + " " + str(self.transformations[-1])


def get_network(taxa, trinets):
    generator_count = np.zeros(len(all_generators))
    for trinet_taxa, trinet in trinets.items():
        generator_number = all_trinets.index(trinet)


def get_network_structure_old(taxa, trinets):
    # Assumes network is biconnected

    """ Make this function by generating al possible generators and check which one it is
    """
    max_ret = 0
    max_trinets = []
    # Checks amount of reticulations in trinet
    for trinet_taxa, trinet in trinets.items():
        taxa_list = [x[1:-1] for x in trinet_taxa[1: -1].split(', ')]
        _, _, reticulations = trinet.number_of_internals_leaves_reticulations()
        if reticulations > max_ret:
            max_ret = reticulations
            max_trinets = [[taxa_list, trinet]]
        elif reticulations == max_ret:
            max_trinets.append([taxa_list, trinet])

    # Based on number of reticulations do:
    if max_ret == 0:
        print("WARNING: taxa_set of size >=4 has no reticulations")
        print(taxa)
        print(trinets)
        return None, None
    elif max_ret == 1:
        net = max_trinets[-1][1]
        reticulations = net.get_reticulations()
        reticulation_leaves = net.leaf_children(set(reticulations))[1]
        return '1A', reticulation_leaves
    elif max_ret > 2:
        print("WARNING: this trinet is level 3")
        print(taxa)
        print(trinets)
        return None, None

    # Check for redundant nodes in max_trinets in order to differentiate between different generators
    # Finds set of minimum number of nodes such that subnet of these nodes still has same amount of reticulations
    min_max_nets = {}
    min_len = 3
    for taxa_list, trinet in max_trinets:
        subset_iterator = all_combinations(taxa_list, 1, 2)
        for subset in subset_iterator:
            net = sub_network_of(trinet, subset, suppress_redundant='none')
            _, _, reticulations = net.number_of_internals_leaves_reticulations()
            if reticulations == max_ret:
                necessary_leaves = subset
                if len(necessary_leaves) < min_len:
                    min_max_nets[str(necessary_leaves)] = sub_network_of(trinet, necessary_leaves, suppress_redundant='none')
                    min_len = len(necessary_leaves)
                elif len(necessary_leaves) == min_len:
                    min_max_nets[str(necessary_leaves)] = sub_network_of(trinet, necessary_leaves, suppress_redundant='none')

    # Based on number of necessary nodes do:
    if min_len == 1:
        net = next(iter(min_max_nets.values()))
        reticulations = net.get_reticulations()
        reticulation_leaves = net.get_leaf_children(set(reticulations))[1]
        return '2A', reticulation_leaves
    elif min_len == 2:
        # If there two or more necessary sets, generator must be of this type
        if len(min_max_nets) > 1:
            net = next(iter(min_max_nets.values()))
            reticulations = net.get_reticulations()
            reticulation_leaves = net.get_leaf_children(set(reticulations))[1]
            return '2D', reticulation_leaves
        else:
            taxa_str = list(min_max_nets.keys())[0]
            binet = min_max_nets[taxa_str]
            taxa_list = [x[1:-1] for x in taxa_str[1: -1].split(', ')]
            leaf = taxa_list[0]
            uninet = sub_network_of(binet, [leaf], suppress_redundant='none', suppress_parallel=False)
            internals, leaves, reticulations = uninet.number_of_internals_leaves_reticulations()
            # Generator D will prune to 2 or to 0 reticulations, B and C to 1
            if reticulations in (0, 2):
                net = next(iter(min_max_nets.values()))
                reticulations = net.get_reticulations()
                reticulation_leaves = net.get_leaf_children(set(reticulations))[1]
                return '2D', reticulation_leaves
            else:
                root = binet.get_root_name()
                children_root = binet.get_children({root}, 2)
                # Generator B has leaf closer to root than C does (2 vs 3 distance respectively)
                if binet.contains_leaf(children_root):
                    # Generator B has asymmetry in reticulations, first one should be the closest to root
                    children_root_1 = binet.get_children({root}, 1)
                    net = next(iter(min_max_nets.values()))
                    reticulations = net.get_reticulations()
                    if reticulations[0] not in children_root_1:
                        reticulations.reverse()
                    reticulation_leaf_0 = net.get_leaf_children(set(reticulations))[1].pop()
                    reticulation_leaf_1 = net.get_leaf_children(set(reticulations))[1].pop()
                    return '2D', [reticulation_leaf_0, reticulation_leaf_1]
                else:
                    net = next(iter(min_max_nets.values()))
                    reticulations = net.get_reticulations()
                    reticulation_leaves = net.get_leaf_children(set(reticulations))[1]
                    net.visualize()
                    return '2C', reticulation_leaves
    else:
        print("WARNING: at least three leaves needed for two reticulations")
        return None, None


