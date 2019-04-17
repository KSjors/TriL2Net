import logging
from datastructures.omega import *
from datastructures.rooted_level_k_network import *
from data.all_trinets import *
from datastructures.trinet_set import TrinetSet
from utils.help_functions import guid
import copy
import time

# regenerate_trinets()
all_generators, all_trinets, all_trinets_gen_sides = get_trinets()


class Solver:
    def __init__(self, trinet_set: TrinetSet):
        self.uid = guid()
        self.logger = logging.getLogger("solver.{}".format(self.uid))
        self.logger.debug("Creating new solver object from trinet set {}.".format(trinet_set.uid))
        self.trinet_sets = [trinet_set]
        self.transformations = {}

    def next_transform(self):
        self.logger.debug("Performing the next transformation.")
        if len(self.trinet_sets[-1].taxa_names) == 1:
            print("WARNING: No more transformations to do")
            return None
        mss = self.get_next_mss()
        return self.transform(mss)

    def get_next_mss(self):
        current_omega = Omega(self.trinet_sets[-1])
        current_mss = current_omega.minimal_sink_sets()
        current_mss.sort(key=len)
        for mss in current_mss:
            if len(mss) >= 2:
                return mss

    def transform(self, mss: list):
        self.logger.debug("Transforming {}".format(mss))
        mss_trinets = dict()
        if len(mss) == 2:
            # TODO: discard all faulty trinets first
            # Find trinet of which mss is part
            current_taxa = list(self.trinet_sets[-1].taxa_names)
            if len(current_taxa) == 2:
                return "Cherry"
            third_leaf = current_taxa.pop(-1)
            while third_leaf in mss:
                third_leaf = current_taxa.pop(-1)
            trinet = sorted(mss + [third_leaf])
            trinet_of_mss = self.trinet_sets[-1].trinet_dict[str(trinet)]
            mss_network = RootedLevelKNetwork.from_network(trinet_of_mss, set(mss))
            # mss_network.visualize()
        elif len(mss) == 3:
            # TODO what if this is the only trinet that disagrees
            trinet = sorted(mss)
            mss_network = self.trinet_sets[-1].trinet_dict[str(trinet)]
            # mss_network.visualize()
        else:
            trinet_iterator = itertools.combinations(mss, 3)
            for trinet in trinet_iterator:
                trinet = sorted(trinet)
                mss_trinets[str(trinet)] = self.trinet_sets[-1].trinet_dict[str(trinet)]
            mss_network = get_network(mss_trinets)

        # TODO: suppress in @classmethod?
        # TODO: save network of mss
        new_trinet_set = TrinetSet.copy_trinet_set(self.trinet_sets[-1])
        new_name = new_trinet_set.suppress_minimal_sink_set(mss)
        self.transformations[new_name] = mss
        self.trinet_sets.append(new_trinet_set)
        return mss_network

    def state(self, number=-1):
        if not -1 <= number < len(self.trinet_sets):
            print("number not in range")
        current_omega = Omega(self.trinet_sets[number])
        current_omega.visualize()

    def __str__(self):
        return str(self.trinet_sets[-1]) + " " + str(self.transformations[-1])

    def __getstate__(self):
        self.logger = 'solver.{}'.format(self.uid)
        result = copy.deepcopy(self.__dict__)
        self.logger = logging.getLogger(self.logger)
        return result

    def __setstate__(self, d):
        self.__dict__ = d
        self.logger = logging.getLogger(self.logger)
        return self.__dict__


def get_network(trinets):
    generator_count = [[] for i in all_generators]
    generator_count.append([])
    c = 0
    for trinet_taxa, trinet in trinets.items():
        c += 1
        try:
            trinet_number = all_trinets.index(trinet)
            generator = all_trinets_gen_sides[trinet_number][1]
            sides = all_trinets_gen_sides[trinet_number][2]
            generator_number = all_generators.index(generator)
            generator_count[generator_number].append(trinet_taxa)
        except ValueError:
            generator_count[-1].append(trinet_taxa)
            try:
                assert not trinet.is_biconnected(leafless=True), "Can not find trinet {} with taxa {} in all trinets even though it is biconnected {}.".format(
                    trinet.uid, trinet.leaf_names, c)
            except:
                trinet.visualize()
                raise AssertionError
    return generator_count
