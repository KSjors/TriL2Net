import logging
from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkSet, NetworkInfo
from datastructures.solver import Solver
from utils.help_functions import *
import typing
import datetime as dt
import time
import json
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
from utils.help_functions import coalesce
import collections
import functools
import operator

"""
                Data:
Input network:
- number of leaves
- evolution probability (recombination, termination)
---> network summary (number of biconnected components, level of each bc, size of each component, generator of each component)
Trinet set noise:
- Noise types: tail move, uniform, deletion --> single / combination
- Noise strength
---> trinet set summary (redundancy, inconsistency, density)
Output network:
- network summary (number of biconnected components, level of each bc, size of each component, generator of each component)
- consistencies (trinet, triplet, cluster, cut-arc set)
Process:
- Logging
- Duration
- Scores

-------------------------------------------------------------
                Experiment:
- For varying networks and varying noise run algorithm
- Save:
    - input network
    - evolution settings
    - noise settings
    - input trinet set
    - output network
    - consistency scores
    - logging
    - duration
    - scores 

-------------------------------------------------------------
                Experiment:
- Take x networks on subsets of leaves of X. Take the union of their trinets as input
- Save:
    - input network
    - evolution settings
    - subset creation settings
    - input sub networks
    - input trinet set
    - output network
    - consistency scores
    - logging
    - duration
    - scores 
     

    
                Analyse:
- relation number of leaves and noise to consistency scores
    Hypothesis:
        1:  n -> n+1: 0.5(n^2 - n) as many trinets, goes to 0.5 n^2 as many for large n. 
            Hence same cut-arc set described by many more trinets and less prone to influence of noise 
        2:  There is a percentage of noise limit, beyond this percentage, inconsistency levels rise steeply
            this limit increases if number of leaves increases, but converges to a certain percentage
            tail move noise has higher limit than random
- Does number of biconnected components influence robustness? 
    Hypothesis
        1: More smaller components has better robustness
- Does level of biconnected components influence robustness?
    Hypothesis
        1:  Lower level probably has better robustness


"""


class NetworkGenerator:
    def __init__(self, n_range, r_range, t_range):
        self.n = n_range
        self.r = r_range
        self.t = t_range

    def iterator(self):
        param_iterator = itertools.product(self.n, self.r, self.t)
        for parameters in param_iterator:
            n, r, t = parameters
            network = RootedLevelKNetwork.random(n, r, level=2)
            if t > 0:
                network.terminate_percentage_leaves(t)
            yield parameters, network

    def __len__(self):
        return len(self.n) * len(self.r) * len(self.t)


class TrinetSetGenerator:
    def __init__(self, t, u, d, all_trinet_list):
        # save parameters
        self.t = t
        self.u = u
        self.d = d

        self.all_trinet_list = all_trinet_list

    def iterator(self, trinet_set):
        param_iterator = itertools.product(self.t, self.u, self.d)

        for parameters in param_iterator:
            t, u, d = parameters
            distorted_trinet_set = copy.deepcopy(trinet_set)
            distorted_trinet_set = NetworkSet.tail_move_distort(distorted_trinet_set, t)
            distorted_trinet_set = NetworkSet.uniform_distort(distorted_trinet_set, u, self.all_trinet_list)
            distorted_trinet_set = NetworkSet.deletion_distort(distorted_trinet_set, d)
            yield parameters, distorted_trinet_set

    def __len__(self):
        return len(self.t) * len(self.u) * len(self.d)


class SolverGenerator:
    def __init__(self
                 , all_trinet_set: NetworkSet
                 , cut_arc_set_count_method, minimal_sink_set_method, leaf_locator_method, level_threshold_method, level_count_method, generator_count_method
                 , symmetric_sides_set_count_method, leaf_order_count_method, leaf_order_method):
        # Set up logging

        # Save config
        self.cut_arc_set_count_method = cut_arc_set_count_method
        self.minimal_sink_set_method = minimal_sink_set_method
        self.leaf_locator_method = leaf_locator_method
        self.level_threshold_method = level_threshold_method
        self.level_count_method = level_count_method
        self.generator_count_method = generator_count_method
        self.symmetric_side_set_count_method = symmetric_sides_set_count_method
        self.leaf_order_count_method = leaf_order_count_method
        self.leaf_order_method = leaf_order_method

        # Save data
        self.all_trinet_set = all_trinet_set

    def iterator(self, problem_trinet_set: NetworkSet):
        param_iterator = itertools.product(self.cut_arc_set_count_method, self.minimal_sink_set_method, self.leaf_locator_method, self.level_threshold_method
                                           , self.level_count_method, self.generator_count_method, self.symmetric_side_set_count_method,
                                           self.leaf_order_count_method, self.leaf_order_method)
        for parameters in param_iterator:
            solver = Solver(self.all_trinet_set, problem_trinet_set, *parameters, is_main=True)
            parameters_dict = {
                'cut_arc_set_count_method'          : parameters[0]
                , 'minimal_sink_set_method'         : parameters[1]
                , 'leaf_locator_method'             : parameters[2]
                , 'level_threshold_method'          : parameters[3]
                , 'level_count_method'              : parameters[4]
                , 'generator_count_method'          : parameters[5]
                , 'symmetric_sides_set_count_method': parameters[6]
                , 'leaf_order_count_method'         : parameters[7]
                , 'leaf_order_method'               : parameters[8]
            }
            yield parameters_dict, solver

    def __len__(self):
        return len(self.cut_arc_set_count_method) * len(self.minimal_sink_set_method) \
               * len(self.leaf_locator_method) * len(self.level_threshold_method) \
               * len(self.level_count_method) * len(self.generator_count_method) \
               * len(self.symmetric_side_set_count_method) * len(self.leaf_order_count_method) \
               * len(self.leaf_order_method)


class ExperimentGenerator:
    def __init__(self
                 , ng: NetworkGenerator
                 , tsg: TrinetSetGenerator
                 , sg: SolverGenerator
                 , name: str
                 , trinet_set_method: int
                 , max_processes: int = 1
                 ):
        # Set up logging

        # Save generators
        self.ng = ng
        self.tsg = tsg
        self.sg = sg

        # Save config
        self.trinet_set_method = trinet_set_method
        self.max_processes = max_processes

        # Description and output
        self.name = name
        self.logger = None

    def iterator(self):
        t0 = time.time()
        for network_params, network in tqdm(self.ng.iterator(), total=len(self)):
            experiment = Experiment(self.name)
            network.visualize()
            try:
                t1 = time.time()
                experiment['input network creation'] = t1 - t0
                experiment['input network parameters'] = network_params
                experiment['input network'] = network.enewick()
                input_trinet_set, t2 = NetworkSet.induced_strict_network_set(network, 3, self.max_processes, False, self.trinet_set_method), time.time()
                experiment['input trinet set creation'] = t2 - t1
                experiment['input trinet set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                                  input_trinet_set.per_network_info()}
                input_triplet_set, t3 = NetworkSet.induced_strict_tree_set(network, 3, 1, False), time.time()  # TODO: from trinets_set?
                experiment['input triplet set creation'] = t3 - t2
                experiment['input triplet set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                                   input_triplet_set.per_network_info()}
                input_cluster_set, t4 = NetworkSet.induced_cluster_set(network, 1, False), time.time()
                experiment['input cluster set creation'] = t4 - t3
                experiment['input cluster set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                                   input_cluster_set.per_network_info()}
                for trinet_set_params, distorted_trinet_set in self.tsg.iterator(input_trinet_set):
                    t5 = time.time()
                    experiment['distorted trinet set creation'] = t5 - t4
                    experiment['distortion parameters'] = trinet_set_params
                    experiment['distorted trinet set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                                          distorted_trinet_set.per_network_info()}
                    for solver_params, solver in self.sg.iterator(distorted_trinet_set):
                        t6 = time.time()
                        experiment['solver creation'] = t6 - t5
                        experiment['solver parameters'] = solver_params
                        experiment['init'] = True
                        output_network, solver_scores = solver.solve()
                        t7 = time.time()
                        experiment['output network'] = output_network.enewick()
                        experiment['solving'] = t7 - t6
                        output_trinet_set, t8 = NetworkSet.induced_strict_network_set(output_network, 3, self.max_processes, False,
                                                                                      self.trinet_set_method), time.time()
                        experiment['output trinet set creation'] = t8 - t7
                        experiment['output trinet set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                                           output_trinet_set.per_network_info()}
                        output_triplet_set, t9 = NetworkSet.induced_strict_tree_set(output_network, 3, 1, False), time.time()  # TODO: from trinets_set?
                        experiment['output triplet set creation'] = t9 - t8
                        experiment['output triplet set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                                            output_triplet_set.per_network_info()}
                        output_cluster_set, t10 = NetworkSet.induced_cluster_set(output_network, 1, False), time.time()
                        experiment['output cluster set creation'] = t10 - t9
                        experiment['output cluster set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                                            output_cluster_set.per_network_info()}
                        IO_tn_cs, t11 = NetworkSet.consistency_score(input_trinet_set, output_trinet_set), time.time()
                        experiment['trinet'] = IO_tn_cs
                        experiment['trinet consistency score computation'] = t11 - t10
                        IO_tp_cs, t12 = NetworkSet.consistency_score(input_triplet_set, output_triplet_set), time.time()
                        experiment['triplet'] = IO_tp_cs
                        experiment['triplet consistency score computation'] = t12 - t11
                        IO_ct_cs, t13 = NetworkSet.consistency_score(input_cluster_set, output_cluster_set), time.time()
                        experiment['cluster'] = IO_ct_cs
                        experiment['cluster consistency score computation'] = t13 - t12
                        IO_cas_cs, t14 = network.cut_arc_set_consistency(output_network), time.time()
                        experiment['cut-arc set'] = IO_cas_cs
                        experiment['cut-arc set consistency score computation'] = t14 - t13
                        equal, t15 = network.equal_structure(output_network, equal_naming=True), time.time()
                        experiment['equal'] = equal[0]
                        experiment['equal computation'] = t15 - t14
                        experiment['full'] = True
                        yield experiment
                    t4 = time.time()
            except Exception as e:
                experiment.save()
                raise e
            t0 = time.time()

    def run(self):
        es = ExperimentSet(self.name)
        for experiment in self.iterator():
            experiment.save()
            es.append(experiment)
        return es

    def run_times(self, times):
        for _ in range(times):
            self.run()

    def __len__(self):
        return len(self.ng) * len(self.tsg) * len(self.sg)


class Experiment:
    def __init__(self, name, data=None):
        if data is None:
            data = {
                'finished'            : {'init': False, 'full': False, 'location': None}
                , 'parameters'        : {
                    'input network parameters': None
                    , 'distortion parameters' : None
                    , 'solver parameters'     : None
                }
                , 'input output'      : {
                    'input network'         : None
                    , 'input trinet set'    : None
                    , 'input triplet set'   : None
                    , 'input cluster set'   : None
                    , 'distorted trinet set': None
                    , 'output network'      : None
                    , 'solver scores'       : None
                    , 'output trinet set'   : None
                    , 'output triplet set'  : None
                    , 'output cluster set'  : None
                }
                , 'consistency scores': {
                    'trinet'       : None
                    , 'triplet'    : None
                    , 'cluster'    : None
                    , 'cut-arc set': None
                    , 'equal'      : None
                }
                , 'run times'         : {
                    'input network creation'                     : None
                    , 'input trinet set creation'                : None
                    , 'input triplet set creation'               : None
                    , 'input cluster set creation'               : None
                    , 'distorted trinet set creation'            : None
                    , 'solver creation'                          : None
                    , 'solving'                                  : None
                    , 'output trinet set creation'               : None
                    , 'output triplet set creation'              : None
                    , 'output cluster set creation'              : None
                    , 'trinet consistency score computation'     : None
                    , 'triplet consistency score computation'    : None
                    , 'cluster consistency score computation'    : None
                    , 'cut-arc set consistency score computation': None
                    , 'equal computation'                        : None
                }
            }
        assert type(data) == dict
        self.data = data
        self.name = name

    @classmethod
    def load(cls, path, name):
        with open(path) as f:
            data = json.load(f)
        return cls(name, data)

    def __setitem__(self, key, value):
        for other_key, values in self.data.items():
            other_key_keys = list(values.keys())
            if other_key == key:
                self.data[other_key] = value
                return
            elif key in other_key_keys:
                self.data[other_key][key] = value
                return
        raise KeyError

    def __getitem__(self, item):
        for key, values in self.data.items():
            if key == item:
                return self.data[key]
            elif item in values.keys():
                return self.data[key][item]
        raise KeyError

    @staticmethod
    def convert(o):
        if isinstance(o, np.int64): return int(o)
        if isinstance(o, np.int32): return int(o)
        raise TypeError

    def save(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        t = dt.datetime.today().strftime('%Y-%m-%d %H%M%S-%f')

        os.makedirs(os.path.dirname(f"{dir_path}\\{self.name}\\"), exist_ok=True)
        self['location'] = f"{dir_path}\\{self.name}\\{t}.txt"
        with open(f"{dir_path}\\{self.name}\\{t}.txt", "w") as f:
            json.dump(obj=self.data, default=self.convert, fp=f)

    @property
    def input_network(self):
        return RootedLevelKNetwork.from_enewick(self['input network'])

    @property
    def output_network(self):
        return RootedLevelKNetwork.from_enewick(self['output network'])

    def visualize_input_output(self):
        self.input_network.visualize()
        self.output_network.visualize()

    def plot_consistency_scores(self):
        data = self['consistency scores']
        plt.bar(x=data.keys(), height=data.values())
        plt.ylim([0, 1])
        plt.xticks(rotation=90)
        plt.show()

    def plot_run_times(self):
        data = self['run times']
        plt.bar(x=data.keys(), height=data.values())
        plt.xticks(rotation=90)
        plt.ylabel("Seconds")
        plt.show()

    def plot_solver_scores(self):
        data = self['solver scores']
        plt.bar(x=data.keys(), height=data.values())
        plt.ylim([0, 1])
        plt.xticks(rotation=90)
        plt.show()

    @property
    def total_run_time(self):
        return sum([t for t in self['run times'].values()])

    @property
    def solving_run_time(self):
        return self['solver creation'] + self['solving']

    def rerun(self, all_trinet_set, trinet_set_method, max_processes):
        try:
            network = RootedLevelKNetwork.from_enewick(self['input network'])
            network.visualize()

            input_trinet_set = NetworkSet(3, network_info_list=[NetworkInfo(RootedLevelKNetwork.from_enewick(enewick), multiplicity) for enewick, multiplicity
                                                                in self['input trinet set'].items()])
            input_triplet_set = NetworkSet(3, network_info_list=[NetworkInfo(RootedLevelKNetwork.from_enewick(enewick), multiplicity) for enewick, multiplicity
                                                                 in self['input triplet set'].items()])
            input_cluster_set = NetworkSet(network.number_of_leaves,
                                           network_info_list=[NetworkInfo(RootedLevelKNetwork.from_enewick(enewick), multiplicity) for enewick, multiplicity
                                                              in self['input cluster set'].items()])
            distorted_trinet_set = NetworkSet(3, network_info_list=[NetworkInfo(network=RootedLevelKNetwork.from_enewick(enewick, multiplicity)) for
                                                                    enewick, multiplicity in self['distorted trinet set'].items()])
            t6 = time.time()
            solver_params = self['solver parameters']
            solver = Solver(all_trinet_set, distorted_trinet_set, **solver_params, is_main=True)
            output_network, solver_scores = solver.solve()
            t7 = time.time()
            self['output network'] = output_network.enewick()
            self['solving'] = t7 - t6
            output_trinet_set, t8 = NetworkSet.induced_strict_network_set(output_network, 3, max_processes, False, trinet_set_method), time.time()
            self['output trinet set creation'] = t8 - t7
            self['output trinet set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                         output_trinet_set.per_network_info()}
            output_triplet_set, t9 = NetworkSet.induced_strict_tree_set(output_network, 3, 1, False), time.time()  # TODO: from trinets_set?
            self['output triplet set creation'] = t9 - t8
            self['output triplet set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                          output_triplet_set.per_network_info()}
            output_cluster_set, t10 = NetworkSet.induced_cluster_set(output_network, 1, False), time.time()
            self['output cluster set creation'] = t10 - t9
            self['output cluster set'] = {network_info.network.enewick(): network_info.multiplicity for network_info in
                                          output_cluster_set.per_network_info()}
            IO_tn_cs, t11 = NetworkSet.consistency_score(input_trinet_set, output_trinet_set), time.time()
            self['trinet'] = IO_tn_cs
            self['trinet consistency score computation'] = t11 - t10
            IO_tp_cs, t12 = NetworkSet.consistency_score(input_triplet_set, output_triplet_set), time.time()
            self['trinet'] = IO_tp_cs
            self['trinet consistency score computation'] = t12 - t11
            IO_ct_cs, t13 = NetworkSet.consistency_score(input_cluster_set, output_cluster_set), time.time()
            self['trinet'] = IO_ct_cs
            self['trinet consistency score computation'] = t13 - t12
            IO_cas_cs, t14 = network.cut_arc_set_consistency(output_network), time.time()
            self['trinet'] = IO_cas_cs
            self['trinet consistency score computation'] = t14 - t13
        except Exception as e:
            self.save()
            raise e
        self.save()


class ExperimentSet:
    def __init__(self, name, experiments=None):
        self.name = name
        self.experiments = coalesce(experiments, [])

    @classmethod
    def load(cls, folder_path, name):
        result = cls(name)
        for filename in os.listdir(folder_path):
            if filename.endswith(".txt"):
                result.experiments.append(Experiment.load(f"{folder_path}\\{filename}", name))
        return result

    def plot_average_run_time(self):
        data = dict(functools.reduce(operator.add, map(collections.Counter, [exp['run times'] for exp in self.experiments])))
        data = {key: value / len(self.experiments) for key, value in data.items()}
        plt.bar(x=data.keys(), height=data.values())
        plt.xticks(rotation=90)
        plt.ylabel("Seconds")
        plt.show()

    def plot_finished(self):
        data = {index: (experiment['init'] + experiment['full']) / 2 for index, experiment in enumerate(self.experiments)}
        plt.bar(x=data.keys(), height=data.values())
        plt.xticks(rotation=90)
        plt.ylim([0, 1])
        plt.ylabel("Finished")
        plt.show()

    def plot_average_consistency_scores(self):
        data = dict(functools.reduce(operator.add, map(collections.Counter, [exp['consistency scores'] for exp in self.experiments])))
        data = {key: value / len(self.experiments) for key, value in data.items()}
        plt.bar(x=data.keys(), height=data.values())
        plt.xticks(rotation=90)
        plt.ylim([0, 1])
        plt.ylabel("Finished")
        plt.show()

    def plot_consistency_scores(self):
        consistency_types = ['trinet', 'triplet', 'cluster', 'cut-arc set']
        x = range(len(self.experiments))
        for i, consistency_type in enumerate(consistency_types):
            plt.subplot(len(consistency_types), 1, i + 1)
            data = [exp[consistency_type] for exp in self.experiments]
            plt.bar(x=x, height=data)
            plt.ylim([0, 1])
            plt.ylabel(f'{consistency_type}')
        plt.show()

    def save(self):
        for experiment in self.experiments:
            experiment.save()

    def append(self, experiment):
        self.experiments.append(experiment)

    def __getitem__(self, item) -> Experiment:
        return self.experiments[item]

    def rerun(self, all_trinet_set, trinet_set_method, max_processes):
        for exp in self.experiments:
            exp.rerun(all_trinet_set, trinet_set_method, max_processes)
