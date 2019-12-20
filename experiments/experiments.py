import logging
from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkSet, NetworkInfo
from datastructures.solver import Solver
from utils.help_functions import *
import pandas as pd
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
import pprint

pp = pprint.PrettyPrinter(indent=4)

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
        # Set up logging
        self.uid = guid()
        self.logger_name = f"network generator.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)
        self.logger.debug("Creating new network generator")

        # parameters
        self.n = n_range
        self.r = r_range
        self.t = t_range

    def iterator(self):
        param_iterator = itertools.product(self.n, self.r, self.t)
        for parameters in param_iterator:
            self.logger.info("Generating next network")
            n, r, t = parameters
            network = RootedLevelKNetwork.random(n, r, level=2)
            if t > 0:
                network.terminate_percentage_leaves(t)
            yield parameters, network

    def __len__(self):
        return len(self.n) * len(self.r) * len(self.t)


class TrinetSetGenerator:
    def __init__(self, t, u, d, all_trinet_list):
        # Set up logging
        self.uid = guid()
        self.logger_name = f"trinet set generator.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)
        self.logger.debug("Creating new trinet set generator")

        # save parameters
        self.t = t
        self.u = u
        self.d = d

        self.all_trinet_list = all_trinet_list

    def iterator(self, trinet_set):
        param_iterator = itertools.product(self.t, self.u, self.d)

        for parameters in param_iterator:
            self.logger.info("Generating next trinet set")
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
        self.uid = guid()
        self.logger_name = f"solver generator.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)
        self.logger.debug("Creating new solver generator")

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
            self.logger.info("Generating next solver")
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
        self.uid = guid()
        self.logger_name = f"experiment generator.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)
        self.logger.debug("Creating new experiment generator")

        # Save generators
        self.ng = ng
        self.tsg = tsg
        self.sg = sg

        # Save config
        self.trinet_set_method = trinet_set_method
        self.max_processes = max_processes

        # Description and output
        self.name = name

    def iterator(self):
        t0 = time.time()
        for network_params, network in tqdm(self.ng.iterator(), total=len(self)):
            self.logger.info("Generating next experiment")
            experiment = Experiment(self.name)
            network.visualize()
            try:
                # Run
                t1 = time.time()
                self.logger.info("Computing input trinet set")
                input_trinet_set, t2 = NetworkSet.induced_strict_network_set(network, 3), time.time()
                self.logger.info("Computing input triplet set")
                input_triplet_set, t3 = NetworkSet.induced_strict_tree_set_of_network_set(input_trinet_set, 3), time.time()
                self.logger.info("Computing input cluster set")
                input_cluster_set, t4 = NetworkSet.induced_cluster_set(network), time.time()

                # Save parameters
                experiment['parameters']['input network'] = network_params

                # Save run times
                experiment['run times']['input network creation'] = t1 - t0
                experiment['run times']['input trinet set creation'] = t2 - t1
                experiment['run times']['input triplet set creation'] = t3 - t2
                experiment['run times']['input cluster set creation'] = t4 - t3

                # Save input and output
                experiment['input output']['input network'] = network.enewick()
                experiment['input output']['input trinet set'] = input_trinet_set.to_enewick()
                experiment['input output']['input triplet set'] = input_triplet_set.to_enewick()
                experiment['input output']['input cluster set'] = input_cluster_set.to_enewick()

                # Save summary
                experiment['summaries']['input network'] = network.summary()

                t4 = time.time()
                for trinet_set_params, distorted_trinet_set in self.tsg.iterator(input_trinet_set):
                    t5 = time.time()

                    # Save consistency score
                    experiment['consistency scores']['distorted trinet set'] = distorted_trinet_set.consistency_score(input_trinet_set)

                    # Save parameters
                    experiment['parameters']['distortion'] = trinet_set_params

                    # Save run time
                    experiment['run times']['distorted trinet set creation'] = t5 - t4

                    # Save input output
                    experiment['input output']['distorted trinet set'] = distorted_trinet_set.to_enewick()

                    # Save summary
                    experiment['summaries']['distorted trinet set'] = distorted_trinet_set.summary()
                    for solver_params, solver in self.sg.iterator(distorted_trinet_set):
                        t6 = time.time()
                        experiment['finished']['init'] = True
                        experiment.save("\\temp\\")
                        output_network, solver_scores = solver.solve()

                        t7 = time.time()
                        self.logger.info("Computing output trinet set")
                        output_trinet_set, t8 = NetworkSet.induced_strict_network_set(output_network, 3), time.time()
                        self.logger.info("Computing output triplet set")
                        output_triplet_set, t9 = NetworkSet.induced_strict_tree_set_of_network_set(output_trinet_set, 3), time.time()
                        self.logger.info("Computing output cluster set")
                        output_cluster_set, t10 = NetworkSet.induced_cluster_set(output_network), time.time()
                        self.logger.info("Computing trinet consistency score")
                        IO_tn_cs, t11 = input_trinet_set.consistency_score(output_trinet_set), time.time()
                        self.logger.info("Computing triplet consistency score")
                        IO_tp_cs, t12 = input_triplet_set.consistency_score(output_triplet_set), time.time()
                        self.logger.info("Computing cluster consistency score")
                        IO_ct_cs, t13 = input_cluster_set.consistency_score(output_cluster_set), time.time()
                        self.logger.info("Computing cut-arc set consistency score")
                        IO_cas_cs, t14 = network.cut_arc_set_consistency(output_network), time.time()
                        self.logger.info("Checking if output network is equal to input network")
                        equal, t15 = network.equal_structure(output_network, equal_naming=True), time.time()

                        # Save consistency scores
                        experiment['consistency scores']['trinet'] = IO_tn_cs
                        experiment['consistency scores']['triplet'] = IO_tp_cs
                        experiment['consistency scores']['cluster'] = IO_ct_cs
                        experiment['consistency scores']['cut-arc set'] = IO_cas_cs
                        experiment['consistency scores']['network'] = equal[0]

                        # Solving parameters
                        experiment['parameters']['solver'] = solver_params

                        # Save run times
                        experiment['run times']['solver creation'] = t6 - t5
                        experiment['run times']['solving'] = t7 - t6
                        experiment['run times']['output trinet set creation'] = t8 - t7
                        experiment['run times']['output triplet set creation'] = t9 - t8
                        experiment['run times']['trinet consistency score'] = t11 - t10
                        experiment['run times']['output cluster set creation'] = t10 - t9
                        experiment['run times']['triplet consistency score'] = t12 - t11
                        experiment['run times']['cluster consistency score'] = t13 - t12
                        experiment['run times']['cut-arc set consistency score'] = t14 - t13
                        experiment['run times']['equal'] = t15 - t14

                        # Save input output
                        experiment['input output']['output trinet set'] = output_trinet_set.to_enewick()
                        experiment['input output']['output triplet set'] = output_triplet_set.to_enewick()
                        experiment['input output']['output cluster set'] = output_cluster_set.to_enewick()
                        experiment['input output']['output network'] = output_network.enewick()

                        # Save summary
                        experiment['summaries']['output network'] = output_network.summary()

                        # Save finsished
                        experiment['finished']['full'] = True

                        yield experiment
                    t4 = time.time()
            except Exception as e:
                self.logger.warning(f"An error occurred: {type(e)}: {str(e)}")
                experiment['finished']['error'] = f"{type(e)}: {str(e)}"
                experiment.save(extra_path="\\failures")
                raise
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
        # Set up logging
        self.uid = guid()
        self.logger_name = f"experiment.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)
        self.logger.debug("Creating new experiment")

        if data is None:
            data = {
                'finished'            : {'init': False, 'full': False, 'location': None, 'error': None}
                , 'parameters'        : {
                    'input network': None
                    , 'distortion' : None
                    , 'solver'     : None
                }
                , 'summaries'         : {
                    'input network'         : None
                    , 'output network'      : None
                    , 'distorted trinet set': None
                }
                , 'input output'      : {
                    'input network'         : None
                    , 'input trinet set'    : None
                    , 'input triplet set'   : None
                    , 'input cluster set'   : None
                    , 'distorted trinet set': None
                    , 'output network'      : None
                    , 'output trinet set'   : None
                    , 'output triplet set'  : None
                    , 'output cluster set'  : None
                }
                , 'consistency scores': {
                    'trinet'       : None
                    , 'triplet'    : None
                    , 'cluster'    : None
                    , 'cut-arc set': None
                    , 'network'    : None
                    , 'solver'     : 0
                    , 'distorted trinet set': None
                }
                , 'run times'         : {
                    'input network creation'         : None
                    , 'input trinet set creation'    : None
                    , 'input triplet set creation'   : None
                    , 'input cluster set creation'   : None
                    , 'distorted trinet set creation': None
                    , 'solver creation'              : None
                    , 'solving'                      : None
                    , 'output trinet set creation'   : None
                    , 'output triplet set creation'  : None
                    , 'output cluster set creation'  : None
                    , 'trinet consistency score'     : None
                    , 'triplet consistency score'    : None
                    , 'cluster consistency score'    : None
                    , 'cut-arc set consistency score': None
                    , 'equal'                        : None
                }
            }
        assert type(data) == dict, "Data is of wrong type"
        self.data = data
        self.name = name

    @classmethod
    def load(cls, path, name):
        with open(path) as f:
            data = json.load(f)
        return cls(name, data)

    def __setitem__(self, key, value):
        self.data[key] = value

    def __getitem__(self, item):
        return self.data[item]

    @staticmethod
    def convert(o):
        if isinstance(o, np.int64): return int(o)
        if isinstance(o, np.int32): return int(o)
        raise TypeError("Can not convert")

    def save(self, extra_path=""):
        # if extra_path != '':
        #     raise NotImplementedError
        dir_path = os.path.dirname(os.path.realpath(__file__))
        t = dt.datetime.today().strftime('%Y-%m-%d %H%M%S-%f')

        os.makedirs(os.path.dirname(f"{dir_path}{extra_path}\\{self.name}\\"), exist_ok=True)
        self['location'] = f"{dir_path}{extra_path}\\{self.name}\\{t}.txt"
        self.logger.info(f"Saving to {dir_path}{extra_path}\\{self.name}\\{t}.txt")
        with open(f"{dir_path}{extra_path}\\{self.name}\\{t}.txt", "w") as f:
            json.dump(obj=self.data, default=self.convert, fp=f)

    def to_df(self):
        data = flatten_dictionary({key: self.data[key] for key in self.data.keys() if key != 'input output'})
        keys = list(data.keys())
        df = pd.DataFrame(columns=keys)
        df.loc[0] = [data[key] for key in keys]
        return df

    @property
    def input_network(self):
        return RootedLevelKNetwork.from_enewick(self['input output']['input network'])

    @property
    def output_network(self):
        return RootedLevelKNetwork.from_enewick(self['input output']['output network'])

    def recreate_network_set(self, category):
        self.logger.info(f"Loading {category}")
        return NetworkSet.from_enewick(self['input output'][category])

    def visualize_input_output(self):
        self.input_network.visualize()
        time.sleep(1)
        self.output_network.visualize()

    def plot_consistency_scores(self, show=True):
        data = self['consistency scores']
        plt.bar(x=data.keys(), height=data.values())
        plt.ylim([0, 1])
        plt.xticks(rotation=90)
        if show:
            plt.show()

    def plot_run_times(self):
        data = self['run times']
        plt.bar(x=data.keys(), height=data.values())
        plt.xticks(rotation=90)
        plt.ylabel("Seconds")
        plt.show()

    def plot_solver_scores(self):
        data = self['consistency scores']['solver']
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
        self.logger.info("Running experiment ...")
        try:
            self.logger.info("Loading input network")
            network = RootedLevelKNetwork.from_enewick(self['input output']['input network'])
            network.visualize()
            input_trinet_set = self.recreate_network_set('input trinet set')
            input_triplet_set = self.recreate_network_set('input triplet set')
            input_cluster_set = self.recreate_network_set('input cluster set')
            distorted_trinet_set = self.recreate_network_set('distorted trinet set')

            t6 = time.time()
            self.logger.info("Creating solver")
            solver_params = self['parameters']['solver']
            solver = Solver(all_trinet_set, distorted_trinet_set, **solver_params, is_main=True)
            output_network, solver_scores = solver.solve()
            t7 = time.time()
            output_trinet_set, t8 = NetworkSet.induced_strict_network_set(output_network, 3), time.time()
            self.logger.info("Computing output triplet set")
            output_triplet_set, t9 = NetworkSet.induced_strict_tree_set_of_network_set(output_trinet_set, 3), time.time()
            self.logger.info("Computing output cluster set")
            output_cluster_set, t10 = NetworkSet.induced_cluster_set(output_network), time.time()
            self.logger.info("Computing trinet consistency")
            IO_tn_cs, t11 = input_trinet_set.consistency_score(input_trinet_set), time.time()
            self.logger.info("Computing triplet consistency")
            IO_tp_cs, t12 = input_triplet_set.consistency_score(output_triplet_set), time.time()
            self.logger.info("Computing cluster consistency")
            IO_ct_cs, t13 = input_cluster_set.consistency_score(output_cluster_set), time.time()
            self.logger.info("Computing cut-arc set consistency")
            IO_cas_cs, t14 = network.cut_arc_set_consistency(output_network), time.time()
            self.logger.info("Checking if output network is equal to input network")
            equal, t15 = network.equal_structure(output_network, equal_naming=True), time.time()

            # Save consistencies
            self['consistency scores']['trinet'] = IO_tn_cs
            self['consistency scores']['triplet'] = IO_tp_cs
            self['consistency scores']['cluster'] = IO_ct_cs
            self['consistency scores']['cut-arc set'] = IO_cas_cs
            self['consistency scores']['network'] = equal[0]

            # Save run times
            self['run times']['solving'] = t7 - t6
            self['run times']['output trinet set creation'] = t8 - t7
            self['run times']['output triplet set creation'] = t9 - t8
            self['run times']['output cluster set creation'] = t10 - t9
            self['run times']['trinet consistency score'] = t11 - t10
            self['run times']['triplet consistency score'] = t12 - t11
            self['run times']['cluster consistency score'] = t13 - t12
            self['run times']['cut-arc set consistency score'] = t14 - t13
            self['run times']['equal'] = t15 - t14

            # Save input output
            self['input output']['output network'] = output_network.enewick()
            self['input output']['output trinet set'] = output_trinet_set.to_enewick()
            self['input output']['output triplet set'] = output_triplet_set.to_enewick()
            self['input output']['output cluster set'] = output_cluster_set.to_enewick()

            # Save summaries
            self['summaries']['output network'] = output_network.summary()

            # Save finished
            self['finished']['full'] = True
        except Exception as e:
            self.logger.warning(f"An error occurred: {type(e)}: {str(e)}")
            self['finished']['error'] = f"{type(e)}: {str(e)}"
            self.save(extra_path="\\failures")
            raise e
        self.save()


class ExperimentSet:
    def __init__(self, name, experiments=None):
        # Set up logging
        self.uid = guid()
        self.logger_name = f"experiment.{self.uid}"
        self.logger = logging.getLogger(self.logger_name)
        self.logger.debug("Creating new experiment set")

        # save config
        self.name = name
        self.experiments = coalesce(experiments, [])

    @classmethod
    def load(cls, folder_path, name):
        result = cls(name)
        for filename in os.listdir(folder_path):
            if filename.endswith(".txt"):
                result.experiments.append(Experiment.load(f"{folder_path}\\{filename}", name))
        return result

    @classmethod
    def filter(cls, experiment_set, category, value, direction: bool = True):
        result = cls(experiment_set.name)
        for experiment in experiment_set.experiments:
            if (experiment[category] == value) is direction:
                result.append(experiment)
        return result

    def to_df(self) -> pd.DataFrame:
        df = pd.concat([e.to_df() for e in self.experiments])
        return df

    def plot_average_run_time(self):
        data = dict(functools.reduce(operator.add, map(collections.Counter, [exp['run times'] for exp in self.experiments])))
        data = {key: value / len(self.experiments) for key, value in data.items()}
        plt.bar(x=data.keys(), height=data.values())
        plt.xticks(rotation=90)
        plt.ylabel("Seconds")
        plt.show()

    def plot_finished(self):
        data = {index: (experiment['finished']['init'] + experiment['finished']['full']) / 2 for index, experiment in enumerate(self.experiments)}
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
        l = len(self.experiments)
        for i, exp in enumerate(self.experiments):
            plt.subplot(1, l, i+1)
            exp.plot_consistency_scores(show=False)
        plt.show()

    def save(self):
        self.logger.info("Saving experiments ... ")
        for experiment in self.experiments:
            experiment.save()

    def append(self, experiment):
        self.experiments.append(experiment)

    def __getitem__(self, item) -> Experiment:
        return self.experiments[item]

    def rerun(self, all_trinet_set, trinet_set_method, max_processes):
        self.logger.info("Running experiments ...")
        for exp in self.experiments:
            exp.rerun(all_trinet_set, trinet_set_method, max_processes)
