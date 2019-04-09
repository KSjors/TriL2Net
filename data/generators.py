# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:48:29 2019

@author: Sjors
"""
from datastructures.rooted_level_k_generator import *
from bidict import bidict
import time

# Level 0
level0 = np.array([[0, 1, 1, 0, 0], [0, 0, 0, 1, 1]])
level0_symmetries = bidict()
level0_necessary = []
generator_level0 = RootedLevelKGenerator(dir_adj_matrix=level0, symmetrical_nodes=level0_symmetries, necessary_edges=level0_necessary)

# Level 1
level1 = np.array([[0, 2, 0], [0, 0, 1]])
level1_symmetries = bidict()
level1_necessary = [['0', '1']]
generator_level1 = RootedLevelKGenerator(dir_adj_matrix=level1, symmetrical_nodes=level1_symmetries, necessary_edges=level1_necessary)

# Level 2
A = np.array([[0, 1, 1, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
A_symmetries = bidict()
A_necessary = []
generator_A = RootedLevelKGenerator(dir_adj_matrix=A, symmetrical_nodes=A_symmetries, necessary_edges=A_necessary)

B = np.array([[0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1]])
B_symmetries = bidict()
B_necessary = []
generator_B = RootedLevelKGenerator(dir_adj_matrix=B, symmetrical_nodes=B_symmetries, necessary_edges=B_necessary)

C = np.array([[0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]])
C_symmetries = bidict({'1': '2', '3': '4'})
C_necessary = []
generator_C = RootedLevelKGenerator(dir_adj_matrix=C, symmetrical_nodes=C_symmetries, necessary_edges=C_necessary)

D = np.array([[0, 1, 1, 0, 0], [0, 0, 0, 2, 0], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0]])
D_symmetries = bidict()
D_necessary = [['1', '3']]
generator_D = RootedLevelKGenerator(dir_adj_matrix=D, symmetrical_nodes=D_symmetries, necessary_edges=D_necessary)
