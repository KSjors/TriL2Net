# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:48:29 2019

@author: Sjors
"""
from datastructures.rooted_level_k_generator import *
import numpy as np
from bidict import bidict
import time

# Level 0
level0 = np.array([[0, 1, 1]])
level0_symmetries = bidict()
level0_necessary = []
generator_level0 = RootedLevelKGenerator(name='0', dir_adj_matrix=level0, symmetrical_nodes=level0_symmetries, level=0)

# Level 1
level1 = np.array([[0, 2, 0], [0, 0, 1]])
level1_symmetries = bidict()
level1_necessary = []
generator_level1 = RootedLevelKGenerator(name='1', dir_adj_matrix=level1, symmetrical_nodes=level1_symmetries, level=1)

# Level 2
A = np.array([[0, 1, 1, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])
A_symmetries = bidict()
A_necessary = []
generator_A = RootedLevelKGenerator(name='2a', dir_adj_matrix=A, symmetrical_nodes=A_symmetries, level=2)

B = np.array([[0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1]])
B_symmetries = bidict()
B_necessary = []
generator_B = RootedLevelKGenerator(name='2b', dir_adj_matrix=B, symmetrical_nodes=B_symmetries, level=2)

C = np.array([[0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]])
C_symmetries = bidict({1: 2, 3: 4})
C_necessary = []
generator_C = RootedLevelKGenerator(name='2c', dir_adj_matrix=C, symmetrical_nodes=C_symmetries, level=2)

D = np.array([[0, 1, 1, 0, 0], [0, 0, 0, 2, 0], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0]])
D_symmetries = bidict()
D_necessary = [['1', '3']]
generator_D = RootedLevelKGenerator(name='2d', dir_adj_matrix=D, symmetrical_nodes=D_symmetries, level=2)
