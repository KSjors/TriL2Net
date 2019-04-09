# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:48:29 2019

@author: Sjors
"""
from datastructures.rooted_level_k_network import *

cherry = np.array([[0, 1, 1, 0, 0], [0, 0, 0, 1, 1]])

cactus = np.array([[0, 2, 0], [0, 0, 1]])

A = np.array([[0, 1, 1, 0, 0], [0, 0, 1, 1, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])

B = np.array([[0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1]])

C = np.array([[0, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 1, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]])

D = np.array([[0, 1, 1, 0, 0], [0, 0, 0, 2, 0], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0]])

# network_A = network_from_dir_adj_matrix(A).standard_form()
# network_B = network_from_dir_adj_matrix(B).standard_form()
# network_C = network_from_dir_adj_matrix(C).standard_form()
# network_D = network_from_dir_adj_matrix(D).standard_form()


