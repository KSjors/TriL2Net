# source: https://www.researchgate.net/figure/Suggested-reticulate-evolution-events-that-could-explain-the-ETS-patterns-in_fig1_37256249

from datastructures.rooted_level_k_network import RootedLevelKNetwork

ETS_NETWORK_dict = {
    0 : [1, 2],
    1 : [3, 4],
    3 : [5, 6],
    5 : [7, 'H. angustifolius'],
    6 : [7, 'H. simulans'],
    7 : ['H. floridanus'],
    4 : [8, 9],
    8 : [10, 'H. californicus'],
    10: [11, 'H. eggertii'],
    11: ['H. hirsutus'],
    12: [11, 'H. schweinitzii'],
    9 : [12, 'H. laevigatus'],
    2 : [13, 14],
    13: [15, 'H. anomalus'],
    15: [16, 'H. pumilus'],
    16: ['H. cusikii', 'H. gracilentus'],
    14: [17, 18],
    17: [19, 20],
    19: ['H. praecox', 'H. debilis'],
    20: [23, 'H. niveus'],
    23: [24, 25],
    24: ['H. petioaris', 'H. deserticola'],
    25: ['H. paradoxis'],
    18: [26, 27],
    26: [25, 28],
    28: [29, 'H. annuus'],
    29: ['H. argophyllus'],
    27: [29, 'H. bolanderi']
}




# ETS_NETWORK.visualize(internal_node_labels=False, rankdir='LR')
#     for leaf in {'H. cusikii', 'H. praecox', 'H. petioaris', 'H. eggertii', 'H. laevigatus', 'H. schweinitzii', 'H. hirsutus', 'H. simulans', 'H. floridanus'}:
#         ETS_NETWORK.terminate_leaf(leaf)
#     ETS_NETWORK.rename_node('H. angustifolius', 'G1')
#     ETS_NETWORK.rename_node('H. californicus', 'G2')
#     ETS_NETWORK.rename_node('H. gracilentus', 'G3')
#     ETS_NETWORK.rename_node('H. debilis', 'G4')
#     ETS_NETWORK.rename_node('H. deserticola', 'G5')
#     ETS_NETWORK = RootedLevelKNetwork.from_network(ETS_NETWORK, node_names=['G1', 'G3', 'G4', 'H. anomalus'])
#     ETS_NETWORK.rename_node('G1', 'G6')
#     ETS_NETWORK.rename_node('G3', 'G7')
#     ETS_NETWORK.rename_node('G4', 'G8')
#     ETS_NETWORK.visualize(file_path='ETS_NETWORK_P2_shrunken', internal_node_labels=False, rankdir='LR')