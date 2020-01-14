# source: https://www.researchgate.net/figure/Suggested-reticulate-evolution-events-that-could-explain-the-ETS-patterns-in_fig1_37256249

from datastructures.rooted_level_k_network import RootedLevelKNetwork, NetworkSet

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


# ETS_NETWORK = RootedLevelKNetwork.from_connections_dict(ETS_NETWORK_dict)
# ETS_NETWORK.terminate_leaves(
#     leaf_names_to_keep={
#         'H. angustifolius', 'H. floridanus', 'H. simulans'
#         , 'H. anomalus', 'H. pumilus', 'H. cusikii', 'H. gracilentus'
#         , 'H. praecox', 'H. debilis', 'H. niveus', 'H. petioaris'
#         , 'H. deserticola', 'H. paradoxis', 'H. annuus', 'H. bolanderi'
#         ,'H. argophyllus'})
#
# T1 = RootedLevelKNetwork.restrict(ETS_NETWORK, ['H. floridanus', 'H. simulans', 'H. niveus'])
# T2 = RootedLevelKNetwork.restrict(ETS_NETWORK, ['H. anomalus', 'H. argophyllus', 'H. deserticola'])
# T3 = RootedLevelKNetwork.restrict(ETS_NETWORK, ['H. angustifolius', 'H. praecox', 'H. deserticola'])
#
# T1.visualize(internal_node_labels=False, rankdir='LR', file_path='ETS_NETWORK_T1', format='pdf')
# T2.visualize(internal_node_labels=False, rankdir='LR', file_path='ETS_NETWORK_T2', format='pdf')
# T3.visualize(internal_node_labels=False, rankdir='LR', file_path='ETS_NETWORK_T3', format='pdf')
#

