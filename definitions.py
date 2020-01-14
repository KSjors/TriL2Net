"""
Generators:
    0A: Cherry
    1A: Cactus
    2A: Zig-zag
    2B: Star-ship
    2C: Kite
    2D: 2-way
"""
#
# for leaf_name in leaves:
#     print()
#     print(f'placing leaf {leaf_name}')
#     leaf_index = leaf_set.inverse[leaf_name]
#
#     # Compute side alignment score
#     print(leaves_per_side)
#     print(leaf_side_alignment_matrix[leaf_index, list(placed_leaf_indicis)] for
#           side_index, placed_leaf_indicis in
#           leaves_per_side.items())
#
#     side_alignment_score = {side_index: sum(leaf_side_alignment_matrix[leaf_index, list(placed_leaf_indicis)]) for
#                             side_index, placed_leaf_indicis in
#                             leaves_per_side.items()}
#     print(f'score is {side_alignment_score}')
#
#     # Normalize this score
#     normalized_side_alignment_score = {}
#     for side_index_0, score_0 in side_alignment_score.items():
#         temp = 2 * score_0 - len(leaves_per_side[side_index_0])
#         for side_index_1, score_1 in side_alignment_score.items():
#             temp -= score_1 - len(leaves_per_side[side_index_1])
#         normalized_side_alignment_score[side_index_0] = temp
#     print(f'Normalized score is {normalized_side_alignment_score}')
#
#     # Find best score
#     best_side = max(normalized_side_alignment_score.items(), key=operator.itemgetter(1))[0]
#     # TODO: in case of ties pick set with least elements
#     leaves_per_side[best_side].add(leaf_index)
#     print(f'best side is {best_side}')