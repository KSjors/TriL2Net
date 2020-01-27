DEFAULT = 0

#   Leaf locator
GREEDY = 0
ILP = 1

#   Counting Methods
MAXIMUM_MULTIPLICITY = 0
WEIGHTED_AVERAGE = 1
WEIGHTED_SUM = 2

#    Thresholds
DEFAULT_THRESHOLD = 0

#    Leaf order methods
DEFAULT_ORDER = 0

#    Minimal sink-set
FIRST_SCC_THAT_IS_MSS = 0
EXPAND_FIRST_SCC = 1

#    Trinet Compute Method
ITERATIVE = 0
RECURSIVE = 1

# All trinet list


FALSE = 0
TRUE = 1

# Formats
FORMAT_eNewick = 0
FORMAT_eNewick_multiplicity = 1
FORMAT_tnets = 2
# FORMAT_trilonet_input = 3

FORMAT_extension = {
    FORMAT_eNewick               : 'eNewick'
    , FORMAT_eNewick_multiplicity: 'eNewickM'
    , FORMAT_tnets               : 'tnets'
    # , FORMAT_trilonet_input: 'tnets'
}
