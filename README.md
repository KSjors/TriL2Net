# TriLoNet-2

## How to use:
Requisites:
- The TriL2Net repository
- Python 3.6+, including the packages: numpy, scipy, mip, inputimeout, multiprocessing, bidict, graphviz, tqdm, tarjan
    (These can be installed using pip: pip install numpy scipy mip inputimeout multiprocessing bidict graphviz tqdm tarjan_


There are currently two ways to let TriL2Net compute a network from a set of trinets:
- from a list of eNewick strings
    + each line represents a trinet according to the eNewick format:
        e.g. ((b,(d)#H1),(c,#H1))root;
    + see https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-532 for information on this format
    + see example\example.eNewick for an example
- from a list of "named" level-1 trinets (this will be extended to level-2)
    + each line represents one of the level-1 trinets:
        e.g. Tr0 = Tr{a, b, c: T1} 
    + the last value (T1) represents the type of level-1 trinet, see Figure 9 of https://academic.oup.com/mbe/article/33/8/2151/2578738 for the list of all these trinets and the corresponding names, the first three values (a,b,c) represent the names of the leaves
    + see example\example.tnet for an example
    
To use any of these two methods, open a terminal and navigate to the TriL2Net folder containing eNewick.py and tnet.py. Then enter:
   "python eNewick.py example\example.eNewick"
 or
   "python tnet.py example\example.tnet"
This instructs TriL2Net to find a network on the given data set and save the resulting network in eNewick format in output.eNewick. A copy of the logs will be saved in output.log. 

Another output file can be specified using 
   "python eNewick.py example\example.eNewick other_output_file"
This will save the network and log in other_output_file.eNewick and other_output_file.log respectively.

Parameters for the solver can be set using:
   "python eNewick.py example\example.eNewick other_output_file parameters_file"
This parameters file must be in json format. An example file is example\default_parameters.json. This contains the default settings used by TriL2Net.
