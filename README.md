# TreeDist
TreeDist is a C package for calculating distances between phylogenetic trees.

TreeDist consists of five programs:
 * `l1_dist` calculates L1-distance (or Node distance, Williams & Clifford, 1971);
 * `l2_dis`t calculates L2-distance (or Path Difference Metric, Penny et al., 1982);
 * `quartet_dist` is a naive and slow implementation of Estabrook quartet distance (Estabrook, 1985);
 * `rf_dist` calculates normalized Robinson-Foulds distance, i.e., the fraction of different splits of two trees;
 * `rfa_dist` calculates a modified version of Robinson-Foulds distance, which is 1 minus average Jaccard measure for pairs of mutually best corresponding splits of two trees (see details in file rfa-algorithm.txt of this repository).

Input of all programs is two trees in Newick format. Trees may be in separate files or in one file.
All five programs are run from the command line:
`*_dist tree1.tre tree2.tre`
or
`*_dist twotrees.tre`
and provide their results to stdout. In the first case (two parameters in the command), first trees from both files are compared. In the second case, first two trees from the input file are regarded as input.

All programs require the set of leaves of one tree to be a subset of the set of leaves of another tree. All implemented distances are defined only for trees with equal sets of leaves. If the leaf set of one tree is a proper subset of the leaf set of another tree, the bigger tree is restricted to the smaller leaf set and the distance is computed between the obtained trees with equal leaf sets. If none of two leaf sets is a subset of another one, a senseless large value is output (it is 1000000 for most programs).

To compile the package on a Unix-like platform, put the Makefile and the src directory somewhere on your computer and run `make`. 
For Windows either use Cygwin or unzip the archive TreeDist_exe.zip.
