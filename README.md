# TreeDist
TreeDist is a C package for calculating distances between phylogenetic trees.

TreeDist consists of five programs:
 * l1_dist calculates L1-distance (or Node distance, Williams & Clifford, 1971);
 * l2_dist calculates L2-distance (or Path Difference Metric, Penny et al., 1982);
 * quartet_dist is a naive and slow implementation of calculation of Estabrook quartet distance (Estabrook, 1985);
 * rf_dist calculates the normalized Robinson-Foulds distance, i.e., the fraction of different splits of two trees;
 * rfa_dist calculated a modified version of Robinson-Foulds distance, which is 1 minus average Jaccard measure for pairs of mutually best corresponding splits of two trees.

Input of all programs is two trees in Newick format. Trees may be in separate files or in one file.
All five programs are run from the command line:
`*_dist tree1.tre tree2.tre`
or
`*_dist twotrees.tre`
and provide their results to stdout. In the first case (two paramters in the command), first trees from both files are compared. In the second case, first two trees from the input files are regarded as input.
