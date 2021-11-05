/*  The program "l1_dist" compares two phylogenetic trees and calculates 
    L1-distance (Node distance of Williams and Clifford), 
    i.e. the mean absolute difference between path lengths for all pairs of leafs. 
    Here "path length" is for the number of branches on the path from one leaf to another.
    This metric is purely topological and does not depend on branch lengths.
    Trees assumed to be unrooted. The sets of leaf labels of two trees
    must either coincide or be embedded into each other, 
    in the latter case the distance between the smaller tree and the constraint 
    of the bigger tree on the set of leaves of the smaller tree is calculated.
    The input file format is Newick, see https://evolution.genetics.washington.edu/phylip/newick_doc.html
    If there is one input file, the distance between two first trees in this file is calculated.
    Otherwise, the distance between the first trees from two input files is calculated.
    The result is output to stdout.

    For compilation, l1_dist requires this file, treedist.c and treedist.h.

    Author: Sergei Spirin, Belozersky Institute of Moscow State University, sas@belozersky.msu.ru

    Copyright 2021 Sergei Spirin 

    This file is part of TreeDist.

    TreeDist is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TreeDist is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TreeDist (a file named "COPYING.txt"). 
    If not, see <https://www.gnu.org/licenses/>.
*/

#include "treedist.h"

int main (int argc, char *argv[])
{
  unsigned i, n;
  int distance;
  unsigned len = 128;
  FILE *inflow;
  struct tree intree1, intree2;
  char *newick;
  char c;

  /* Checking command line */
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <input tree 1> [<input tree 2>]\n", argv[0]);
    fprintf(stderr, "Example: %s tree1.tre tree2.tre\n", argv[0]);
    return 1;
  }
  if (strcmp(argv[1], "-h") == 0 ||
      strcmp(argv[1], "-help") == 0 ||
      strcmp(argv[1], "--help") == 0 ) {
    fprintf(stderr, "l1_dist is a program computing L1 distance ");
    fprintf(stderr, "between two phylogenetic trees.\n");
    fprintf(stderr, "Trees should be in Newick format.\n");
    fprintf(stderr, "Both trees can be in one file or in separate files.\n");
    fprintf(stderr, "Usage: %s <input file> [<input file 2>]\n", argv[0]);
    fprintf(stderr, "Example: %s tree1.tre tree2.tre\n", argv[0]);
    return 0;
  }

  inflow = fopen(argv[1], "r");
  if ( inflow == NULL ) {
    fprintf(stderr, "Can not open input file \"%s\"!\n", argv[1]);
    exit(1);
  }

  newick = (char*)malloc(sizeof(char) * len);
  i = 0;
  c = '\0';
  while (!feof(inflow) && (c != ';') ) {
    c = fgetc(inflow);
    if ( !feof(inflow) && c != '\n' && c != '\r') {
      i++;
      if (i > len - 2) {
        len = 2*len;
        newick = (char*)realloc(newick, sizeof(char) * len);
      }
      newick[i-1] = c;
    }
  }
  if ( feof(inflow) ) {
    fprintf(stderr, "Wrong tree format in \"%s\"!\n", argv[1]); 
    exit(1);
  }
  newick[i] = '\0';
  intree1 = readbrackets(newick);
  free(newick);

  if (argc > 2) { /* second tree is in a separate file */
    fclose(inflow);
    inflow = fopen(argv[2], "r");
    if ( inflow == NULL ) {
      fprintf(stderr, "Can not open input file \"%s\"!\n", argv[2]);
      exit(1);
    }
  }
  newick = (char*) malloc(sizeof(char)*len);
  i = 0;
  c = '\0';
  while (!feof(inflow) && (c != ';') ) {
    c = fgetc(inflow);
    if ( !feof(inflow) && c != '\n' && c != '\r') {
      i++;
      if (i > len - 2) {
        len = 2*len;
        newick = (char*)realloc(newick, sizeof(char)*len);
      }
      newick[i-1] = c;
    }
  }
  if ( feof(inflow) ) {
    if (argc > 2) {
      fprintf(stderr, "Wrong tree format in \"%s\"!\n", argv[2]);
    }
    else {
      fprintf(stderr, "Only one tree in \"%s\"!\n", argv[1]);
    }
    exit(1);
  }
  fclose(inflow);
  newick[i] = '\0';
  intree2 = readbrackets(newick);
  free(newick);

  n = intree1.leavesnum;
  if ( n == intree2.leavesnum ) {
    distance = treedist2(intree1, intree2, 1);
  }
  if ( n < intree2.leavesnum ) {
    distance = treedist2(intree1, subtree(intree2, intree1.leaf, n), 1);
  }
  if ( n > intree2.leavesnum ) {
    n = intree2.leavesnum;
    distance = treedist2(subtree(intree1, intree2.leaf, n), intree2, 1);
  }
  if ( distance >= 0 ) printf("%.4f\n", (float)(distance*2)/(n*(n-1)) );
  else puts("1000000");
  return 0;
} /* main */
