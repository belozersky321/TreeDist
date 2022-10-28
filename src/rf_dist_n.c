/*  The program "rf_dist" compares two phylogenetic trees and prints 
    the normalized Robinson-Foulds distance, i.e. the fraction of different branches (splits), to the stdout.
    Trees assumed to be unrooted. The sets of leaf labels of two trees
    must either coincide or be embedded into each other, 
    in the latter case the distance between the smaller tree and the constraint 
    of the bigger tree on the set of leaves of the smaller tree is calculated.
    The input file format is Newick, see https://evolution.genetics.washington.edu/phylip/newick_doc.html
    If there is one input file, the distance between two first trees in this file is calculated.
    Otherwise, the distance between the first trees from two input files is calculated.

    For compilation, rf_dist requires this file, treedist.c and treedist.h.

    Author: Sergei Spirin, Belozersky Institute of Moscow State University, sas@belozersky.msu.ru

    Copyright 2021 Sergei Spirin 

    This file is part of Treedist.

    Treedist is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Treedist is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Treedist (a file named "COPYING.txt"). 
    If not, see <https://www.gnu.org/licenses/>.
*/

#include "treedist.h"

int main (int argc, char *argv[])
{
  unsigned i, n;
  unsigned distance;
  unsigned len = 128;
  FILE *inflow;
  struct tree intree1, intree2, tmptree;
  char *newick;
  char c;

  /* Checking command line */
  if (argc < 2)
  {
    fprintf(stderr, "Usage: %s <input tree 1> [<input tree 2>]\n", argv[0]);
    fprintf(stderr, "Example: %s tree1.tre tree2.tre\n", argv[0]);
    return 1;
  }
  if (strcmp(argv[1], "-h") == 0 ||
      strcmp(argv[1], "-help") == 0 ||
      strcmp(argv[1], "--help") == 0 ) {
    fprintf(stderr, "rf_dist is a program computing the split (Robinson-Foulds) ");
    fprintf(stderr, "distance between two phylogenetic trees.\n");
    fprintf(stderr, "Trees should be in Newick format.\n");
    fprintf(stderr, "Both trees can be in one file or in separate files.\n");
    fprintf(stderr, "Usage: %s <input file> [<input file 2>]\n", argv[0]);
    fprintf(stderr, "Example: %s tree1.tre tree2.tre\n", argv[0]);
    fprintf(stderr, "Example: %s twotrees.tre\n", argv[0]);
    return 0;
  }

  inflow = fopen(argv[1], "r");
  if ( inflow == NULL ) {
    fprintf(stderr, "Can not open input file \"%s\"!\n", argv[1]);
    exit(1);
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
    fprintf(stderr, "Wrong tree format in \"%s\"!\n", argv[1]);
    exit(1);
  }
  newick[i] = '\0';
  intree1 = readbrackets(newick);
  free(newick);

  if (argc > 2) { /* second tree is in separate file */
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
    return 1;
  }
  fclose(inflow);
  newick[i] = '\0';
  intree2 = readbrackets(newick);
  free(newick);
  if (intree1.leavesnum == intree2.leavesnum) {
    distance = branchdist_n(intree1, intree2);
    if(intree1.branchnum < intree2.branchnum) n = intree1.branchnum - intree1.leavesnum;
    else n = intree2.branchnum - intree2.leavesnum;
  }
  if (intree1.leavesnum < intree2.leavesnum) {
    tmptree = subtree(intree2, intree1.leaf, intree1.leavesnum);
    distance = branchdist(intree1, tmptree);
    if(intree1.branchnum < tmptree.branchnum) n = intree1.branchnum - intree1.leavesnum;
    else n = tmptree.branchnum - tmptree.leavesnum;
  }
  if (intree1.leavesnum > intree2.leavesnum) {
    tmptree = subtree(intree1, intree2.leaf, intree2.leavesnum);
    distance = branchdist(tmptree, intree2);
    if(tmptree.branchnum < intree2.branchnum) n = tmptree.branchnum - tmptree.leavesnum;
    else n = intree2.branchnum - intree2.leavesnum;
  }
  
  if(n == 0) printf("0.0000\n");
  else printf("%.4f\n", (float)distance/n );
  return 0;
} /* main */
