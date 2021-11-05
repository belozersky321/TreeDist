/*  treedist.h is a header file for TreeDist package.
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

#ifndef _HEADER_H_
#define _HEADER_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <ctype.h>

/* Structure for tree */
struct tree {
  unsigned leavesnum; /* number of leaves */
  unsigned branchnum; /* number of branches, including trivial ones */
  char **leaf; /* leaves' names */
  char **branch; /* branches as strings of 0 and 1 of length leavesnum */
  char rooted; /* 1 if rooted, 0 if unrooted */
  unsigned root; /* index of the root branch */
  char phylogram; /* 1 if lengths are really lengths, 0 otherwise */
  float *length; /* branch lengths */
  float rootlocation; /* distance from the root to the node from the side of leaf #1 (index 0) */
};

struct tree readbrackets(char *brackets); /* converting string with Newick to rooted tree */

unsigned combdistance(struct tree intree, unsigned leaf1, unsigned leaf2); 
/* combinatorial distance (number of branches in path) between two leaves  */
long treedist2 (struct tree intree1, struct tree intree2, char p);

unsigned branchdist(struct tree tree1, struct tree tree2);

float aligndist(struct tree tree1, struct tree tree2);
float ffminf(float a, float b);
float ffmaxf(float a, float b);
float jaccard(char *br1, char *br2, unsigned *corresp, unsigned n);

char whichsplittree(unsigned a, unsigned b, unsigned c, unsigned d, 
                        struct tree intree);
long treedist4 (struct tree intree1, struct tree intree2);

struct tree subtree(struct tree intree, char **leaflist, unsigned listlen);

#endif
