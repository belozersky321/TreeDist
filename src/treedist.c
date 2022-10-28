/*  treedist.c is a library of subroutines for Treedist package.
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
#define MAXNAME 40
#define MAXNUMLEN 25

/****************************************************************
* readbrackets: converting a string with Newick to a rooted tree 
*****************************************************************/
struct tree readbrackets(char *brackets) {
  struct tree result;
  char c;
  char tmpstr[2];
  char lenstr[MAXNUMLEN];
  char flag;
  unsigned i, j, k, s;
  unsigned nname;
  unsigned stacklen = 0, maxstacklen = 3;
  unsigned leafi, branchi, stacki, rooti, rootj; 
  char **newbranchstack;

  result.leavesnum = 0;
  result.branchnum = 0;
  result.leaf = (char**)malloc(sizeof(char*));
  result.branch = (char**)malloc(sizeof(char*));
  result.length = (float*)malloc(sizeof(float));
  newbranchstack = (char**)malloc(sizeof(char*) * maxstacklen);

  result.rooted = 0;
  result.phylogram = 0;

  i = 0;
  c = brackets[i];
  flag = 0;
  while ( (c != ';') && (c != '\0') ) {

    if ( c == ',' ) {
      if ( flag == 2 && strlen(lenstr) > 0 ) { /* branch length is ready */
        result.length[branchi] = atof(lenstr);
      }
      flag = 0; /* wait for a new leaf or new branch */
    }

    else if ( c == '(' ) { /* opening new branch */
      if (flag != 0) {
        fprintf(stderr, "Wrong Newick format\n");
        exit(1);
      }
      stacklen++;
      if (stacklen > maxstacklen) {
        maxstacklen = stacklen;
        newbranchstack = (char**)realloc(newbranchstack, sizeof(char*) * maxstacklen);
      }
      stacki = stacklen - 1;
      newbranchstack[stacki] = (char*)malloc(sizeof(char) * (result.leavesnum + 1));
      for (k = 0; k < result.leavesnum; k++) {
        newbranchstack[stacki][k] = 0;
      }
    } /* if c == '(' */

    else if ( c == ')' ) { /* new branch is ready */
      if ( flag == 2 && strlen(lenstr) > 0 ) { /* branch length is ready */
        result.length[branchi] = atof(lenstr);
      }
      flag = 3;
      if ( stacklen == 0 ) {
        fprintf(stderr, "Unbalanced brackets in Newick\n");
        exit(1);
      }
      stacklen--;
      result.branchnum++;
      branchi = result.branchnum - 1;
      result.branch = (char**)realloc(result.branch, sizeof(char*) * result.branchnum);
      result.branch[branchi] = (char*)malloc(sizeof(char) * result.leavesnum);
      for (k = 0; k < result.leavesnum; k++) {
        result.branch[branchi][k] = newbranchstack[stacklen][k];
      } /* for */
      free(newbranchstack[stacklen]);
      result.length = (float*) realloc(result.length, sizeof(float) * result.branchnum);
      result.length[branchi] = 1.0;
    }  /* if c == ')' */

    else if ( c == ':' ) {
      if ( flag == 0 || flag == 2 ) {
        fprintf(stderr, "Wrong Newick format\n");
        exit(1);
      }
      flag = 2; /* wait for a branch length */
      strcpy(lenstr, "");
    }

    else if ( flag == 0 && !isspace(c) ) { /* regard c as the first symbol of a new leaf name */
      result.leavesnum++;     /* new leaf */
      result.branchnum++;     /* new trivial (one-leaf) branch */
      leafi = result.leavesnum - 1; /* index of the current leaf */
      branchi = result.branchnum - 1; /* index of the current branch */
      result.leaf = (char **) realloc(result.leaf, sizeof(char *)*result.leavesnum);
      result.leaf[leafi] = (char *) malloc(sizeof(char) * MAXNAME);  /* memory for the name of the current leaf */
      result.leaf[leafi][0] = c;
      result.leaf[leafi][1] = '\0';
      
      result.branch = (char**)realloc(result.branch, sizeof(char*) * result.branchnum);
      result.length = (float*)realloc(result.length, sizeof(float) * result.branchnum);
      result.branch[branchi] = (char*)malloc(sizeof(char) * result.leavesnum);
      for (j = 0; j < branchi; j++) {
        result.branch[j] = (char*)realloc(result.branch[j], sizeof(char) * result.leavesnum);
        result.branch[j][leafi] = 0; /* increasing leaf array for all branches */
      }
      for ( k = 0; k < leafi; k++) {
        result.branch[branchi][k] = 0; 
      }
      result.branch[branchi][leafi] = 1; /* new branch discriminate new leaf from all old leaves */
      result.length = (float*)realloc (result.length, sizeof(float) * result.branchnum);
      result.length[branchi] = 1.0; /* length is not read yet */

      for (j = 0; j < stacklen; j++) {
        newbranchstack[j] = (char*)realloc(newbranchstack[j], sizeof(char) * result.leavesnum);
        newbranchstack[j][leafi] = 1;
      }

      flag = 1; /* in a leaf name */
      nname = 1; /* name length */
    }  /* if c is the first symbol of new leaf name */

    else if ( flag == 1 && !isspace(c) ) { /* continuation of leaf name */
      nname++;
      if (nname < MAXNAME) {
        tmpstr[0] = c;
        tmpstr[1] = '\0';
        strcat (result.leaf[leafi], tmpstr);
      } /* if */
      else {
        fprintf(stderr, "Warning: the name %s... is too long!\n", result.leaf[leafi]); 
      }
    } /* if c is a part of the name */

    else if ( flag == 2 ) { /* branch length */
      if ( c == '.' || isdigit(c) ) {
        tmpstr[0] = c;
        tmpstr[1] = '\0';
        if (strlen(lenstr) + 2 < MAXNUMLEN) {
          strcat(lenstr,tmpstr);
        }
        else {
          fprintf(stderr, "Too long number: %s...!", lenstr);
        }
      } /* if c is a part of a number (is a digit or the decimal point) */
    } /* if flag == 2 */

    i++;
    c = brackets[i];
  } /* while */
  if ( c == '\0' ) {
    fprintf(stderr, "No \';\' at the end of Newick string!\n");
    exit(1);
  }
  if ( stacklen > 0 ) {
    fprintf(stderr, "Error: unbalanced brackets in Newick\n");
    exit(1);
  }

  /* removing equal branches and the trivial branch */
  s = 0;
  i = result.branchnum - 1;
  for (k = 0; k < result.leavesnum; k++) {
    s += result.branch[i][k];
  }
  if (s == result.leavesnum) {
    result.branchnum--;
  }
  else {
    fprintf(stderr, "No outer bracket pair!\n");
    exit(1);
  }
  rooti = 0;
  rootj = 0;
  for (i = 1; i < result.branchnum && rooti == 0; i++) {
    for (j = 0; j < i && rooti == 0; j++) {
      s = 0;
      for (k = 0; k < result.leavesnum; k++) {
        if ( result.branch[i][k] == result.branch[j][k] ) {
          s++;
        }
      }
      if (s == 0) { /* branches rooti and rootj are equal <=> coincide with the root */
        rooti = i;
        rootj = j;
      }
    }
  }
  if (rooti != 0) { /* the root was found */
    result.rooted = 1;
    result.root = rooti;
    result.length[rooti] += result.length[rootj];
    result.rootlocation = result.length[rootj];
    result.branchnum--;
    for (j = rootj; j < result.branchnum; j++) { /* elimination of the branch rootj */
      result.branch[j] = result.branch[j + 1];
      result.length[j] = result.length[j + 1];
      if ( result.root == j + 1 ) result.root = j;
    }
  }

  return result;
} /* readbrackets */

/*******************************************************************
* combdistance: the number of branches (the combinatorial distance)
* between two leaves of a tree 
********************************************************************/
unsigned combdistance(struct tree intree, unsigned leaf1, unsigned leaf2) { 
  unsigned result;
  unsigned j;

  if ( leaf1 < intree.leavesnum && leaf2 < intree.leavesnum ) {
    result = 0;
    for ( j = 0; j < intree.branchnum; j++ ) {
      if ( intree.branch[j][leaf1] != intree.branch[j][leaf2] ) {
        result++;
      }
    }
  }
  else { /* at least one of indices is out of range */
    result = 0;
    fprintf(stderr, "Combdistance: input index is out of range!\n");
  }

  return result;
} /* combdistance */

/*************************************************************
*  treedist2 returns the sum of square differences, if p is 2,
*  and the sum of absulute differences, if p is 1, 
*  between combinatorial leaf-to-leaf distances.
**************************************************************/
long treedist2 (struct tree intree1, struct tree intree2, char p) {
  long result = -1;
  unsigned *corresp;
  unsigned cd1, cd2, diff;
  unsigned i, j;
  unsigned a, b;
  char flag;

  if (intree1.leavesnum == intree2.leavesnum) {
    corresp = (unsigned*)malloc(sizeof(unsigned) * intree1.leavesnum);
    for ( i = 0; i < intree1.leavesnum; i++ ) {
      flag = 1;
      for ( j = 0; j < intree2.leavesnum && flag; j++ ) {
        if ( strcmp(intree2.leaf[j], intree1.leaf[i]) == 0 ) {
          flag = 0;
          corresp[i] = j;
        }
      }
      if ( flag ) return -1;
    }
    result = 0;
    for ( a = 0; a < intree1.leavesnum - 1; a++ )
    for ( b = a + 1; b < intree1.leavesnum; b++ ) { /* for all pairs of species */
      cd1 = combdistance(intree1, a, b);
      cd2 = combdistance(intree2, corresp[a], corresp[b]);
      if (cd1 > cd2) diff = cd1 - cd2;
      else diff = cd2 - cd1;
      if ( p == 1) result += diff;
      else result += diff * diff;
    } /* for all pairs */
  } /* if */
  return result;
}  /* treedist2 */

/*******************************************************************************
* whichsplittree returns 2 if there is a split separating a and b from c and d,
* 3 if there is a split separating a and c from b and d,
* and 4 if there is a split separating a and d from b and c.
********************************************************************************/
char whichsplittree(unsigned a, unsigned b, unsigned c, unsigned d, struct tree intree) {
  char result = 0;
  unsigned i;

  if ( a >= intree.leavesnum ) {
    fprintf(stderr, "whichsplittree: a is out of range");
    exit(1);
  }
  if ( b >= intree.leavesnum ) {
    fprintf(stderr, "whichsplittree: b is out of range");
    exit(1);
  }
  if ( c >= intree.leavesnum ) {
    fprintf(stderr, "whichsplittree: c is out of range");
    exit(1);
  }
  if ( d >= intree.leavesnum ) {
    fprintf(stderr, "whichsplittree: d is out of range");
    exit(1);
  }

  for ( i = 0; i < intree.branchnum && result == 0; i++ ) {
    if ( intree.branch[i][a] == intree.branch[i][b] 
      && intree.branch[i][a] != intree.branch[i][c] 
      && intree.branch[i][c] == intree.branch[i][d] ) {
      result = 2;
    }
    else if ( intree.branch[i][a] == intree.branch[i][c] 
           && intree.branch[i][a] != intree.branch[i][b] 
           && intree.branch[i][b] == intree.branch[i][d] ) {
      result = 3;
    }
    else if ( intree.branch[i][a] == intree.branch[i][d] 
           && intree.branch[i][a] != intree.branch[i][b] 
           && intree.branch[i][b] == intree.branch[i][c] ) {
      result = 4;
    }
  }
 
  return result;
} /* whichsplittree */

/*************************************************
*  treedist4 returns the number of common fourths  
**************************************************/
long treedist4 (struct tree intree1, struct tree intree2) {
  long result;
  unsigned *corresp;
  unsigned i, j;
  unsigned a, b, c, d;
  char flag;

  if ( (intree1.leavesnum == intree2.leavesnum) && (intree1.leavesnum > 3) ) {
    corresp = (unsigned*)malloc(sizeof(unsigned) * intree1.leavesnum);
    for ( i = 0; i < intree1.leavesnum; i++ ) {
      flag = 1;
      for ( j = 0; j < intree2.leavesnum && flag; j++ ) {
        if ( strcmp(intree2.leaf[j], intree1.leaf[i]) == 0 ) {
          flag = 0;
          corresp[i] = j;
        }
      }
      if ( flag ) {
        fprintf(stderr, "Warning: no leaf with name \"%s\" in tree #2\n", intree1.leaf[i]);
        return -1;
      }
    }
    result = 0;
    for ( a = 0; a < intree1.leavesnum - 3; a++ )
    for ( b = a + 1; b < intree1.leavesnum - 2; b++ )
    for ( c = b + 1; c < intree1.leavesnum - 1; c++ )
    for ( d = c + 1; d < intree1.leavesnum; d++ ) { /* for all 4-ths of species */
      if ( whichsplittree(a, b, c, d, intree1) 
           == whichsplittree(corresp[a], corresp[b], corresp[c],  corresp[d], intree2) ) 
        result++;
    } /* for all 4-ths */
  } /* if */
  else {
    if (intree1.leavesnum <= 3) result = 0;
    else result = -1;
  }
  return result;
} /* treedist4 */

/****************************************************************
*  branchdist returns the number of different branches (splits)  
*****************************************************************/
unsigned branchdist (struct tree tree1, struct tree tree2) {
  unsigned result;
  unsigned common = 0;
  unsigned *corresp;
  unsigned i, j, k;
  char flag, iflag;

  if ( tree1.leavesnum == tree2.leavesnum ) {
    corresp = (unsigned*)malloc(sizeof(unsigned) * tree1.leavesnum);
    for ( i = 0; i < tree1.leavesnum; i++ ) {
      flag = 1;
      for ( j = 0; j < tree1.leavesnum && flag; j++ ) {
        if ( strcmp(tree1.leaf[i], tree2.leaf[j]) == 0 ) {
          flag = 0;
          corresp[i] = j;
        }
      }
      if ( flag ) {
        fprintf(stderr, "The leaf %s has no correspondence in the tree 2\n", 
                tree1.leaf[i]);
        return tree1.branchnum + tree2.branchnum - tree1.leavesnum - tree2.leavesnum;
      }
    }
 
    for ( i = 0; i < tree1.branchnum; i++ ) {
      iflag = 0;
      for ( j = 0; j < tree2.branchnum && (!iflag); j++ ) {
        flag = 0;
        for ( k = 0; k < tree1.leavesnum && flag == 0; k++ ) {
          if ( tree1.branch[i][k] != tree2.branch[j][corresp[k]] ) {
            flag = 1; /* branches are either different or of opposite orientation */
          }
        }
        if (flag == 1) { /* equal branches mean opposite orientation = completely different values of .branch[i][k] */
          for ( k = 0; k < tree1.leavesnum && flag == 1; k++ ) {
            if ( tree1.branch[i][k] == tree2.branch[j][corresp[k]] ) {
              flag = 2; /* branches are definetely different */
            }
          }
        }
        if ( flag < 2 ) { /* i and j are equal branches */
          common++;
          iflag = 1;
        }
      }
    }
    free(corresp);
    result = tree1.branchnum + tree2.branchnum - 2 * common; 
  } /* if */
  else result = tree1.branchnum + tree2.branchnum - tree1.leavesnum - tree2.leavesnum;

  return result;
} /* branchdist */

/****************************************************************
*  branchdist_n returns the number of branches (splits) 
*  in the tree with less number of branches that are not
*  presented in the tree with more branches
*****************************************************************/
unsigned branchdist_n (struct tree tree1, struct tree tree2) {
  unsigned result;
  unsigned common = 0;
  unsigned *corresp;
  unsigned i, j, k, n;
  char flag, iflag;

  if ( tree1.leavesnum == tree2.leavesnum ) {
    corresp = (unsigned*)malloc(sizeof(unsigned) * tree1.leavesnum);
    for ( i = 0; i < tree1.leavesnum; i++ ) {
      flag = 1;
      for ( j = 0; j < tree1.leavesnum && flag; j++ ) {
        if ( strcmp(tree1.leaf[i], tree2.leaf[j]) == 0 ) {
          flag = 0;
          corresp[i] = j;
        }
      }
      if ( flag ) {
        fprintf(stderr, "The leaf %s has no correspondence in the tree 2\n", 
                tree1.leaf[i]);
        return tree1.branchnum + tree2.branchnum - tree1.leavesnum - tree2.leavesnum;
      }
    }
 
    for ( i = 0; i < tree1.branchnum; i++ ) {
      iflag = 0;
      for ( j = 0; j < tree2.branchnum && (!iflag); j++ ) {
        flag = 0;
        for ( k = 0; k < tree1.leavesnum && flag == 0; k++ ) {
          if ( tree1.branch[i][k] != tree2.branch[j][corresp[k]] ) {
            flag = 1; /* branches are either different or of opposite orientation */
          }
        }
        if (flag == 1) { /* equal branches mean opposite orientation = completely different values of .branch[i][k] */
          for ( k = 0; k < tree1.leavesnum && flag == 1; k++ ) {
            if ( tree1.branch[i][k] == tree2.branch[j][corresp[k]] ) {
              flag = 2; /* branches are definetely different */
            }
          }
        }
        if ( flag < 2 ) { /* i and j are equal branches */
          common++;
          iflag = 1;
        }
      }
    }
    free(corresp);
    if(tree1.branchnum < tree2.branchnum) n = tree1.branchnum; else n = tree2.branchnum;
    result = n - common; 
  } /* if */
  else result = tree1.branchnum + tree2.branchnum - tree1.leavesnum - tree2.leavesnum;

  return result;
} /* branchdist_n */

/************************************************************************
*  aligndist: find best bidirectional hits between branches of two trees,
*  return the sum of Jaccard measures for all BBH
*************************************************************************/ 
float aligndist (struct tree tree1, struct tree tree2) {
  float result;
  float common;
  float best, curr;
  unsigned *corresp;
  unsigned **correspbranches1;
  unsigned **correspbranches2;
  unsigned *numcorrbr1;
  unsigned *numcorrbr2; 
  float *maxjacc;
  char *chosen;
  unsigned i, j, k, m;
  char flag;

  if ( tree1.leavesnum == tree2.leavesnum ) {
    corresp = (unsigned*)malloc(sizeof(unsigned) * tree1.leavesnum);
    for ( i = 0; i < tree1.leavesnum; i++ ) {
      flag = 1;
      for ( j = 0; j < tree1.leavesnum && flag; j++ ) {
        if ( strcmp(tree1.leaf[i], tree2.leaf[j]) == 0 ) {
          corresp[i] = j;   /* leaf i in tree 1 corresponds to leaf j in tree 2 */
          flag = 0;
        }
      }
      if ( flag ) {
        fprintf(stderr, "The leaf \"%s\" of tree #1 has no correspondence in the tree #2\n", tree1.leaf[i]);
        return 1.0;
      }
    }

    numcorrbr1 = (unsigned*)calloc(tree1.branchnum, sizeof(unsigned));
    numcorrbr2 = (unsigned*)calloc(tree2.branchnum, sizeof(unsigned));
    correspbranches1 = (unsigned**)malloc(sizeof(unsigned*) * tree1.branchnum);
    for ( i = 0; i < tree1.branchnum; i++ ) {
      correspbranches1[i] = (unsigned*)malloc(sizeof(unsigned));
    }
    correspbranches2 = (unsigned**)malloc(sizeof(unsigned*) * tree2.branchnum);
    for ( j = 0; j < tree2.branchnum; j++ ) {
      correspbranches2[j] = (unsigned*)malloc(sizeof(unsigned));
    }
    maxjacc = (float*)calloc(tree2.branchnum, sizeof(float));

    for ( i = 0; i < tree1.branchnum; i++ ) {
      best = 0.0;
      for ( j = 0; j < tree2.branchnum; j++ ) {
        curr =  jaccard(tree1.branch[i], tree2.branch[j], corresp, tree1.leavesnum);
        if ( curr > best) {
          best = curr;
          numcorrbr1[i] = 1;
          correspbranches1[i][0] = j;
        }
        if ( curr == best ) {
          numcorrbr1[i]++;
          correspbranches1[i] = (unsigned*)realloc(correspbranches1[i], sizeof(unsigned) * numcorrbr1[i]);
          correspbranches1[i][numcorrbr1[i] - 1] = j;
        }
        if ( curr > maxjacc[j] ) {
          maxjacc[j] = curr;
          numcorrbr2[j] = 1;
          correspbranches2[j][0] = i;
        }
        if ( curr == maxjacc[j] ) {
          numcorrbr2[j]++;
          correspbranches2[j] = (unsigned*)realloc(correspbranches2[j], sizeof(unsigned) * numcorrbr2[j]);
          correspbranches2[j][numcorrbr2[j] - 1] = i;
        }
      }
    }

    chosen = (char*)calloc(tree1.branchnum, sizeof(char));
    common = 0.0;
    for ( j = 0; j < tree2.branchnum; j++ ) {
      flag = 1;
      for ( k = 0; k < numcorrbr2[j] && flag; k++ ) {
        i = correspbranches2[j][k];
        if ( !chosen[i] ) {
          for ( m = 0; m < numcorrbr1[i] && flag; m++ ) {
            if ( j == correspbranches1[i][m] ) {
              flag = 0;
              chosen[i] = 1; /* BBH for i is found */
              common += maxjacc[j]; 
            }
          }
        }
      }
    }
    free(corresp);
    free(numcorrbr1); 
    free(numcorrbr2);
    for ( i = 0; i < tree1.branchnum; i++ ) free(correspbranches1[i]);
    for ( j = 0; j < tree1.branchnum; j++ ) free(correspbranches2[j]);
    free(correspbranches1); 
    free(correspbranches2); 
    free(maxjacc);
    result = tree1.branchnum + tree2.branchnum - 2 * common; 
  } /* if ( tree1.leavesnum == tree2.leavesnum ) */
  else result = -1.0;
  return result;
} /* aligndist */

/****************************************************
* ffminf and ffmaxf are min and max for float
* added for compatibility between compilators
*****************************************************/
float ffminf(float a, float b) {
  if (a > b) return b;
  return a;
}
float ffmaxf(float a, float b) {
  if (a > b) return a;
  return b;
}

/*********************************************************
* jaccard: return the Jaccard measure of two splits
**********************************************************/
float jaccard(char *br1, char *br2, unsigned *corresp, unsigned n) {
  float res00, res01, res10, res11;
  unsigned is00 = 0, is01 = 0, is10 = 0, is11 = 0;
  unsigned un00 = 0, un01 = 0, un10 = 0, un11 = 0;
  unsigned i, j;

  for ( i = 0; i < n; i++ ) {
    j = corresp[i];
    if ( br1[i] == 0 && br2[j] == 0) is00++;
    if ( br1[i] == 0 && br2[j] == 1) is01++;
    if ( br1[i] == 1 && br2[j] == 0) is10++;
    if ( br1[i] == 1 && br2[j] == 1) is11++;
    if ( br1[i] == 0 || br2[j] == 0) un00++;
    if ( br1[i] == 0 || br2[j] == 1) un01++;
    if ( br1[i] == 1 || br2[j] == 0) un10++;
    if ( br1[i] == 1 || br2[j] == 1) un11++;
  }
  res00 = ((float)is00)/un00;
  res01 = ((float)is01)/un01;
  res10 = ((float)is10)/un10;
  res11 = ((float)is11)/un11;

  return ffmaxf(ffminf(res00, res11), ffminf(res01, res10));
} /* jaccard */

/*************************************************************
* subtree:
*  for a given tree and a list of names create the subtree
*  on species contained in the list
**************************************************************/
struct tree subtree(struct tree intree, char **leaflist, unsigned int listlen) {
  unsigned i, j, k;
  char flag;
  char *newbranches;
  char *newleaves;
  unsigned newbranchnum, countleaf;
  unsigned *correspleaf;
  unsigned *correspbranch;
  struct tree result;

  result.rooted = intree.rooted;
  result.leavesnum = 0;
  result.leaf = (char**)malloc(sizeof(char*) * listlen);
  correspleaf = (unsigned*)malloc(sizeof(unsigned) * listlen);
  correspbranch = (unsigned*)malloc(sizeof(unsigned) * intree.branchnum);
  newleaves = (char*)calloc(intree.leavesnum, sizeof(char));
  i = j = 0;
  while ( i < listlen) {
    flag = 0;
    for ( k = 0; k < intree.leavesnum && flag == 0; k++) {
      if ( strcmp(intree.leaf[k], leaflist[i]) == 0 ) {
        flag = 1;
        newleaves[k] = 1;
        result.leavesnum++;
        result.leaf[j] = (char*)malloc(sizeof(char) * (strlen(leaflist[i]) + 1));
        strcpy(result.leaf[j], leaflist[i]);
        correspleaf[j] = k;
        j++;
      } /* if */
    } /* for k */
    i++;
  } /* while i < listlen */

  fprintf(stderr, "Restricting tree from %u to %u leafs\n", intree.leavesnum, listlen);

  newbranches = (char*)calloc(intree.branchnum, sizeof(char));
  newbranchnum = 0;
  for ( i = 0; i < intree.branchnum; i++ ) {
    newbranches[i] = 1; 
    countleaf = 0;
    for ( k = 0; k < intree.leavesnum; k++ ) {
      if ( newleaves[k] ) countleaf += intree.branch[i][k];
    }
    if ( countleaf > 0 && countleaf < result.leavesnum ) {
      flag = 3;
      for ( j = 0; j < i && flag == 3; j++ ) {
        if ( newbranches[j] ) {
          flag = 0;
          for ( k = 0; k < intree.leavesnum && flag < 3; k++ ) {
            if ( newleaves[k] ) {
              if ( flag == 0 && intree.branch[i][k] == intree.branch[j][k] ) flag = 1;
              else if ( flag == 0 && intree.branch[i][k] != intree.branch[j][k] ) flag = 2;
              else if ( flag == 1 && intree.branch[i][k] != intree.branch[j][k] ) flag = 3;
              else if ( flag == 2 && intree.branch[i][k] == intree.branch[j][k] ) flag = 3;
            } /* if */
          } /* for k */
          if ( flag < 3 ) {
            newbranches[i] = 0;
            correspbranch[i] = correspbranch[j];
          }
        } /* if newbranches[j] == 1 */
      } /* for j */
    } /* if */
    else { 
      newbranches[i] = 0;
      correspbranch[i] = intree.branchnum; 
    }
    if ( newbranches[i] ) {
      correspbranch[i] = newbranchnum;
      newbranchnum++;
    }
  }
  result.branchnum = newbranchnum;
  result.branch = (char**)malloc(sizeof(char*) * newbranchnum);

  i = j = 0;
  while ( i < newbranchnum ) {
    if ( newbranches[j] ) {
      result.branch[i] = (char*)calloc(result.leavesnum, sizeof(char));
      for ( k = 0; k < result.leavesnum; k++ ) {
        result.branch[i][k] = intree.branch[j][correspleaf[k]];
      }
      i++;
    }
    j++;
  }

  if ( intree.rooted ) {
    if ( correspbranch[intree.root] < intree.branchnum ) {
      result.root = correspbranch[intree.root];
    }
    else {
      result.root = 0; /* !!! Subject to change !!! */
    }
  }

  result.phylogram = intree.phylogram;
  result.length = (float*)calloc(result.branchnum, sizeof(float));
  if ( intree.phylogram ) {
    for ( j = 0; j < intree.branchnum; j++ ) {
      if ( correspbranch[j] < intree.branchnum ) result.length[correspbranch[j]] += intree.length[j];
    }
    if ( intree.rooted ) {
      result.rootlocation = intree.rootlocation; /* !!! Subject to change !!! */
    }
  }
  else {
    for ( i = 0; i < result.branchnum; i++ ) {
      result.length[i] = 1.0;
    }
    result.rootlocation = 1.0;
  }

  return result;
} /* subtree */
