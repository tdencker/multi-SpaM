/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef _WIN32 // change to WIN32 if error
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>


#include "axml_extract.h"



extern int Thorough;
extern infoList iList;
extern char seq_file[1024];
extern char resultFileName[1024];
extern char tree_file[1024];
extern char workdir[1024];
extern char run_id[128];
extern FILE *INFILE;
extern double masterTime;



fakeboolean initrav (tree *tr, nodeptr p)
{ 
  nodeptr  q;
  
  if (!isTip(p->number, tr->rdta->numsp)) 
    {      
      q = p->next;
      
      do 
	{	   
	  if (! initrav(tr, q->back))  return FALSE;		   
	  q = q->next;	
	} 
      while (q != p);  
      
      newviewGeneric(tr, p);
    }
  
  return TRUE;
}


fakeboolean smooth (tree *tr, nodeptr p)
{
  nodeptr  q;
  
  if (! update(tr, p))               return FALSE; /*  Adjust branch */
  if (! isTip(p->number, tr->rdta->numsp)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  if (! smooth(tr, q->back))   return FALSE;
	  q = q->next;
	}	
      
      if(tr->multiBranch)		  
	newviewGenericMasked(tr, p);	
      else
	newviewGeneric(tr, p);     
    }
  
  return TRUE;
} 

static fakeboolean allSmoothed(tree *tr)
{
  int i;
  fakeboolean result = TRUE;
  
  for(i = 0; i < tr->numBranches; i++)
    {
      if(tr->partitionSmoothed[i] == FALSE)
	result = FALSE;
      else
	tr->partitionConverged[i] = TRUE;
    }

  return result;
}



fakeboolean smoothTree (tree *tr, int maxtimes)
{
  nodeptr  p, q;   
  int i, count = 0;
   
  p = tr->start;
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;

  while (--maxtimes >= 0) 
    {    
      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;		

      if (! smooth(tr, p->back))       return FALSE;
      if (!isTip(p->number, tr->rdta->numsp)) 
	{
	  q = p->next;
	  while (q != p) 
	    {
	      if (! smooth(tr, q->back))   return FALSE;
	      q = q->next;
	    }
	}
         
      count++;

      if (allSmoothed(tr)) 
	break;      
    }

  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;



  return TRUE;
} 



fakeboolean localSmooth (tree *tr, nodeptr p, int maxtimes)
{ 
  nodeptr  q;
  int i;
  
  if (isTip(p->number, tr->rdta->numsp)) return FALSE;
  
   for(i = 0; i < tr->numBranches; i++)	
     tr->partitionConverged[i] = FALSE;	

  while (--maxtimes >= 0) 
    {     
      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;
	 	
      q = p;
      do 
	{
	  if (! update(tr, q)) return FALSE;
	  q = q->next;
        } 
      while (q != p);
      
      if (allSmoothed(tr)) 
	break;
    }

  for(i = 0; i < tr->numBranches; i++)
    {
      tr->partitionSmoothed[i] = FALSE; 
      tr->partitionConverged[i] = FALSE;
    }

  return TRUE;
}



void insertInfoList(nodeptr node, double likelihood)
{
  int i;
  int min = 0;
  double min_l =  iList.list[0].likelihood;

  for(i = 1; i < iList.n; i++)
    {
      if(iList.list[i].likelihood < min_l)
	{
	  min = i;
	  min_l = iList.list[i].likelihood;
	}
    }

  if(likelihood > min_l)
    {
      iList.list[min].likelihood = likelihood;
      iList.list[min].node = node;
      iList.valid += 1;
    }

  if(iList.valid > iList.n)
    iList.valid = iList.n;
}



fakeboolean treeEvaluate (tree *tr, double smoothFactor)       
{
  fakeboolean result;
 
  if(tr->useBrLenScaler)
    assert(0);
 
  result = smoothTree(tr, (int)((double)smoothings * smoothFactor));
  
  assert(result); 

  evaluateGeneric(tr, tr->start);         

  return TRUE;
}


/*void treeEvaluateRandom (tree *tr, double smoothFactor, analdef *adef)       */
/*{*/
/* */
/*  smoothTreeRandom(tr, (int)((double)smoothings * smoothFactor), adef);*/
/*  */

/*  evaluateGeneric(tr, tr->start);       */
/*}*/

/*void treeEvaluateProgressive(tree *tr)*/
/*{  */
/*  int */
/*    i, k;*/

/*  tr->branchCounter = 0;*/
/*  tr->numberOfBranches = 2 * tr->mxtips - 3;*/
/*  */
/*  tr->bInf = (branchInfo*)rax_malloc(tr->numberOfBranches * sizeof(branchInfo));  */

/*  setupBranches(tr, tr->start->back, tr->bInf);*/

/*  assert(tr->branchCounter == tr->numberOfBranches);*/
/*  */
/*  for(i = 0; i < tr->numBranches; i++)*/
/*    tr->partitionConverged[i] = FALSE;*/
/*  */
/*  for(i = 0; i < 10; i++)*/
/*    {*/
/*      for(k = 0; k < tr->numberOfBranches; k++)*/
/*	{      */
/*	  update(tr, tr->bInf[k].oP);*/
/*	  newviewGeneric(tr, tr->bInf[k].oP);     */
/*	}*/
/*      evaluateGenericInitrav(tr, tr->start);*/
/*      printf("It %d %f \n", i, tr->likelihood);*/
/*    }*/
/*}*/
/*   */

  



