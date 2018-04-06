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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses
 *  with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#include "axml_extract.h"
#include "globalVariables.h"
#include <string.h>
#include <float.h>

unsigned char getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}

void errorExit(int e)
{
  exit(e);
}

void set_custom_options(analdef * adef, tree * tr, int nbr_runs, fakeboolean gamma_model)
{
    /*********** tr inits **************/
    #ifdef _USE_PTHREADS
      NumberOfThreads = 0;
    #endif
      tr->doSubtreeEPA = FALSE;
      tr->useFastScaling = TRUE; 
      tr->bootStopCriterion = -1;
      tr->wcThreshold = 0.03;
      tr->doCutoff = TRUE;
      tr->secondaryStructureModel = SEC_16; /* default setting */
      tr->searchConvergenceCriterion = FALSE;
      tr->catOnly = FALSE;
      tr->useEpaHeuristics = FALSE;
      tr->fastEPAthreshold = -1.0;
      tr->multiStateModel  = GTR_MULTI_STATE;
      tr->saveMemory = FALSE;
      tr->useGammaMedian = FALSE;
      tr->noRateHet = FALSE;
      tr->perPartitionEPA = FALSE;
      tr->useBrLenScaler = FALSE;
      tr->ascertainmentCorrectionType = NOT_DEFINED;
      tr->autoProteinSelectionType = AUTO_ML;
      //EPA related stuff 
      tr->numberOfEPAEntries = 7;
      tr->accumulatedEPACutoff = 0.95;
      tr->useAccumulatedEPACutoff = FALSE;
      tr->probThresholdEPA = 0.01;
      //JC and K80 and HKY85
      tr->useK80 = FALSE;
      tr->useJC69 = FALSE;
      tr->useHKY85 = FALSE;

    #ifdef _BASTIEN
      tr->doBastienStuff = FALSE;
    #endif
      /********* tr inits end*************/
      
      // no warnings and checks for identical sequences and undetermined sites
      adef->silent = TRUE; 
      adef->noSequenceCheck = TRUE;
      adef->sampleQuartetsWithoutReplacement = TRUE;
      
      #ifdef _USE_PTHREADS
      printf("Error: Trying to use pthreads within raxml\n");
      exit(1);
      #endif
      
      adef->gapyness = 0;
      
      adef->parsimonySeed = 12345;
      
      // -f q -# 1000
      adef->mode = QUARTET_CALCULATION;
      adef->multipleRuns = nbr_runs;
      
      // model -m GTRGAMMA
      if(gamma_model)
        adef->model = M_GTRGAMMA;
      else
        adef->model = M_GTRCAT;
      adef->useInvariant = FALSE;
}

void initAdef(analdef *adef)
{  
  adef->useSecondaryStructure  = FALSE;
  adef->bootstrapBranchLengths = FALSE;
  adef->model                  = M_GTRCAT;
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = adef->bestTrav = 10;
  adef->initialSet             = FALSE;
  adef->restart                = FALSE;
/*  adef->mode                   = BIG_RAPID_MODE;*/
  adef->categories             = 25;
  adef->boot                   = 0;
  adef->rapidBoot              = 0;
  adef->useWeightFile          = FALSE;
  adef->checkpoints            = 0;
  adef->startingTreeOnly       = 0;
  adef->multipleRuns           = 1;
  adef->useMultipleModel       = FALSE;
  adef->likelihoodEpsilon      = 0.1;
  adef->constraint             = FALSE;
  adef->grouping               = FALSE;
  adef->randomStartingTree     = FALSE;
  adef->parsimonySeed          = 0;
  adef->constraintSeed         = 0;
  adef->proteinMatrix          = JTT;
  adef->protEmpiricalFreqs     = 0;
  adef->outgroup               = FALSE;
  adef->useInvariant           = FALSE;
  adef->permuteTreeoptimize    = FALSE;
  adef->useInvariant           = FALSE;
  adef->allInOne               = FALSE;
  adef->likelihoodTest         = FALSE;
  adef->perGeneBranchLengths   = FALSE;
  adef->generateBS             = FALSE;
  adef->bootStopping           = FALSE;
  adef->gapyness               = 0.0;
  adef->similarityFilterMode   = 0;
  adef->useExcludeFile         = FALSE;
  adef->userProteinModel       = FALSE;
  adef->computeELW             = FALSE;
  adef->computeDistance        = FALSE;
  adef->compressPatterns       = TRUE; 
  adef->readTaxaOnly           = FALSE; 
  adef->useBinaryModelFile     = FALSE;
  adef->leaveDropMode          = FALSE;
  adef->slidingWindowSize      = 100;
  adef->checkForUndeterminedSequences = TRUE;
  adef->useQuartetGrouping = FALSE;
  adef->alignmentFileType = FASTA;
  adef->calculateIC = FALSE;
  adef->verboseIC = FALSE;
  adef->stepwiseAdditionOnly = FALSE;
  adef->optimizeBaseFrequencies = FALSE;
  adef->ascertainmentBias = FALSE;
  adef->rellBootstrap = FALSE;
  adef->mesquite = FALSE;
  adef->silent = FALSE;
  adef->noSequenceCheck = FALSE;
  adef->useBFGS = TRUE;
  adef->setThreadAffinity = FALSE;
  adef->bootstopPermutations = 100;
  adef->fcThreshold = 99;
  adef->sampleQuartetsWithoutReplacement = FALSE;
}

static int iterated_bitcount(unsigned int n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}

static char bits_in_16bits [0x1u << 16];
#define BIT_COUNT(x)  precomputed16_bitcount(x)
unsigned int precomputed16_bitcount (unsigned int n)
{
  /* works only for 32-bit int*/
    
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}

static void compute_bits_in_16bits(void)
{
    unsigned int i;    

    assert(sizeof(unsigned int) == 4);

    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);
    
    return ;
}

double randum (int64_t  *seed)
{
  int64_t  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
    y = 362436069,
    z = 21288629,
    w = 14921776,
    c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}

static fakeboolean setupTree (tree *tr, analdef *adef)
{
  nodeptr  p0, p, q;
  int
    i,
    j,  
    tips,
    inter; 
  
  
  
  tr->storedBrLens = (double*)NULL;

  if(!adef->readTaxaOnly)
    {
      tr->bigCutoff = FALSE;

      tr->patternPosition = (int*)NULL;
      tr->columnPosition = (int*)NULL;

      tr->maxCategories = MAX(4, adef->categories);

      tr->partitionContributions = (double *)rax_malloc(sizeof(double) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->partitionContributions[i] = -1.0;

      tr->perPartitionLH = (double *)rax_malloc(sizeof(double) * tr->NumberOfModels);
      tr->storedPerPartitionLH = (double *)rax_malloc(sizeof(double) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  tr->perPartitionLH[i] = 0.0;
	  tr->storedPerPartitionLH[i] = 0.0;
	}

      if(adef->grouping)
	tr->grouped = TRUE;
      else
	tr->grouped = FALSE;

      if(adef->constraint)
	tr->constrained = TRUE;
      else
	tr->constrained = FALSE;

      tr->treeID = 0;
    }

  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  if(!adef->readTaxaOnly)
    {
      tr->yVector      = (unsigned char **)  rax_malloc((tr->mxtips + 1) * sizeof(unsigned char *));     

      tr->likelihoods  = (double *)rax_malloc(adef->multipleRuns * sizeof(double));
    }

  tr->numberOfTrees = -1;

  tr->treeStringLength = 
    2 * (size_t)tr->mxtips + //parentheses
    2 * (size_t)tr->mxtips * 64 + //branche lengths with : and . and branch labels and
    (size_t)tr->mxtips + //commas
    1 + //closing semicolon 
    (size_t)tr->mxtips * nmlngth; //taxon names

  //tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  //printf("tips %d Tree String Length %d old length %d\n", tr->mxtips, tr->treeStringLength,tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2 );

  tr->tree_string  = (char*)rax_calloc(tr->treeStringLength, sizeof(char)); 

  if(!adef->readTaxaOnly)
    {
           
      tr->td[0].count = 0;
      tr->td[0].ti    = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * tr->mxtips);	

     
      tr->constraintVector = (int *)rax_malloc((2 * tr->mxtips) * sizeof(int));

      tr->nameList = (char **)rax_malloc(sizeof(char *) * (tips + 1));
    }

  if (!(p0 = (nodeptr) rax_malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }

  if (!(tr->nodep = (nodeptr *) rax_malloc((2*tr->mxtips) * sizeof(nodeptr))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
    {
      p = p0++;

      p->hash   =  KISS32(); /* hast table stuff */
      p->x      =  0;
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL;
      p->bInf   = (branchInfo *)NULL;

      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++)
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++)
	{	 
	  p = p0++;
	  if(j == 1)
	    p->x = 1;
	  else
	    p->x =  0;
	  p->number = i;
	  p->next   = q;
	  p->bInf   = (branchInfo *)NULL;
	  p->back   = (node *) NULL;
	  p->hash   = 0;

	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;

  

  tr->ntips       = 0;
  tr->nextnode    = 0;

  if(!adef->readTaxaOnly)
    {
      for(i = 0; i < tr->numBranches; i++)
	tr->partitionSmoothed[i] = FALSE;
    }

  return TRUE;
}

void getyspace (rawdata *rdta)
{
  size_t size = 4 * ((size_t)(rdta->sites / 4 + 1));
  int    i;
  unsigned char *y0;

  rdta->y = (unsigned char **) rax_malloc((rdta->numsp + 1) * sizeof(unsigned char *));
  assert(rdta->y);   

  y0 = (unsigned char *) rax_malloc(((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char));
  assert(y0);   

  rdta->y0 = y0;

  for (i = 0; i <= rdta->numsp; i++)
    {
      rdta->y[i] = y0;
      y0 += size;
    }

  return;
}

stringHashtable *initStringHashTable(hashNumberType n)
{
  /* 
     init with primes 
  */
    
  static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
 

  /* init with powers of two

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};
  */
  
  stringHashtable *h = (stringHashtable*)rax_malloc(sizeof(stringHashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];  

  h->table = (stringEntry**)rax_calloc(tableSize, sizeof(stringEntry*));
  h->tableSize = tableSize;    

  return h;
}

static hashNumberType  hashString(char *p, hashNumberType tableSize)
{
  hashNumberType h = 0;
  
  for(; *p; p++)
    h = 31 * h + *p;
  
  return (h % tableSize);
}

void addword(char *s, stringHashtable *h, int nodeNumber)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  

  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return;	  	
    }

  p = (stringEntry *)rax_malloc(sizeof(stringEntry));

  assert(p);
  
  p->nodeNumber = nodeNumber;
  p->word = (char *)rax_malloc((strlen(s) + 1) * sizeof(char));

  strcpy(p->word, s);
  
  p->next =  h->table[position];
  
  h->table[position] = p;
}

static void sitesort(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  
    gap, 
    i, 
    j, 
    jj, 
    jg, 
    k, 
    n, 
    nsp,  
    *index, 
    *category = (int*)NULL;

  fakeboolean  
    flip, 
    tied;
  
  unsigned char  
    **data;

/*  if(adef->useSecondaryStructure)*/
/*    {*/
/*      assert(tr->NumberOfModels > 1 && adef->useMultipleModel);*/

/*      adaptRdataToSecondary(tr, rdta);*/
/*    }*/

  if(adef->useMultipleModel)    
    category      = tr->model;
  

  index    = cdta->alias;
  data     = rdta->y;
  n        = rdta->sites;
  nsp      = rdta->numsp;
  index[0] = -1;


  if(adef->compressPatterns)
    {
      for (gap = n / 2; gap > 0; gap /= 2)
	{
	  for (i = gap + 1; i <= n; i++)
	    {
	      j = i - gap;

	      do
		{
		  jj = index[j];
		  jg = index[j+gap];
/*		  if(adef->useMultipleModel)*/
/*		    {		     		      */
/*		      assert(category[jj] != -1 &&*/
/*			     category[jg] != -1);*/
/*		     */
/*		      flip = (category[jj] > category[jg]);*/
/*		      tied = (category[jj] == category[jg]);		     */
/*		    }*/
/*		  else*/
/*		    {*/
		      flip = 0;
		      tied = 1;
/*		    }*/
		  for (k = 1; (k <= nsp) && tied; k++)
		    {
		      flip = (data[k][jj] >  data[k][jg]);
		      tied = (data[k][jj] == data[k][jg]);
		    }

		  if (flip)
		    {
		      index[j]     = jg;
		      index[j+gap] = jj;
		      j -= gap;
		    }
		}
	      while (flip && (j > 0));
	    }
	}
    }
}


static void sitecombcrunch (rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef, int countAscBias)
{
  fakeboolean  
    tied;
  
  int   
    i, 
    sitei, 
    j, 
    sitej, 
    k,
    *aliasModel = (int*)NULL,
    *aliasSuperModel = (int*)NULL,
    undeterminedSites = 0;

 

  if(adef->useMultipleModel)
    {
      aliasSuperModel = (int*)rax_malloc(sizeof(int) * (rdta->sites + 1));
      aliasModel      = (int*)rax_malloc(sizeof(int) * (rdta->sites + 1));
    } 

  i = 0;
  cdta->alias[0]    = cdta->alias[1];
  cdta->aliaswgt[0] = 0;
  
 //if(adef->mode == PER_SITE_LL || adef->mode == ANCESTRAL_STATES) // not commented out by me
  {
    int i;
    
    tr->patternPosition = (int*)rax_malloc(sizeof(int) * rdta->sites);
    tr->columnPosition  = (int*)rax_malloc(sizeof(int) * rdta->sites);
    
    for(i = 0; i < rdta->sites; i++)
      {
	tr->patternPosition[i] = -1;
	tr->columnPosition[i]  = -1;
      }
  }

  i = 0;

  for (j = 1; j <= rdta->sites; j++)
    {
      int 
	allGap = TRUE;

      unsigned char 
	undetermined;
      
      sitei = cdta->alias[i];
      sitej = cdta->alias[j];

      undetermined = getUndetermined(tr->dataVector[sitej]);

      for(k = 1; k <= rdta->numsp; k++)
	{	 
	  if(rdta->y[k][sitej] != undetermined)
	    {
	      allGap = FALSE;
	      break;
	    }
	}

      if(allGap)      
	undeterminedSites++;

      if(!adef->compressPatterns)
	    tied = 0;
      else
	{
	  if(adef->useMultipleModel)
	    {	     
	      tied = (tr->model[sitei] == tr->model[sitej]);
	      if(tied)
		assert(tr->dataVector[sitei] == tr->dataVector[sitej]);
	    }
	  else
	  {
	    tied = 1;
	  }
	}

      for (k = 1; tied && (k <= rdta->numsp); k++)
	{
	tied = (rdta->y[k][sitei] == rdta->y[k][sitej]);
	}

      assert(!(tied && allGap));
      
      if(tied && !allGap)
	{	  
	  tr->patternPosition[j - 1] = i;
	  tr->columnPosition[j - 1] = sitej;
	  
	  cdta->aliaswgt[i] += rdta->wgt[sitej];

	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
      else
	{
	  if(!allGap)
	    {
	      if(cdta->aliaswgt[i] > 0) 
		i++;
	      	      
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	
	      cdta->aliaswgt[i] = rdta->wgt[sitej];
	      cdta->alias[i] = sitej;
	      
	      if(adef->useMultipleModel)
		{
		  aliasModel[i]      = tr->model[sitej];
		  aliasSuperModel[i] = tr->dataVector[sitej];
		}
	    }	
	}
    }

  cdta->endsite = i;

  if (cdta->aliaswgt[i] > 0) 
    cdta->endsite++;
    
  if( (countAscBias > 0)) // adef->mode == PER_SITE_LL || adef->mode == ANCESTRAL_STATES ||
    {
        assert(undeterminedSites == 0);
/*      if(undeterminedSites > 0)*/
/*	{*/
/*	  printBothOpen("You are trying to infer per site likelihoods or ancestral states or\n");*/
/*	  printBothOpen("do calculations with an ascertainment bias correction\n");*/
/*	  printBothOpen("on an alignment containing %d sites consisting only of undetermined\n", undeterminedSites);*/
/*	  printBothOpen("characters. Please remove them first and then re-run RAxML!\n");*/

/*	  errorExit(-1);*/
/*	}*/

      for(i = 0; i < rdta->sites; i++)
	{
	  int 
	    p  = tr->patternPosition[i],
	    c  = tr->columnPosition[i];
	  	  
	  assert(p >= 0 && p < cdta->endsite);
	  assert(c >= 1 && c <= rdta->sites);
	}
    }

 
 

  if(adef->useMultipleModel)
    {
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->model[i]      = aliasModel[i];
	  tr->dataVector[i] = aliasSuperModel[i];
	}    
    }

  if(adef->useMultipleModel)
    {
      rax_free(aliasModel);
      rax_free(aliasSuperModel);
    }     

    assert(undeterminedSites == 0);
/*  if(undeterminedSites > 0)    */
/*    printBothOpen("\nAlignment has %d completely undetermined sites that will be automatically removed from the input data\n\n", undeterminedSites);*/

  //exit(-1);
}

static fakeboolean makeweights (analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr, int countAscBias)
{
  int  i;

  for (i = 1; i <= rdta->sites; i++)
    cdta->alias[i] = i;

  sitesort(rdta, cdta, tr, adef);
  sitecombcrunch(rdta, cdta, tr, adef, countAscBias);

  return TRUE;
}

static fakeboolean makevalues(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  i, j, model, fullSites = 0, modelCounter;

  unsigned char
    *y    = (unsigned char *)rax_malloc(((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char)),
    *yBUF = (unsigned char *)rax_malloc( ((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));

  for (i = 1; i <= rdta->numsp; i++)
    for (j = 0; j < cdta->endsite; j++)
	{
      y[(((size_t)(i - 1)) * ((size_t)cdta->endsite)) + j] = rdta->y[i][cdta->alias[j]];
	}

  rax_free(rdta->y0);
  rax_free(rdta->y);

  rdta->y0 = y;
  memcpy(yBUF, y, ((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));
  rdta->yBUF = yBUF;

  if(!adef->useMultipleModel)
    tr->NumberOfModels = 1;

  if(adef->useMultipleModel)
    {
      tr->partitionData[0].lower = 0;

      model        = tr->model[0];
      modelCounter = 0;
     
      i            = 1;

      while(i <  cdta->endsite)
	{	  
	  if(tr->model[i] != model)
	    {	     
	      tr->partitionData[modelCounter].upper     = i;
	      tr->partitionData[modelCounter + 1].lower = i;

	      model = tr->model[i];	     
	      modelCounter++;
	    }
	  i++;
	}

      if(modelCounter <  tr->NumberOfModels - 1)
	{
	  printf("\nYou specified %d partitions, but after parsing and pre-processing ExaML only found %d partitions\n", tr->NumberOfModels, modelCounter + 1);
	  printf("Presumably one or more partitions vanished because they consisted entirely of undetermined characters.\n");
	  printf("Please fix your data!\n\n");
	  exit(-1);
	}

      tr->partitionData[tr->NumberOfModels - 1].upper = cdta->endsite;      
    
      for(i = 0; i < tr->NumberOfModels; i++)		  
	tr->partitionData[i].width      = tr->partitionData[i].upper -  tr->partitionData[i].lower;
	 
      model        = tr->model[0];
      modelCounter = 0;
      tr->model[0] = modelCounter;
      i            = 1;
	
      while(i < cdta->endsite)
	{	 
	  if(tr->model[i] != model)
	    {
	      model = tr->model[i];
	      modelCounter++;
	      tr->model[i] = modelCounter;
	    }
	  else
	    tr->model[i] = modelCounter;
	  i++;
	}      
    }
  else
    {
      tr->partitionData[0].lower = 0;
      tr->partitionData[0].upper = cdta->endsite;
      tr->partitionData[0].width =  tr->partitionData[0].upper -  tr->partitionData[0].lower;
    }

  tr->rdta       = rdta;
  tr->cdta       = cdta;

  tr->invariant          = (int *)rax_malloc(cdta->endsite * sizeof(int));
  tr->originalDataVector = (int *)rax_malloc(cdta->endsite * sizeof(int));
  tr->originalModel      = (int *)rax_malloc(cdta->endsite * sizeof(int));
  tr->originalWeights    = (int *)rax_malloc(cdta->endsite * sizeof(int));

  memcpy(tr->originalModel, tr->model,            cdta->endsite * sizeof(int));
  memcpy(tr->originalDataVector, tr->dataVector,  cdta->endsite * sizeof(int));
  memcpy(tr->originalWeights, tr->cdta->aliaswgt, cdta->endsite * sizeof(int));


  tr->originalCrunchedLength = tr->cdta->endsite;
  for(i = 0; i < tr->cdta->endsite; i++)
    fullSites += tr->cdta->aliaswgt[i];

  tr->fullSites = fullSites;

  for(i = 0; i < rdta->numsp; i++)
    tr->yVector[i + 1] = &(rdta->y0[((size_t)tr->originalCrunchedLength) * ((size_t)i)]);

  return TRUE;
}

static void setupPresenceMask(tree *tr)
{
  int 
    model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {     
      int 
	j;

      for(j = 1; j <= tr->mxtips; j++)
	{
	  unsigned int 
	    presenceMask = 0,	  
	    i;
	    
	  for(i = 0; i < tr->partitionData[model].width; i++)	      	   	 	             	  		  	    	 		
	    presenceMask = presenceMask | mask32[tr->partitionData[model].yVector[j][i]];       
	  
	  tr->partitionData[model].presenceMap[j] = presenceMask;
	}
    }  
}

int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}

static fakeboolean isCat(analdef *adef)
{
  if(adef->model == M_PROTCAT || adef->model == M_GTRCAT || adef->model == M_BINCAT || adef->model == M_32CAT || adef->model == M_64CAT)
    return TRUE;
  else
    return FALSE;
}

static void setRateHetAndDataIncrement(tree *tr, analdef *adef)
{
  int model;

  if(isCat(adef))
    tr->rateHetModel = CAT;
  else
    {
      if(adef->useInvariant)
	tr->rateHetModel = GAMMA_I;
      else
	tr->rateHetModel = GAMMA;
    }

  switch(tr->rateHetModel)
    {
    case GAMMA:
    case GAMMA_I:
      tr->discreteRateCategories = 4;      
      break;
    case CAT:
      if((adef->boot && !adef->bootstrapBranchLengths) || (adef->mode == CLASSIFY_ML) || (tr->catOnly))	
	tr->discreteRateCategories = 1; 	
      else
	tr->discreteRateCategories = 4;
      break;
    default:
      assert(0);
    }

  if(adef->bootstrapBranchLengths)
    assert(tr->discreteRateCategories == 4);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	states = -1,
	maxTipStates = getUndetermined(tr->partitionData[model].dataType) + 1;
    assert(tr->partitionData[model].dataType == DNA_DATA);
	  states = getStates(tr->partitionData[model].dataType);	 

      tr->partitionData[model].states       = states;
      tr->partitionData[model].maxTipStates = maxTipStates;
    }
}

partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
    states    = p->states,
    tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  pLength.leftLength = pLength.rightLength = states * states;
  pLength.eignLength = states -1;
  pLength.evLength   = states * states;
  pLength.eiLength   = states * states - states;
  pLength.substRatesLength = (states * states - states) / 2;
  pLength.frequenciesLength = states;
  pLength.tipVectorLength   = tipLength * states;
  pLength.symmetryVectorLength = (states * states - states) / 2;
  pLength.frequencyGroupingLength = states;
  pLength.nonGTR = FALSE;
  //pLength.optimizeBaseFrequencies = FALSE;

  return (&pLengths[dataType]); 
}

static void allocPartitions(tree *tr)
{
  int
    i,
    maxCategories = tr->maxCategories;

  for(i = 0; i < tr->NumberOfModels; i++)
    {
      const partitionLengths 
	*pl = getPartitionLengths(&(tr->partitionData[i]));
            
      if(tr->useFastScaling)	
	tr->partitionData[i].globalScaler    = (unsigned int *)rax_calloc(2 * tr->mxtips, sizeof(unsigned int));  	         

      
      tr->partitionData[i].left              = (double *)rax_malloc(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[i].right             = (double *)rax_malloc(pl->rightLength * (maxCategories + 1) * sizeof(double));      
      tr->partitionData[i].EIGN              = (double*)rax_malloc(pl->eignLength * sizeof(double));
      tr->partitionData[i].EV                = (double*)rax_malloc(pl->evLength * sizeof(double));
      tr->partitionData[i].EI                = (double*)rax_malloc(pl->eiLength * sizeof(double));
      tr->partitionData[i].substRates        = (double *)rax_malloc(pl->substRatesLength * sizeof(double));
      tr->partitionData[i].frequencies       = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[i].freqExponents     = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[i].tipVector         = (double *)rax_malloc(pl->tipVectorLength * sizeof(double));

      tr->partitionData[i].invariableFrequencies  = (double *)rax_malloc(pl->states * sizeof(double));      
      

#ifdef _HET
      
      tr->partitionData[i].EIGN_TIP          = (double*)rax_malloc(pl->eignLength * sizeof(double));
      tr->partitionData[i].EV_TIP            = (double*)rax_malloc(pl->evLength * sizeof(double));
      tr->partitionData[i].EI_TIP            = (double*)rax_malloc(pl->eiLength * sizeof(double));
      tr->partitionData[i].substRates_TIP    = (double *)rax_malloc(pl->substRatesLength * sizeof(double));      
      tr->partitionData[i].tipVector_TIP     = (double *)rax_malloc(pl->tipVectorLength * sizeof(double));

#endif
      


/*      if(tr->partitionData[i].protModels == LG4 || tr->partitionData[i].protModels == LG4X)      */
/*	{	  	  */
/*	  int */
/*	    k;*/
/*	  */
/*	  for(k = 0; k < 4; k++)*/
/*	    {	    */
/*	      tr->partitionData[i].EIGN_LG4[k]              = (double*)rax_malloc(pl->eignLength * sizeof(double));*/
/*	      tr->partitionData[i].rawEIGN_LG4[k]              = (double*)rax_malloc(pl->eignLength * sizeof(double));	      */
/*	      tr->partitionData[i].EV_LG4[k]                = (double*)rax_malloc(pl->evLength * sizeof(double));*/
/*	      tr->partitionData[i].EI_LG4[k]                = (double*)rax_malloc(pl->eiLength * sizeof(double));*/
/*	      tr->partitionData[i].substRates_LG4[k]        = (double *)rax_malloc(pl->substRatesLength * sizeof(double));*/
/*	      tr->partitionData[i].frequencies_LG4[k]       = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));*/
/*	      tr->partitionData[i].tipVector_LG4[k]         = (double *)rax_malloc(pl->tipVectorLength * sizeof(double));*/
/*	    }*/
/*	}*/

      

      tr->partitionData[i].symmetryVector    = (int *)rax_malloc(pl->symmetryVectorLength  * sizeof(int));
      tr->partitionData[i].frequencyGrouping = (int *)rax_malloc(pl->frequencyGroupingLength  * sizeof(int));
      tr->partitionData[i].perSiteRates      = (double *)rax_malloc(sizeof(double) * tr->maxCategories);
      tr->partitionData[i].unscaled_perSiteRates = (double *)rax_malloc(sizeof(double) * tr->maxCategories);
      
      
      tr->partitionData[i].nonGTR = FALSE;     
      

      tr->partitionData[i].gammaRates = (double*)rax_malloc(sizeof(double) * 4);
      tr->partitionData[i].yVector = (unsigned char **)rax_malloc(sizeof(unsigned char*) * (tr->mxtips + 1));

           
      tr->partitionData[i].xVector = (double **)rax_malloc(sizeof(double*) * tr->innerNodes);     
      tr->partitionData[i].xSpaceVector = (size_t *)rax_calloc(tr->innerNodes, sizeof(size_t));	
           
      tr->partitionData[i].expVector      = (int **)rax_malloc(sizeof(int*) * tr->innerNodes);
      tr->partitionData[i].expSpaceVector = (size_t *)rax_calloc(tr->innerNodes, sizeof(size_t));

      tr->partitionData[i].mxtips  = tr->mxtips;

     
      //andre-opt
      tr->partitionData[i].presenceMap = (unsigned int *)rax_calloc((size_t)tr->mxtips + 1 , sizeof(unsigned int));

/*#ifndef _USE_PTHREADS    */
      {
	int j;

	for(j = 1; j <= tr->mxtips; j++)
	  tr->partitionData[i].yVector[j] = &(tr->yVector[j][tr->partitionData[i].lower]);
      }
/*#endif*/

    }
}

static void allocNodex (tree *tr)
{
  size_t
    i,   
    model,
    offset,
    memoryRequirements = 0;

  allocPartitions(tr);

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      size_t 
	width = tr->partitionData[model].upper - tr->partitionData[model].lower;

      int 
	undetermined, 
	j;

      memoryRequirements += (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states) * width;              
	
      //asc
          
     assert(!tr->partitionData[model].ascBias);
/*      if(tr->partitionData[model].ascBias)*/
/*	{	 */
/*	  tr->partitionData[model].ascOffset = 4 * tr->partitionData[model].states * tr->partitionData[model].states;*/

/*	  tr->partitionData[model].ascVector = (double *)rax_malloc(((size_t)tr->innerNodes) **/
/*								    ((size_t)tr->partitionData[model].ascOffset) * */
/*								    sizeof(double));*/
/*	  	 */
/*	  tr->partitionData[model].ascExpVector = (int *)rax_calloc(((size_t)tr->innerNodes) * ((size_t)tr->partitionData[model].states),*/
/*								 sizeof(int));*/
/*	 */
/*	  */
/*	  tr->partitionData[model].ascSumBuffer = (double *)rax_malloc(((size_t)tr->partitionData[model].ascOffset) **/
/*								       sizeof(double));*/
/*	  */
/*	}*/
      
      //asc

      tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
      
      tr->partitionData[model].gapVector = (unsigned int*)rax_calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));


      tr->partitionData[model].initialGapVectorSize = tr->partitionData[model].gapVectorLength * 2 * tr->mxtips * sizeof(int);
	
      /* always multiply by 4 due to frequent switching between CAT and GAMMA in standard RAxML */
      
      tr->partitionData[model].gapColumn = (double *)rax_malloc(((size_t)tr->innerNodes) *
								    ((size_t)4) * 
								    ((size_t)(tr->partitionData[model].states)) *
								    sizeof(double));		  		
	
      

      undetermined = getUndetermined(tr->partitionData[model].dataType);

      for(j = 1; j <= tr->mxtips; j++)
	for(i = 0; i < width; i++)
	  if(tr->partitionData[model].yVector[j][i] == undetermined)
	    tr->partitionData[model].gapVector[tr->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];      
    }

  tr->perSiteLL       = (double *)rax_malloc((size_t)tr->cdta->endsite * sizeof(double));
  assert(tr->perSiteLL != NULL);

  tr->sumBuffer  = (double *)rax_malloc(memoryRequirements * sizeof(double));
  assert(tr->sumBuffer != NULL);
 
  offset = 0;

  /* C-OPT for initial testing tr->NumberOfModels will be 1 */

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      size_t 
	lower = tr->partitionData[model].lower,
	width = tr->partitionData[model].upper - lower;

      /* TODO all of this must be reset/adapted when fixModelIndices is called ! */

      
      tr->partitionData[model].sumBuffer       = &tr->sumBuffer[offset];
     

      tr->partitionData[model].perSiteLL    = &tr->perSiteLL[lower];        


      tr->partitionData[model].wgt          = &tr->cdta->aliaswgt[lower];
      tr->partitionData[model].invariant    = &tr->invariant[lower];
      tr->partitionData[model].rateCategory = &tr->cdta->rateCategory[lower];

      offset += (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states) * width;      
    }

  for(i = 0; i < tr->innerNodes; i++)
    {     
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	{	 	  
	  tr->partitionData[model].expVector[i] = (int*)NULL;
	  tr->partitionData[model].xVector[i]   = (double*)NULL;		  			      		  		      		      		  	    	    	  	 
	}
    }
}

void getinput_begin(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  int i;
  assert(rdta->numsp > 0 && rdta->sites > 0);
  tr->mxtips            = rdta->numsp;
  rdta->wgt             = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int)); 
  cdta->alias           = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));
  cdta->aliaswgt        = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));
  cdta->rateCategory    = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));
  tr->model             = (int *)    rax_calloc((rdta->sites + 1), sizeof(int));
  tr->initialDataVector  = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));
  tr->extendedDataVector = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));     
  cdta->patrat          = (double *) rax_malloc((rdta->sites + 1) * sizeof(double));
  cdta->patratStored    = (double *) rax_malloc((rdta->sites + 1) * sizeof(double));      
  
  compute_bits_in_16bits();

  for (i = 1; i <= rdta->sites; i++)
    rdta->wgt[i] = 1;

  tr->multiBranch = 0;
  tr->numBranches = 1;

    assert(!adef->useMultipleModel);
/*      if(adef->useMultipleModel)*/
/*	{*/
/*	  int ref;*/
/*	  */
/*	  parsePartitions(adef, rdta, tr);	  */
/*	  */
/*	  for(i = 1; i <= rdta->sites; i++)*/
/*	    {*/
/*	      ref = tr->model[i];*/
/*	      tr->initialDataVector[i] = tr->initialPartitionData[ref].dataType;*/
/*	    }*/
/*	}*/
/*      else*/
	{
	  int 
	    dataType = -1;
	  
	  tr->initialPartitionData  = (pInfo*)rax_malloc(sizeof(pInfo));
	  tr->initialPartitionData[0].partitionName = (char*)rax_malloc(128 * sizeof(char));
	  strcpy(tr->initialPartitionData[0].partitionName, "No Name Provided");
	  
	  tr->initialPartitionData[0].protModels = adef->proteinMatrix;
	  if(adef->protEmpiricalFreqs)
	    tr->initialPartitionData[0].usePredefinedProtFreqs  = FALSE;
	  else
	    tr->initialPartitionData[0].usePredefinedProtFreqs  = TRUE;
	  
	  if(adef->optimizeBaseFrequencies)
	    {	     
	      tr->initialPartitionData[0].optimizeBaseFrequencies  = TRUE;
	      tr->initialPartitionData[0].usePredefinedProtFreqs  = FALSE;
	    }
	  else
	    tr->initialPartitionData[0].optimizeBaseFrequencies  = FALSE;

	  if(adef->ascertainmentBias)
	    tr->initialPartitionData[0].ascBias = TRUE;
	  else
	    tr->initialPartitionData[0].ascBias = FALSE;


	  tr->NumberOfModels = 1;
	  
	  if(adef->model == M_PROTCAT || adef->model == M_PROTGAMMA)
	    dataType = AA_DATA;
	  if(adef->model == M_GTRCAT || adef->model == M_GTRGAMMA)
	    dataType = DNA_DATA;
	  if(adef->model == M_BINCAT || adef->model == M_BINGAMMA)
	    dataType = BINARY_DATA;
	  if(adef->model == M_32CAT || adef->model == M_32GAMMA)
	    dataType = GENERIC_32;
	  if(adef->model == M_64CAT || adef->model == M_64GAMMA)
	    dataType = GENERIC_64;
	     
	     
/*	  assert(dataType == BINARY_DATA || dataType == DNA_DATA || dataType == AA_DATA || */
/*		 dataType == GENERIC_32  || dataType == GENERIC_64);*/

        assert(dataType == DNA_DATA); /* aa stuff cut out */

	  tr->initialPartitionData[0].dataType = dataType;
	  
	  for(i = 0; i <= rdta->sites; i++)
	    {
	      tr->initialDataVector[i] = dataType;
	      tr->model[i]      = 0;
	    }
	}
	  tr->dataVector    = tr->initialDataVector;
	  tr->partitionData = tr->initialPartitionData;
	  
      tr->executeModel   = (fakeboolean *)rax_malloc(sizeof(fakeboolean) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->executeModel[i] = TRUE;

      getyspace(rdta);

  setupTree(tr, adef);

      assert(adef->alignmentFileType == FASTA); /* removed PHYLIP stuff*/
	  //parseFasta(adef, rdta, tr);
      
	// -> end
}

void getinput_end(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
    int countAscBias = 0, i = 0;

	for(i = 0; i < tr->NumberOfModels; i++)	  
	  if(tr->partitionData[i].ascBias)
	    countAscBias++;

	makeweights(adef, rdta, cdta, tr, countAscBias);
	makevalues(rdta, cdta, tr, adef); 
	
	tr->innerNodes = tr->mxtips;
	setRateHetAndDataIncrement(tr, adef);
	
	assert(countAscBias == 0);
/*	checkAscBias(tr);*/
	
	if(!adef->readTaxaOnly)  
      allocNodex(tr);
      
    assert(tr->ascertainmentCorrectionType == NOT_DEFINED);
/*  readAscFiles(tr);*/
  
  setupPresenceMask(tr);
}

fakeboolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}

static void markNodesInTree(nodeptr p, tree *tr, unsigned char *nodesInTree)
{
  if(isTip(p->number, tr->mxtips))
    nodesInTree[p->number] = 1;
  else
    {
      markNodesInTree(p->next->back, tr, nodesInTree);
      markNodesInTree(p->next->next->back, tr, nodesInTree);
    }

}

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}

static nodeptr findAnyTipFast(nodeptr p, int numsp)
{ 
  return  (p->number <= numsp)? p : findAnyTipFast(p->next->back, numsp);
} 

static void checkSeed(analdef *adef)
{ 
  static fakeboolean seedChecked = FALSE;

  if(!seedChecked) 
    {
      /*printf("Checking seed\n");*/

      if(adef->parsimonySeed <= 0)
	{
	  printf("Error: you need to specify a random number seed with \"-p\" for the randomized stepwise addition\n");
	  printf("parsimony algorithm or random tree building algorithm such that runs can be reproduced and debugged ... exiting\n");      
	}
  
      assert(adef->parsimonySeed > 0);
      seedChecked = TRUE;
    }
}

static void getxnodeLocal (nodeptr p)
{
  nodeptr  s;

  if((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }
}

void getxnode (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }

  assert(p->x);
}

static void computeTraversalInfoParsimony(nodeptr p, int *ti, int *counter, int maxTips, fakeboolean full)
{        
  nodeptr 
    q = p->next->back,
    r = p->next->next->back;
  
  if(! p->x)
    getxnodeLocal(p);  
  
  if(full)
    {
       if(q->number > maxTips) 
	 computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  else
    {
      if(q->number > maxTips && !q->x) 
	computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips && !r->x) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  
  
  ti[*counter]     = p->number;
  ti[*counter + 1] = q->number;
  ti[*counter + 2] = r->number;
  *counter = *counter + 4;
}

static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q, tr->numBranches);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 

static void insertParsimony (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q, tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches); 
   
  newviewParsimony(tr, p);     
} 

static void buildSimpleTree (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupDefault(p, tr->nodep[iq], tr->numBranches);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertParsimony(tr, s, p);
}

void makePermutation(int *perm, int lower, int n, analdef *adef)
{    
  int  i, j, k;

  checkSeed(adef);          

  for (i = lower; i <= n; i++)    
    perm[i] = i;               

  for (i = lower; i <= n; i++) 
    {    
      k =  (int)((double)(n + 1 - i) * randum(&adef->parsimonySeed));

      assert(i + k <= n);
      
      j        = perm[i];
      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}



#if (defined(__SIM_SSE3) || defined(__AVX))

static inline unsigned int populationCount(INT_TYPE v_N)
{
  unsigned int
    res[INTS_PER_VECTOR] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  
  unsigned int 
    i,
    a = 0;
    
  VECTOR_STORE((CAST)res, v_N);
    
  for(i = 0; i < INTS_PER_VECTOR; i++)
    a += BIT_COUNT(res[i]);
    
  return a;	   
}

#else

static inline unsigned int populationCount(unsigned int n)
{
  return BIT_COUNT(n);
}

#endif

#if (defined(__SIM_SSE3) || defined(__AVX))

void newviewParsimonyIterativeFast(tree *tr)
{    
  INT_TYPE
    allOne = SET_ALL_BITS_ONE;

  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 

  for(index = 4; index < count; index += 4)
    {      
      unsigned int
	totalScore = 0;

      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  size_t
	    k,
	    states = tr->partitionData[model].states,
	    width = tr->partitionData[model].parsimonyLength;	 
            
	  unsigned int	
	    i;      
                 
	  switch(states)
	    {
	    case 2:       
	      {
		parsimonyNumber
		  *left[2],
		  *right[2],
		  *thisOne[2];

		for(k = 0; k < 2; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 2 * rNumber) + width * k]);
		    thisOne[k]  = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    INT_TYPE
		      s_r, s_l, v_N,
		      l_A, l_C,
		      v_A, v_C;	    	 
		    
		    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
		    l_A = VECTOR_BIT_AND(s_l, s_r);
		    v_A = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
		    l_C = VECTOR_BIT_AND(s_l, s_r);
		    v_C = VECTOR_BIT_OR(s_l, s_r);		  		  		  		  
		    
		    v_N = VECTOR_BIT_OR(l_A, l_C);
		    
		    VECTOR_STORE((CAST)(&thisOne[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
		    VECTOR_STORE((CAST)(&thisOne[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);		  
		  }
	      }
	      break;
	    case 4:
	      {
		parsimonyNumber
		  *left[4],
		  *right[4],
		  *thisOne[4];

		for(k = 0; k < 4; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 4 * rNumber) + width * k]);
		    thisOne[k]  = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    INT_TYPE
		      s_r, s_l, v_N,
		      l_A, l_C, l_G, l_T,
		      v_A, v_C, v_G, v_T;	    	 
		    
		    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
		    l_A = VECTOR_BIT_AND(s_l, s_r);
		    v_A = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
		    l_C = VECTOR_BIT_AND(s_l, s_r);
		    v_C = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[2][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[2][i]));
		    l_G = VECTOR_BIT_AND(s_l, s_r);
		    v_G = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[3][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[3][i]));
		    l_T = VECTOR_BIT_AND(s_l, s_r);
		    v_T = VECTOR_BIT_OR(s_l, s_r);
		    
		    v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));	  	 	    	  
		    
		    VECTOR_STORE((CAST)(&thisOne[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
		    VECTOR_STORE((CAST)(&thisOne[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));
		    VECTOR_STORE((CAST)(&thisOne[2][i]), VECTOR_BIT_OR(l_G, VECTOR_AND_NOT(v_N, v_G)));
		    VECTOR_STORE((CAST)(&thisOne[3][i]), VECTOR_BIT_OR(l_T, VECTOR_AND_NOT(v_N, v_T)));	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);	
		  }
	      }
	      break;
	    case 20:
	      {
		parsimonyNumber
		  *left[20],
		  *right[20],
		  *thisOne[20];

		for(k = 0; k < 20; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 20 * rNumber) + width * k]);
		    thisOne[k]  = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    size_t j;
		    
		    INT_TYPE
		      s_r, s_l, 
		      v_N = SET_ALL_BITS_ZERO,
		      l_A[20], 
		      v_A[20];	    	 
		    
		    for(j = 0; j < 20; j++)
		      {
			s_l = VECTOR_LOAD((CAST)(&left[j][i]));
			s_r = VECTOR_LOAD((CAST)(&right[j][i]));
			l_A[j] = VECTOR_BIT_AND(s_l, s_r);
			v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
			v_N = VECTOR_BIT_OR(v_N, l_A[j]);
		      }
		    
		    for(j = 0; j < 20; j++)		    
		      VECTOR_STORE((CAST)(&thisOne[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);
		  }
	      }
	      break;
	    default:
	      {
		parsimonyNumber
		  *left[32], 
		  *right[32],
		  *thisOne[32];

		assert(states <= 32);
		
		for(k = 0; k < states; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * states * rNumber) + width * k]);
		    thisOne[k]  = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    size_t j;
		    
		    INT_TYPE
		      s_r, s_l, 
		      v_N = SET_ALL_BITS_ZERO,
		      l_A[32], 
		      v_A[32];	    	 
		    
		    for(j = 0; j < states; j++)
		      {
			s_l = VECTOR_LOAD((CAST)(&left[j][i]));
			s_r = VECTOR_LOAD((CAST)(&right[j][i]));
			l_A[j] = VECTOR_BIT_AND(s_l, s_r);
			v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
			v_N = VECTOR_BIT_OR(v_N, l_A[j]);
		      }
		    
		    for(j = 0; j < states; j++)		    
		      VECTOR_STORE((CAST)(&thisOne[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);
		  }	  			
	      }
	    }	  	 
	}

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}



unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  INT_TYPE 
    allOne = SET_ALL_BITS_ONE;

  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr); 

  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = tr->partitionData[model].states,
	width = tr->partitionData[model].parsimonyLength, 
	i;

       switch(states)
	 {
	 case 2:
	   {
	     parsimonyNumber
	       *left[2],
	       *right[2];
	     
	     for(k = 0; k < 2; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
	       }     
	     
	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	                       
		 INT_TYPE      
		   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
		   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),		 
		   v_N = VECTOR_BIT_OR(l_A, l_C);
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += populationCount(v_N);
		 
		 if(sum >= bestScore)
		   return sum;		   	       
	       }
	   }
	   break;
	 case 4:
	   {
	     parsimonyNumber
	       *left[4],
	       *right[4];
      
	     for(k = 0; k < 4; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
	       }        

	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	                        
		 INT_TYPE      
		   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
		   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),
		   l_G = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[2][i])), VECTOR_LOAD((CAST)(&right[2][i]))),
		   l_T = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[3][i])), VECTOR_LOAD((CAST)(&right[3][i]))),
		   v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));     
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += populationCount(v_N);
		 
		 if(sum >= bestScore)		 
		   return sum;	        
	       }	   	 
	   }
	   break;
	 case 20:
	   {
	     parsimonyNumber
	       *left[20],
	       *right[20];
	     
	      for(k = 0; k < 20; k++)
		{
		  left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		  right[k] = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		}  
	   
	      for(i = 0; i < width; i += INTS_PER_VECTOR)
		{                	       
		  int 
		    j;
		  
		  INT_TYPE      
		    l_A,
		    v_N = SET_ALL_BITS_ZERO;     
		  
		  for(j = 0; j < 20; j++)
		    {
		      l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
		      v_N = VECTOR_BIT_OR(l_A, v_N);
		    }
		  
		  v_N = VECTOR_AND_NOT(v_N, allOne);
		  
		  sum += populationCount(v_N);	       
		  
		  if(sum >= bestScore)	    
		    return sum;		    	       
		}
	   }
	   break;
	 default:
	   {
	     parsimonyNumber
	       *left[32],  
	       *right[32]; 

	     assert(states <= 32);

	     for(k = 0; k < states; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
	       }  
	   
	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	       
		 size_t
		   j;
		 
		 INT_TYPE      
		   l_A,
		   v_N = SET_ALL_BITS_ZERO;     
		 
		 for(j = 0; j < states; j++)
		   {
		     l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
		     v_N = VECTOR_BIT_OR(l_A, v_N);
		   }
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += populationCount(v_N);	       
		 
		 if(sum >= bestScore)	      
		   return sum;		       
	       }
	   }
	 }
    }
  
  return sum;
}


#else

void newviewParsimonyIterativeFast(tree *tr)
{    
  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 

  for(index = 4; index < count; index += 4)
    {      
      unsigned int
	totalScore = 0;

      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  size_t
	    k,
	    states = tr->partitionData[model].states,
	    width = tr->partitionData[model].parsimonyLength;	 
            
	  unsigned int	
	    i;      
                 
	  switch(states)
	    {
	    case 2:       
	      {
		parsimonyNumber
		  *left[2],
		  *right[2],
		  *thisOne[2];
		
		parsimonyNumber
		   o_A,
		   o_C,
		   t_A,
		   t_C,	
		   t_N;
		
		for(k = 0; k < 2; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 2 * rNumber) + width * k]);
		    thisOne[k]  = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i++)
		  {	 	  
		    t_A = left[0][i] & right[0][i];
		    t_C = left[1][i] & right[1][i];		   

		    o_A = left[0][i] | right[0][i];
		    o_C = left[1][i] | right[1][i];
		  
		    t_N = ~(t_A | t_C);	  

		    thisOne[0][i] = t_A | (t_N & o_A);
		    thisOne[1][i] = t_C | (t_N & o_C);		   
		    
		    totalScore += populationCount(t_N);   
		  }
	      }
	      break;
	    case 4:
	      {
		parsimonyNumber
		  *left[4],
		  *right[4],
		  *thisOne[4];

		for(k = 0; k < 4; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 4 * rNumber) + width * k]);
		    thisOne[k]  = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
		  }

		parsimonyNumber
		   o_A,
		   o_C,
		   o_G,
		   o_T,
		   t_A,
		   t_C,
		   t_G,
		   t_T,	
		   t_N;

		for(i = 0; i < width; i++)
		  {	 	  
		    t_A = left[0][i] & right[0][i];
		    t_C = left[1][i] & right[1][i];
		    t_G = left[2][i] & right[2][i];	  
		    t_T = left[3][i] & right[3][i];

		    o_A = left[0][i] | right[0][i];
		    o_C = left[1][i] | right[1][i];
		    o_G = left[2][i] | right[2][i];	  
		    o_T = left[3][i] | right[3][i];

		    t_N = ~(t_A | t_C | t_G | t_T);	  

		    thisOne[0][i] = t_A | (t_N & o_A);
		    thisOne[1][i] = t_C | (t_N & o_C);
		    thisOne[2][i] = t_G | (t_N & o_G);
		    thisOne[3][i] = t_T | (t_N & o_T); 
		    
		    totalScore += populationCount(t_N);   
		  }
	      }
	      break;
	    case 20:
	      {
		parsimonyNumber
		  *left[20],
		  *right[20],
		  *thisOne[20];

		parsimonyNumber
		  o_A[20],
		  t_A[20],	  
		  t_N;

		for(k = 0; k < 20; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 20 * rNumber) + width * k]);
		    thisOne[k]  = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i++)
		  {	 	  
		    size_t k;
		    
		    t_N = 0;

		    for(k = 0; k < 20; k++)
		      {
			t_A[k] = left[k][i] & right[k][i];
			o_A[k] = left[k][i] | right[k][i];
			t_N = t_N | t_A[k];
		      }
		    
		    t_N = ~t_N;

		    for(k = 0; k < 20; k++)		      
		      thisOne[k][i] = t_A[k] | (t_N & o_A[k]);		   
		    
		    totalScore += populationCount(t_N); 
		  }
	      }
	      break;
	    default:
	      {		
		parsimonyNumber
		  *left[32],
		  *right[32],
		  *thisOne[32];
		
		parsimonyNumber
		  o_A[32],
		  t_A[32],	  
		  t_N;
		
		assert(states <= 32);
		
		for(k = 0; k < states; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * states * rNumber) + width * k]);
		    thisOne[k]  = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
		  }
		
		for(i = 0; i < width; i++)
		  {	 	  
		    t_N = 0;
		    
		    for(k = 0; k < states; k++)
		      {
			t_A[k] = left[k][i] & right[k][i];
			o_A[k] = left[k][i] | right[k][i];
			t_N = t_N | t_A[k];
		      }
		    
		    t_N = ~t_N;
		    
		    for(k = 0; k < states; k++)		      
		      thisOne[k][i] = t_A[k] | (t_N & o_A[k]);		   
		    
		    totalScore += populationCount(t_N); 
		  }
	      }			      
	    } 
	}

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}



unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr); 

  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = tr->partitionData[model].states,
	width = tr->partitionData[model].parsimonyLength, 
	i;

       switch(states)
	 {
	 case 2:
	   {
	     parsimonyNumber 
	       t_A,
	       t_C,	      
	       t_N,
	       *left[2],
	       *right[2];
	     
	     for(k = 0; k < 2; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
	       }     
	     
	     for(i = 0; i < width; i++)
	       {                	                       
		 t_A = left[0][i] & right[0][i];
		 t_C = left[1][i] & right[1][i];
		 
		  t_N = ~(t_A | t_C);

		  sum += populationCount(t_N);    
		 
		 if(sum >= bestScore)
		   return sum;		   	       
	       }
	   }
	   break;
	 case 4:
	   {
	     parsimonyNumber
	       t_A,
	       t_C,
	       t_G,
	       t_T,
	       t_N,
	       *left[4],
	       *right[4];
      
	     for(k = 0; k < 4; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
	       }        

	     for(i = 0; i < width; i++)
	       {                	                        
		  t_A = left[0][i] & right[0][i];
		  t_C = left[1][i] & right[1][i];
		  t_G = left[2][i] & right[2][i];	  
		  t_T = left[3][i] & right[3][i];

		  t_N = ~(t_A | t_C | t_G | t_T);

		  sum += populationCount(t_N);     
		 
		 if(sum >= bestScore)		 
		   return sum;	        
	       }	   	 
	   }
	   break;
	 case 20:
	   {
	     parsimonyNumber
	       t_A,
	       t_N,
	       *left[20],
	       *right[20];
	     
	      for(k = 0; k < 20; k++)
		{
		  left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		  right[k] = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		}  
	   
	      for(i = 0; i < width; i++)
		{ 
		  t_N = 0;
		  
		  for(k = 0; k < 20; k++)
		    {
		      t_A = left[k][i] & right[k][i];
		      t_N = t_N | t_A;
		    }
  	       
		  t_N = ~t_N;

		  sum += populationCount(t_N);      
		  
		  if(sum >= bestScore)	    
		    return sum;		    	       
		}
	   }
	   break;
	 default:
	   {
	     parsimonyNumber
	       t_A,
	       t_N,
	       *left[32], 
	       *right[32];  

	     assert(states <= 32);

	     for(k = 0; k < states; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
	       }  
	   
	     for(i = 0; i < width; i++)
	       {                	       
		 t_N = 0;
		  
		 for(k = 0; k < states; k++)
		   {
		     t_A = left[k][i] & right[k][i];
		     t_N = t_N | t_A;
		   }
  	       
		  t_N = ~t_N;

		  sum += populationCount(t_N);      
		  		  		 
		 if(sum >= bestScore)			  
		   return sum;			   
	       }	     	     
	   }
	 }
    }
  
  return sum;
}

#endif






static unsigned int evaluateParsimony(tree *tr, nodeptr p, fakeboolean full)
{
  volatile unsigned int result;
  nodeptr q = p->back;
  int
    *ti = tr->ti,
    counter = 4;
  
  ti[1] = p->number;
  ti[2] = q->number;

  if(full)
    {
      if(p->number > tr->mxtips)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }
  else
    {
      if(p->number > tr->mxtips && !p->x)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips && !q->x)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }

  ti[0] = counter;

  result = evaluateParsimonyIterativeFast(tr);

  return result;
}


void newviewParsimony(tree *tr, nodeptr  p)
{     
  if(p->number <= tr->mxtips)
    return;

  {
    int 
      counter = 4;     
           
    computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
    tr->ti[0] = counter;            
    
    newviewParsimonyIterativeFast(tr);      
  }
}

static void compressDNA(tree *tr, int *informative, fakeboolean saveMemory)
{
  size_t
    totalNodes,
    i,
    model;
  
  if(saveMemory)
    totalNodes = (size_t)tr->innerNodes + 1 + (size_t)tr->mxtips;
  else
    totalNodes = 2 * (size_t)tr->mxtips;

 

  for(model = 0; model < (size_t) tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = (size_t)tr->partitionData[model].states,       
	compressedEntries,
	compressedEntriesPadded,
	entries = 0, 
	lower = tr->partitionData[model].lower,
	upper = tr->partitionData[model].upper;

      parsimonyNumber 
	**compressedTips = (parsimonyNumber **)rax_malloc(states * sizeof(parsimonyNumber*)),
	*compressedValues = (parsimonyNumber *)rax_malloc(states * sizeof(parsimonyNumber));
      
      for(i = lower; i < upper; i++)    
	if(informative[i])
	  entries += (size_t)tr->cdta->aliaswgt[i];     
  
      compressedEntries = entries / PCF;

      if(entries % PCF != 0)
	compressedEntries++;

#if (defined(__SIM_SSE3) || defined(__AVX))
      if(compressedEntries % INTS_PER_VECTOR != 0)
	compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
      else
	compressedEntriesPadded = compressedEntries;
#else
      compressedEntriesPadded = compressedEntries;
#endif     

      
      tr->partitionData[model].parsVect = (parsimonyNumber *)rax_malloc((size_t)compressedEntriesPadded * states * totalNodes * sizeof(parsimonyNumber));
     
      for(i = 0; i < compressedEntriesPadded * states * totalNodes; i++)      
	tr->partitionData[model].parsVect[i] = 0;          

      for(i = 0; i < (size_t)tr->mxtips; i++)
	{
	  size_t
	    w = 0,
	    compressedIndex = 0,
	    compressedCounter = 0,
	    index = 0;

	  for(k = 0; k < states; k++)
	    {
	      compressedTips[k] = &(tr->partitionData[model].parsVect[(compressedEntriesPadded * states * (i + 1)) + (compressedEntriesPadded * k)]);
	      compressedValues[k] = 0;
	    }                
	      
	  for(index = lower; index < (size_t)upper; index++)
	    {
	      if(informative[index])
		{
		  const unsigned int 
		    *bitValue = getBitVector(tr->partitionData[model].dataType);

		  parsimonyNumber 
		    value = bitValue[tr->yVector[i + 1][index]];	  
	      
		  for(w = 0; w < (size_t)tr->cdta->aliaswgt[index]; w++)
		    {	   
		      for(k = 0; k < states; k++)
			{
			  if(value & mask32[k])
			    compressedValues[k] |= mask32[compressedCounter];
			}
		     
		      compressedCounter++;
		  
		      if(compressedCounter == PCF)
			{
			  for(k = 0; k < states; k++)
			    {
			      compressedTips[k][compressedIndex] = compressedValues[k];
			      compressedValues[k] = 0;
			    }			 
			  
			  compressedCounter = 0;
			  compressedIndex++;
			}
		    }
		}
	    }
                           
	  for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	    {	
	      for(;compressedCounter < PCF; compressedCounter++)	      
		for(k = 0; k < states; k++)
		  compressedValues[k] |= mask32[compressedCounter];		  
	  
	      for(k = 0; k < states; k++)
		{
		  compressedTips[k][compressedIndex] = compressedValues[k];
		  compressedValues[k] = 0;
		}	      	      
	      
	      compressedCounter = 0;
	    }	 	
	}               
  
      tr->partitionData[model].parsimonyLength = compressedEntriesPadded;   

      rax_free(compressedTips);
      rax_free(compressedValues);
    }
  
  tr->parsimonyScore = (unsigned int*)rax_malloc(sizeof(unsigned int) * totalNodes);  
          
  for(i = 0; i < totalNodes; i++) 
    tr->parsimonyScore[i] = 0;
}

fakeboolean tipHomogeneityChecker(tree *tr, nodeptr p, int grouping)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(tr->constraintVector[p->number] != grouping) 
	return FALSE;
      else 
	return TRUE;
    }
  else
    {   
      return  (tipHomogeneityChecker(tr, p->next->back, grouping) && tipHomogeneityChecker(tr, p->next->next->back,grouping));      
    }
}

void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];

#ifdef _BASTIEN
  for(i = 0; i < numBranches; i++)
    p->secondDerivativeValid[i] = q->secondDerivativeValid[i] = FALSE;
#endif
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;

#ifdef _BASTIEN
  for(i = 0; i < numBranches; i++)
    p->secondDerivativeValid[i] = q->secondDerivativeValid[i] = FALSE;
#endif
}


static nodeptr  removeNodeParsimony (nodeptr p, tree *tr)
{ 
  nodeptr  q, r;         

  q = p->next->back;
  r = p->next->next->back;   
    
  hookupDefault(q, r, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;
  
  return  q;
}

static void restoreTreeParsimony(tree *tr, nodeptr p, nodeptr q)
{ 
  nodeptr
    r = q->back;
  
  int counter = 4;
  
  hookupDefault(p->next,       q, tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches);
  
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
    
  newviewParsimonyIterativeFast(tr); 
}

static void restoreTreeRearrangeParsimony(tree *tr)
{    
  removeNodeParsimony(tr->removeNode, tr);  
  restoreTreeParsimony(tr, tr->removeNode, tr->insertNode);  
}

int checker(tree *tr, nodeptr p)
{
  int group = tr->constraintVector[p->number];

  if(isTip(p->number, tr->mxtips))
    {
      group = tr->constraintVector[p->number];
      return group;
    }
  else
    {
      if(group != -9) 
	return group;

      group = checker(tr, p->next->back);
      if(group != -9) 
	return group;

      group = checker(tr, p->next->next->back);
      if(group != -9) 
	return group;

      return -9;
    }
}


static void testInsertParsimony (tree *tr, nodeptr p, nodeptr q)
{ 
  unsigned int 
    mp;
 
  nodeptr  
    r = q->back;   

  fakeboolean 
    doIt = TRUE;
    
  if(tr->grouped)
    {
      int 
	rNumber = tr->constraintVector[r->number],
	qNumber = tr->constraintVector[q->number],
	pNumber = tr->constraintVector[p->number];

      doIt = FALSE;
     
      if(pNumber == -9)
	pNumber = checker(tr, p->back);
      if(pNumber == -9)
	doIt = TRUE;
      else
	{
	  if(qNumber == -9)
	    qNumber = checker(tr, q);

	  if(rNumber == -9)
	    rNumber = checker(tr, r);

	  if(pNumber == rNumber || pNumber == qNumber)
	    doIt = TRUE;       
	}
    }

  if(doIt)
    {
      insertParsimony(tr, p, q);   
  
      mp = evaluateParsimony(tr, p->next->next, FALSE);          
      
      if(mp < tr->bestParsimony)
	{
	  tr->bestParsimony = mp;
	  tr->insertNode = q;
	  tr->removeNode = p;
	}
  
      hookupDefault(q, r, tr->numBranches);
      p->next->next->back = p->next->back = (nodeptr) NULL;
    }
       
  return;
} 


static void addTraverseParsimony (tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav, fakeboolean doAll)
{        
  if (doAll || (--mintrav <= 0))               
    testInsertParsimony(tr, p, q);	                 

  if (((q->number > tr->mxtips)) && ((--maxtrav > 0) || doAll))
    {	      
      addTraverseParsimony(tr, p, q->next->back, mintrav, maxtrav, doAll);	      
      addTraverseParsimony(tr, p, q->next->next->back, mintrav, maxtrav, doAll);              	     
    }
}

static int rearrangeParsimony(tree *tr, nodeptr p, int mintrav, int maxtrav, fakeboolean doAll)  
{   
  nodeptr  
    p1, 
    p2, 
    q, 
    q1, 
    q2;
  
  int      
    mintrav2; 

  fakeboolean 
    doP = TRUE,
    doQ = TRUE;
           
  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3; 

  assert(mintrav == 1);

  if(maxtrav < mintrav)
    return 0;

  q = p->back;

  if(tr->constrained)
    {    
      if(! tipHomogeneityChecker(tr, p->back, 0))
	doP = FALSE;
	
      if(! tipHomogeneityChecker(tr, q->back, 0))
	doQ = FALSE;
		        
      if(doQ == FALSE && doP == FALSE)
	return 0;
    }  

  if((p->number > tr->mxtips) && doP) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
      if ((p1->number > tr->mxtips) || (p2->number > tr->mxtips)) 
	{	  	  
	  removeNodeParsimony(p, tr);	  	 

	  if ((p1->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, p, p1->next->back, mintrav, maxtrav, doAll);         
	      addTraverseParsimony(tr, p, p1->next->next->back, mintrav, maxtrav, doAll);          
	    }
	 
	  if ((p2->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, p, p2->next->back, mintrav, maxtrav, doAll);
	      addTraverseParsimony(tr, p, p2->next->next->back, mintrav, maxtrav, doAll);          
	    }
	    
	   
	  hookupDefault(p->next,       p1, tr->numBranches); 
	  hookupDefault(p->next->next, p2, tr->numBranches);	   	    	    

	  newviewParsimony(tr, p);
	}
    }  
       
  if ((q->number > tr->mxtips) && (maxtrav > 0) && doQ) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;

      if (
	  (
	   (q1->number > tr->mxtips) && 
	   ((q1->next->back->number > tr->mxtips) || (q1->next->next->back->number > tr->mxtips))
	   )
	  ||
	  (
	   (q2->number > tr->mxtips) && 
	   ((q2->next->back->number > tr->mxtips) || (q2->next->next->back->number > tr->mxtips))
	   )
	  )
	{	   

	  removeNodeParsimony(q, tr);
	  
	  mintrav2 = mintrav > 2 ? mintrav : 2;
	  
	  if ((q1->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, q, q1->next->back, mintrav2 , maxtrav, doAll);
	      addTraverseParsimony(tr, q, q1->next->next->back, mintrav2 , maxtrav, doAll);         
	    }
	 
	  if ((q2->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, q, q2->next->back, mintrav2 , maxtrav, doAll);
	      addTraverseParsimony(tr, q, q2->next->next->back, mintrav2 , maxtrav, doAll);          
	    }	   
	   
	  hookupDefault(q->next,       q1, tr->numBranches); 
	  hookupDefault(q->next->next, q2, tr->numBranches);
	   
	  newviewParsimony(tr, q);
	}
    }

  return 1;
} 

static void reorderNodes(tree *tr, nodeptr *np, nodeptr p, int *count)
{
  int i, found = 0;

  if((p->number <= tr->mxtips))    
    return;
  else
    {              
      for(i = tr->mxtips + 1; (i <= (tr->mxtips + tr->mxtips - 1)) && (found == 0); i++)
	{
	  if (p == np[i] || p == np[i]->next || p == np[i]->next->next)
	    {
	      if(p == np[i])			       
		tr->nodep[*count + tr->mxtips + 1] = np[i];		 		
	      else
		{
		  if(p == np[i]->next)		  
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next;		     	   
		  else		   
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next->next;		    		    
		}

	      found = 1;	      	     
	      *count = *count + 1;
	    }
	}            
     
      assert(found != 0);

      reorderNodes(tr, np, p->next->back, count);     
      reorderNodes(tr, np, p->next->next->back, count);                
    }
}

void nodeRectifier(tree *tr)
{
  nodeptr *np = (nodeptr *)rax_malloc(2 * tr->mxtips * sizeof(nodeptr));
  int i;
  int count = 0;
  
  tr->start       = tr->nodep[1];
  tr->rooted      = FALSE;

  /* TODO why is tr->rooted set to FALSE here ?*/
  
  for(i = tr->mxtips + 1; i <= (tr->mxtips + tr->mxtips - 1); i++)
    np[i] = tr->nodep[i];           
  
  reorderNodes(tr, np, tr->start->back, &count); 

 
  rax_free(np);
}

static void stepwiseAddition(tree *tr, nodeptr p, nodeptr q)
{            
  nodeptr 
    r = q->back;

  unsigned int 
    mp;
  
  int 
    counter = 4;
  
  p->next->back = q;
  q->back = p->next;

  p->next->next->back = r;
  r->back = p->next->next;
   
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
  tr->ti[1] = p->number;
  tr->ti[2] = p->back->number;
    
  mp = evaluateParsimonyIterativeFast(tr);
  
  if(mp < tr->bestParsimony)
    {    
      tr->bestParsimony = mp;
      tr->insertNode = q;     
    }
 
  q->back = r;
  r->back = q;
   
  if(q->number > tr->mxtips && tr->parsimonyScore[q->number] > 0)
    {	      
      stepwiseAddition(tr, p, q->next->back);	      
      stepwiseAddition(tr, p, q->next->next->back);              	     
    }
}

static fakeboolean isInformative(tree *tr, int dataType, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = getUndetermined(dataType);

  const unsigned int
    *bitVector = getBitVector(dataType);

  unsigned char
    nucleotide;
  
	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;
      assert(bitVector[nucleotide] > 0);	           
    }
  
  for(j = 0; j < undetermined; j++)
    {
      if(check[j] > 0)
	informativeCounter++;    
    } 
	  
  if(informativeCounter <= 1)
    return FALSE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(check[j] > 1)
	    return TRUE;
	} 
    }
     
  return FALSE;	     
}

static void determineUninformativeSites(tree *tr, int *informative)
{
  int 
    i,
    number = 0;

  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.

     Alternative definition: informative columns contain at least two types
     of nucleotides, and each nucleotide must appear at least twice in each 
     column. Kind of a pain if we intend to check for this when using, e.g.,
     amibiguous DNA encoding.
  */

  for(i = 0; i < tr->cdta->endsite; i++)
    {
      if(isInformative(tr, tr->dataVector[i], i))
	informative[i] = 1;
      else
	{
	  informative[i] = 0;
	  number++;
	}            
    }
  
 
  /* printf("Uninformative Patterns: %d\n", number); */
}

fakeboolean update(tree *tr, nodeptr p)
{       
  nodeptr  q; 
  fakeboolean smoothedPartitions[NUM_BRANCHES];
  int i;
  double   z[NUM_BRANCHES], z0[NUM_BRANCHES];
  double _deltaz;

  q = p->back;   

  for(i = 0; i < tr->numBranches; i++)
    z0[i] = q->z[i];    

  if(tr->numBranches > 1)
    makenewzGeneric(tr, p, q, z0, newzpercycle, z, TRUE);  
  else
    makenewzGeneric(tr, p, q, z0, newzpercycle, z, FALSE);
  
  for(i = 0; i < tr->numBranches; i++)    
    smoothedPartitions[i]  = tr->partitionSmoothed[i];
      
  for(i = 0; i < tr->numBranches; i++)
    {         
      if(!tr->partitionConverged[i])
	{	  
	    _deltaz = deltaz;
	    
	  if(ABS(z[i] - z0[i]) > _deltaz)  
	    {	    
	      smoothedPartitions[i] = FALSE;       
	    }	 
	  	  
	  p->z[i] = q->z[i] = z[i];	 
#ifdef _BASTIEN
	  //printf("update %f\n", tr->secondDerivative[i]);
	  p->secondDerivative[i] = q->secondDerivative[i] = tr->secondDerivative[i];
	  p->secondDerivativeValid[i] = q->secondDerivativeValid[i] = TRUE;
#endif
	}
    }
  
#ifdef _DEBUG_UPDATE
  evaluateGeneric(tr, p);

  if(tr->likelihood <= startLH)
    {
      if(fabs(tr->likelihood - startLH) > 0.01)
	{
	  printf("%f %f\n", startLH, tr->likelihood);	  
	  assert(0);             
	}
    }
#endif

  for(i = 0; i < tr->numBranches; i++)    
    tr->partitionSmoothed[i]  = smoothedPartitions[i];
  
  return TRUE;
}


double FABS(double x)
{
  /*  if(x < -1.0E-10)
      assert(0);*/
  
  /* if(x < 0.0)
     printf("%1.40f\n", x); */

  return fabs(x);
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

void nniSmooth(tree *tr, nodeptr p, int maxtimes)
{
  int
    i;

  for(i = 0; i < tr->numBranches; i++)	
    tr->partitionConverged[i] = FALSE;	

 

  while (--maxtimes >= 0) 
    {     
      

      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;
      
      

      assert(!isTip(p->number, tr->mxtips)); 	

     

      assert(!isTip(p->back->number, tr->mxtips));  
      
      update(tr, p);
     
      update(tr, p->next);
     
      update(tr, p->next->next);
      
      update(tr, p->back->next);
      
      update(tr, p->back->next->next);           
     
      if (allSmoothed(tr)) 
	break;
      
    }

  

  for(i = 0; i < tr->numBranches; i++)
    {
      tr->partitionSmoothed[i] = FALSE; 
      tr->partitionConverged[i] = FALSE;
    }

  
}

void makeParsimonyTreeFast(tree *tr, analdef *adef, fakeboolean full)
{   
  nodeptr  
    p, 
    f;    

  size_t
    model;

  int 
    i, 
    nextsp,
    *perm        = (int *)rax_malloc((size_t)(tr->mxtips + 1) * sizeof(int)),
    *informative = (int *)rax_malloc(sizeof(int) * (size_t)tr->cdta->endsite);  

  unsigned int 
    randomMP, 
    startMP;        

  /* double t; */

  determineUninformativeSites(tr, informative);     

  compressDNA(tr, informative, FALSE);

  rax_free(informative); 

  tr->ti = (int*)rax_malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  
 
  /*t = gettime();*/

  if(!full)
    {                	
      int 
	j = 0;

      unsigned char
	*nodesInTree = (unsigned char*)rax_calloc((size_t)(tr->mxtips + 1), sizeof(unsigned char));	      
	
      tr->start = findAnyTipFast(tr->start, tr->rdta->numsp);
	
      tr->bestParsimony = INT_MAX;

      evaluateParsimony(tr, tr->start->back, TRUE);
		
      assert(tr->start);
      
      checkSeed(adef);

      markNodesInTree(tr->start, tr, nodesInTree);
      markNodesInTree(tr->start->back, tr, nodesInTree);
	
      j = tr->ntips + 1;
	
      if(tr->grouped)
	{
	  for(i = 1; i <= tr->mxtips; i++)      
	    {
	      if(tr->constraintVector[i] == -1) 
		{
		  perm[j++] = i;		
		  tr->constraintVector[i] = -9;
		}
	    }
	}
      else
	{
	  if(tr->constrained)
	    { 
	      for(i = 1; i <= tr->mxtips; i++)
		tr->constraintVector[i] = 0;
	      
	      for(i = 1; i <= tr->mxtips; i++)
		{		  
		  if(nodesInTree[i] == 0) 	      
		    perm[j++] = i;
		  else
		    tr->constraintVector[i] = 1;		    
		}
	    }
	  else
	    {
	      for(i = 1; i <= tr->mxtips; i++)      
		if(nodesInTree[i] == 0) 
		  perm[j++] = i;	  
	    }
	}
	
      for(i = tr->ntips + 1; i <= tr->mxtips; i++) 
	{	     
	  int k, j;
	  
	  k =  (int)((double)(tr->mxtips + 1 - i) * randum(&adef->parsimonySeed));
	  
	  assert(i + k <= tr->mxtips);
	  j        = perm[i];
	  perm[i]     = perm[i + k];
	  perm[i + k] = j;
	}    
	
      f = tr->start;     

      rax_free(nodesInTree);
    }
  else
    {
      assert(!tr->constrained);

      makePermutation(perm, 1, tr->mxtips, adef);
      
      tr->ntips = 0;    
      
      tr->nextnode = tr->mxtips + 1;       
      
      buildSimpleTree(tr, perm[1], perm[2], perm[3]);      
      
      f = tr->start;
    }     
  
  while(tr->ntips < tr->mxtips) 
    {	
      nodeptr q;
      
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];                 
      q = tr->nodep[(tr->nextnode)++];
      p->back = q;
      q->back = p;
        
      if(tr->grouped && !full)
	{
	  int 
	    number = p->back->number;	  	 

	  tr->constraintVector[number] = -9;
	}
          
      stepwiseAddition(tr, q, f->back);      	  	 
      
      {
	nodeptr	  
	  r = tr->insertNode->back;
	
	int counter = 4;
	
	hookupDefault(q->next,       tr->insertNode, tr->numBranches);
	hookupDefault(q->next->next, r, tr->numBranches);
	
	computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips, FALSE);              
	tr->ti[0] = counter;
	
	newviewParsimonyIterativeFast(tr);	
      }
    }    
  
  //printf("ADD: %d\n", tr->bestParsimony); 
  
  nodeRectifier(tr);
  
  if(adef->stepwiseAdditionOnly == FALSE)
    {
      randomMP = tr->bestParsimony;        
      
      do
	{
	  startMP = randomMP;
	  nodeRectifier(tr);
	  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	    {
	      rearrangeParsimony(tr, tr->nodep[i], 1, 20, FALSE);
	      if(tr->bestParsimony < randomMP)
		{		
		  restoreTreeRearrangeParsimony(tr);
		  randomMP = tr->bestParsimony;
		}
	    }      		  	   
	}
      while(randomMP < startMP);
    }
  
  //printf("OPT: %d\n", tr->bestParsimony);

  
     
  rax_free(perm);  
  rax_free(tr->parsimonyScore);
  
  for(model = 0; model < (size_t) tr->NumberOfModels; model++)
    rax_free(tr->partitionData[model].parsVect);
  
  rax_free(tr->ti);
} 

void getStartingTree(tree *tr, analdef *adef)
{
  tr->likelihood = unlikely;
  
  assert(!adef->restart); 

      assert(adef->mode != PARSIMONY_ADDITION &&
	     adef->mode != MORPH_CALIBRATOR   &&
	     adef->mode != ANCESTRAL_STATES   &&
	     adef->mode != OPTIMIZE_BR_LEN_SCALER);

/*      if(adef->randomStartingTree)	  */
/*	makeRandomTree(tr, adef);       	   	 	   	  */
/*      else*/
    assert(!adef->randomStartingTree);
	makeParsimonyTreeFast(tr, adef, TRUE);	   	    	      		      	
      
/*      if(adef->startingTreeOnly)*/
/*	{*/
/*	  printStartingTree(tr, adef, TRUE);*/
/*	  exit(0);*/
/*	}*/
/*      else   	         */
/*	printStartingTree(tr, adef, FALSE);     	         */
            
      
      evaluateGenericInitrav(tr, tr->start);   

     
      
      treeEvaluate(tr, 1);        	 
        

  tr->start = tr->nodep[1];
}

fakeboolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}

static uint64_t f2(int n, int a) 
{
  return ((n - a) * (n - 1 - a) * (n - 2 - a ) / 6);
};

static uint64_t f3(int n, int b) 
{
  return ((n - b) * (n - 1 - b) / 2);
};

static uint64_t f4(int n, int c) 
{
  return (n-c);
};

static void preprocessQuartetPrefix(int numberOfTaxa, uint64_t *prefixSumF2, uint64_t *prefixSumF3, uint64_t *prefixSumF4)
{
  int 
    i,
    n = numberOfTaxa;
  
  

  prefixSumF2[0] = 1;
  prefixSumF3[0] = 1;
  prefixSumF4[0] = 1;
  
  for (i = 1; i < n - 3; ++i) 
    {
      prefixSumF2[i] = prefixSumF2[i - 1] + f2(n, i);
      prefixSumF3[i] = prefixSumF3[i - 1] + f3(n, i+1);
      prefixSumF4[i] = prefixSumF4[i - 1] + f4(n, i+2);
  }
}

static double quartetLikelihood(tree *tr, nodeptr p1, nodeptr p2, nodeptr p3, nodeptr p4, nodeptr q1, nodeptr q2, analdef *adef, fakeboolean firstQuartet)
{
  /* 
     build a quartet tree, where q1 and q2 are the inner nodes and p1, p2, p3, p4
     are the tips of the quartet where the sequence data is located.

     initially set all branch lengths to the default value.
  */

  /* 
     for the tree and node data structure used, please see one of the last chapter's of Joe 
     Felsensteins book. 
  */

  hookupDefault(q1, q2, tr->numBranches);
  
  hookupDefault(q1->next,       p1, tr->numBranches);
  hookupDefault(q1->next->next, p2, tr->numBranches);
  
  hookupDefault(q2->next,       p3, tr->numBranches);
  hookupDefault(q2->next->next, p4, tr->numBranches);
  
  /* now compute the likelihood vectors at the two inner nodes of the tree,
     here the virtual root is located between the two inner nodes q1 and q2.
  */

  newviewGeneric(tr, q1);
  newviewGeneric(tr, q2);

  
#ifdef __BLACKRIM 
  if(firstQuartet)
    {
      tr->start = q1->next->back;
  
      if(modOpt(tr, adef, TRUE, adef->likelihoodEpsilon) != 0)
        return -DBL_MAX;
    }
#endif
  
  /* call a function that is also used for NNIs that iteratively optimizes all 
     5 branch lengths in the tree.

     Note that 16 is an important tuning parameter, this integer value determines 
     how many times we visit all branches until we give up further optimizing the branch length 
     configuration.
  */

  nniSmooth(tr, q1, 16);

  /* now compute the log likelihood of the tree for the virtual root located between inner nodes q1 and q2 */
  
  /* debugging code 
     {
    double l;
  */
  
  evaluateGeneric(tr, q1->back->next->next);
  
  /* debugging code 
     
     l = tr->likelihood;

     newviewGeneric(tr, q1);
     newviewGeneric(tr, q2);
     evaluateGeneric(tr, q1);
     
   
     assert(ABS(l - tr->likelihood) < 0.00001);
     }
  */

  return (tr->likelihood);
}

static quartetResult computeAllThreeQuartets(tree *tr, nodeptr q1, nodeptr q2, int t1, int t2, int t3, int t4, analdef *adef)
{
  /* set the tip nodes to different sequences 
     with the tip indices t1, t2, t3, t4 */
	       
  nodeptr 
    p1 = tr->nodep[t1],
    p2 = tr->nodep[t2],
    p3 = tr->nodep[t3], 
    p4 = tr->nodep[t4];
  
  double 
    l; 

  quartetResult qr;
  
  /* first quartet */	    
  
  /* compute the likelihood of tree ((p1, p2), (p3, p4)) */
  
  l = quartetLikelihood(tr, p1, p2, p3, p4, q1, q2, adef, TRUE);
  
  qr.a1 = p1->number;
  qr.b1 = p2->number;
  qr.c1 = p3->number;
  qr.d1 = p4->number;
  qr.l1 = l;
  /* second quartet */	    
  
  /* compute the likelihood of tree ((p1, p3), (p2, p4)) */
  
  l = quartetLikelihood(tr, p1, p3, p2, p4, q1, q2, adef, FALSE);

  qr.a2 = p1->number;
  qr.b2 = p3->number;
  qr.c2 = p2->number;
  qr.d2 = p4->number;
  qr.l2 = l;
  /* third quartet */	    
  
  /* compute the likelihood of tree ((p1, p4), (p2, p3)) */
  
  l = quartetLikelihood(tr, p1, p4, p2, p3, q1, q2, adef, FALSE);
  
  qr.a3 = p1->number;
  qr.b3 = p4->number;
  qr.c3 = p2->number;
  qr.d3 = p3->number;
  qr.l3 = l;

    return qr;
}

static unsigned int binarySearch(uint64_t* array, uint64_t z, int n)
{
  unsigned int 
    first = 0,
    last = n-3,
    middle = (first + last) / 2, 
    lastSmaller = 0;
  
  while(first <= last)
    {
      if(array[middle] < z)
	{
	  first = middle + 1;
	  lastSmaller = middle;
	}
      else 
	{
	  if (array[middle] > z)	  
	    last = middle-1;	 
	  else 
	    { 
	      // array[middle] == z
	      lastSmaller = middle;
	      break;
	    }
	}
      
      middle = (first + last)/2;
    }

  return lastSmaller;
}

static void mapNumberToQuartet(int numberOfTaxa, uint64_t z, int *t1, int *t2, int *t3, int *t4, uint64_t *prefixSumF2, uint64_t *prefixSumF3, uint64_t *prefixSumF4)
{
  uint64_t    
    wantedT1 = z;

  *t1 = binarySearch(prefixSumF2, z, numberOfTaxa) + 1;

  uint64_t 
    foundT1 = prefixSumF2[*t1 - 1];
  
  if(wantedT1 == foundT1) 
    {
      *t2 = *t1+1;
      *t3 = *t1+2;
      *t4 = *t1+3;
      return;
    }
  
  uint64_t 
    wantedT2 = (prefixSumF3[*t1 - 1]) + (wantedT1 - foundT1);
  
  *t2 = binarySearch(prefixSumF3, wantedT2, numberOfTaxa) + 2;

  uint64_t 
    foundT2 = prefixSumF3[*t2 - 2];
  
  if(wantedT2 == foundT2) 
    {
      *t3 = *t2 + 1;
      *t4 = *t2 + 2;
      return;
    }
  
  uint64_t 
    wantedT3 = (prefixSumF4[*t2 - 2]) + (wantedT2 - foundT2);
  
  *t3 = binarySearch(prefixSumF4, wantedT3, numberOfTaxa) + 3;

  uint64_t 
    foundT3 = prefixSumF4[*t3 - 3];
  
  if (wantedT3 == foundT3) 
    {
      *t4 = *t3 + 1;
      return;
    }

  *t4 = wantedT3 - foundT3 + *t3 + 1;
}

static void sampleQuartetsWithoutReplacementA(quartetResult * result_vec, tree *tr, int numberOfTaxa, int64_t seed, uint64_t numberOfQuartets, uint64_t randomQuartets, nodeptr q1, nodeptr q2, uint64_t *prefixSumF2, uint64_t *prefixSumF3, uint64_t *prefixSumF4, analdef *adef, uint64_t actVal, uint64_t initialQuartetCounter)
{
  int64_t 
    myseed = seed;

  uint64_t    
    sampleSize = randomQuartets,
    quartetCounter = initialQuartetCounter,
    top = numberOfQuartets - sampleSize,
    s;
  
  int 
    t1,
    t2,
    t3,
    t4;

  double 
    NReal = (double)numberOfQuartets, 
    v, 
    quot; 
  
  while(sampleSize >= 2)
    {
      v = randum(&myseed);
      s = 0;
      quot = top / NReal;
    
      while (quot > v)
	{
	  s++; 
	  top--; 
	  NReal--;
	  quot = (quot * top) / NReal;
	}
    // Skip over the next s records and select the following one for the sample
      actVal += s+1;
      mapNumberToQuartet(numberOfTaxa, actVal, &t1, &t2, &t3, &t4, prefixSumF2, prefixSumF3, prefixSumF4);
      result_vec[quartetCounter] = computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, adef);
      quartetCounter++;
      
      NReal--;
      sampleSize--;
    }
  
  // Special case sampleSize == 1
  s = trunc(round(NReal) * randum(&myseed));
  // Skip over the next s records and select the following one for the sample
  actVal += s+1;
  
  mapNumberToQuartet(numberOfTaxa, actVal, &t1, &t2, &t3, &t4, prefixSumF2, prefixSumF3, prefixSumF4);
  result_vec[quartetCounter] = computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, adef);
  quartetCounter++;

  assert(quartetCounter == randomQuartets);
}

static void sampleQuartetsWithoutReplacementD(quartetResult * result_vec, tree *tr, int numberOfTaxa, int64_t seed, uint64_t numberOfQuartets, uint64_t randomQuartets, nodeptr q1, nodeptr q2, 
					      uint64_t *prefixSumF2, uint64_t *prefixSumF3, uint64_t *prefixSumF4, analdef *adef, uint64_t actVal)
{
  int64_t       
    myseed = seed;
  
  uint64_t
    sampleSize = randomQuartets,
    quartetCounter = 0,
    s,
    qu1,
    threshold,
    t,
    limit;
    
  int 
    t1,
    t2,
    t3,
    t4;
    
  double
    negalphainv = -1.0/13,
    nreal = sampleSize,
    ninv = 1.0 / nreal,
    Nreal = numberOfQuartets,
    vprime = exp(log(randum(&myseed)) * ninv),
    qu1real,
    nmin1inv,
    x,
    u, 
    negSreal,
    y1,
    y2,
    top,
    bottom;
    
  qu1 = -sampleSize + 1 + numberOfQuartets;
  qu1real = -nreal + 1.0 + Nreal;
  threshold = -negalphainv * sampleSize;

  while((sampleSize > 1) && (threshold < numberOfQuartets))
  {
    nmin1inv = 1.0 / (-1.0 + nreal);
    while(TRUE)
      {
	while (TRUE)
	  // step D2: Generate U and X
	  {
	    x = Nreal * (-vprime + 1.0);
	    s = trunc(x);
	    if (s < qu1) break;
	    vprime = exp(log(randum(&myseed)) * ninv);
	  }
	u = randum(&myseed);
	negSreal = (double) s * (-1);
	// step D3: Accept?
	y1 = exp(log(u * Nreal / qu1real) * nmin1inv);
	vprime = y1 * (-x / Nreal + 1.0) * (qu1real / (negSreal + qu1real));
	if (vprime <= 1.0) break; // Accept! test (2.8) is true
	// step D4: Accept?
	y2 = 1.0;
	top = -1.0 + Nreal;
	if(-1 + sampleSize > s)
	  {
	    bottom = -nreal + Nreal;
	    limit = -s + numberOfQuartets;
	  }
	else
	  {
	    bottom = -1.0 + negSreal + Nreal;
	    limit = qu1;
	  }
	for (t = -1 + numberOfQuartets; t >= limit; t--)
	  {
	    y2 = (y2 * top)/bottom;
	    top = -1.0 + top;
	    bottom = -1.0 + bottom;
	  }
	
	if(Nreal / (-x + Nreal) >= y1 * exp(log(y2) * nmin1inv))
	  {
	    // Accept!
	    vprime = exp(log(randum(&myseed)) * nmin1inv);
	    break;
	  }
	vprime = exp(log(randum(&myseed)) * ninv);
      }
    // Step D5: Select the (s+1)st record
    // Skip over the next s records and select the following one for the sample
    actVal += s+1;
    mapNumberToQuartet(numberOfTaxa, actVal, &t1, &t2, &t3, &t4, prefixSumF2, prefixSumF3, prefixSumF4);
    result_vec[quartetCounter] = computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, adef);
    quartetCounter++;
    
    numberOfQuartets = -s + (-1 + numberOfQuartets);
    Nreal = negSreal + (-1.0 + Nreal);
    sampleSize--;
    nreal = nreal - 1.0;
    ninv = nmin1inv;
    qu1 = qu1 - s;
    qu1real += negSreal;
    threshold += negalphainv;
  }
  if (sampleSize > 1)
    {
      // Use Method A to finish the sampling
      assert(quartetCounter == randomQuartets - sampleSize);
      sampleQuartetsWithoutReplacementA(result_vec, tr, numberOfTaxa, seed, numberOfQuartets, sampleSize, q1, q2, prefixSumF2, prefixSumF3, prefixSumF4, adef, actVal, quartetCounter);
    }
  else // Special case sampleSize == 1
    {
      s = trunc(numberOfQuartets * vprime);
      // Skip over the next s records and select the following one for the sample
      actVal += s+1;
      mapNumberToQuartet(numberOfTaxa, actVal, &t1, &t2, &t3, &t4, prefixSumF2, prefixSumF3, prefixSumF4);
      result_vec[quartetCounter] = computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, adef);
      quartetCounter++;
      assert(quartetCounter == randomQuartets);
    }
}

quartetResult * computeQuartets(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  /* some indices for generating quartets in an arbitrary way */

  int
    flavor = ALL_QUARTETS, //type of quartet calculation 
    i, 
    t1, 
    t2, 
    t3, 
    t4, 
    *groups[4],
    groupSize[4];    

  double
    fraction = 0.0, //fraction of random quartets to compute
    t;

  uint64_t
    randomQuartets = (uint64_t)(adef->multipleRuns), //number of random quartets to compute 
    quartetCounter = 0, 
    //total number of possible quartets, note that we count the following ((A,B),(C,D)), ((A,C),(B,D)), ((A,D),(B,C)) as one quartet here 
    numberOfQuartets = ((uint64_t)tr->mxtips * ((uint64_t)tr->mxtips - 1) * ((uint64_t)tr->mxtips - 2) * ((uint64_t)tr->mxtips - 3)) / 24; 
  
  /* use two inner tree nodes for building quartet trees */

  nodeptr 	
    q1 = tr->nodep[tr->mxtips + 1],
    q2 = tr->nodep[tr->mxtips + 2];

  char 
    quartetFileName[1024];
    
    
    /* my changes */
    uint64_t size = randomQuartets > numberOfQuartets ? numberOfQuartets : randomQuartets;
    quartetResult * result_vec = (quartetResult *)rax_malloc(sizeof(quartetResult) * (size + 1));
    result_vec[size].a1 = -1;
/*    memset(result_vec,0,sizeof(quartetResult) * randomQuartets);*/
    /****/

/* initialize model parameters */

  initModel(tr, rdta, cdta, adef);

  

  if(!adef->useBinaryModelFile)
    {
#ifdef _QUARTET_MPI
      //the parallel version requires a pre-computed model parameter file as input!
      assert(0);
#endif

      /* get a starting tree on which we optimize the likelihood model parameters: either reads in a tree or computes a randomized stepwise addition parsimony tree */

      getStartingTree(tr, adef);
   
      /* optimize model parameters on that comprehensive tree that can subsequently be used for evaluation of quartet likelihoods */

#ifndef __BLACKRIM //if BLACKRIM is defined, the model parameters will be optimized for each quartet individually
      if(modOpt(tr, adef, TRUE, adef->likelihoodEpsilon) != 0)
      {
        rax_free(result_vec);
        return NULL;
      }
#endif
    }else{assert(0);}


  /* figure out which flavor of quartets we want to compute */

  assert(!adef->useQuartetGrouping);
  
    {
      //if the user specified more random quartets to sample than there actually 
      //exist for the number of taxa, then fix this.
      if(randomQuartets > numberOfQuartets)
	randomQuartets = 1;
  
      if(randomQuartets == 1)   
	//change flavor if randomQuartets > possibleQuartets
	flavor = ALL_QUARTETS;
      else
	{      
	  //compute the fraction of random quartets to sample 
	  //there may be an issue here with the unit64_t <-> double cast
	  fraction = (double)randomQuartets / (double)numberOfQuartets;      
	  flavor = RANDOM_QUARTETS;
	}
    }

  /* print taxon name to taxon number correspondance table to output file */
/*#ifdef _QUARTET_MPI*/
/*  if(processID == 0)   */
/*#endif*/
/*    {*/
/*      fprintf(f, "Taxon names and indices:\n\n");*/

/*      for(i = 1; i <= tr->mxtips; i++)*/
/*	{*/
/*	  fprintf(f, "%s %d\n", tr->nameList[i], i);*/
/*	  assert(tr->nodep[i]->number == i);*/
/*	}*/
/*      */
/*      fprintf(f, "\n\n");*/
/*    }*/
  
  /* do a loop to generate some quartets to test.
     note that tip nodes/sequences in RAxML are indexed from 1,...,n
     and not from 0,...,n-1 as one might expect 
     
     tr->mxtips is the maximum number of tips in the alignment/tree
  */


  //now do the respective quartet evaluations by switching over the three distinct flavors 

#ifdef _QUARTET_MPI
  if(processID > 0)   
#endif
    {
      switch(flavor)
	{
	case ALL_QUARTETS:
	  {
	    assert(randomQuartets == 1);
	    
	    /* compute all possible quartets */
	    
	    for(t1 = 1; t1 <= tr->mxtips; t1++)
	      for(t2 = t1 + 1; t2 <= tr->mxtips; t2++)
		for(t3 = t2 + 1; t3 <= tr->mxtips; t3++)
		  for(t4 = t3 + 1; t4 <= tr->mxtips; t4++)
		    {
#ifdef _QUARTET_MPI
		      if((quartetCounter % (uint64_t)(processes - 1)) == (uint64_t)(processID - 1))
#endif
			result_vec[quartetCounter] = computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, adef);
		      quartetCounter++;
		    }
	    
	    assert(quartetCounter == numberOfQuartets);
	  }
	  break;
	case RANDOM_QUARTETS:
	  {	 
	    //code contributed by Sarah for drawing quartets without replacement :-) 
	    
	    if(adef->sampleQuartetsWithoutReplacement)
	      {
		uint64_t
		  *prefixSumF2 = (uint64_t*)rax_malloc(sizeof(uint64_t) * (size_t)(tr->mxtips - 2)),
		  *prefixSumF3 = (uint64_t*)rax_malloc(sizeof(uint64_t) * (size_t)(tr->mxtips - 2)),
		  *prefixSumF4 = (uint64_t*)rax_malloc(sizeof(uint64_t) * (size_t)(tr->mxtips - 2));

		preprocessQuartetPrefix(tr->mxtips, prefixSumF2, prefixSumF3, prefixSumF4);
		
		if(randomQuartets >= numberOfQuartets / 13)		
		  sampleQuartetsWithoutReplacementA(result_vec, tr, tr->mxtips, adef->parsimonySeed, numberOfQuartets, randomQuartets, q1, q2, prefixSumF2, prefixSumF3, prefixSumF4, adef, 0, 0);		
		else		
		  sampleQuartetsWithoutReplacementD(result_vec, tr, tr->mxtips, adef->parsimonySeed, numberOfQuartets, randomQuartets, q1, q2, prefixSumF2, prefixSumF3, prefixSumF4, adef, 0);	       

		rax_free(prefixSumF2);
		rax_free(prefixSumF3);
		rax_free(prefixSumF4);
	      }
	    else
	      {
		//endless loop ta make sure we randomly sub-sample exactly as many quartets as the user specified

		//This is not very elegant, but it works, note however, that especially when the number of 
		//random quartets to be sampled is large, that is, close to the total number of quartets 
		//some quartets may be sampled twice by pure chance. To randomly sample unique quartets 
		//using hashes or bitmaps to store which quartets have already been sampled is not memory efficient.
		//Insetad, we need to use a random number generator that can generate a unique series of random numbers 
		//and then have a function f() that maps those random numbers to the corresponding index quartet (t1, t2, t3, t4),
		//see above 
		
		do
		  {	      
		    //loop over all quartets 
		    for(t1 = 1; t1 <= tr->mxtips; t1++)
		      for(t2 = t1 + 1; t2 <= tr->mxtips; t2++)
			for(t3 = t2 + 1; t3 <= tr->mxtips; t3++)
			  for(t4 = t3 + 1; t4 <= tr->mxtips; t4++)
			    {
			      //chose a random number
			      double
				r = randum(&adef->parsimonySeed);
			      
			      //if the random number is smaller than the fraction of quartets to subsample
			      //evaluate the likelihood of the current quartet
			      if(r < fraction)
				{
#ifdef _QUARTET_MPI
				  //MPI version very simple and naive way to determine which processor 
				  //is goingt to do the likelihood calculations for this quartet
				  if((quartetCounter % (uint64_t)(processes - 1)) == (uint64_t)(processID - 1))
#endif
				    //function that computes the likelihood for all three possible unrooted trees 
				    //defined by the given quartet of taxa 
				    result_vec[quartetCounter] = computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, adef);
				  //increment quartet counter that counts how many quartets we have evaluated
				  quartetCounter++;
				}
			      
			      //exit endless loop if we have randomly sub-sampled as many quartets as the user specified
			      if(quartetCounter == randomQuartets)
				goto DONE;
			    }
		  }
		while(1);
		
	      DONE:
		assert(quartetCounter == randomQuartets);	  
	      }
	  }
	  break;
	default:
	  assert(0);
	}
    }
    
    return result_vec;
}

