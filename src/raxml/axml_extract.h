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

#ifndef _H_AXML_EXTRACT
#define _H_AXML_EXTRACT

#include <stdlib.h>
#include <inttypes.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>

#ifdef __AVX
#define BYTE_ALIGNMENT 32
#else
#define BYTE_ALIGNMENT 16
#endif

#define NUM_BRANCHES   128
#define M_GTRCAT         1
#define M_GTRGAMMA       2
#define M_BINCAT         3
#define M_BINGAMMA       4
#define M_PROTCAT        5
#define M_PROTGAMMA      6
#define M_32CAT          7
#define M_32GAMMA        8
#define M_64CAT          9
#define M_64GAMMA        10
#define GENERIC_32       6
#define GENERIC_64       7
#define TRUE            1
#define FALSE            0
#define AUTO_ML   0
#define AUTO_BIC  1
#define AUTO_AIC  2
#define AUTO_AICC 3
#define JTT          2
#define SEC_16   11
#define GTR_MULTI_STATE     2
#define NOT_DEFINED              0
#define PHYLIP 0
#define FASTA  1

#define ORDERED_MULTI_STATE 0
#define MK_MULTI_STATE      1
#define GTR_MULTI_STATE     2

#define AA_SCALE 10.0
#define AA_SCALE_PLUS_EPSILON 10.001

#define ITMAX 100

#define LEWIS_CORRECTION         1
#define FELSENSTEIN_CORRECTION   2
#define STAMATAKIS_CORRECTION    3
#define GOLDMAN_CORRECTION_1     4
#define GOLDMAN_CORRECTION_2     5
#define GOLDMAN_CORRECTION_3     6


#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))

#define ABS(x)    (((x)<0)   ?  (-(x)) : (x))
#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))
#define MAX(x,y)  (((x)>(y)) ?    (x)  : (y))
#define NINT(x)   ((int) ((x)>0 ? ((x)+0.5) : ((x)-0.5)))


#define LOG(x)  log(x)
#define EXP(x)  exp(x)

#define  TREE_EVALUATION                 0
#define  BIG_RAPID_MODE                  1
#define  CALC_BIPARTITIONS               2
#define  SPLIT_MULTI_GENE                3
#define  CHECK_ALIGNMENT                 4
#define  PER_SITE_LL                     5
#define  PARSIMONY_ADDITION              6
#define  CLASSIFY_ML                     7
#define  DISTANCE_MODE                   8
#define  GENERATE_BS                     9
#define  COMPUTE_ELW                     10
#define  BOOTSTOP_ONLY                   11
#define  COMPUTE_LHS                     12
#define  COMPUTE_BIPARTITION_CORRELATION 13
#define  COMPUTE_RF_DISTANCE             14
#define  MORPH_CALIBRATOR                15
#define  CONSENSUS_ONLY                  16
#define  FAST_SEARCH                     17        
#define  EPA_SITE_SPECIFIC_BIAS          18
#define  SH_LIKE_SUPPORTS                19
#define  CLASSIFY_MP                     20
#define  ANCESTRAL_STATES                21
#define  QUARTET_CALCULATION             22
#define  THOROUGH_OPTIMIZATION           23
#define  OPTIMIZE_BR_LEN_SCALER          24
#define  ANCESTRAL_SEQUENCE_TEST         25
#define  PLAUSIBILITY_CHECKER            26
#define  CALC_BIPARTITIONS_IC            27
#define  ROOT_TREE                       28
#define  STEAL_BRANCH_LENGTHS            29
#define  SUBTREE_EPA                     30

#define DAYHOFF      0
#define DCMUT        1
#define JTT          2
#define MTREV        3
#define WAG          4
#define RTREV        5
#define CPREV        6
#define VT           7
#define BLOSUM62     8
#define MTMAM        9
#define LG           10
#define MTART        11
#define MTZOA        12
#define PMB          13
#define HIVB         14
#define HIVW         15
#define JTTDCMUT     16
#define FLU          17 
#define STMTREV      18
#define DUMMY        19
#define DUMMY2       20
#define AUTO         21
#define LG4          22
#define LG4X         23
#define PROT_FILE    24
#define GTR_UNLINKED 25
#define GTR          26  /* GTR always needs to be the last one */

#define NUM_PROT_MODELS 27

#define ALPHA_MIN    0.02
#define ALPHA_MAX    1000.0

#define RATE_MIN     0.0001
#define RATE_MAX     1000000.0

#define INVAR_MIN    0.0001
#define INVAR_MAX    0.9999

#define TT_MIN       0.0000001
#define TT_MAX       1000000.0

#define FREQ_MIN     0.001

#define NUM_RELL_BOOTSTRAPS 1000

#define badRear         -1

#define TIP_TIP     0
#define TIP_INNER   1
#define INNER_INNER 2

#define BIPARTITIONS_ALL       0
#define GET_BIPARTITIONS_BEST  1
#define DRAW_BIPARTITIONS_BEST 2
#define BIPARTITIONS_BOOTSTOP  3
#define BIPARTITIONS_RF  4
#define GATHER_BIPARTITIONS_IC 5
#define FIND_BIPARTITIONS_IC 6
#define BIPARTITIONS_PARTIAL_TC 7

#define LG4X_RATE_MIN 0.0000001
#define LG4X_RATE_MAX 1000.0

#define SEC_6_A 0
#define SEC_6_B 1
#define SEC_6_C 2
#define SEC_6_D 3
#define SEC_6_E 4

#define SEC_7_A 5
#define SEC_7_B 6
#define SEC_7_C 7
#define SEC_7_D 8
#define SEC_7_E 9
#define SEC_7_F 10

#define SEC_16   11
#define SEC_16_A 12
#define SEC_16_B 13
#define SEC_16_C 14
#define SEC_16_D 15
#define SEC_16_E 16
#define SEC_16_F 17
#define SEC_16_I 18
#define SEC_16_J 19
#define SEC_16_K 20

#define MIN_MODEL        -1
#define BINARY_DATA      0
#define DNA_DATA         1
#define AA_DATA          2
#define SECONDARY_DATA   3
#define SECONDARY_DATA_6 4
#define SECONDARY_DATA_7 5
#define GENERIC_32       6
#define GENERIC_64       7
#define MAX_MODEL        8
#define PROT_FILE    24
#define MAX_TIP_EV     0.999999999 /* max tip vector value, sum of EVs needs to be smaller than 1.0, otherwise the numerics break down */
#define smoothings     32          /* maximum smoothing passes through tree */
#define iterations     10          /* maximum iterations of iterations per insert */
#define newzpercycle   1           /* iterations of makenewz per tree traversal */
#define nmlngth        256         /* number of characters in species name */
#define deltaz         0.00001     /* test of net branch length change in update */
#define defaultz       0.9         /* value of z assigned as starting point */
#define unlikely       -1.0E300    /* low likelihood for initialization */

#define CAT         0
#define GAMMA       1
#define GAMMA_I     2

#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))

#define ABS(x)    (((x)<0)   ?  (-(x)) : (x))
#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))
#define MAX(x,y)  (((x)>(y)) ?    (x)  : (y))
#define NINT(x)   ((int) ((x)>0 ? ((x)+0.5) : ((x)-0.5)))


#define LOG(x)  log(x)
#define EXP(x)  exp(x)

#define ALL_QUARTETS 0
#define RANDOM_QUARTETS 1
#define GROUPED_QUARTETS 2

#define PCF 32

#define zmin       1.0E-15  /* max branch prop. to -log(zmin) (= 34) */
#define zmax (1.0 - 1.0E-6) /* min branch prop. to 1.0-zmax (= 1.0E-6) */

#define twotothe256  \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0
                                                     /*  2**256 (exactly)  */

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood

#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))

#ifdef __AVX

#ifdef __SIM_SSE3

#define _SSE3_WAS_DEFINED

#undef __SIM_SSE3

#endif

#endif


#ifdef __SIM_SSE3

#include <xmmintrin.h>
#include <pmmintrin.h>
  
#endif

#ifdef __AVX

#include <xmmintrin.h>
#include <immintrin.h>

#endif

#ifdef __SIM_SSE3

#define INTS_PER_VECTOR 4
#define INT_TYPE __m128i
#define CAST __m128i*
#define SET_ALL_BITS_ONE _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm_load_si128
#define VECTOR_BIT_AND _mm_and_si128
#define VECTOR_BIT_OR  _mm_or_si128
#define VECTOR_STORE  _mm_store_si128
#define VECTOR_AND_NOT _mm_andnot_si128

#endif

#ifdef __AVX

#define INTS_PER_VECTOR 8
#define INT_TYPE __m256d
#define CAST double*
#define SET_ALL_BITS_ONE (__m256d)_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO (__m256d)_mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm256_load_pd
#define VECTOR_BIT_AND _mm256_and_pd
#define VECTOR_BIT_OR  _mm256_or_pd
#define VECTOR_STORE  _mm256_store_pd
#define VECTOR_AND_NOT _mm256_andnot_pd

#endif

typedef  int fakeboolean;

typedef struct {
  double lh;
  int tree;
  double weight;
} elw;

struct ent
{
  unsigned int *bitVector;
  unsigned int *treeVector;
  unsigned int amountTips;
  int *supportVector;
  unsigned int bipNumber;
  unsigned int bipNumber2;
  unsigned int supportFromTreeset[2]; 
  
  //added by Kassian for TC/IC correction on partial gene trees
  unsigned int *taxonMask;
  unsigned int   bLink;
  double         adjustedSupport;
  double         tempSupport;
  int            tempSupportFrom;
  unsigned int   coveredNumber;
  fakeboolean        covered;
  //Kassian modif end 

  struct ent *next;
};

typedef struct ent entry;

typedef unsigned int hashNumberType;

typedef unsigned int parsimonyNumber;

typedef struct
{
  hashNumberType tableSize;
  entry **table;
  hashNumberType entryCount;
}
  hashtable;


struct stringEnt
{
  int nodeNumber;
  char *word;
  struct stringEnt *next;
};

typedef struct stringEnt stringEntry;
 
typedef struct
{
  hashNumberType tableSize;
  stringEntry **table;
}
  stringHashtable;




typedef struct ratec
{
  double accumulatedSiteLikelihood;
  double rate;
}
  rateCategorize;


typedef struct
{
  int tipCase;
#ifdef _HET
  fakeboolean parentIsTip;
#endif
#ifdef _BASTIEN
  double secondDerivativeQ[NUM_BRANCHES];
  double secondDerivativeR[NUM_BRANCHES];
  double secondDerivativeP[NUM_BRANCHES];
#endif
  int pNumber;
  int qNumber;
  int rNumber;
  double qz[NUM_BRANCHES];
  double rz[NUM_BRANCHES];
} traversalInfo;

typedef struct
{
  traversalInfo *ti;
  int count;
} traversalData;


struct noderec;

typedef struct epBrData
{
  int    *countThem;
  int    *executeThem;
  unsigned int *parsimonyScore;
  double *branches;
  double *distalBranches; 
  double *likelihoods;
  double originalBranchLength;
  char branchLabel[64];
  int leftNodeNumber;
  int rightNodeNumber;
  int *leftScaling;
  int *rightScaling;
  double branchLengths[NUM_BRANCHES];
  double *left;
  double *right;
  int branchNumber;
  int jointLabel;
} epaBranchData;

typedef struct
{
  epaBranchData *epa;
  unsigned int *vector; 
  int support;
  int *supports;
  double ic;
  double icAll;
  struct noderec *oP;
  struct noderec *oQ;
} branchInfo;

typedef struct
{
  fakeboolean valid;
  int partitions;
  int *partitionList;
}
  linkageData;

typedef struct
{
  int entries;
  linkageData* ld;
}
  linkageList;


typedef  struct noderec
{  
  branchInfo      *bInf;
  double           z[NUM_BRANCHES];
#ifdef _BASTIEN
  double           secondDerivative[NUM_BRANCHES];
  fakeboolean          secondDerivativeValid[NUM_BRANCHES];
#endif
  struct noderec  *next;
  struct noderec  *back;
  hashNumberType   hash;
  int              support;
  int              number;
  char             x;
}
  node, *nodeptr;

typedef struct
  {
    double lh;
    double pendantBranch;
    double distalBranch;    
    int number;
  }
  info;

typedef struct bInf {
  double likelihood;
  nodeptr node;
} bestInfo;

typedef struct iL {
  bestInfo *list;
  int n;
  int valid;
} infoList;


typedef  struct
{
  int              numsp;
  int              sites;
  unsigned char             **y;
  unsigned char             *y0;
  unsigned char             *yBUF;
  int              *wgt;
} rawdata;

typedef  struct {
  int             *alias;       /* site representing a pattern */
  int             *aliaswgt;    /* weight by pattern */
  int             *rateCategory;
  int              endsite;     /* # of sequence patterns */
  double          *patrat;      /* rates per pattern */
  double          *patratStored; 
} cruncheddata;


typedef struct {
  int     states;
  int     maxTipStates;
  size_t    lower;
  size_t     upper;
  size_t     width;
  int     dataType;
  int     protModels;
  int     autoProtModels;
  fakeboolean usePredefinedProtFreqs;
  int     mxtips;
  fakeboolean optimizeBaseFrequencies;
  int     numberOfCategories;
  int             **expVector;
  double          **xVector;
  size_t             *xSpaceVector;
  size_t             *expSpaceVector;
 
  unsigned char            **yVector;
 

  //asc bias
  fakeboolean ascBias;  
  int     ascOffset;
  int     *ascExpVector;
  double  *ascSumBuffer;
  double  *ascVector;
  double ascScaler[64];
  //asc bias end


  char   *partitionName;
  char   proteinSubstitutionFileName[2048];
  char   ascFileName[2048];
  double externalAAMatrix[420];

  double *sumBuffer;
   double *gammaRates;

  double *EIGN;
  double *EV;
  double *EI;  

  double *left;
  double *right;

  double    *invariableFrequencies;
  double    invariableWeight;

#ifdef _HET
  /* heterotachy */
 
  double *EIGN_TIP;
  double *EV_TIP;
  double *EI_TIP;  
  double *tipVector_TIP;
 
  double *substRates_TIP;
#endif


  /* LG4 */

  double *EIGN_LG4[4];
  double *rawEIGN_LG4[4];
  double *EV_LG4[4];
  double *EI_LG4[4];  

 

  double *frequencies_LG4[4];
  double *tipVector_LG4[4];
  double *substRates_LG4[4];
  
  /* LG4X */

  double weights[4];
  double weightExponents[4];

  double weightsBuffer[4];
  double weightExponentsBuffer[4];

  /* LG4 */

  double *frequencies;
  double *freqExponents;
  double *tipVector;
 
  double *substRates;
  double *perSiteLL;
  
  double *perSiteRates;
  double *unscaled_perSiteRates;

  unsigned int    *globalScaler;
 
  int    *wgt;
  int    *invariant;
  int    *rateCategory;
  int    *symmetryVector;
  int    *frequencyGrouping;
  fakeboolean nonGTR;
  double alpha;
  double propInvariant;

  int gapVectorLength;
  unsigned int *gapVector;
  double *gapColumn;

  size_t initialGapVectorSize;

  size_t parsimonyLength;
  parsimonyNumber *parsVect; 

  double brLenScaler;
  //andre opt
  unsigned int *presenceMap;

} pInfo;



typedef struct 
{
  int left;
  int right;
  double likelihood;
} lhEntry;


typedef struct 
{
  int count;
  int size;
  lhEntry *entries;
} lhList;



typedef struct idlist
{
  int value; 
  struct idlist *next; 
} IdList;   

typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;


/***************************************************************/

typedef struct
{
  double z[NUM_BRANCHES];
  nodeptr p, q;
  int cp, cq;
}
  connectRELL, *connptrRELL;

typedef  struct
{
  connectRELL     *connect; 
  int             start;
  double          likelihood;
}
  topolRELL;


typedef  struct
{
  int max;
  topolRELL **t;
}
  topolRELL_LIST;


/* simple tree structure */

typedef struct
{
  nodeptr p, q;
  //int cp, cq;
}
  connectTree, *connptrTree;

typedef  struct
{
  connectTree     *connectt; 
  int             start;
  double          likelihood;
}
  topolTree;


typedef  struct
{
  int max;
  topolTree **t;
}
  treeList;

  typedef struct conntyp {
      double           z[NUM_BRANCHES];           /* branch length */
      node            *p, *q;       /* parent and child sectors */
      void            *valptr;      /* pointer to value of subtree */
      int              descend;     /* pointer to first connect of child */
      int              sibling;     /* next connect from same parent */
  } connectt, *connptr;

typedef  struct {
    double           likelihood;
  int              initialTreeNumber;
    connectt         *links;       /* pointer to first connect (start) */
    node            *start;
    int              nextlink;    /* index of next available connect */
                                  /* tr->start = tpl->links->p */
    int              ntips;
    int              nextnode;
    int              scrNum;      /* position in sorted list of scores */
    int              tplNum;      /* position in sorted list of trees */

    } topol;



typedef struct {
    double           best;        /* highest score saved */
    double           worst;       /* lowest score saved */
    topol           *start;       /* starting tree for optimization */
    topol          **byScore;
    topol          **byTopol;
    int              nkeep;       /* maximum topologies to save */
    int              nvalid;      /* number of topologies saved */
    int              ninit;       /* number of topologies initialized */
    int              numtrees;    /* number of alternatives tested */
    fakeboolean          improved;
    } bestlist;

typedef  struct  {
  fakeboolean optimizeAllTrees;
  fakeboolean saveMemory;
  int    *resample;
  treeList *rellTrees;
  int numberOfBranches;
  int    numberOfTipsForInsertion;
  int    *readPartition;
  fakeboolean perPartitionEPA;
  int    *inserts;
  int    branchCounter;
  int *ti; 
  int numberOfTrees; 
  stringHashtable  *nameHash;
  pInfo            *partitionData;
  pInfo            *initialPartitionData;
  pInfo            *extendedPartitionData;
  int              *dataVector;
  int              *initialDataVector;
  int              *extendedDataVector;
  int              *patternPosition;
  int              *columnPosition;
  char             *secondaryStructureInput;
  fakeboolean          *executeModel;
  double           *perPartitionLH;
  double           *storedPerPartitionLH;
  traversalData td[1];
  unsigned int *parsimonyScore;
  int              maxCategories;
  double           *sumBuffer;
  double           *perSiteLL;  
  double           coreLZ[NUM_BRANCHES];
  int              modelNumber;
  int              multiBranch;
  int              numBranches;
  int              maxNodes;
  int              bootStopCriterion;
  int              consensusType;
  int              consensusUserThreshold;
  double           wcThreshold;
  double          *storedBrLens;
  fakeboolean         useBrLenScaler;
  fakeboolean          useFastScaling;
  branchInfo	   *bInf;
  int              multiStateModel;
  size_t innerNodes;
  fakeboolean curvatOK[NUM_BRANCHES];
  /* the stuff below is shared among DNA and AA, span does
     not change depending on datatype */
  double           *invariants;
  /* model stuff end */
  unsigned char             **yVector;
  int              secondaryStructureModel;
  int              discreteRateCategories;
  int              originalCrunchedLength;
  int              fullSites;
  int              *originalModel;
  int              *originalDataVector;
  int              *originalWeights;
  int              *secondaryStructurePairs;
  double            *partitionContributions; 
  int               ascertainmentCorrectionType;
  int               autoProteinSelectionType;
  unsigned int      numberOfEPAEntries;
  double            accumulatedEPACutoff;
  fakeboolean           useAccumulatedEPACutoff;
  double            probThresholdEPA;
#ifdef _BASTIEN
  double           secondDerivative[NUM_BRANCHES];
  fakeboolean          doBastienStuff;
#endif
  double            lhCutoff;
  double            lhAVG;
  uint64_t          lhDEC;
  uint64_t          itCount;
  int               numberOfInvariableColumns;
  int               weightOfInvariableColumns;
  int               rateHetModel;
  double           startLH;
  double           endLH;
  double           likelihood;
  double          *likelihoods;
  int             *invariant;
  node           **nodep;
  node            *start;
  int              mxtips;
  int              mxtipsVector[NUM_BRANCHES];
  int              *model;
  int              *constraintVector;
  int              numberOfSecondaryColumns;
  fakeboolean          searchConvergenceCriterion;
  int              branchLabelCounter;
  int              ntips;
  int              binaryFile_ntips;
  int              nextnode;
  int              NumberOfModels;
  int              parsimonyLength;
  int              checkPointCounter;
  int              treeID;
  int              numberOfOutgroups;
  int             *outgroupNums;
  char           **outgroups;
  fakeboolean          useEpaHeuristics;
  double           fastEPAthreshold;
  fakeboolean          bigCutoff;
  fakeboolean          partitionSmoothed[NUM_BRANCHES];
  fakeboolean          partitionConverged[NUM_BRANCHES];
  fakeboolean          rooted;
  fakeboolean          grouped;
  fakeboolean          constrained;
  fakeboolean          doCutoff;
  fakeboolean          catOnly;
  rawdata         *rdta;
  cruncheddata    *cdta;
  char **nameList;
  char *tree_string;
  size_t treeStringLength;
  unsigned int bestParsimony;
  double bestOfNode;
  nodeptr removeNode;
  nodeptr insertNode;
  double zqr[NUM_BRANCHES];
  double currentZQR[NUM_BRANCHES];
  double currentLZR[NUM_BRANCHES];
  double currentLZQ[NUM_BRANCHES];
  double currentLZS[NUM_BRANCHES];
  double currentLZI[NUM_BRANCHES];
  double lzs[NUM_BRANCHES];
  double lzq[NUM_BRANCHES];
  double lzr[NUM_BRANCHES];
  double lzi[NUM_BRANCHES];
  int mr_thresh;
  fakeboolean wasRooted;
  nodeptr leftRootNode;
  nodeptr rightRootNode;
  int rootLabel;
  fakeboolean useGammaMedian;
  fakeboolean noRateHet;
  fakeboolean corrected_IC_Score;
  fakeboolean useK80;
  fakeboolean useHKY85;
  fakeboolean useJC69;  
  fakeboolean doSubtreeEPA;
} tree;

typedef  struct {
  int              categories;
  int              model;
  int              bestTrav;
  int              max_rearrange;
  int              stepwidth;
  int              initial;
  fakeboolean          initialSet;
  int              mode;
  int64_t             boot;
  int64_t             rapidBoot;
  fakeboolean          bootstrapBranchLengths;
  fakeboolean          restart;
  fakeboolean          useWeightFile;
  fakeboolean          useMultipleModel;
  fakeboolean          constraint;
  fakeboolean          grouping;
  fakeboolean          randomStartingTree;
  fakeboolean          useInvariant;
  int            protEmpiricalFreqs;
  int            proteinMatrix;
  int            checkpoints;
  int            startingTreeOnly;
  int            multipleRuns;
  int64_t           parsimonySeed;
  int64_t           constraintSeed;
  fakeboolean        perGeneBranchLengths;
  fakeboolean        likelihoodTest;
  fakeboolean        outgroup;
  fakeboolean        permuteTreeoptimize;
  fakeboolean        allInOne;
  fakeboolean        generateBS;
  fakeboolean        bootStopping;
  fakeboolean        useExcludeFile;
  fakeboolean        userProteinModel;
  fakeboolean        computeELW;
  fakeboolean        computeDistance;
  fakeboolean        compressPatterns;
  fakeboolean        useSecondaryStructure; 
  double         likelihoodEpsilon;
  double         gapyness;
  int            similarityFilterMode;
  double        externalAAMatrix[420];
  fakeboolean       useFloat;
  fakeboolean       readTaxaOnly;
  fakeboolean       veryFast;
  fakeboolean       useBinaryModelFile;
  fakeboolean       leaveDropMode;
  int           slidingWindowSize;
  fakeboolean       checkForUndeterminedSequences;
  fakeboolean       useQuartetGrouping;
  int           alignmentFileType;
  fakeboolean       calculateIC;
  fakeboolean       verboseIC;
  fakeboolean       stepwiseAdditionOnly;
  fakeboolean       optimizeBaseFrequencies;
  fakeboolean       ascertainmentBias;
  fakeboolean       rellBootstrap;
  fakeboolean       mesquite;
  fakeboolean       silent;
  fakeboolean       noSequenceCheck;
  fakeboolean       useBFGS;
  fakeboolean       setThreadAffinity;
  int           bootstopPermutations;
  int           fcThreshold; 
  fakeboolean       sampleQuartetsWithoutReplacement;
} analdef;

typedef struct 
{
  int leftLength;
  int rightLength;
  int eignLength;
  int evLength;
  int eiLength;
  int substRatesLength;
  int frequenciesLength;
  int tipVectorLength;
  int symmetryVectorLength;
  int frequencyGroupingLength;
  fakeboolean nonGTR;
  int undetermined;
  const char *inverseMeaning;
  int states;
  fakeboolean smoothFrequencies;
  const unsigned  int *bitVector;
} partitionLengths;

typedef struct 
{
  int a1;
  int b1;
  int c1; 
  int d1;
  
  int a2;
  int b2;
  int c2; 
  int d2;

  int a3;
  int b3;
  int c3; 
  int d3;

  double l1;
  double l2;
  double l3;
} quartetResult;

// functions 

void errorExit(int e);

extern void *rax_malloc( size_t size );
extern void *rax_realloc(void *p, size_t size, fakeboolean needsMemoryAlignment);
extern void rax_free(void *p);
extern void *rax_calloc(size_t n, size_t size);
extern double evaluateGenericInitrav (tree *tr, nodeptr p);
extern void initModel ( tree *tr, rawdata *rdta, cruncheddata *cdta, analdef *adef );
extern fakeboolean treeEvaluate ( tree *tr, double smoothFactor );
extern int modOpt ( tree *tr, analdef *adef , fakeboolean resetModel, double likelihoodEpsilon);
extern void newviewGeneric (tree *tr, nodeptr p);
extern void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, fakeboolean mask);
extern double evaluateGeneric (tree *tr, nodeptr p);

void newviewParsimony(tree *tr, nodeptr  p);
double FABS(double x);
void hookup ( nodeptr p, nodeptr q, double *z, int numBranches);
void hookupDefault ( nodeptr p, nodeptr q, int numBranches);
extern double PointChi2 ( double prob, double v );
extern void scaleLG4X_EIGN(tree *tr, int model);
extern unsigned int precomputed16_bitcount(unsigned int n);

extern void computeBootStopOnly(tree *tr, char *bootStrapFileName, analdef *adef);
extern fakeboolean bootStop(tree *tr, hashtable *h, int numberOfTrees, double *pearsonAverage, unsigned int **bitVectors, int treeVectorLength, unsigned int vectorLength, analdef *adef);
extern void computeConsensusOnly(tree *tr, char* treeSetFileName, analdef *adef, fakeboolean computeIC);
extern double evaluatePartialGeneric (tree *, int i, double ki, int _model);
extern double evaluateGeneric (tree *tr, nodeptr p);
extern void newviewGeneric (tree *tr, nodeptr p);
extern void newviewGenericMulti (tree *tr, nodeptr p, int model);
extern void newviewGenericMasked (tree *tr, nodeptr p);
extern void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, fakeboolean mask);
extern void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2);
extern double evaluatePartitionGeneric (tree *tr, nodeptr p, int model);
extern void newviewPartitionGeneric (tree *tr, nodeptr p, int model);
extern double evaluateGenericVector (tree *tr, nodeptr p);
extern void categorizeGeneric (tree *tr, nodeptr p);
extern double makenewzPartitionGeneric(tree *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model);
extern fakeboolean isTip(int number, int maxTips);
extern void computeTraversalInfo(tree *tr, nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches);

extern fakeboolean initrav ( tree *tr, nodeptr p );
extern void initravPartition ( tree *tr, nodeptr p, int model );
extern fakeboolean update ( tree *tr, nodeptr p );
extern fakeboolean smooth ( tree *tr, nodeptr p );
extern fakeboolean smoothTree ( tree *tr, int maxtimes );
extern fakeboolean localSmooth ( tree *tr, nodeptr p, int maxtimes );

extern void scaleLG4X_EIGN(tree *tr, int model);
fakeboolean getSmoothFreqs(int dataType);
void initAdef(analdef *);
void set_custom_options(analdef * , tree *, int, fakeboolean);
void getinput_begin(analdef *, rawdata *, cruncheddata *, tree *);
void getinput_end(analdef *, rawdata *, cruncheddata *, tree *);
unsigned char getUndetermined(int dataType);
const unsigned int *getBitVector(int dataType);
extern void updatePerSiteRates(tree *tr, fakeboolean scaleRates);

void newviewIterative (tree *tr);
extern void determineFullTraversal(nodeptr p, tree *tr);
extern void ascertainmentBiasSequence(unsigned char tip[32], int numStates, int dataType, int nodeNumber);
void getxnode (nodeptr p);
quartetResult * computeQuartets(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta);
double randum (int64_t  *seed);

extern void makeGammaCats(int rateHetModel, double alpha, double *gammaRates, int K, fakeboolean useMedian, double propInvariant);
extern void initReversibleGTR(tree *tr, int model);
extern void initRateMatrix(tree *tr);
extern void onlyInitrav(tree *tr, nodeptr p);

void getyspace (rawdata *rdta);
void addword(char *s, stringHashtable *h, int nodeNumber);
stringHashtable *initStringHashTable(hashNumberType n);

#ifdef __AVX

void newviewGTRGAMMAPROT_AVX_LG4(int tipCase,
				 double *x1, double *x2, double *x3, double *extEV[4], double *tipVector[4],
				 int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
				 double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling);

void newviewGTRCAT_AVX_GAPPED_SAVE(int tipCase,  double *EV,  int *cptr,
				   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling,
				   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

void newviewGTRCATPROT_AVX_GAPPED_SAVE(int tipCase, double *extEV,
				       int *cptr,
				       double *x1, double *x2, double *x3, double *tipVector,
				       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling,
				       unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				       double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats);

void  newviewGTRGAMMA_AVX_GAPPED_SAVE(int tipCase,
				      double *x1_start, double *x2_start, double *x3_start,
				      double *extEV, double *tipVector,
				      int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				      const int n, double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling,
				      unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
				      double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn,
				      const unsigned int x1_presenceMap,
				      const unsigned int x2_presenceMap
				      );

void newviewGTRGAMMAPROT_AVX_GAPPED_SAVE(int tipCase,
					 double *x1_start, double *x2_start, double *x3_start, double *extEV, double *tipVector,
					 int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
					 double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling,
					 unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
					 double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn); 

void newviewGTRCAT_AVX(int tipCase,  double *EV,  int *cptr,
    double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    int n,  double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling);


void newviewGenericCATPROT_AVX(int tipCase, double *extEV,
    int *cptr,
    double *x1, double *x2, double *x3, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    int n, double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling);


void newviewGTRGAMMA_AVX(int tipCase,
			 double *x1_start, double *x2_start, double *x3_start,
			 double *EV, double *tipVector,
			 int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			 const int n, double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling,
			 const unsigned int x1_presenceMap,
			 const unsigned int x2_presenceMap
			 );

void newviewGTRGAMMAPROT_AVX(int tipCase,
			     double *x1, double *x2, double *x3, double *extEV, double *tipVector,
			     int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
			     double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling);

void newviewGTRCATPROT_AVX(int tipCase, double *extEV,
			       int *cptr,
			       double *x1, double *x2, double *x3, double *tipVector,
			       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			   int n, double *left, double *right, int *wgt, int *scalerIncrement, const fakeboolean useFastScaling);

#endif

fakeboolean isGap(unsigned int *x, int pos);
fakeboolean noGap(unsigned int *x, int pos);

#endif
