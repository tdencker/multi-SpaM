/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#ifndef RAXMLWRAPPER_HPP_
#define RAXMLWRAPPER_HPP_

extern "C"
{
#include "raxml/axml_extract.h"
}
#include "pseudoalignment.hpp"
#include <vector>
#include <iostream>

inline int idx2code(PseudoAlignment & pa, int idx)
{
    return pa.getSeq(idx - 1);
}

/**
* @brief RAxML does not free most of its large structs, so each run of RAxML would lead to memory leaks. This function frees
* all the things that RAxML does not so that it can be used multiple times.
**/

void freeAllTheThings(analdef * adef, tree * tr, cruncheddata * cdta, rawdata * rdta)
{
    // setuptree
    rax_free(tr->partitionContributions);
    rax_free(tr->perPartitionLH);
    rax_free(tr->storedPerPartitionLH);
    rax_free(tr->yVector);
    rax_free(tr->likelihoods);
    rax_free(tr->tree_string);
    rax_free(tr->td[0].ti);
    rax_free(tr->constraintVector);
    rax_free(tr->nameList);
    rax_free(tr->nodep[1]);
    rax_free(tr->nodep);
    
    // sitecombcrunch
    rax_free(tr->patternPosition);
    rax_free(tr->columnPosition);
    
    // makevalues
    rax_free(tr->invariant);
    rax_free(tr->originalDataVector);
    rax_free(tr->originalModel);
    rax_free(tr->originalWeights);
    rax_free(rdta->y0);
    rax_free(rdta->yBUF);
    
    // getinput_begin
    rax_free(rdta->wgt);
    rax_free(cdta->alias);
    rax_free(cdta->aliaswgt);
    rax_free(cdta->rateCategory);
    rax_free(tr->model);
    rax_free(tr->initialDataVector);
    rax_free(tr->extendedDataVector);    
    rax_free(cdta->patrat);
    rax_free(cdta->patratStored);
    rax_free(tr->executeModel);
    
    //allocnodex
    rax_free(tr->sumBuffer);
    rax_free(tr->perSiteLL);
    
    //allocPartitions + allocnodex
    for(int i = 0; i < tr->NumberOfModels; i++)
    {
        assert(!tr->partitionData[i].ascBias);
        rax_free(tr->partitionData[i].gapVector);
        rax_free(tr->partitionData[i].gapColumn	);
        if(tr->useFastScaling)	
            rax_free(tr->partitionData[i].globalScaler); 	         
        rax_free(tr->partitionData[i].left);
        rax_free(tr->partitionData[i].right);    
        rax_free(tr->partitionData[i].EIGN);
        rax_free(tr->partitionData[i].EV);
        rax_free(tr->partitionData[i].EI);
        rax_free(tr->partitionData[i].substRates);
        rax_free(tr->partitionData[i].frequencies);
        rax_free(tr->partitionData[i].freqExponents);
        rax_free(tr->partitionData[i].tipVector);
        rax_free(tr->partitionData[i].invariableFrequencies);    
        rax_free(tr->partitionData[i].symmetryVector);
        rax_free(tr->partitionData[i].frequencyGrouping);
        rax_free(tr->partitionData[i].perSiteRates);
        rax_free(tr->partitionData[i].unscaled_perSiteRates);
        rax_free(tr->partitionData[i].gammaRates);
        rax_free(tr->partitionData[i].yVector);
        for(unsigned j = 0; j < tr->innerNodes; ++j)
        {
            if(tr->partitionData[i].xVector[j] != nullptr)
                rax_free(tr->partitionData[i].xVector[j]); 
            if(tr->partitionData[i].expVector[j] != nullptr)
                rax_free(tr->partitionData[i].expVector[j]);
        }
        rax_free(tr->partitionData[i].xVector);   
        rax_free(tr->partitionData[i].xSpaceVector);
        rax_free(tr->partitionData[i].expVector);
        rax_free(tr->partitionData[i].expSpaceVector);
        rax_free(tr->partitionData[i].presenceMap);
    }
    
    rax_free(tr->initialPartitionData[0].partitionName);
    rax_free(tr->initialPartitionData);
        
    rax_free(adef);
    rax_free(tr);
    rax_free(cdta);
    rax_free(rdta);
}

/**
* @brief This function runs RAxML to find the best quartet trees. It creates multiple structs required by RAxML and sets 
* all the necessary options equal to the "-f q -m GTRGAMMA -p 12345" flags. It checks for multiple errors such as not having
* all nucleotides in the quartet block or bad model errors. In the latter case, it restarts RAxML with the GTRCAT model which
* requires a little more runtime. The results are written to the ostream out in the format 1,2|3,4 where 1-4 are indices of
* the sequences used by this program.
* @return the number of quartets written to file
**/

uint64_t computeAndPrintBestQuartets(PseudoAlignment & pa, std::ostream & out = std::cout, bool gamma_model = true)
{
    analdef * adef          = (analdef *) rax_malloc(sizeof(analdef));
    tree * tr               = (tree *) rax_malloc(sizeof(tree));
    rawdata * rdta          = (rawdata *) rax_malloc(sizeof(rawdata));
    cruncheddata * cdta     = (cruncheddata *) rax_malloc(sizeof(cruncheddata));
    rdta->numsp = pa.size();
    rdta->sites = pa.getLength();
    initAdef(adef);
    set_custom_options(adef, tr, 1000, gamma_model);
    getinput_begin(adef, rdta, cdta, tr);
    
    int meaningDNA[256];
    meaningDNA[0] =  1;
    meaningDNA[1] =  2;
    meaningDNA[2] =  4;
    meaningDNA[3] =  8;

    std::vector<int> count_vec(256,0);
    constexpr int target_count = 4;
    
    for(unsigned i = 0; i < pa.size(); ++i)
    {
        if(pa.getRevComp(i))
        {
            std::reverse_copy(pa.getPos(i) - pa.getLength(), pa.getPos(i), &rdta->y[i+1][1]);
        }
        else
        {
            std::copy(pa.getPos(i), pa.getPos(i) + pa.getLength(), &rdta->y[i+1][1]);
        }
        for(auto ptr = &rdta->y[i+1][1]; ptr <= &rdta->y[i+1][pa.getLength()]; ++ptr)
        {
            if (*ptr == (std::numeric_limits<char>::max)())
            {
                mspamstats::bad_character_error++;
                return 0; // illegal character at the don't care positions
            }
            *ptr = meaningDNA[*ptr];
            count_vec[*ptr]++;
        }
    }
    if(std::count_if(count_vec.begin(), count_vec.end(), [&](int x){return x > 0;}) != target_count)
    {
        mspamstats::not_all_nucleotides_error++;
        return 0;
    }
    getinput_end(adef, rdta, cdta, tr);
    
    quartetResult * result_vec = computeQuartets(tr, adef, rdta, cdta);
    if(result_vec == nullptr)
    {
        mspamstats::bad_model_errors++;
        freeAllTheThings(adef, tr, cdta, rdta);
        return computeAndPrintBestQuartets(pa, out, false);
    }
    
    uint64_t num_quartets = 0;
    for(size_t i = 0; result_vec[i].a1 != -1; ++i)
    {
        num_quartets++;
        int bestQuartet = -1;
        double bestVal = 0;
        if(result_vec[i].l1 > result_vec[i].l2 && result_vec[i].l1 > result_vec[i].l3)
        {
            bestVal = result_vec[i].l1;
            bestQuartet = 1;
        }
        if(result_vec[i].l2 > result_vec[i].l1 && result_vec[i].l2 > result_vec[i].l3)
        {
            bestVal = result_vec[i].l2;
            bestQuartet = 2;
        }
        if(result_vec[i].l3 > result_vec[i].l1 && result_vec[i].l3 > result_vec[i].l2)
        {
            bestVal = result_vec[i].l3;
            bestQuartet = 3;
        }
        float score = bestVal;
        #pragma omp critical
        switch(bestQuartet)
        {
            case 1: out << idx2code(pa, result_vec[i].a1) << "," 
                << idx2code(pa, result_vec[i].b1) << "|" 
                << idx2code(pa, result_vec[i].c1) << "," 
                << idx2code(pa, result_vec[i].d1) << ":" << score << std::endl; break;
            case 2: out << idx2code(pa, result_vec[i].a2) << "," 
                << idx2code(pa, result_vec[i].b2) << "|" 
                << idx2code(pa, result_vec[i].c2) << "," 
                << idx2code(pa, result_vec[i].d2) << ":" << score << std::endl; break;
            case 3: out << idx2code(pa, result_vec[i].a3) << "," 
                << idx2code(pa, result_vec[i].b3) << "|" 
                << idx2code(pa, result_vec[i].c3) << "," 
                << idx2code(pa, result_vec[i].d3) << ":" << score << std::endl; break;
            case -1: mspamstats::same_likelihood_errors++; num_quartets--; break;
            default: std::cerr << "There is no best quartet!" << std::endl;
        }
    }
    
    rax_free(result_vec);
    freeAllTheThings(adef, tr, cdta, rdta);
    return num_quartets;
}
#endif
