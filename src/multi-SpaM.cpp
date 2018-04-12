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

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <map>
#include <memory>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <climits>
#include <chrono>
#include <unordered_set>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "pattern.hpp"
#include "sequence.hpp"
#include "word.hpp"
#include "matchfinder.hpp"
#include "pseudoalignment.hpp"
#include "rasbhari/variance.h"
#include "options.hpp"
// #include "mergesort.hpp"
#include "stats.hpp"
// #include "wordarray.hpp"
#include "raxmlwrapper.hpp"
#include "randommatchfinder.hpp"

#include <numeric>

constexpr int steps = 5;

void RunRAxML(std::vector<PseudoAlignment> & pa_vec)
{
    auto start = std::chrono::steady_clock::now();
    std::cout << "[Step 5 / " << steps << "] Calculating optimal quartet trees for block: " << std::flush;
    std::ofstream out_file(options::output_file);
    assert(out_file.is_open());
    
    // compute all quartets
    
    uint64_t nbr_quartets = 0;
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < (int)pa_vec.size(); ++i)
    {
        #ifdef _OPENMP
        if(omp_get_thread_num() == 0)
        #endif
        std::cout << "\r[Step 5 / " << steps << "] Calculating optimal quartet trees for block: " << i << " / " << pa_vec.size() << " ( Length: " << pa_vec[i].getLength() << " )              " << std::flush;

        #pragma omp atomic
        nbr_quartets += computeAndPrintBestQuartets(pa_vec[i],  out_file);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end-start;
    #pragma omp master
    std::cout << "\r[Step 5 / " << steps << "] Calculating optimal quartet trees in " << diff.count() << " seconds.                            " << std::endl;
}

std::vector<PseudoAlignment> samplingPseudoAlignmentsMemSave(std::vector<Word> &words, Pattern & current_pattern, size_t nbr_sequences, unsigned & progress, int thread_id, int thread_num)
{
    #pragma omp master
    std::cout << "[Step 4 / " << steps << "] Sampling blocks ..." << std::flush;
    auto start = std::chrono::steady_clock::now();
    
    RandomMatchFinder rmf(words, thread_id, thread_num);
    std::vector<PseudoAlignment> pa_vec;
    progress = 0; // TODO: change the progress to something useful
    while(progress < (options::nbr_samples / (options::patterns * 256)) ) // TODO: bucket stuff
    {
        #pragma omp master
        if(progress % 100 == 0) // TODO: doesn't work well for many threads
            std::cout << "\r[Step 4 / " << steps << "] Sampling blocks ... Progress: " << progress << " / " << options::nbr_samples << std::flush;

        try
        {
            pa_vec.push_back(rmf.next(current_pattern, nbr_sequences));
        }
        catch(const std::exception & e)
        {
            break;
        }

        #pragma omp atomic
        progress++;
    }
    
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end-start;
    #pragma omp master
    std::cout << "\r[Step 4 / " << steps << "] Sampling blocks in " << diff.count() << " seconds.                            " << std::endl;
    return pa_vec;
}

std::vector<PseudoAlignment> samplingPseudoAlignments(std::vector<Word> &words, Pattern & current_pattern, size_t nbr_sequences, unsigned & progress, int thread_id, int thread_num)
{
    #pragma omp master
    std::cout << "[Step 4 / " << steps << "] Sampling blocks ..." << std::flush;
	auto start = std::chrono::steady_clock::now();
    
    RandomMatchFinder rmf(words, thread_id, thread_num);
    std::vector<PseudoAlignment> pa_vec;

    while(progress < options::nbr_samples)
    {
        #pragma omp master
        if(progress % 100 == 0) // TODO: doesn't work well for many threads
            std::cout << "\r[Step 4 / " << steps << "] Sampling blocks ... Progress: " << progress << " / " << options::nbr_samples << std::flush;

        try
        {
            pa_vec.push_back(rmf.next(current_pattern, nbr_sequences));
        }
        catch(const std::exception & e)
        {
            break;
        }

        #pragma omp atomic
        progress++;
    }
    
    auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;
    #pragma omp master
    std::cout << "\r[Step 4 / " << steps << "] Sampling blocks in " << diff.count() << " seconds.                            " << std::endl;
    return pa_vec;
}

std::vector<PseudoAlignment> findAllPseudoAlignments(std::vector<Word> & words, Pattern & current_pattern, size_t nbr_sequences, unsigned & progress, int thread_id, int thread_num)
{
    #pragma omp master
    std::cout << "[Step 4 / " << steps << "] Finding all blocks ... " << std::flush;
	auto start = std::chrono::steady_clock::now();
    
    MatchFinder finder(words, thread_id, thread_num);
    std::vector<PseudoAlignment> pa_vec;
    unsigned total_mil = words.size() / 1000000;
    unsigned current_mil = -1;
    
    while(finder.next()) // while there is a match
    {
        #pragma omp atomic
        progress += finder.progress();
        #pragma omp master
        {
        if( ( progress / 1000000 ) != current_mil )
        {
            current_mil++;
            std::cout << "\r[Step 4 / " << steps << "] Finding all blocks ... Progress: " << current_mil << " / " << total_mil << " million words" << std::flush;
        }
        }
        std::vector<Component> components = finder.getCurrentComponents(current_pattern, nbr_sequences); // get all matches with positive score
        for(auto & comp : components)
        {
            pa_vec.emplace_back(current_pattern.size());
            for(auto & w : comp)
            {
                pa_vec.back().push_back(*w);
            }
            pa_vec.back().sort();
        }
    }
    
    auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;
    #pragma omp master
    std::cout << "\r[Step 4 / " << steps << "] Finding all blocks in " << diff.count() << " seconds.                       " << std::endl;
    return pa_vec;
}

void sortSpacedWords(std::vector <Word> &words)
{
    std::cout << "[Step 3 / " << steps << "] Sorting spaced words ..." << std::flush;
	auto start = std::chrono::steady_clock::now();
	std::sort(words.begin(), words.end());
    // parallel merge sort sometimes allocates more space than necessary for some reason...
    // words.sort() // parallel merge sort
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;
    std::cout << "\r[Step 3 / " << steps << "] Sorted spaced words in " << diff.count() << " seconds." << std::endl;
}

std::vector<Word> createSpacedWordsMemSave(std::vector<Sequence> & sequences, Pattern & current_pattern, bool compute_rev_comp, uint64_t bucket)
{
    auto start = std::chrono::steady_clock::now();
    std::vector<Word> words;
    int finished_sequences = 0;
    std::cout << "\r[Step 2 / " << steps << "] Creating spaced words: " << finished_sequences << " / " << sequences.size() << " sequences for bucket: " << bucket << std::flush;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < (int)sequences.size(); ++i)
    {
        std::vector<Word> localWords;
        auto & seq = sequences[i].content;
        for (auto it = seq.begin(); it != seq.end() - current_pattern.size() + 1; ++it)
        {
            const unsigned shift = options::symbol_bits;
            uint64_t key = 0;
            for (unsigned i = 0; i < 4; i++)
            {
                unsigned nuc = 0;
                int offset = current_pattern[i];
                nuc = *(it + offset);

                if (nuc == (std::numeric_limits<char>::max)())
                {
                    continue;
                }

                key <<= shift;
                key |= nuc;
            }
            if (key != bucket)
            {
                continue;
            }
            try {
                localWords.push_back( Word(current_pattern, it, i, false) );
            }
            catch (std::exception & e) {}
        }
        if (compute_rev_comp == true)
        {
            for (auto it = seq.end() - 1; it != seq.begin() + current_pattern.size() - 2; --it)
            {
                const unsigned mask = options::mask;
                const unsigned shift = options::symbol_bits;
                uint64_t key = 0;
                for (unsigned i = 0; i < 4; i++)
                {
                    unsigned nuc = 0;
                    int offset = -current_pattern[i];

                    nuc = mask - *(it + offset);

                    if (nuc == (std::numeric_limits<char>::max)())
                    {
                        continue;
                    }

                    key <<= shift;
                    key |= nuc;
                }
                if (key != bucket)
                {
                    continue;
                }
                try {
                    localWords.push_back(Word(current_pattern, it, i, true));
                }
                catch (std::exception & e) {}
            }
        }
#pragma omp critical
        {
            words.insert(words.end(), localWords.begin(), localWords.end());
            std::cout << "\r[Step 2 / " << steps << "] Creating spaced words: " << ++finished_sequences << " / " << sequences.size() << " sequences for bucket: " << bucket << std::flush;
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "\r[Step 2 / " << steps << "] Created spaced words in " << diff.count() << " seconds." << std::endl;
    return words;
}

std::vector<Word> createSpacedWords(std::vector<Sequence> & sequences, Pattern & current_pattern, bool compute_rev_comp = true)
{
    auto start = std::chrono::steady_clock::now();
    std::vector<size_t> start_points(sequences.size());
    start_points[0] = 0;
    int factor = compute_rev_comp ? 2 : 1;
    for(unsigned i = 1; i < sequences.size(); ++i)
    {
        start_points[i] = start_points[i - 1] + (sequences[i - 1].content.size() - current_pattern.size() + 1) * factor;
    }
    size_t size = start_points[sequences.size() - 1] + (sequences.back().content.size() - current_pattern.size() + 1) * factor;
    std::vector<Word> words;
    words.resize(size);
	int finished_sequences = 0;
    std::cout << "\r[Step 2 / " << steps << "] Creating spaced words: " <<  finished_sequences << " / " << sequences.size() << " sequences" << std::flush;
	#pragma omp parallel for schedule(static)
	for(int i = 0; i < (int)sequences.size(); ++i)
	{
		auto & seq = sequences[i].content;
 		for(auto it = seq.begin(); it != seq.end() - current_pattern.size() + 1; ++it)
		{
            try {
                words[start_points[i]++] = Word(current_pattern, it, i, false);
            }
            catch (std::exception & e) {}
		}
		if(compute_rev_comp == true)
		{
		    for(auto it = seq.end() - 1; it != seq.begin() + current_pattern.size() - 2; --it)
		    {
                try {
                    words[start_points[i]++] = Word(current_pattern, it, i, true);
                }
                catch (std::exception & e) {}
		    }
		}
		#pragma omp critical
		std::cout << "\r[Step 2 / " << steps << "] Creating spaced words: " <<  ++finished_sequences << " / " << sequences.size() << " sequences" << std::flush;
	}
    auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;
	std::cout << "\r[Step 2 / " << steps << "] Created spaced words in " << diff.count() << " seconds." << std::endl;
    return words;
}

std::vector<Sequence> readSequences()
{
	auto start = std::chrono::steady_clock::now();
	std::cout << "[Step 1 / " << steps << "] Reading sequences." << std::flush;
	std::vector<Sequence> sequences = Sequence::read(options::input_file, false);
	if(options::all_sequences)
		options::min_sequences = sequences.size();
	assert(sequences.size() < (std::numeric_limits<uint16_t>::max)());
	assert(sequences.size() >= options::min_sequences);
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> diff = end-start;
	std::cout << "\r[Step 1 / " << steps << "] Read " << sequences.size() << " sequences in " << diff.count() << " seconds." << std::endl;
    return sequences;
}

inline void initOMP()
{
    #ifdef _OPENMP
	int threads = options::threads;
	omp_set_dynamic(0);
	omp_set_num_threads(threads);
	#endif
}

// TODO: neues rasbhari
inline std::vector<Pattern> getPatternSet()
{
    std::vector<Pattern> pattern_set;
	variance var(options::patterns, options::weight, options::dontcare, options::dontcare);
	var.Init(true, true, true, true, false, NULL);
	var.Improve(500);
	for(auto & e : var.GetPattern())
	{
		std::string str(e.begin(), e.end());
		pattern_set.emplace_back(str);
	}
	return pattern_set;
}

std::vector<PseudoAlignment> runMemSave(std::vector<Sequence> & sequences, std::vector<Pattern> & pattern_set)
{
    std::vector<PseudoAlignment> pa_vec;
    unsigned progress = 0;
    
    for(auto current_pattern : pattern_set)
    {
        uint64_t bucketCount = 256;
        for (uint64_t bucket = 0; bucket < bucketCount; ++bucket)
        {
            std::vector<Word> all_spaced_words = createSpacedWordsMemSave(sequences, current_pattern, true, bucket);
            sortSpacedWords(all_spaced_words);

#pragma omp parallel
            {
                int thread_id = 0;
                int thread_num = 1;

#ifdef _OPENMP
                thread_id = omp_get_thread_num();
                thread_num = omp_get_num_threads();
#endif

                std::vector<PseudoAlignment> thread_pa_vec = samplingPseudoAlignmentsMemSave(all_spaced_words, current_pattern, sequences.size(), progress, thread_id, thread_num);

                // reduction of the index (if necessary)
#pragma omp critical
                {
                    if (thread_num > 1 || pattern_set.size() > 1)
                    {
                        pa_vec.insert(pa_vec.end(), thread_pa_vec.begin(), thread_pa_vec.end());
                    }
                    else
                    {
                        std::swap(thread_pa_vec, pa_vec);
                    }
                } // end of critical
            } // end of parallel
        }
        
    }
    assert(pa_vec.size() > 0);
    return pa_vec;
}

std::vector<PseudoAlignment> runStandard(std::vector<Sequence> & sequences, std::vector<Pattern> & pattern_set)
{
    std::vector<PseudoAlignment> pa_vec;
    unsigned progress = 0;
    
    for(auto current_pattern : pattern_set)
    {
        std::vector<Word> all_spaced_words = createSpacedWords(sequences, current_pattern);
        sortSpacedWords(all_spaced_words);
        
        #pragma omp parallel
        {
        int thread_id = 0;
        int thread_num = 1;
        
        #ifdef _OPENMP
        thread_id = omp_get_thread_num();
        thread_num = omp_get_num_threads();
        #endif
        
        std::vector<PseudoAlignment> thread_pa_vec = samplingPseudoAlignments(all_spaced_words, current_pattern, sequences.size(), progress, thread_id, thread_num);

        // reduction of the index (if necessary)
        #pragma omp critical
        {
        if(thread_num > 1 || pattern_set.size() > 1)
        {   
            pa_vec.insert(pa_vec.end(), thread_pa_vec.begin(), thread_pa_vec.end());
        }else
        {
            std::swap(thread_pa_vec, pa_vec);
        }
        } // end of critical
        } // end of parallel
    }
    assert(pa_vec.size() > 0);
    return pa_vec;
}

// TODO: try filling an array before calculation of match score
// TODO: try using matchpos instead of ismatch when calculating the score
// TODO: is there a way to use chris' bucket like sorting to be faster?
// TODO: make a hashmap of all words -> components that were removed because of too many revcomp words, just reverse the components and use them automatically without computing the scores again
// or even better: just reverse the words in the component and mark the reverse complement word to be skipped (does this work?)
// TODO: test if it works without linking openmp
// TODO: maybe check which version of the component ( more or less revcomp ) has less conflicts (defaulting to the non revcomp version) and insert that one
// TODO: maybe make an option for "aligning" / mapping instead of removing uncertainties, but it needs to be rethought
// TODO: also on option: compute all match scores for more accuracy and full connectivity can be enforced
// TODO: finally test different score computation strategies (with array etc.)
int main(int argc, char** argv)
{
	options::parseParameters(argc, argv);
    options::printParameters();
	initOMP();
	std::vector<Sequence> sequences = readSequences();
    std::vector<Pattern> pattern_set = getPatternSet();
    
    std::vector<PseudoAlignment> pa_vec = options::mem_save_mode ? runMemSave(sequences, pattern_set) : runStandard(sequences, pattern_set);

    // ungappedExtension(pa_vec, sequences);
    
    // quartet calculation using raxml
    
    RunRAxML(pa_vec);

    stats::printStats();
}

