/**
 * @file
 *
 * This is the main file. It contains all top level functions for
 * the game logic.
 *
 * @brief the main file
 * @author Thomas Dencker
 *
 * @section License
 *
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

#include <algorithm>
#include <cassert>
#include <chrono>
#include <climits>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "mergesort.hpp"
#include "options.hpp"
#include "pattern.hpp"
#include "quartetblock.hpp"
#include "randommatchfinder.hpp"
#include "rasbhari/variance.h"
#include "raxmlwrapper.hpp"
#include "sequence.hpp"
#include "stats.hpp"
#include "word.hpp"

constexpr int num_steps = 5;
constexpr int num_buckets = 256;
constexpr int left_column = 40;
constexpr int right_column = 13;

void printQuartets( std::vector<Sequence> & sequences, std::vector<QuartetBlock> & qb_vec, unsigned length)
{
    std::ofstream out_file( mspamoptions::output_file );
    assert( out_file.is_open() );
    for(auto & quartet : qb_vec)
    {
        for(auto & w_it : quartet)
        {
            out_file << ">" << sequences[w_it.getSeq()].id << " (Pos: " << std::distance(sequences[w_it.getSeq()].content.begin(), w_it.getPos()) << ")" << std::endl;
            out_file << w_it.toString(length) << std::endl;
        }
        out_file << std::endl << std::endl;
    }
}

/**
* @brief Runs RAxML and writes the resulting quartet trees in the outfile. For
* more detail, see raxmlwrapper.hpp.
**/

void RunRAxML( std::vector<QuartetBlock> & qb_vec )
{
    auto start = std::chrono::steady_clock::now();
    std::ofstream out_file( mspamoptions::output_file );
    assert( out_file.is_open() );

// compute all quartets

#pragma omp parallel for schedule( dynamic )
    for ( int i = 0; i < (int) qb_vec.size(); ++i )
    {
#ifdef _OPENMP
        if ( omp_get_thread_num() == 0 )
#endif
        {
            std::cout << "\r[Step 5 / " << num_steps << "] RAxML: " << i << " / " << qb_vec.size() << std::flush;
        }
#pragma omp atomic
        mspamstats::num_quartet_trees += computeAndPrintBestQuartets( qb_vec[i], out_file );
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << std::setw( left_column ) << std::left
              << ( "\r[Step 5 / " + std::to_string( num_steps ) + "] RaxML step" ) << std::setw( right_column )
              << std::right << std::setprecision( 2 ) << std::fixed << diff.count() << " seconds" << std::endl;
}

/**
* @brief This function samples quartet blocks from the word vector until
* mspamoptions::num_samples quartet blocks have been found. The sampling is done
* with the RandomMatchFinder.
**/

std::vector<QuartetBlock> samplingQuartetBlocks( std::vector<Word> & words, Pattern & current_pattern,
                                                 size_t num_sequences, unsigned & progress, int thread_id,
                                                 int thread_num, bool mem_save = false )
{
    unsigned limit = mspamoptions::num_samples / mspamoptions::num_patterns;

    if ( mem_save == true )
    {
        progress = 0;
        limit = std::ceil( limit / num_buckets );
    }

    auto start = std::chrono::steady_clock::now();

    RandomMatchFinder rmf( words, thread_id, thread_num );
    std::vector<QuartetBlock> qb_vec;

    size_t local_progress = 0;
    while ( local_progress < (limit/thread_num) )
    {
#pragma omp critical
        if ( mem_save == false && progress % 100 == 0 )
            std::cout << "\r[Step 4 / " << num_steps << "] Sampling: " << progress << " / " << mspamoptions::num_samples
                      << std::flush;

        try
        {
            qb_vec.push_back( rmf.next( current_pattern, num_sequences ) );
        }
        catch ( const std::exception & e )
        {
            break;
        }

        local_progress++;
#pragma omp atomic
        progress++;
    }

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

#pragma omp master
    if ( mem_save == false )
    {
        std::cout << std::setw( left_column ) << std::left
                  << ( "\r[Step 4 / " + std::to_string( num_steps ) + "] Sampled quartet blocks" )
                  << std::setw( right_column ) << std::right << std::setprecision( 2 ) << std::fixed << diff.count()
                  << " seconds" << std::endl;
    }

    return qb_vec;
}

/**
* @brief sorts the spaced words vector
**/

void sortSpacedWords( std::vector<Word> & words, bool mem_save = false )
{
    if ( mem_save == false )
    {
        std::cout << "[Step 3 / " << num_steps << "] Sorting spaced words ..." << std::flush;
    }
    auto start = std::chrono::steady_clock::now();
    // parallel merge sort sometimes allocates more space than necessary for
    // some reason...
    mergeSort( words.begin(), words.end() );
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

    if ( mem_save == false )
    {
        std::cout << std::setw( left_column ) << std::left
                  << ( "\r[Step 3 / " + std::to_string( num_steps ) + "] Sorted spaced words" )
                  << std::setw( right_column ) << std::right << std::setprecision( 2 ) << std::fixed << diff.count()
                  << " seconds" << std::endl;
    }
}

/**
* @brief See createSpacedWords for the basic information. The size of the vector
* is not precalculated here. The hash for the first 4 nucleotides has to be
* identical to the bucket unsigned integer, otherwise the word is not created.
**/

std::vector<Word> createSpacedWordsMemSave( std::vector<Sequence> & sequences, Pattern & current_pattern,
                                            bool compute_rev_comp, uint64_t bucket )
{
    std::vector<Word> words;

#pragma omp parallel for schedule( static )
    for ( int i = 0; i < (int) sequences.size(); ++i )
    {
        std::vector<Word> localWords;
        auto & seq = sequences[i].content;
        for ( auto it = seq.begin(); it != seq.end() - current_pattern.size() + 1; ++it )
        {
            const unsigned shift = mspamoptions::symbol_bits;
            uint64_t key = 0;
            for ( unsigned i = 0; i < 4; i++ )
            {
                unsigned nuc = 0;
                int offset = current_pattern[i];
                nuc = *( it + offset );

                if ( nuc == ( std::numeric_limits<char>::max )() )
                {
                    continue;
                }

                key <<= shift;
                key |= nuc;
            }
            if ( key != bucket )
            {
                continue;
            }
            localWords.push_back( Word( current_pattern, it, i, false ) );
        }
        if ( compute_rev_comp == true )
        {
            for ( auto it = seq.end() - 1; it != seq.begin() + current_pattern.size() - 2; --it )
            {
                const unsigned mask = mspamoptions::mask;
                const unsigned shift = mspamoptions::symbol_bits;
                uint64_t key = 0;
                for ( unsigned i = 0; i < 4; i++ )
                {
                    unsigned nuc = 0;
                    int offset = -current_pattern[i];

                    nuc = mask - *( it + offset );

                    if ( nuc == ( std::numeric_limits<char>::max )() )
                    {
                        continue;
                    }

                    key <<= shift;
                    key |= nuc;
                }
                if ( key != bucket )
                {
                    continue;
                }
                localWords.push_back( Word( current_pattern, it, i, true ) );
            }
        }
#pragma omp critical
        {
            words.insert( words.end(), localWords.begin(), localWords.end() );
        }
    }

    return words;
}

/**
* @brief The spaced words are created for all sequences and a given pattern. The
* size of the vector is precalculated. By default, the spaced words for the
* reverse complement are also created.
**/

std::vector<Word> createSpacedWords( std::vector<Sequence> & sequences, Pattern & current_pattern,
                                     bool compute_rev_comp = mspamoptions::use_rev_comp )
{
    auto start = std::chrono::steady_clock::now();
    std::vector<size_t> start_points( sequences.size() );
    start_points[0] = 0;
    int factor = compute_rev_comp ? 2 : 1;
    for ( unsigned i = 1; i < sequences.size(); ++i )
    {
        size_t effective_size = sequences[i-1].content.size() > current_pattern.size() ? ( sequences[i - 1].content.size() - current_pattern.size() + 1 ) : 0;
        start_points[i] = start_points[i - 1] + effective_size * factor;
    }
    size_t effective_size = sequences.back().content.size() > current_pattern.size() ? ( sequences.back().content.size() - current_pattern.size() + 1 ) : 0;
    size_t size = start_points[sequences.size() - 1] + effective_size * factor;
    std::vector<Word> words;
    words.resize( size );
    int finished_sequences = 0;

    std::cout << "\r[Step 2 / " << num_steps << "] Creating spaced words: " << finished_sequences << " / "
              << sequences.size() << " sequences" << std::flush;

#pragma omp parallel for schedule( static )
    for ( int i = 0; i < (int) sequences.size(); ++i )
    {
        auto & seq = sequences[i].content;
        if(seq.size() < current_pattern.size())
        {
            continue;
        }
        for ( auto it = seq.begin(); it != seq.end() - current_pattern.size() + 1; ++it )
        {
            words[start_points[i]++] = Word( current_pattern, it, i, false );
        }
        if ( compute_rev_comp == true )
        {
            for ( auto it = seq.end() - 1; it != seq.begin() + current_pattern.size() - 2; --it )
            {
                words[start_points[i]++] = Word( current_pattern, it, i, true );
            }
        }
#pragma omp critical
        std::cout << "\r[Step 2 / " << num_steps << "] Creating spaced words: " << ++finished_sequences << " / "
                  << sequences.size() << " sequences" << std::flush;
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << std::setw( left_column ) << std::left
              << ( "\r[Step 2 / " + std::to_string( num_steps ) + "] Created spaced words" )
              << std::setw( right_column ) << std::right << std::setprecision( 2 ) << std::fixed << diff.count()
              << " seconds" << std::endl;

    return words;
}

/**
* @brief See runStandard for main functionality. The memory saving mode only
* stores words that start with the same 4 nucleotides at the don't care
* positions. Thus, there are 256 iterations per pattern. The number of samples
* are divided by 256, so they are not gonna be exactly the same amount of blocks
* as for the standard run.
**/

std::vector<QuartetBlock> runMemSave( std::vector<Sequence> & sequences, std::vector<Pattern> & pattern_set )
{
    auto start = std::chrono::steady_clock::now();
    std::vector<QuartetBlock> qb_vec;
    unsigned progress = 0;

    for ( auto current_pattern : pattern_set )
    {
        for ( uint64_t bucket = 0; bucket < num_buckets; ++bucket )
        {
            std::cout << "\r[Step 2-4 / " << num_steps << "] Partition: " << bucket << " / " << num_buckets
                      << std::flush;
            std::vector<Word> all_spaced_words = createSpacedWordsMemSave( sequences, current_pattern, true, bucket );

            sortSpacedWords( all_spaced_words, true );

#pragma omp parallel
            {
                int thread_id = 0;
                int thread_num = 1;

#ifdef _OPENMP
                thread_id = omp_get_thread_num();
                thread_num = omp_get_num_threads();
#endif

                std::vector<QuartetBlock> thread_qb_vec = samplingQuartetBlocks(
                    all_spaced_words, current_pattern, sequences.size(), progress, thread_id, thread_num, true );

// reduction of the index (if necessary)
#pragma omp critical
                {
                    if ( thread_num > 1 || pattern_set.size() > 1 )
                    {
                        qb_vec.insert( qb_vec.end(), thread_qb_vec.begin(), thread_qb_vec.end() );
                    }
                    else
                    {
                        std::swap( thread_qb_vec, qb_vec );
                    }
                } // end of critical
            }     // end of parallel
        }
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << std::left << std::setw( left_column )
              << "\r[Step 2-4 / " + std::to_string( num_steps ) + "] Memory saving mode" << std::fixed
              << std::setw( right_column ) << std::right << std::setprecision( 2 ) << diff.count() << " seconds"
              << std::endl;

    assert( qb_vec.size() > 0 );
    return qb_vec;
}

/**
* @brief Creates a list of spaced words for every sequence in the sequence
* vector. The words are then sorted and passed to the sampling step. This
* process is repeated for every pattern in the pattern set.
**/

std::vector<QuartetBlock> runStandard( std::vector<Sequence> & sequences, std::vector<Pattern> & pattern_set )
{
    std::vector<QuartetBlock> qb_vec;
    unsigned progress = 0;

    for ( auto current_pattern : pattern_set )
    {
        std::vector<Word> all_spaced_words = createSpacedWords( sequences, current_pattern );
        sortSpacedWords( all_spaced_words );

#pragma omp parallel
        {
            int thread_id = 0;
            int thread_num = 1;

#ifdef _OPENMP
            thread_id = omp_get_thread_num();
            thread_num = omp_get_num_threads();
#endif

            std::vector<QuartetBlock> thread_qb_vec = samplingQuartetBlocks(
                all_spaced_words, current_pattern, sequences.size(), progress, thread_id, thread_num );

// reduction of the index (if necessary)
#pragma omp critical
            {
                if ( thread_num > 1 || pattern_set.size() > 1 )
                {
                    qb_vec.insert( qb_vec.end(), thread_qb_vec.begin(), thread_qb_vec.end() );
                }
                else
                {
                    std::swap( thread_qb_vec, qb_vec );
                }
            } // end of critical
        }     // end of parallel
    }
    if(qb_vec.size() == 0)
    {
        std::cerr << "No quartetblocks could be found!" << std::endl;
        exit(-1);
    }
    return qb_vec;
}

/**
* @brief Read the sequences from the FASTA file provided with the "-i" flag.
**/

std::vector<Sequence> readSequences()
{
    auto start = std::chrono::steady_clock::now();
    std::cout << "[Step 1 / " << num_steps << "] Reading sequences ..." << std::flush;
    std::vector<Sequence> sequences = Sequence::read( mspamoptions::input_file, false );
    if ( mspamoptions::all_sequences )
        mspamoptions::min_sequences = sequences.size();
    assert( sequences.size() < ( std::numeric_limits<uint16_t>::max )() );
    if(sequences.size() < mspamoptions::min_sequences)
    {
        std::cerr << "Only " << sequences.size() << " sequences were found. At least " << mspamoptions::min_sequences << " are necessary for multi-SpaM to run." << std::endl;
        exit(-1);
    }
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

    std::cout << std::left << std::setw( left_column ) << ( "\r[Step 1 / " + std::to_string( num_steps ) + "] Read " +
                                                            std::to_string( sequences.size() ) + " sequences" )
              << std::fixed << std::setw( right_column ) << std::right << std::setprecision( 2 ) << diff.count()
              << " seconds" << std::endl;

    return sequences;
}

/**
* @brief Initiatlize OMP with the number of threads.
**/

inline void initOMP()
{
#ifdef _OPENMP
    int threads = mspamoptions::num_threads;
    omp_set_dynamic( 0 );
    omp_set_num_threads( threads );
#endif
}

/**
* @brief Create one or multiple pattern with rasbhari.
**/

// TODO: use newer version of rasbhari
std::vector<Pattern> getPatternSet()
{
    if( mspamoptions::use_seed == true )
    {
        setRasbhariSeed( mspamoptions::seed );
    }
    std::vector<Pattern> pattern_set;
    variance var( mspamoptions::num_patterns, mspamoptions::weight, mspamoptions::dontcare, mspamoptions::dontcare );
    var.Init( true, true, true, true, false, NULL );
    var.Improve( 500 );
    for ( auto & e : var.GetPattern() )
    {
        std::string str( e.begin(), e.end() );
        pattern_set.emplace_back( str );
    }
    return pattern_set;
}

/**
* @brief main function
*
* The main function read paramaters, sequence files, creates random pattern(s)
* and calls the run function either with or without memory saving mode
**/

int main( int argc, char ** argv )
{
    mspamoptions::parseParameters( argc, argv );
    mspamoptions::printParameters();
    initOMP();
    std::vector<Sequence> sequences = readSequences();
    std::vector<Pattern> pattern_set = getPatternSet();

    std::vector<QuartetBlock> qb_vec =
        mspamoptions::mem_save_mode ? runMemSave( sequences, pattern_set ) : runStandard( sequences, pattern_set );
    mspamstats::num_quartet_blocks = qb_vec.size();

    // quartet calculation using raxml

    if ( mspamoptions::print_only == true)
    {
        printQuartets( sequences, qb_vec, pattern_set[0].size() );
    } else 
    {
        RunRAxML( qb_vec );
    }

    if ( mspamoptions::show_stats == true )
        mspamstats::printStats();
}
