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

#include <iostream>
#include <iomanip>
#include "options.hpp"

namespace mspamstats
{
    uint64_t bad_characters;
    uint64_t random_matches;
    uint64_t ambigious_sequences;
    uint64_t total_iterations;
    uint64_t num_quartet_blocks;
    uint64_t bad_model_errors;
    uint64_t bad_character_error;
    uint64_t not_all_nucleotides_error;
    uint64_t same_likelihood_errors;
    uint64_t num_quartet_trees;
    
    void printStats()
    {
        std::cout << std::endl << std::endl;
        std::cout << "########################## STATS ###########################" << std::endl;
        std::cout << std::setw(50) << std::left << "Bad characters in input sequences: " << std::setw(10) 
        << std::right << bad_characters << std::endl;
        std::cout << std::setw(50) << std::left << "Random matches (score below threshold): " << std::setw(10) 
        << std::right << random_matches << std::endl;
        std::cout << std::setw(50) << std::left << "Multiple spaced words in the same sequence: " << std::setw(10) 
        << std::right << ambigious_sequences << std::endl;
        std::cout << std::setw(50) << std::left << "Total number of iterations: " << std::setw(10) 
        << std::right << total_iterations << std::endl;
        std::cout << std::setw(50) << std::left << "Number of quartet blocks: " << std::setw(10) 
        << std::right << num_quartet_blocks << std::endl;
        std::cout << std::setw(50) << std::left << "Bad model errors (recoverable): " << std::setw(10) 
        << std::right << bad_model_errors << std::endl;
        std::cout << std::setw(50) << std::left << "Bad character (at don't care positions) errors: " << std::setw(10) 
        << std::right << bad_character_error << std::endl;
        std::cout << std::setw(50) << std::left << "Not all 4 nucleotides errors: " << std::setw(10) 
        << std::right << not_all_nucleotides_error << std::endl;
        std::cout << std::setw(50) << std::left << "Same likelihood errors: " << std::setw(10) 
        << std::right << same_likelihood_errors << std::endl;
        std::cout << std::setw(50) << std::left << "Number of quartet trees: " << std::setw(10) 
        << std::right << num_quartet_trees << std::endl;
        std::cout << "############################################################" << std::endl << std::endl;
    }
}
