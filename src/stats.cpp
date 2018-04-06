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
 */#include <iostream>
#include <iomanip>
#include "options.hpp"

namespace stats
{
    uint64_t matches;
    uint64_t ambigious_sequences;
    uint64_t accepted_components;
    uint64_t multiple_components;
    uint64_t inserted_pas;
    uint64_t fully_merged_pas;
    uint64_t partially_merged_pas;
    uint64_t removed_pas;
    uint64_t removed_sequences;
    uint64_t componentless_words;
    uint64_t score_computations;
    uint64_t bad_model_errors;
    uint64_t same_likelihood_errors;
    
    void printStats()
    {
        std::cout << std::endl << std::endl;
        std::cout << "### STATS ###" << std::endl;
        std::cout << std::setw(50) << std::left << ("Matches with at least " + std::to_string(options::min_sequences) + " spaced words: ") << std::setw(10) << std::right << matches << std::endl;
        std::cout << std::setw(50) << std::left << "Sequences with repeated spaced word: " << std::setw(10) << std::right << ambigious_sequences << std::endl;
        std::cout << std::setw(50) << std::left << ("Components with at least " + std::to_string(options::min_sequences) + " sequences: ") << std::setw(10) << std::right << accepted_components << std::endl;
        std::cout << std::setw(50) << std::left << "Multiple components for the same spaced word: " << std::setw(10) << std::right << multiple_components << std::endl;
        std::cout << std::setw(50) << std::left << "Inserted PseudoAlignments without overlaps: " << std::setw(10) << std::right << inserted_pas << std::endl;
        std::cout << std::setw(50) << std::left << "Fully merged PseudoAlignments: " << std::setw(10) << std::right << fully_merged_pas << std::endl;
        std::cout << std::setw(50) << std::left << "Partially merged PseudoAlignments: " << std::setw(10) << std::right << partially_merged_pas << std::endl;
        std::cout << std::setw(50) << std::left << "PseudoAlignments removed due to few sequences: " << std::setw(10) << std::right << removed_pas << std::endl;
        std::cout << std::setw(50) << std::left << "Sequences removed due to overlaps: " << std::setw(10) << std::right << removed_sequences << std::endl;
        std::cout << std::setw(50) << std::left << "Random matches (no component): " << std::setw(10) << std::right << componentless_words << std::endl;
        std::cout << std::setw(50) << std::left << "Score computations: " << std::setw(10) << std::right << score_computations << std::endl;
        std::cout << std::setw(50) << std::left << "Bad model errors: " << std::setw(10) << std::right << bad_model_errors << std::endl;
        std::cout << std::setw(50) << std::left << "Same likelihood errors: " << std::setw(10) << std::right << same_likelihood_errors << std::endl;
        std::cout << "#############" << std::endl << std::endl;
    }
}
