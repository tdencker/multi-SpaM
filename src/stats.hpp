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

#ifndef STATS_HPP_
#define STATS_HPP_

/**
* @brief Stats that can be used to identify problems with runtime or
* the number of quartet trees for the final tree
**/

namespace mspamstats
{
extern uint64_t bad_characters;
extern uint64_t random_matches;
extern uint64_t ambigious_sequences;
extern uint64_t total_iterations;
extern uint64_t num_quartet_blocks;
extern uint64_t bad_model_errors;
extern uint64_t bad_character_error;
extern uint64_t not_all_nucleotides_error;
extern uint64_t same_likelihood_errors;
extern uint64_t num_quartet_trees;
void printStats();
}
#endif
