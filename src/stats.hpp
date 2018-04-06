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
 */#ifndef STATS_HPP_
#define STATS_HPP_

namespace stats
{
    extern uint64_t matches;
    extern uint64_t ambigious_sequences;
    extern uint64_t accepted_components;
    extern uint64_t multiple_components;
    extern uint64_t inserted_pas;
    extern uint64_t fully_merged_pas;
    extern uint64_t partially_merged_pas;
    extern uint64_t removed_pas;
    extern uint64_t removed_sequences;
    extern uint64_t componentless_words;
    extern uint64_t score_computations;
    extern uint64_t bad_model_errors;
    extern uint64_t same_likelihood_errors;
    void printStats();
}
#endif
