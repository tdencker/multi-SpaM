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

#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <string>

namespace mspamoptions
{
extern int weight;
extern int dontcare;
extern int num_patterns;
extern unsigned mask;
extern unsigned symbol_bits;
extern int min_score;
extern unsigned min_sequences;
extern int num_threads;
extern unsigned num_samples;
extern uint64_t seed;
extern std::string input_file;
extern std::string output_file;
extern bool all_sequences;
extern bool mem_save_mode;
extern bool show_stats;
extern bool print_only;
extern bool use_seed;
void parseParameters( int argc, char * argv[] );
void printParameters();
}

#endif
