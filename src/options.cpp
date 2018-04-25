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
#include "options.hpp"

namespace mspamoptions
{
int weight = 10;
int dontcare = 100;
int num_patterns = 1;
unsigned mask = 3;
unsigned symbol_bits = 2;
int min_score = 0;
unsigned min_sequences = 4;
unsigned num_samples = 1000000;
int num_threads = 1;
std::string input_file = "";
std::string output_file = "outfile";
bool all_sequences = false;
bool mem_save_mode = false;
bool show_stats = false;

void printHelp()
{
    // TODO: show all options
    std::cout << "-h, --help          show this help message and exit"
              << std::endl
              << "-i <string>         [Required] Input file (multi-fasta file) "
                 "(default: None)"
              << std::endl
              << "-o <string>         Output file (tree file in newick-format) "
                 "(default:"
              << std::endl
              << "                      outtree)" << std::endl
              << "-k <int>, -w <int>  weight of the pattern (the number of "
                 "matching positions)"
              << std::endl
              << "                      (default: 10)" << std::endl
              << "-d <int>            number of don't care positions in the "
                 "pattern (default:"
              << std::endl
              << "                      100)" << std::endl
              << "-n <int>            number of sampled quartet blocks "
                 "(default: 1000000)"
              << std::endl
              << "-t <int>            number of threads (default: 1)"
              << std::endl
              << "--mem-save          memory saving mode (default: False)"
              << std::endl
              << "--show-stats        additional stats (mostly for debugging) "
                 "(default: False)"
              << std::endl
              << "-v, --version       show program's version number and exit"
              << std::endl;
}

void printVersion()
{
    std::cout << "Multi-SpaM version 1.0" << std::endl
              << "Copyright (C) 2018 Thomas Dencker" << std::endl
              << "License GPLv3+: GNU GPL version 3 or later" << std::endl;
}

void printError( std::string message )
{
    std::cerr << message << std::endl;
    exit( 1 );
}

void parseParameters( int argc, char * argv[] )
{
    struct option long_options[] = {
        {"version", no_argument, nullptr, 'v'},
        {"mem-save", no_argument, nullptr, 0},
        {"show-stats", no_argument, nullptr, 0},
        {"all-sequences", no_argument, nullptr, 'a'},
        {"input", required_argument, nullptr, 'i'},
        {"help", no_argument, nullptr, 'h'},
        {"output", required_argument, nullptr, 'o'},
        {"threads", required_argument, nullptr, 't'},
        {"weight", required_argument, nullptr, 'w'},
        {"threads", required_argument, nullptr, 't'},
        {"dontcare", required_argument, nullptr, 'd'},
        {"min-sequences", required_argument, nullptr, 'm'},
        {"num-patterns", required_argument, nullptr, 'p'},
        {"num-samples", required_argument, nullptr, 'n'},
        {"min-score", required_argument, nullptr, 's'},
        {nullptr, 0, nullptr, 0}};

    while ( 1 )
    {

        int option_index = 0;
        int c = getopt_long( argc, argv, "vai:ho:t:w:t:d:m:p:n:s:",
                             long_options, &option_index );

        if ( c == -1 )
        {
            break;
        }

        switch ( c )
        {
        case 'v':
            printVersion();
            exit( 0 );
            break;
        case 'p':
            mspamoptions::num_patterns = std::atoi( optarg );
            if ( mspamoptions::num_patterns < 1 )
                printError( "Number of patterns (-p) must be an integer larger "
                            "than 0" );
            break;
        case 'd':
            mspamoptions::dontcare = atoi( optarg );
            if ( !isdigit( optarg[0] ) || mspamoptions::dontcare < 0 )
                printError( "Number of don't care positions (-d) must be a "
                            "positive integer" );
            break;
        case 'w':
            mspamoptions::weight = std::atoi( optarg );
            if ( mspamoptions::weight < 1 )
                printError( "Weight (-k) must be an integer larger than 1" );
            if ( mspamoptions::weight > 16 )
            {
                std::cout << "Weight is too large for a 32 bit integer. "
                             "Setting weight to 16."
                          << std::endl;
                mspamoptions::weight = 16;
            }
            break;
        case 't':
            mspamoptions::num_threads = std::atoi( optarg );
            if ( mspamoptions::num_threads < 1 )
                printError( "Threads (-t) must be an integer larger than 0" );
            break;
        case 'n':
            mspamoptions::num_samples = std::atoi( optarg );
            if ( mspamoptions::num_samples < 1 )
                printError(
                    "The number of samples (-n) must be larger than 0" );
            break;
        case 'o':
            mspamoptions::output_file = optarg;
            break;
        case 'i':
            mspamoptions::input_file = optarg;
            break;
        case 's':
            mspamoptions::min_score = std::atoi( optarg );
            break;
        case 'm':
            mspamoptions::min_sequences = std::atoi( optarg );
            if ( mspamoptions::min_sequences < 1 )
                printError(
                    "Min_sequences (-x) must be an integer larger than 1" );
            break;
        case 'a':
            all_sequences = true;
            break;
        case 0:
            if ( long_options[option_index].name == std::string( "mem-save" ) )
                mem_save_mode = true;
            else
                show_stats = true;
            break;
        case 'h':
            printHelp();
            exit( 0 );
            break;
        case '?':
            exit( 1 );
        }
    }

    if ( input_file == "" )
    {
        printError( "No input file specified! Exiting early ..." );
    }
}

void printParameters()
{
    constexpr int text_width = 50;
    constexpr int par_width = 10;
    std::cout << std::endl
              << "######################## Parameters ########################"
              << std::endl
              << std::setw( text_width ) << std::left
              << "Match positions (weight): " << std::setw( par_width )
              << std::right << weight << std::endl
              << std::setw( text_width ) << std::left
              << "Don't care positions: " << std::setw( par_width )
              << std::right << dontcare << std::endl
              << std::setw( text_width ) << std::left
              << "Number of patterns: " << std::setw( par_width ) << std::right
              << num_patterns << std::endl
              << std::setw( text_width ) << std::left
              << "Threads: " << std::setw( par_width ) << std::right
              << num_threads << std::endl
              << std::setw( text_width ) << std::left
              << "Number of samples: " << std::setw( par_width ) << std::right
              << std::to_string( num_samples ) << std::endl
              << "############################################################"
              << std::endl
              << std::endl;
}
}
