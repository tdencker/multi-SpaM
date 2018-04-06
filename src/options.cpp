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
 */#include "options.hpp"

namespace options
{
	int weight = 10;
	int dontcare = 100;
	int patterns = 1;
	unsigned mask = 3;
	unsigned symbol_bits = 2; 
	int min_score = 0;
	unsigned min_sequences = 4;
	unsigned nbr_samples = 1000000;
	int threads = 1;
	std::string input_file = "";
	std::string output_file = "outfile";
	bool all_sequences = false;
    bool mem_save_mode = false;
	
	void printHelp()
	{
		std::cout << "TODO" << std::endl;
	}
	
	void parseParameters(int argc, char *argv[])
	{
		int option_char;
		while ((option_char = getopt (argc, argv, "hp:n:d:k:t:o:s:i:m:ac")) != -1)
		{ 
			int i = 0;
			switch (option_char)
			{  
			case 'p': 
				options::patterns = std::atoi (optarg); 
				if(options::patterns<1) options::printError("Number of patterns (-p) must be an integer larger than 0");
				break;
			case 'd': 
				options::dontcare = atoi (optarg); 
				if(!isdigit(optarg[0]) || options::dontcare<0) options::printError("Number of don't care positions (-d) must be a positive integer");
				break;
			case 'k': 
				options::weight = std::atoi (optarg); 
				if(options::weight < 1) options::printError("Weight (-k) must be an integer larger than 1");
				if(options::weight > 16) { std::cout << "Weight is too large for a 32 bit integer. Setting weight to 16." << std::endl; options::weight = 16;}
				break;
			case 't': 
				options::threads = std::atoi (optarg); 
				if(options::threads<1) options::printError("Threads (-t) must be an integer larger than 0");
				break;
			case 'n':
				i = std::atoi(optarg);
				if(i < 1) options::printError("The number of samples (-n) must be larger than 0");
				else options::nbr_samples = i;
				break;
			case 'o': 
				options::output_file = optarg; 
				break;
			case 'i':
				options::input_file = optarg;
				break;
			case 's':
				options::min_score = std::atoi(optarg);
				break;
			case 'm':
				i = std::atoi(optarg);
				if(i < 1) options::printError("Min_sequences (-m) must be an integer larger than 1");
				options::min_sequences = (unsigned) i;
				break;
			case 'a':
				all_sequences = true;
				break;
            case 'c':
                mem_save_mode = true;
                break;
			case 'h': 
					printHelp();
					exit (0);
		
				break;
			case '?': 		
				exit (1);
			}
		}
        
        if(input_file == "")
        {
            printError("No input file specified! Exiting early ...");
        }
	}
    
    void printParameters()
    {
        constexpr int text_width = 30;
        constexpr int par_width = 10;
        std::cout << std::endl << "############## Parameters ##############" << std::endl
        << std::setw(text_width) << std::left << "Match positions (weight): " << std::setw(par_width) << std::right << weight << std::endl
        << std::setw(text_width) << std::left << "Don't care positions: " << std::setw(par_width) << std::right << dontcare << std::endl
        << std::setw(text_width) << std::left << "Number of patterns: " << std::setw(par_width) << std::right << patterns << std::endl
        << std::setw(text_width) << std::left << "Threads: " << std::setw(par_width) << std::right << threads << std::endl
        << std::setw(text_width) << std::left << "Min score: " << std::setw(par_width) << std::right << min_score << std::endl
        << std::setw(text_width) << std::left << "Min sequences: " << std::setw(par_width) << std::right << min_sequences << std::endl
        << std::setw(text_width) << std::left << "Number of samples: " << std::setw(par_width) << std::right << std::to_string(nbr_samples) << std::endl << std::endl;
    }
	
	void printError(std::string message)
	{
		std::cerr << message << std::endl;
		exit(1);
	}
}
