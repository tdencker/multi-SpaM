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

#include <algorithm>
#include <exception>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <limits>
#include <cassert>
#include "sequence.hpp"

std::vector<Sequence> Sequence::read(std::string & file_name, bool verbose)
{
	std::vector<Sequence> sequences;
	std::ifstream infile(file_name.c_str());
	assert(infile.is_open());
	std::string line;
	std::string header;
	std::getline(infile, line, '>');
	while (!infile.eof())
	{
		std::getline(infile, header);
		std::getline(infile, line, '>');
		sequences.emplace_back(header, line);
		if(verbose)
			std::cout << sequences.back().id << ": " << sequences.back().content.size() << " residues." << std::endl;
    }
	return sequences;
}

Sequence::Sequence(std::string header, std::string & seq) : id(header)
{
	std::unordered_map<char,char> map;
	map['A'] = 0;
	map['C'] = 1;
	map['G'] = 2;
	map['T'] = 3;
	content.reserve(seq.size());
	for (auto it = seq.begin(); it != seq.end(); it++) 
	{	
		if(!std::isspace(*it))
		{
			char c = std::toupper(*it);
			if(map.count(c) > 0)
				content.push_back(map[c]);
            else
                content.push_back(std::numeric_limits<char>::max());
		}
	}
	content.shrink_to_fit();

	size_t test_length = std::min(content.size(), (size_t) 100);
	if( std::count_if(content.begin(), content.begin() + test_length, 
		[](char c)
		{
			return c == std::numeric_limits<char>::max();
		}) 
		> 0.5 * test_length )
	{
		throw std::runtime_error("The sequence " + id + " is not a DNA sequence!");
	}
}
