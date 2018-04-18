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

#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>
#include "quartetblock.hpp"
#include "word.hpp"
#include "stats.hpp"

//test (mb remove)
#include <limits>
#include <unordered_map>

#include "options.hpp" // remove me
#include <cassert>
#include <iomanip>


static inline void reset(std::unordered_map<char,int> & map)
{
    // assert(map.size() == 4);
	for(auto & e : map)
		e.second = 0;
}

std::string QuartetBlock::to_string(std::vector<Sequence> & sequences)
{
    std::string str;
    str += "##################### Length: " + std::to_string(_length) + " ######################\n";
    for(auto & w : _words)
    {
        str += std::to_string(w.getSeq()) + " : " + std::to_string(std::distance(sequences[w.getSeq()].content.begin(), w.getPos())) + "\n";
    }
    return str;
}

inline void QuartetBlock::setSequenceKey()
{
    _seq_key = 0;
    for(auto & e : _words)
    {
        _seq_key |= 1 << e.getSeq();
    }
}

QuartetBlock::QuartetBlock(std::vector< std::vector<Word>::iterator > & vec, int pattern_length, uint64_t seq_key) 
    : _length(pattern_length), _seq_key(seq_key)
{
	assert(vec.size() >= mspamoptions::min_sequences);
	assert(vec.size() <= 32 );
	_words.reserve(vec.size());
	for(auto & it : vec)
		_words.push_back(*it);
	std::sort(_words.begin(), _words.end(), [](const Word & a, const Word & b)
    {
        return a.getSeq() < b.getSeq();
    });
}

QuartetBlock::QuartetBlock(std::vector<Word> && vec, int length) : _words(vec), _length(length)
{
    setSequenceKey();
}

QuartetBlock::QuartetBlock(int length) : _length(length), _seq_key(0)
{}