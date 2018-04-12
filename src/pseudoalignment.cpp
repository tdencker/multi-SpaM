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
#include "pseudoalignment.hpp"
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

std::string PseudoAlignment::to_string(std::vector<Sequence> & sequences)
{
    std::string str;
    str += "##################### Length: " + std::to_string(_length) + " ######################\n";
    for(auto & w : _words)
    {
        str += std::to_string(w.getSeq()) + " : " + std::to_string(std::distance(sequences[w.getSeq()].content.begin(), w.getPos())) + "\n";
    }
    return str;
}

inline void PseudoAlignment::setSequenceKey()
{
    _seq_key = 0;
    for(auto & e : _words)
    {
        _seq_key |= 1 << e.getSeq();
    }
    
}

// TODO: which constructors are actually used?
PseudoAlignment::PseudoAlignment(std::vector< std::vector<Word>::iterator > & vec, int pattern_length, uint64_t seq_key, bool aligned) : _length(pattern_length), _seq_key(seq_key), _merge_count(0), _was_aligned(aligned)
{
	assert(vec.size() >= options::min_sequences);
	assert(vec.size() <= 32 );
	_words.reserve(vec.size());
	for(auto & it : vec)
		_words.push_back(*it); // TODO: move?
	std::sort(_words.begin(), _words.end(), [](const Word & a, const Word & b){return a.getSeq() < b.getSeq();});
//	setSequenceKey();
}

PseudoAlignment::PseudoAlignment(std::vector<Word> && vec, int length, bool aligned) : _words(vec), _length(length), _merge_count(0), _was_aligned(aligned)
{
    setSequenceKey();
}

bool PseudoAlignment::checkBoundaries(std::vector<Sequence> & sequences)
{
        return checkBoundariesLeft(0, sequences) && checkBoundariesRight(_length - 1, sequences);
        
//        bool b = checkBoundariesLeft(0, sequences) && checkBoundariesRight(_length - 1, sequences);
//        if(b == false)
//        {
//                for(unsigned i = 0; i < _words.size(); ++i)
//                {
//                        std::cout << i << ": " << std::boolalpha << checkBoundariesSingle(i, sequences,0,0) << std::endl;
//                }
//        }
//        return b;
}

bool PseudoAlignment::checkBoundariesSingle(int seq, std::vector<Sequence> & sequences, int add_to_left, int length)
{
    if(length < 0)
        length = _length;
    if(_words[seq].revComp())
	{
		if(_words[seq].getPos() + add_to_left >= sequences[_words[seq].getSeq()].content.end() || _words[seq].getPos() + add_to_left - length + 1 < sequences[_words[seq].getSeq()].content.begin())
			return false;
	}
	else
	{
		if(_words[seq].getPos() - add_to_left < sequences[_words[seq].getSeq()].content.begin() || _words[seq].getPos() - add_to_left + length - 1 >= sequences[_words[seq].getSeq()].content.end())
			return false;
	}
	return true;
}

bool PseudoAlignment::checkBoundariesLeft(int pos, std::vector<Sequence> & sequences)
{
	assert(pos <= 0);
	for(unsigned i = 0; i < _words.size(); ++i)
	{
		if(_words[i].revComp())
		{
			if(_words[i].getPos() - pos >= sequences[_words[i].getSeq()].content.end())
				return false;
		}
		else
		{
			if(_words[i].getPos() + pos < sequences[_words[i].getSeq()].content.begin())
				return false;
		}
	}
	return true;
}

bool PseudoAlignment::checkBoundariesRight(int pos, std::vector<Sequence> & sequences)
{
	assert(pos >= 0);
	for(unsigned i = 0; i < _words.size(); ++i)
	{
		if(_words[i].revComp())
		{
			if(_words[i].getPos() - pos < sequences[_words[i].getSeq()].content.begin())
				return false;
		}
		else
		{
			if(_words[i].getPos() + pos >= sequences[_words[i].getSeq()].content.end())
				return false;
		}
	}
	return true;
}