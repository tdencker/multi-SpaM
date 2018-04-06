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
 */#include <iomanip>
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

int64_t PseudoAlignment::getDiagonalHash(std::vector<Sequence> & sequences) const
{
	int64_t diagonal_key = 0;
	for(unsigned i = 1; i < _words.size(); ++i)
	{
		diagonal_key *= 37;
		auto d1 = _words[i-1].getPos(sequences); // are iterators really better than int for position
		auto d2 = _words[i].getPos(sequences);
		diagonal_key += d1 - d2;
	}
	return diagonal_key;
}

void PseudoAlignment::print(std::vector<Sequence> & sequences, std::ostream & stream) // remove this whole thing
{
	for(auto & w : _words)
	{
		stream << std::setw(3) << w.getSeq() << "\t" << std::setw(8) << w.getPos(sequences) << "\t" << w.to_string( (std::min)(_length, 250)) << "\t" << w.revComp() << std::endl;
	}
	stream << std::endl;
}

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

// requires that both pseudoalignments involve the same set of sequences
bool PseudoAlignment::hasOverlap(PseudoAlignment & other)
{
	if(_words[0].revComp() != other._words[0].revComp())
	{
		return false;
	}
	int dist = _words[0].getDistance(other._words[0]);
	if(dist > _length || dist < -other._length)
	{
		return false;
	}
	for(unsigned i = 1; i < _words.size(); ++i)
	{
		if(_words[i].revComp() != other._words[i].revComp())
		{
			return false;
		}
		int newdist = _words[i].getDistance(other._words[i]);
		if(newdist != dist)
		{
			return false;
		}
	}
	return true;
}

int PseudoAlignment::getPartialOverlap(PseudoAlignment & other) 
{	
	int overlap_count = 0;
	auto it = _words.begin();
    auto it2 = other._words.begin();
    while( it != _words.end() && it2 != other._words.end() )
    {
        if( it->revComp() == it2->revComp() && it->getSeq() == it2->getSeq() )
        {
            int dist = std::distance(it->getPos(), it2->getPos());
            if(dist <= _length && dist >= -other._length)
                overlap_count++;
        }
        if( it->getSeq() == it2->getSeq() )
        {
            it++;
            it2++;
            continue;
        }
        if( it->getSeq() > it2->getSeq() )
        {
            it2++;
        }
        else
        {
            it++;
        }
    }
    return overlap_count;
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

// just for testing oles stuff
void PseudoAlignment::ungappedExtension(int xdrop, std::vector<Sequence> & sequences, int max_mismatches)
{
    std::unordered_map<char,int> nuc_map;
	nuc_map.insert(std::make_pair(0, 0));
	nuc_map.insert(std::make_pair(1, 0));
	nuc_map.insert(std::make_pair(2, 0));
	nuc_map.insert(std::make_pair(3, 0));
	int mm_left = 0, mm_right = 0;
	int cdc = 0;
	
	while(mm_left + mm_right < max_mismatches)
	{
	    bool check_left = mm_left - mm_right != 0 ? mm_left < mm_right : cdc++ % 2 == 0;
	    while(1)
	    {
	        //TODO: be less lazy
	        bool left_boundaries = checkBoundariesLeft(-1, sequences);
	        bool right_boundaries = checkBoundariesRight( _length + 1, sequences);
	        if( left_boundaries == false && right_boundaries == false )
	            return;
	        if( check_left == true && left_boundaries == false )
	            check_left = false;
	        if( check_left == false && right_boundaries == false )
	            check_left = true;
	        
	        reset(nuc_map);
	        if(check_left == false)
	            _length++;
	        for(auto & w : _words)
	        {
	            if(check_left == true)
	                w.moveLeft(1);
	            int pos  = check_left ? 0 : _length;
	            nuc_map[w[pos]]++;
	        }
	        auto frequent_nuc = std::max_element(nuc_map.begin(), nuc_map.end(), 
	        [](const std::pair<char, int>& p1, const std::pair<char, int>& p2) {
                return p1.second < p2.second; });
	        if(_words.size() - frequent_nuc->second != 0) // there were mismatches
	        {
	            if(check_left == true)
	                mm_left += _words.size() - frequent_nuc->second;
	            else
	                mm_right += _words.size() - frequent_nuc->second;
	            // possibly go for the other direction
	            break;
	        }
        }
    }
}