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
#ifndef PSEUDOALIGNMENT_HPP_
#define PSEUDOALIGNMENT_HPP_

#include <vector>
#include <string>
#include <iostream>
#include "sequence.hpp"
#include "word.hpp"

#include <bitset>

#include <cassert> // remove me

class PseudoAlignment
{
	public:
	typedef std::vector<Word>::iterator iterator;

	PseudoAlignment(std::vector< iterator > &, int, uint64_t);
	PseudoAlignment(std::vector<Word> &&, int);
    PseudoAlignment(int);

    bool operator<(const PseudoAlignment & other) const { return _seq_key < other._seq_key; }
    
	size_t size() const { return _words.size(); }
    
    uint64_t getSequenceKey() const;
	
	std::string getHeader(int idx, std::vector<Sequence> & sequences) { return sequences[_words[idx].getSeq()].id; }
	
	std::string to_string(std::vector<Sequence> &);
	
    void printFull(std::vector<Sequence> & sequences, std::ostream & out = std::cout)
    {
        for(auto & w : _words)
        {
            out << ">" << sequences[w.getSeq()].id << std::endl;
            out << w.toString(_length) << std::endl;
        }
    }
	std::vector<char>::iterator getPos(int idx) const { return _words[idx].getPos(); }
	int getLength() const { return _length;}
	int getSeq(int idx) const { return _words[idx].getSeq(); }
	bool getRevComp(int idx) const { return _words[idx].revComp(); }
	
    void clear() { _words.clear(); } 
    inline void setSequenceKey();
    
    PseudoAlignment::iterator begin() { return _words.begin(); }
    PseudoAlignment::iterator end() { return _words.end(); }
    
    
    
    void push_back(const Word & w) { _seq_key |= 1 << w.getSeq(); _words.push_back(w); }
    
	private:
	std::vector<Word> _words;
	int _length;
	uint64_t _seq_key;
};

#endif
