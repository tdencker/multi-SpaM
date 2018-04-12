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

class OverlapInfo;

class PseudoAlignment
{
	public:
	typedef std::vector<Word>::iterator iterator;
	PseudoAlignment(std::vector< std::vector<Word>::iterator > &, int, uint64_t, bool = false);
	PseudoAlignment(std::vector<Word> &&, int, bool = false);
//	bool operator<(const PseudoAlignment & other) const { return _words[0].getPos() < other._words[0].getPos(); } // TODO: this might not be enough if revcomp
    bool operator<(const PseudoAlignment & other) const { return _seq_key < other._seq_key; }
    
	size_t size() const { return _words.size(); }
	
	bool checkBoundaries(std::vector<Sequence> &);
	bool checkBoundariesSingle(int, std::vector<Sequence> &, int = 0, int = -1);
	bool checkBoundariesLeft(int, std::vector<Sequence> &);
	bool checkBoundariesRight(int, std::vector<Sequence> &);
	
	//static bool hashCompare(const PseudoAlignment & pa1, const PseudoAlignment & pa2, std::vector<Sequence> & sequences) { return pa1.getDiagonalHash(sequences) < pa2.getDiagonalHash(sequences); } // TODO: make this better if it works
	
	uint64_t getSequenceKey() const;
	
	std::string getHeader(int idx, std::vector<Sequence> & sequences) { return sequences[_words[idx].getSeq()].id; }
	
	std::string to_string(std::vector<Sequence> &);
	
	std::string findSeqAndPrint(int seq, std::vector<Sequence> & sequences) // TODO: clean all these functions up
	{
	    for(auto & e : _words)
	        if(e.getSeq() == seq)
	        {
	            assert(checkBoundaries(sequences));
	            return e.toString(_length);
	        }
	    for(auto & e : _words)
	    {
	        std::cout << e.getSeq() << " ";
	    }
	    std::cout << " | " << seq << std::endl;
	    exit(1);
	}
	
    void printFull(std::vector<Sequence> & sequences, std::ostream & out = std::cout)
    {
        for(auto & w : _words)
        {
            out << ">" << sequences[w.getSeq()].id << std::endl;
            out << w.toString(_length) << std::endl;
        }
    }
	//remove maybe
	std::vector<char>::iterator getPos(int idx) const { return _words[idx].getPos(); }
	int getLength() const { return _length;}
	int getSeq(int idx) const { return _words[idx].getSeq(); }
	bool getRevComp(int idx) const { return _words[idx].revComp(); }
	unsigned getMergeCount() const { return _merge_count; }
	bool wasAligned() const { return _was_aligned; }
	
	//for paindex
//	std::vector<Word>::iterator begin() { return _words.begin(); }
//	std::vector<Word>::iterator end() { return _words.end(); }
//    Word & operator[](std::size_t idx) { return _words[idx]; }
    void extend(int, int, int, int);
    void removeSequences(std::vector<unsigned> &);
    void clear() { _words.clear(); } //_words.shrink_to_fit(); }
    inline void setSequenceKey();
    
    PseudoAlignment::iterator begin() { return _words.begin(); }
    PseudoAlignment::iterator end() { return _words.end(); }
    
    PseudoAlignment(int length) : _length(length), _seq_key(0) {}
    
    void push_back(const Word & w) { _seq_key |= 1 << w.getSeq(); _words.push_back(w); }
    
    void sort() { std::sort(_words.begin(), _words.end(), [](const Word & a, const Word & b) { return a.getSeq() < b.getSeq(); }); };
    
    void extend(const OverlapInfo & info);
    
    
    void erase(iterator pos) { _words.erase(pos); }
    
    
    bool hasIterator(std::vector<char>::iterator it)
    {
        for(auto & e : _words)
            if(e.getPos() == it)
                return true;
        return false;
    }
	
	private:
	std::vector<Word> _words;
	int _length;
	uint64_t _seq_key;
	// TODO: make sure these are actually used for anything useful
	unsigned _merge_count;
	bool _was_aligned;
};

#endif
