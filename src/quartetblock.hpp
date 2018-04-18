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
#ifndef QUARTETBLOCK_HPP_
#define QUARTETBLOCK_HPP_

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "sequence.hpp"
#include "stats.hpp"
#include "word.hpp"


/**
* @brief This class contains the words that make up a quartet block. It is assumed that there is only one
* Word per sequence. This class was previously used for more tasks and has, therefore, a bunch of convinience
* functions. The words are copied from the original vector and can be deleted after pushing them into the quartet
* block.
**/

class QuartetBlock
{
	public:
	typedef std::vector<Word>::iterator iterator;

	QuartetBlock(std::vector< iterator > &, int, uint64_t);
	QuartetBlock(std::vector<Word> &&, int);
    QuartetBlock(int);

    void push_back(const Word & w);

    // convenience functions
    bool operator<(const QuartetBlock & other) const;
	size_t size() const;
    uint64_t getSequenceKey() const;
	std::string getHeader(int idx, std::vector<Sequence> & );
	std::string to_string(std::vector<Sequence> &);
	std::vector<char>::iterator getPos(int idx) const;
	int getLength() const;
	int getSeq(int idx) const;
	bool getRevComp(int idx) const;
    void clear();
    void setSequenceKey();
    
    QuartetBlock::iterator begin();
    QuartetBlock::iterator end();

	private:
	std::vector<Word> m_words;
	int m_length;
	uint64_t m_seq_key;
};

#endif
