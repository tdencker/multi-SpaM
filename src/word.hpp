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

#ifndef WORD_HPP_
#define WORD_HPP_

#include "options.hpp"
#include "pattern.hpp"
#include "sequence.hpp"

#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

class Word
{
	public:
    Word();
	Word(Pattern & p, std::vector<char>::iterator , int16_t , bool ); 
	bool revComp() const;
    bool isDummy() const;
	char operator[](const int & ) const;
	bool operator<(const Word & ) const;
	bool operator==(const Word & ) const;
	int getSeq() const;
	uint32_t getKey() const;
	std::vector<char>::iterator getPos() const;
	std::string toString(unsigned);
    
	private:
	std::vector<char>::iterator m_pos;
	uint32_t m_key;
	int16_t m_seq;
	bool m_rev_comp;
    static std::vector<char> dummy_vec;
};

#endif
