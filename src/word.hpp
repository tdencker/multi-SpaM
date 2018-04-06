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
 */#ifndef WORD_HPP_
#define WORD_HPP_

#include "options.hpp"
#include "pattern.hpp"
#include "sequence.hpp"
#include <vector>
#include <limits>
#include <cassert>
#include <string>
#include <iostream>

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
	void reverseIterator(int);
	void moveLeft(int);
	void moveRight(int);
	std::string to_string(unsigned);
	int getPos( std::vector<Sequence> & ) const; 
	void swapPos(Word &);
	int getDistance(Word &) const;
    
	private:
	std::vector<char>::iterator _pos;
	uint32_t _key;
	int16_t _seq;
	bool _rev_comp;
    static std::vector<char> dummy_vec;
};

#endif
