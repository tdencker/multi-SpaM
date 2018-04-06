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
 */#include "word.hpp"

static inline char reverse(char & c)
{
	return ~c & 3;
}

std::vector<char> Word::dummy_vec;

int Word::getPos( std::vector<Sequence> & sequences ) const
{
	auto dist = std::distance( sequences[_seq].content.begin(), _pos);
	if(_rev_comp)
		dist = sequences[_seq].content.size() - dist;
	return dist;
}

int Word::getDistance(Word & other) const
{
	int dist = std::distance(_pos, other._pos);
	if(_rev_comp)
		dist *= -1;
	return dist;
}

std::string Word::to_string(unsigned length)
{
    std::string buf;
    buf.reserve(length);
	if(_rev_comp)
	{
		for(unsigned i = 0; i < length; ++i)
		{
			buf[i] = reverse(*(_pos - i));
		}
	}else
	{
		for(unsigned i = 0; i < length; ++i)
		{
			buf[i] = *(_pos + i);
		}
	}
	for(unsigned i = 0; i < length; ++i)
	{
		switch(buf[i])
		{
			case 0: buf[i] = 'A';break;
			case 1: buf[i] = 'C';break;
			case 2: buf[i] = 'G';break;
			case 3: buf[i] = 'T';break;
			default: std::cerr << "Illegal character: " << buf[i] << std::endl; //exit(1); //TODO: uncomment when errors are fixed
		}
	}
	return buf;
}

void Word::reverseIterator(int pattern_length)
{
	int factor = _rev_comp ? -1 : 1;
	_pos += (pattern_length - 1) * factor;
	_rev_comp = ! _rev_comp;
}

Word::Word(Pattern & p, std::vector<char>::iterator it, int16_t seq, bool rev_comp) :  _pos(it), _key(0), _seq(seq), _rev_comp(rev_comp)
{
    const unsigned mask = options::mask;
    const unsigned shift = options::symbol_bits;
	for(unsigned i = 0; i < p.weight; i++)
	{
		unsigned nuc = 0;
		int offset = p[i];
		if(_rev_comp)
			offset *= -1;
			
		nuc = * (it + offset);
        
        if(nuc == (std::numeric_limits<char>::max)())
        {
            throw std::exception(); // illegal character
        }

		if(_rev_comp)
			nuc = mask - nuc;
		_key <<= shift;
		_key |= nuc;
	}
}

Word::Word() : _pos(dummy_vec.begin()), _key(0), _seq(0), _rev_comp(false)
{}

char Word::operator[](const int & idx) const 
{ 
	return _rev_comp ? reverse(*(_pos - idx)) : *(_pos + idx); 
}

bool Word::revComp() const 
{ 
    return _rev_comp;
}

bool Word::isDummy() const 
{
    return _pos == dummy_vec.begin();
}

bool Word::operator<(const Word & other) const
{ 
    return _key < other._key; 
}

bool Word::operator==(const Word & other) const 
{
    return _pos == other._pos;
}

int Word::getSeq() const 
{
    return _seq;
}

uint32_t Word::getKey() const 
{
    return _key;
}

std::vector<char>::iterator Word::getPos() const 
{
    return _pos;
}

void Word::moveLeft(int distance) 
{
    assert(distance >= 0);
    if(_rev_comp) 
        _pos += distance; 
    else 
        _pos -= distance;
}
void Word::moveRight(int distance) 
{
    assert(distance >= 0);
    if(_rev_comp) 
        _pos -= distance; 
    else 
        _pos += distance;
}

void Word::swapPos(Word & other) 
{ 
    std::swap(_pos, other._pos);
}
