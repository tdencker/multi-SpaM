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

#include "pattern.hpp"

Pattern::Pattern(std::string && str)
{
	_boolPattern.reserve(str.size());
	for(unsigned k = 0; k < str.size(); ++k)
	{
		if(str[k] == '1')
			_matchPos.push_back(k);
		_boolPattern.push_back(str[k] == '1');
	}
}

Pattern::Pattern(std::string & str)
{
	_boolPattern.reserve(str.size());
	for(unsigned k = 0; k < str.size(); ++k)
	{
		if(str[k] == '1')
			_matchPos.push_back(k);
		_boolPattern.push_back(str[k] == '1');
	}
}

const uint8_t & Pattern::operator[] (size_t idx) const 
{
	return _matchPos[idx];
}

bool Pattern::isMatch(const int & pos) const 
{
	return _boolPattern[pos];
}

size_t Pattern::size() const 
{ 
	return _boolPattern.size();
}