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

#ifndef PATTERN_HPP_
#define PATTERN_HPP_

#include <vector>
#include <string>

class Pattern
{
	public:
	unsigned weight; // kinda superfluous
	unsigned dontCare;
	
	Pattern(std::string && str)
	{
//		_matchPos.reserve(str.size());
		_boolPattern.reserve(str.size());
		for(unsigned k = 0; k < str.size(); ++k)
		{
			if(str[k] == '1')
				_matchPos.push_back(k);
			_boolPattern.push_back(str[k] == '1');
		}
		weight = _matchPos.size();
		dontCare = _boolPattern.size() - weight;
	}
	
	Pattern(std::string & str)
	{
//		_matchPos.reserve(str.size());
		_boolPattern.reserve(str.size());
		for(unsigned k = 0; k < str.size(); ++k)
		{
			if(str[k] == '1')
				_matchPos.push_back(k);
			_boolPattern.push_back(str[k] == '1');
		}
		weight = _matchPos.size();
		dontCare = _boolPattern.size() - weight;
	}
	
	const uint8_t & operator[] (size_t idx) const { return _matchPos[idx];}
	
	bool isMatch(const int & pos) const { return _boolPattern[pos];}
	
	int size() const { return weight + dontCare; }
	
	private:
	std::vector<uint8_t> _matchPos;
	std::vector<bool> _boolPattern;
};

#endif
