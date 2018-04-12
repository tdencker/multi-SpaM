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
#ifndef MATCHFINDER_HPP_
#define MATCHFINDER_HPP_

#include <vector>
#include "component.hpp"
#include "options.hpp"
#include "pattern.hpp"
#include "stats.hpp"
#include "word.hpp"
//#include "wordarray.hpp"

#include <iostream>

class MatchFinder
{
	public:
	MatchFinder(std::vector<Word> &, int, int);
	bool next();
	size_t size() { return std::distance(_start,_end); }
	size_t progress();
	std::vector<Component> getCurrentComponents(const Pattern &, int);
	uint64_t getCurrentKey();
	private:
    typedef std::vector<Word> ::iterator iterator;
	iterator _vec_start;
    iterator _start;
    iterator _end;
    iterator _vec_end;
	size_t last_progress;
};

#endif
