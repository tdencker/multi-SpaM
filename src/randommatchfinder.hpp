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

#ifndef RANDOMMATCHFINDER_HPP_
#define RANDOMMATCHFINDER_HPP_

//TODO:check
#include <vector>
#include "component.hpp"
#include "options.hpp"
#include "pattern.hpp"
#include "quartetblock.hpp"
#include "stats.hpp"
#include "word.hpp"

/**
* @brief This class is used to sample QuartetBlocks at random. It will divide the (sorted) vector
* in <number threads> parts. The constructor needs to be called by every thread. The next function
* will provide QuartetBlocks until it can't find more in which case it will throw an exception.
* Words that were used for QuartetBlocks are overwritten by "dummy" words.
**/

class RandomMatchFinder
{
	public:
	RandomMatchFinder(std::vector<Word> &, int, int);
	QuartetBlock next(const Pattern &, int);
	private:
	std::vector<Word>::iterator _vec_start;
	std::vector<Word>::iterator _start;
	std::vector<Word>::iterator _end;
	std::vector<Word>::iterator _vec_end;
	std::vector<char> _dummy_vec;
};

#endif