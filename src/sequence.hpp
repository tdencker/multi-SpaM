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

#ifndef SEQUENCE_HPP_
#define SEQUENCE_HPP_

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "stats.hpp"

/**
* @brief Sequence container class. ID is the header of the sequence as described
* in the FASTA file and CONTENT is the DNA sequence converted to 0-3.
**/

class Sequence
{
  public:
    std::vector<char> content;
    std::string id;
    static std::vector<Sequence> read( std::string &, bool = false );
    Sequence( std::string, std::string & );
};

#endif
