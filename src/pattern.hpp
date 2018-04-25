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

#include <string>
#include <vector>

/**
* @brief Class for the binary pattern. A `1` denotes a matching position and a
*`0` denotes
* a don't care positions. At these positions, mismatches are allowed.
**/

class Pattern
{
  public:
    Pattern( std::string & );
    Pattern( std::string && );
    const uint8_t & operator[]( size_t ) const;
    bool isMatch( const int & pos ) const;
    size_t size() const;

  private:
    std::vector<uint8_t> m_matchPos;
    std::vector<bool> m_boolPattern;
};

#endif
