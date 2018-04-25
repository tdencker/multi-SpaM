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

Pattern::Pattern( std::string && str )
{
    m_boolPattern.reserve( str.size() );
    for ( unsigned k = 0; k < str.size(); ++k )
    {
        if ( str[k] == '1' )
            m_matchPos.push_back( k );
        m_boolPattern.push_back( str[k] == '1' );
    }
}

Pattern::Pattern( std::string & str )
{
    m_boolPattern.reserve( str.size() );
    for ( unsigned k = 0; k < str.size(); ++k )
    {
        if ( str[k] == '1' )
            m_matchPos.push_back( k );
        m_boolPattern.push_back( str[k] == '1' );
    }
}

/**
* @brief Returns the idx-th match position, i.e. `1` in the pattern.
* This function does not check wether idx is larger than the weight
* of the pattern, so be careful because idx is not the idx-th position
* in the pattern.
**/

const uint8_t & Pattern::operator[]( size_t idx ) const
{
    return m_matchPos[idx];
}

/**
* @brief Checks wether the pos-th position is a match position, i.e. a `1`
* in the binary pattern.
**/

bool Pattern::isMatch( const int & pos ) const
{
    return m_boolPattern[pos];
}

/**
* @brief Returns the size of the pattern which should be #match positions +
* #don't care positions.
**/

size_t Pattern::size() const
{
    return m_boolPattern.size();
}