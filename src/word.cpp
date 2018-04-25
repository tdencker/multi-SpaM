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

/**
* @file
* @brief Definitions of the Word class.
**/

#include "word.hpp"

/**
* @brief Computes the reverse nucleotide (encoded as 2 bit)
**/

static inline char reverse( char & c )
{
    return ~c & 3;
}

/**
* @brief Vector to mark words as "removed"
**/

std::vector<char> Word::dummy_vec;

/**
* @brief String representation of the Word (also accounting for reverse
*complement)
**/

std::string Word::toString( unsigned length )
{
    std::string buf;
    buf.reserve( length );
    if ( m_rev_comp )
    {
        for ( unsigned i = 0; i < length; ++i )
        {
            buf[i] = reverse( *( m_pos - i ) );
        }
    }
    else
    {
        for ( unsigned i = 0; i < length; ++i )
        {
            buf[i] = *( m_pos + i );
        }
    }
    for ( unsigned i = 0; i < length; ++i )
    {
        switch ( buf[i] )
        {
        case 0:
            buf[i] = 'A';
            break;
        case 1:
            buf[i] = 'C';
            break;
        case 2:
            buf[i] = 'G';
            break;
        case 3:
            buf[i] = 'T';
            break;
        default:
            std::cerr << "Illegal character: " << buf[i] << std::endl;
        }
    }
    return buf;
}

/**
* @brief Constructor of Word: Encodes the nucleotides at the match positions
* as indicated by the pattern into a 32 bit unsigned int (2 bit each). Throws if
* there is a non-nucleotide symbol at one of the matching positions.
*
* @param p - Pattern (e.g. 10001011011) indicating match and don't care
*positions
* @param it - Starting position of the spaced word
* @param seq - Id of the sequence of the word
* @param rev_comp - Reverse complement word?
**/

Word::Word( Pattern & p, std::vector<char>::iterator it, int16_t seq,
            bool rev_comp )
    : m_pos( it ), m_key( 0 ), m_seq( seq ), m_rev_comp( rev_comp )
{
    const unsigned mask = mspamoptions::mask;
    const unsigned shift = mspamoptions::symbol_bits;
    for ( int i = 0; i < mspamoptions::weight; i++ )
    {
        unsigned nuc = 0;
        int offset = p[i];
        if ( m_rev_comp )
            offset *= -1;

        nuc = *( it + offset );

        if ( nuc == ( std::numeric_limits<char>::max )() )
        {
            throw std::exception(); // illegal character
        }

        if ( m_rev_comp )
            nuc = mask - nuc;
        m_key <<= shift;
        m_key |= nuc;
    }
}

/**
* @brief Constructor of a "removed" word. Only useful for isDummy().
**/

Word::Word()
    : m_pos( dummy_vec.begin() ), m_key( 0 ), m_seq( 0 ), m_rev_comp( false )
{
}

char Word::operator[]( const int & idx ) const
{
    return m_rev_comp ? reverse( *( m_pos - idx ) ) : *( m_pos + idx );
}

bool Word::revComp() const
{
    return m_rev_comp;
}

bool Word::isDummy() const
{
    return m_pos == dummy_vec.begin();
}

bool Word::operator<( const Word & other ) const
{
    return m_key < other.m_key;
}

bool Word::operator==( const Word & other ) const
{
    return m_pos == other.m_pos;
}

int Word::getSeq() const
{
    return m_seq;
}

uint32_t Word::getKey() const
{
    return m_key;
}

std::vector<char>::iterator Word::getPos() const
{
    return m_pos;
}