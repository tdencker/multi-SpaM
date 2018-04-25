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

#include "quartetblock.hpp"

/**
* @brief First way: provide the Word iterators to create this block. If you do
* not have a sequence key, you can just call setSequenceKey() afterwards and
* pass 0.
**/

QuartetBlock::QuartetBlock( std::vector<std::vector<Word>::iterator> & vec, int pattern_length, uint64_t seq_key )
    : m_length( pattern_length ), m_seq_key( seq_key )
{
    assert( vec.size() >= mspamoptions::min_sequences );
    assert( vec.size() <= 32 );
    m_words.reserve( vec.size() );
    for ( auto & it : vec )
        m_words.push_back( *it );
    std::sort( m_words.begin(), m_words.end(),
               []( const Word & a, const Word & b ) { return a.getSeq() < b.getSeq(); } );
}

QuartetBlock::QuartetBlock( std::vector<Word> && vec, int length ) : m_words( vec ), m_length( length )
{
    setSequenceKey();
}

/**
* @brief Second way: create an empty QuartetBlock and only pass the length
*(should be the pattern length)
* and then call push_back to add the words.
**/

QuartetBlock::QuartetBlock( int length ) : m_length( length ), m_seq_key( 0 )
{
}

void QuartetBlock::push_back( const Word & w )
{
    m_seq_key |= 1 << w.getSeq();
    m_words.push_back( w );
}

bool QuartetBlock::operator<( const QuartetBlock & other ) const
{
    return m_seq_key < other.m_seq_key;
}

// convenience functions

size_t QuartetBlock::size() const
{
    return m_words.size();
}

std::string QuartetBlock::getHeader( int idx, std::vector<Sequence> & sequences )
{
    return sequences[m_words[idx].getSeq()].id;
}

std::string QuartetBlock::to_string( std::vector<Sequence> & sequences )
{
    std::string str;
    str += "##################### Length: " + std::to_string( m_length ) + " ######################\n";
    for ( auto & w : m_words )
    {
        str += std::to_string( w.getSeq() ) + " : " +
               std::to_string( std::distance( sequences[w.getSeq()].content.begin(), w.getPos() ) ) + "\n";
    }
    return str;
}

void QuartetBlock::setSequenceKey()
{
    m_seq_key = 0;
    for ( auto & e : m_words )
    {
        m_seq_key |= 1 << e.getSeq();
    }
}

uint64_t QuartetBlock::getSequenceKey() const
{
    return m_seq_key;
}

std::vector<char>::iterator QuartetBlock::getPos( int idx ) const
{
    return m_words[idx].getPos();
}

int QuartetBlock::getLength() const
{
    return m_length;
}

int QuartetBlock::getSeq( int idx ) const
{
    return m_words[idx].getSeq();
}

bool QuartetBlock::getRevComp( int idx ) const
{
    return m_words[idx].revComp();
}

void QuartetBlock::clear()
{
    m_words.clear();
}

QuartetBlock::iterator QuartetBlock::begin()
{
    return m_words.begin();
}

QuartetBlock::iterator QuartetBlock::end()
{
    return m_words.end();
}