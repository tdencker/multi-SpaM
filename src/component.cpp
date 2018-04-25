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
#include "component.hpp"

Component::Component( const std::vector<Word>::iterator & first, int n ) : m_total( n, 0 ), m_component( this )
{
    m_words.push_back( first );
    m_total[first->getSeq()]++;
}

/**
* @brief Removes all words that are from sequences which have more than 1 word
* in this instance.
**/

void Component::removeUncertainties()
{
    for ( auto it = m_words.begin(); it != m_words.end(); )
    {
        auto & w = *it;
        if ( m_total[w->getSeq()] > 1 || m_total[w->getSeq()] == ambigious )
        {
            m_words.erase( it );
            m_total[w->getSeq()] = ambigious;
#pragma omp atomic
            mspamstats::ambigious_sequences++;
        }
        else
        {
            ++it;
        }
    }
}

unsigned Component::countSequences() const
{
    return std::count_if( m_total.begin(), m_total.end(), []( int x ) { return x > 0; } );
}

size_t Component::size() const
{
    return m_words.size();
}

Component::component_iter Component::begin()
{
    return m_words.begin();
}

Component::component_iter Component::end()
{
    return m_words.end();
}

void Component::erase( component_iter pos )
{
    m_words.erase( pos );
}

/**
* @brief Find function for the union find structure.
**/

Component & Component::getComponent()
{
    while ( m_component != m_component->m_component )
    {
        m_component = m_component->m_component;
    }
    return *m_component;
}

/**
* @brief Union function for the union find structure.
**/

void Component::merge( Component & other )
{
    assert( other.m_words.empty() == false );
    assert( this->m_words.empty() == false );
    // already in the same component
    if ( &other == this )
        return;
    // keep the larger component
    Component & merging_from = this->m_words.size() >= other.m_words.size() ? other : *this;
    Component & merging_into = this->m_words.size() >= other.m_words.size() ? *this : other;
    merging_into.m_words.insert( merging_into.m_words.end(), merging_from.m_words.begin(), merging_from.m_words.end() );
    for ( unsigned i = 0; i < m_total.size(); ++i )
        merging_into.m_total[i] += merging_from.m_total[i];
    merging_from.m_component = merging_into.m_component;
    merging_from.m_words.clear();
}
