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
#ifndef COMPONENT_HPP_
#define COMPONENT_HPP_

#include <algorithm>
#include <cassert>
#include <vector>

#include "options.hpp"
#include "stats.hpp"
#include "word.hpp"

/**
* @brief Mostly a legacy class. Mostly just a QuartetBlock light. Used in the
* RandomMatchFinder class. It is a union find datastructure but no longer used
* that way.
**/

class Component
{
  public:
    typedef std::vector<std::vector<Word>::iterator>::iterator component_iter;
    Component( const std::vector<Word>::iterator &, int );
    unsigned countSequences() const;
    void merge( Component & );
    Component & getComponent();
    void removeUncertainties();
    size_t size() const;
    component_iter begin();
    component_iter end();
    void erase( component_iter );

  private:
    std::vector<unsigned> m_total;
    std::vector<std::vector<Word>::iterator> m_words;
    Component * m_component;
    static constexpr int ambigious = 0;
};

#endif
