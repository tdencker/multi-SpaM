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
 */#ifndef COMPONENT_HPP_
#define COMPONENT_HPP_

#include <algorithm>
#include <cassert>
#include <vector>
#include "options.hpp"
#include "stats.hpp"
#include "word.hpp"

class Component
{
    public:
    typedef std::vector< std::vector<Word>::iterator >::iterator component_iter;
	Component(const std::vector<Word>::iterator &, int);
	Component(const std::vector<Word>::iterator & first, int count, int n); // remove after testing
	unsigned countSequences() const;
	void merge(Component &);
	Component & getComponent();
	
	// testing
	
	void setComponent(Component & comp)
	{
	    _component = &comp;
	}

	void erase(component_iter first, component_iter last)
	{
		assert(first >= _words.begin() && last <= _words.end());
		_words.erase(first, last);
	}
	
	
	void removeUncertainties();
	size_t size() const;
	component_iter begin();
	component_iter end();
	void erase(component_iter);
	private:
	std::vector<unsigned> _total;
	std::vector< std::vector<Word>::iterator > _words;
	Component * _component;
	static constexpr int ambigious = 0;
};

#endif
