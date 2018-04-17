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

Component::Component(const std::vector<Word>::iterator & first, int n)
    : _total(n,0),
    _component(this)
{
    _words.push_back(first);
    _total[first->getSeq()]++;
}

Component::Component(const std::vector<Word>::iterator & first, int count, int n)
    : _total(n,0),
    _component(this)
{
    auto it = first;
    for(int i = 0; i < count; ++i)
    {
        _words.push_back(it);
        _total[it->getSeq()]++;
        it++;
    }
}

void Component::removeUncertainties()
{
    for(auto it = _words.begin(); it != _words.end();)
    {
        auto & w = *it;
        if(_total[w->getSeq()] > 1 || _total[w->getSeq()] == ambigious)
        {
            _words.erase(it);
            _total[w->getSeq()] = ambigious;
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
	return std::count_if(_total.begin(), _total.end(), [](int x){return x > 0;});
}

size_t Component::size() const
{
    return _words.size();
}

Component::component_iter Component::begin()
{
    return _words.begin();
}

Component::component_iter Component::end()
{
    return _words.end();
}

void Component::erase(component_iter pos)
{
    _words.erase(pos);
}

Component & Component::getComponent()
{
	while( _component != _component->_component)
	{
		_component = _component->_component;
	} 
	return *_component; 
}

void Component::merge(Component & other)
{
    assert(other._words.empty() == false);
    assert(this->_words.empty() == false);
    // already in the same component
	if(&other == this)
		return;
	// keep the larger component
	Component & merging_from = this->_words.size() >= other._words.size() ? other : *this;
	Component & merging_into = this->_words.size() >= other._words.size() ? *this : other;
	merging_into._words.insert(merging_into._words.end(), merging_from._words.begin(), merging_from._words.end());
	for(unsigned i = 0; i < _total.size(); ++i)
		merging_into._total[i] += merging_from._total[i];
	merging_from._component = merging_into._component;
	merging_from._words.clear();
}
