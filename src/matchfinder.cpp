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
 */#include "matchfinder.hpp"

size_t MatchFinder::progress()
{
	auto dist = std::distance(_vec_start, _end);
	auto new_progress = dist - last_progress;
	last_progress = dist;
	return new_progress;
}

MatchFinder::MatchFinder(std::vector<Word> & sorted_array, int thread_id, int thread_num) : last_progress(0)
{
    float x = (float) thread_id / thread_num;
//    #pragma omp critical
//    std::cout << thread_id << " " << thread_num << " " << x << std::endl;
	_end = sorted_array.begin() + x * sorted_array.size();
	if(thread_id != 0)
	{
		auto key = _end->getKey();
		while( ++_end != sorted_array.end() && _end->getKey() == key)
			;
	}
	_vec_start = _end;
    x = (float) (thread_id + 1) / thread_num;
//    #pragma omp critical
//    std::cout << thread_id << " " << thread_num << " " << x << std::endl;
	_vec_end = sorted_array.begin() + x * sorted_array.size();
	if(thread_id != thread_num - 1)
	{
		auto key = _vec_end->getKey();
		while( ++_vec_end != sorted_array.end() && _vec_end->getKey() == key)
			;
	}
//	#pragma omp critical
//	std::cout << thread_id << ": " << std::distance(sorted_array.begin(), _vec_start) << " " << std::distance(sorted_array.begin(), _vec_end) << " " << sorted_array.size() << std::endl;
}

bool MatchFinder::next()
{
	while(_end != _vec_end)
	{
		_start = _end;
		auto current_key = _start->getKey();
		while(++_end != _vec_end && _end->getKey() == current_key)
			;
		if(std::distance(_start, _end) >= options::min_sequences)
		{
			return true;
		}
	}
	return false;
}

uint64_t MatchFinder::getCurrentKey()
{
    return _start->getKey();
}

int8_t score_vec[16] = { 91, -114, -31, -123, -114, 100, -125, -31, -31, -125, 100, -114, -123, -31, -114, 91};

// TODO: compute all match scores again?
std::vector<Component> MatchFinder::getCurrentComponents(const Pattern & p, int nbr_sequences)
{
	std::vector<Component> components;
	int size = std::distance(_start, _end);
	components.reserve(size);
	for(auto it = _start; it != _end; ++it)
	{
		components.emplace_back(it, nbr_sequences);
	}
	
	for(int i = 0; i < size; ++i)
	{
	    for(int j = i + 1; j < size; ++j)
	    {
	        // Same sequences = match not possible
	        // Same component = even if score is negative, it will be in the same pseudoalignment, so no need to compute
	        if( (_start + i)->getSeq() == (_start + j)->getSeq() 
	            || &components[i].getComponent() == &components[j].getComponent() )
	        {
	            continue;
	        }
            
	        if((_start + i)->isDummy() || (_start + j)->isDummy())
	        {
	        	continue;
	        }
	        
	        #pragma omp atomic
	        stats::score_computations++;
	        double score = 0;
	        auto pos_seq1 = (_start + i)->getPos();
	        auto pos_seq2 = (_start + j)->getPos();
	        const int step_seq1 = (_start + i)->revComp() ? -1 : 1;
	        const int step_seq2 = (_start + j)->revComp() ? -1 : 1;
	        constexpr int alphabet_size = 4;
	        for(int k = 0; k < p.size(); ++k, std::advance(pos_seq1, step_seq1), std::advance(pos_seq2, step_seq2))
	        {
//	            if(p.isMatch(k))
//	            {
//	                continue;
//	            }
	            int idx = *pos_seq1 * alphabet_size + *pos_seq2;
	            score += score_vec[ idx ];
	        }
	        
	        if(score >= options::min_score)
	        {
	            components[i].getComponent().merge(components[j].getComponent());
	        }
	    }
	}
	int cnt = 0;
	for(auto & e : components)
	{
	    e.removeUncertainties();
		if(e.size() >= options::min_sequences)
			cnt++;
	}
	if(cnt > 1)
	    #pragma omp atomic
		stats::multiple_components++;
		
	components.erase(std::remove_if(components.begin(), components.end(), [](Component & c){ return c.size() < options::min_sequences; }), components.end());
	return components;
}
