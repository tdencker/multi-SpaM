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

#include "randommatchfinder.hpp"

/**
* @brief The vector has to be sorted. It will divide the range according to the
* thread id/num.
**/

RandomMatchFinder::RandomMatchFinder( std::vector<Word> & sorted_array, int thread_id, int thread_num )
{
    float x = (float) thread_id / thread_num;
    _vec_start = sorted_array.begin() + x * sorted_array.size();
    if ( thread_id != 0 )
    {
        auto key = _vec_start->getKey();
        while ( ++_vec_start != sorted_array.end() && _vec_start->getKey() == key )
            ;
    }
    x = (float) ( thread_id + 1 ) / thread_num;
    _vec_end = sorted_array.begin() + x * sorted_array.size();
    if ( thread_id != thread_num - 1 )
    {
        auto key = _vec_end->getKey();
        while ( ++_vec_end != sorted_array.end() && _vec_end->getKey() == key )
            ;
    }
}

int8_t score_mat[16] = {91, -114, -31, -123, -114, 100, -125, -31, -31, -125, 100, -114, -123, -31, -114, 91};
constexpr int max_iterations = 10000;

/**
* @brief This function selects a Word at random and checks for the range of
*Words that have the same hash/key.
* Among those words, random words are selected and a score is computed. If there
*are 3 words with a positive
* score, the QuartetBlock is returned. The words used are replaced with a
*"dummy" word.
*
* It will throw if max_iterations iterations have happened, i.e. either a Dummy
*word has been hit, there are
* less than 4 words left with that key or there were not enough words in 4
*sequences that have a positive score
* with the initial word.
**/

QuartetBlock RandomMatchFinder::next( const Pattern & p, int nbr_sequences )
{
    // choose a spaced word in the array
    // set _start and _end accordingly
    // restart if there are not even 4 words
    int iterations = 0;
backup:
    bool repeat = true;
    while ( repeat == true )
    {
        int x = rand() % std::distance( _vec_start, _vec_end );
        _start = _end = _vec_start + x;
        auto key = _end->getKey();
        while ( ++_end != _vec_end && _end->getKey() == key )
            ;
        while ( _start != _vec_start && ( --_start )->getKey() == key )
            ;
        if ( _start->getKey() != key )
            _start++;
        if ( std::count_if( _start, _end, [&]( Word & w ) { return w.getPos() != _dummy_vec.end(); } ) >= 4 )
        {
            repeat = false;
        }
#pragma omp atomic
        mspamstats::total_iterations++;
        if ( iterations++ > max_iterations )
            throw std::runtime_error( "maximum number of iterations reached." );
    }

    std::vector<Component> components;
    int size = std::distance( _start, _end );
    components.reserve( size );
    for ( auto it = _start; it != _end; ++it )
    {
        components.emplace_back( it, nbr_sequences );
    }

    for ( int i = 0; i < size; ++i )
    {
        for ( int j = i + 1; j < size; ++j )
        {
            // Same sequences = match not possible
            // Same component = even if score is negative, it will be in the
            // same QuartetBlock, so no need to compute
            if ( ( _start + i )->getSeq() == ( _start + j )->getSeq() ||
                 &components[i].getComponent() == &components[j].getComponent() )
            {
                continue;
            }

            double score = 0;
            auto pos_seq1 = ( _start + i )->getPos();
            auto pos_seq2 = ( _start + j )->getPos();

            // either word has previously been "removed"
            if ( pos_seq1 == _dummy_vec.end() || pos_seq2 == _dummy_vec.end() )
            {
                continue;
            }
            const int step_seq1 = ( _start + i )->revComp() ? -1 : 1;
            const int step_seq2 = ( _start + j )->revComp() ? -1 : 1;
            constexpr int alphabet_size = 4;
            for ( size_t k = 0; k < p.size();
                  ++k, std::advance( pos_seq1, step_seq1 ), std::advance( pos_seq2, step_seq2 ) )
            {
                //	            if(p.isMatch(k))
                //	            {
                //	                continue;
                //	            }
                int idx = *pos_seq1 * alphabet_size + *pos_seq2;
                score += score_mat[idx];
            }
            if ( score >= mspamoptions::min_score )
            {
                components[i].getComponent().merge( components[j].getComponent() );
            }
            else
            {
#pragma omp atomic
                mspamstats::random_matches++;
            }
        }
    }
    int cnt = 0;
    for ( auto & e : components )
    {
        e.removeUncertainties();
        if ( e.size() >= mspamoptions::min_sequences )
            cnt++;
    }

    components.erase( std::remove_if( components.begin(), components.end(),
                                      []( Component & c ) { return c.size() < mspamoptions::min_sequences; } ),
                      components.end() );

    // choose component and extract 4 words
    // then return QuartetBlock and "remove" the words from the array
    // if there are no components, "remove" all words and rerun the functions to
    // check another spaced word

    if ( components.size() == 0 )
    {
        // return next(p, nbr_sequences);
        goto backup;
    }

    QuartetBlock qb( p.size() );
    int idx = rand() % components.size();
    Component & comp = components[idx];
    idx = rand() % ( comp.size() - 3 );
    for ( int i = 0; i < 4; ++i )
    {
        qb.push_back( **( comp.begin() + idx + i ) );
        **( comp.begin() + idx + i ) = Word(); // "remove" to prevent realloc
    }
    return qb;
}