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
    m_vec_start = sorted_array.begin() + x * sorted_array.size();
    if ( thread_id != 0 )
    {
        auto key = m_vec_start->getKey();
        while ( ++m_vec_start != sorted_array.end() && m_vec_start->getKey() == key )
            ;
    }
    x = (float) ( thread_id + 1 ) / thread_num;
    m_vec_end = sorted_array.begin() + x * sorted_array.size();
    if ( thread_id != thread_num - 1 )
    {
        auto key = m_vec_end->getKey();
        while ( ++m_vec_end != sorted_array.end() && m_vec_end->getKey() == key )
            ;
    }
    if(mspamoptions::use_seed == true)
    {
        m_gen = std::mt19937_64(mspamoptions::seed + thread_id);
    }else
    {
        std::random_device rd;
	    m_gen = std::mt19937_64(rd());
    }
	m_distr = std::uniform_int_distribution<>(0, std::distance(m_vec_start, m_vec_end) - 1);
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
    // set m_start and m_end accordingly
    // restart if there are not even 4 words
    int iterations = 0;
backup:
    bool repeat = true;
    int i;
    while ( repeat == true )
    {
        i = m_distr(m_gen);
        m_start = m_end = m_vec_start + i;
        if (m_start->isDummy() == true)
        {
            continue;
        }
        auto key = m_end->getKey();
        while ( ++m_end != m_vec_end && m_end->getKey() == key )
            ;
        while ( m_start != m_vec_start && ( --m_start )->getKey() == key )
            ;
        if ( m_start->getKey() != key )
            m_start++;
        repeat = std::count_if( m_start, m_end, [&]( Word & w ) { return w.isDummy() == false; } ) < 4;

#pragma omp atomic
        mspamstats::total_iterations++;
        if ( iterations++ > max_iterations )
            throw std::runtime_error( "maximum number of iterations reached." );
    }

    i -= std::distance(m_vec_start, m_start);

    std::vector<Component> components;
    int size = std::distance( m_start, m_end );
    components.reserve( size );
    for ( auto it = m_start; it != m_end; ++it )
    {
        components.emplace_back( it, nbr_sequences );
    }

    std::vector<int> shuf_vec(size, 0);
    std::iota(shuf_vec.begin(), shuf_vec.end(), 0);
    std::shuffle(shuf_vec.begin(), shuf_vec.end(), m_gen);

    for (int k = 0; k < size; ++k)
    {
        if (components[i].size() == mspamoptions::min_sequences)
        {
            break;
        }
        int j = shuf_vec[k];

        // Same sequences = match not possible
        // Same component = even if score is negative, it will be in the same pseudoalignment, so no need to compute
        if ((m_start + i)->getSeq() == (m_start + j)->getSeq()
            || &components[i].getComponent() == &components[j].getComponent())
        {
            continue;
        }

        // the word has previously been "removed"
        if ((m_start + j)->isDummy())
        {
            continue;
        }

        double score = 0;
        auto pos_seq1 = ( m_start + i )->getPos();
        auto pos_seq2 = ( m_start + j )->getPos();

        const int step_seq1 = ( m_start + i )->revComp() ? -1 : 1;
        const int step_seq2 = ( m_start + j )->revComp() ? -1 : 1;
        constexpr int alphabet_size = 4;
        for ( size_t l = 0; l < p.size();
                ++l, std::advance( pos_seq1, step_seq1 ), std::advance( pos_seq2, step_seq2 ) )
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

    components[i].removeUncertainties();

    // choose component and extract 4 words
    // then return QuartetBlock and "remove" the words from the array
    // if there are no components, "remove" all words and rerun the functions to
    // check another spaced word

    if ( components[i].size() < mspamoptions::min_sequences )
    {
        // return next(p, nbr_sequences);
        goto backup;
    }

    QuartetBlock qb( p.size() );
    Component & comp = components[i];

    for ( size_t j = 0; j < comp.size(); ++j )
    {
        qb.push_back( **( comp.begin() + j ) );
        **( comp.begin() + j ) = Word(); // "remove" to prevent realloc
    }
    return qb;
}