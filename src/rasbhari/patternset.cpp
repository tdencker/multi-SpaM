/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * patternset object file
 *
 * For theory please have a look at:
 *
 * - L. Hahn, C.-A. Leimeister, B. Morgenstern (2015)
 * RasBhari: optimizing spaced seeds for database searching, read mapping and
 * alignment-free sequence comparison
 * arXiv:1511.04001[q-bio.GN]
 *
 * - B. Morgenstern, B. Zhu, S. Horwege, C.-A Leimeister (2015)
 * Estimating evolutionary distances between genomic sequences from spaced-word
 * matches
 * Algorithms for Molecular Biology 10, 5.
 * (http://www.almob.org/content/10/1/5/abstract)
 *
 *
 * @author: Lars Hahn - 12.05.2016, Georg-August-Universitaet Goettingen
 * @version: 1.2.0 05/2016
 */
#include "patternset.h"

/*---Variables---------------------------------------------------------------*/
std::default_random_engine generator( std::random_device{}() ); // time(0)

/*===Main-Part===============================================================*/
/*---Constructor-------------------------------------------------------------*/
/**
 * Default constructor, sets the default vaulues.
 */
patternset::patternset()
{
    size = 10;
    weight = 8;
    max_dontcare = 6;
    min_dontcare = 6;
    silent = true;
    isInFile = false;
    bitmode = true;
    random_length = false;
    improve = true;
    update = false;
    initialized = false;
}

/**
 * Constructor, sets some default vaulues, just pattern dimension is set.
 *
 * @param size			Number of patterns for the patternset.
 *
 * @param weigth		Number of match positions ('1') in each pattern.
 *
 * @param min_dontcare	Number of minimum don't care positions ('0') in each
 * pattern.
 *
 * @param max_dontcare	Number of maximum don't care positions ('0') in each
 * pattern.
 *
 * @param silent 		Boolean which (de)activates std::output on
 * commandline (exclusive errors)
 */
patternset::patternset( uint32_t size, uint32_t weight, uint32_t min_dontcare, uint32_t max_dontcare, bool silent )
{
    this->size = size;
    this->weight = weight;
    this->max_dontcare = max_dontcare;
    this->min_dontcare = min_dontcare;
    this->silent = silent;
    isInFile = false;
    bitmode = true;
    random_length = false;
    improve = true;
    update = false;
    size = 10;
    weight = 8;
    max_dontcare = 6;
    min_dontcare = 6;
    initialized = false;
}

/**
 * File constructor, sets values only from files.
 *
 * @param inFile		The file which contains a patternset.
 *
 * @param silent 		Boolean which (de)activates std::output on
 * commandline (exclusive errors)
 */
patternset::patternset( const char * inFile, bool silent )
{
    this->silent = silent;
    this->inFile = inFile;
    isInFile = true;
    improve = true;
    update = false;
    bitmode = true;
    size = 10;
    weight = 8;
    max_dontcare = 6;
    min_dontcare = 6;
    initialized = false;
}

/**
 * Default destructor, deletes all vectors and matrices in the object.
 */
patternset::~patternset()
{
    pattern_set.clear();
    bit_pattern_set.clear();
    pattern_dontcare.clear();
    pattern_improve.clear();
}

/*---Initialization----------------------------------------------------------*/
/**
 * Creates for the submitted or default values a set of pattern; if errors
 * occur, these will be fixed.
 * It is possible to get all patterns for a specific combination.
 * Size and lenght can be reduced; each pattern is unique.
 *
 * @param random_length	Boolean which (de)activates random generated pattern
 * lengths.
 */
void patternset::Init( bool random_length )
{
    bool read_error;

    read_error = false;
    this->random_length = random_length;

    if ( isInFile )
    {
        read_error = ReadPattern();
        if ( read_error )
        {
            SecureMessage( "nkl", -1 );
            pattern_set.clear();
            update = true;
        }
        else
        {
            CheckParams();
            improve = SetImprove();
            if ( bitmode )
            {
                PatternToBit();
            }
        }
    }

    if ( !isInFile || read_error )
    {
        if ( weight != 1 )
        {
            CheckParams();
            pattern_dontcare = CreateLengths();
            pattern_set = CreatePattern();
            improve = SetImprove();
            if ( bitmode )
            {
                PatternToBit();
            }
        }
        else
        {
            pattern_dontcare.push_back( 0 );
            pattern_set = CreatePattern();
            improve = false;
            pattern_improve.push_back( improve );
            PatternToBit();
        }
    }
    initialized = true;
}

/**
 * Standard Reinitialization of variance options, creates a new patternset.
 */
void patternset::ReInit()
{
    pattern_set.clear();
    bit_pattern_set.clear();
    pattern_dontcare.clear();
    pattern_improve.clear();
    bitmode = true;
    improve = true;
    update = false;
    initialized = false;
    Init( random_length );
}

/**
 * Checks the first conditions of the patternset, in special cases resets them.
 * Checks, if bitmode is possible.
 */
void patternset::CheckParams()
{
    // if(size < 0){
    //	size = std::abs(size);
    //	if(size <= 0){
    //		size = 10;
    //		update = true;
    //	}
    //}
    // if(max_dontcare < 0 || min_dontcare < 0){
    //	max_dontcare = std::abs(max_dontcare);
    //	min_dontcare = std::abs(min_dontcare);
    //}
    // if(max_dontcare < min_dontcare){
    //	std::swap(max_dontcare,min_dontcare);
    //}
    if ( weight < 2 )
    {
        if ( weight < 0 )
        {
            // weight = std::abs(weight);
        }
        else
        {
            if ( size != 1 )
            {
                size = 1;
                update = true;
            }
            if ( max_dontcare != 0 || max_dontcare != 0 )
            {
                min_dontcare = 0;
                max_dontcare = 0;
                update = true;
            }
        }
    }
    if ( weight == 2 )
    {
        improve = false;
    }
    if ( max_dontcare + weight > 63 || min_dontcare + weight > 63 )
    {
        bitmode = false;
    }
    if ( update && !silent )
    {
        SecureMessage( "update", -1 );
    }
}

/*---Functions---------------------------------------------------------------*/
/*---Func-Create-------------------------------------------------------------*/
/**
 * Generates automatically patternlengths; depends on min_dontcare and
 *max_dontcare.
 * Returns a vector which contains in descending order the pattern lenghts.
 * Two options: 1) nearly linear distributed lenghts
 *				2) random distributed lenghts
 *
 * @return 				All lenghts for the patternset (type:
 *std::vector<uint32_t>)
 */
std::vector<uint32_t> patternset::CreateLengths()
{
    std::vector<uint32_t> leng_vec;
    double steps, leng, leng_tmp, max_size;
    uint32_t intervals;

    if ( random_length )
    {
        std::uniform_int_distribution<uint32_t> lengths( min_dontcare, max_dontcare );
        for ( uint32_t i = 0; i < size; i++ )
        {
            leng_vec.push_back( lengths( generator ) );
        }
    }
    else
    {
        intervals = max_dontcare - min_dontcare + 1;
        steps = ( (double) intervals / (double) size );
        leng = (double) min_dontcare;
        for ( uint32_t i = 0; i < size; i++ )
        {
            leng_vec.push_back( (uint32_t) leng );
            leng += steps;
        }
    }

    std::sort( leng_vec.begin(), leng_vec.end() );

    leng_tmp = leng_vec[0];
    leng = 0;

    max_size = MaxNumberPattern( weight - 2, leng_tmp + weight - 2 );

    for ( uint32_t i = 0; i < size; i++ )
    {
        if ( leng_vec[i] <= leng_tmp )
        {
            if ( leng_vec[i] < leng_tmp )
            {
                leng_vec[i] = leng_tmp;
            }
            if ( (double) leng >= max_size )
            {
                if ( leng_tmp != leng_vec[leng_vec.size() - 1] )
                {
                    leng = 0;
                    leng_vec[i]++;
                    leng_tmp = leng_vec[i];
                    max_size = MaxNumberPattern( weight - 2, leng_tmp + weight - 2 );
                }
                else
                {
                    update = true;
                    leng_vec.erase( leng_vec.begin() + i, leng_vec.end() );
                    size = leng_vec.size();
                    SecureMessage( "max_number_pattern", i + 1 );
                }
            }
        }
        else
        {
            leng_tmp = leng_vec[i];
            leng = 0;
            max_size = MaxNumberPattern( weight - 2, leng_tmp + weight - 2 );
        }
        leng++;
    }
    std::sort( leng_vec.rbegin(), leng_vec.rend() );

    return leng_vec;
}

/**
 * Generates for submitted (and corrected) values a patternset as std::vector<
 * std::vector<char> >.
 *
 * @return 				Patternset (type: std::vector<
 * std::vector<char>
 * >)
 */
std::vector<std::vector<char>> patternset::CreatePattern()
{
    std::vector<std::vector<char>> patterns;
    std::vector<char> pat;
    uint32_t tmp_weight, pos;
    if ( !silent )
    {
        std::cout << "Generating pattern automatically ...\n" << std::endl;
    }

    if ( weight != 1 )
    {
        for ( uint32_t i = 0; i < size; i++ )
        {
            for ( uint32_t j = 0; j < pattern_dontcare[i] + weight; j++ )
            {
                pat.push_back( '0' );
            }
            pat[0] = '1';
            pat[pat.size() - 1] = '1';
            tmp_weight = weight - 2;
            std::uniform_int_distribution<uint32_t> positions( 0, pattern_dontcare[i] + weight - 1 );
            while ( tmp_weight > 0 )
            {
                pos = positions( generator );
                if ( pat[pos] == '0' )
                {
                    tmp_weight--;
                    pat[pos] = '1';
                }
            }
            if ( UniqPattern( patterns, pat ) )
            {
                patterns.push_back( pat );
                if ( !silent )
                {
                    for ( uint32_t i = 0; i < pat.size(); i++ )
                    {
                        std::cout << pat[i];
                    }
                    std::cout << std::endl;
                }
            }
            else
            {
                i--;
            }
            pat.clear();
        }
    }
    else
    {
        pat.push_back( '1' );
        patterns.push_back( pat );
        pat.clear();
    }

    if ( !silent )
    {
        std::cout << "\nDone." << std::endl << std::endl;
    }

    return patterns;
}

/**
 * Reads pattern from a submitted file which (may) contains binary strings.
 * Each pattern will be verified, if its really a binary pattern with the same
 * weight.
 * If an error occurs (file not found, etc.), it returns false, otherwise true.
 *
 * @return 				Boolean, if an error occured or not
 * while
 * pattern
 * processing.
 */
bool patternset::ReadPattern()
{
    std::ifstream fasta;
    std::vector<char> pattern;
    uint32_t p_size;
    char tokens[6] = {'.', ' ', ',', ';', '\n', '\t'};
    char c;
    bool token, read_error, bad;

    read_error = false;
    bad = false;
    fasta.open( inFile );
    if ( fasta.fail() )
    {
        read_error = true;
        SecureMessage( "file", -1 );
    }
    else if ( fasta.peek() == std::ifstream::traits_type::eof() )
    {
        read_error = true;
        SecureMessage( "empty", -1 );
    }
    else
    {
        while ( !fasta.eof() )
        {
            c = fasta.get();
            token = false;
            for ( uint32_t i = 0; i < 6; i++ )
            {
                if ( c == tokens[i] )
                {
                    token = true;
                }
            }
            if ( token || -1 == (int) c )
            {
                if ( !bad && pattern.size() != 0 )
                {
                    pattern_set.push_back( pattern );
                }
                bad = false;
                pattern.clear();
            }
            else if ( c == '1' || c == '0' )
            {
                pattern.push_back( c );
            }
            else
            {
                bad = true;
                SecureMessage( "format", -1 );
            }
        }
    }
    fasta.close();
    size = (uint32_t) pattern_set.size();
    if ( size == 0 )
    {
        if ( !read_error )
        {
            SecureMessage( "empty", -1 );
        }
        read_error = true;
        size = 10;
    }

    for ( uint32_t i = 0; i < pattern_set.size(); i++ )
    {
        p_size = pattern_set[i].size();
        if ( p_size != 0 )
        {
            p_size--;
        }
        if ( pattern_set[i][0] != '1' || pattern_set[i][p_size] != '1' )
        {
            SecureMessage( "startend", i );
            read_error = true;
        }
    }
    if ( !read_error )
    {
        read_error = VerifyPatternCondition();
        if ( !read_error )
        {
            for ( uint32_t i = 0; i < size; i++ )
            {
                pattern_dontcare.push_back( pattern_set[i].size() - weight );
            }
        }
    }
    return read_error;
}

/**
* Checks, if the submitted pattern file is in correct format; if it
* contains only binary pattern with the same weight.
* Sets the wight and max_ and min_dontcare number, if no error occured.
*/
bool patternset::VerifyPatternCondition()
{
    uint32_t first_weight, p_size;
    bool pattern_error;

    pattern_error = false;

    min_dontcare = pattern_set[0].size();
    max_dontcare = min_dontcare;

    first_weight = GetWeight( pattern_set[0] );

    for ( uint32_t i = 0; i < size; i++ )
    {
        p_size = pattern_set[i].size();
        min_dontcare = std::min( min_dontcare, p_size );
        max_dontcare = std::max( max_dontcare, p_size );
        if ( first_weight != GetWeight( pattern_set[i] ) )
        {
            pattern_error = true;
            SecureMessage( "weight", i );
        }
        if ( !IsBinary( pattern_set[i] ) )
        {
            pattern_error = true;
            SecureMessage( "pattern", i );
        }
    }

    weight = first_weight;

    min_dontcare -= weight;
    max_dontcare -= weight;

    if ( min_dontcare < 0 || max_dontcare < 0 )
    {
        pattern_error = true;
    }

    return pattern_error;
}

/**
 * Checks, if there is at least one pattern, which can be improved.
 * Returns true, if there is at least one, otherwise false.
 *
 * Checks and saves for each pattern individually, if it can be improved.
 *
  * @return 			Boolean, if at least one pattern can be changed.
 */
bool patternset::SetImprove()
{
    uint32_t leng, p_size, offset;
    bool isImprovable, tmp_impro;

    isImprovable = false;

    leng = pattern_dontcare[0];
    p_size = 0;
    offset = 0;

    for ( uint32_t i = 0; i < size; i++ )
    {
        if ( leng <= pattern_dontcare[i] )
        {
            p_size++;
        }
        else
        {
            if ( p_size < MaxNumberPattern( weight - 2, leng + weight - 2 ) )
            {
                isImprovable = true;
                tmp_impro = true;
            }
            else
            {
                tmp_impro = false;
            }
            for ( uint32_t j = offset; j < i; j++ )
            {
                pattern_improve.push_back( tmp_impro );
            }
            p_size = 1;
            offset = i;
            leng = pattern_dontcare[i];
        }
    }
    if ( p_size < MaxNumberPattern( weight - 2, leng + weight - 2 ) )
    {
        isImprovable = true;
        tmp_impro = true;
    }
    else
    {
        tmp_impro = false;
    }
    for ( uint32_t j = offset; j < size; j++ )
    {
        pattern_improve.push_back( tmp_impro );
    }

    return isImprovable;
}

/*--Func-Change--------------------------------------------------------------*/
/**
 * Changes two different positions ('1' and '0') in a specific bit-pattern
 * Start and end are excluded
 *
 * @param number 		The pattern which has to be modified
 */
void patternset::ChangeBitsRandom( uint32_t number )
{
    uint64_t zero, one;
    bool change = false;
    while ( !change )
    {
        zero = GetSymbolRandPos( number, '0' );
        one = GetSymbolRandPos( number, '1' );
        change = ChangeBitPos( number, zero, one );
    }
}

/**
 * Exchange of two positions; changing match to don't care and vice versa.
 *
 * @param pos			The pattern index of the modifying pattern
 *
 * @param pos_one		Matchposition which has to be a don't care
 *
 * @param pos_zero		Don'tCare position which has to be a match
 *
 * @return				Boolean, if two positions could be
 * changed
 * or
 * not.
 */
bool patternset::ChangeBitPos( uint32_t number, uint64_t zero, uint64_t one )
{
    bool change;
    change = false;
    if ( initialized )
    {
        if ( pattern_improve[number % size] )
        {
            if ( bitmode )
            {
                uint64_t bit_tmp, shift;
                bit_tmp = bit_pattern_set[number % size];
                shift = 1;
                if ( bit_tmp > zero && bit_tmp > one )
                {
                    bit_tmp -= ( shift << one );
                    bit_tmp += ( shift << zero );
                    if ( UniqBit( bit_tmp ) )
                    {
                        bit_pattern_set[number] = bit_tmp;
                        change = true;
                    }
                }
            }
            else
            {
                std::vector<char> pattern_tmp( pattern_set[number] );
                if ( pattern_tmp.size() > zero && pattern_tmp.size() > one )
                {
                    std::swap( pattern_tmp[zero], pattern_tmp[one] );
                    if ( UniqPattern( pattern_set, pattern_tmp ) )
                    {
                        std::swap( pattern_set[number], pattern_tmp );
                        change = true;
                    }
                    pattern_tmp.clear();
                }
            }
        }
    }
    return change;
}

/**
 * Investigates for a specific symbol all positions of a pattern.
 * For a match, the first and last positions are excluded!
 * Chooses a position randomly and returns the position
 *
 * @param number		The pattern index of the modifying pattern
 *
 * @param symb			The specific symbol, either '1' or '0'
 *
 * @return 				A random chosen position of all
 * positions
 * matching the symbol.
 */
uint64_t patternset::GetSymbolRandPos( uint32_t number, char symb )
{
    std::vector<uint32_t> positions;
    int pos = -1;
    if ( pattern_improve[number % size] )
    {
        pos = GetLength( number ) - 1;
        if ( bitmode )
        {
            uint64_t tmp_pat, one;
            one = 1;
            for ( int i = 1; i < pos; i++ )
            {
                tmp_pat = bit_pattern_set[number % size] & ( one << ( i ) );
                if ( tmp_pat != 0 && symb == '1' )
                {
                    positions.push_back( i );
                }
                if ( tmp_pat == 0 && symb == '0' )
                {
                    positions.push_back( i );
                }
            }
        }
        else
        {
            for ( int i = 1; i < pos; i++ )
            {
                if ( pattern_set[number][i] == symb )
                {
                    positions.push_back( i );
                }
            }
        }

        std::uniform_int_distribution<uint32_t> posits( 0, positions.size() - 1 );
        pos = posits( generator );
        pos = positions[pos];
        positions.clear();
    }
    return pos;
}

/*---stuff-------------------------------------------------------------------*/
/**
 * A Method to collect all errormessages. Just easier for programmer to change
 *  	the text or extend.
 *
 * @param errmsg			Due to a few possible errormessages,
 * this
 * is
 * the
 * option, which has to be printed.
 *
 * @param pos			The position of the incorrect patterns.
 */
void patternset::SecureMessage( std::string errmsg, int pos )
{
    if ( errmsg == "file" )
    {
        printf( "%c[1;31mERROR! ", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "Pattern file \'" << inFile << "\' could not be found!" << std::endl;
        std::cerr << "Return to submitted or default values\n" << std::endl;
        return;
    }
    if ( errmsg == "empty" )
    {
        printf( "%c[1;31mERROR! ", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "File \'" << inFile << "\' is an empty file!" << std::endl;
        std::cerr << "Return to submitted or default values\n" << std::endl;
        return;
    }
    if ( errmsg == "pattern" )
    {
        printf( "%c[1;31mERROR! ", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "Patternconditions from pattern " << pos + 1 << " file were not correct (different weight, etc.)!"
                  << std::endl;
        std::cerr << "Go on to next pattern.\n" << std::endl;
        return;
    }
    if ( errmsg == "startend" )
    {
        printf( "%c[1;31mFORMAT-ERROR: ", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "Pattern " << pos + 1 << " has to start and end with a match position '1' !\n" << std::endl;
        return;
    }
    if ( errmsg == "format" )
    {
        printf( "%c[1;31mFORMAT-ERROR: ", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "Pattern containes illegal characters!" << std::endl;
        std::cerr << "Allowed characters: '0','1' for pattern; ','|'.'|';'|' ' "
                     "to seperate patterns"
                  << std::endl;
        std::cerr << "Go on to next pattern.\n" << std::endl;
        return;
    }
    if ( errmsg == "nkl" )
    {
        printf( "%c[1;31mERROR! ", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "Wrong values for weight, pattern number or pattern length!" << std::endl;
        std::cerr << "Return to default values\n" << std::endl;
        return;
    }
    if ( errmsg == "max_number_pattern" )
    {
        printf( "%c[1;33m", 27 );
        std::cerr << "The number of patterns is too high for your "
                     "configuration and will be reset!"
                  << std::endl;
        printf( "%c[0m", 27 );
        return;
    }
    if ( errmsg == "weight" )
    {
        printf( "%c[1;31mERROR! ", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "By comparing with the first pattern, the " << pos + 1 << ". pattern has a different weight!\n"
                  << std::endl;
        return;
    }
    if ( errmsg == "update" )
    {
        printf( "%c[1;31m\n####IMPROTANT####\n", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "Due to some configuration errors, your submitted "
                     "parameters have been updated!"
                  << std::endl;
        return;
    }
}

/**
 *	Estimates the number of matches ('1') in a std::vector<char>.
 *
 * @param pattern 		Pattern as std::vector<char>
 *
 * @return 				Number of matches ('1').
 */
uint32_t patternset::GetWeight( std::vector<char> pattern )
{
    uint32_t weight_p = 0;
    for ( uint32_t i = 0; i < pattern.size(); i++ )
    {
        if ( pattern[i] == '1' )
        {
            weight_p++;
        }
    }
    return weight_p;
}

/**
 *	Estimates if a std::vector<char> pattern is in binary format ('1' or
 *'0').
 *
 * @param pattern 		Pattern as std::vector<char>
 *
 * @return 				Boolean, if the submitted
 *std::vector<char>
 *pattern is binary.
 */
bool patternset::IsBinary( std::vector<char> pattern )
{
    bool binary = true;
    for ( uint32_t i = 0; i < pattern.size(); i++ )
    {
        if ( pattern[i] != '0' && pattern[i] != '1' )
        {
            binary = false;
        }
    }
    return binary;
}

/**
 * Creates a bitpattern (64-Bit Integer) from a char vector pattern.
 *
 *@param pattern 		Pattern as std::vector<char>
 *
 *@return 				Pattern as uint64_t (64-Bit Integer)
 */
uint64_t patternset::StringToBit( std::vector<char> pattern )
{
    uint64_t bitpat, tmp, tmp_size;

    bitpat = 0;
    tmp = 1;
    tmp_size = (uint64_t) pattern.size() - 1;

    for ( uint32_t i = 0; i < pattern.size(); i++ )
    {
        if ( pattern[i] == '1' )
        {
            bitpat += ( tmp << ( tmp_size - i ) );
        }
    }
    return bitpat;
}

/**
 * Creates a char vector pattern. from a bitpattern (64-Bit Integer).
 *
 *@param pattern 		Pattern as uint64_t (64-Bit Integer)
 *
 *@return 				Pattern as std::vector<char>
 */
std::vector<char> patternset::BitToString( uint64_t bitpat )
{
    std::vector<char> pattern;
    uint64_t one, tmp;

    one = 1;

    while ( bitpat > 0 )
    {
        tmp = bitpat & one;
        if ( tmp == 0 )
        {
            pattern.push_back( '0' );
        }
        else
        {
            pattern.push_back( '1' );
        }
        bitpat >>= one;
    }
    std::reverse( pattern.begin(), pattern.end() );
    return pattern;
}

/**
 *	Changes the current patternset into a bitpattern set.
 */
void patternset::PatternToBit()
{
    uint64_t tmp;
    for ( uint32_t i = 0; i < size; i++ )
    {
        tmp = StringToBit( pattern_set[i] );
        bit_pattern_set.push_back( tmp );
    }
}

/**
 *	Changes the current bitpattern set into a patternset (vector of char
 *vectors).
 */
void patternset::BitToPattern()
{
    std::vector<char> tmp;
    for ( uint32_t i = 0; i < size; i++ )
    {
        pattern_set[i] = BitToString( bit_pattern_set[i] );
    }
}

/*---------------------------------------------------------------------------*/
/**
 * Determines the maximum number of patterns with weight and length
 *
 *@param weight			Complete pattern weight-2, start and end have
 *to be match and do not change
 *
 *@param length			Complete pattern lengt-2, start and end have to
 *be match and do not change
 *
 *@return 			max number of possible pattern
 */
double patternset::MaxNumberPattern( uint32_t k, uint32_t l )
{
    double max_number;

    max_number = 1;

    for ( uint32_t i = 1; i <= k; i++ )
    {
        max_number *= ( (double) ( l - i + 1 ) ) / (double) i;
    }

    return max_number;
}

/**
 * Scans if there is another pattern in the same format in a patternset
 *
 * @param pats 			The patternset, in which we are looking for a
 * similar pattern.
 *
 * @param pat 			The pattern, the query, we are looking for to be
 * unique.
 *
 * @return 				Returns a boolean, if the submitted
 * pattern
 * is
 * not existing in the patternset.
 */
bool patternset::UniqPattern( std::vector<std::vector<char>> pats, std::vector<char> pat )
{
    bool uniq, no_common;
    uniq = true;

    for ( uint32_t i = 0; i < pats.size(); i++ )
    {
        if ( pats[i].size() == pat.size() )
        {
            no_common = false;
            for ( uint32_t j = 0; j < pat.size(); j++ )
            {
                if ( pats[i][j] != pat[j] )
                {
                    no_common = true;
                    j = pat.size();
                }
            }
            if ( !no_common )
            {
                uniq = false;
            }
        }
    }

    return uniq;
}

/**
 * Scans if there is another bitpattern in the same format in a bitpattern set
 *
 * @param pats 			The bitpattern set, in which we are looking for
 * a
 * similar bitpattern.
 *
 * @param pat 			The bitpattern, the query, we are looking for to
 * be
 * unique.
 *
 * @return 				Returns a boolean, if the submitted
 * bitpattern
 * is
 * not existing in the bitpattern set.
 */
bool patternset::UniqBit( uint64_t bit_pat )
{
    return std::find( bit_pattern_set.begin(), bit_pattern_set.end(), bit_pat ) == bit_pattern_set.end();
}

/**
 * Prints the current patternset to the commandline.
 */
void patternset::Print()
{
    if ( bitmode )
    {
        BitToPattern();
    }
    for ( uint32_t i = 0; i < size; i++ )
    {
        for ( uint32_t j = 0; j < pattern_set[i].size(); j++ )
        {
            std::cout << pattern_set[i][j];
        }
        std::cout << std::endl;
    }
}

/*---------------------------------------------------------------------------*/
/**
 * Returns the patternset as std::vector< std::vector<char> >
 *
 * @return 				Copy of the current patternset vector
 * (type:
 * std::vector< std::vector<char> >)
 */
std::vector<std::vector<char>> patternset::GetPattern()
{
    std::vector<std::vector<char>> tmp;
    if ( bitmode )
    {
        BitToPattern();
    }
    return pattern_set;
}

/**
 * Returns the patternset as bitpattern as std::vector<uint64_t>
 *
 * @return 				Copy of the current bit-patternset
 * vector
 * (type:
 * std::vector<uint64_t>)
 */
std::vector<uint64_t> patternset::GetBitPattern()
{
    std::vector<uint64_t> tmp;
    if ( initialized )
    {
        if ( bitmode )
        {
            tmp = bit_pattern_set;
        }
    }
    return tmp;
}

/**
 * Returns a pattern of the patternset as std::vector<char>
 *
 * @return 				Copy of the 'number'th pattern of the
 * current
 * patternset vector (type: std::vector<char>)
 */
std::vector<char> patternset::GetPattern( uint32_t number )
{
    std::vector<char> tmp;
    if ( bitmode )
    {
        tmp = BitToString( bit_pattern_set[number % size] );
    }
    else
    {
        tmp = pattern_set[number % size];
    }
    return tmp;
}

/**
 * Returns a pattern of the patternset as uint64_t.
 *
 * @return 				Copy of the 'number'th pattern of the
 * current
 * patternset vector (type: uint64_t)
 */
uint64_t patternset::GetBitPattern( uint32_t number )
{
    if ( !bitmode )
    {
        return 0;
    }
    return bit_pattern_set[number % size];
}

/**
 * Returns a vector in which (in pattern order) for each pattern there is a
 * boolean, if the pattern
 * can be changed or not.
 *
 * @return 				Vector<bool> which contains for each
 * pattern
 * the
 * boolean of improvement.
 */
std::vector<uint32_t> patternset::GetDontCareVec()
{
    std::vector<uint32_t> tmp;
    if ( initialized )
    {
        tmp = pattern_dontcare;
    }
    return tmp;
}

/**
 * Returns a boolean, if the 'number'th patterncan be improved/changed, or not.
 *
 * @return 				Boolean, if the 'number'th patterncan be
 * improved/changed.
 */
bool patternset::GetImprove( uint32_t number )
{
    bool tmp = false;
    if ( initialized )
    {
        tmp = pattern_improve[number % size];
    }
    return tmp;
}

/**
 * Returns a boolean, if there is at least one pattern, which can be
 * improved/changed, or not.
 *
 * @return 				Boolean, if at least one pattern can be
 * changed.
 */
bool patternset::GetImprove()
{
    return improve;
}

/**
 * Returns the current patternset weight.
 *
 * @return 				Current pattern weight.
 */
uint32_t patternset::GetSize()
{
    return size;
}

/**
 * Returns the current patternset weight.
 *
 * @return 				Current pattern weight.
 */
uint32_t patternset::GetWeight()
{
    return weight;
}

/**
 * Returns the current lenght of the 'number'th pattern.
 * Weight + don't care positions
 *
 *@return 				Current lenght of the 'number'th
 *pattern.
 */
uint32_t patternset::GetLength( uint32_t number )
{
    uint32_t leng = 0;
    if ( initialized )
    {
        leng = pattern_dontcare[number % size] + weight;
    }
    return leng;
}

/**
 * Returns boolean if the current patternset modus is bitmode oder not.
 *
 * @return 				Current bitmode boolean (true, if
 * bitmode).
 */
bool patternset::GetBitMode()
{
    return bitmode;
}

/**
 * Returns boolean if the current patternset modus has been updated
 *
 * @return 				Current patternmode has been updated.
 */
bool patternset::GetUpdate()
{
    return update;
}

/**
 * Sets a bitpattern at the 'number'th positions of the bitpattern set.
 *
 * @param number 		The position where the bitpattern should be
 * saved.
 *
 * @param pat 			The bitpattern, that should be saved.
 */
void patternset::SetPattern( uint32_t number, uint64_t pat )
{
    if ( (uint32_t) std::log2( pat ) == (uint32_t) std::log2( bit_pattern_set[number % size] ) )
    {
        bit_pattern_set[number % size] = pat;
    }
}

/**
 * Sets a pattern at the 'number'th positions of the patternset.
 *
 * @param number 		The position where the pattern should be saved.
 *
 * @param pat 			The pattern, that should be saved.
 */
void patternset::SetPattern( uint32_t number, std::vector<char> pat )
{
    uint64_t tmp;
    if ( bitmode && pat.size() == ( uint32_t )( std::log2( bit_pattern_set[number % size] ) + 1 ) )
    {
        tmp = StringToBit( pat );
        bit_pattern_set[number % size] = tmp;
    }
    else if ( pat.size() == pattern_set[number % size].size() )
    {
        pattern_set[number % size] = pat;
    }
}