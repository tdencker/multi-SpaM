/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of
 * pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * sensitivity object file
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
#include "sensitivity.h"

/*---Variables---------------------------------------------------------------*/

/*===Main-Part===============================================================*/
/*---Constructor-------------------------------------------------------------*/
/**
 * Default constructor, sets the default vaulues for std. output on cmd-line.
 */
sensitivity::sensitivity()
{
    p = 0.75;
    q = 0.25;
    seq_leng = 16000;
    H = 64;
    size = 10;
    weight = 8;
    max_dontcare = 6;
    min_dontcare = 6;
    isInFile = false;
    init = false;
}

/**
 * File constructor, where patternset will be read from a file.
 *
 * @param inFile 		File, which contains the patternset.
 *
 * @param p				Match probability ( = #matches /
 * #positions),
 * used for sensitivity and variance.
 *
 * @param q				Background probability, used for
 * variance.
 *
 * @param seq_leng		Length of sequences in database, used for
 * variance.
 *
 * @param H				Length of a homolgue random region in a
 * Dataset.
 */
sensitivity::sensitivity( const char * inFile, double p, double q, uint32_t H,
                          int seq_leng )
{
    this->p = p;
    this->H = H;
    this->inFile = inFile;
    this->q = q;
    this->seq_leng = seq_leng;
    isInFile = true;
    init = false;
}

/**
 * Long constructor, where the complete pattern and sensitivity options can be
 * set.
 *
 * @param size 			Number of patterns for the patternset.
 *
 * @param weight		Number of match positions ('1') in each pattern.
 *
 * @param min_dontcare	Number of minimum don't care positions ('0') in each
 * pattern.
 *
 * @param max_dontcare	Number of maximum don't care positions ('0') in each
 * pattern.
 *
 * @param p				Match probability ( = #matches /
 * #positions),
 * used for sensitivity and variance.
 *
 * @param q				Background probability, used for
 * variance.
 *
 * @param H				Length of a homolgue random region in a
 * Dataset.
 *
 * @param seq_leng		Length of sequences in database, used for
 * variance.
 */
sensitivity::sensitivity( uint32_t size, uint32_t weight, uint32_t min_dontcare,
                          uint32_t max_dontcare, double p, double q, double H,
                          uint32_t seq_leng )
{
    this->p = p;
    this->q = q;
    this->H = H;
    this->size = size;
    this->weight = weight;
    this->max_dontcare = max_dontcare;
    this->min_dontcare = min_dontcare;
    this->seq_leng = seq_leng;
    isInFile = false;
    init = false;
}

/*---Initialization----------------------------------------------------------*/
/**
 * Standard Initialization of sensitivity options.
 */
void sensitivity::Init()
{
    if ( !init )
    {
        std::string str = "outpat.pat";
        sens = false;
        oc = true;
        silent = false;
        quiet = true;
        randpatleng = false;
        bitmode = true;
        improve = false;
        outFile = str.c_str();
        Init( sens, oc, improve, randpatleng, quiet, silent, outFile );
    }
}

/**
 * Initialization of sensitivity options with specific parameter.
 *
 * @param sens 			Boolean which (de)activates the sensitivity
 * calculation.
 *
 * @param oc 			Boolean which (de)activates the overlap
 * complexity/variance calculation.
 *
 * @param improve 		Boolean which (de)activates the improvement of
 * the
 * patternset.
 *
 * @param randpatleng	Boolean which (de)activates random generated pattern
 * lengths.
 *
 * @param quiet			Boolean which (de)activates just a short output
 * of
 * the improvement.
 *
 * @param silent 		Boolean which (de)activates std::output on
 * commandline (exclusive errors).
 *
 * @param outFile 		Char array which contains name and path to
 * outputfile (commandline parameter).
 */
void sensitivity::Init( bool sens, bool oc, bool improve, bool randpatleng,
                        bool quiet, bool silent, const char * outFile )
{
    if ( !init )
    {
        this->sens = sens;
        this->oc = oc;
        this->improve = improve;
        this->randpatleng = randpatleng;
        this->quiet = quiet;
        this->silent = silent;
        this->outFile = outFile;
        opt_sens = 1;
        opt_oc = 1;
        ReInit();
    }
}

/**
 * Standard Reinitialization of sensitivity options, with a new patternset.
 */
void sensitivity::ReInit()
{
    init = false;
    if ( isInFile )
    {
        var_pattern = variance( inFile, p, q, seq_leng );
    }
    else
    {
        var_pattern = variance( size, weight, min_dontcare, max_dontcare,
                                seq_leng, p, q );
    }
    var_pattern.Init( oc, improve, quiet, true, randpatleng, NULL );
    ReInit( var_pattern );
}

/**
 * Standard Reinitialization of sensitivity options, with a given patternset.
 */
void sensitivity::ReInit( variance & var_pattern )
{
    this->var_pattern = var_pattern;
    init = false;
    if ( outFile != NULL && outFile[0] != '\0' )
    {
        isOutFile = true;
    }
    else
    {
        isOutFile = false;
    }
    bitmode = var_pattern.GetBitMode();
    size = var_pattern.GetSize();
    weight = var_pattern.GetWeight();
    max_dontcare = var_pattern.GetMaxDontcare();
    min_dontcare = var_pattern.GetMinDontcare();
    update = var_pattern.GetUpdate();
    if ( improve )
    {
        improve = var_pattern.GetImprove();
    }
    if ( !bitmode )
    {
        if ( sens )
        {
            SecureMessage( "bitmode", -1 );
            update = true;
        }
        sens = false;
    }
    if ( sens )
    {
        sens_value =
            CalculateSensitivity( var_pattern.GetPattern(), p, size, H );
    }
    else
    {
        sens_value = 0;
    }
    init = true;
}

/*---Improvement-------------------------------------------------------------*/
/**
 * Improvement Method, which starts the improvement method from the variance
 * object.
 *
 * @param limits 		Number of tries of random permutation to
 * optimize
 * the
 * ov/variance of the patternset
 *
 * @param opt_oc		Number of reinitializations of the patternset,
 * to
 * find
 * the best oc/variance of the permutated patternsets.
 */
void sensitivity::Improve( uint32_t limits, uint32_t opt_oc )
{
    this->opt_oc = opt_oc;
    var_pattern.Improve( limits, opt_oc );
    if ( sens )
    {
        sens_value =
            CalculateSensitivity( var_pattern.GetPattern(), p, size, H );
    }
}

/**
 * Improvement Method, which improves the sensitivity improvement incl. the
 * variance/oc improvement.
 *
 * @param limits 		Number of tries of random permutation to
 * optimize
 * the
 * ov/variance of the patternset
 *
 * @param opt_oc		Number of reinitializations of the patternset,
 * to
 * find
 * the best oc/variance of the permutated patternsets.
 *
 * @param opt_sens		Number of reinitializations of the patternset,
 * to
 * find the best sensitivity of the best oc/variance optimized permutated
 * patternsets.
 */
void sensitivity::Improve( uint32_t limits, uint32_t opt_oc, uint32_t opt_sens )
{
    variance best_var_pattern;
    std::vector<char> out_pat;
    std::ofstream pattern_out;
    std::string out_str;
    double time_begin, best_time, to_second, best_sens_value;
    uint32_t better_pattern;

    this->opt_sens = opt_sens;
    this->opt_oc = opt_oc;
    to_second = 1 / 1000000.0;
    time_begin = clock();
    best_time = 0;

    if ( opt_sens == 1 || !sens )
    {
        var_pattern.ImproSilent( true );
    }
    if ( !silent )
    {
        std::cout << "\n ===== First patternset =====" << std::endl;
        var_pattern.Print();
        std::cout << var_pattern.GetFormat() << var_pattern.GetVariance()
                  << std::endl;
        std::cout << "norm_" << var_pattern.GetFormat()
                  << var_pattern.GetNormVariance() << std::endl;
        if ( sens )
        {
            std::cout << "sensitivity\t\t\t= " << sens_value << std::endl;
        }
        std::cout << "time\t\t\t\t= " << best_time * to_second << " sec"
                  << std::endl
                  << std::endl;
    }
    if ( sens && improve )
    {
        best_var_pattern = var_pattern;
        best_sens_value = sens_value;
        better_pattern = 0;

        for ( uint32_t i = 1; i <= opt_sens; i++ )
        {
            Improve( limits, opt_oc );
            best_time = clock() - time_begin;
            if ( sens_value > best_sens_value )
            {
                best_sens_value = sens_value;
                better_pattern++;
                if ( !silent && opt_sens != 1 )
                {
                    if ( !quiet )
                    {
                        std::cout << "\n*** BETTER PATTERN " << better_pattern
                                  << " *** \t(sensitivity optimization)"
                                  << std::endl;
                        std::cout << "Step " << i << " / " << opt_sens
                                  << std::endl
                                  << "Patternset: \n";
                        var_pattern.Print();
                        std::cout << var_pattern.GetFormat()
                                  << var_pattern.GetVariance() << std::endl;
                        std::cout << "norm_" << var_pattern.GetFormat()
                                  << var_pattern.GetNormVariance() << std::endl;
                        if ( sens )
                        {
                            std::cout << "sensitivity\t\t\t= " << sens_value
                                      << std::endl;
                        }
                        std::cout << "time\t\t\t\t= " << best_time * to_second
                                  << " sec" << std::endl
                                  << std::endl;
                    }
                    else
                    {
                        std::cout << "\r*** BETTER PATTERN " << better_pattern
                                  << " *** \t(sensitivity optimization)";
                        std::cout.flush();
                    }
                }
                best_var_pattern = var_pattern;
            }
            ReInit();
        }
        var_pattern = best_var_pattern;
        ReInit( best_var_pattern );
    }
    if ( !sens && improve )
    {
        Improve( limits, opt_oc );
        best_time = clock() - time_begin;
    }
    if ( !silent && improve )
    {
        std::cout << "\n\n ===== Best patternset ======" << std::endl;
        var_pattern.Print();
        std::cout << var_pattern.GetFormat() << var_pattern.GetVariance()
                  << std::endl;
        std::cout << "norm_" << var_pattern.GetFormat()
                  << var_pattern.GetNormVariance() << std::endl;
        if ( sens )
        {
            std::cout << "sensitivity\t\t\t= " << sens_value << std::endl;
        }
        std::cout << "time\t\t\t\t= " << best_time * to_second << " sec"
                  << std::endl
                  << std::endl;
    }
    if ( isOutFile )
    {
        std::cout << "=======================================\nSaving pattern "
                     "into file '"
                  << outFile << "' ..." << std::endl;
        pattern_out.open( outFile );
        for ( uint32_t i = 0; i < size; i++ )
        {
            out_pat = var_pattern.GetPattern( i );
            out_str = std::string( out_pat.begin(), out_pat.end() );
            pattern_out << out_str << std::endl;
        }
        pattern_out << var_pattern.GetFormat() << var_pattern.GetVariance()
                    << std::endl;
        pattern_out << "norm_" << var_pattern.GetFormat()
                    << var_pattern.GetNormVariance() << std::endl;
        if ( sens )
        {
            pattern_out << "sensitivity\t\t\t= " << sens_value << std::endl;
        }
        pattern_out << "time\t\t\t\t= " << best_time * to_second << " sec"
                    << std::endl
                    << std::endl;
        std::cout << "... Done!\n======================================="
                  << std::endl;
    }
    if ( !silent && update )
    {
        SecureMessage( "update", -1 );
    }
}

/*---stuff-------------------------------------------------------------------*/
/**
 * Prints the current patternset to the commandline.
 */
void sensitivity::Print()
{
    var_pattern.Print();
}

/**
 * Returns the current patternset as vector of char vectors.
 *
 * @return Patternset in std::vector< std::vector<char> > format.
 */
std::vector<std::vector<char>> sensitivity::GetPattern()
{
    return var_pattern.GetPattern();
}

/**
 * Returns the current sensitvity (no calculation: 0).
 *
 * @return 				The sensitivity of the current set (0 if
 * not
 * calculated)
 */
double sensitivity::GetSensitivity()
{
    return sens_value;
}

/**
 * Returns the current patternset size.
 *
 * @return 				Current pattern size.
 */
uint32_t sensitivity::GetSize()
{
    return size;
}

/**
 * Returns the current patternset weight.
 *
 * @return 				Current pattern weight.
 */
uint32_t sensitivity::GetWeight()
{
    return weight;
}

/**
 * Returns the current patternset minimum don't care number.
 *
 * @return 				Current pattern minimum don't care
 * number.
 */
uint32_t sensitivity::GetMaxDontcare()
{
    return max_dontcare;
}

/**
 * Returns the current patternset maximum don't care number.
 *
 * @return 				Current pattern maximum don't care
 * number.
 */
uint32_t sensitivity::GetMinDontcare()
{
    return min_dontcare;
}

/**
 * Returns the current dataset sequence length for variance calculation.
 *
 * @return 				Current dataset sequence length for
 * variance
 * calculation.
 */
uint32_t sensitivity::GetSequenceLength()
{
    return seq_leng;
}

/**
 * Returns the current match probability for variance calculation.
 *
 * @return 				Current match probability for variance
 * calculation.
 */
double sensitivity::GetP()
{
    return p;
}

/**
 * Returns the current background probability for variance calculation.
 *
 * @return 				Current background probability for
 * variance
 * calculation.
 */
double sensitivity::GetQ()
{
    return q;
}

/**
 * Returns the current length of a random homolgue region in a possible datasat.
 *
 * @return 				Current length of a random homolgue
 * region
 * in
 * a
 * possible datasat.
 */
uint32_t sensitivity::GetH()
{
    return H;
}

void sensitivity::SecureMessage( std::string errmsg, int pos )
{
    if ( errmsg == "noimprove" )
    {
        printf( "%c[1;33m", 27 );
        std::cerr << "Using your pattern conditions it is not sensible to "
                     "improve your pattern, sorry!"
                  << std::endl;
        std::cerr << "Deactivating improve mode\n" << std::endl;
        printf( "%c[0m", 27 );
        return;
    }
    if ( errmsg == "update" )
    {
        printf( "%c[1;31m\n####IMPROTANT####\n", 27 );
        printf( "%c[0m", 27 );
        std::cerr << "Due to some configuration errors (see above), your "
                     "submitted parameters have been updated!\n"
                  << std::endl;
        return;
    }
    if ( errmsg == "bitmode" )
    {
        printf( "%c[1;33m", 27 );
        std::cerr << "A patternlength is over 63, leaving bitmode ..."
                  << std::endl;
        std::cerr << "Using your pattern conditions it is not possible to "
                     "calculate the sensitivity!"
                  << std::endl;
        std::cerr << "Deactivating sensitivity calculation\n" << std::endl;
        printf( "%c[0m", 27 );
        return;
    }
}

/*---------------------------------------------------------------------------*/
/**
 * Interface function to create the right pattern format for the speed
 * functions.
 *
 * @param pattern	Contains the patternset.
 *
 * @param p		The match probability ( = #matches / #l_hom).
 *
 * @param n		Number of Patterns.
 *
 * @param R		The length of random homolog region in a dataset (=H).
 */
double
sensitivity::CalculateSensitivity( std::vector<std::vector<char>> pattern,
                                   double p, int n, int R )
{
    double sens;
    uint32_t p_size = pattern.size(), length;

    char ** pats = new char *[p_size]; // set of seeds
    for ( uint32_t i = 0; i < p_size; i++ )
    {
        length = pattern[i].size();
        pats[i] = new char[length + 1];
    }

    for ( uint32_t i = 0; i < p_size; i++ )
    {
        length = pattern[i].size();
        for ( uint32_t j = 0; j < length; j++ )
        {
            if ( pattern[i][j] == '1' )
            {
                pats[i][j] = '1';
            }
            else
            {
                pats[i][j] = '0';
            }
        }
        pats[i][length] = '\0';
    }
    sens = MULTIPLE_SENSITIVITY2( pats, n, R, p );

    return sens;
}

/*===SPEED===================================================================*/
/*---Acknowledgements--------------------------------------------------------*/
/**
 * We would like to thank Ilie and Ilie for using there programm part of
 * SpEED to calculate the sensitivity for a patternset off:
 *
 * Lucian Ilie, Silvana Ilie, and Anahita M. Bigvand. SpEED: fast computation
 * of sensitive spaced seeds. Bioinformatics, 27:2433–2434, 2011.
 */
inline long long sensitivity::BIN_REVERSED_TO_INT2(
    char * s ) // converts the reversed of the binary string s into integer
{              // works also with * instead of 0
    long long val = 0;
    long long l = strlen( s ), i = 0, temp = 1;
    for ( i = 0; i <= l - 1; i++ )
    {
        if ( s[i] == '1' )
            val += temp;
        temp *= 2;
    }
    return ( val );
}
/**
* Computing sensitivity of a set of SEEDS with the given parameters
* using the dynamic programming of (Li et al., 2004)
*/
double sensitivity::MULTIPLE_SENSITIVITY2( char ** SEEDS, int NO_SEEDS,
                                           long long N, double P )
{
    long long i = 0, j = 0, b = 0, pos = 0, MAX_L = 0, level = 0,
              prev_level_start = 0, prev_level_end = 0, compatible = 0, hit = 0,
              suffix_link = 0, zero_link = 0, new_i = 0, tmp = 0;
    long long b_zero = 0, b_one = 0;
    double f0 = 0, f1 = 0;
    // compute the lengths of the seeds and MAX_L = the length of the longest
    // seed

    long long * seed_length = new long long[NO_SEEDS];
    for ( i = 0; i <= NO_SEEDS - 1; i++ )
    {
        seed_length[i] = strlen( SEEDS[i] );
        // cout << seed_length[i] << endl;
        if ( MAX_L < seed_length[i] )
            MAX_L = seed_length[i];
    }
    // compute the integer values of the reversed seeds INTeger REVersed SEEDS
    long long * INT_REV_SEEDS = new long long[NO_SEEDS];
    for ( i = 0; i <= NO_SEEDS - 1; i++ )
        INT_REV_SEEDS[i] =
            BIN_REVERSED_TO_INT2( SEEDS[i] ); // !!! this works like * = 0

    // create the tree of BS --- 1..NO_BS-1
    // *********************************************************
    // BS[i][0] = the integer value of b^r (except for epsilon, any b starts
    // with 1)
    // BS[i][1] = index j in BS of left son: BS[j][0] = integer value of b^r0 =
    // (0b)^r (-1 if it doesn't exist)
    // BS[i][2] = index j in BS of right son: BS[j][0] = integer value of b^r1 =
    // (1b)^r (-1 if it doesn't exist)
    // BS[i][3] = (suffix link) index j in BS of (B(b^r))^r i.e. BS[j][0] =
    // integer value of (B(b^r))^r
    // BS[i][4] = 1 if b is a hit and 0 otherwise			//B(x)
    // is
    // the
    // longest
    // prefix
    // of x that is in B
    // B = set of compatible but not hits b's
    // BS[i][5] = its level = the length of the string
    // BS[i][6] = the longest prefix of 0b which reversed means the longest
    // suffix followed by a 0

    long long MAX_NO_BS = NO_SEEDS;
    for ( i = 0; i <= NO_SEEDS - 1; i++ )
    { // compute maximum possible no of b's
        tmp = 1;
        for ( j = strlen( SEEDS[i] ) - 1; j >= 0; j-- )
        {
            if ( SEEDS[i][j] != '1' )
                tmp *= 2;
            MAX_NO_BS += tmp; // add previous value if 1 in seed or double (tmp
                              // *= 2 above) if a * in seed
        }
    }
    long long NO_BS = MAX_NO_BS;
    // cout << NO_BS << endl;
    // bound for computing sensitivity (not allocate more than 120GB)

    long long ** BS = new long long *[NO_BS];
    for ( i = 0; i <= NO_BS - 1; i++ )
    {
        BS[i] = new long long[7];
        BS[i][0] = BS[i][3] = BS[i][4] = BS[i][5] = BS[i][6] = 0; // initialize
        BS[i][1] = BS[i][2] = -1;
    };
    // create the tree by levels: all b's of length i are on level i
    BS[0][0] = 0;  // epsilon
    BS[0][1] = -1; // no left son since b=0 is not compatible -- seeds end with
                   // 1
    BS[0][2] = 1;  // right son is BS[1][0] = 1
    BS[0][3] = 0;  // suffix link to itself
    BS[0][4] = 0;  // epsilon is not hit
    BS[0][5] = 0;
    BS[0][6] = 0;
    BS[1][0] = 1;
    BS[1][3] = 0;
    BS[1][4] = 0; // assume 1 is not a hit
    BS[1][5] = 1;
    prev_level_start = 1;
    prev_level_end = 1; // indices in BS between which previous level is found
    pos = 2;            // first empty position in BS
    for ( level = 2; level <= MAX_L; level++ )
    { // complete level "level"
        for ( i = prev_level_start; i <= prev_level_end; i++ )
            if ( BS[i][4] != 1 )
            {                   // not a hit
                b = BS[i][0];   // integer value
                b_zero = 2 * b; // try b0
                compatible = 0;
                hit = 0;
                for ( j = 0; j <= NO_SEEDS - 1; j++ ) // check long enough seeds
                                                      // to seee if b0 is
                                                      // compat/hit
                    if ( seed_length[j] >= level )
                        if ( ( ( INT_REV_SEEDS[j] >>
                                 ( seed_length[j] - level ) ) &
                               ( ~b_zero ) ) == 0 )
                        {
                            compatible = 1;
                            if ( level == seed_length[j] )
                                hit = 1;
                        }
                if ( compatible )
                {
                    BS[i][1] = pos;
                    BS[pos][0] = 2 * b;
                    BS[pos][4] = hit; // hit = 1 if it is hit by a seed
                    BS[pos][5] = level;
                    suffix_link = BS[i][3];
                    while ( ( suffix_link != 0 ) &&
                            ( BS[suffix_link][1] == -1 ) )
                    {
                        suffix_link = BS[suffix_link][3];
                    }
                    if ( suffix_link != 0 )
                    {
                        BS[pos][3] = BS[suffix_link][1];
                        if ( BS[BS[pos][3]][4] == 1 ) // if suffix link is hit
                                                      // then also itself is hit
                            BS[pos][4] = 1;
                    }
                    pos++;
                }
                b_one = 2 * b + 1; // try b1
                compatible = 0;
                hit = 0;
                for ( j = 0; j <= NO_SEEDS - 1; j++ ) // check all long enough
                                                      // seed to seee if b0 is
                                                      // compat/hit
                    if ( seed_length[j] >= level )
                        if ( ( ( INT_REV_SEEDS[j] >>
                                 ( seed_length[j] - level ) ) &
                               ( ~b_one ) ) == 0 )
                        {
                            compatible = 1;
                            if ( level == seed_length[j] )
                                hit = 1;
                        }
                if ( compatible )
                {
                    BS[i][2] = pos;
                    BS[pos][0] = 2 * b + 1;
                    BS[pos][4] = hit; // hit = 1 if it is hit by a seed
                    BS[pos][5] = level;
                    suffix_link = BS[i][3];
                    while ( BS[suffix_link][2] == -1 )
                    {
                        suffix_link = BS[suffix_link][3];
                    }
                    BS[pos][3] = BS[suffix_link][2];
                    if ( BS[BS[pos][3]][4] ==
                         1 ) // if suffix link is hit then also itself is hit
                        BS[pos][4] = 1;

                    pos++;
                }
            }
        prev_level_start = prev_level_end + 1;
        prev_level_end = pos - 1;
    }
    // zero_links -- longest suffix of b0 in the tree
    for ( i = 1; i <= NO_BS - 1; i++ )
        if ( BS[i][1] != -1 ) // has left son, that is, 0-son
            BS[i][6] = BS[i][1];
        else
        {
            zero_link = BS[i][3];
            while ( ( zero_link != 0 ) && ( BS[zero_link][1] == -1 ) )
            {
                zero_link = BS[zero_link][3];
            }
            if ( zero_link != 0 )
                BS[i][6] = BS[zero_link][1];
        }
    // compute the f's  f[i][j] = probab to hit a prefix of length i that ends
    // with INT_TO_BIN_REVERSED[BS[j][0]]
    double ** f;
    f = new double *[N + 1];
    for ( i = 0; i <= N; i++ )
    {
        f[i] = new double[NO_BS];
        for ( j = 0; j <= NO_BS - 1; j++ )
            f[i][j] = 0; // initialize
    }

    for ( i = 0; i <= N; i++ )
    {
        for ( j = NO_BS - 1; j >= 0; j-- )
        {
            if ( i == 0 )
                f[i][j] = 0; // empty prefix of random region cannot be hit
            else if ( i < BS[j][5] )
                f[i][j] = 0; // too short
            else if ( BS[j][4] == 1 )
                f[i][j] = 1; // hit
            else
            {
                new_i = i - BS[j][5] + BS[BS[j][6]][5] - 1;
                if ( new_i < 0 )
                    new_i = 0;
                f0 = f[new_i][BS[j][6]];
                if ( BS[j][2] < 0 )
                    f1 = 1;
                else
                    f1 = f[i][BS[j][2]];
                f[i][j] = ( 1 - P ) * f0 + P * f1;
            }
        }
    }
    double result = f[N][0];
    // free memory
    delete[] seed_length;
    delete[] INT_REV_SEEDS;
    for ( i = 0; i <= N; i++ )
        delete[] f[i];
    delete[] f;
    for ( i = 0; i <= MAX_NO_BS - 1; i++ )
        delete[] BS[i];
    delete[] BS;

    return ( result );
}
