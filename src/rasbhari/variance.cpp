/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * variance object file
 *
 * For theory please have a look at:
 *
 * - L. Hahn, C.-A. Leimeister, B. Morgenstern (2015)
 * RasBhari: optimizing spaced seeds for database searching, read mapping and alignment-free sequence comparison
 * arXiv:1511.04001[q-bio.GN]
 *
 * - B. Morgenstern, B. Zhu, S. Horwege, C.-A Leimeister (2015)
 * Estimating evolutionary distances between genomic sequences from spaced-word matches
 * Algorithms for Molecular Biology 10, 5. (http://www.almob.org/content/10/1/5/abstract)
 *
 *
 * @author: Lars Hahn - 12.05.2016, Georg-August-Universitaet Goettingen
 * @version: 1.2.0 05/2016
 */
#include "variance.h"

/*---Variables---------------------------------------------------------------*/



/*===Main-Part===============================================================*/
/*---Constructor-------------------------------------------------------------*/
/**
 * Default constructor, sets the default vaulues, pattern will be generated automatically.
 */
variance::variance(){
	size = 10;
	weight = 8;
	min_dontcare = 6;
	max_dontcare = 6;
	seq_leng = 16000;
	p = 0.75;
	q = 0.25;
	isInFile = false;
	initialized = false;
}


/**
 * Short constructor, sets some default vaulues, just pattern dimension is set.
 *
 * @param size			Number of patterns for the patternset.
 *
 * @param weigth		Number of match positions ('1') in each pattern.
 *
 * @param min_dontcare	Number of minimum don't care positions ('0') in each pattern.
 *
 * @param max_dontcare	Number of maximum don't care positions ('0') in each pattern.
 */
variance::variance(uint32_t size, uint32_t weight, uint32_t min_dontcare, uint32_t max_dontcare){
	this->size = size;
	this->weight = weight;
	this->min_dontcare = min_dontcare;
	this->max_dontcare = max_dontcare;
	seq_leng = 16000;
	p = 0.75;
	q = 0.25;
	isInFile = false;
	initialized = false;
}


/**
 * Long constructor, sets the values; resets automatically, if there are problems.
 *
 * @param size 			Number of patterns for the patternset.
 *
 * @param weight		Number of match positions ('1') in each pattern.
 * 
 * @param min_dontcare	Number of minimum don't care positions ('0') in each pattern.
 *
 * @param max_dontcare	Number of maximum don't care positions ('0') in each pattern.
 *
 * @param seq_leng		Length of sequences in database, used for variance.
 *
 * @param p				Match probability ( = #matches / #positions), used for sensitivity and variance.
 *
 * @param q				Background probability, used for variance.
 */
variance::variance(uint32_t size, uint32_t weight, uint32_t min_dontcare, uint32_t max_dontcare, uint32_t seq_leng, double p, double q){
	this->size = size;
	this->weight = weight;
	this->min_dontcare = min_dontcare;
	this->max_dontcare = max_dontcare;
	this->seq_leng = seq_leng;
	this->p = p;
	this->q = q;
	isInFile = false;
	initialized = false;
}


/**
 * File constructor, sets values only from files; resets automatically if there are problems.
 *
 * @param inFile	The file which contains a patternset.
 */
variance::variance(const char *inFile, double p, double q, uint32_t seq_leng){
	this->inFile = inFile;
	this->p = p;
	this->q = q;
	this->seq_leng = seq_leng;
	size = 10;
	weight = 8;
	min_dontcare = 6;
	max_dontcare = 6;
	isInFile = true;
	initialized = false;
}


/**
 * Default destructor, deletes all vectors and matrices in the object.
 */
variance::~variance(){
	Clear();
}


/*---Initialization----------------------------------------------------------*/
/**
 * Standard Initialization of variance options.
 */
void variance::Init(){
	if(!initialized){
		std::string str = "outfile.pat";
		oc = true;
		improve = true;
		quiet = true;
		silent = false;
		isOutFile = false;
		random_leng = false;
		outFile = str.c_str();
		Init(oc, improve, quiet, silent, random_leng, outFile);
	}
}

/**
 * Initialization of variance options with specific parameter.
 *
 * @param oc 			Boolean which (de)activates the overlap complexity/variance calculation
 *
 * @param improve 		Boolean which (de)activates the improvement of the patternset.
 *
 * @param quiet			Boolean which (de)activates just a short output of the improvement
 * 
 * @param silent 		Boolean which (de)activates std::output on commandline (exclusive errors)
 *
 * @param randpatleng	Boolean which (de)activates random generated pattern lengths
 *
 * @param outFile 		Char array which contains name and path to outputfile (commandline parameter)
 */
void variance::Init(bool oc, bool improve, bool quiet, bool silent, bool random_leng, const char *outFile){
	if(!initialized){
		this->oc = oc;
		this->improve = improve;
		this->quiet = quiet;
		this->silent = silent;
		this->outFile = outFile;
		this->random_leng = random_leng;
		best_silent = false;
		impro_silent = false;
		loops = 1;
		if(outFile != NULL && outFile[0] != '\0'){
			isOutFile = true;
		}
		else{
			isOutFile = false;
		}
		if(oc){
			outvar = "overlap-complexity\t\t= ";
		}
		else{
			outvar = "variance\t\t\t= ";
		}
		if(silent){
			quiet = true;
		}
		ReInit();
	}
}


/**
 * Standard Reinitialization of variance options, with a new patternset.
 */
void variance::ReInit(){
	if(isInFile){
		pattern_set = patternset(inFile, true);
	}
	else{	
		pattern_set = patternset(size, weight, min_dontcare, max_dontcare, true);
	}
	pattern_set.Init(random_leng);
	ReInit(pattern_set);
}


/**
 * Standard Reinitialization of variance options, with a given patternset.
 */
void variance::ReInit(patternset &pat){
	if(initialized){
		Clear();
	}
	initialized = false;
	pattern_set = pat;
	
	size = pattern_set.GetSize();
	weight = pattern_set.GetWeight();
	min_dontcare = pattern_set.GetLength(size-1)- weight;
	max_dontcare = pattern_set.GetLength(0) - weight;
	bitmode = pattern_set.GetBitMode();
	if(improve){
		improve = pattern_set.GetImprove();
		if(!improve){
			SecureMessage("noimprove",-1);
		}
	}
	norm_size = (size*(size+1))/2;
	p = std::abs(p);
	q = std::abs(q);
	if(p < q || p > 1 || q > 1){
		p = 0.75;
		q = 0.25;
		SecureMessage("pq",-1);
	}	
	InitVariance();
	initialized = true;
}


/**
 * First Initializing of the variance/oc where each value of the matrix has 
 * to be created, instead of updated variance, where it is not necessary to
 * recalculate each value/entry.
 */
void variance::InitVariance(){
	std::pair<double,uint32_t> tmp_pair;
	double var_hom, var_bac, hom, back;
	uint64_t oc_hom;
	uint32_t lower_length, upper_length, length_mean;
	int shift;
	oc_hom = 1;
	variance_val = 0;

	variance_mat = std::vector< std::vector<double> >(size,std::vector<double>(size,0));
	pattern_variance = std::vector<double>(size,0);
	pattern_current_contribute = std::vector<double>(size,0);
	pattern_last_contribute = std::vector<double>(size,0);
	for(uint32_t i = 0; i < size; i++){
		lower_length = pattern_set.GetLength(i)-1;
		for(uint32_t j = i; j < size; j++){
			upper_length = pattern_set.GetLength(j);
	
			if(i == j && !oc){
				upper_length = 1;
			}
			var_hom = 0.0;
			var_bac = 0.0;
			for (int s = -1 * lower_length; s < (int)upper_length; s++) {
				shift = ShiftPos(i, j, s);
				if (!oc) {
					var_hom += (pow(p, shift) - pow(p, 2 * weight));
					var_bac += (pow(q, shift) - pow(q, 2 * weight));
				}
				else{	
					var_hom += (double) (oc_hom << (shift));				
				}
			}
			if(!oc){
				length_mean = (lower_length + 1 + upper_length) / 2;
				hom = (seq_leng - length_mean + 1);
				hom *= var_hom;
	
				back = (seq_leng - length_mean + 1);
				back *= (seq_leng - length_mean);
				back *= var_bac;
			}
			else{
				hom = var_hom;
				back = 0.0; 
			}
			variance_mat[i][j] = hom + back;
			variance_mat[j][i] = hom + back;
			variance_val += variance_mat[i][j];
			pattern_variance[i] += variance_mat[i][j];
			pattern_variance[j] += variance_mat[i][j];
		}
	}

	for(uint32_t i = 0; i < size; i++){
		if(pattern_set.GetImprove(i)){
			tmp_pair = std::make_pair(pattern_variance[i],i);
			variance_set.insert(tmp_pair);
		}
	}
}

/**
 * Part of the destructor; could be used in other cases for 
 * reinitializing the variance/oc new.
 */
void variance::Clear(){
	variance_mat.clear();
	variance_set.clear();
	pattern_variance.clear();
	pattern_current_contribute.clear();
	pattern_last_contribute.clear();
}


/*---Variance-----------------------------------------------------------------*/
/**
 * Instead of recalculating each value, only the matrix entries which
 * belong to the modified pattern have to be recalculated. 
 * 
 * @param number		Patternindex of the modified pattern 
 */
void variance::UpdateVariance(uint32_t number){
	double var_hom, var_bac, hom, back;
	uint64_t oc_hom;
	uint32_t lower_length, upper_length, length_mean;
	int shift;

	if(initialized){
		oc_hom = 1;
		lower_length = pattern_set.GetLength(number)-1;
		for(uint32_t i = 0; i < size; i++){
			pattern_last_contribute[i] = variance_mat[number][i];
			variance_val -= pattern_last_contribute[i];
			pattern_variance[number] -= pattern_last_contribute[i];
			pattern_variance[i] -= pattern_last_contribute[i];
			upper_length = pattern_set.GetLength(i);

			if(number == i && !oc){
				upper_length = 1;
			}
			var_hom = 0.0;
			var_bac = 0.0;
			for (int s = -1 * lower_length; s < (int)upper_length; s++) {
				shift = ShiftPos(number, i, s);
				if (!oc) {
					var_hom += (pow(p, shift) - pow(p, 2 * weight));
					var_bac += (pow(q, shift) - pow(q, 2 * weight));
				}
				else{	
					var_hom += (double) (oc_hom << (shift));				
				}
			}
			if(!oc){
				length_mean = (lower_length + 1 + upper_length) / 2;
				hom = (seq_leng - length_mean + 1);
				hom *= var_hom;
	
				back = (seq_leng - length_mean + 1);
				back *= (seq_leng - length_mean);
				back *= var_bac;
			}
			else{
				hom = var_hom;
				back = 0.0; 
			} 
			variance_mat[number][i] = hom + back;
			variance_mat[i][number] = hom + back;
			pattern_current_contribute[i] = variance_mat[number][i];
			variance_val += pattern_current_contribute[i];
			pattern_variance[number] += pattern_current_contribute[i];
			pattern_variance[i] += pattern_current_contribute[i];

		}
	}
}

/**
 * Shifts pattern and counts the number of common match positions.
 *
 * @param p1 			Position of the first used pattern of the pattern set.
 *
 * @param p2			Position of the second used pattern of the pattern set.
 *							NOTE: possible is p1 = p2
 *
 * @param s				The shift of the second pattern, s < 0 := shift left 
 *							pattern 2, s > := shift right pattern 2.
 * 
 * @return 				Calculates and returns current variance.
 */
uint32_t variance::ShiftPos(uint32_t number, uint32_t number2, int shift){
	std::vector<char> stra, strb;
	uint64_t pata, patb, pos, tmp;
	uint32_t counter, maxi, maxa, maxb;
	counter = 0;
	pos = 1;
	if(bitmode){
		pata = pattern_set.GetBitPattern(number);
		patb = pattern_set.GetBitPattern(number2);
		if(shift < 0){
			shift= std::abs(shift);
			pata = pata >> shift;
		}
		else{
			patb = patb >> shift;
		}
		maxa = (uint32_t) std::log2(pata)+1;
		maxb = (uint32_t) std::log2(patb)+1;
		maxi = std::min(maxa,maxb);
		tmp = pata & patb;
		for(uint32_t i = 0; i < maxi; i++){
			if((tmp & (pos << i)) != 0){
				counter ++;
			}
		}
		if(!oc){
			counter = 2*weight-counter;
		}
	}
	else{
		stra = pattern_set.GetPattern(number);
		strb = pattern_set.GetPattern(number2);

		maxa = pattern_set.GetLength(number);
		maxb = pattern_set.GetLength(number2);
		if(shift < 0){
			shift = std::abs(shift);
			maxa -= shift;
		}
		else{
			maxb -= shift;
			std::swap(stra,strb);
		}
		maxi = std::min(maxa,maxb);
		for(uint32_t i = 0; i < maxi; i++){
			if(stra[i] == '1' && strb[strb.size() - maxi + i] == '1'){
				counter++;
			}
		}
	}
	return counter;
}


/**
 * If the random swap leads to a worse oc/variance, the old values have
 * to be restored.
 */
void variance::ResetContribute(uint32_t number){
	for(uint32_t i = 0; i < size; i++){
		variance_val -= pattern_current_contribute[i];
		pattern_variance[number] -= pattern_current_contribute[i];
		pattern_variance[i] -= pattern_current_contribute[i];

		variance_val += pattern_last_contribute[i];
		pattern_variance[number] += pattern_last_contribute[i];
		pattern_variance[i] += pattern_last_contribute[i];
		variance_mat[number][i] = pattern_last_contribute[i];
		variance_mat[i][number] = pattern_last_contribute[i];
	}
}


/*---Improvment--------------------------------------------------------------*/
/**
 * The improvement method.
 *
 * @param limits		The number of patterns which have to be selected and mayby modified.
 */
void variance::Improve(uint32_t limits){
	std::set<std::pair<double,uint32_t> >::iterator var_iterator;
	std::vector<char> last_pat;
	std::pair<double,uint32_t> tmp_pair;
	uint64_t last_bit_pat;
	double best_variance_val;
	uint32_t worst_pat, better_pattern;

	if(initialized && improve){
		best_variance_val = variance_val;
		var_iterator = variance_set.end();
		var_iterator--;
		better_pattern = 0;

		if(bitmode){
			last_bit_pat = 0;
		}

		for(uint32_t i = 1; i <= limits; i++){
			worst_pat = var_iterator->second;
			if(bitmode){
				last_bit_pat = pattern_set.GetBitPattern(worst_pat);
			}
			else{
				last_pat = pattern_set.GetPattern(worst_pat);
			}

			pattern_set.ChangeBitsRandom(worst_pat);
			UpdateVariance(worst_pat);

			if(variance_val <= best_variance_val){
				if(best_variance_val == variance_val){
					if(var_iterator == variance_set.begin()){
						var_iterator = variance_set.end();
					}
				}
				else{
					best_variance_val = variance_val;
					variance_set.erase(var_iterator);
					tmp_pair = std::make_pair(pattern_variance[worst_pat],worst_pat);
					variance_set.insert(tmp_pair);
					var_iterator = variance_set.end();
					better_pattern++;
					
					if(loops == 1 && (!silent || impro_silent)){
						if(!quiet){
							std::cout << "*** BETTER PATTERN " << better_pattern << " *** \t(random permutation)" << std::endl;
							std::cout << "Step " << i << " / " << limits << std::endl << "Patternset: \n";
							pattern_set.Print();
							std::cout << outvar << GetVariance() << std::endl;
							std::cout << "norm_" << outvar << GetNormVariance() << std::endl << std::endl;
						}
						else{
							std::cout << "\r*** BETTER PATTERN " << better_pattern << " *** \t(random permutation)";
							std::cout.flush();
						}
					}
				}
				if(!bitmode){
					last_pat.clear();
				}
			}
			else{
				if(var_iterator == variance_set.begin()){
					var_iterator = variance_set.end();
				}
				if(bitmode){
					pattern_set.SetPattern(worst_pat, last_bit_pat);
				}
				else{
					pattern_set.SetPattern(worst_pat, last_pat);
				}
				ResetContribute(worst_pat);
			}
			var_iterator--;
		}
		if(!silent && quiet){
			std::cout << std::endl;
		}
	}
	else if(!improve && !silent){
		SecureMessage("noimprove",-1);
	}
}

/**
 * Improvement with reinitialization. Improves a patternset, resets it 'loops' times and
 * and sets the best patternset with the best detected oc/variance.
 *
 * @param limits		The number of patterns which have to be selected and mayby modified.
 *
 * @param loops			The number of reinitializations to find better oc/variance patternsets.
 */
void variance::Improve(uint32_t limits, uint32_t loops){
	if(initialized && improve){
		patternset best_pattern_set;
		std::ofstream pattern_out;
		std::string out_str;
		double best_variance_val;
		uint32_t better_pattern;

		this->loops = loops;

		if(!silent || best_silent){
			std::cout << "\n===== First patternset =====" << std::endl;
			pattern_set.Print();
			std::cout << outvar << GetVariance() << std::endl;
			std::cout << "norm_" << outvar << GetNormVariance() << std::endl << std::endl;
		}

		best_variance_val = variance_val;
		better_pattern = 0;
		best_pattern_set = pattern_set;
		for(uint32_t i = 1; i <= loops; i++){
			Improve(limits);
			if(best_variance_val > variance_val){
				better_pattern++;
				 if(loops != 1 && (!silent || impro_silent)){
				 		if(!quiet){
				 			std::cout << "\n*** BETTER PATTERN " << better_pattern << " *** \t(oc/variance optimization)" << std::endl;
				 			std::cout << "Step " << i << " / " << loops << std::endl << "Patternset: \n";
				 			pattern_set.Print();
				 			std::cout << outvar << GetVariance() << std::endl;
				 			std::cout << "norm_" << outvar << GetNormVariance() << std::endl << std::endl;
				 		}
				 		else{
				 			std::cout << "\r*** BETTER PATTERN " << better_pattern << " *** \t(oc/variance optimization)";
				 			std::cout.flush();
				 		}
				 	}
				best_variance_val = variance_val;
				best_pattern_set = pattern_set;	
			}
			ReInit();
		}
		if(!silent || best_silent){
			std::cout << std::endl;
		}
		pattern_set = best_pattern_set;
		ReInit(pattern_set);

		if(!silent || best_silent){
			std::cout << "\n===== Best patternset ======" << std::endl;
			pattern_set.Print();
			std::cout << outvar << GetVariance() << std::endl;
			std::cout << "norm_" << outvar << GetNormVariance() << std::endl << std::endl;
		}
		if(isOutFile && best_silent){
			pattern_out.open(outFile);
			for(uint32_t i = 0; i < size; i++){
				out_str = std::string(pattern_set.GetPattern(i).begin(),pattern_set.GetPattern(i).end());
				pattern_out << out_str << std::endl;
			}
			pattern_out << GetFormat() << GetVariance() << std::endl;
			pattern_out << "norm_" << GetFormat() << GetNormVariance() << std::endl;
		}
	}
	else if(!improve && !silent){
		SecureMessage("noimprove",-1);
	}
}


/*---stuff-------------------------------------------------------------------*/
/**
 * Prints the current patternset to the commandline.
 */
void variance::PrintPatternSet(){
	pattern_set.Print();
}


/**
 * A Method to collect all errormessages. Just easier for programmer to change
 *  	the text or extend.
 *
 * @param errmsg			Due to a few possible errormessages, this is the option, which has to be printed.
 *
 * @param pos			The position of the incorrect patterns.
 */
void variance::SecureMessage(std::string errmsg, int pos){
	if (errmsg == "noimprove" && !silent) {
		printf("%c[1;33m", 27);
		std::cerr << "Using your pattern conditions it is not sensible to improve your pattern, sorry!" << std::endl;
		std::cerr << "Deactivating improve mode\n" << std::endl;
		printf("%c[0m", 27);
		return;
	}
	if (errmsg == "pq") {
		printf("%c[1;32m\n\n--> IMPORTANT <--\n", 27);
		printf("%c[0m", 27);
		printf("%c[1;31mError ", 27);
		printf("%c[0m", 27);
		std::cerr << "while parsing your p and/or q value: \t0 < q <= p <= 1!" << std::endl;
		std::cerr << "Return to default values:\tp = 0.75 \tq=0.25\n" << std::endl;
		return;
	}
}


/*---Set&Get-----------------------------------------------------------------*/
/**
 * Returns the patternset object
 * @return 				Copy of the current patternset object (type: patternset)
 */
patternset variance::GetPatternSet(){
	return pattern_set;
}


/**
 * Returns the patternset as std::vector< std::vector<char> >
 *
 * @return 				Copy of the current patternset vector (type: std::vector< std::vector<char> >)
 */
std::vector< std::vector<char> > variance::GetPattern(){
	return pattern_set.GetPattern();
}

/**
 * Returns a pattern of the patternset as std::vector<char>
 *
 * @return 				Copy of the 'number'th pattern of the current patternset vector (type: std::vector<char>)
 */
std::vector<char> variance::GetPattern(uint32_t number){
	return pattern_set.GetPattern(number);
}


/**
 * Returns the patternset as bitpattern as std::vector<uint64_t>
 *
 * @return 				Copy of the current bit-patternset vector (type: std::vector<uint64_t>)
 */
std::vector<uint64_t> variance::GetBitPattern(){
	return pattern_set.GetBitPattern();
}


/**
 * Returns a pattern of the patternset as uint64_t.
 *
 * @return 				Copy of the 'number'th pattern of the current patternset vector (type: uint64_t)
 */
uint64_t variance::GetBitPattern(uint32_t number){
	return pattern_set.GetBitPattern(number);
}


/**
 * Returns the current variance/oc value.
 *
 * @return 				Current pattern variance/oc.
 */
double variance::GetVariance(){
	return variance_val;
}


/**
 * Returns the current norm variance/oc value.
 *
 * @return 				Current pattern norm variance/oc.
 */
double variance::GetNormVariance(){
	return variance_val/norm_size;
}


/**
 * Returns the current patternset size.
 *
 * @return 				Current pattern size.
 */
uint32_t variance::GetSize(){
	return size;
}


/**
 * Returns the current patternset weight.
 *
 * @return 				Current pattern weight.
 */
uint32_t variance::GetWeight(){
	return weight;
}


/**
 * Returns the current patternset minimum don't care number.
 *
 * @return 				Current pattern minimum don't care number.
 */
uint32_t variance::GetMaxDontcare(){
	return max_dontcare;
}


/**
 * Returns the current patternset maximum don't care number.
 *
 * @return 				Current pattern maximum don't care number.
 */
uint32_t variance::GetMinDontcare(){
	return min_dontcare;
}


/**
 * Returns the current dataset sequence length for variance calculation.
 *
 * @return 				Current dataset sequence length for variance calculation.
 */
uint32_t variance::GetSequenceLength(){
	return seq_leng;
}


/**
 * Returns the current match probability for variance calculation.
 *
 * @return 				Current match probability for variance calculation.
 */
double variance::GetP(){
	return p;
}


/**
 * Returns the current background probability for variance calculation.
 *
 * @return 				Current background probability for variance calculation.
 */
double variance::GetQ(){
	return q;
}


/**
 * Returns boolean if the current patternset modus is bitmode oder not.
 *
 * @return 				Current bitmode boolean (true, if bitmode).
 */
bool variance::GetBitMode(){
	return pattern_set.GetBitMode();
}


/**
 * Returns boolean if the current patternset modus has been updated
 *
 * @return 				Current patternmode has been updated.
 */
bool variance::GetUpdate(){
	return pattern_set.GetUpdate();
}

/**
 * A special option, if silent == true but, first and best pattern should be printed.
 *
 * @param sil 			Boolean, if output first and best pattern, although silent == true.
 */
void variance::BestSilent(bool sil){
	best_silent = sil;
}


/**
 * A special option, if silent == true but the number of improved pattern should be printed.
 *
 * @param sil 			Boolean, if print number of improvements, although silent == true.
 */
void variance::ImproSilent(bool sil){
	impro_silent = sil;
}

/**
 * Returns a string containing a speacial output format, if its overlap complexity or not.
 *
 * @return 				String which contains current pattern mode overlap complexity or variance calculation. 
 */
std::string variance::GetFormat(){
	return outvar;
}


/**
 * Returns a boolean, if there is at least one pattern, which can be improved/changed, or not.
 *
 * @return 				Boolean, if at least one pattern can be changed.
 */
bool variance::GetImprove(){
	return improve;
}


/**
 * Prints the current patternset to the commandline.
 */
void variance::Print(){
	pattern_set.Print();
}