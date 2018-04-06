/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * variance object header
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
#ifndef VARIANCE_H_
#define VARIANCE_H_

#include <cmath>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>
#include <string>
#include "patternset.h"

class variance{
public:
	variance();
	variance(uint32_t size, uint32_t weight, uint32_t min_dontcare, uint32_t max_dontcare);
	variance(uint32_t size, uint32_t weight, uint32_t min_dontcare, uint32_t max_dontcare, uint32_t seq_leng, double p, double q);
	variance(const char *inFile, double p, double q, uint32_t seq_leng);
	~variance();

	void Init();
	void Init(bool oc, bool improve, bool quiet, bool silent, bool random_leng, const char *outFile);
	void ReInit();

	void Improve(uint32_t limits);
	void Improve(uint32_t limits, uint32_t loops);
	void PrintPatternSet();

	/*-----------------------------------------------------------------------*/
	patternset GetPatternSet();
	std::vector< std::vector<char> > GetPattern();
	std::vector<char> GetPattern(uint32_t number);
	std::vector<uint64_t> GetBitPattern();
	uint64_t GetBitPattern(uint32_t number);
	
	double GetVariance();
	double GetNormVariance();

	uint32_t GetSize();
	uint32_t GetWeight();
	uint32_t GetMaxDontcare();
	uint32_t GetMinDontcare();
	uint32_t GetSequenceLength();
	double GetP();
	double GetQ();
	bool GetBitMode();
	bool GetUpdate();
	void BestSilent(bool sil);
	void ImproSilent(bool sil);
	std::string GetFormat();
	bool GetImprove();
	void Print();

protected:
	void Clear();
	void ReInit(patternset &pat);
	void InitVariance();
	void SetVariance();
	void UpdateVariance(uint32_t number);
	uint32_t ShiftPos(uint32_t number, uint32_t number2, int shift);
	void ResetContribute(uint32_t number);
	void SecureMessage(std::string errmsg, int pos);


private:
	patternset pattern_set;
	std::vector< std::vector<double> > variance_mat;
	std::set< std::pair<double, uint32_t> > variance_set;
	std::vector<double> pattern_variance;
	std::vector<double> pattern_current_contribute;
	std::vector<double> pattern_last_contribute;

	std::string outvar;
	double variance_val;
	double p;
	double q;

	uint32_t size;
	uint32_t weight;
	uint32_t min_dontcare;
	uint32_t max_dontcare;
	uint32_t seq_leng;
	uint32_t norm_size;
	uint32_t loops;

	const char *inFile;
	const char *outFile;

	bool isInFile;
	bool isOutFile;
	bool bitmode;
	bool initialized;
	bool random_leng;
	bool improve;
	bool quiet;
	bool silent;
	bool best_silent;
	bool impro_silent;
	bool oc;
};
#endif