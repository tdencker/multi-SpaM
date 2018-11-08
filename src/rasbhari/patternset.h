/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * patternset object header
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
#ifndef PATTERNSET_H_
#define PATTERNSET_H_

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>

void setRasbhariSeed(uint64_t seed);

class patternset{
public:
	patternset();
	patternset(uint32_t size, uint32_t weight, uint32_t min_dontcare, uint32_t max_dontcare, bool silent);
	patternset(const char *patternfile, bool silent);
	~patternset();

	void Init(bool random_length);
	void ReInit();

	void ChangeBitsRandom(uint32_t number);
	bool ChangeBitPos(uint32_t number, uint64_t zero, uint64_t one);

	void Print();

	static uint64_t StringToBit(std::vector<char> pattern);
	static std::vector<char> BitToString(uint64_t bitpat);

	/*-----------------------------------------------------------------------*/
	std::vector<std::vector<char> > GetPattern();
	std::vector<uint64_t> GetBitPattern();
	std::vector<char> GetPattern(uint32_t number);
	uint64_t GetBitPattern(uint32_t number);

	std::vector<uint32_t> GetDontCareVec();
	bool GetImprove(uint32_t number);
	bool GetImprove();
	uint32_t GetSize();
	uint32_t GetWeight();
	uint32_t GetLength(uint32_t number);
	bool GetBitMode();	
	bool GetUpdate();

	void SetPattern(uint32_t number, uint64_t pat);
	void SetPattern(uint32_t number, std::vector<char> pat);

protected:
	void CheckParams();
	std::vector<uint32_t> CreateLengths();
	std::vector<std::vector<char> > CreatePattern();
	bool ReadPattern();
	bool VerifyPatternCondition();
	bool SetImprove();
	void PatternToBit();
	void BitToPattern();

	uint64_t GetSymbolRandPos(uint32_t number, char symb);

	bool UniqPattern(std::vector<std::vector<char> > pats, std::vector<char> pat);
	bool UniqBit(uint64_t bit_pat);

	static uint32_t GetWeight(std::vector<char> pattern);
	static bool IsBinary(std::vector<char> pattern);

	static double MaxNumberPattern(uint32_t k, uint32_t l);
	static double Faculty(uint32_t value);

	void SecureMessage(std::string errmsg, int pos);

private:
	std::vector<std::vector<char> > pattern_set;
	std::vector<uint64_t> bit_pattern_set;
	std::vector<uint32_t> pattern_dontcare;
	std::vector<bool> pattern_improve;
	std::vector<double> pattern_uniq;

	uint32_t size;
	uint32_t weight;
	uint32_t min_dontcare;
	uint32_t max_dontcare;

	const char *inFile;

	bool bitmode;
	bool random_length;
	bool isInFile;
	bool improve;
	bool update;
	bool silent;
	bool initialized;
};
#endif