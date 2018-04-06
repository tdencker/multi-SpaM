/**
 * This programm calculates the variance/OC and/or the sensitivity of a set of pattern with the same weight.
 * It is possible to improve your patternset and read patterns from a file.
 *
 * sensitivity object header
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
#ifndef SENSITIVITY_H_
#define SENSITIVITY_H_

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "variance.h"

class sensitivity{
public:
	sensitivity();
	sensitivity(const char *inFile, double p, double q, uint32_t H, int seq_leng);
	sensitivity(uint32_t size, uint32_t weight,uint32_t min_dontcare, uint32_t max_dontcare, double p, double q, double H, uint32_t seq_leng);

	void Init();
	void Init(bool sens, bool oc, bool improve, bool randpatleng, bool quiet, bool silent, const char *outFile);
	void ReInit();

	void Improve(uint32_t limits, uint32_t opt_oc);
	void Improve(uint32_t limits, uint32_t opt_oc, uint32_t opt_sens);

	void Print();
	std::vector< std::vector<char> > GetPattern();
	double GetSensitivity();

	uint32_t GetSize();
	uint32_t GetWeight();
	uint32_t GetMaxDontcare();
	uint32_t GetMinDontcare();
	uint32_t GetSequenceLength();
	uint32_t GetH();
	double GetP();
	double GetQ();

	double CalculateSensitivity(std::vector< std::vector<char> > pattern, double p, int n, int R);

protected:
	void ReInit(variance &var_pattern);

	inline long long BIN_REVERSED_TO_INT2(char *s);
	double MULTIPLE_SENSITIVITY2(char** SEEDS, int NO_SEEDS, long long N, double P);
	void SecureMessage(std::string errmsg, int pos);

private:
	variance var_pattern;
	double sens_value;
	double p;
	double q;
	uint32_t size;
	uint32_t weight;
	uint32_t max_dontcare;
	uint32_t min_dontcare;
	uint32_t H;
	uint32_t seq_leng;

	uint32_t opt_sens;
	uint32_t opt_oc;

	const char *inFile;
	const char *outFile;

	bool sens;
	bool oc;
	bool improve;
	bool silent;
	bool quiet;
	bool randpatleng;
	bool isInFile;
	bool isOutFile;
	bool init;
	bool bitmode;
	bool update;
};
#endif