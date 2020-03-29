#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
#include "common.h"
#ifndef __NVCC__
	#include "../xavier/xavier.h"
#endif
#include <omp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

#ifdef __NVCC__
#include "../loganGPU/logan.cuh"
#endif

#define BATCH_SIZE 30000

using namespace seqan;
using namespace std;

char
complementbase(char n) {
	switch(n)
	{
	case 'A':
		return 'T';
	case 'T':
		return 'A';
	case 'G':
		return 'C';
	case 'C':
		return 'G';
	}
	assert(false);
	return ' ';
}

std::string
reversecomplement(const std::string& seq) {

	std::string cpyseq = seq;
	std::reverse(cpyseq.begin(), cpyseq.end());

	std::transform(
		std::begin(cpyseq),
		std::end  (cpyseq),
		std::begin(cpyseq),
	complementbase);

	return cpyseq;
}

double slope(double error)
{
	double p_mat = pow(1-error,2);  // match
	double p_mis = 1-p_mat;         // mismatch/gap
	double alpha = 1;               // match penalty
	double beta  = 1;                // mismatch/gap penalty

	return alpha*p_mat - beta*p_mis;
}

/**
 * @brief alignSeqAn does the seed-and-extend alignment
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param i is the starting position of the k-mer on the first read
 * @param j is the starting position of the k-mer on the second read
 * @param xDrop
 * @return alignment score and extended seed
 */
seqAnResult alignSeqAn(const std::string & row, const std::string & col, int rlen, int i, int j, int xDrop, int kmerSize) {
	// printLog("SeqAn");
	Score<int, Simple> scoringScheme(1,-1,-1);

	Dna5String seqH(row); 
	Dna5String seqV(col); 
	Dna5String seedH;
	Dna5String seedV;
	string strand;
	int longestExtensionTemp;
	seqAnResult longestExtensionScore;


	TSeed seed(i, j, i + kmerSize, j + kmerSize);
	seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
	seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));

	/* we are reversing the "row", "col" is always on the forward strand */
	Dna5StringReverseComplement twin(seedH);

	if(twin == seedV)
	{
		strand = 'c';
		Dna5StringReverseComplement twinRead(seqH);
		i = rlen - i - kmerSize;

		setBeginPositionH(seed, i);
		setBeginPositionV(seed, j);
		setEndPositionH(seed, i + kmerSize);
		setEndPositionV(seed, j + kmerSize);

		/* Perform match extension */
		longestExtensionTemp = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, xDrop, kmerSize, GappedXDrop());

	} else
	{
		strand = 'n';
		longestExtensionTemp = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, xDrop,  kmerSize, GappedXDrop());
	} 

	longestExtensionScore.score = longestExtensionTemp;
	longestExtensionScore.seed = seed;
	longestExtensionScore.strand = strand;
	return longestExtensionScore;
}

//////////////////////////////////////////
//	A:				//	B:				//
//	>------->		//	<-------<		//
//		<-------<	//		>------->	//
//////////////////////////////////////////
//	D:				//	E:				//
//		>------->	//		<-------<	//
//	<-------<		//	>------->		//
//////////////////////////////////////////
//	C:				//	F:				//
//	>------->		//	<-------<		//
//		>------->	//		<-------<	//
//////////////////////////////////////////
//	G:				//	H:				//
//		>------->	//		<-------<	//
//	>------->		//	<-------<		//
//////////////////////////////////////////

#ifndef __NVCC__
/**
 * @brief alignLogan does the seed-and-extend alignment
 * @param row
 * @param col
 * @param rowLen is the length of the row sequence
 * @param i is the starting position of the k-mer on the first read
 * @param j is the starting position of the k-mer on the second read
 * @param xDrop
 * @return alignment score and extended seed
 */
xavierResult xavierAlign(const std::string& row, const std::string& col, int rowLen, PairType read_i, PairType read_j, int xDrop, int kmerSize)
{
	// GGGG: result.first = best score, result.second = exit score when (if) x-drop termination is satified
	std::pair<int, int> tmp;
	xavierResult result;

	short match    =  1;
	short mismatch = -1;
	short gap 	   = -1;

	// GGGG: initialize scoring scheme
	ScoringSchemeX scoringScheme(match, mismatch, gap);	// enalties (LOGAN currently supports only linear gap penalty and penalty within +/- 3)

	std::string read_rc;
	SeedX seed(read_i.first, read_j.first, kmerSize);

	// GGGG: reads are on opposite strand
	// GGGG: for miniasm I need only one option but for combter I might use both?
	if(read_i.second != read_j.second)
	{
		//	* B: >---< | j.E ---> i.B | ~A
		//	* E: <---> | i.B ---> j.E | ~D
		if(read_i.second)
		{
			read_rc = row;
			std::reverse  (std::begin(read_rc), std::end(read_rc));
			std::transform(std::begin(read_rc), std::end(read_rc), std::begin(read_rc), complementbase);

			setBeginPositionH(seed, rowLen - read_i.first - kmerSize);
			setEndPositionH  (seed, rowLen - read_i.first);

			tmp = XavierXDrop(seed, XAVIER_EXTEND_BOTH, read_rc, col, scoringScheme, xDrop);

			result.type = "B|E";
		}
		//	* A: <---> | i.E ---> j.B | ~B
		//	* D: >---< | j.B ---> i.E | ~E
		else if(read_j.second)
		{
			read_rc = col;
			std::reverse  (std::begin(read_rc), std::end(read_rc));
			std::transform(std::begin(read_rc), std::end(read_rc), std::begin(read_rc), complementbase);

			setBeginPositionH(seed, rowLen - read_j.first - kmerSize);
			setEndPositionH  (seed, rowLen - read_j.first);

			tmp = XavierXDrop(seed, XAVIER_EXTEND_BOTH, row, read_rc, scoringScheme, xDrop);
			result.type = "A|D";			
		}
		// GGGG: DO NOT output contained read for combter for now
		else
		{
			std::cout << "err" << std::endl;
		}	
	}
	else
	{
		//	* C: >--> | j.E --> i.E | ~F
		//	* G: >--> | i.E --> j.E | ~H
		if(!read_i.second)
		{
			tmp = XavierXDrop(seed, XAVIER_EXTEND_BOTH, row, col, scoringScheme, xDrop);
			result.type = "C|G";
		}
		//	* F: <--< | i.B --> j.B	| ~C
		//	* H: <--< | j.B --> i.B	| ~G
		if(read_i.second)
		{
			tmp = XavierXDrop(seed, XAVIER_EXTEND_BOTH, row, col, scoringScheme, xDrop);
			result.type = "F|H";
		}
		// GGGG: DO NOT output contained read for combter for now
		else
		{
			std::cout << "err" << std::endl;
		}
	}

	result.score = tmp.first; 	// best score

	setBeginPositionH(result.seed, getBeginPositionH(seed));	// updated extension
	setBeginPositionV(result.seed, getBeginPositionV(seed));	// updated extension

	setEndPositionH(result.seed, getEndPositionH(seed));		// updated extension
	setEndPositionV(result.seed, getEndPositionV(seed));		// updated extension

	return result;
}

#else

// ======================================= //
// 				GPU Functions			   //
// ======================================= //

void alignLogan(vector<string>&	target, vector<string>&	query, vector<SeedL>& seeds, 
	const BELLApars& b_pars, vector<loganResult>& longestExtensionScore)
{

	ScoringSchemeL sscheme(1, -1, -1, -1);
	std::vector<ScoringSchemeL> scoring;
	scoring.push_back(sscheme);

	int AlignmentsToBePerformed = seeds.size();
	printLog(AlignmentsToBePerformed);
	//int* res = (int*)malloc(BATCH_SIZE*sizeof(int));
	int numAlignmentsLocal = BATCH_SIZE * b_pars.numGPU; 
	cout <<"///////////////////////////////////////////////" <<b_pars.numGPU << endl;
	

	//	Divide the alignment in batches of 100K alignments
	for(int i = 0; i < AlignmentsToBePerformed; i += BATCH_SIZE * b_pars.numGPU)
	{
		if(AlignmentsToBePerformed < (i + BATCH_SIZE * b_pars.numGPU))
			numAlignmentsLocal = AlignmentsToBePerformed % (BATCH_SIZE * b_pars.numGPU);

		int* res = (int*)malloc(numAlignmentsLocal * sizeof(int));	

		std::vector<string>::const_iterator first_t = target.begin() + i;
		std::vector<string>::const_iterator last_t  = target.begin() + i + numAlignmentsLocal;
		std::vector<string> target_b(first_t, last_t);

		std::vector<string>::const_iterator first_q = query.begin() + i;
		std::vector<string>::const_iterator last_q  = query.begin() + i + numAlignmentsLocal;
		std::vector<string> query_b(first_q, last_q);

		std::vector<SeedL>::const_iterator first_s = seeds.begin() + i;
		std::vector<SeedL>::const_iterator last_s  = seeds.begin() + i + numAlignmentsLocal;
		std::vector<SeedL> seeds_b(first_s, last_s);

		extendSeedL(seeds_b, EXTEND_BOTHL, target_b, query_b, scoring, b_pars.xDrop, b_pars.kmerSize, res, numAlignmentsLocal, b_pars.numGPU);

		for(int j=0; j<numAlignmentsLocal; j++)
		{
			longestExtensionScore[j+i].score = res[j];
			longestExtensionScore[j+i].seed = seeds_b[j];
		}

		free(res);
	}
}
#endif // __NVCC__
#endif // _ALIGNMENT_H_
