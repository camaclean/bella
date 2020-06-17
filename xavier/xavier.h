//===========================================================================
// Title:  Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment
// Author: G. Guidi, E. Younis
// Date:   29 April 2019
//===========================================================================

#ifndef XAVIER_H
#define XAVIER_H

#include<vector>
#include<iostream>
#include<omp.h>
#include<algorithm>
#include<inttypes.h>
#include<assert.h>
#include<iterator>
#include<x86intrin.h>
#include"simdutils.h"


vectorUnionType left_shift(const vectorType& v)
{
	signed char tmp[33];
	_mm256_storeu_si256((__m256i*) &tmp[0], v);
	vectorUnionType ret;
	ret.simd = _mm256_lddqu_si256((__m256i*) &tmp[1]);
	return ret;
}


void
XavierPhase1(XavierState& state)
{
	myLog("Phase1");
	// we need one more space for the off-grid values and one more space for antiDiag2
	int DPmatrix[LOGICALWIDTH + 2][LOGICALWIDTH + 2];

	// DPmatrix initialization
	DPmatrix[0][0] = 0;
	for(int i = 1; i < LOGICALWIDTH + 2; i++ )
	{
		DPmatrix[0][i] = -i;
		DPmatrix[i][0] = -i;
	}

	// DPmax tracks maximum value in DPmatrix for xdrop condition
	int DPmax = 0;

	// Dynamic programming loop to fill DPmatrix
	for(int i = 1; i < LOGICALWIDTH + 2; i++)
	{
		for (int j = 1; j <= LOGICALWIDTH + 2 - i; j++) // we only need the upper-left triangular matrix
		{
			int oneF = DPmatrix[i-1][j-1];

			// Comparing bases
			if(state.queryh[i-1] == state.queryv[j-1])
				oneF += state.get_match_cost();
			else
				oneF += state.get_mismatch_cost();

			int twoF = std::max(DPmatrix[i-1][j], DPmatrix[i][j-1]);
			twoF += state.get_gap_cost();

			DPmatrix[i][j] = std::max(oneF, twoF);

			// Heuristic to keep track of the max in phase1
			if(DPmatrix[i][j] > DPmax)
				DPmax = DPmatrix[i][j];
		}
	}

	for(int i = 0; i < LOGICALWIDTH; ++i)
	{
		state.update_vqueryh(i, state.queryh[i + 1]);
		state.update_vqueryv(i, state.queryv[LOGICALWIDTH - i]);
	}

	state.update_vqueryh(LOGICALWIDTH, NINF);
	state.update_vqueryv(LOGICALWIDTH, NINF);

	int antiDiagMax = std::numeric_limits<int8_t>::min();

	// Load DPmatrix into antiDiag1 and antiDiag2 vector and find max elem at the end of phase1 in antiDiag1
	for(int i = 1; i < LOGICALWIDTH + 1; ++i)
	{
		int value1 = DPmatrix[i][LOGICALWIDTH - i + 1];
		int value2 = DPmatrix[i + 1][LOGICALWIDTH - i + 1];

		state.update_antiDiag1(i - 1, value1);
		state.update_antiDiag2(i, value2);

		if(value1 > antiDiagMax)
			antiDiagMax = value1;
	}

	state.update_antiDiag1(LOGICALWIDTH, NINF);
	state.update_antiDiag2(0, NINF);
	printf("antidiag1\n");
	printVectorD(state.antiDiag1.simd);
	printf("antidiag2\n");
	printVectorD(state.antiDiag2.simd);
	//state.antiDiag1 = shiftLeft(state.antiDiag1.simd);
	/*
	{
		auto a = state.antiDiag1;
		auto start = std::chrono::system_clock::now();
		for (size_t i = 0; i < 100000000000; ++i);
			a = shiftLeft(a.simd);
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> sec = end-start;
		std::cout << "Microseconds: " << std::chrono::duration_cast<std::chrono::microseconds>(sec).count() << std::endl;
	}
	{
		auto a = state.antiDiag1;
		auto start = std::chrono::system_clock::now();
		for (size_t i = 0; i < 100000000000; ++i);
			a = left_shift(a.simd);
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> sec = end-start;
		std::cout << "Microseconds: " << std::chrono::duration_cast<std::chrono::microseconds>(sec).count() << std::endl;
	}
	*/
	//auto b = shiftLeft(state.antiDiag2.simd);
	//printVectorD(a.simd);
	//printVectorD(b.simd);

	state.broadcast_antiDiag3(NINF);

	state.set_best_score(DPmax);
	state.set_curr_score(antiDiagMax);

	if(antiDiagMax < DPmax - state.get_score_dropoff())
	{
		state.xDropCond = true;

		setEndPositionH(state.seed, state.hoffset);
		setEndPositionV(state.seed, state.voffset);

		return;
	}
	// antiDiag2 going right, first computation of antiDiag3 is going down
}

void
XavierPhase2(XavierState& state)
{
	myLog("Phase2");
	int count = 0;
	while(state.hoffset < state.hlength && state.voffset < state.vlength)
	{
		printf("-------- %03d --------\n",count++);
		printVectorC(state.get_vqueryh());
		printVectorC(state.get_vqueryv());
		printVectorD(state.get_antiDiag1());
		printVectorD(state.get_antiDiag2());
		// antiDiag1F (final)
		// NOTE: -1 for a match and 0 for a mismatch
		vectorType match = cmpeqOp(state.get_vqueryh(), state.get_vqueryv());
		match = blendvOp(state.get_vmismatchCost(), state.get_vmatchCost(), match);
		printVectorD(match);
		vectorType antiDiag1F = addOp(match, state.get_antiDiag1());
		printVectorD(antiDiag1F);

		// antiDiag2S (shift)
		vectorUnionType antiDiag2S = shiftLeft(state.get_antiDiag2());
		printf("tmp: ");
		printVectorD(antiDiag2S.simd);

		// antiDiag2M (pairwise max)
		vectorType antiDiag2M = maxOp(antiDiag2S.simd, state.get_antiDiag2());

		// antiDiag2F (final)
		vectorType antiDiag2F = addOp(antiDiag2M, state.get_vgapCost());

		// Compute antiDiag3
		state.set_antiDiag3(maxOp(antiDiag1F, antiDiag2F));

		// we need to have always antiDiag3 left-aligned
		state.update_antiDiag3(LOGICALWIDTH, NINF);

		// TODO: x-drop termination
		// Note: Don't need to check x-drop every time
		// Create custom max_element that also returns position to save computation
		int8_t antiDiagBest = *std::max_element(state.antiDiag3.elem, state.antiDiag3.elem + VECTORWIDTH);
		state.set_curr_score(antiDiagBest + state.get_score_offset());

		int64_t scoreThreshold = state.get_best_score() - state.get_score_dropoff();
		if (state.get_curr_score() < scoreThreshold)
		{
			state.xDropCond = true;

			setBeginPositionH(state.seed, 0);
			setBeginPositionV(state.seed, 0);

			setEndPositionH(state.seed, state.hoffset);
			setEndPositionV(state.seed, state.voffset);

			return; // GG: it's a void function and the values are saved in XavierState object
		}

		if (antiDiagBest > CUTOFF)
		{
			int8_t min = *std::min_element(state.antiDiag3.elem, state.antiDiag3.elem + LOGICALWIDTH);
			state.set_antiDiag2(subOp(state.get_antiDiag2(), setOp(min)));
			state.set_antiDiag3(subOp(state.get_antiDiag3(), setOp(min)));
			state.set_score_offset(state.get_score_offset() + min);
		}

		// Update best
		if (state.get_curr_score() > state.get_best_score())
			state.set_best_score(state.get_curr_score());

		// TODO: optimize this
		int maxpos, max = 0;
		for (int i = 0; i < VECTORWIDTH; ++i)
		{
			if (state.antiDiag3.elem[i] > max)
			{
				maxpos = i;
				max = state.antiDiag3.elem[i];
			}
		}

		setEndPositionH(state.seed, state.hoffset);
		setEndPositionV(state.seed, state.voffset);

		if (maxpos > MIDDLE) {
			printf("moving right %d\n", maxpos);
			state.moveRight();
		}
		else
		{
			printf("moving down %d\n", maxpos);
			state.moveDown();
		}
		printVectorD(state.get_antiDiag1());
		printVectorD(state.get_antiDiag2());
	}
}

void
XavierPhase4(XavierState& state)
{
	myLog("Phase4");
	int dir = state.hoffset >= state.hlength ? goDOWN : goRIGHT;

	for (int i = 0; i < (LOGICALWIDTH - 3); i++)
	{
		// antiDiag1F (final)
		// NOTE: -1 for a match and 0 for a mismatch
		vectorType match = cmpeqOp(state.get_vqueryh(), state.get_vqueryv());
		match = blendvOp(state.get_vmismatchCost(), state.get_vmatchCost(), match);
		vectorType antiDiag1F = addOp(match, state.get_antiDiag1());

		// antiDiag2S (shift)
		vectorUnionType antiDiag2S = shiftLeft(state.get_antiDiag2());

		// antiDiag2M (pairwise max)
		vectorType antiDiag2M = maxOp(antiDiag2S.simd, state.get_antiDiag2());

		// antiDiag2F (final)
		vectorType antiDiag2F = addOp(antiDiag2M, state.get_vgapCost());

		// Compute antiDiag3
		state.set_antiDiag3(maxOp(antiDiag1F, antiDiag2F));

		// we need to have always antiDiag3 left-aligned
		state.update_antiDiag3(LOGICALWIDTH, NINF);

		// TODO: x-drop termination
		// note: Don't need to check x drop every time
		// Create custom max_element that also returns position to save computation
		int8_t antiDiagBest = *std::max_element(state.antiDiag3.elem, state.antiDiag3.elem + VECTORWIDTH);
		state.set_curr_score(antiDiagBest + state.get_score_offset());

		int64_t scoreThreshold = state.get_best_score() - state.get_score_dropoff();

		if (state.get_curr_score() < scoreThreshold)
		{
			state.xDropCond = true;
			return; // GG: it's a void function and the values are saved in XavierState object
		}

		if (antiDiagBest > CUTOFF)
		{
			int8_t min = *std::min_element(state.antiDiag3.elem, state.antiDiag3.elem + LOGICALWIDTH);
			state.set_antiDiag2(subOp(state.get_antiDiag2(), setOp(min)));
			state.set_antiDiag3(subOp(state.get_antiDiag3(), setOp(min)));
			state.set_score_offset(state.get_score_offset() + min);
		}

		// Update best
		if (state.get_curr_score() > state.get_best_score())
			state.set_best_score(state.get_curr_score());

		// antiDiag swap, offset updates, and new base load
		short nextDir = dir ^ 1;

		if (nextDir == goRIGHT)
			state.moveRight();
		else
			state.moveDown();

		// Update direction
		dir = nextDir;
	}
}

//======================================================================================
// X-DROP ADAPTIVE BANDED ALIGNMENT
//======================================================================================

void
XavierOneDirection (XavierState& state) {
	printf("queryh: %s\n", state.queryh);
	printf("queryv: %s\n", state.queryv);

	// PHASE 1 (initial values load using dynamic programming)
	XavierPhase1(state);
	if(state.xDropCond)  return;

	// PHASE 2 (core vectorized computation)
	XavierPhase2(state);
	if(state.xDropCond) return;

	// PHASE 3 (align on one edge)
	// GG: Phase3 removed to read to code easily (can be recovered from simd/ folder or older commits)

	// PHASE 4 (reaching end of sequences)
	XavierPhase4(state);
	if(state.xDropCond) return;
}

std::pair<int, int>
XavierXDrop
(
	SeedX& seed,
	ExtDirectionL direction,
	std::string const& target,
	std::string const& query,
	ScoringSchemeX& scoringScheme,
	int const &scoreDropOff
)
{
	// TODO: add check scoring scheme correctness/input parameters

	if (direction == XAVIER_EXTEND_LEFT)
	{
		SeedX _seed = seed; // need temporary datastruct

		std::string targetPrefix = target.substr (0, getEndPositionH(seed));	// from read start til start seed (seed included)
		std::string queryPrefix  = query.substr  (0, getEndPositionV(seed));	// from read start til start seed (seed included)
		std::reverse (targetPrefix.begin(), targetPrefix.end());
		std::reverse (queryPrefix.begin(),  queryPrefix.end());

		XavierState result (_seed, targetPrefix, queryPrefix, scoringScheme, scoreDropOff);

		if (targetPrefix.length() >= VECTORWIDTH || queryPrefix.length() >= VECTORWIDTH)
			XavierOneDirection (result);

		setBeginPositionH(seed, getEndPositionH(seed) - getEndPositionH(result.seed));
		setBeginPositionV(seed, getEndPositionV(seed) - getEndPositionV(result.seed));

		return std::make_pair(result.get_best_score(), result.get_curr_score());
	}
	else if (direction == XAVIER_EXTEND_RIGHT)
	{
		SeedX _seed = seed; // need temporary datastruct

		std::string targetSuffix = target.substr (getBeginPositionH(seed), target.length()); 	// from end seed until the end (seed included)
		std::string querySuffix  = query.substr  (getBeginPositionV(seed), query.length());		// from end seed until the end (seed included)

		XavierState result (_seed, targetSuffix, querySuffix, scoringScheme, scoreDropOff);

		if (targetSuffix.length() >= VECTORWIDTH || querySuffix.length() >= VECTORWIDTH)
			XavierOneDirection (result);

		setEndPositionH (seed, getBeginPositionH(seed) + getEndPositionH(result.seed));
		setEndPositionV (seed, getBeginPositionV(seed) + getEndPositionV(result.seed));

		return std::make_pair(result.get_best_score(), result.get_curr_score());
	}
	else
	{
		SeedX _seed1 = seed; // need temporary datastruct
		SeedX _seed2 = seed; // need temporary datastruct

		std::string targetPrefix = target.substr (0, getEndPositionH(seed));	// from read start til start seed (seed not included)
		std::string queryPrefix  = query.substr  (0, getEndPositionV(seed));	// from read start til start seed (seed not included)

		std::reverse (targetPrefix.begin(), targetPrefix.end());
		std::reverse (queryPrefix.begin(),  queryPrefix.end());

		std::cout << targetPrefix << std::endl << std::endl;
		std::cout << queryPrefix << std::endl;

		XavierState result1(_seed1, targetPrefix, queryPrefix, scoringScheme, scoreDropOff);

		if (targetPrefix.length() < VECTORWIDTH || queryPrefix.length() < VECTORWIDTH)
		{
			setBeginPositionH (seed, getEndPositionH(seed) - targetPrefix.length());
			setBeginPositionV (seed, getEndPositionV(seed) - queryPrefix.length());
		}
		else
		{
			XavierOneDirection (result1);

			setBeginPositionH (seed, getEndPositionH(seed) - getEndPositionH(result1.seed));
			setBeginPositionV (seed, getEndPositionV(seed) - getEndPositionV(result1.seed));
		}
		std::cout << "here" << std::endl;

		std::string targetSuffix = target.substr (getEndPositionH(seed), target.length()); 	// from end seed until the end (seed included)
		std::string querySuffix  = query.substr  (getEndPositionV(seed), query.length());	// from end seed until the end (seed included)

		std::cout << targetSuffix << std::endl << std::endl;
		std::cout << querySuffix << std::endl;

		XavierState result2(_seed2, targetSuffix, querySuffix, scoringScheme, scoreDropOff);

		if (targetSuffix.length() < VECTORWIDTH || querySuffix.length() < VECTORWIDTH)
		{
			std::cout << "here2" << std::endl;
			setBeginPositionH (seed, getEndPositionH(seed) + targetSuffix.length());
			setBeginPositionV (seed, getEndPositionV(seed) + querySuffix.length());
		}
		else
		{
			XavierOneDirection (result2);

			setEndPositionH (seed, getEndPositionH(seed) + getEndPositionH(result2.seed));
			setEndPositionV (seed, getEndPositionV(seed) + getEndPositionV(result2.seed));
		}

		// seed already updated and saved in result1
		// this operation sums up best and exit scores for result1 and result2 and stores them in result1
		result1 += result2;
		return std::make_pair(result1.get_best_score(), result1.get_curr_score());
	}
}
#endif
