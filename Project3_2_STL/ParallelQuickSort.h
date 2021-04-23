#pragma once
#include <omp.h>
#include "math.h"
#include <algorithm>
#include <list>
#include <iostream>
#include <omp.h>
#include <cstdlib>

class ParallelQuickSort
{
public: 
	
	double* pData;
	int first, last,Size;
	

public:
	
	void serial_quick_sort(double* pData, int first, int last);
	void parallel_quick_sort(double* pData, int Size);
	void set_block_pairs(int* BlockPairs, int Iter);
	void compare_split_blocks(double* pFirstBlock, int& FirstBlockSize, double* pSecondBlock, int& SecondBlockSize, double Pivot);
	int find_my_pair(int* BlockPairs, int ThreadID, int Iter);
	
};