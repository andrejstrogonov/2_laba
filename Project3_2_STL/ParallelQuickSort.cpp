#include "ParallelQuickSort.h"
#include <algorithm>
#include <list>
#include <iostream>
#include <omp.h>
#include <cstdlib>
#include <vector>
#include <chrono>
int ThreadNum = omp_get_num_threads();
static size_t minParallelSize = 128 * 1024;
static unsigned maxParallelDepth = 4;
int DimSize, MAX_VALUE, MIN_VALUE;
void ParallelQuickSort::serial_quick_sort(double*pData, int first, int last) 
{
	if (first >= last)
	{
		return;
	}
	int PivotPos = first;
	double Pivot = pData[first];
	for (int i = first + 1; i <= last; i++) 
	{
		if (pData[i] < Pivot) 
		{
			if (i != PivotPos + 1)
			{
				std::swap(pData[i], pData[PivotPos + 1]);
				PivotPos++;
			}
		}
	}
	std::swap(pData[first], pData[PivotPos]);
	serial_quick_sort( pData, first, PivotPos - 1);
	serial_quick_sort( pData, PivotPos + 1, last);
}

void ParallelQuickSort::parallel_quick_sort(double* pData, int Size)
{
	int ThreadID = omp_get_thread_num();
	double** pTempData = new double* [2 * ThreadNum];
	int* BlockSize = new int[2 * ThreadNum];
	double* Pivots = new double[ThreadNum];
	int* BlockPairs = new int[2 * ThreadNum];
	for (int i = 0; i < 2 * ThreadNum; i++) 
	{
		pTempData[i] = new double[Size];
		BlockSize[i] = Size / (2 * ThreadNum);
	}
	for (int j = 0; j < Size; j++)
	{
		pTempData[2 * j * ThreadNum / Size]	[j % (Size / (2 * ThreadNum))] = pData[j];
	}
	#pragma omp parallel
	{
		#pragma omp parallel
		for (int i = 0; i < DimSize; i++)
		{
			#pragma omp sections
			{
				#pragma omp section
				{
					for (int j = 0; j < ThreadNum; j++)
					{
						Pivots[j] = (MAX_VALUE + MIN_VALUE) / 2;
					}
				}
				#pragma omp section
				for (int iter = 1; iter <= i; iter++)
				{
					for (int j = 0; j < ThreadNum; j++)
					{
						Pivots[j] = Pivots[j] - pow(-1.0f, j / ((2 * ThreadNum) >> (iter + 1))) * (MAX_VALUE - MIN_VALUE) / (2 << iter);
					}
				}
			}
			set_block_pairs(BlockPairs, i);
#pragma omp parallel 
			{
				int MyPair = find_my_pair(BlockPairs, ThreadID, i);
				int FirstBlock = BlockPairs[2 * MyPair];
				int SecondBlock = BlockPairs[2 * MyPair + 1];
				compare_split_blocks(pTempData[FirstBlock], BlockSize[FirstBlock], pTempData[SecondBlock], BlockSize[SecondBlock], Pivots[ThreadID]);
			}
#pragma omp parallel if (0 <maxParallelDepth && (last - first) > minParallelSize)
			if (BlockSize[2 * ThreadID] > 0)
			{
				serial_quick_sort(pData, BlockSize[2 * ThreadID], BlockSize[2 * ThreadID]);
			}
			if (BlockSize[2 * ThreadID + 1] > 0)
			{
				serial_quick_sort(pData, BlockSize[2 * ThreadID + 1], BlockSize[2 * ThreadID + 1]);
			}
		}
	}
	int curr = 0;
	for (int i = 0; i < 2 * ThreadNum; i++)
		for (int j = 0; (j < BlockSize[i]) && (curr < Size); j++)
			pData[curr++] = pTempData[i][j];
	for (int i = 0; i < ThreadNum; i++)
		delete[] pTempData[i];
	delete[] pTempData;
	delete[] BlockSize;
	delete[] Pivots;
	delete[] BlockPairs;
}
void ParallelQuickSort::set_block_pairs(int* BlockPairs, int Iter)
{
	int PairNum = 0, FirstValue, SecondValue, ThreadID;
	bool Exist;
	for (int i = 0; i < 2 * ThreadNum; i++) 
	{
		FirstValue = i;
		Exist = false;
		for (int j = 0; (j < PairNum) && (!Exist); j++)
		{
			if (BlockPairs[2 * j + 1] == FirstValue)
			{
				Exist = true;
			}
		}
		if (!Exist) 
		{
			SecondValue = FirstValue ^ (1 <<(DimSize - Iter - 1));
			BlockPairs[2 * PairNum] = FirstValue;
			BlockPairs[2 * PairNum + 1] = SecondValue;
			PairNum++;
		} 
	}
}
void ParallelQuickSort::compare_split_blocks(double* pFirstBlock, int& FirstBlockSize, double* pSecondBlock, int& SecondBlockSize, double Pivot)
{
	int TotalSize = FirstBlockSize + SecondBlockSize;
	double* pTempBlock = new double[TotalSize];
	int LastMin = 0, FirstMax = TotalSize - 1;
	for (int i = 0; i < FirstBlockSize; i++) 
	{
		if (pFirstBlock[i] < Pivot)
		{
			pTempBlock[LastMin++] = pFirstBlock[i];
		}
		else
		{
			pTempBlock[FirstMax--] = pFirstBlock[i];
		}
	}
	for (int i = 0; i < SecondBlockSize; i++) 
	{
		if (pSecondBlock[i] < Pivot)
		{
			pTempBlock[LastMin++] = pSecondBlock[i];
		}
		else
		{
			pTempBlock[FirstMax--] = pSecondBlock[i];
		}
	}
	FirstBlockSize = LastMin;
	SecondBlockSize = TotalSize - LastMin;
	for (int i = 0; i < FirstBlockSize; i++)
	{
		pFirstBlock[i] = pTempBlock[i];
	}
	for (int i = 0; i < SecondBlockSize; i++)
	{
		pSecondBlock[i] = pTempBlock[FirstBlockSize + i];
	}
	delete[] pTempBlock;
}
int ParallelQuickSort::find_my_pair(int* BlockPairs, int ThreadID,int Iter) {
	int DimSize = int(log10(double(ThreadNum))/log10(2.0));
	int BlockID = 0, index, result;
	for (int i = 0; i < ThreadNum; i++) 
	{
		BlockID = BlockPairs[2 * i];
		if (Iter == 0)
		{
			index = BlockID % (1 << DimSize - Iter - 1);
		}
		if ((Iter > 0) && (Iter < DimSize - 1))
		{
			index = ((BlockID >> (DimSize - Iter)) << (DimSize - Iter - 1)) | (BlockID % (1 << (DimSize - Iter - 1)));
		}
		if (Iter == DimSize - 1)
		{
			index = BlockID >> 1;
		}
		if (index == ThreadID) 
		{
			result = i;
			break;
		}
	}
	return result;
}