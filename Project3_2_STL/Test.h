#pragma once
#include <iostream>
#include <vector>
#include "ParallelQuickSort.h"

class Test:public ParallelQuickSort
{
public:
	std::vector<unsigned long long> a{n};
	size_t n = 1024*1024*128;
public:
	
	int Testirovanie(unsigned long long n);
	
};