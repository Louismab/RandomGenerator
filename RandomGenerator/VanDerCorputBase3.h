#pragma once
#include "QuasiGenerator.h"

class VanDerCorputBase3 : public QuasiGenerator
{
public:
	VanDerCorputBase3();
	virtual double Generate();
	
};

std::vector<int> IntToTernary(int n);
