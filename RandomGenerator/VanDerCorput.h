#pragma once
#include "QuasiGenerator.h"

class VanDerCorput : public QuasiGenerator
{
public:
	VanDerCorput();
	virtual double Generate();
	

};

std::vector<int> IntToInverseBinary(int n);