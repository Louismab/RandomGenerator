#pragma once
#include "QuasiGenerator.h"

class VanDerCorput : public QuasiGenerator
{
public:
	VanDerCorput(int _base=2);
	virtual double Generate();

private:
	int base;

};

std::vector<int> IntToInversePAdicExpansion(int n, int base);