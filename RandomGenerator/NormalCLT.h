#pragma once
#include "Normal.h"

class NormalCLT : public Normal
{
public:
	NormalCLT();
	NormalCLT(double _mean, double _var, UniformGenerator* _gen);
	virtual double Generate();
};

