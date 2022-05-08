#pragma once
#include "ContinuousGenerator.h" 

class Exponential :public ContinuousGenerator
{
public:
	Exponential();
	Exponential(double _lambda, UniformGenerator* _gen);

protected:
	double lambda;
};

