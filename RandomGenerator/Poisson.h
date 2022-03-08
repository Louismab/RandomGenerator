#pragma once
#include "DiscreteGenerator.h"

class Poisson : public DiscreteGenerator
{
public:
	Poisson();
	Poisson(double _lambda, UniformGenerator* _gen);

protected:
	double lambda;
};

