#pragma once
#include "Poisson.h"

class PoissonAlgo1 : public Poisson
{

public:
	PoissonAlgo1();
	PoissonAlgo1(double _lambda, UniformGenerator* _gen);
	virtual double Generate();

};

