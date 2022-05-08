#pragma once
#include "DiscreteGenerator.h"

class Bernoulli : public DiscreteGenerator
{
public:
	Bernoulli();
	Bernoulli(UniformGenerator* _gen, double _p);
	virtual double Generate();

private:
	double p;
};

