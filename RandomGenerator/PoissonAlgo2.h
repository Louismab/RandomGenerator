#pragma once
#include "Poisson.h"
#include "Exponential.h"

class PoissonAlgo2 : public Poisson
{

public:
	PoissonAlgo2();
	PoissonAlgo2(Exponential* _gen);
	virtual double Generate();

private:
	Exponential* generatorPtr;


};

