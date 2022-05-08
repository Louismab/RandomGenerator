#include "pch.h"
#include "ExponentialID.h"
#include <math.h>

ExponentialID::ExponentialID()
{
}

ExponentialID::ExponentialID(double _lambda, UniformGenerator* _gen)
	:Exponential(_lambda,_gen)
{
}

double ExponentialID::Generate()
{
	double u = generatorPtr->Generate();

	return -log(u) / lambda;
}