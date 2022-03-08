#include "pch.h"
#include "HeadTail.h"
#include "LinearCongruential.h"

HeadTail::HeadTail()
{
}

HeadTail::HeadTail(UniformGenerator* _gen)
	: DiscreteGenerator(_gen)
{
	// generatorPtr = _gen;  -> equiv to call the constructor of discreteGenerator
}

double HeadTail::Generate()
{
	double u = generatorPtr->Generate();
	if (u <= 0.5)
		return 1.;
	else
		return 0.;

	//return (generatorPtr->Generate() <= 0.5);
}