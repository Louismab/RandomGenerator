#include "pch.h"
#include "DiscreteGenerator.h"


DiscreteGenerator::DiscreteGenerator()
{
}

DiscreteGenerator::DiscreteGenerator(UniformGenerator* _gen)
{
	generatorPtr = _gen;
}