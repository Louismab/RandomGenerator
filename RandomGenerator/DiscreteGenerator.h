#pragma once
//#include "RandomGenerator.h" no need anymore because uniformGenerator include RandomGenerator
#include "UniformGenerator.h"

class DiscreteGenerator : public RandomGenerator
{
public:
	DiscreteGenerator();
	DiscreteGenerator(UniformGenerator* _gen); //si on passe pas un pointeur mais direct la classe, on pourra pas passer un linearcongruential ou ecuyer par exemple mais juste un objet uniformgenerator

protected:
	UniformGenerator* generatorPtr;
};

