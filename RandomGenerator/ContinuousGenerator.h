#pragma once
#include "RandomGenerator.h"
#include "UniformGenerator.h"

class ContinuousGenerator : public RandomGenerator
{
public:
	ContinuousGenerator();
	ContinuousGenerator(UniformGenerator* _gen); //si on passe pas un pointeur mais direct la classe, on pourra pas passer un linearcongruential ou ecuyer par exemple mais juste un objet uniformgenerator

protected:
	UniformGenerator* generatorPtr;
};

