#pragma once
#include "UniformGenerator.h"

class QuasiGenerator : public UniformGenerator
{
	public:
		QuasiGenerator();

	protected:
		myLong current_n;
};

