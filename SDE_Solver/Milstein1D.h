#pragma once
#include "BS1D.h"

class Milstein1D : public BS1D
{
public:
	Milstein1D();
	Milstein1D(RandomGenerator* _gen, double _s, double _r, double _vol);
	void Simulate(double start_time, double end_time, size_t nb_steps);
	void Simulate_Antithetic(double start_time, double end_time, size_t nb_steps);

};

