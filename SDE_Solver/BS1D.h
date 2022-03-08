#pragma once
#include "RandomProcess.h"

class BS1D : public RandomProcess
{
	public:
		BS1D();
		BS1D(RandomGenerator* _gen, double _s, double _r, double _vol);
		//const double Get_Value(double time);

	protected:
		
		double s;
		double r;
		double vol;


};

