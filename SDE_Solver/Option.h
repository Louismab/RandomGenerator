#pragma once
#include "RandomGenerator.h"
#include "RandomProcess.h"

class Option
{
	public:
		Option();
		Option(RandomProcess* _process,  double _s, double _K, double _r, double _vol, double _T);
		virtual double ComputePrice(int NbSim) = 0;

	protected:
		double s;
		double K;
		double r;
		double vol;
		double T;
		RandomProcess* process;



};

