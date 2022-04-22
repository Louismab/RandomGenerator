#pragma once
#include "Option.h"


class EUCall : public Option
{
	public:
		EUCall();
		EUCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T);

		double ComputePrice(int NbSim, bool antithetic = false);
		double ComputePrice_ControlVariate(int NbSim);


};

