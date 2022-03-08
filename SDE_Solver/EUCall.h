#pragma once
#include "Option.h"


class EUCall : public Option
{
	public:
		EUCall();
		EUCall(RandomProcess* _process, double _s, double _K, double _r, double _vol, double _T);

		double ComputePrice(int NbSim);

};

