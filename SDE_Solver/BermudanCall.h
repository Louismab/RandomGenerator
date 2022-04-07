#pragma once
#include "Option.h"
#include <eigen-3.4.0/Eigen/Dense>

class BermudanCall : public Option
{

public:
	BermudanCall();
	BermudanCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _c);

	double ComputePrice(int NbSim, bool antithetic = false);
	double ComputePrice_ControlVariate(int NbSim);

private:
	std::vector<double> c;

};

