#pragma once
#include "Option.h"
#include <eigen-3.4.0/Eigen/Dense>

class BermudanBasketCall : public Option
{
public:
	BermudanBasketCall();
	BermudanBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _exeDates,  int _L = 5);
	BermudanBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _exeDates, std::vector<double> _S, Eigen::MatrixXd _VCV, int _L = 5);

	double ComputePrice(int NbSim, bool antithetic = false);
	double ComputePrice_ControlVariate(int NbSim);
	double ComputePrice_VDC(int NbSim);


private:
	std::vector<double> exeDates;
	int L;
	std::vector<double> weights;
	std::vector<double> S;
	Eigen::MatrixXd VCV;
};

double Compute_E_Ybis(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd VCV, double T);

