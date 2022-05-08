#pragma once
#include "Option.h"
#include <eigen-3.4.0/Eigen/Dense>

class EUBasketCall :public Option
{
public:
	EUBasketCall();
	EUBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights);
	EUBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _S, Eigen::MatrixXd _VCV);

	double ComputePrice(int NbSim, bool antithetic = false);
	double ComputePrice_ControlVariate(int NbSim);
	double ComputePrice_VDC(int NbSim);
	//double ComputePrice_Antithetic(int NbSim);

private:
	std::vector<double> weights;
	std::vector<double> S;
	Eigen::MatrixXd VCV;
};

void print_vector(const std::vector<double>& v);
double Compute_E_Y(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd Vol, double T);