#include "pch.h"
#include "EUBasketCall.h"
#include <numeric>
#include <iterator>
#include "BS_ClosedForm.h"

EUBasketCall::EUBasketCall()
{
}

EUBasketCall::EUBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights)
	: Option(_process, _K, _r, _T), weights(_weights)
{

}
EUBasketCall::EUBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _S, Eigen::MatrixXd _Vol)
	: Option(_process, _K, _r, _T), weights(_weights), S(_S), Vol(_Vol)
{
}

double EUBasketCall::ComputePrice(int NbSim,bool antithetic)
{
	double somme = 0.;
	double last_value;
	double WS_T;

	for (int n = 0; n < NbSim; ++n)
	{
		
		if (antithetic)
		{
			process->Simulate_Antithetic(0, T, T * 365);
			std::vector<double> S_T = process->Get_ValueND(T, antithetic);
			//print_vector(S_T);
			int dim = S_T.size() / 2;
			std::vector<double> S_T_anti(dim);
			for (int i = 0;i < dim;i++)
			{
				S_T_anti[i] = (S_T[i] + S_T[i + dim]) / 2;
			}
			//print_vector(S_T_anti);
			WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T_anti), 0.0);
			//std::cout << "WS_T: " << WS_T << std::endl;
		}
		else
		{
			process->Simulate(0, T, T * 365);
			std::vector<double> S_T = process->Get_ValueND(T, antithetic);
			//print_vector(S_T);
			WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T), 0.0);
			//std::cout << "WS_T normal: " << WS_T << std::endl;
		}
		
		last_value = std::max(WS_T - K, 0.);
		somme = somme + last_value;
	}

	double price = std::exp(-r[0] * T) * (somme / NbSim); //mettre r en double
	return price;

}

double EUBasketCall::ComputePrice_ControlVariate(int NbSim)
{
	double somme = 0.;
	double last_value;
	double WS_T;
	double WS_T_L;
	double X;
	double Y;
	double X_prime;

	double E_Y = Compute_E_Y(S, weights, K, r[0], Vol, T);
	
	for (int n = 0; n < NbSim; ++n)
	{
		process->Simulate(0, T, T * 365);
		std::vector<double> S_T = process->Get_ValueND(T);
		std::vector<double> S_T_L(S.size());
		for (int j = 0; j < S.size();j++)
		{
			S_T_L[j] = log(S_T[j]);
		}
		WS_T = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T), 0.0);
		WS_T_L = std::inner_product(std::begin(weights), std::end(weights), std::begin(S_T_L), 0.0);

		X = std::max(WS_T - K, 0.);
		Y= std::max(exp(WS_T_L) - K, 0.);
		X_prime = X - Y + E_Y;

		somme = somme + X_prime;
	}

	double price = std::exp(-r[0] * T) * (somme / NbSim);
	return price;
}



void print_vector(const std::vector<double>& v)
{
	std::ostream_iterator<double> out_it(std::cout, ", ");
	std::copy(v.begin(), v.end(), out_it);
	std::cout << std::endl;
}

double Compute_E_Y(std::vector<double> S, std::vector<double> weights, double K, double r, Eigen::MatrixXd Vol, double T)
{
	
	double S_Y=1;
	double R_Y;
	double vol_Y;

	Eigen::VectorXd weights_M = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());
	Eigen::MatrixXd B;

	Eigen::EigenSolver<Eigen::MatrixXd> es(Vol);

	if (Vol.determinant() == 0) //if vol not definite positive
	{
		Eigen::MatrixXd D = es.pseudoEigenvalueMatrix();
		Eigen::MatrixXd P = es.pseudoEigenvectors();
		B = P * D;
	}
	else //vol is definite positive
	{
		B = Vol.llt().matrixL();
	}

	for (int i = 0;i < S.size(); i++)
	{
		S_Y *= pow(S[i], weights[i]);	
	}

	Eigen::MatrixXd sigma2 = B * B;
	Eigen::VectorXd sigmai = sigma2.colwise().sum();
	double R1 = 0.5* weights_M.transpose() * sigmai;
	double R2 = 0.5 * weights_M.transpose() * B * B.transpose() * weights_M;
	R_Y = r - R1 + R2;


	vol_Y= pow(weights_M.transpose() * B * B.transpose() * weights_M,0.5);

	double price = bs_price(S_Y, K, vol_Y, T, R_Y, true);

	return price;
	
}

