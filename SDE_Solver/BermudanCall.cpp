#include "pch.h"
#include "BermudanCall.h"

BermudanCall::BermudanCall()
{
}

BermudanCall::BermudanCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _c)
	: Option(_process, _K, _r, _T),c(_c)
{
}

double BermudanCall::ComputePrice(int NbSim, bool antithetic)
{

	Eigen::MatrixXd S_exe_M(NbSim, c.size());
	Eigen::MatrixXd values_M(NbSim, c.size()+1);
	

	double somme = 0.;
	double last_value;

	for (int n = 0; n < NbSim; ++n)
	{
		std::cout << "ok" << std::endl;
		process->Simulate(0, T, T * 365);
		
		for (int j = 0; j < c.size();j++)
		{
			S_exe_M(n, j) = process->Get_Value(c[j]);
		}
		
	}

	std::cout << S_exe_M << std::endl;


	//double price = std::exp(-r[0] * T) * (somme / NbSim);
	return 0;

}

double BermudanCall::ComputePrice_ControlVariate(int NbSim)
{
	return 0.;
}