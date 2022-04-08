#include "pch.h"
#include "BermudanBasketCall.h"


BermudanBasketCall::BermudanBasketCall()
{
}

BermudanBasketCall::BermudanBasketCall(RandomProcess* _process, double _K, std::vector<double> _r, double _T, std::vector<double> _weights, std::vector<double> _exeDates, int _L)
	: Option(_process, _K, _r, _T), exeDates(_exeDates), L(_L), weights(_weights)
{
}

double BermudanBasketCall::ComputePrice(int NbSim, bool antithetic)
{

	int nb_exe_dates = exeDates.size();
	Eigen::MatrixXd S_exe_M(NbSim, nb_exe_dates);
	Eigen::MatrixXd Tau(NbSim, nb_exe_dates);
	Eigen::VectorXd weights_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

	double tk;
	double somme = 0.;
	double last_value;

	for (int n = 0; n < NbSim; ++n)
	{
		Tau(n, nb_exe_dates - 1) = T;
		process->Simulate(0, T, T * 365);

		for (int k = 0; k < nb_exe_dates;k++)
		{
			std::vector<double> Snk = process->Get_ValueND(exeDates[k]);
			Eigen::VectorXd Snk_V = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Snk.data(), Snk.size());
			S_exe_M(n, k) = weights_V.transpose()* Snk_V;
		}

	}

	for (int k = nb_exe_dates - 2; k >= 0;k--)
	{
		tk = exeDates[k];
		Eigen::MatrixXd A(NbSim, L);
		Eigen::MatrixXd B(NbSim, 1);

		for (int n = 0; n < NbSim; n++)
		{
			double Tk1 = Tau(n, k + 1);
			int j = nb_exe_dates - 1;
			while (Tk1 != exeDates[j])
			{
				j--;
			}

			B(n, 0) = std::exp(-r[0] * (Tau(n, k + 1) - tk)) * std::max(S_exe_M(n, j) - K, 0.);


			for (int z = 0;z < L;z++)
			{
				A(n, z) = fp(S_exe_M(n, k), z);
			}

		}

		Eigen::VectorXd a(A.colPivHouseholderQr().solve(B));

		for (int n = 0; n < NbSim; n++)
		{
			double payoff_actual = std::max(S_exe_M(n, k) - K, 0.);
			double expected_payoff = A.row(n) * a;
			//std::cout << expected_payoff << std::endl;

			if (payoff_actual >= expected_payoff)
			{
				Tau(n, k) = tk;
			}
			else
			{
				Tau(n, k) = Tau(n, k + 1);
			}
		}

	}

	for (int n = 0; n < NbSim; n++)
	{
		double T0 = Tau(n, 0);
		int j = nb_exe_dates - 1;
		while (T0 != exeDates[j])
		{
			j--;
		}
		somme += std::exp(-r[0] * T0) * std::max(S_exe_M(n, j) - K, 0.);

	}

	double price = (somme / NbSim);
	return price;


}

double BermudanBasketCall::ComputePrice_ControlVariate(int NbSim)
{
	return 0.;
}