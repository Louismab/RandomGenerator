#pragma once
#include "BSND.h"

class MilsteinND : public BSND
{
public:
	MilsteinND();
	MilsteinND(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _VCV, int dim);
	void Simulate(double start_time, double end_time, size_t nb_steps);
	void Simulate_Antithetic(double start_time, double end_time, size_t nb_steps);


private:
	Eigen::MatrixXd B;
};

