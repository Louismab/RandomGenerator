#pragma once
#include "BSND.h"

class BSEULERND : public BSND
{
public:
	BSEULERND();
	BSEULERND(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _vol, int dim);
	void Simulate(double start_time, double end_time, size_t nb_steps);
	//int get_dim();

private:
	Eigen::MatrixXd B;
};

std::vector<double> operator*(std::vector<double> lhs, double rhs);

