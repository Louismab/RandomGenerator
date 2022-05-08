#pragma once
#include "RandomProcess.h"
#include <vector>
#include <iostream>
#include <eigen-3.4.0/Eigen/Dense>

class BSND : public RandomProcess
{
public:
	BSND();
	BSND(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _VCV,int dim);
	//const double Get_Value(double time);

protected:

	std::vector<double> s;
	std::vector<double> r;
	Eigen::MatrixXd VCV;
};

