#include "pch.h"
#include "BSND.h"

BSND::BSND()
{
}

BSND::BSND(RandomGenerator* _gen, std::vector<double> _s, std::vector<double> _r, Eigen::MatrixXd _VCV,int dim)
	: RandomProcess(_gen, dim), s(_s), r(_r), VCV(_VCV)
{
}


