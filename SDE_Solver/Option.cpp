#include "pch.h"
#include "Option.h"


Option::Option()
{
}

Option::Option(RandomProcess* _process,  double _K, std::vector<double> _r,  double _T)
	: process(_process), K(_K), r(_r), T(_T)
{

}