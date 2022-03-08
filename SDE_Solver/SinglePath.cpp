#include "pch.h"
#include "SinglePath.h"
#include "iostream"

SinglePath::SinglePath()
{
}

SinglePath::SinglePath(double _start, double _end, size_t _nbSteps)
	:start(_start),end(_end), nbSteps(_nbSteps)
{
	timeStep = (end - start) / nbSteps;
	Times.push_back(_start);
}

void SinglePath::AddValue(double val)
{
	Values.push_back(val);
	
	if (Values.size() > 1)
	{
		//std::cout << "j'ajoute: " << val << "au temps: " << Times[Times.size() - 1] + timeStep << std::endl;
		Times.push_back(Times[Times.size() - 1] + timeStep);
	}
	
}

const double SinglePath::GetValue(double time)
{
	size_t i = 1;
	if (time +1< Times.size())
	{
		while (Times[i] <= time)
		{
			i++;
		}
	}
	else
	{
		i = Times.size();
	}
	
	if (Times[i - 1] == time)
	{
		return Values[i - 1];
	}
	else
	{
		//std::cout << "je rentre dans la condition 2" << std::endl;
		return ((time - Times[i - 1]) * Values[i] + (Times[i] - time) * Values[i - 1]) / (timeStep);
		//g�rer les cas o� i < start ou > end
	}

}

void SinglePath::print_vector()
{
	for (size_t i = 0; i < Times.size(); ++i)
	{
		std::cout << "t = " << Times[i] << " ; value = " << GetValue(Times[i]) << std::endl;
	}
}