#include "pch.h"
#include "BSEULER1D.h"
#include "SinglePath.h"
#include "iostream"


BSEULER1D::BSEULER1D()
{
}

BSEULER1D::BSEULER1D(RandomGenerator* _gen, double _s, double _r, double _vol)
	: BS1D(_gen, _s, _r, _vol)
{
}

void BSEULER1D::Simulate(double start_time, double end_time, size_t nb_steps)
{
	SinglePath* path=new SinglePath(start_time, end_time, nb_steps);
	path->AddValue(s);
	double last = s;
	double dt = (end_time- start_time)/nb_steps;
	
	for (size_t i=0; i < nb_steps; i++)
	{
		double dW = pow(dt,0.5) * gen->Generate();
		double next = last + last * (r * dt + vol * dW);
		//std::cout << "simulation " << i << ": " << next << std::endl;
		path->AddValue(next);
		last = next; 
		//std::cout << i << std::endl;
	}

	if (paths[0])
	{
		delete paths[0];
	}
	paths[0]=path;

}

void BSEULER1D::Simulate_Antithetic(double start_time, double end_time, size_t nb_steps)
{
	SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
	SinglePath* path_anti = new SinglePath(start_time, end_time, nb_steps);

	path->AddValue(s);
	path_anti->AddValue(s);

	double last = s;
	double last_anti = s;
	double dt = (end_time - start_time) / nb_steps;

	for (size_t i = 0; i < nb_steps; i++)
	{
		double dW = pow(dt, 0.5) * gen->Generate();

		double next = last + last * (r * dt + vol * dW);
		double next_anti = last_anti + last_anti * (r * dt - vol * dW);

		path->AddValue(next);
		path_anti->AddValue(next_anti);

		last = next;
		last_anti = next_anti;
	}

	if (paths[0])
	{
		delete paths[0];
	}
	paths[0] = path;

	if (paths_antithetic[0])
	{
		delete paths_antithetic[0];
	}
	paths_antithetic[0] = path_anti;

}