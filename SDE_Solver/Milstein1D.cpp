#include "pch.h"
#include "Milstein1D.h"
#include "SinglePath.h"
#include <cmath>

Milstein1D::Milstein1D()
{
}

Milstein1D::Milstein1D(RandomGenerator* _gen, double _s, double _r, double _vol)
	: BS1D(_gen, _s, _r, _vol)
{
}

void Milstein1D::Simulate(double start_time, double end_time, size_t nb_steps)
{
	SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
	path->AddValue(s);
	double last = s;
	double dt = (end_time - start_time) / nb_steps;

	for (size_t i{}; i < nb_steps; i++)
	{
		double dW = pow(dt, 0.5) * gen->Generate();
		double next = last + last * ( (r - 0.5 * pow(vol, 2.)) * dt + vol * dW + 0.5 * pow(vol, 2) * pow(dW, 2) );
		path ->AddValue(next);
		last = next;
	}

	if (paths[0])
		delete paths[0];
	paths[0] = path;

}

void Milstein1D::Simulate_Antithetic(double start_time, double end_time, size_t nb_steps)
{
	//create two singlepaths : a normal path and the antithetic path

	SinglePath* path = new SinglePath(start_time, end_time, nb_steps);
	SinglePath* path_anti = new SinglePath(start_time, end_time, nb_steps);

	path->AddValue(s);
	path_anti->AddValue(s);

	double last = s;
	double last_anti = s;
	double dt = (end_time - start_time) / nb_steps;

	for (size_t i{}; i < nb_steps; i++)
	{
		double dW = pow(dt, 0.5) * gen->Generate();
		double next = last + last * ((r - 0.5 * pow(vol, 2.)) * dt + vol * dW + 0.5 * pow(vol, 2) * pow(dW, 2));
		double next_anti = last + last * ((r - 0.5 * pow(vol, 2.)) * dt - vol * dW + 0.5 * pow(vol, 2) * pow(dW, 2));

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