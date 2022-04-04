#pragma once
#include "RandomGenerator.h"
#include "SinglePath.h"

class RandomProcess
{
	public:
		RandomProcess();
		RandomProcess(RandomGenerator* _gen, int _dim);
		virtual void Simulate(double start_time, double end_time, size_t nb_steps) = 0;
		virtual void Simulate_Antithetic(double start_time, double end_time, size_t nb_steps) = 0;
		void add_path(SinglePath* Path);
		SinglePath* GetPath(int dimension = 0);
		const double Get_Value(double time,int dim=0);
		const std::vector<double> Get_ValueND(double time, bool antithetic=false);
	
	protected:
		std::vector<SinglePath*> paths;
		RandomGenerator* gen;
		int dim;

};

