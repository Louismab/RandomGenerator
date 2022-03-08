#include "pch.h"
#include "RandomProcess.h"
#include <vector>
#include "SinglePath.h"


RandomProcess::RandomProcess()
{
}

RandomProcess::RandomProcess(RandomGenerator* _gen, int _dim)
	:gen(_gen), dim(_dim),paths(_dim,new SinglePath())
{
	//std::vector<SinglePath*> paths;
}

void RandomProcess::add_path(SinglePath* Path)
{
	paths.push_back(Path);
}

SinglePath* RandomProcess::GetPath(int dimension)
{
	return paths[0];
}

const double RandomProcess::Get_Value(double time,int dim)
{
	return paths[dim]->GetValue(time);
}


