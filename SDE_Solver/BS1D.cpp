#include "pch.h"
#include "BS1D.h"


BS1D::BS1D()
{
}

BS1D::BS1D(RandomGenerator* _gen, double _s, double _r, double _vol)
	: RandomProcess(_gen, 1), s(_s), r(_r), vol(_vol)
{
}

//const double BS1D::Get_Value(double time)
//{
//	return paths[0]->GetValue(time);
//}