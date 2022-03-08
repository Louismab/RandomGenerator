#include "pch.h"
#include "Option.h"


Option::Option()
{
}

Option::Option(RandomProcess* _process, double _s, double _K, double _r, double _vol, double _T)
	: process(_process), s(_s), K(_K), r(_r), vol(_vol), T(_T)
{

}