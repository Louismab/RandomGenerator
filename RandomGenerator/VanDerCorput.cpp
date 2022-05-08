#include "pch.h"
#include "VanDerCorput.h"
#include <vector>

VanDerCorput::VanDerCorput()
	: QuasiGenerator()
{
}

double VanDerCorput::Generate()
{
    
    std::vector<int> b = IntToInverseBinary(current_n);
    double phi = 0;
    for (int k = 0; k < b.size(); k++)
    {
        phi += b[k] / pow(2, (k + 1));
    }
    current_n = current_n + 1;
    return phi;
    
}

std::vector<int> IntToInverseBinary(int n)
{
    std::vector<int> b;
    
    while (n > 0)
    {
        int d = n % 2;
        b.push_back(d);
        n = n / 2;
    }

    return b;
}


