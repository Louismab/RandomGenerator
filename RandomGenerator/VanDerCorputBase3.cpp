#include "pch.h"
#include "VanDerCorputBase3.h"

VanDerCorputBase3::VanDerCorputBase3()
	: QuasiGenerator(2)
{
}

double VanDerCorputBase3::Generate()
{
    std::vector<int> c = IntToTernary(current_n);
    double phi = 0;
    for (int k = 0;k < c.size();k++)
    {
        phi += c[k] / pow(3, (k + 1));
    }
    current_n = current_n + 1;
    return phi;
}


std::vector<int> IntToTernary(int n)
{
    std::vector<int> c;

    while (n > 0)
    {
        int d = n % 3;
        c.push_back(d);
        n = n / 3;
    }

    return c;
}