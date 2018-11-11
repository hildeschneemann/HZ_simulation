#include <iostream>
#include <cmath>
#include "MersenneTwister.h"
using namespace std;

extern MTRand rnd;

void brownian_bridge(double * Brown, int lv, int nv)
{
	int i =0;
	double gasdev();
	double ii =1.0;
	double * times = new double [lv];
	double * distW = new double [lv];
	double * cumW = new double [lv];
	//double * Brown = new double [nbSv+1];

        times[0] = distW[0] = cumW[0] = Brown[0] = 0.0;

        for (i = 1; i < lv; i++)
        {
                times[i] = ii / (lv-1);
                distW[i] = gasdev() / sqrt(lv * nv);
                cumW[i] = cumW[i-1] + distW[i];
                ii += 1.0;
        }

        for (i = 1; i < lv; i++)
        {
                Brown[i] = cumW[i] - times[i] * cumW[lv-1];
	//	cout << i << "\t" << Brown[i] << "\n";
	//cout << times[i] << "\t" << distW[i] << "\t" << cumW[i] << "\t" << Brown[i] << "\n";
         }
//	Brown[nbSv] = 0.0;
}
