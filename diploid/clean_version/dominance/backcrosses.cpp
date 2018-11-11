#include "fisher.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <bitset>
#include <string>
#include <random>
#include <boost/dynamic_bitset.hpp>
using namespace std;

#define nb_bc 1000


void backcrosses(chr * tmp, double * mutations, int lv, int nv, int Lv, double kd2, double av, int twoND,int nb4, bool withrec)
{
	chr * bcs = new chr [twoND / 2];
	double * bcs_sz2 = new double [nb_bc];
	int i, loc, k, loc1;
	double d, sz2, w, bcs_mean, bcs_var;
	vector<int> loci;

	auto rng = std::default_random_engine {};

	for (loc=0; loc < lv; loc++)
	{
		loci.push_back(loc);
	}


	for (i = 0; i < (nb_bc); i++)
	{
		bcs[i].sel.resize(lv);	
	}
	
	
	for (i = 0; i < nb_bc; i++)
	{
		std::shuffle(std::begin(loci), std::end(loci),rng);
		

		// create recombined gametes produced by F1 individuals
		
		sz2 = 0;

		for (k = 0; k < nv; k++)  // phenotypes of the individual
		{
			d = 0;
			for (loc = 0; loc < (lv / 2); loc++)
			{
				loc1 = loci[loc];
				bcs[i].sel[loc1] = 1;
				d += mutations[lv * k + loci[loc1]];
			}
			
			sz2 += d * d;
		}
		
		//print genotype
		//for (loc=0; loc < lv; loc++)
		//{	
		//	cout << bcs[i].sel[loc];	
		//}

		bcs_sz2[i] = sz2;
		bcs_mean += sz2;
		//cout << " breakdown ind " << i << " is " << sz2 << "\n";
	}
	
	bcs_mean /= nb_bc;
	
	for (i =0; i < nb_bc; i++)
	{
		bcs_var += (bcs_sz2[i] - bcs_mean) * (bcs_sz2[i] - bcs_mean);
	}
	bcs_var /= (nb_bc - 1);
	cout << bcs_mean << ", " << bcs_var << "\n";
}
