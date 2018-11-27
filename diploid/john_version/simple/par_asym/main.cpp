// main() function: reads parameter values in input file,
// and runs the simulation.

#include "fisher.h"
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// input and output files:
FILE * fichierE;

// random number generator (Mersenne Twister):
MTRand rnd;

int main(int argc, char * argv[])
{
	// definitions of variables:
		
	int D, N, l, R, t, reps;
	double m, L, k, a, f, fp1, fp2, r;
	
	// opens input and output files:

	bool fin;
	ouvrirFichierE(argv[1]);
	fin = false;

	do
	{
		// reads parameter values;
		fin = lireFichier( D,  N,  m, l,  L, k,  a,  f,  fp1, fp2,  r, R, t, reps);
		if (!fin)
		{
			// runs the simulation:
			for(int i = 0; i<reps; i++)
				recursion( D,  N,  m, l,  L, k,  a,  f,  fp1, fp2,  r, R,  t, i+1);
		}
	} while (!fin);
	// closes files:
	fclose(fichierE);
	return 0 ;
}
