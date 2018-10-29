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
		
	int D, N, l, t, pas;
	double m, L, k, a, f, p, r;
	
	// opens input and output files:

	bool fin;
	ouvrirFichierE(argv[1]);
	fin = false;

	int no = 1;
	do
	{
		// reads parameter values;

		fin = lireFichier( D,  N,  m, l,  L, k,  a,  f,  p,  r, t,  pas);

		if (!fin)
		{
			// writes parameter values in output file:

			ecrireParametres( D,  N,  m, l,  L, k,  a,  f,  p,  r, t,  pas);

			// runs the simulation:

			recursion( D,  N,  m, l,  L, k,  a,  f,  p,  r, t,  pas, no);

			no++;
		}
	} while (!fin);

	// closes files:
	fclose(fichierE);
	return 0 ;
}
