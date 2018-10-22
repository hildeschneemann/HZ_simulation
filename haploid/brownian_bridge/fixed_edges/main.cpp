// main() function: reads parameter values in input file,
// and runs the simulation.

#include "fisher.h"
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// input and output files:

FILE * fichierE;
FILE * fichierS;

// random number generator (Mersenne Twister):

MTRand rnd;

int main(int argc, char * argv[])
{
	// definitions of variables:

	int d, Nt, b, n, m, nbS, T1, T2, T3, pas;
	double mig, sig, a, diff, Q, U, L;

	// opens input and output files:

	bool fin;
	ouvrirFichierE(argv[1]);
	ouvrirFichierS();
	fin = false;

	int no = 1;
	do
	{
		// reads parameter values;

		fin = lireFichier(d, Nt, mig, b, n, m, sig, a, diff, Q, U, nbS, L, T1, T2, T3, pas);

		if (!fin)
		{
			// writes parameter values in output file:

			ecrireParametres(d, Nt, mig, b, n, m, sig, a, diff, Q, U, nbS, L, T1, T2, T3, pas);

			// runs the simulation:

			recursion(d, Nt, mig, b, n, m, sig, a, diff, Q, U, nbS, L, T1, T2, T3, pas, no);

			no++;
		}
	} while (!fin);

	// closes files:

	fclose(fichierE);
	fclose(fichierS);

	return 0 ;
}
