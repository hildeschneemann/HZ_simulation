// main() function: reads parameter values in input file,
// and runs the simulation.

#include "fisher.h"
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// input and output files:

FILE * fichierE;
FILE * fichierS;
FILE * fichierM;

// random number generator (Mersenne Twister):

MTRand rnd;

int main(int argc, char * argv[])
{
	// definitions of variables:

	int d, Nt, b, n, m, nbS, T1, T2, T3, pas;
	double mig, sig, a, diff, Q, U, L;

	// opens input and output files:

	bool fin;
	bool newmut = true;
	ouvrirFichierE(argv[1]);
	

	//read mutational effects from input file if this is provided
	if(argc > 2)
	{
		cout << "mutational effects read from file\n";
		newmut = false;
		ouvrierFichierM(argv[2]);
	}
	else
		cout << "mutational effects freshly generated\n";
	
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

			recursion(d, Nt, mig, b, n, m, sig, a, diff, Q, U, nbS, L, T1, T2, T3, pas, no, newmut);

			no++;
		}
	} while (!fin);

	// closes files:

	fclose(fichierE);
	fclose(fichierS);
	fclose(fichierM);

	return 0 ;
}
