// main() function: reads parameter values in input file,
// and runs the simulation.

#include "fisher.h"
#include <iostream>
#include "MersenneTwister.h"
using namespace std;

// input and output files:

FILE * fichierE;
FILE * fichierM;

// random number generator (Mersenne Twister):

MTRand rnd;

int main(int argc, char * argv[])
{
	// definitions of variables:
		
	int N, t, n, R, reps;
	double L, U, k, a, F, q;
	
	// opens input and output files:

	bool fin;
	bool newmut = true;
	
	//read mutational effects from input file if this is provided
        if(argc > 2)
        {
                cout << "mutational effects read from file\n";
                newmut = false;
                ouvrirFichierM(argv[2]);
        }
        else
                cout << "mutational effects freshly generated\n";



	ouvrirFichierE(argv[1]);
	

	fin = false;

	do
	{
		// reads parameter values;
		fin = lireFichier( N, L, U, k,  a,  n, F, q, R, t, reps);
		if (!fin)
		{
			// runs the simulation:
			for(int i=0; i<reps; i++)
				recursion( N, L, U, k,  a,  n, F, q, R, t, i+1, newmut);
		}
	} while (!fin);

	// closes files:
	fclose(fichierE);
	//fclose(fichierM);
	return 0 ;
}
