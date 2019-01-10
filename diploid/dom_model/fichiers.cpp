// Functions to open input and output files,
// read parameter values from input file and 
// write them in output file.

#include "fisher.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierE;
extern FILE * fichierM;

// opens input file:
void ouvrirFichierE(char * param)    
{						 
	fichierE = fopen(param,"r");
}

//open input file with mutational effects and dominance values
void ouvrirFichierM(char * muteff)
{
        fichierM = fopen(muteff, "r");
}


// reads parameter values from input file,
// returns 1 if end of input file, else returns 0
bool lireFichier(int & Nr, 
				double & Lr, double & Ur,
				double & kr, double & ar, int & nr, double & Fr, double & qr,
				int & Rr, int & tr, int & repsr)
{					 
	int x;
	bool term;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
		// each parameter set must start with *
	if (x == EOF)
		term = true;
	else
	{
        fscanf(fichierE,"%d ",&Nr);
		fscanf(fichierE,"%lf ",&Lr);
		fscanf(fichierE,"%lf ",&Ur);	
		fscanf(fichierE,"%lf ",&kr);
		fscanf(fichierE,"%lf ",&ar);
		fscanf(fichierE,"%d ",&nr);
		fscanf(fichierE,"%lf ",&Fr);
		fscanf(fichierE,"%lf ",&qr);		

		fscanf(fichierE,"%d ",&Rr);
		fscanf(fichierE,"%d ",&tr);
		fscanf(fichierE,"%d ",&repsr);
		
		term = false;
	} 
	return term;
}
