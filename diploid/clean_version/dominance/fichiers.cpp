// Functions to open input and output files,
// read parameter values from input file and 
// write them in output file.

#include "fisher.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierE;

// opens input file:
void ouvrirFichierE(char * param)    
{						 
	fichierE = fopen(param,"r");
}





// reads parameter values from input file,
// returns 1 if end of input file, else returns 0
bool lireFichier(int & Dr, int & Nr, double & mr, 
				int & lr, double & Lr,
				double & kr, double & ar, int & nr, double & Fr, double & qr,
				int & tr, int & pasr)
{					 
	int x;
	bool term;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
		// each parameter set must start with *
	if (x == EOF)
		term = true;
	else
	{
		fscanf(fichierE,"%d ",&Dr);
        fscanf(fichierE,"%d ",&Nr);
		fscanf(fichierE,"%lf ",&mr);
        
        fscanf(fichierE,"%d ",&lr);
		fscanf(fichierE,"%lf ",&Lr);
		
		fscanf(fichierE,"%lf ",&kr);
		fscanf(fichierE,"%lf ",&ar);
		fscanf(fichierE,"%d ",&nr);
		fscanf(fichierE,"%lf ",&Fr);
		fscanf(fichierE,"%lf ",&qr);		

		fscanf(fichierE,"%d ",&tr);
		fscanf(fichierE,"%d ",&pasr);
		
		term = false;
	} 
	return term;
}
