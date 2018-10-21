// Functions to open input and output files,
// read parameter values from input file and 
// write them in output file.

#include "fisher.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;

// opens input file:

void ouvrirFichierE(char * param)    
{						 
	fichierE = fopen(param,"r");
}


// opens output file:

void ouvrirFichierS()   
{
	fichierS = fopen(fichierEcriture,"a");
}


// reads parameter values from input file,
// returns 1 if end of input file, else returns 0

bool lireFichier(int &dr, int &Nr, double &migr, int &br, int &nr, int &mr, double &sigr, double &ar, double &diffr, double &Qr, double &Ur, int &nbSr, double &Lr,
                 int &T1r, int &T2r, int &T3r, int &pasr)
{					 
	int x;
	bool term;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
		// each parameter set must start with *
	if (x == EOF)
		term = true;
	else
	{
		fscanf(fichierE,"%d ",&dr);
        fscanf(fichierE,"%d ",&Nr);
		fscanf(fichierE,"%lf ",&migr);
        fscanf(fichierE,"%d ",&br);
		fscanf(fichierE,"%d ",&nr);
		fscanf(fichierE,"%d ",&mr);
		fscanf(fichierE,"%lf ",&sigr);
		fscanf(fichierE,"%lf ",&ar);
		fscanf(fichierE,"%lf ",&diffr);
		fscanf(fichierE,"%lf ",&Qr);
		fscanf(fichierE,"%lf ",&Ur);
		fscanf(fichierE,"%d ",&nbSr);
		fscanf(fichierE,"%lf ",&Lr);
		fscanf(fichierE,"%d ",&T1r);
		fscanf(fichierE,"%d ",&T2r);
		fscanf(fichierE,"%d ",&T3r);
		fscanf(fichierE,"%d ",&pasr);
		
		term = false;
	} 
	return term;
}


// writes parameter values in output file:

void ecrireParametres(int dv, int Nv, double migv, int bv, int nv, int mv, double sigv, double av, double diffv, double Qv, double Uv, int nbSv, double Lv,
                      int T1v, int T2v, int T3v, int pasv)
{
	fprintf(fichierS,"\n_________________________________________\n");
	fprintf(fichierS,"\nd = %d", dv);
    fprintf(fichierS,", N = %d", Nv);
	fprintf(fichierS,", mig = %g", migv);
    fprintf(fichierS,", b = %d", bv);
	fprintf(fichierS,", n = %d", nv);
	fprintf(fichierS,", m = %d", mv);
	fprintf(fichierS,", sigma = %g", sigv);
	fprintf(fichierS,", alpha = %g", av);
	fprintf(fichierS,", diff = %g", diffv);
	fprintf(fichierS,"\nQ = %g", Qv);
	fprintf(fichierS,", U = %g", Uv);
	fprintf(fichierS,", nbS = %d", nbSv);
	fprintf(fichierS,", L = %g", Lv);
	fprintf(fichierS,", T1 = %d", T1v);
	fprintf(fichierS,", T2 = %d", T2v);
	fprintf(fichierS,", T3 = %d", T3v);
	fprintf(fichierS,", pas = %d", pasv);
}
