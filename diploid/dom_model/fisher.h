// Header file: definitions of global variables, function prototypes

#ifndef FISHER_H
#define FISHER_H

#include <iostream>
#include <dynamic_bitset.hpp>	 /* for dynamic_bitset objects (www.boost.org)*/
#include "MersenneTwister.h"
using namespace std;

// Global variables:
// "chr": represents a chromosome:

struct chr
{
	boost::dynamic_bitset<> sel; // selected loci (chain of 0 and 1)
};

// Function prototypes:

void ouvrirFichierE(char * param);
void ouvrirFichierM(char * muteff);
bool lireFichier(int & Nr, 
				double & Lr, double & Ur,
				double & kr, double & ar, int & nr, double & Fr, double & qr,
				int & Rr, int & tr, int & repsv);
void recursion(int Nv,
				double Lv, double Uv,
				double kv, double av, int nv, double Fv, double qv,
				int Rv, int tv, int repv, bool newmut);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(chr &res, chr &c1, chr &c2, double R, int nS);
boost::dynamic_bitset<> RandomMask(int N);
void freerec(chr &res, chr &c1, chr &c2, int nS);
//void backcrosses(chr * tmp, double * mutations, int lv, int nv, int Lv, double kd2, double av, int twoND, int nb4, bool withrec);
#endif
