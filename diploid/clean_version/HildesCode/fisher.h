// Header file: definitions of global variables, function prototypes

#ifndef FISHER_H
#define FISHER_H

#include <iostream>
#include <dynamic_bitset.hpp>	 /* for dynamic_bitset objects (www.boost.org)*/
#include "MersenneTwister.h"
using namespace std;

// Global variables:
//#define fichierLecture "parametres.txt"     // names of input
// "chr": represents a chromosome:

struct chr
{
	boost::dynamic_bitset<> sel; // selected loci (chain of 0 and 1)
};

// Function prototypes:

void ouvrirFichierE(char * param);
bool lireFichier(int & Dr, int & Nr, double & mr, 
				int & lr, double & Lr,
				double & kr, double & ar, double & fr, double & pr, double & rr,
				int & tr, int & pasr);
void recursion(int Dv, int Nv, double mv, 
				int lv, double Lv,
				double kv, double av, double fv, double pv, double rv,
				int tv, int pasv);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(chr &res, chr &c1, chr &c2, double R, int nS);
boost::dynamic_bitset<> RandomMask(int N);
void freerec(chr &res, chr &c1, chr &c2, int nS);
#endif
