// Header file: definitions of global variables, function prototypes

#ifndef FISHER_H
#define FISHER_H

#include <iostream>
#include <dynamic_bitset.hpp>	 /* for dynamic_bitset objects (www.boost.org)*/
#include "MersenneTwister.h"
#include "H5Cpp.h"
using namespace std;


// Global variables:

//#define fichierLecture "parametres.txt"     // names of input
#define fichierEcriture "resultats.txt"		// and output files

// "chr": represents a chromosome:

struct chr
{
	boost::dynamic_bitset<> sel; // selected loci (chain of 0 and 1)
};

// Function prototypes:

void ouvrirFichierE(char * param);
void ouvrirFichierS();
void ecrireParametres(int dv, int Nv, double migv, int bv, int nv, int mv, double sigv, double av, double diffv, double Qv, double Uv, int nbSv, double Lv,
                      int T1v, int T2v, int T3v, int pasv);
bool lireFichier(int &dr, int &Nr, double &migr, int &br, int &nr, int &mr, double &sigr, double &ar, double &diffr, double &Qr, double &Ur, int &nbSr, double &Lr,
                 int &T1r, int &T2r, int &T3r, int &pasr);
void recursion(int dv, int Nv, double migv, int bv, int nv, int mv, double sigv, double av, double diffv, double Qv, double Uv, int nbSv, double Lv,
               int T1v, int T2v, int T3v, int pasv, int nov);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(chr &res, chr &c1, chr &c2, double R, int nS);
boost::dynamic_bitset<> RandomMask(int N);
void freerec(chr &res, chr &c1, chr &c2, int nS);
void createHDF5(const H5std_string &FILE_NAME,const H5std_string &DSET_FREQ_NAME,
    const int &RANK_FREQ, const int &DIM0, const int &DIM1, const int &DIM2,
    const H5std_string &DSET_W_NAME, const int &RANK_W, const int &DIM_W0, const int &DIM_W1,
    const H5std_string &DSET_GEN_NAME);
void writeTimeStepHDF5(H5::H5File &file, H5::DataSet &dataset,
    const int &RANK, hsize_t dimsm[],
    double sdata[], int &indexGen);
void writeGenSaved(H5::H5File &file, H5::DataSet &dataset, int data[]);
void writeAttributeHDF5(H5::H5File &file, H5::DataSet &dataset,
    const char *ATTR_NAME, int attr_data);
void writeAttributeHDF5(H5::H5File &file, H5::DataSet &dataset,
        const char *ATTR_NAME, double attr_data);
void brownian_bridge(double * Brown, int nbSv, int nv);
#endif
