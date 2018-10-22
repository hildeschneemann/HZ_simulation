#include "fisher.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
#include "H5Cpp.h"
using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;
extern MTRand rnd;


/*----------------------------------------------------------
Function recursion: simulates a 1D stepping-stone, with a barrier
appearing after T1 generations, followed
by a secondary contact after T1 + T2 generations.

Parameters:

dv: number of demes
Nv: size of the initial population
migv: migration rate (between adjacent demes)
bv: position of the barrier
propv: number of individuals in the first deme (after the split)
nv: number of phenotypic dimensions (complexity)
mv: number of dimensions affected by a mutation (pleiotropy)
sigv: size of mutational steps along each phenotypic axis
diffv: change in optimum after split
av: strength of selection (steepness of fitness function)
Qv: curvature of fitness function (cf Gros et al 2009 Genetics)
Uv: total mutation rate per genome (on selected loci)
nbSv: total number of selected loci per genome
Lv: genome map length (mean nb of cross-overs per meiosis)
	set to -1 for free recombination
-----------------------------------------------------------*/


void recursion(int dv, int Nv, double migv, int bv, int nv, int mv, double sigv, double av, double diffv, double Qv, double Uv, int nbSv, double Lv,
			   int T1v, int T2v, int T3v, int pasv, int nov)
{
	// variables:

	int i, j, k, loc, gen, mut, chr1, chr2, ind, nb, nb1, nb2, nb3, nb4, nbMig, ns, part, nbCo;
	double w, wbar, varw, rd, pp, d, x, sz2, delta_p, p_old;
	vector<int> store;
	bool withrec;

	// various fixed quantities:

    int Nd = Nv * dv;
    int nd = nv * dv;
    int twoN = 2*Nv;
	int N_1 = Nv - 1;
	int n_1 = nv - 1;
	int nbS_1 = nbSv - 1;
	int nS = nbSv * nv;
	int Nn = Nv * nv;
	double hQ = Qv / 2.0;
    int NbGen = T1v + T2v + T3v; // total number of generations
	int NbGen_1 = NbGen - 1;
    int Tcontact = T1v + T2v; // time of secondary contact
	withrec = Lv == -1 ? false : true;
	int Nd1 = bv * Nv;
	bool last = 0;
	bool equi = false;
	bool sign = true;
	int nbSign = 0;
	withrec = Lv == -1 ? false : true;	
	int round=0;

	// HDF5 constants and variables

	string fileName;
	stringstream nameF;
	nameF << "result_d" << dv << "_N" << Nv << "_mig" << migv << "_b" << bv
		<< "_n" << nv << "_m" << mv << "_sig" << sigv << "_diff" << diffv
        << "_a" << av << "_Q" << Qv << "_U" << Uv << "_nbS" << nbSv
		<< "_L" << Lv << "_Ts" << T1v << "-" << T2v << "-" << T3v
		<< "_" << nov << ".h5";
	nameF >> fileName;

	const H5std_string FILE_NAME(fileName);
	const H5std_string DSET_FREQ_NAME("alleleFreq");
	const H5std_string DSET_W_NAME("w");
	const H5std_string DSET_GEN_NAME("savedGen");

	const int RANK_FREQ = 3; // dimensions of the dataset
	const int RANK_W = 2;
	const int DIM0_SUB = 1;	// subset dimensions
	const int DIM1_SUB = nbSv;
	const int DIM2_SUB = dv;
	const int DIM0 = NbGen / pasv + 1; // size of dataset frequencies
	const int DIM1 = nbSv;
	const int DIM2 = dv;
	const int DIM_W0 = DIM0; // size of dataset fitness
	const int DIM_W1 = dv;
	hsize_t dim_sub_freq[RANK_FREQ];
	hsize_t dim_sub_w[RANK_W];
	int indexGen = 0; // where to write in HDF5 dataset
	int accGen=0;
	double sdata_freq[DIM0_SUB*DIM1_SUB*DIM2_SUB]; // subset to write to HDF5 dataset
	double sdata_w[1*DIM_W1];
	int savedGen[DIM0];

	// population: table of N*d chromosomes:

    chr * pop = new chr [Nd];
	chr * temp = new chr [Nd]; // used to generate next generation
    chr * cp;

	// "mutations" will hold the effect of the 1 allele at each locus
	// on each phenotypic axis:

	double * mutations = new double [nS];

	// "Wtot" will hold the fitness of each individual:

	double * Wtot = new double [Nd];

	double Wmean[dv]; // mean fitness in each deme

    // maximal fitness in each deme:

    double * Wmax = new double [dv];

	// tables for means and variances of phenotypic traits:

	double * m = new double [nv];
	double * v = new double [nv];

	//brownian bridge
	
	double * Brown = new double [nbSv+1];	

	// create HDF5 file and open it

	createHDF5(FILE_NAME, DSET_FREQ_NAME, RANK_FREQ, DIM0, DIM1, DIM2,
		DSET_W_NAME, RANK_W, DIM_W0, DIM_W1, DSET_GEN_NAME);

	H5::H5File file;
	H5::DataSet dset_freq;
	H5::DataSet dset_w;
	H5::DataSet dset_gen;
	file.openFile(FILE_NAME, H5F_ACC_RDWR);
	dset_freq = file.openDataSet(DSET_FREQ_NAME);
	dset_w = file.openDataSet(DSET_W_NAME);
	dset_gen = file.openDataSet(DSET_GEN_NAME);

	dim_sub_freq[0] = DIM0_SUB;
	dim_sub_freq[1] = DIM1_SUB;
	dim_sub_freq[2] = DIM2_SUB;
	dim_sub_w[0] = 1;
	dim_sub_w[1] = DIM_W1;

	// write parameters as attributes of alleleFreq dataset
	//(int dv, int Nv, double migv, int bv, int nv, int mv, double sigv, double av, double diffv, double Qv, double Uv, int nbSv, double Lv, int T1v, int T2v, int T3v, int pasv, int nov)
	writeAttributeHDF5(file, dset_freq, "d", dv);
	writeAttributeHDF5(file, dset_freq, "N", Nv);
	writeAttributeHDF5(file, dset_freq, "mig", migv);
	writeAttributeHDF5(file, dset_freq, "b", bv);
	writeAttributeHDF5(file, dset_freq, "n", nv);
	writeAttributeHDF5(file, dset_freq, "m", mv);
	writeAttributeHDF5(file, dset_freq, "sig", sigv);
	writeAttributeHDF5(file, dset_freq, "a", av);
	writeAttributeHDF5(file, dset_freq, "diff", diffv);
	writeAttributeHDF5(file, dset_freq, "Q", Qv);
	writeAttributeHDF5(file, dset_freq, "U", Uv);
	writeAttributeHDF5(file, dset_freq, "nbS", nbSv);
	writeAttributeHDF5(file, dset_freq, "L", Lv);
	writeAttributeHDF5(file, dset_freq, "T1", T1v);
	writeAttributeHDF5(file, dset_freq, "T2", T2v);
	writeAttributeHDF5(file, dset_freq, "T3", T3v);
	writeAttributeHDF5(file, dset_freq, "pas", pasv);
	writeAttributeHDF5(file, dset_freq, "no", nov);

	// for time length measure:

	time_t debut, fin;
	struct tm *ptr;
	debut = time(0);

	// initialization: allele 0 is fixed at all selected loci:

	for (i = 0; i < Nd; i++)
    {
		pop[i].sel.resize(nbSv);
        temp[i].sel.resize(nbSv);
		if(i <= Nd1) //complete divergence
			pop[i].sel.flip();
			temp[i].sel.flip();
	}	


	// samples the effect of allele 1 at each locus on each m (randomly sampled) phenotypic axes.
	// stores the result in "mutations" table:

	for (i = 0; i < nS; i++)
		mutations[i] = 0;
	for (i = 0; i < nv; i++)
	{
		nb = nbSv * i;
		brownian_bridge(Brown, nbSv+1);
		for (j = 0; j < nbSv; j++)
		{
		//	cout << "ok before mut?\n";
			mutations[nb + j] = sigv * (Brown[j+1] - Brown[j]);
			cout << nb << "\t" << j << "\t" << mutations[nb+j] << "\n";
		}
	}
    // generations:
	while(equi == false & accGen < 10000)
	{
        cout << "last: "<< last << ", nbSign" << nbSign << ", equi: "<< equi << ", indexGen: "<< indexGen << "\n";
		for (gen = 0; gen < NbGen; gen++)
		{
			// fitness of each individual, maximal fitnesses,
			// mean fitness and variance in fitness:

			wbar = 0;
			varw = 0;
		for (k = 0; k < nv; k++)
		{
		    m[k] = 0;
		    v[k] = 0;
		}

			for (i = 0; i < dv; i++)  // for each deme
		{
			    nb = i * Nv;
		    Wmax[i] = 0;

		    for (j = 0; j < Nv; j++)  // for each individual
		    {
			nb2 = nb + j;
			sz2 = 0;

			for (k = 0; k < nv; k++)  // phenotypes of the individual
			{
			    d = 0;
			    for (loc = 0; loc < nbSv; loc++)
				if (pop[nb2].sel[loc] == 1)
				    d += mutations[nbSv * k + loc];

			    if ((gen > T1v) && (k == 0)) // change in optimum along first axis after T1 generations
				d -= diffv;

			    m[k] += d;
			    v[k] += d * d;
			    sz2 += d * d; // "sz2" is the square of the distance to the optimum
			}

			//cout << "generation" << gen << ", ok after d+=?\n";
			// fitness:

			w = exp(-av * pow(sz2, hQ));
			Wtot[nb2] = w;

			wbar += w;
			varw += w * w;
			if (Wmax[i] < w)
			    Wmax[i] = w;
		    }
		}
			wbar /= Nd;
			varw /= Nd;
		for (k = 0; k < nv; k++)
		{
		    m[k] /= Nd;
		    v[k] /= Nd;
		    v[k] -= (m[k] * m[k]);
		}

			//cout << "generation" << gen << ", ok after d+=?\n";
		// sampling the next generation:

		for (i = 0; i < dv; i++) // for each deme i
		{
		    nb3 = Nv*i;

		    // number of immigrants in deme i:

		    nbMig = int(binldev(migv, Nv));

		    // migrant individuals:

		    for (ind = 0; ind < nbMig; ind++)
		    {
			nb4 = nb3 + ind;

			// selection of a deme of origin (j):

			if ((gen < T1v) || (gen > Tcontact))  // no barrier
			{
			    if ((i > 0) && (i < dv-1))
			    {
				if (rnd.rand() < 0.5)
				    j = i - 1;
				else
				    j = i + 1;
			    }
			    else
			    {
				if (i == 0)
				    j = 1;
				else if (i == dv - 1)
				    j = dv - 2;
			    }
			}

			   
			nb2 = Nv*j;

			// sampling parents:

			do
			{
			    chr1 = nb2 + int(rnd.randInt(N_1));

			} while (rnd.rand()> (Wtot[chr1]/Wmax[j]));

			do
			{
			    chr2 = nb2 + int(rnd.randInt(N_1));

			} while (rnd.rand()> (Wtot[chr2] / Wmax[j]));

			// recombination

					if (withrec)
						rec(temp[nb4], pop[chr1], pop[chr2], Lv, nbSv);
					else
						freerec(temp[nb4], pop[chr1], pop[chr2], nbSv);
		    }

		    // philopatric individuals:

		    for (ind = nbMig; ind < Nv; ind++)
		    {
			nb4 = nb3 + ind;

			// sampling parents:

			do
			{
			    chr1 = nb3 + int(rnd.randInt(N_1)) ;

			} while (rnd.rand()> (Wtot[chr1] / Wmax[i]));

			do
			{
			    chr2 = nb3 + int(rnd.randInt(N_1));

			} while (rnd.rand()> (Wtot[chr2] / Wmax[i]));

			// recombination

					if (withrec)
						rec(temp[nb4], pop[chr1], pop[chr2], Lv, nbSv);
					else
						freerec(temp[nb4], pop[chr1], pop[chr2], nbSv);
		    }
		}

		// mutation

	       // for (i = 0; i < Nd; i++)
		//{
		 //   mut = int(poisdev(Uv)); // number of new mutations
		  //  for (j = 0; j < mut; j++)
		   //     temp[i].sel.flip(int(rnd.randInt(nbS_1)));
		//}

		// update population:

		cp = pop;
		pop = temp;
		temp = cp;


			// write result in HDF5 file
		if (gen % pasv == 0)
		{
		indexGen = gen / pasv;
				// Allele frequency data
				for (loc = 0; loc < nbSv; loc++)
				{
					for (i = 0; i < dv; i++)
					{
						nb = Nv * i;
						d = 0;
						for (j = 0; j < Nv; j++)
						{
							if (pop[nb + j].sel[loc] == 1)
							{
								d += 1;
								sdata_freq[loc * dv + i] = d / Nv;
							}
						}
					}
				}
				writeTimeStepHDF5(file, dset_freq, RANK_FREQ,
					dim_sub_freq, sdata_freq, indexGen);

		



			delta_p = sdata_freq[bv] - p_old;
                        p_old = sdata_freq[bv];
                        //cout << ", delta_p:" << delta_p << ", p:" << p_old;

                        if (accGen > 200) //burnin period
                        {
                                if (delta_p ==0)
                                        nbSign += 1;
                                else if (delta_p < 0 && sign ==true)
                                        {
                                        nbSign +=1;
                                        sign = false;
                                        }
                                else if (delta_p > 0 && sign==false)
                                        {
                                        nbSign+=1;
                                        sign=true;
                                        }
                         }

				// mean fitness of demes
				for (i = 0; i < dv; i++) // demes
				{
					sdata_w[i] = 0;
					for (j = 0; j < Nv; j++) // individuals
					{
						 sdata_w[i] += Wtot[i * Nv + j];
					}
					sdata_w[i] /= Nv;
				}
				writeTimeStepHDF5(file, dset_w, RANK_W,
					dim_sub_w, sdata_w, indexGen);

		//	cout << " ok after writing fitness in HDF5?\n";

				// which generation was saved
				savedGen[gen] = accGen;
			}
			accGen +=	1;

		} // end gen loop
	//	cout << "ok after end gen loop?\n";
	
		if (nbSign > 200)
		{
			if (last ==0)
				last = 1;
			else
				equi=true;
		 }
	}
		// write savedGen to file
	writeGenSaved(file, dset_gen, savedGen);
//	cout << "ok after writing gen in HDF5?\n";

	dset_freq.close();
	dset_w.close();
	dset_gen.close();
	file.close();

	fin = time(0);


	// time length:
	int temps = int(difftime(fin, debut));
	cout << NbGen << " generations ont pris " << temps << " secondes\n";	


	delete [] pop;
	delete [] temp;
	delete [] mutations;
	delete [] Wtot;
    delete [] Wmax;
	delete [] m;
	delete [] v;
}
