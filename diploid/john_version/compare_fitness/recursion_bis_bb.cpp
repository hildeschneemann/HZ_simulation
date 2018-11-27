#include "fisher.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <boost/dynamic_bitset.hpp>
#include <boost/math/distributions.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

#define BURNIN 2

extern FILE * fichierE;
extern FILE * fichierM;
extern MTRand rnd;

double readmut(double * mutations, int j)
{
        double mut = 0.0;
        fscanf(fichierM, "%lf", &mut);
        fscanf(fichierM, "%lf", &mut);
        fscanf(fichierM, "%lf", &mutations[j]);
	fscanf(fichierM, "%lf", &mut);
	cout << "muteff: " << mutations[j] << "\n";
}

double getfitness(double sz2, double kd2, double a)
{
	//double sz2 = d * d;
	return(exp(-a * pow(sz2,kd2)));
}

double getfitness_simple(double HI, int p12, int lv, double kd2, double a, double f, double p, double r)
{
        double sz2;
        double p12prop = double(p12) / double(lv);
        sz2 = p + (4.0 - 2.0*p)*HI*(1.0-HI)+(f-1.0)*p12prop+(f-1.0-r)*p12prop*(1.0-p12prop);
        return(exp(-a * pow(sz2,kd2)));
}

int getfitness_index(int p2, int p12, int lv)
{
        int notp2 = lv - p2;
        return((notp2 * (notp2+1))/ 2 + p12);
}

double getdominance(double qv, double Fv)
{
    if (Fv == 0)
		return(qv);
		
	if (Fv == 1)
		return(rnd.rand() < qv ? 1 : 0);
	
	if (Fv == -1)
		return(rnd.rand());
	
	if(Fv > 0 & Fv < 1)
	{
		double alpha = qv * (1.0-Fv)/Fv;
		double beta = (1.0-Fv)*(1.0-qv)/Fv;
		boost::math::beta_distribution<double> betadist(alpha, beta);
		return(quantile(betadist, rnd.rand()));
	}
   
    exit(0);   
}

/*----------------------------------------------------------

* Demographic Parameters:
Dv: number of demes
Nv: size of the population in each deme
mv: migration rate between adjacent demes

* Genomic parameters
lv: total number of selected loci per genome
Lv: genome map length (mean nb of cross-overs per meiosis)
	set to -1 for free recombination

* FGM parameters
kv: curvature of fitness function (cf Gros et al 2009 Genetics)
av: strength of selection (steepness of fitness function)
nv: number of phenotypic dimensions (complexity)
Fv: parameter of distribution of phenotypic dominance distribution (variance within/between loci)
qv: parameter of phenotypic dominance distribution

* Simulation parameters
Rv: number of records
tv: number of generations between each record
----------------------------------------------------------*/

void recursion(	int Dv, int Nv, double mv, 
				int lv, double Lv,
				double kv, double av, int nv, double Fv, double qv,
				int Rv, int tv, int repv, bool newmut)

{
	// variables:

	int i, j, k, loc, gen, par1, par2, ind, nb, nb2, nb3, nb4, nbMig, p2, p12, geno;
	double w, wbar, d, HI, sz2, w2;
	bool withrec;
	vector<int> store;

	// various fixed quantities:
	int MaxGen = (Rv-1)*tv + BURNIN + 1;
	int ND = Nv * Dv;
    int twoN = 2*Nv;
    int twoND = 2*Nv*Dv;
	int N_1 = Nv - 1;
	double kd2 = kv / 2.0;
	int tv_1 = tv - 1;
	int bv = Dv / 2;
	int ND1 = twoN * bv;
	withrec = Lv == -1 ? false : true;
	int nbW = (lv + 2) * (lv + 1) / 2.0;
	int indexW = 0;
	double fv = Fv;
	double pv = 0;
	double rv = 1.0;

	boost::dynamic_bitset<> tmp1;
	boost::dynamic_bitset<> tmp2;
	boost::dynamic_bitset<> tmp3;
	boost::dynamic_bitset<> tmp4;


	string fileName;
	stringstream nameF;
	nameF << "res_D" << Dv << "_N" << Nv << "_m" << mv << 
		"_l" << lv << "_L" << Lv << 
		"_k" << kv << "_a" << av << "_n" << nv << "_F" << Fv << "_q" << qv <<
		"_rep" << repv << ".txt";
	nameF >> fileName;

	string fileName2;
	stringstream nameF2;
	nameF2 << "mutEff_D" << Dv << "_N" << Nv << "_m" << mv << 
		"_l" << lv << "_L" << Lv << 
		"_k" << kv << "_a" << av << "_n" << nv << "_F" << Fv << "_q" << qv <<
		"_rep" << repv << ".txt";
	nameF2 >> fileName2;


	int accGen = 0; //which generations are recorded
	int poplength = 2 * pow(3,lv);	
    //chr * pop = new chr [poplength];
	chr * pop = new chr [twoND];
	chr * temp = new chr [twoND];
	//chr * temp = new chr [poplength]; // used to generate next generation
    chr * cp;


	double * mutations = new double [lv * nv];
	double * delta = new double [lv * nv];

	// "Wtot" will hold the fitness of each individual:
	double * Wtot = new double [ND];
    
    // maximal fitness in each deme:
    double * Wmax = new double [Dv];
	
	// "Wtable" will be the lookup table for the fitness corresponding to a genotype
        double * Wtable = new double [nbW];

	double * Brown = new double [lv + 1];
	
	time_t debut, fin;
	struct tm *ptr;
	debut = time(0);

       // create fitness lookup table
        for (i = lv; i > -1; i--)
        {
                p2 = i;
                for (j=0; j < (lv - i + 1); j++)
                {
                        p12 = j;
                        HI = (p2 + p12 / 2.0) / lv;
                        indexW = getfitness_index(p2, p12, lv);
                        Wtable[indexW] = getfitness_simple(HI, p12, lv, kd2, av,  fv,  pv, rv);
                 }
         }
	ofstream fout;
	ofstream fout2;
	
	// initialization: allele 0 is fixed at all selected loci:

	//for (i = 0; i < poplength; i++)
    	for (i = 0; i < twoND; i++)
	{
		pop[i].sel.resize(lv);
        temp[i].sel.resize(lv);
		if(i < ND1) //complete divergence
		{
			pop[i].sel.flip();
			temp[i].sel.flip();
		}
    }

	//create mutational effects matrix and dominance matrix
	//initiate delta
	fout2.open(fileName2);	
	for (i = 0; i < nv; i++)
	{
		nb = lv * i;
		brownian_bridge(Brown, lv+1, nv);
		for (j = 0; j < lv; j++)
		{
			if(newmut == true)
				mutations[nb + j] = Brown[j+1] - Brown[j];
			else
				readmut(mutations, (nb+j));
			
			delta[nb + j] = getdominance(qv,Fv);
			fout2 << i << "\t" << j << "\t" << mutations[nb + j]  << "\t" << delta[nb + j] << "\n";
		}
	}
	fout2.close();	
	
	
	// Open output file and add header
	fout.open(fileName);
	fout << "gen" << ", ";
	for(i= 0; i < lv; i++)
	{
		for(j = 0; j < Dv; j++)
		{	
			fout << "d" << j << "l" << i << ", ";
		}
	}
	for (i = 0; i < Dv; i++)
	{
		fout << "d" << i << "W" << ", ";
	}
	fout << endl;

//generate systematic genotypes to test fitness functions

        for (int loc = 0; loc < lv; loc++)
        {
                geno = 0;
                for (int t = 0; t < pow(3,loc); t++)
                {
                        for (int i = 0; i < 3; i++)//3 states: 0 homo P1, 1 hetero, 2 homo P2
                        {
                                for (int j = 0; j < pow(3, (lv - loc - 1)); j++)
                                {
                                        //cout << i;
                	                if (i == 1)
                        	        {
                                	        pop[2*geno].sel[loc].flip();
                               		}
                                	if (i == 2)
                                	{
                                        	pop[2 * geno].sel[loc].flip();
                                        	pop[2 * geno + 1].sel[loc].flip();
                                	}
                                        geno += 1;
                                }
                        }
                }
                //cout << "\n";
        }

for (int ind = 0; ind < pow(3,lv); ind++) //for each individual
{
	d = 0.0;
	p12 = 0;
	p2 = 0; 
	sz2 = 0.0;
	cout << "ind: " << ind << ", genotype: ";
	for (int n = 0; n < lv; n++) //for each trait
	{
		d = 0.0;
		p2 = 0;
		p12 = 0;
		for(i = 0; i < lv; i++) // for each locus
		{
		
			if(pop[2 * ind].sel[i] ^ pop[2 * ind + 1].sel[i])
			{
				d += mutations[lv * n + i] * delta[lv * n + i];
				p12 += 1;
			}
			if(pop[2 * ind].sel[i] & pop[2 * ind + 1].sel[i])
			{
				d += mutations[lv * n + i];
				p2 += 1;
			}
			if (n == 0)
				cout << pop[2 * ind].sel[i] << pop[2 * ind + 1].sel[i];
		}
		sz2 += d * d;
	}
	w = getfitness(sz2, kd2, av);
	indexW = getfitness_index(p2,p12, lv);
	w2 = Wtable[indexW];
	cout << ", full: " << w << ", simple: " << w2 << ", p2: " << p2 << ", p12: " << p12 << "sz2: " << sz2 << "\n";
}
	

	//backcrosses(temp, mutations, lv, nv, Lv, kd2, av, twoND, nb4, withrec);
	// generations:
	for(gen=0; gen < MaxGen; gen++)
	{
		//cout << "round: " << round << "\tnbSign: " << nbSign << "\tlast?: " << last << "\taccGen: " << accGen << "\n";
		for (i = 0; i < Dv; i++)  // for each deme
		{
		    nb = i * Nv;
		    Wmax[i] = 0;
		    for (j = 0; j < Nv; j++)  // for each individual
			{
				nb2 = 2 * (nb + j);
				sz2 = 0.0;
				for (k = 0; k < nv; k++)  // phenotypes of the individual
				{
					d = 0.0;
					for (loc = 0; loc < lv; loc++)
					{
						if (pop[nb2].sel[loc] == 1 & pop[nb2+1].sel[loc] == 1)
							d += mutations[lv * k + loc];
						else if(pop[nb2].sel[loc] == 1 ^ pop[nb2+1].sel[loc] == 1)
							d += mutations[lv *k + loc] * delta[lv * k + loc];
						if (k == 0)
							cout << pop[nb2].sel[loc] << pop[nb2+1].sel[loc];
					}
					sz2 += d * d; // "sz2" is the square of the distance to the optimum
				}

				// fitness:
				w = getfitness(sz2, kd2, av);
				Wtot[nb2/2] = w;
				if (Wmax[i] < w)
					Wmax[i] = w;

				//simple fitness
				tmp1 = (pop[nb2].sel ^ pop[nb2+1].sel);
                                p12 = tmp1.count();
				tmp2 = (pop[nb2].sel & pop[nb2+1].sel);
                                p2 = tmp2.count();
				indexW = getfitness_index(p2, p12, lv);
                                w2 = Wtable[indexW];
				cout << ", full: " << w << ", simple: " << w2 << ", p12: " << p12 << ", p2: " << p2 << ", gen: " << gen << "\n";
		    }
		}
			
		// sampling the next generation:
		for (i = 0; i < Dv; i++) // for each deme i
		{
			nb3 = Nv*i;
			// number of immigrants in deme i:
			nbMig = int(binldev(mv, Nv));
			// migrant individuals:
			for (ind = 0; ind < nbMig; ind++)
			{
				nb4 = nb3 + ind;
				// selection of a deme of origin (j):
				if ((i > 0) && (i < Dv-1))
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
					else if (i == Dv - 1)
						j = Dv - 2;
				}
				nb2 = Nv*j;
				// sampling first parent:
				do
				{
					par1 = nb2 + int(rnd.randInt(N_1));
				} while (rnd.rand()> (Wtot[par1]/Wmax[j]));
				// recombination
				if (withrec)
					rec(temp[2*nb4], pop[2*par1], pop[2*par1+1], Lv, lv);
				else
					freerec(temp[2*nb4], pop[2*par1], pop[2*par1+1], lv);
				// sampling second parent:
				do
				{
					par2 = nb2 + int(rnd.randInt(N_1));
				} while (rnd.rand()> (Wtot[par2] / Wmax[j]));
				// recombination
				if (withrec) 
					rec(temp[2*nb4+1], pop[2*par2], pop[2*par2+1], Lv, lv);
				else
					freerec(temp[2*nb4+1], pop[2*par2], pop[2*par2+1], lv);
			}	
			// philopatric individuals:
			for (ind = nbMig; ind < Nv; ind++)
			{               
				nb4 = nb3 + ind;
				// sampling first parent:
				do
				{
					par1 = nb3 + int(rnd.randInt(N_1)) ;
				} while (rnd.rand()> (Wtot[par1] / Wmax[i]));
				// recombination
				if (withrec)
					rec(temp[2*nb4], pop[2*par1], pop[2*par1+1], Lv, lv);
				else
					freerec(temp[2*nb4], pop[2*par1], pop[2*par1+1], lv);
				do
				{
					par2 = nb3 + int(rnd.randInt(N_1));
				} while (rnd.rand()> (Wtot[par2] / Wmax[i]));
				// recombination
				if (withrec)
					rec(temp[2*nb4+1], pop[2*par2], pop[2*par2+1], Lv, lv);
				else
					freerec(temp[2*nb4+1], pop[2*par2], pop[2*par2+1], lv);
			}
		} // Sampling next generation
			
		// update population:
		cp = pop;
		pop = temp;
		temp = cp;
		
		// Record data to output file..
		if (gen >= BURNIN & (gen-BURNIN) % tv == 0)
		{
			fout << gen << ", ";	
			// Allele frequency data
			for (loc = 0; loc < lv; loc++)
			{
				for (i = 0; i < Dv; i++)
				{
					nb = twoN * i;
					d = 0;
					for (j = 0; j < twoN; j++)
					{
						if (pop[nb + j].sel[loc] == 1)
							d += 1;
					}
					fout << d / twoN << ",";
				}
			}
			// Number of individuals with given fitness value per deme
			for (i = 0; i < Dv; i++)
			{
				wbar = 0;
				for ( j =0; j < Nv; j++)
				{	
					wbar += Wtot[j];				
				}	
				wbar /= Nv;
				fout << wbar << ",";	
			}
			fout << endl;
		} // End record loop
		fprintf(stderr,"%d %d\n",gen,MaxGen);	
	} // for(gen )... The main gen loop


	fin = time(0);
	int temps = int(difftime(fin, debut));
	cout << gen << " generations ont pris " << temps << " secondes\n";
	
	fout.close();
	delete [] pop;
	delete [] temp;
	delete [] Wtot;
    delete [] Wmax;
}
