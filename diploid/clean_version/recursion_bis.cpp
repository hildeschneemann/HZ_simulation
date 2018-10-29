#include "fisher.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

#define MAXGEN 100000

extern FILE * fichierE;
extern FILE * fichierS;
extern MTRand rnd;

double getfitness(double HI, double p12, double kd2, double a, double f, double p, double r)
{
	double sz2;
	sz2 = p2 + (4.0 - 2.0*p2)*4.0*HI*(1.0-HI)+(f-1.0)*p12+(f-1.0+r)*p12*(1.0-p12)
	return(exp(-a * pow(sz2,kd2)));
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
fv: fitness breakdown in F1 hybrids
pv: fitness breakdown in parents
rv: contribution of heterozygosity to curvature

* Simulation parameters
tv: number of generations that are saved and after which the
	equilibrium state is checked
----------------------------------------------------------*/


void recursion(	int Dv, int Nv, double mv, 
				int lv, double Lv,
				double kv, double av, double fv, double pv, double rv,
				int tv)

{
	// variables:

	int i, j, k, loc, gen, par1, par2, ind, nb, nb2, nb3, nb4, nbMig;
	double w, wbar, varw, d, HI, p12, p2, delta_p, p_old;
	bool withrec;
	vector<int> store;

	// various fixed quantities:
	
    int ND = Nv * Dv;
    int twoN = 2*Nv;
    int twoND = 2*Nv*Dv;
	int N_1 = Nv - 1;
	double kd2 = kv / 2.0;
	int tv_1 = tv - 1;
	int ND1 = twoN * bv;
	bool sign = true;
	bool equi = false;
	bool last = 0;
	int nbSign = 0;
	withrec = Lv == -1 ? false : true;
	int round = 0;	
	int nbW = lv * (lv + 1) / 2.0;
	int indexW = 0;

	boost::dynamic_bitset<> tmp1;
	boost::dynamic_bitset<> tmp2;
	boost::dynamic_bitset<> tmp3;
	boost::dynamic_bitset<> tmp4;

	
	// HDF5 constants and variables

	string fileName;
	stringstream nameF;
	nameF << "res_D" << Dv << "_N" << Nv << "_m" << mv << 
		"_l" << lv << "_L" << Lv << 
		"_k" << kv << "_a" << av << "_f" << fv << "_p" << pv << "_r" << rv <<
		"_t" << tv << ".txt";
	nameF >> fileName;

	int accGen = 0; //which generations are recorded
	
    chr * pop = new chr [twoND];
	chr * temp = new chr [twoND]; // used to generate next generation
    chr * cp;

	// "Wtot" will hold the fitness of each individual:
	double * Wtot = new double [ND];
    
    // maximal fitness in each deme:
    double * Wmax = new double [Dv];
	
	// "HIend" will hold the hybrid index of each individual at the last generation
	
	double * HIend = new double [ND];
	
	// "Wtable" will be the lookup table for the fitness corresponding to a genotype
	
	double * Wtable = new double [nbW];
	
	// "Wsave" will save the number of individuals with a certain fitness value to output
	
	double * Wsave = new double [nbW][Dv];
	
	// create HDF5 file and open it

	time_t debut, fin;
	struct tm *ptr;
	debut = time(0);
	
	// initialization: allele 0 is fixed at all selected loci:

	for (i = 0; i < twoND; i++)
    {
		pop[i].sel.resize(lv);
        temp[i].sel.resize(lv);
		if(i <= ND1) //complete divergence
			pop[i].sel.flip();
			temp[i].sel.flip();
    }	
	
	
	// create fitness lookup table
	for (i =0; i < nbW; i++)
	{
		p2 = nbW - i;
		for (j=0; j < (nbW - i); j++)
		{	
			p12 = j;
			Wtable[indexW] = getfitness( p2, p12, kd2, a,  f,  p, r);
			indexW += 1;
		}
	}
	
	
	
   
    // generations:
	while(equi == false & accGen < MAXGEN)
	{
		for (gen=0; gen < tv; gen++)
		{
		
		  
			for (i = 0; i < Dv; i++)  // for each deme
				{
					nb = i * Nv;
					Wmax[i] = 0;
				
					for (j = 0; j < Nv; j++)  // for each individual
					{
						nb2 = 2 * (nb + j);
				  
						// determine hybrid index and heterozygosity
						tmp1 = (pop[nb2].sel ^ pop[nb2+1].sel);
						p12 = tmp1.count(); 
						//p12 /= lv; //creates a vector with a 0 for homozygote and 1 for heterozygote sites
						tmp2 = (pop[nb2].sel & pop[nb2+1].sel);
						p2 = tmp2.count();
						//p2 /= lv;
						//HI = p2 + p12/2.0; 
							
						// fitness
					
						w = W[ ((p2^2 - p2) / 2) + p12 ];
						
						
						Wtot[nb2/2] = w;
							
						if (gen == tv_1)
							HIend[nb2/2] = HI;

						if (Wmax[i] < w)
							Wmax[i] = w;
						
					}
				}
			
			// sampling the next generation:
			
			for (i = 0; i < Dv; i++) // for each deme i
			{
				nb3 = Nv*i;
				
				// number of immigrants in deme i:
				
				nbMig = int(binldev(migv, Nv));
				
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
			}
			

			// update population:

			cp = pop;
			pop = temp;
			temp = cp;
				
			
			if (gen % pasv == 0)
			{
				indexGen = gen / pasv; 
				
				Wtot = sort(Wtot.begin(), Wtot.end()),
			
				for (i = 0; i < Dv; i++)
				{
					sumk =0;
					
					for ( j =0; j < nbW; j++)
					{
						k = 0;
						do
						{
							k += 1;
						}while (Wtot[sumk+1] == Wtot[sumk])
							
						Wsave[j][i] = k;	
					}
				}
		
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
							{	
								d += 1;
								sdata_freq[loc * Dv + i] = d / twoN;
							}
						}

					}					

				}

				
				//count how often the sign of the change in allele frequency changes
				delta_p = sdata_freq[bv] - p_old;
				p_old = sdata_freq[bv];

	 		
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
							last = 1;
					else 
					{
						nbSign+=1;
						sign=true;
					}
				}

				
				// which generation was saved
				
				accGen;
			}	

			accGen +=	1;
				
		} // end gen loop

		if (nbSign > 100)
		{
			if (last ==0) 
				last = 1;
			else 
				equi=true;			
		}

	round +=1;
	}



	fin = time(0);
	int temps = int(difftime(fin, debut));
	cout << tv << " generations ont pris" << temps << "secondes\n";


	
	delete [] pop;
	delete [] temp;
	delete [] Wtot;
    delete [] Wmax;
}
