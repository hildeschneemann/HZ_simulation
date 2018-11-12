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
#define BURNIN 1000
#define SIGNCHANGES 100


extern FILE * fichierE;
extern FILE * fichierS;
extern MTRand rnd;

double getfitness(double HI, int p12, int lv, double kd2, double a, double f, double p, double r)
{
	double sz2;
	double p12prop = double(p12) / double(lv);
	sz2 = p + (4.0 - 2.0*p)*HI*(1.0-HI)+(f-1.0)*p12prop+(f-1.0+r)*p12prop*(1.0-p12prop);
	return(exp(-a * pow(sz2,kd2)));
}

int getfitness_index(int p2, int p12, int lv)
{
	int notp2 = lv - p2;
	return((notp2 * (notp2+1))/ 2 + p12);
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
pasv: generation steps
----------------------------------------------------------*/


void recursion(	int Dv, int Nv, double mv, 
				int lv, double Lv,
				double kv, double av, double fv, double pv, double rv,
				int tv, int pasv)

{
	// variables:

	int i, j, k, loc, gen, par1, par2, ind, nb, nb2, nb3, nb4, nbMig, p2, p12;
	double w, wbar, varw, d, delta_p, p_old, HI;
	bool withrec;
	vector<int> store;

	// various fixed quantities:
	
    int ND = Nv * Dv;
    int twoN = 2*Nv;
    int twoND = 2*Nv*Dv;
	int N_1 = Nv - 1;
	double kd2 = kv / 2.0;
	int tv_1 = tv - 1;
	int bv = Dv / 2;
	int ND1 = twoN * bv;
	bool sign = true;
	bool equi = false;
	bool last = 0;
	int nbSign = 0;
	withrec = Lv == -1 ? false : true;
	int round = 0;	
	int nbW = (lv + 2) * (lv + 1) / 2.0;
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
	
	double * Wsave = new double [nbW];
	
	// create HDF5 file and open it

	time_t debut, fin;
	struct tm *ptr;
	debut = time(0);


	ofstream fout;
	
	// initialization: allele 0 is fixed at all selected loci:

	for (i = 0; i < twoND; i++)
    {
		pop[i].sel.resize(lv);
        temp[i].sel.resize(lv);
		if(i < ND1) //complete divergence
			pop[i].sel.flip();
			temp[i].sel.flip();
    }	
	
	
	// create fitness lookup table
	for (i = lv; i > -1; i--)
	{
		p2 = i;
		for (j=0; j < (lv - i + 1); j++)
		{	
			p12 = j;
			HI = (p2 + p12 / 2.0) / lv;
			indexW = getfitness_index(p2, p12, lv);
			Wtable[indexW] = getfitness(HI, p12, lv, kd2, av,  fv,  pv, rv);
			cout << indexW << "\tp2: " << p2 << "\tp12: " << p12 << "\tHI: " << HI << "\tfitness: " << Wtable[indexW] << "\n";
			//indexW += 1;
		}
	}
	
	
	
   
    // generations:
	while(equi == false & accGen < MAXGEN)
	{
		if (last == 1)
		{
			fout.open(fileName);

			// Make header output file
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
				for (j = 0; j < nbW; j++)
				{
					fout << "d" << i << "W" << Wtable[j] << ", ";
				}
			}
			
			fout << endl;
		}

		cout << "round: " << round << "\tnbSign: " << nbSign << "\tlast?: " << last << "\taccGen: " << accGen << "\n";
	
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
						indexW = getfitness_index(p2, p12, lv);	
						w = Wtable[indexW];
						
						
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
			}
			

			// update population:

			cp = pop;
			pop = temp;
			temp = cp;
				
			
			if (gen % pasv == 0)
			{
				if (last == 1)
				{	
				
					fout << gen << ", ";	
		
					// Allele frequency data

					for (loc = 0; loc < lv; loc++)
					{
						//fout << gen << "," <<loc << ",";

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
						//fout << endl;

					}
					
					// Number of individuals with given fitness value per deme
					for (i = 0; i < Dv; i++)
					{
					
						for ( j =0; j < nbW; j++)
						{	
							int histW = 0;

							for ( ind = 0; ind < Nv; ind++)
							{
								if ( Wtot[ind + i * Nv] == Wtable[j])
								{
									histW += 1;
								}
							} 	
							
							fout << histW << ",";
						}	
					}
					fout << endl;
				}

				
				else if (last == 0)
				{
					d=0;
					for (j =0; j < twoN; j++)
					{
						if (pop[ND1 + j].sel[0] == 1)
							d += 1;
					}
					delta_p = d / twoN - p_old;
					p_old = d / twoN;
				}	

	 		
				if (accGen > BURNIN) //burnin period
				{
					if (delta_p ==0)
						nbSign += 1;
					else if (delta_p < 0 && sign ==true)
						{
						nbSign +=1;
						sign = false;
						}
					else if (delta_p > 0 && sign==false)
						nbSign +=1;
						sign=true;
				}
				
			}	

			accGen +=1;
				
		} // end gen loop

		if (nbSign > SIGNCHANGES)
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
