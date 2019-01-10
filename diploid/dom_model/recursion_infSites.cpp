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


double getfitness(double sz2, double kd2, double a)
{
        //double sz2 = d * d;
        return(exp(-a * pow(sz2,kd2)));
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
Nv: size of the population

* Genomic parameters
lv: total number of selected loci per genome
Lv: genome map length (mean nb of cross-overs per meiosis)
        set to -1 for free recombination
Uv: mutation rate

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

void recursion( int Nv,
                                double Lv, double Uv,
                                double kv, double av, int nv, double Fv, double qv,
                                int Rv, int tv, int repv, bool newmut)

{
        // variables:

        int i, j, k, loc, gen, par1, par2, ind, nb, nb2, nb3, nb4, nbMig, p2, p12, mut, mut_j, indmut;
        double w, wbar, d, HI, sz2;
        bool withrec;
        vector<int> store;

        // various fixed quantities:
        int lv = 0;
		int MaxGen = (Rv-1)*tv + BURNIN + 1;
        int twoN = 2*Nv;
        int N_1 = Nv - 1;
        double kd2 = kv / 2.0;
        int tv_1 = tv - 1;
        withrec = Lv == -1 ? false : true;
        int l1 = lv - 1;

        boost::dynamic_bitset<> tmp1;
        boost::dynamic_bitset<> tmp2;
        boost::dynamic_bitset<> tmp3;
        boost::dynamic_bitset<> tmp4;


        string fileName;
        stringstream nameF;
        nameF << "res_N" << Nv <<
                        "_L" << Lv << "_U" << Uv <<
                        "_k" << kv << "_a" << av << "_n" << nv << "_F" << Fv << "_q" << qv <<
                        "_rep" << repv << ".txt";
        nameF >> fileName;

        string fileName2;
        stringstream nameF2;
        nameF2 << "mutEff_N" << Nv <<
                        "_L" << Lv << "_U" << Uv <<
                        "_k" << kv << "_a" << av << "_n" << nv << "_F" << Fv << "_q" << qv <<
                        "_rep" << repv << ".txt";
        nameF2 >> fileName2;


        int accGen = 0; //which generations are recorded

        chr * pop = new chr [twoN];
        chr * temp = new chr [twoN]; // used to generate next generation
        chr * cp;

        // "Wtot" will hold the fitness of each individual:
        double * Wtot = new double [Nv];

        vector<double> mutations;
        vector<double> delta;
	vector<int> muttimes;
	vector<int> fixtimes;	

        time_t debut, fin;
        struct tm *ptr;
        debut = time(0);


        ofstream fout;
        ofstream fout2;

        // initialization: allele 0 is fixed at all selected loci:

        for (i = 0; i < twoN; i++)
        {
                pop[i].sel.resize(lv);
                temp[i].sel.resize(lv);
        }

        // generations:
        for(gen=0; gen < MaxGen; gen++)
        {

                double Wmax = 0;
                for (j = 0; j < Nv; j++)  // for each individual
                {
                                nb2 = 2 * j;
                                sz2 = 0;
                                for (loc = 0; loc < lv; loc++)
                                {
                                                //cout << pop[nb2].sel[loc] << pop[nb2+1].sel[loc];
						d = 0;
                                                for (k = 0; k < nv; k++)  // phenotypes of the individual
                                                {
                                                                if (pop[nb2].sel[loc] == 1 & pop[nb2+1].sel[loc] == 1)
                                                                                d += mutations[nv * loc + k];
                                                                else if(pop[nb2].sel[loc] == 1 || pop[nb2+1].sel[loc] == 1)
                                                                                d += mutations[nv * loc + k] * delta[nv * loc + k];
                                                }
                                                sz2 += d * d; // "sz2" is the square of the distance to the optimum
                                }
				//cout << endl;
                                // fitness:
                                w = getfitness(sz2, kd2, av);
                                Wtot[nb2/2] = w;
                                if (Wmax < w)
                                                Wmax = w;
                }


                // sampling the next generation:

                for (ind = 0; ind < Nv; ind++)
                {
                                // sampling first parent:
                                do
                                {
                                                par1 = int(rnd.randInt(N_1)) ;
                                } while (rnd.rand()> (Wtot[par1] / Wmax));
                                // recombination
                                if (withrec)
                                                rec(temp[2*ind], pop[2*par1], pop[2*par1+1], Lv, lv);
                                else
                                                freerec(temp[2*ind], pop[2*par1], pop[2*par1+1], lv);
                                do
                                {
                                                par2 = int(rnd.randInt(N_1));
                                } while (rnd.rand()> (Wtot[par2] / Wmax));
                                // recombination
                                if (withrec)
                                                rec(temp[2*ind+1], pop[2*par2], pop[2*par2+1], Lv, lv);
                                else
                                                freerec(temp[2*ind+1], pop[2*par2], pop[2*par2+1], lv);
                }

        // mutation
        	//cout << "mut, trait, mutational effect, delta\n";
		mut = int(poisdev(Uv)) * twoN; // number of new mutations
		for (j = 0; j < mut; j++)
		{
			indmut = int(rnd.randInt(twoN));
			//cout << "mutated individual: " << indmut << "\n";
			for (i = 0; i < twoN; i++) //loop through individuals
			{
                		if ( i == indmut)
				{
					temp[i].sel.push_back(1);
				}
				else
					temp[i].sel.push_back(0);

            }
			
			for (k = 0; k < nv; k++) // generate mutational effect and dominance of new mutation for each trait
			{
				mutations.push_back(gasdev());
				delta.push_back(getdominance(qv, Fv));
				//cout << lv + j << " , " << k << " , " << mutations[nv * (lv+j) + k] << " , "<< delta[nv * (lv + j) + k] << '\n';
			}
			muttimes.push_back(gen);
			fixtimes.push_back(0);
        }
		lv += mut;
		//cout << "new lv:" << lv << "\n";


                // update population:
                cp = pop;
                pop = temp;
                temp = cp;

		// Allele frequency data and fixation times
		for (loc = 0; loc < lv; loc++)
                {

                     	d = 0;
                     	for (j = 0; j < twoN; j++)
                     	{	
				if (pop[j].sel[loc] == 1)
				d += 1;
			}
			if ((d / twoN) == 1 & fixtimes[loc] == 0)
				fixtimes[loc] = gen;
	         }
		
		if (gen >= BURNIN & (gen-BURNIN) % tv == 0)    //display progress of simulation           
 			fprintf(stderr,"%d %d\n",gen,MaxGen);

        } // for(gen )... The main gen loop
	
	// output mutational effects, dominance, fixation time and position of FIXED mutations to file
	fout.open(fileName);
	fout << "position, mut_time, fix_time, ";
	for (i = 0; i < nv; i++)
	{
		fout << "mut_effect" << i << ", ";
	}
	fout << "dom\n";

	for (loc = 0; loc < lv; loc++)
	{

		d = 0;
		for (j = 0; j < twoN; j++)
		{
			if (pop[j].sel[loc] == 1)
			d += 1;
		}
		if ((d / twoN) == 1 )
		{
			if (fixtimes[loc] == 0)
				fixtimes[loc] = gen;
			fout << loc << ", " << muttimes[loc] << ", " << fixtimes[loc] << ", ";
			for (i = 0; i < nv; i++)
			{
				fout << mutations[loc * nv + i] << ", ";
			}
			fout << delta[loc] << endl;
		}
	 }	
	fout.close();

	fout2.open(fileName2); // output mutational effects of ALL mutations to file
        fout2 << "trait\tloc\tmut_effect\tdom\n";
	for (j = 0; j < lv; j++)
        {
		for (i = 0; i < nv; i++)
		{
                	nb = nv * j;
			fout2 << i << "\t" << j << "\t" << mutations[nb + i]  << "\t" << delta[nb + i] << "\n";
           	}
        }
        fout2.close();
		
        fin = time(0);
        int temps = int(difftime(fin, debut));
        cout << gen << " generations ont pris " << temps << " secondes\n";

        delete [] pop;
        delete [] temp;
        delete [] Wtot;

}
